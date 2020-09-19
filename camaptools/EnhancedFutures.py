#!/usr/bin/env python3

import functools

class EnhancedFuture(object):
    def __init__(self, executor, kind):
        self.executor = executor
        self.kind = kind

    def submit(self, func, *args):
        if self.kind == 'thread':
            p = functools.partial(func, *args)
        elif self.kind == 'process':
            p = self.executor.pex.submit(func, *args)
        self.process = p

    def result(self):
        if self.kind == 'thread':
            if self.process in self.executor.results:
                return self.executor.results[self.process]
            else:
                res = self.process()
                self.executor.results[self.process] = res
                return res
        elif self.kind == 'process':
            return self.process.result()


class EnhancedProcessPoolExecutor(object):
    def __init__(self, max_workers=0, use_threads=True, mpi=False):

        self.workers = max_workers
        self.use_threads = use_threads
        try:
            assert self.workers or self.use_threads
        except AssertionError:
            print('WARNING: Specifying no workers and no threads is not allowed, forcing use_threads=True.')
            self.use_threads = True
        #self.tex = ThreadPoolExecutor(max_workers=1) # MIXING threads and processes DOESN'T WORK

        self.processes = []
        self.processes_count = 0
        self.threads = []
        self.threads_count = 0
        self.results = {}

        self.pex = None
        if self.workers:
            if mpi:
                from mpi4py.futures import MPIPoolExecutor
                self.pex = MPIPoolExecutor(max_workers=self.workers)
            else:
                from concurrent.futures import ProcessPoolExecutor
                self.pex = ProcessPoolExecutor(max_workers=self.workers)

    def __enter__(self, max_workers=0):
        return self

    def __exit__(self, *exc):
        self.shutdown()
        return False

    def __str__(self):
        return  "<EnhancedExecutor spawned: %d, processes: %d, threads: %d>" % (
            self.workers, len(self.processes), len(self.threads))

    def submit(self, func, *args):
        # Check if threading is to be used
        if not self.use_threads:
            thread = False
        elif not self.workers or not self.threads:  # make the first job is assigned to main
            thread = True
        else:
            quotient = self.processes_count // self.workers
            remainder = self.processes_count % self.workers
            if remainder == 0 and quotient > self.threads_count:
                thread = True
            else:
                thread = False

        # Launch thread or process
        if thread:
            future = EnhancedFuture(self, 'thread')
            future.submit(func, *args)
            self.threads.append(future)
            self.threads_count += 1
        else:
            future = EnhancedFuture(self, 'process')
            future.submit(func, *args)
            self.processes.append(future)
            self.processes_count += 1

        return future

    def map(self, func, *args):
        assert args[0] and len(set([len(x) for x in args])) == 1
        futures = [self.submit(func, *a) for a in zip(*args)]
        print(self)
        self.run()
        return [self.result(f) for f in futures]

    def result(self, future):
        return future.result()

    def run(self):
        for future in self.threads:
            _ = future.result()

    def shutdown(self):
        self.run()
        if not self.pex is None:
            self.pex.shutdown()


class EnhancedMPIPoolExecutor(EnhancedProcessPoolExecutor):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs, mpi=True)


def as_completed(futures):
    ptypes = [type(f.executor) for f in futures]
    assert len(set(ptypes)) == 1
    ptype = ptypes[0]

    threads = [f for f in futures if f.kind == 'thread']
    processes = [f for f in futures if f.kind == 'process']

    if processes:
        proc_dct = {f.process: f for f in futures}
        if ptype == EnhancedMPIPoolExecutor:
            from mpi4py.futures import as_completed
        else:
            from concurrent.futures import as_completed
        for future in as_completed(proc_dct):
            yield proc_dct[future]

    if threads:
        for future in threads:
            yield future
