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
    def __init__(self, max_workers=0):
        from concurrent.futures import ProcessPoolExecutor#, ThreadPoolExecutor
        self.workers = max_workers
        self.pex = ProcessPoolExecutor(max_workers=self.workers) if self.workers else None
        #self.tex = ThreadPoolExecutor(max_workers=1) # MIXING threads and processes DOESN'T WORK
        self.processes = []
        self.processes_count = 0
        self.threads = []
        self.threads_count = 0
        self.results = {}

    def __enter__(self, max_workers=0):
        return self

    def __exit__(self, *exc):
        self.shutdown()
        return False

    def __str__(self):
        return  "<EnhancedExecutor spawned: %d, processes: %d, threads: %d>" % (
            self.workers, len(self.processes), len(self.threads))

    def submit(self, func, *args):
        if not self.workers or not self.threads:
            thread = True
        else:
            quotient = self.processes_count // self.workers
            remainder = self.processes_count % self.workers
            if remainder == 0 and quotient > self.threads_count:
                thread = True
            else:
                thread = False

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

    def result(self, future):
        return future.result()

    def run(self):
        for future in self.threads:
            _ = future.result()

    def shutdown(self):
        if self.pex is not None:
            self.pex.shutdown()


class EnhancedMPIPoolExecutor(EnhancedProcessPoolExecutor):
    def __init__(self, max_workers=0):
        super().__init__(max_workers)
        from mpi4py.futures import MPIPoolExecutor
        self.pex = MPIPoolExecutor(max_workers=self.workers) if self.workers else None
