#!/usr/bin/env bash

wget https://pubs.acs.org/doi/suppl/10.1021/acs.jproteome.6b00971/suppl_file/pr6b00971_si_001.xlsx
wget https://pubs.acs.org/doi/suppl/10.1021/acs.jproteome.6b00971/suppl_file/pr6b00971_si_002.xlsx

xlsx2csv -d tab MS/pr6b00971_si_001.xlsx | sed '/^\t*$/d' > MS/pr6b00971_si_001_EL4.txt

