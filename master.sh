#!/bin/bash

# preliminar.sh                       contains qualcheck etc using bwameth
# bowtie_based.sh                     same with bowtie2
# stranded_dna_meth_caller.sh         meth calls using methyldackel using really convoluted sam flags parsing
# simple_stranded_dna_meth_caller.sh  kiss version of the upper
# simple_stranded_dna_meth_caller.R   analyzes the upper
# paired_end_workflow                 detected that a sample did not map as pe, writing the whole flow


# this is getting messy, let's clarify
# bowtie_based.sh                         bismark stuff
# paired_end_workflow_bulk.sh             bwa meth paired end stuff
# deptools_comparison.sh                  test on deptools representations just in case there were
                                             # bwa and bowtie mapping differences
# diagnostic_plots.R                      preliminar integration
