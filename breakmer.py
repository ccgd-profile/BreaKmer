#! /usr/bin/local/python

import sys
import utils
import time
from sv_processor import runner
import argparse


'''
Author: Ryan Abo
12/23/2013

Requires a configuration file containing key=value style formatting.

This is the main script that initiates the analysis. 

Submodules:
utils.py - utility functions.
 classes:
 - fq_read
 - FastqFile
 - anno
 - params
sv_processor.py - classes to manage analysis, process, analyze, write.
 classes:
 - runner - analysis-wide information and analysis (loop through each)
 - target - target-level information and analysis
 - contig - contig-level information and analysis
sv_assembly.py - classes to take sample kmers and produce contigs.
 classes:
 - kmers
 - kmer 
 - buffer
 - contig
sv_caller.py - call SVs from blat result.
 classes:
 - align_manager
 - blat_manager
 - blat_res
 - sv_event
'''


def parse_config_f(config_fn, opts):
    param_opts = {}
    for line in open(config_fn, 'rU'):
        line = line.strip()
        linesplit = line.split("=")
        if len(linesplit) == 1 : 
            print 'Configuration line', line, 'not set correctly. Exiting.'
            sys.exit(2)
        else :
            k,v = linesplit
            param_opts[k] = v
    for opt in vars(opts): 
        param_opts[opt] = vars(opts)[opt]
    return param_opts


parser = argparse.ArgumentParser(description="Program to identify structural variants within targeted locations.", usage='%(prog)s [options]', add_help=True)
parser.add_argument('config_fn', help='Path to the configuration file. Required')
parser.add_argument('-l', '--log_level', dest='log_level', default='DEBUG', help='Log level [default: DEBUG]')
parser.add_argument('-a', '--keep_repeat_regions', dest='keep_repeat_regions', default=False, action='store_true', help='Keep indels in repeat regions. No repeat mask bed file required if set. [default: False]')
parser.add_argument('-p', '--preset_ref_data', dest='preset_ref_data', default=False, action='store_true', help='Preset all the reference data for all the targets before running analysis. [default: False]')
parser.add_argument('-s', '--indel_size', dest='indel_size', default=15, type=int, help='Indel size filter [default: %(default)s]')
parser.add_argument('-c', '--trl_sr_thresh', dest='trl_sr_thresh', default=2, type=int, help='Split read support threshold for translocations [default: %(default)s]')
parser.add_argument('-d', '--indel_sr_thresh', dest='indel_sr_thresh', default=5, type=int, help='Split read support threshold for indels [default: %(default)s]')
parser.add_argument('-r', '--rearr_sr_thresh', dest='rearr_sr_thresh', default=3, type=int, help='Split read support threshold for rearrangements [default: %(default)s]')
parser.add_argument('-g', '--gene_list', dest='gene_list', default=None, help='Gene list to consider for analysis [default: %(default)s]')
parser.add_argument('-k', '--keep_intron_vars', dest='keep_intron_vars', default=False, action='store_true', help='Keep intronic indels or rearrangements [default: %(default)s]')
parser.add_argument('-v', '--var_filter', dest='var_filter', default='all', help='Variant types to report (all, indel, trl, rearrangment) [default: %(default)s]')
parser.add_argument('-m', '--rearr_min_seg_len', dest='rearr_minseg_len', default=30, type=int, help='Threshold for minimum segment to be rearranged [default: %(default)s]')
parser.add_argument('-n', '--trl_min_seg_len', dest='trl_minseg_len', default=25, type=int, help='Threshold for minimum length of translocation segment [default: %(default)s]')
parser.add_argument('-t', '--align_thresh', dest='align_thresh', default=.90, type=int, help='Threshold for minimum read alignment for assembly [default: %(default)s]')
parser.add_argument('-z', '--no_output_header', dest='no_output_header', default=False, action='store_true', help='Suppress output headers [default: %(default)s]')
parser.add_argument('-f', '--filter_list', dest='filter_list', default=None, help='Input a set of events to filter out. [default: %(default)s]')


if __name__ == '__main__':
    try:
        args = sys.argv[1:]
        opts = parser.parse_args()
        if len(args) != 1:
            parser.error('Requires input of configuration file name')
            parser.usage()
        else:
            config_fn = args[0]
    except parser.error:
        parser.usage()
        sys.exit(2)

    tic = time.clock()
    config_d = parse_config_f(config_fn, opts)
    utils.setup_logger(config_d, 'root')
    r = runner(config_d)
    r.run(tic)