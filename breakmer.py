#! /usr/bin/local/python

# Fix get_svs_result()
# Fix filter_by_feature

import sys
import os
from utils import *
from sv_processor import *
import logging
from optparse import OptionParser

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

#-----------------------------------------------------------
# TODO: Set required options, exit if these are not set.
#-----------------------------------------------------------
def parse_config_file_and_command_line_options(config_file_name, command_line_options) :
  param_opts = {}
  config_file = open(config_file_name, 'rU')
  flines = config_file.readlines()
  for line in flines :
    line = line.strip()

    # Allow blank lines
    if (len (line) == 0):
      continue

    # Allow comment lines
    if (line.startswith('#')):
      continue

    linesplit = line.split("=")
    if len(linesplit) == 1 : 
      print 'Config line', line, ' not set correctly. Exiting.'
      sys.exit()
    else :
      key,value = linesplit
      param_opts[key] = value

  # This pulls the key-value pairs out of the command line options
  # and adds them to the dictionary created from the config file
  for opt in vars(command_line_options) :
    param_opts[opt] = vars(command_line_options)[opt]

  return param_opts
#-----------------------------------------------------------

#````````````````````````````````````````````````````````````
usage = '%prog [options] <config file name>'
desc = """Script to identify structural variants within targeted locations."""
parser = OptionParser(usage=usage,description=desc)
parser.add_option('-l', '--log_level', dest='log_level', default='DEBUG', type='string', help='Log level [default: %default]')
parser.add_option('-a', '--keep_repeat_regions', dest='keep_repeat_regions', default=False, action='store_true', help='Keep indels in repeat regions. Requires a repeat mask bed file. [default: %default]')
parser.add_option('-p', '--preset_ref_data', dest='preset_ref_data', default=False, action="store_true", help='Preset all the reference data for all the targets before running analysis. [default: %default]')
parser.add_option('-s', '--indel_size', dest='indel_size', default=15, type='int', help='Indel size filter [default: %default]')
parser.add_option('-c', '--trl_sr_thresh', dest='trl_sr_thresh', default=2, type='int', help='Split read support threshold for translocations [default: %default]')
parser.add_option('-d', '--indel_sr_thresh', dest='indel_sr_thresh', default=5, type='int', help='Split read support threshold for indels [default: %default]')
parser.add_option('-r', '--rearr_sr_thresh', dest='rearr_sr_thresh', default=3, type='int', help='Split read support threshold for rearrangements [default: %default]')
parser.add_option('-g', '--gene_list', dest='gene_list', default=None, type='string', help='Gene list to consider for analysis [default: %default]')
parser.add_option('-k', '--keep_intron_vars', dest='keep_intron_vars', default=False, action='store_true', help='Keep intronic indels or rearrangements [default: %default]')
parser.add_option('-v', '--var_filter', dest='var_filter', default='all', type='string', help='Variant types to report (all, indel, trl, rearrangment) [default: %default]')
parser.add_option('-m', '--rearr_min_seg_len', dest='rearr_minseg_len', default=30, type='int', help='Threshold for minimum segment to be rearranged [default: %default]')
parser.add_option('-n', '--trl_min_seg_len', dest='trl_minseg_len', default=25, type='int', help='Threshold for minimum length of translocation segment [default: %default]')
parser.add_option('-t', '--align_thresh', dest='align_thresh', default=.90, type='int', help='Threshold for minimum read alignment for assembly [default: %default]')
parser.add_option('-z', '--no_output_header', dest='no_output_header', default=False, action='store_true', help='Suppress output headers [default: %default]')


parser.add_option('-o', '--tar_gzip_output', dest='tar_gzip_output', default=False, action='store_true', help='Create a tar.gz file of the final output.  Useful for running on the cloud. [default: %default]')


if __name__ == '__main__' :
  command_line_options, args = parser.parse_args(sys.argv[1:])
  config_file_name = args[0]

  start_time = time.clock()
  config_dictionary = parse_config_file_and_command_line_options(config_file_name,command_line_options)
  setup_logger(config_dictionary,'root')
  r = runner(config_dictionary)
  r.run(start_time)
#````````````````````````````````````````````````````````````
