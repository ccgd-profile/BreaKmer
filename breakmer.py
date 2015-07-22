#! /usr/bin/local/python
# -*- coding: utf-8 -*-

import sys
import argparse
import breakmer.params as params
import breakmer.processor.analysis as breakmer_analysis

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"

args = sys.argv
if len(args) < 2:
    print 'No arguments provided. At least one command required (run, start_blat_server, prepare_reference_kmers).'
    sys.exit(1)
else:
    PARSER = argparse.ArgumentParser(description="Program to identify structural variants within targeted locations.", usage='%(prog)s [options]', add_help=True)
    fncCmd = args[1]
    if fncCmd == 'run':
        # PARSER.add_argument('funcCommand', help='Command for BreaKmer function to execute. Required')
        # PARSER.add_argument('configFn', help='Path to the configuration file. Required')
        # PARSER.add_argument('-a', '--keep_repeat_regions', dest='keep_repeat_regions', default=False, action='store_true', help='Keep indels in repeat regions. No repeat mask bed file required if set. [default: False]')
        # PARSER.add_argument('-p', '--preset_ref_data', dest='preset_ref_data', default=False, action='store_true', help='Preset all the reference data for all the targets before running analysis. [default: False]')
        # PARSER.add_argument('-k', '--keep_intron_vars', dest='keep_intron_vars', default=False, action='store_true', help='Keep intronic indels or rearrangements [default: %(default)s]')
        # PARSER.add_argument('-v', '--var_filter', dest='var_filter', default='all', help='Variant types to report (all, indel, trl, rearrangment) [default: %(default)s]')
        PARSER.add_argument('--log_level', dest='log_level', default='DEBUG', help='Log level [default: DEBUG]')
        PARSER.add_argument('--indel_size', dest='indel_size', default=15, type=int, help='Indel size filter. [default: %(default)s]')
        PARSER.add_argument('--trl_sr_thresh', dest='trl_sr_thresh', default=2, type=int, help='Split read support threshold for translocations. [default: %(default)s]')
        PARSER.add_argument('--indel_sr_thresh', dest='indel_sr_thresh', default=5, type=int, help='Split read support threshold for indels. [default: %(default)s]')
        PARSER.add_argument('--rearr_sr_thresh', dest='rearr_sr_thresh', default=2, type=int, help='Split read support threshold for rearrangements. [default: %(default)s]')
        PARSER.add_argument('--rearr_min_seg_len', dest='rearr_minseg_len', default=30, type=int, help='Threshold for minimum segment to be rearranged. [default: %(default)s]')
        PARSER.add_argument('--trl_min_seg_len', dest='trl_minseg_len', default=25, type=int, help='Threshold for minimum length of translocation segment. [default: %(default)s]')
        PARSER.add_argument('--align_thresh', dest='align_thresh', default=.90, type=int, help='Threshold for minimum read alignment for assembly. [default: %(default)s]')
        PARSER.add_argument('--no_output_header', dest='no_output_header', default=False, action='store_true', help='Suppress output headers. [default: %(default)s]')
        PARSER.add_argument('--discread_only_thresh', dest='discread_only_thresh', default=2, type=int, help='The number of discordant read pairs in a cluster to output without evidence from a split read event. [default: %(default)s]')
        PARSER.add_argument('--generate_image', dest='generate_image', default=False, action='store_true', help='Generate pileup image for events. [default: %(default)s]')
        PARSER.add_argument('-g', '--gene_list', dest='gene_list', default=None, help='Gene list to consider for analysis. [default: %(default)s]')
        PARSER.add_argument('-f', '--filter_list', dest='filterList', default=None, help='Input a set of events to filter out. [default: %(default)s]')
        PARSER.add_argument('-n', '--nprocessors', dest='nprocs', default=1, type=int, help='The number of processors to use for analysis. [default: %(default)s]')
        PARSER.add_argument('-s', '--start_blat_server', dest='start_blat_server', default=False, action='store_true', help='Start the blat server. Random port number and localhost will be used if neither specified. [default: %(default)s]')
        PARSER.add_argument('-k', '--keep_blat_server', dest='keep_blat_server', default=False, action='store_true', help='Keep the blat server alive. [default: %(default)s]')
        PARSER.add_argument('-p', '--port_number', dest='blat_port', default=None, type=int, help='The port number for the blat server. A random port number (8000-9500) will be used if not specified. [default: %(default)s]')
        PARSER.add_argument('--hostname', dest='blat_hostname', default='localhost', help='The hostname for the blat server. Localhost will be used if not specified. [default: %(default)s]')
        PARSER.add_argument('-c', '--config', dest='config_fn', default=None, required=True, help='The configuration filename that contains additional parameters. [default: %(default)s]')

    elif fncCmd == 'start_blat_server':
        PARSER.add_argument('-p', '--port_number', dest='blat_port', default=None, type=int, help='The port number for the blat server. A random port number (8000-9500) will be used if not specified. [default: %(default)s]')
        PARSER.add_argument('--hostname', dest='blat_hostname', default='localhost', help='The hostname for the blat server. Localhost will be used if not specified. [default: %(default)s]')
        PARSER.add_argument('-c', '--config', dest='config_fn', default=None, required=True, help='The configuration filename that contains additional parameters. [default: %(default)s]')

    elif fncCmd == 'prepare_reference_data':
        PARSER.add_argument('-g', '--gene_list', dest='gene_list', default=None, help='Gene list to consider for analysis. [default: %(default)s]')
        PARSER.add_argument('-c', '--config', dest='config_fn', default=None, required=True, help='The configuration filename that contains additional parameters. [default: %(default)s]')
        PARSER.add_argument('-n', '--nprocessors', dest='nprocs', default=1, type=int, help='The number of processors to use for analysis. [default: %(default)s]')
    else:
        print 'The first argument %s is not one of the options.' % fncCmd
        sys.exit(1)

    RUN_TRACKER = breakmer_analysis.RunTracker(params.ParamManager(fncCmd, PARSER.parse_args(args[2:])))
    RUN_TRACKER.run()
