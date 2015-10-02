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

'''

Main script that initiates the BreaKmer analysis or auxiliary functions to setup BreaKmer for analysis.

There are three functions provided:
1. run                    = perform analysis to detect structural variation.
2. start_blat_server      = start the blat server in the background for analysis.
3. prepare_reference_data = prepare the reference data for the target regions that are specified in the input files.


The blat server provides a challenge in workflow. The best method is to:
1. prepare reference data using 'prepare_reference_data' function
2. start the blat server using 'start_blat_server' function
    breakmer.py start_blat_server -p <port_number> --hostname <hostname> -c <config file>
3. run the analysis and keep the blat server alive in the background for use in other analyses.
    breakmer.py run -k -p <port_number> --hostname <hostname> -c <config file> -n <nprocessors> -g <gene_list>

'''

args = sys.argv

PARSER = argparse.ArgumentParser(description='Program to identify structural variants within targeted locations.', usage='%(prog)s [options]', add_help=False)
SUBPARSERS = PARSER.add_subparsers(help='Program mode (run, start_blat_server, prepare_reference_data).', dest='fncCmd')

RUN_PARSER = SUBPARSERS.add_parser('run', help='Run analysis to detect structural variants.')
SERVER_PARSER = SUBPARSERS.add_parser('start_blat_server', help='Start the blat server prior to performing the analysis.')
REF_PARSER = SUBPARSERS.add_parser('prepare_reference_data', help='Prepare the reference sequence data for target regions prior to analysis.')

RUN_PARSER.add_argument('--log_level', dest='log_level', default='DEBUG', help='Log level [default: DEBUG]')
RUN_PARSER.add_argument('--indel_size', dest='indel_size', default=15, type=int, help='Indel size filter. [default: %(default)s]')
RUN_PARSER.add_argument('--trl_sr_thresh', dest='trl_sr_thresh', default=2, type=int, help='Split read support threshold for translocations. [default: %(default)s]')
RUN_PARSER.add_argument('--indel_sr_thresh', dest='indel_sr_thresh', default=5, type=int, help='Split read support threshold for indels. [default: %(default)s]')
RUN_PARSER.add_argument('--rearr_sr_thresh', dest='rearr_sr_thresh', default=2, type=int, help='Split read support threshold for rearrangements. [default: %(default)s]')
RUN_PARSER.add_argument('--rearr_min_seg_len', dest='rearr_minseg_len', default=30, type=int, help='Threshold for minimum segment to be rearranged. [default: %(default)s]')
RUN_PARSER.add_argument('--trl_min_seg_len', dest='trl_minseg_len', default=25, type=int, help='Threshold for minimum length of translocation segment. [default: %(default)s]')
RUN_PARSER.add_argument('--align_thresh', dest='align_thresh', default=.90, type=int, help='Threshold for minimum read alignment for assembly. [default: %(default)s]')
RUN_PARSER.add_argument('--no_output_header', dest='no_output_header', default=False, action='store_true', help='Suppress output headers. [default: %(default)s]')
RUN_PARSER.add_argument('--discread_only_thresh', dest='discread_only_thresh', default=2, type=int, help='The number of discordant read pairs in a cluster to output without evidence from a split read event. [default: %(default)s]')
RUN_PARSER.add_argument('--generate_image', dest='generate_image', default=False, action='store_true', help='Generate pileup image for events. [default: %(default)s]')
RUN_PARSER.add_argument('--hostname', dest='blat_hostname', default='localhost', help='The hostname for the blat server. Localhost will be used if not specified. [default: %(default)s]')
RUN_PARSER.add_argument('-g', '--gene_list', dest='gene_list', default=None, help='Gene list to consider for analysis. [default: %(default)s]')
RUN_PARSER.add_argument('-f', '--filter_list', dest='filterList', default=None, help='Input a set of events to filter out. [default: %(default)s]')
RUN_PARSER.add_argument('-n', '--nprocessors', dest='nprocs', default=1, type=int, help='The number of processors to use for analysis. [default: %(default)s]')
RUN_PARSER.add_argument('-s', '--start_blat_server', dest='start_blat_server', default=False, action='store_true', help='Start the blat server. Random port number and localhost will be used if neither specified. [default: %(default)s]')
RUN_PARSER.add_argument('-k', '--keep_blat_server', dest='keep_blat_server', default=False, action='store_true', help='Keep the blat server alive. [default: %(default)s]')
RUN_PARSER.add_argument('-p', '--port_number', dest='blat_port', default=None, type=int, help='The port number for the blat server. A random port number (8000-9500) will be used if not specified. [default: %(default)s]')
RUN_PARSER.add_argument('-c', '--config', dest='config_fn', default=None, required=True, help='The configuration filename that contains additional parameters. [default: %(default)s]')

SERVER_PARSER.add_argument('-p', '--port_number', dest='blat_port', default=None, type=int, help='The port number for the blat server. A random port number (8000-9500) will be used if not specified. [default: %(default)s]')
SERVER_PARSER.add_argument('--hostname', dest='blat_hostname', default='localhost', help='The hostname for the blat server. Localhost will be used if not specified. [default: %(default)s]')
SERVER_PARSER.add_argument('-c', '--config', dest='config_fn', default=None, required=True, help='The configuration filename that contains additional parameters. [default: %(default)s]')

REF_PARSER.add_argument('-g', '--gene_list', dest='gene_list', default=None, help='Gene list to consider for analysis. [default: %(default)s]')
REF_PARSER.add_argument('-c', '--config', dest='config_fn', default=None, required=True, help='The configuration filename that contains additional parameters. [default: %(default)s]')
REF_PARSER.add_argument('-n', '--nprocessors', dest='nprocs', default=1, type=int, help='The number of processors to use for analysis. [default: %(default)s]')

RUN_TRACKER = breakmer_analysis.RunTracker(params.ParamManager(PARSER.parse_args()))
RUN_TRACKER.run()
