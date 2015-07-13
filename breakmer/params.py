#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import logging
import random
import subprocess
import time
import pysam
import breakmer.utils as utils
import breakmer.caller.filter as resultfilter

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


class ParamManager:
    """ParamManager class stores all the input specifications provided to the program to run. These include
    file paths, thresholds, directories, etc...
    Attributes:
        opts (dict):                   Containing parameter options and input values as key-values.
        gene_annotations (Annotation): Tracks the annotation information.
        targets (dict):                Target region coordinates, key-values.
        paths (dict):                  Dictionary containing the top level directories for the analysis output.
        logging_name (str):            Logging string name object for logging messages.
        repeat_mask: Dictionary containing chromosome IDs as keys with a list of tuples for each repeat
                     feature (i.e., chr, left_bp, right_bp, repeat_name).
    """

    def __init__(self, fncCmd, arguments):
        """Initializes params class.
        Args:
            arguments: The argparse dictionary object from the command line parameters.
        Returns:
            None
        Raises:
            None
        """
        self.opts = None
        self.gene_annotations = None
        self.filter = None
        self.targets = {}
        self.paths = {}
        self.repeat_mask = None
        self.loggingName = 'breakmer.params'
        self.fncCmd = fncCmd
        self.set_params(arguments)

    def set_params(self, arguments):
        """Organize and format all input parameters into class variables to access
        later. Specific instances of parameters are checked and set. All parameters that are
        set are logged. The target objects are set along with the paths.
        Args:
            arguments: The argparse dictionary object from the command line options.
        Returns:
            No returns, it sets the class variable self.opts.
        """
        # Parse config file and command line options store in self.opts
        self.parse_opts(arguments)
        # Create logging object.
        utils.setup_logger(self.get_param('analysis_dir', True), 'breakmer')
        utils.log(self.loggingName, 'info', 'Setting up parameters')

        """DEPRECATED
        # Set the output filter
        var_filter = self.opts['var_filter']
        if var_filter == 'all':
            self.opts['var_filter'] = ['indel', 'rearrangement', 'trl']
        else:
            self.opts['var_filter'] = var_filter.split(",")
            if 'indel' in self.opts['var_filter'] or 'rearrangement' in self.opts['var_filter'] or 'trl' in self.opts['var_filter']:
                utils.log(self.loggingName, 'info', 'Variant filters %s set, only reporting these variants' % ','.join(self.opts['var_filter']))
            else:
                utils.log(self.loggingName, 'debug', 'Variant filter options %s are not valid, using all' % ','.join(self.opts['var_filter']))
                self.opts['var_filter'] = ['indel', 'rearrangement', 'trl']
        """

        # Log all parameters passed in, warn for poor paths
        for param in self.opts:
            value = self.opts[param]
            utils.log(self.loggingName, 'info', '%s = %s' % (param, value))

        self.set_targets(self.get_param('gene_list'))
        self.paths['ref_data'] = os.path.abspath(os.path.normpath(self.opts['reference_data_dir']))
        self.opts['reference_fasta_dir'] = os.path.split(self.opts['reference_fasta'])[0]

        """DEPRECATED
        # Setup alternative reference sequence files if they were passed in.
        if 'alternate_reference_fastas' in self.opts:
            utils.log(self.loggingName, 'info', 'Alternate reference fastas listed in configuration %s' % self.opts['alternate_reference_fastas'])
            self.opts['alternate_reference_fastas'] = self.opts['alternate_reference_fastas'].split(',')
            # Check for altref 2bit files
            for fn in self.opts['alternate_reference_fastas']:
                twobit_fn = os.path.splitext(fn)[0] + '.2bit'
                if not os.path.exists(twobit_fn):
                    utils.log(self.loggingName, 'info', 'Creating 2bit from %s alternate reference fasta' % fn)
                    curdir = os.getcwd()
                    os.chdir(os.path.dirname(fn))
                    cmd = '%s %s %s' % (self.opts['fatotwobit'], fn, twobit_fn)
                    utils.log(self.loggingName, 'info', 'fatotwobit command %s' % cmd)
                    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                    output, errors = p.communicate()
                    os.chdir(curdir)
        else:
            utils.log(self.loggingName, 'info', 'No alternate reference fasta files listed in configuration.')
        """

        # If only preseting the reference data, then move on.
        if self.fncCmd == 'prepare_reference_data':
            utils.log(self.loggingName, 'info', 'Preset reference data option set! Only the reference data directory will be setup.')
            return

        # TODO - set this as an optional parameter.
        # self.gene_annotations = utils.Annotation()
        # self.gene_annotations.add_genes(self.opts['gene_annotation_file'])

        # if 'other_regions_file' in self.opts:
        #     self.gene_annotations.add_regions(self.opts['other_regions_file'])

        # Setup analysis directories
        self.paths['analysis'] = os.path.abspath(os.path.normpath(self.opts['analysis_dir']))
        self.paths['output'] = os.path.join(self.paths['analysis'], 'output')
        if 'targets_dir' in self.opts:
            self.paths['targets'] = os.path.abspath(os.path.normpath(self.opts['targets_dir']))
        else:
            self.paths['targets'] = os.path.join(self.paths['analysis'], 'targets')

        for path in self.paths:
            utils.log(self.loggingName, 'info', 'Creating %s directory (%s)' % (path, self.paths[path]))
            if not os.path.exists(self.paths[path]):
                os.makedirs(self.paths[path])

        if self.fncCmd == 'start_blat_server':
            utils.log(self.loggingName, 'info', 'Starting the blat server.')
            return

        # Check if Jellyfish and Cutadapt work.
        self.check_binaries()
        self.filter = resultfilter.ResultFilter(self.get_param('filterList'), self)

        # Set the expected insert size threshold from the properly mapped read pairs.
        self.set_insertsize_thresh()

        # Set repeats if specified
        # if not self.opts['keep_repeat_regions'] and 'repeat_mask_file' in self.opts:
        #     utils.log(self.loggingName, 'info', 'Storing all repeats by chrom from file %s' % self.opts['repeat_mask_file'])
        #     self.repeat_mask = utils.setup_rmask_all(self.opts['repeat_mask_file'])

    def set_insertsize_thresh(self):
        """Store the insert sizes for a small number of "properly mapped" reads
        and determine an upperbound cutoff to use to determine discordantly mapped read
        pairs.
        """
        nSampleReads = 100000
        bamF = pysam.Samfile(self.get_param('sample_bam_file'), 'rb')
        testReads = bamF.fetch()
        insertSizes = []
        readIter = 0
        for read in testReads:
            proper_map = read.flag == 83 or read.flag == 99 and read.mapq > 0
            if read.is_read1 and proper_map:
                readIter += 1
                insertSizes.append(abs(read.tlen))
            if readIter == nSampleReads:
                break
        isMedian = utils.median(insertSizes)
        isSD = utils.stddev(utils.remove_outliers(insertSizes))
        self.opts['insertsize_thresh'] = isMedian + (5 * isSD)

    def parse_opts(self, arguments):
        """Formats input parameters into dictionary, self.opts
        It first parses the configuration file and stores the key, values in the self.opts dictionary.
        It will exit with error if the configuration file does not have lines in the proper format (i.e.,
        key=value). It will also iterate through the command line paramaters and store the keys and values
        in the opts dictionary. A final check if performed for the required parameters depending on the
        parameters that have been passed.
        Args:
            arguments: The argparse dictionary object from the command line options.
        Returns:
            No returns, it sets the class variable self.opts.
        """

        self.opts = {}
        for line in open(arguments.config_fn, 'rU'):
            line = line.strip()
            # Allow for blank lines and comments
            if line == '' or line.find('#') > -1:
                continue
            linesplit = line.split("=")
            if len(linesplit) == 1:
                err_msg = 'Config line', line, ' not set correctly. Exiting.'
                print err_msg
                utils.log(self.loggingName, 'error', err_msg)
                sys.exit(1)
            else:
                key, value = linesplit
                self.opts[key] = value

        for opt in vars(arguments):
            self.opts[opt] = vars(arguments)[opt]

        # Sanity check for required params
        # Required when preset_ref_data = True:
        # - reference_data_dir
        # - reference_fasta
        # - targets_bed_file

        # Required when preset_ref_data = False
        # - analysis_name
        # - targets_bed_file
        # - sample_bam_file
        # - analysis_dir
        # - reference data_dir
        # - cutadapt_config_file
        # - reference_fasta
        # - gene_annotation_file

        required = ['analysis_name', 'targets_bed_file', 'sample_bam_file', 'analysis_dir', 'reference_data_dir', 'cutadapt_config_file', 'reference_fasta', 'gene_annotation_file']
        if self.get_param('preset_ref_data'):
            required = ['reference_data_dir', 'reference_fasta', 'targets_bed_file']

        for req in required:
            self.get_param(req, True)

    def check_binaries(self):
        """Checks if the six required binaries are available to use.
        There are six required binaries to perform the complete analysis (blat, gfserver,
        gfclient, fatotwobit, cutadapt, jellyfish). Each binary is checked whether
        the path provided in the configuration file has an executable file attached or
        if no path was provided that the binary is on the path. Cutadapt and Jellyfish
        are also tested using hardcoded data.
        Args:
            None
        Returns:
            None
        """

        binaries = ('blat', 'gfserver', 'gfclient', 'fatotwobit', 'cutadapt', 'jellyfish')
        for bin in binaries:
            # Check provided binary
            if bin in self.opts:
                bin_path = self.opts[bin]
                bin_check = utils.which(bin_path)
            # Check if binary in path
            else:
                bin_check = utils.which(bin)
                self.opts[bin] = bin_check
            if not bin_check:
                print 'Missing path/executable for', bin
                utils.log(self.loggingName, 'error', 'Missing path/executable for %s' % bin)
                sys.exit(1)
            utils.log(self.loggingName, 'info', '%s path = %s' % (bin, bin_check))
        utils.log(self.loggingName, 'info', 'All the required binaries have been check successfully!')

        # Test cutadapt and jellyfish binaries
        test_dir = os.path.join(self.paths['analysis'], 'bin_test')
        test_fq = os.path.join(test_dir, 'test.fq')
        if not os.path.exists(test_dir):
            os.makedirs(test_dir)
        utils.write_test_fq(test_fq)
        clean_fq, rc = utils.test_cutadapt(test_fq, self.opts['cutadapt'], self.opts['cutadapt_config_file'])
        if clean_fq:
            utils.log(self.loggingName, 'info', 'Test cutadapt ran successfully')
            jfish_prgm, rc = utils.test_jellyfish(self.opts['jellyfish'], clean_fq, test_dir)
            if rc != 0:
                utils.log(self.loggingName, 'error', '%s unable to run successfully, exit code %s. Check installation and correct version.' % (jfish_prgm, str(rc)))
                sys.exit(1)
            else:
                utils.log(self.loggingName, 'info', 'Test jellyfish ran successfully')
        else:
            utils.log(self.loggingName, 'error', 'Cutadapt failed to run, exit code %s. Check installation and version.' % str(rc))
            sys.exit(1)

    def set_targets(self, geneList):
        """Parse the targets bed file and store them in a dictionary. Limit to a gene
        list if input.
        A list of genes can be passed in by the user to limit the analysis. This will
        limit which targets are stored in the dictionary as the target bed file is parsed.
        The target bed file is a tab-delimited text file that should have at minimum,
        four columns (chromosome, start, end, name) with an optional fourth column
        containing a coding feature (i.e., exon or intron). Each row is either a tiled
        region with sequencing coverage or it is just a region to analyze by BreaKmer.
        The name can be applied to multiple rows, and if multiple tiled regions are input
        with the same name they are aggregated together under the same key.
        Args:
            gene_list: A list containing gene names that match those in the target bed
                       file.
        Returns:
            None
        """

        regionList = None
        if geneList:
            regionList = []
            for line in open(geneList, 'r'):
                line = line.strip()
                regionList.append(line.upper())

        utils.log(self.loggingName, 'info', 'Parsing target list')
        # TODO: Check to make sure there aren't duplicate genes.
        cur_region = ['', []]
        for target in open(self.opts['targets_bed_file'], 'rU'):
            # Each target is formatted like a bed, chr bp1 bp2 name
            target = target.strip()
            targetsplit = target.split()
            chrm, bp1, bp2, name = targetsplit[0:4]
            if regionList:
                if name.upper() not in regionList:
                    continue
            feature = None
            # Allow a fifth column containing indication of what type of region it is.
            # Typically exon/intron designation. This will be deprecated.
            if len(targetsplit) > 4:
                feature = targetsplit[4]
            # All the names are converted to uppercase.
            if name.upper() not in self.targets:
                self.targets[name.upper()] = []
            self.targets[name.upper()].append((chrm, int(bp1), int(bp2), name, feature))
        utils.log(self.loggingName, 'info', '%d targets' % len(self.targets))

    def check_blat_server(self):
        """

        """
        test_dir = os.path.join(self.paths['analysis'], 'blatserver_test')
        test_fa_fn = os.path.join(test_dir, 'test.fa')
        if not os.path.exists(test_dir):
            os.makedirs(test_dir)
        test_fa = open(test_fa_fn, 'w')
        test_fa.write('>test\nCCAAGGGAGACTTCAAGCAGAAAATCTTTAAGGGACCCTTGCATAGCCAGAAGTCCTTTTCAGGCTGATGTACATAAAATATTTAGTAGCCAGGACAGTAGAAGGACTGAAGAGTGAGAGGAGCTCCCAGGGCCTGGAAAGGCCACTTTGTAAGCTCATTCTTG')
        test_fa.close()

        resultFn = os.path.join(test_dir, 'blatserver_test.psl')
        cmd = '%s -t=dna -q=dna -out=psl -minScore=20 -nohead %s %d %s %s %s' % (self.get_param('gfclient'), self.get_param('blat_hostname'), self.get_param('blat_port'), self.get_param('reference_fasta_dir'), test_fa_fn, resultFn)
        utils.log(self.loggingName, 'info', 'Blat server test system command %s' % cmd)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, errors = p.communicate()
        utils.log(self.loggingName, 'info', 'Realignment output file %s' % resultFn)
        serverSuccess = True
        print 'Blat errors', errors
        if errors != '':
            serverSuccess = False
            utils.log(self.loggingName, 'info', 'Realignment errors %s' % errors)
        print 'Sever success', serverSuccess
        return serverSuccess

    def start_blat_server(self):
        """Fire up a blat server instance using a random port number and localhost.
        The required files to start a blat server are first checked and created, if
        necessary. These include a genome-wide reference fasta file and a 2bit
        file generated from that fasta file. The faToTwoBit program is used if the
        2bit file needs to be generated on the fly. The gfServer is started and
        we wait while the server is successfully started.
        Args:
            None
        Return:
            None
        """

        if self.fncCmd == 'prepare_reference_data':
            return
        elif self.fncCmd == 'start_blat_server':
            port = self.get_param('blat_port')
            hostname = self.get_param('blat_hostname')
            self.opts['blat_hostname'] = hostname
            self.opts['blat_port'] = port
            if port is None:
                self.opts['blat_port'] = random.randint(8000, 9500)
                utils.log(self.loggingName, 'info', 'Starting blat server on port %d on host %s.' % (self.opts['blat_port'], self.opts['blat_hostname']))    
        elif self.fncCmd == 'run':
            if not self.get_param('start_blat_server'):
                port = self.get_param('blat_port')
                hostname = self.get_param('blat_hostname')
                self.opts['blat_hostname'] = hostname
                if port is None:
                    utils.log(self.loggingName, 'debug', 'BreaKmer set to run and start_blat_server is set to False, but no blat server port is specified. Setting blat port to random value and starting blat server.')
                    self.opts['blat_port'] = random.randint(8000, 9500)
                else:
                    # Blat server is already running in this instance. Check it to make sure with a test blat.
                    self.opts['blat_port'] = int(self.get_param('blat_port'))
                    if self.check_blat_server():
                        return
                    else:
                        utils.log(self.loggingName, 'debug', 'Blat server with port %d and hostname %s did not pass test query. Please check specifications.' % (self.opts['blat_port'], self.opts['blat_hostname']))

        self.opts['reference_fasta_dir'] = os.path.split(self.opts['reference_fasta'])[0]

        # This makes the assumption that an existing blat server exists at this port.
        # Typically nice for debugging and --keep_blat_server was used.
        # TODO - make this more stable
        # if self.get_param('blat_port'):
        #     return

        ref_fasta_name = os.path.basename(self.opts['reference_fasta']).split(".fa")[0]

        # Check if 2bit is there.
        self.opts['blat_2bit'] = os.path.join(self.opts['reference_fasta_dir'], ref_fasta_name + ".2bit")
        if not os.path.exists(self.opts['blat_2bit']):
            utils.log(self.loggingName, 'info', 'Creating 2bit from %s reference fasta' % ref_fasta_name + ".fa")
            # Create 2bit requires faToTwoBit
            curdir = os.getcwd()
            os.chdir(self.opts['reference_fasta_dir'])
            cmd = '%s %s %s' % (self.opts['fatotwobit'], ref_fasta_name + ".fa", ref_fasta_name + ".2bit")
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            output, errors = p.communicate()
            os.chdir(curdir)

        curdir = os.getcwd()
        os.chdir(self.opts['reference_fasta_dir'])
        # Start gfServer, change dir to 2bit file, gfServer start localhost 8000 .2bit
        # self.opts['blat_port'] = random.randint(8000, 9500)
        # self.opts['blat_hostname'] = 'localhost'
        self.opts['gfserver_log'] = os.path.join(self.paths['output'], 'gfserver_%d.log' % self.opts['blat_port'])
        cmd = '%s -canStop -log=%s -stepSize=5 start %s %d %s &' % (self.opts['gfserver'], self.opts['gfserver_log'], self.opts['blat_hostname'], self.opts['blat_port'], ref_fasta_name + ".2bit")
        utils.log(self.loggingName, 'info', "Starting gfServer %s" % cmd)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        start_time = time.time()
        while not utils.server_ready(self.opts['gfserver_log']):
            new_time = time.time()
            wait_time = new_time - start_time
            if wait_time > 1000:
                utils.log(self.loggingName, 'error', 'gfServer wait time exceeded ~15 minutes, exiting')
                sys.exit(1)
            utils.log(self.loggingName, 'info', 'Waiting for blat gfServer to load reference seq')
            time.sleep(60)
        utils.log(self.loggingName, 'info', 'Server ready!')
        os.chdir(curdir)

    def get_kmer_size(self):
        return int(self.opts['kmer_size'])

    def get_min_segment_length(self, type):
        return int(self.opts[type + '_minseg_len'])

    def get_sr_thresh(self, type):
        """Get the threshold input for the number of reads that are required to
        support a structural variant event.
        """
        if type == 'min':
            return min(self.get_sr_thresh('trl'), self.get_sr_thresh('rearrangement'), self.get_sr_thresh('indel'))
        else:
            if type == 'trl':
                return int(self.opts['trl_sr_thresh'])
            elif type == 'rearrangement':
                return int(self.opts['rearr_sr_thresh'])
            elif type == 'indel':
                return int(self.opts['indel_sr_thresh'])

    def get_param(self, key, required=False):
        """Get the parameter value.

        If the parameer is required to be availale, then exit the program
        and throw an error.
        Args:
            key: The key in the opts dictionary to access the parameter value.
            required: Boolean value to indicate if the key should be required to
            be in the dictionary or not.
        Returns:
            The value of the parameter if it is found. If the parameter is
            required and not found the program will exit. If the parameter is
            not required and not found, it will return None.
        """
        value = None
        if key in self.opts:
            value = self.opts[key]
        elif required:
            utils.log(self.loggingName, 'error', 'Missing required parameter %s, exiting.' % key)
            sys.exit(1)
        return value

    def set_param(self, key, value):
        """Set the parameter value"""
        self.opts[key] = value
