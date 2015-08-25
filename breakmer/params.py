#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import logging
import random
import subprocess
import time
import pysam
import shutil
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
    """

    def __init__(self, fncCmd, arguments):
        """Initialize ParamManager class.

        Args:
            fncCmd (str):       The command to execute - run / prepare_reference_data / start_blat_server.
            arguments (dict):   The argparse dictionary object from the command line parameters.
        Returns:
            None
        Raises:
            None
        """

        self.__loggingName = 'breakmer.params'
        self.__opts = {}
        self.filter = None
        self.targets = {}
        self.paths = {}
        self.fncCmd = fncCmd
        self.__set_params(arguments)

    def __set_params(self, arguments):
        """Organize and format all input parameters into class variables to access
        later. Specific instances of parameters are checked and set. All parameters that are
        set are logged. The target objects are set along with the paths.

        Args:
            arguments (dict): The argparse dictionary object from the command line options.
        Returns:
            None
        Raises:
            None
        """

        self.__parse_opts(arguments)  # Parse the config file and command line parameters into the self.opts dictionary.
        utils.setup_logger(self.get_param('analysis_dir', True), 'breakmer')  # Create logging object.
        utils.log(self.loggingName, 'info', 'Setting up parameters')

        # Log all parameters passed in, warn for poor paths
        for paramKey, paramValue in self.opts.items():
            utils.log(self.loggingName, 'info', '%s = %s' % (paramKey, paramValue))

        self.__set_targets()
        self.paths['ref_data'] = os.path.abspath(os.path.normpath(self.opts['reference_data_dir'])) # Path to target reference sequence fast files.
        self.opts['reference_fasta_dir'] = os.path.split(self.opts['reference_fasta'])[0] # Path to genome fasta file.

        # If only preseting the reference data no need to continue.
        if self.fncCmd == 'prepare_reference_data':
            utils.log(self.loggingName, 'info', 'Preset reference data option set! Only the reference data directory will be setup.')
            return

        # Setup directories
        self.paths['analysis'] = os.path.abspath(os.path.normpath(self.opts['analysis_dir']))
        self.paths['output'] = os.path.join(self.paths['analysis'], 'output')
        if 'targets_dir' in self.opts:
            self.paths['targets'] = os.path.abspath(os.path.normpath(self.opts['targets_dir']))
        else:
            self.paths['targets'] = os.path.join(self.paths['analysis'], 'targets')

        # Create all the paths.
        for path in self.paths:
            utils.log(self.loggingName, 'info', 'Creating %s directory (%s)' % (path, self.paths[path]))
            if not os.path.exists(self.paths[path]):
                os.makedirs(self.paths[path])

        # If starting the blat server then return.
        if self.fncCmd == 'start_blat_server':
            utils.log(self.loggingName, 'info', 'Starting the blat server.')
            return

        self.__check_binaries()  # Check if Jellyfish and Cutadapt work.
        self.filter = resultfilter.ResultFilter(self.get_param('filterList'), self)  # Instantiate the filter class.
        self.__set_insertsize_thresh()  # Set the expected insert size threshold from the properly mapped read pairs.

    def __parse_opts(self, arguments):
        """Formats input parameters into self.opts dictionary. It first parses the configuration file and stores the key, values in the self.opts dictionary.
        It will exit with an error if the configuration file does not have lines in the proper format (i.e., key=value).
        It will also iterate through the command line paramaters and store the keys and values in the opts dictionary.
        A final check is performed for the required parameters depending on the parameters that have been passed.

        Args:
            arguments (dict):   The argparse dictionary object from the command line options.
        Returns:
            None
        Raises:
            None
        """

        for line in open(arguments.config_fn, 'rU'):
            line = line.strip()
            if line == '' or line.find('#') > -1:  # Allow for blank lines and comments
                continue
            linesplit = line.split("=")
            if len(linesplit) == 1:  # Make sure the lines in the configuration file are set properly.
                err_msg = 'Config line', line, ' not set correctly. Exiting.'
                print err_msg
                utils.log(self.loggingName, 'error', err_msg)
                sys.exit(1)
            else:
                key, value = linesplit
                self.__opts[key] = value  # Store key-value in opts dictionary.

        for opt in vars(arguments):
            self.__opts[opt] = vars(arguments)[opt]

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
        if self.fncCmd == 'prepare_reference_data':
            required = ['reference_data_dir', 'reference_fasta', 'targets_bed_file']

        for req in required:
            self.get_param(req, True)

    def __check_binaries(self):
        """Check the required binaries.
        There are six required binaries to perform the complete analysis (blat, gfserver,
        gfclient, fatotwobit, cutadapt, jellyfish). Each binary is checked whether
        the path provided in the configuration file has an executable file attached or
        if no path was provided that the binary is on the path. Cutadapt and Jellyfish
        are also tested using hardcoded data.

        Args:
            None
        Returns:
            None
        Raises:
            None
        """

        binaries = ('blat', 'gfserver', 'gfclient', 'fatotwobit', 'cutadapt', 'jellyfish')
        for binary in binaries:
            binaryPath = self.get_param(binary)
            if binaryPath is not None:
                binaryCheck = utils.which(binaryPath)  # Use the binary path specified in the config file.
            else:
                binaryCheck = utils.which(binary)  # Perform a which on the server to see if the binary is in the path.
                self.set_param(binary, binaryCheck)  # Store the result in the opts dictionary.
            if not binaryCheck:  # No binary found or specified. Throw an error.
                print 'Missing path/executable for', binary
                utils.log(self.loggingName, 'error', 'Missing path/executable for %s' % binary)
                sys.exit(1)
            utils.log(self.loggingName, 'info', '%s path = %s' % (binary, binaryCheck))
        utils.log(self.loggingName, 'info', 'All the required binaries have been checked successfully!')

        # Test cutadapt and jellyfish binaries
        testDir = os.path.join(self.paths['analysis'], 'bin_test')
        testFq = os.path.join(testDir, 'test.fq')
        if not os.path.exists(testDir):
            os.makedirs(testDir)

        fq_f = open(testFq, 'w')
        fq_f.write("@H91H9ADXX140327:1:2102:19465:23489/2\nCACCCCCACTGAAAAAGATGAGTATGCCTGCCGTGTGAACCATGTGACTTTACAATCTGCATATTGGGATTGTCAGGGAATGTTCTTAAAGATC\n+\n69EEEFBAFBFABCCFFBEFFFDDEEHHDGH@FEFEFCAGGCDEEEBGEEBCGBCCGDFGCBBECFFEBDCDCEDEEEAABCCAEC@>>BB?@C\n@H91H9ADXX140327:2:2212:12198:89759/2\nTCTTGTACTACACTGAATTCACCCCCACTGAAAAAGATGAGTATGCCTGCCGTGTGAACCATGTGACTTTACAATCTGCATATTGGGATTGTCAGGGA\n+\nA@C>C;?AB@BBACDBCAABBDDCDDCDEFCDDDDEBBFCEABCGDBDEEF>@GBGCEDGEDGCGFECAACFEGDFFGFECB@DFGCBABFAECEB?=")
        fq_f.close()

        cleanFq, rc = utils.test_cutadapt(testFq, self.opts['cutadapt'], self.opts['cutadapt_config_file'])
        if cleanFq:
            utils.log(self.loggingName, 'info', 'Test cutadapt ran successfully')
            jfish_prgm, rc = utils.test_jellyfish(self.opts['jellyfish'], cleanFq, testDir)
            if rc != 0:
                utils.log(self.loggingName, 'error', '%s unable to run successfully, exit code %s. Check installation and correct version.' % (jfish_prgm, str(rc)))
                sys.exit(1)
            else:
                utils.log(self.loggingName, 'info', 'Test jellyfish ran successfully')
        else:
            utils.log(self.loggingName, 'error', 'Cutadapt failed to run, exit code %s. Check installation and version.' % str(rc))
            sys.exit(1)
        shutil.rmtree(testDir)  # Remove the test directory.

    def __set_insertsize_thresh(self):
        """Store the insert sizes for a small number of "properly mapped" reads
        and determine an upperbound cutoff to use to determine discordantly mapped read
        pairs.

        Args:
            None
        Returns:
            None
        Raises:
            None
        """

        nSampleReads = 100000
        bamF = pysam.Samfile(self.get_param('sample_bam_file'), 'rb')
        testReads = bamF.fetch()
        insertSizes = []
        readIter = 0
        for read in testReads:
            proper_map = read.flag == 83 or read.flag == 99 and read.mapq > 0
            if read.is_read1 and proper_map:  # Sample the read and store the insert size to its partner.
                readIter += 1
                insertSizes.append(abs(read.tlen))
                if 'readLen' not in self.opts:  # Store the read length if it is not already stored.
                    self.opts['readLen'] = read.rlen
            if readIter == nSampleReads:
                break
        isMedian = utils.median(insertSizes)
        isSD = utils.stddev(utils.remove_outliers(insertSizes))  # Calculate the standard deviation of the sample read pairs insert sizes.
        self.opts['insertsize_thresh'] = isMedian + (5 * isSD)  # Set the threshold to be median + 5 standard deviations.

    def __set_targets(self):
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

        Store the target information in the self.target dictionary with the name as the key
        and a list of tuples of interval genomic locations as the values.
        self.target[gene_name] = [(chrom, start_bp, end_bp, name, feature),...]

        Args:
            None
        Returns:
            None
        Raises:
            None
        """

        geneList = self.get_param('gene_list')
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
            if name.upper() not in self.targets:  # All the names are converted to uppercase.
                self.targets[name.upper()] = []
            self.targets[name.upper()].append((chrm, int(bp1), int(bp2), name, feature))
        utils.log(self.loggingName, 'info', '%d targets' % len(self.targets))

    def __check_blat_server(self):
        """Run a test query on the specified blat server to make sure it is running. 

        Args:
            None
        Returns:
            serverSuccess (boolean): Indicates whether the test ran without errors.
        Raises:
            None
        """

        testDir = os.path.join(self.paths['analysis'], 'blatserver_test')
        test_fa_fn = os.path.join(testDir, 'test.fa')
        if not os.path.exists(testDir):
            os.makedirs(testDir)
        test_fa = open(test_fa_fn, 'w')
        test_fa.write('>test\nCCAAGGGAGACTTCAAGCAGAAAATCTTTAAGGGACCCTTGCATAGCCAGAAGTCCTTTTCAGGCTGATGTACATAAAATATTTAGTAGCCAGGACAGTAGAAGGACTGAAGAGTGAGAGGAGCTCCCAGGGCCTGGAAAGGCCACTTTGTAAGCTCATTCTTG')
        test_fa.close()

        resultFn = os.path.join(testDir, 'blatserver_test.psl')
        cmd = '%s -t=dna -q=dna -out=psl -minScore=20 -nohead %s %d %s %s %s' % (self.get_param('gfclient'), self.get_param('blat_hostname'), self.get_param('blat_port'), self.get_param('reference_fasta_dir'), test_fa_fn, resultFn)
        utils.log(self.loggingName, 'info', 'Blat server test system command %s' % cmd)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, errors = p.communicate()
        utils.log(self.loggingName, 'info', 'Realignment output file %s' % resultFn)
        serverSuccess = True
        if errors != '':
            serverSuccess = False
            utils.log(self.loggingName, 'info', 'Realignment errors %s' % errors)
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

        if self.fncCmd == 'prepare_reference_data':  # Do not start blat server for this function.
            return
        elif self.fncCmd == 'start_blat_server':
            port = self.get_param('blat_port')
            hostname = self.get_param('blat_hostname')
            self.opts['blat_hostname'] = hostname
            self.opts['blat_port'] = port
            # If no port is specified for this function, then randomly select a port between 8000-9500.
            if port is None:
                self.opts['blat_port'] = random.randint(8000, 9500)
                utils.log(self.loggingName, 'info', 'Starting blat server on port %d on host %s.' % (self.opts['blat_port'], self.opts['blat_hostname']))    
        elif self.fncCmd == 'run':  # Start the blat server if it is not already running.
            if not self.get_param('start_blat_server'):  # Start blat server option is not set. Check that one is running, if not, start it.
                port = self.get_param('blat_port')
                hostname = self.get_param('blat_hostname')
                self.opts['blat_hostname'] = hostname
                if port is None:  # No port is specified for a server that should be running. It will start a new one on a random numbered port.
                    utils.log(self.loggingName, 'debug', 'BreaKmer set to run and start_blat_server is set to False, but no blat server port is specified. Setting blat port to random value and starting blat server.')
                    self.opts['blat_port'] = random.randint(8000, 9500)
                else:  # Blat server is already running in this instance. Check it to make sure with a test blat.
                    self.opts['blat_port'] = int(self.get_param('blat_port'))
                    if self.__check_blat_server():  # Both port and hostname are specified. Check that the server is running.
                        return
                    else:
                        utils.log(self.loggingName, 'debug', 'Blat server with port %d and hostname %s did not pass test query. Please check specifications.' % (self.opts['blat_port'], self.opts['blat_hostname']))

        self.opts['reference_fasta_dir'] = os.path.split(self.opts['reference_fasta'])[0]
        ref_fasta_name = os.path.basename(self.opts['reference_fasta']).split(".fa")[0]

        self.opts['blat_2bit'] = os.path.join(self.opts['reference_fasta_dir'], ref_fasta_name + ".2bit")
        if not os.path.exists(self.opts['blat_2bit']):  # Create 2bit file to use for running the blat server.
            utils.log(self.loggingName, 'info', 'Creating 2bit from %s reference fasta' % ref_fasta_name + ".fa")
            curdir = os.getcwd()
            os.chdir(self.opts['reference_fasta_dir'])
            cmd = '%s %s %s' % (self.opts['fatotwobit'], ref_fasta_name + ".fa", ref_fasta_name + ".2bit")
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            output, errors = p.communicate()
            os.chdir(curdir)

        curdir = os.getcwd()
        os.chdir(self.opts['reference_fasta_dir'])
        # Start gfServer, change dir to 2bit file, gfServer start localhost 8000 .2bit
        self.opts['gfserver_log'] = os.path.join(self.paths['output'], 'gfserver_%d.log' % self.opts['blat_port'])
        cmd = '%s -canStop -log=%s -stepSize=5 start %s %d %s &' % (self.opts['gfserver'], self.opts['gfserver_log'], self.opts['blat_hostname'], self.opts['blat_port'], ref_fasta_name + ".2bit")
        utils.log(self.loggingName, 'info', "Starting gfServer %s" % cmd)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        start_time = time.time()

        while not utils.server_ready(self.opts['gfserver_log']):  # Wait for the blat server to initiate. Timeout if it has not started in 15 minutes.
            new_time = time.time()
            wait_time = new_time - start_time
            if wait_time > 1000:
                utils.log(self.loggingName, 'error', 'gfServer wait time exceeded ~15 minutes, exiting')
                sys.exit(1)
            utils.log(self.loggingName, 'info', 'Waiting for blat gfServer to load reference seq')
            time.sleep(60)
        utils.log(self.loggingName, 'info', 'Server ready!')
        os.chdir(curdir)

    def get_target_names(self):
        """Get a list of target names.

        Args:
            None
        Returns:
            A list of the target region names that were defined from the input bed file (list).
        """

        return self.targets.keys()

    def get_target_intervals(self, targetName):
        """Return the stored intervals for a specific target.
        """

        if targetName in self.targets:
            return self.targets[targetName]
        else:
            utils.log(self.loggingName, 'debug', '%s target name not in target dictionary.' % targetName)
            sys.exit(1)

    def get_kmer_size(self):
        """Get the input kmer size. This should be an integer value.

        Args:
            None
        Returns:
            kmer size (int): Kmer size that was input to use.
        Raises:
            TypeError when the kmer_size is not an integer.
        """

        try:
            int(self.opts['kmer_size'])
        except ValueError:
            print 'The specified kmer size is not an integer.'
            raise
        else:
            return int(self.opts['kmer_size'])

    def get_min_segment_length(self, type):
        """Get the input segment length limit. This should be an integer value.

        Args:
            type (str): The variant type to get the minimum segment length - trl / rearr
        Returns:
            min_seg_len (int)
        Raises:
            TypeError when the value is not an integer.
        """

        try:
            int(self.opts[type + '_minseg_len'])
        except ValueError:
            print 'The specified minsegment limit is not an integer.'
            raise
        else:
            return int(self.opts[type + '_minseg_len'])

    def get_sr_thresh(self, type):
        """Get the threshold input for the number of reads that are required to
        support a structural variant event.

        Args:
            type (str): The variant type to get the read support threshold.
        Returns:
            Integer of the split read threshold for specific events.
        Raises:
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
        """Get the parameter value in the self.opts dictionary.

        If the parameer is required to be availale, then exit the program
        and throw an error.
        Args:
            key (str): The key in the opts dictionary to access the parameter value.
                       required: Boolean value to indicate if the key should be required to
                       be in the dictionary or not.
        Returns:
            value (int, str, boolean): The value of the parameter if it is found. If the parameter is
                                       required and not found the program will exit with error. If the parameter is
                                       not required and not found, it will return None.
        Raises:
            None
        """

        value = None
        if key in self.__opts:
            value = self.__opts[key]
        elif required:
            utils.log(self.loggingName, 'error', 'Missing required parameter %s, exiting.' % key)
            sys.exit(1)
        return value

    def set_param(self, key, value):
        """Set the parameter value in the self.opts dict.

        Args:
            key (str):              Dictionary key
            value (int/str/boolean):  Value to store
        Returns:
            None
        Raises:
            None
        """

        self.opts[key] = value
