#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
import pysam
import shutil
import subprocess
import breakmer.utils as utils
import breakmer.processor.bam_handler as bam_handler
import breakmer.assembly.assembler as assembly

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


'''
target.py module contains classes and functions to handle target-level analysis.

For each target region:
1. Set reference data
2. Extract structural variant reads and clean them
3. Kmer subtraction
4. Assemble extracted reads and realign.
    - breakmer.assembly.assembler.init_assembly()
5. Make calls and output results.
    1. breakmer.assembly.contig.make_calls()

Classes:
    - TargetManager - Interface for target region functions and values.
    - Variation - A "subclass" to TargetManager. It is designed to manage the extracted read data for a target.

Functions:
    - load_kmers
'''


def load_kmers(fns, kmerDict):
    """Iterate through the Jellyfish kmer flat files and store them in the kmerDict dictionary.

    Store the kmer sequence string as the key and the count of the number of reads
    containing it as the value. The kmerDict is passed in by reference and any changes
    should persist even though it is not returned.

    Args:
        fns (str):       Filenames of the kmer flat files.
        kmerDict (dict): kmer sequence as key, frequency of sequence as values.

    Returns:
        None
    """

    if not fns:
        return kmerDict
    fns = fns.split(',')
    for fn in fns:  # Iterate through all the jellyfish kmer files and store the kmer as key and count as value.
        for line in open(fn, 'rU'):
            line = line.strip()
            mer, count = line.split()
            if mer not in kmerDict:
                kmerDict[mer] = 0
            kmerDict[mer] += int(count)


class Variation:
    """This class handles the storage and interaction of all the variant reads that could
    be contributing to the support of a structural variant.

    Attributes:
        params (ParamManager):    Parameters for breakmer analysis.
        loggingName (str):        Module name for logging file purposes.
        var_reads (dict):         Tumor/normal sample variation read objects (breakmer.process.bam_handler.VariantReadTracker).
        cleaned_read_recs (dict): Cleaned reads for tumor/normal samples.
        files (dict):             File paths for tumor/normal analysis.
        kmer_clusters (list):     A list of Contig objects that are assembled from tumor extracted reads.  
        kmers (dict):             A set of processed kmer sequences from tumor/normal/reference sequences.
        results (list):           
        discReadClusters (dict):
        discReadFormatted (list):
    """

    def __init__(self, params):
        self.loggingName = 'breakmer.processor.target'
        self.params = params
        self.var_reads = {}
        self.cleaned_read_recs = None
        self.kmer_clusters = []
        self.kmers = {}
        self.results = []
        self.files = {}
        # self.svs = {}
        self.discReadClusters = {}
        self.discReadFormatted = []

    def setup_cleaned_reads(self, sampleType):
        """Initiate the cleaned_read_recs dictionary for sample or normal data.

        Args:
            sampleType (str): Indicates the sample type - sv or normal

        Returns:
            None
        """

        if not self.cleaned_read_recs:
            self.cleaned_read_recs = {}
        self.cleaned_read_recs[sampleType] = None

    def get_var_reads(self, sampleType):
        """
        """

        return self.var_reads[sampleType]

    def clear_sv_reads(self, sampleType):
        """Wrapper function to call the clear_sv_reads() function in bam_handler.VariationReadTracker.

        Clears the data structure containing the extracted reads for this target.

        Args:
            sampleType (str): A string to indicate if the data being processed is from
                              tumor/normal.

        Returns:
            None
        """

        self.var_reads[sampleType].clear_sv_reads()

    def clear_cleaned_reads(self):
        """
        """

        self.cleaned_read_recs = None

    def continue_analysis_check(self, type):
        """
        """

        check = True
        if len(self.cleaned_read_recs[type]) == 0:
            check = False
        return check

    def get_sv_reads(self, sampleType):
        """
        """

        return self.var_reads[sampleType].sv

    def add_result(self, result):
        """
        """

        self.results.append(result)

    def set_var_reads(self, sampleType, bamFile, chrom, start, end, regionBuffer):
        """Pass in the bam file and the coordinates for a target region and extract the
        reads that indicate a structural variant is present.

        The region is extended by the size specified in the regionBuffer, both upstream and downstream
        from the start and end inputs, respectively.

        A VariantReadTracker object is returned and stored in a dictionary, where the key is the sampleType (sv or norm).
        The extracted sequences are written to multiple files:
        1. Fastq
        2. fasta
        3. Bam

        Args:
            sampleType (str):   'sv' or 'norm' to indicate tumor or normal sample being processed.
            bamFile (str):      The full path to the bam file containing the reads to be extracted.
            chrom (str):        The chromosome of the target region to extract.
            start (int):        The start coordinate of the target region.
            end (int):          The end coordinate of the target region.
            regionBuffer (int): The additional base pairs to extract upstream and downstream of the target region.

        Returns:
            None
        """

        # Get VariantReadTracker object from bam_handler module and extract reads.
        self.var_reads[sampleType] = bam_handler.get_variant_reads(bamFile, chrom, start - regionBuffer, end - regionBuffer, self.params.get_param('insertsize_thresh'))
        # Iterate through reads that are not perfectly aligned and store necessary information for downstream analysis.
        # Store the reads with softclipped sequences that are high quality in VariantReadTracker.sv dictionary.
        self.var_reads[sampleType].check_clippings(self.params.get_kmer_size(), start, end)

        # Write the bam, fastq, and fasta files with the extracted reads.
        svBam = None
        if sampleType == 'sv':
            svBam = pysam.Samfile(self.files['sv_bam'], 'wb', template=pysam.Samfile(bamFile, 'rb'))
        readsFq = open(self.files['%s_fq' % sampleType], 'w')
        scFa = open(self.files['%s_sc_unmapped_fa' % sampleType], 'w')
        # Write all the stored sequences into files.
        self.var_reads[sampleType].write_seqs(scFa, readsFq, svBam, self.params.get_kmer_size())
        readsFq.close()
        scFa.close()

        # Close the bam file, sort and index.
        if sampleType == 'sv':
            svBam.close()
            utils.log(self.loggingName, 'info', 'Sorting bam file %s to %s' % (self.files['sv_bam'], self.files['sv_bam_sorted']))
            pysam.sort(self.files['sv_bam'], self.files['sv_bam_sorted'].replace('.bam', ''))
            utils.log(self.loggingName, 'info', 'Indexing sorted bam file %s' % self.files['sv_bam_sorted'])
            pysam.index(self.files['sv_bam_sorted'])

    def setup_read_extraction_files(self, sampleType, dataPath, name):
        """Create file names to store the extracted reads.

        This creates four files (for tumor samples):
        1. fastq with extracted reads = sv_fq or normal_fq
        2. fasta file with softclipped sequences = sv_sc_unmapped_fa
        3. bam file with extracted reads = sv_bam
        4. sorted bam file with extracted reads = sv_bam_sorted

        Args:
            sampleType (str): The type of input data - sv / normal
            dataPath (str):   The path to the data files for this target.
            name (str):       The target name.

        Returns:
            None
        """

        # Store extracted reads in <data_path>/<target_name>_<type>_reads.fastq
        self.files['%s_fq' % sampleType] = os.path.join(dataPath, name + '_%s_reads.fastq' % sampleType)
        # Store softclipped sequences in a fasta file <data_path>/<target_name>_<type>_sc_seqs.fa
        self.files['%s_sc_unmapped_fa' % sampleType] = os.path.join(dataPath, name + '_%s_sc_seqs.fa' % sampleType)

        if sampleType == 'sv':
            # Store variant reads in bam formatted file <data_path>/<target_name>_sv_reads.bam
            self.files['sv_bam'] = os.path.join(dataPath, name + '_sv_reads.bam')
            # Store variant reads in sorted bam file
            self.files['sv_bam_sorted'] = os.path.join(dataPath, name + '_sv_reads.sorted.bam')

    def clean_reads(self, dataPath, name, sampleType):
        """Trim adapter sequences from the extracted reads, format and organize
        the cleaned reads into new files.

        Cutadapt is run to trim the adapter sequences from the sequence reads to
        remove any 'noise' from the assembly process. The cleaned reads output
        from cutadapt are then reprocessed to determine if the softclipped sequences
        were trimmed off or not to further filter out reads.

        The softclipped sequences that remain are stored and a new fastq file is written.

        Args:
            dataPath (str):  The path to the data files for this target.
            name (str):      The target name.
            type (str):      A string indicating a tumor ('sv') or normal ('norm') sample being processed.

        Returns:
            check (boolean): A boolean to indicate whether the are any reads left after
                             cleaning is complete.
        """

        cutadapt = self.params.get_param('cutadapt')  # Cutadapt binary
        cutadaptConfigFn = self.params.get_param('cutadapt_config_file')  # Get the config file containing the adapter sequences to remove.
        utils.log(self.loggingName, 'info', 'Cleaning reads using %s with configuration file %s' % (cutadapt, cutadaptConfigFn))
        self.files['%s_cleaned_fq' % sampleType] = os.path.join(dataPath, name + '_%s_reads_cleaned.fastq' % sampleType)
        utils.log(self.loggingName, 'info', 'Writing clean reads to %s' % self.files['%s_cleaned_fq' % sampleType])
        output, errors = utils.run_cutadapt(cutadapt, cutadaptConfigFn, self.files['%s_fq' % sampleType], self.files['%s_cleaned_fq' % sampleType], self.loggingName)

        self.setup_cleaned_reads(sampleType)
        self.files['%s_cleaned_fq' % sampleType], self.cleaned_read_recs[sampleType] = utils.get_fastq_reads(self.files['%s_cleaned_fq' % sampleType], self.get_sv_reads(sampleType))
        self.clear_sv_reads(sampleType)
        check = self.continue_analysis_check(sampleType)
        utils.log(self.loggingName, 'info', 'Clean reads exist %s' % check)
        return check

    def set_reference_kmers(self, targetRefFns):
        """Populate the self.kmers['ref'] dictionary with the kmer sequences in the 
        target reference sequence file.

        Args:
            targetRefFns (list): A list of strings, which are the paths to the fasta files
                                 containing the forward and reverse sequences of the target reference.

        Returns:
            None
        """

        self.kmers['ref'] = {}
        for i in range(len(targetRefFns)):
            utils.log(self.loggingName, 'info', 'Indexing kmers for reference sequence %s' % targetRefFns[i])
            self.get_kmers(targetRefFns[i], self.kmers['ref'])

    def set_sample_kmers(self):
        """Populate the self.kmers dictionary for 'case' and 'case_sc' keys with the kmers
        within the extracted, cleaned tumor sequences and softclipped and unmapped sequences.

        Args:
            None

        Returns:
            None
        """

        utils.log(self.loggingName, 'info', 'Indexing kmers for sample sequence %s' % self.files['sv_cleaned_fq'])
        self.kmers['case'] = {}
        self.kmers['case_sc'] = {}
        self.get_kmers(self.files['sv_cleaned_fq'], self.kmers['case'])
        self.get_kmers(self.files['sv_sc_unmapped_fa'], self.kmers['case_sc'])

    def get_kmers(self, seqFn, kmerDict):
        """Generic function to run jellyfish on a set of sequences.

        The kmerDict is passed in by reference, which means it will be modified
        even when it is not returned by the function.

        1. Run jellyfish to generate kmers of specified size (kmer_size) and get the
           the flat file containing all the kmers with their frequencies of occurance in the
           reads that were input into Jellyfish.
        2. Store the kmer sequence as keys in the kmerDict and the frequency of the sequence
           as the value.

        Args:
            seqFn (str):     Full path to the fastq or fasta file containing sequences to
                             be split into kmer-sized sequences.
            kmerDict (dict): A dictionary with kmer sequence as key, frequency of sequence as values.

        Returns:
            None
        """

        # Load the kmers into the kmer dictionary based on keyStr value.
        load_kmers(utils.run_jellyfish(seqFn, self.params.get_param('jellyfish'), self.params.get_kmer_size()), kmerDict)

    def compare_kmers(self, kmerPath, name, readLen, targetRefFns):
        """

        Args:
            kmerPath (str):      Path to the kmer directory to store the files.
            name (str):          Target name.
            readLen (int):       Read length.
            targetRefFns (list): A list of the target reference fasta sequence files.

        Returns:
            None
        """

        # Store the kmers from the reference sequence in the target region in
        # self.kmers['ref'].
        self.set_reference_kmers(targetRefFns)

        # Store the extracted sequences in the tumor sample into
        # self.kmers['case'] and self.kmers['case_sc'] dictionaries.
        self.set_sample_kmers()

        # Merge the kmers from the cleaned sample sequences and the unmapped and softclipped sequences.
        scKmers = set(self.kmers['case'].keys()) & set(self.kmers['case_sc'].keys())
        # Take the difference from the reference kmers.
        sampleOnlyKmers = list(scKmers.difference(set(self.kmers['ref'].keys())))
        # Add normal sample kmers if available.
        if self.params.get_param('normal_bam_file'):
            normKmers = {}
            self.get_kmers(self.files['norm_cleaned_fq'], normKmers)
            sampleOnlyKmers = list(set(sampleOnlyKmers).difference(set(normKmers.keys())))

        # Write case only kmers out to file.
        self.files['sample_kmers'] = os.path.join(kmerPath, name + "_sample_kmers.out")
        sample_kmer_fout = open(self.files['sample_kmers'], 'w')
        kmer_counter = 1
        self.kmers['case_only'] = {}  # Populate the case_only dictioanry with the kmer subtraction.
        for mer in sampleOnlyKmers:
            sample_kmer_fout.write("\t".join([str(x) for x in [mer, str(self.kmers['case'][mer])]]) + "\n")
            self.kmers['case_only'][mer] = self.kmers['case'][mer]
        sample_kmer_fout.close()

        # Clean out data structures.
        self.kmers['ref'] = {}
        self.kmers['case'] = {}
        self.kmers['case_sc'] = {}

        utils.log(self.loggingName, 'info', 'Writing %d sample-only kmers to file %s' % (len(self.kmers['case_only']), self.files['sample_kmers']))
        self.files['kmer_clusters'] = os.path.join(kmerPath, name + "_sample_kmers_merged.out")
        utils.log(self.loggingName, 'info', 'Writing kmer clusters to file %s' % self.files['kmer_clusters'])

        # Initialize the assembly using all the extracted reads and tumor only kmers.
        # Pass in:
        #      1. self.kmers['case_only'] - dictionary kmer sequence as key, count as value.
        #      2. self.cleaned_read_recs['sv'] - dictionary read sequence as key, list of fq_read objects as value.
        #      3. kmer_size
        #      4. readcount threshold for contig
        #      5. read length
        # Returns a list of Contig objects
        self.kmers['clusters'] = assembly.init_assembly(self.kmers['case_only'], self.cleaned_read_recs['sv'], self.params.get_kmer_size(), self.params.get_sr_thresh('min'), readLen)

        # Clear the data structures
        self.clear_cleaned_reads()
        self.kmers['case_only'] = {}

    def get_disc_reads(self):
        """
        """

        return self.var_reads['sv'].get_disc_reads()

    def write_results(self, outputPath, targetName):
        """
        """

        if len(self.results) > 0:
            resultFn = os.path.join(outputPath, targetName + "_svs.out")
            utils.log(self.loggingName, 'info', 'Writing %s result file %s' % (targetName, resultFn))
            resultFile = open(resultFn, 'w')
            for i, result in enumerate(self.results):
                headerStr, formattedResultValuesStr = result.get_formatted_output_values()
                if i == 0:
                    resultFile.write(headerStr + '\n')
                resultFile.write(formattedResultValuesStr + '\n')
            resultFile.close()
        if len(self.discReadClusters) > 0:
            resultFn = os.path.join(outputPath, targetName + "_discreads.out")
            utils.log(self.loggingName, 'info', 'Writing %s discordant read cluster result file %s' % (targetName, resultFn))
            resultFile = open(resultFn, 'w')
            for i, discReadRes in enumerate(self.discReadFormatted):
                headerStr, outStr = discReadRes
                if i == 0:
                    resultFile.write(headerStr + '\n')
                resultFile.write(outStr + '\n')
            resultFile.close()

    def get_formatted_output(self):
        """
        """

        formattedResultsDict = {'contigs': [], 'discreads': []}
        if len(self.results) > 0:
            for i, result in enumerate(self.results):
                formattedResultsDict['contigs'].append(result.get_formatted_output_values())
        if len(self.discReadClusters) > 0:
            for i, discReadRes in enumerate(self.discReadFormatted):
                formattedResultsDict['discreads'].append(discReadRes)
        return formattedResultsDict

    def cluster_discreads(self, targetName, targetChrom):
        """
        """

        self.discReadClusters = self.var_reads['sv'].cluster_discreads()
        self.discReadFormatted = []
        headerStr = '\t'.join(['Target_name', 'sv_type', 'cluster_id', 'left_breakpoint_estimate', 'right_breakpoint_estimate', 'strands', 'discordant_readpair_count', 'cluster_distance'])
        for key in self.discReadClusters:
            readCount = self.discReadClusters[key]['readCount']
            k1, k2, k3, c1, c2 = key.split('|')
            checkCount = readCount
            clusterDist = '0'
            if k1 == 'inter':
                checkCount = self.discReadClusters[key]['interClusterCount']
            if checkCount < self.params.get_param('discread_only_thresh'):
                continue
            svType = 'inter-chromosomal'
            lChrom = 'chr' + targetChrom.replace('chr', '')
            if k1 == 'inter':
                rChrom = 'chr' + k2.replace('chr', '')
            elif k1 == 'intra':
                svType = 'intra-chromosomal_' + k2
                rChrom = lChrom
                clusterDist = abs(self.discReadClusters[key]['leftBrkpt'] - self.discReadClusters[key]['rightBrkpt'])
            lStrand, rStrand = k3.split(':')
            lBrkpt = self.discReadClusters[key]['leftBrkpt']
            rBrkpt = self.discReadClusters[key]['rightBrkpt']
            clusterId = targetName + '_' + str(self.discReadClusters[key]['clusterId'])
            outStr = '\t'.join([targetName, svType, clusterId, lChrom + ':' + str(lBrkpt), rChrom + ':' + str(rBrkpt), lStrand + ',' + rStrand, str(readCount), str(clusterDist)])
            self.discReadFormatted.append((headerStr, outStr))


class TargetManager(object):
    """TargetManager class handles all the high level information relating to a target.
    The analysis is peformed at the target level, so this class contains all the information
    necessary to perform an independent analysis.

    Call the setup() function from the init() function. This will iterate over the input interval
    regions and sets the chromosome of the target region and the max(interval coordinates) as the end
    and the min(interval coordinates) as the start.

    This will setup all the directories and files for this target.

    Attributes:
        params (ParamManager):      Parameters for breakmer analysis.
        loggingName (str):          Module name for logging file purposes.
        name (str):                 Target name specified in the input bed file.
        chrom (str):                Chromosome ID as specified in the input bed file.
        start (int):                Genomic position for the target region (minimum value among all intervals).
        end (int):                  Genomic position for the target region (maximum value among all intervals).
        paths (dict):               Contains the analysis paths for this target.
        files (dict):               Contains the paths to file names needed for analysis.
        read_len (int):             Length of a single read.
        variation (Variation):      Stores data for variants identified within the target.
        regionBuffer (int):         Base pairs to add or subtract from the target region end and start locations.
    """

    def __init__(self, name, params):
        self.loggingName = 'breakmer.processor.target'
        self.params = params
        self.name = name
        self.chrom = None
        self.start = None
        self.end = None
        self.paths = {}
        self.files = {}
        self.readLen = int(params.get_param('readLen'))
        self.variation = Variation(params)
        self.regionBuffer = 200
        self.setup()

    @property
    def values(self):
        """Return the defined features of this target
        """

        return (self.chrom, self.start, self.end, self.name, self.get_target_intervals(), self.regionBuffer)

    @property
    def fnc(self):
        """Return the function of the program.
        """

        return self.params.fncCmd

    def setup(self):
        """Setup the TargetManager object with the input params.

        Define the location (chrom, start, end), file paths, directory paths, and name.

        Args:
            None

        Returns:
            None
        """

        # Define the target boundaries based on the intervals input.
        # The target start is the minimum start of the intervals and the end
        # is the maximum end of the intervals.
        intervals = self.params.get_target_intervals(self.name)
        for values in intervals:
            chrom, start, end = values[0], int(values[1]), int(values[2])
            if self.chrom is None:
                self.chrom = chrom
            if self.start is None:
                self.start = start
            elif start < self.start:
                self.start = start
            if self.end is None:
                self.end = end
            elif end > self.end:
                self.end = end

        # Create the proper paths for the target analysis.
        '''
        Each target analyzed has a set of directories associated with it.
        targets/
            <target name>/
                data/
                contigs/
                kmers/

        There is separate directory for each target in the output directory.
        output/
            <target name>/
        '''
        self.add_path('ref_data', os.path.join(self.params.paths['ref_data'], self.name))
        if self.params.fncCmd == 'run':
            self.add_path('base', os.path.join(self.params.paths['targets'], self.name))
            self.add_path('data', os.path.join(self.paths['base'], 'data'))
            self.add_path('contigs', os.path.join(self.paths['base'], 'contigs'))
            self.add_path('kmers', os.path.join(self.paths['base'], 'kmers'))
            self.add_path('output', os.path.join(self.params.paths['output'], self.name))

        '''
        Each target has reference files associated with it.
        <ref_data_dir>/
            <target_name>/
                <target_name>_forward_refseq.fa
                <target_name>_reverse_refseq.fa
                <target_name>_forward_refseq.fa_dump
                <target_name>_reverse_refseq.fa_dump
        '''
        self.files['target_ref_fn'] = [os.path.join(self.paths['ref_data'], self.name + '_forward_refseq.fa'), os.path.join(self.paths['ref_data'], self.name + '_reverse_refseq.fa')]
        # ref_fa_marker_f = open(os.path.join(self.paths['ref_data'], '.reference_fasta'), 'w')
        # ref_fa_marker_f.write(self.params.get_param('reference_fasta'))
        # ref_fa_marker_f.close()
        self.files['ref_kmer_dump_fn'] = [os.path.join(self.paths['ref_data'], self.name + '_forward_refseq.fa_dump'), os.path.join(self.paths['ref_data'], self.name + '_reverse_refseq.fa_dump')]

    def add_path(self, key, path):
        """Utility function to create all the output directories.

        Args:
            key (str):  String value to store the file path value.
            path (str): File path value.

        Returns:
            None
        """

        utils.log(self.loggingName, 'info', 'Creating %s %s path (%s)' % (self.name, key, path))
        self.paths[key] = path
        if not os.path.exists(self.paths[key]):
            os.makedirs(self.paths[key])

    def set_ref_data(self):
        """Write the reference sequence to a fasta file for this specific target if it does not
        exist.

        This will write the forward sequence, which should be what is contained in the reference
        sequence file, as well as the reverse sequence.

        If blastn is specified in the configuration file, it will setup the necessary blastn files as
        well. This is currently under development and needs more testing.

        The target reference files stored in:
        <ref_data_dir>/
            <target_name>/
                <target_name>_forward_refseq.fa
                <target_name>_reverse_refseq.fa
                <target_name>_forward_refseq.fa_dump
                <target_name>_reverse_refseq.fa_dump

        Args:
            None

        Returns:
            None
        """

        # Write reference fasta file if needed.
        for i in range(len(self.files['target_ref_fn'])):
            fn = self.files['target_ref_fn'][i]
            direction = "forward" if fn.find("forward") != -1 else "reverse"
            utils.log(self.loggingName, 'info', 'Extracting refseq sequence and writing %s' % fn)
            utils.extract_refseq_fa(self.values, self.paths['ref_data'], self.params.get_param('reference_fasta'), direction, fn)

        blastn = self.params.get_param('blast')  # If using blatn for target realignment, the db must be available.
        if blastn is not None:
            # Check if blast db files are available for each target.
            if not os.path.isfile(self.files['target_ref_fn'][0] + '.nin'):
                makedb = os.path.join(os.path.split(blastn)[0], 'makeblastdb')  # Create blast db
                cmd = "%s -in %s -dbtype 'nucl' -out %s" % (makedb, self.files['target_ref_fn'][0], self.files['target_ref_fn'][0])
                utils.log(self.loggingName, 'info', 'Creating blast db files for target %s with reference file %s' % (self.name, self.files['target_ref_fn'][0]))
                p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                output, errors = p.communicate()
                if errors != '':
                    utils.log(self.loggingName, 'debug', 'Failed to make blast db files using reference file %s' % self.files['target_ref_fn'][0])

    def find_sv_reads(self):
        """Entry function to extract sequence reads from sample or normal bam file.
        It extracts and cleans the sample reads from the target region that may
        be used to build a variant contig.

        1. Extract bam reads for tumor and normal (if specified in config).
        2. Clean reads

        Args:
            None

        Returns:
            check (boolean): Variable to determine if the analysis should continue. It is
                             False when there are no reads extracted or left after cleaning
                             and True when there are.
        """

        self.extract_bam_reads('sv')  # Extract variant reads.
        if self.params.get_param('normal_bam_file'):  # Extract reads from normal sample, if input.
            self.extract_bam_reads('norm')
            self.clean_reads('norm')
        check = True
        if not self.clean_reads('sv'):  # Check if there are any reads left to analyze after cleaning.
            shutil.rmtree(self.paths['output'])  # Remove the output directory since there is nothing to analyze
            check = False
        return check

    def extract_bam_reads(self, sampleType):
        """Wrapper for Variation extract_bam_reads function for both tumor and normal
        samples if both are input.

        Args:
            sampleType (str): Indicates a tumor ('sv') or normal ('norm') sample being processed.

        Returns:
            None
        """

        # Create the file paths for the files that will be created from the read extraction.
        self.variation.setup_read_extraction_files(sampleType, self.paths['data'], self.name)
        bamType = 'sample'
        if sampleType == 'norm':
            bamType = 'normal'
        bamFile = self.params.get_param('%s_bam_file' % bamType)
        utils.log(self.loggingName, 'info', 'Extracting bam reads from %s to %s' % (bamFile, self.variation.files['%s_fq' % sampleType]))
        self.variation.set_var_reads(sampleType, bamFile, self.chrom, self.start, self.end, self.regionBuffer)

    def clean_reads(self, sampleType):
        """Wrapper for Variation clean_reads function.

        Args:
            sampleType (str): A string indicating a tumor ('sv') or normal ('norm') sample being processed.

        Returns:
            boolean:          A boolean to indicate whether the are any reads left after
                              cleaning is complete.
        """

        return self.variation.clean_reads(self.paths['data'], self.name, sampleType)

    def compare_kmers(self):
        """Wrapper for Variation compare_kmers function.

        Args:
            None

        Returns:
            None
        """

        self.variation.compare_kmers(self.paths['kmers'], self.name, self.readLen, self.files['target_ref_fn'])

    def resolve_sv(self):
        """Perform variant calling on the Contig object that was generated from the split reads in the target.

        All the assembled contigs are iterated through and processed.
        1. Realign to reference sequence.
        2. Process the realignment.
        3. Make variant call based on realignment.
        4. Check filters
        5. Annotate call.
        6. Output call.
        7. Perform discordant read clustering.

        Args:
            None

        Returns:
            None
        """

        iter = 1  # Label the contigs with an index.
        contigs = self.variation.kmers['clusters']  # Get the assembled Contig objects.
        utils.log(self.loggingName, 'info', 'Resolving structural variants from %d kmer clusters' % len(contigs))

        # Iterate through all the contigs, realign and make variant calls
        for contig in contigs:
            contigId = self.name + '_contig' + str(iter)
            utils.log(self.loggingName, 'info', 'Assessing contig %s, %s' % (contigId, contig.seq))
            contig.set_meta_information(contigId, self.params, self.values, self.paths['contigs'], self.variation.files['kmer_clusters'], self.variation)
            contig.query_ref(self.files['target_ref_fn'])
            contig.make_calls()
            if contig.svEventResult:
                contig.filter_calls()
                contig.annotate_calls()
                contig.output_calls(self.paths['output'], self.variation.files['sv_bam_sorted'])
                if contig.svEventResult:
                    self.variation.add_result(contig.svEventResult)
            else:
                utils.log(self.loggingName, 'info', '%s has no structural variant result.' % contigId)
            iter += 1
        self.variation.cluster_discreads(self.name, self.chrom)  # Cluster discordant reads.

    def complete_analysis(self):
        """
        """

        if len(self.variation.results) > 0 or len(self.variation.discReadFormatted) > 0:
            self.variation.write_results(self.paths['output'], self.name)
        else:
            shutil.rmtree(self.paths['output'])

    def get_target_intervals(self):
        """Return the list of tuples defining intervals for this target
        """

        return self.params.targets[self.name]

    def get_sv_reads(self, type):
        """ """

        return self.variation.get_sv_reads(type)

    def clear_sv_reads(self, type):
        """
        """

        self.variation.clear_sv_reads(type)

    def clear_cleaned_reads(self):
        """ """

        self.variation.clear_cleaned_reads()

    def has_results(self):
        """ """

        if len(self.variation.results) > 0 or len(self.variation.discReadFormatted) > 0:
            return True
        else:
            return False

    def get_results(self):
        """ """

        return self.variation.results

    def get_formatted_output(self):
        """ """

        return self.variation.get_formatted_output()
