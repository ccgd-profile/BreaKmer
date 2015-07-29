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


def load_kmers(fns, kmers):
    """Iterate through the kmer flat files and store them in the kmers dictionary.
    Store the kmer sequence string as the key and the count of the number of reads
    containing it as the value.

    Args:
        fns (str):      Filename of the kmer flat files.
        kmers (dict):   Dictionary of the kmer, count values,
    Returns:
        None
    Raises:
        None
    """

    if not fns:
        return kmers
    fns = fns.split(',')
    for fn in fns:  # Iterate through all the jellyfish kmer files and store the kmer as key and count as value.
        for line in open(fn, 'rU'):
            line = line.strip()
            mer, count = line.split()
            if mer not in kmers:
                kmers[mer] = 0
            kmers[mer] += int(count)


class Variation:
    """This class handles the storage and interaction of all the variant reads that could
    be contributing to the support of a structural variant.


    """

    def __init__(self, params):
        self.params = params
        self.var_reads = {}
        self.sv_reads = None
        self.cleaned_read_recs = None
        self.kmer_clusters = []
        self.kmers = {}
        self.results = []
        self.files = {}
        self.svs = {}
        self.discReadClusters = {}
        self.discReadFormatted = []
        self.loggingName = 'breakmer.processor.target'

    def setup_cleaned_reads(self, sampleType):
        """Initiate the cleaned_read_recs dictionary for sample or normal data.

        Args:
            sampleType (str): String indicating the sample type - sv / normal
        Returns:
            None
        Raises:
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
        """
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

    def get_sv_reads(self, type):
        """
        """
        return self.var_reads[type].sv

    def add_result(self, result):
        """
        """
        self.results.append(result)

    def set_var_reads(self, sampleType, bamFile, chrom, start, end, regionBuffer):
        """

        Args:
            sampleType ():
            bamFile ():
            chrom ():
            start ():
            end ():
            regionBuffer ():
        Returns:
            None
        Raises:
            None
        """

        # Get VariantReadTracker object from bam_handler module.
        self.var_reads[sampleType] = bam_handler.get_variant_reads(bamFile, chrom, start - regionBuffer, end - regionBuffer, self.params.get_param('insertsize_thresh'))
        # Iterate through reads that are not perfectly aligned and store necessary information for downstream analysis.
        self.var_reads[sampleType].check_clippings(self.params.get_kmer_size(), start, end)

        svBam = None
        if sampleType == 'sv':
            svBam = pysam.Samfile(self.files['sv_bam'], 'wb', template=pysam.Samfile(bamFile, 'rb'))
        readsFq = open(self.files['%s_fq' % sampleType], 'w')
        scFa = open(self.files['%s_sc_unmapped_fa' % sampleType], 'w')
        # Write all the stored sequences into files.
        self.var_reads[sampleType].write_seqs(scFa, readsFq, svBam, self.params.get_kmer_size())
        readsFq.close()
        scFa.close()

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
            sampleType (str):   The type of input data - sv / normal
            dataPath (str):     The path to the data files for this target.
            name (str):         The target name.
        Returns:
            None
        Raises:
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
        remove any 'noise' that bogs down the assembly process or analysis. The
        cleaned reads output from cutadapt are then re-processed to determine
        if the soft-clipped sequences were trimmed off or not to further filter
        out reads. The soft-clipped sequences that remain are stored and a new
        fastq file is written.

        Args:
            dataPath (str): The path to the data files for this target.
            name (str):     The target name.
            type (str):     A string indicating a tumor ('sv') or normal ('norm') sample being processed.
        Return:
            check (boolean): A boolean to indicate whether the are any reads left after
                             cleaning is complete.
        """

        cutadapt = self.params.get_param('cutadapt')  # Cutadapt binary
        cutadaptConfigFn = self.params.get_param('cutadapt_config_file')
        utils.log(self.loggingName, 'info', 'Cleaning reads using %s with configuration file %s' % (cutadapt, cutadaptConfigFn))
        self.files['%s_cleaned_fq' % sampleType] = os.path.join(dataPath, name + '_%s_reads_cleaned.fastq' % sampleType)
        utils.log(self.loggingName, 'info', 'Writing clean reads to %s' % self.files['%s_cleaned_fq' % sampleType])
        output, errors = utils.run_cutadapt(cutadapt, cutadaptConfigFn, self.files['%s_fq' % sampleType], self.files['%s_cleaned_fq' % sampleType], self.loggingName)

        self.setup_cleaned_reads(sampleType)
        self.files['%s_cleaned_fq' % sampleType], self.cleaned_read_recs[sampleType], self.read_len = utils.get_fastq_reads(self.files['%s_cleaned_fq' % sampleType], self.get_sv_reads(sampleType))
        self.clear_sv_reads(sampleType)
        check = self.continue_analysis_check(sampleType)
        utils.log(self.loggingName, 'info', 'Clean reads exist %s' % check)
        return check

    def set_reference_kmers(self, targetRefFns):
        """Set the reference sequence kmers"""
        self.kmers['ref'] = {}
        for i in range(len(targetRefFns)):
            utils.log(self.loggingName, 'info', 'Indexing kmers for reference sequence %s' % targetRefFns[i])
            self.get_kmers(targetRefFns[i], self.kmers['ref'])

    def set_sample_kmers(self):
        """Set the sample kmers
        """

        utils.log(self.loggingName, 'info', 'Indexing kmers for sample sequence %s' % self.files['sv_cleaned_fq'])
        self.kmers['case'] = {}
        self.kmers['case_sc'] = {}
        self.get_kmers(self.files['sv_cleaned_fq'], self.kmers['case'])
        self.get_kmers(self.files['sv_sc_unmapped_fa'], self.kmers['case_sc'])

    def get_kmers(self, seqFn, kmerDict):
        """Generic function to run jellyfish on a set of sequences
        """

        jellyfish = self.params.get_param('jellyfish')
        kmer_size = self.params.get_kmer_size()
        # Load the kmers into the kmer dictionary based on keyStr value.
        load_kmers(utils.run_jellyfish(seqFn, jellyfish, kmer_size), kmerDict)

    def compare_kmers(self, kmerPath, name, readLen, targetRefFns):
        """
        """

        # Set the reference sequence kmers.
        self.set_reference_kmers(targetRefFns)

        # Set sample kmers.
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
        self.kmers['case_only'] = {}
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

        self.kmers['clusters'] = assembly.init_assembly(self.kmers['case_only'], self.cleaned_read_recs['sv'], self.params.get_kmer_size(), self.params.get_sr_thresh('min'), readLen)
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
        headerStr = '\t'.join(['Target_name', 'sv_type', 'left_breakpoint_estimate', 'right_breakpoint_estimate', 'strands', 'discordant_readpair_count'])
        for key in self.discReadClusters:
            readCount = self.discReadClusters[key]['readCount']
            if readCount < self.params.get_param('discread_only_thresh'):
                continue
            k1, k2, k3, c1, c2 = key.split('|')
            svType = 'inter-chromosomal'
            lChrom = 'chr' + targetChrom.replace('chr', '')
            if k1 == 'inter':
                rChrom = 'chr' + k2.replace('chr', '')
            elif k1 == 'intra':
                svType = 'intra-chromosomal_' + k2
                rChrom = lChrom
            lStrand, rStrand = k3.split(':')
            lBrkpt = self.discReadClusters[key]['leftBrkpt']
            rBrkpt = self.discReadClusters[key]['rightBrkpt']
            outStr = '\t'.join([targetName, svType, lChrom + ':' + str(lBrkpt), rChrom + ':' + str(rBrkpt), lStrand + ',' + rStrand, str(readCount)])
            self.discReadFormatted.append((headerStr, outStr))


class TargetManager:
    """TargetManager class handles all the high level information relating to a target.
    The analysis is peformed at the target level, so this class contains all the information
    necessary to perform an independent analysis.

    Attributes:
        params (ParamManager):      Parameters for breakmer analysis.
        loggingName (str):          Module name for logging file purposes.
        name (str):                 Target name specified in the input bed file.
        chrom (str):                Chromosome ID as specified in the input bed file.
        start (int):                Genomic position for the target region (minimum value
                                    among all intervals).
        end (int):                  Genomic position for the target region (maximum value among
                                    all intervals).
        paths (dict):               Contains the analysis paths for this target.
        files (dict):               Dicionary containing paths to file names needed for analysis.
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
        self.read_len = params.get_param('readLen')
        # self.repeat_mask = None
        self.variation = Variation(params)
        self.regionBuffer = 200
        self.setup()

    def setup(self):
        """Setup the TargetManager object with the input params.

        Args:
            None
        Returns:
            None
        """

        intervals = self.params.targets[self.name]
        for values in intervals:
            chrom, start, end = values[0], int(values[1]), int(values[2])
            if not self.chrom:
                self.chrom = chrom
            if not self.start:
                self.start = start
            if not self.end:
                self.end = end
            if start < self.start:
                self.start = start
            if end > self.end:
                self.end = end

        # Create the proper paths for the target analysis.
        self.add_path('ref_data', os.path.join(self.params.paths['ref_data'], self.name))
        if self.params.fncCmd == 'run':
            self.add_path('base', os.path.join(self.params.paths['targets'], self.name))
            self.add_path('data', os.path.join(self.paths['base'], 'data'))
            self.add_path('contigs', os.path.join(self.paths['base'], 'contigs'))
            self.add_path('kmers', os.path.join(self.paths['base'], 'kmers'))
            self.add_path('output', os.path.join(self.params.paths['output'], self.name))
        self.files['target_ref_fn'] = [os.path.join(self.paths['ref_data'], self.name + '_forward_refseq.fa'), os.path.join(self.paths['ref_data'], self.name + '_reverse_refseq.fa')]
        ref_fa_marker_f = open(os.path.join(self.paths['ref_data'], '.reference_fasta'), 'w')
        ref_fa_marker_f.write(self.params.opts['reference_fasta'])
        ref_fa_marker_f.close()
        self.files['ref_kmer_dump_fn'] = [os.path.join(self.paths['ref_data'], self.name + '_forward_refseq.fa_dump'), os.path.join(self.paths['ref_data'], self.name + '_reverse_refseq.fa_dump')]

    def add_path(self, key, path):
        """Utility function to create all the output directories.

        Args:
            key (str):  String value to store the file path value.
            path (str): File path value.
        Returns:
            None
        Raises:
            None
        """

        utils.log(self.loggingName, 'info', 'Creating %s %s path (%s)' % (self.name, key, path))
        self.paths[key] = path
        if not os.path.exists(self.paths[key]):
            os.makedirs(self.paths[key])

    def set_ref_data(self):
        """
        """

        # Write reference fasta file if needed.
        for i in range(len(self.files['target_ref_fn'])):
            fn = self.files['target_ref_fn'][i]
            direction = "forward"
            if fn.find("forward") == -1:
                direction = "reverse"
            utils.log(self.loggingName, 'info', 'Extracting refseq sequence and writing %s' % fn)
            utils.extract_refseq_fa(self.get_values(), self.paths['ref_data'], self.params.opts['reference_fasta'], direction, fn)

        blastn = self.params.get_param('blast')
        if blastn is not None:
            # Check if blast db files are available for each target.
            if not os.path.isfile(self.files['target_ref_fn'][0] + '.nin'):
                # Use makedb
                makedb = os.path.join(os.path.split(blastn)[0], 'makeblastdb')
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

        Args:
            None
        Returns:
            check (boolean):    Variable to determine if the analysis should continue. It is
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

    def clean_reads(self, sampleType):
        """Wrapper for Variation clean_reads function.
        Args:
            type (str): A string indicating a tumor ('sv') or normal ('norm') sample being processed.
        Return:
            check (boolean): A boolean to indicate whether the are any reads left after
                             cleaning is complete.
        """
        return self.variation.clean_reads(self.paths['data'], self.name, sampleType)

    def extract_bam_reads(self, sampleType):
        """
        Args:
            sampleType (str): Indicates a tumor ('sv') or normal ('norm') sample being processed.
        Return:
            None
        """

        # Create the file paths for the files that will be created from the read extraction.
        self.variation.setup_read_extraction_files(sampleType, self.paths['data'], self.name)
        bamType = 'sample'
        if sampleType == 'norm':
            bamType = 'normal'
        bamFile = self.params.opts['%s_bam_file' % bamType]
        utils.log(self.loggingName, 'info', 'Extracting bam reads from %s to %s' % (bamFile, self.variation.files['%s_fq' % sampleType]))
        self.variation.set_var_reads(sampleType, bamFile, self.chrom, self.start, self.end, self.regionBuffer)

    def compare_kmers(self):
        """Obtain the sample only kmers and initiate assembly of reads with these kmers.
        """

        self.variation.compare_kmers(self.paths['kmers'], self.name, self.read_len, self.files['target_ref_fn'])

    def resolve_sv(self):
        """
        """

        iter = 1
        contigs = self.variation.kmers['clusters']
        utils.log(self.loggingName, 'info', 'Resolving structural variants from %d kmer clusters' % len(contigs))
        for contig in contigs:
            contigId = self.name + '_contig' + str(iter)
            utils.log(self.loggingName, 'info', 'Assessing contig %s, %s' % (contigId, contig.seq))
            contig.set_meta_information(contigId, self.params, self.get_values(), self.paths['contigs'], self.variation.files['kmer_clusters'], self.variation)
            contig.query_ref(self.files['target_ref_fn'])
            contig.make_calls()
            if contig.svEventResult:
                contig.filter_calls()
                contig.annotate_calls()
                contig.output_calls(self.paths['output'], self.variation.files['sv_bam_sorted'])
                self.add_result(contig.svEventResult)
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

    def get_values(self):
        """Return the defined features of this target
        """

        return(self.chrom, self.start, self.end, self.name, self.get_target_intervals(), self.regionBuffer)

    def get_sv_reads(self, type):
        """ """

        return self.variation.get_sv_reads(type)

    def clear_sv_reads(self, type):
        """ """

        self.variation.clear_sv_reads(type)

    def clear_cleaned_reads(self):
        """ """

        self.variation.clear_cleaned_reads()

    def add_result(self, result):
        """ """

        if result:
            self.variation.add_result(result)

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
        """
        """

        return self.variation.get_formatted_output()
