#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
import logging
import shutil
import pysam
import breakmer.assembly.olc as olcAssembly
import breakmer.assembly.utils as assemblyUtils
import breakmer.realignment.realigner as realigner
import breakmer.caller.sv_caller as sv_caller
import breakmer.utils as utils
import breakmer.annotation.sv_annotation as annotator
import breakmer.plotting.sv_viz as svplotter

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


def get_read_kmers(new_seq, kmerLen, kmer_seqs, order='for'):
    """Return new sample kmers from the existing contig sequence that can help extend
    the contig sequence.
    All the k-length mers are determined from the new_seq. These kmer sequences are put
    into a set and intersected with the kmer sequences in the sample, ordered according
    to the position of the kmer in the new_seq string and returned.
    Args:
        new_seq: String of the contig sequence to create kmers from.
        kmerLen: Integer of the kmer length
        kmer_seqs: The set of kmer sequences from the pool of extracted reads.
        order: String for the direction to order the news set of kmer sequences. A None
               value indicates no ordering.

    Return:
        kmers: List of tuples containing:
               1. String kmer seq
               2. Integer kmer position
               3. Boolean if kmer seq is in the first half of the sequence
               4. Integer of position distance to middle of sequence
               5. String of how to order tuples in the list
    """
    m = len(new_seq) / 2
    kmers = map(lambda x: (new_seq[x:x + kmerLen], x, int(x < m), abs(x - m), order), range(0, (len(new_seq) - kmerLen)))
    ks = set(map(lambda x: x[0], kmers))
    ss = ks & kmer_seqs
    kmers = filter(lambda x: x[0] in ss, kmers)
    if order == 'rev':
        kmers.reverse()
    elif order == 'mid':
        kmers = sorted(kmers, key=lambda x: (x[2], x[3]))
    else:
        kmers = list(ss)
    return kmers


class AssemblyRead:
    """Wrapper class for a sequence read used in a contig assembly. This will
    track meta information about the sequence read.
    Attributes:
        read:           fq_read object
        redundant:      Boolean to indicate whether the read is duplicated.
        alignChecked:   Boolean to indicate if the read has been checked against
                        the contig sequence.
        aligned:        Boolean to indicate if the read aligned to the contig sequence.
    """
    def __init__(self, read, redundant, checked, aligned):
        self.read = read
        self.redundant = redundant
        self.alignChecked = checked
        self.aligned = aligned


class ReadBatch:
    """A class to track the reads that are being considered for building a contig
    sequence.
    Attributes:
        delete:     Set of fq_read objects to remove from further analysis.
        alt:        List of tuples containing (fq_read object, integer of nreads with the same sequence)
        reads:      List of AssemblyRead objects containing the reads used to build a contig.
        mer_pos_d:  Dictionary containing kmer position information. DEPRECATED
    """
    def __init__(self, read, mer_pos):
        self.delete = set()
        self.alt = []
        self.reads = [AssemblyRead(read, False, True, True)]
        # self.mer_pos_d = {mer_pos: [0]} DEPRECATED

    def check_kmer_read(self, kmer_read_align_pos, read):
        """Adds AssemblyRead to reads list. Note that the check for add_to_pos_d is deprecated.
        Args:
            kmer_read_align_pos: Integer of the position the kmer sequence found
                                 in the read sequence.
            read:                fq_read object.
        Return: None
        """
        check = True
        redund_read = False
        add_read = True
        add_to_pos_d = False

        """
        # Deprecated code.
        if add_read :
            if add_to_pos_d :
                if pos not in self.mer_pos_d :
                    self.mer_pos_d[pos] = []
                self.mer_pos_d[pos].append(len(self.reads))
            self.reads.append(AssemblyRead(read, False, check, False))
        return check
        """
        self.reads.append(AssemblyRead(read, False, check, False))

    def set_last_read_aligned(self):
        """Sets the last read added to the reads list as aligned."""
        self.reads[-1].aligned = True

    def clean(self, fq_reads, contigBuffer, last_keep_read):
        """Remove all data from data structures.
        Iterate through reads in delete set and delete them from the fq dictionary.
        Check if the delete reads are in the contigBuffer contig dictionary. If the
        contig associated with the read is not setup then delete the read from the dictionary.
        Args:
            fq_reads: Dictionary containing the extracted reads.
            contigBuffer: ContigBuffer object.
            last_keep_read: fq_read object kept for further use.
        Return: None
        """

        map(fq_reads.__delitem__, map(lambda x: x.seq, list(self.delete)))
        for read_id in filter(lambda x: x in contigBuffer.contigs, list(self.delete)):
            if not contigBuffer.contigs[read_id].setup:
                del contigBuffer.contigs[read_id]
        self.delete = set()
        self.alt = []
        self.reads = [last_keep_read]
        self.mer_pos_d = {}


class ContigCounts:
    """A class to track the number of read sequences that support a consensus sequence.

    Initially set counts for the first read in the contig.
    Attributes:
        indel_only: List of integers, providing the count for number of indel only reads support
                    the given position in the consensus sequence.
        others: List of integers, providing the count of non indel only reads that are assembled
                at a given position of the consensus sequence.
    """

    def __init__(self, read, nreads):
        self.indel_only = [0] * len(read.seq)
        self.others = [0] * len(read.seq)
        self.set_counts(0, len(read.seq), nreads, read.indel_only)

    def get_counts(self, p1, p2, sv_type):
        """Return the counts for a range of positions in the consensus sequence.
        If the positions are the same, then return the counts for the single position.
        Args:
            p1: Integer indicating the first position.
            p2: Integer indicating the second position.
            sv_type: String indicating what kind of event the count is intended to support.
        Return:
            counts: List of integers for counts of reads assembled at the provided range.
        """
        counts = []
        if sv_type == 'indel' or sv_type == 'rearr':
            if p1 == p2:
                counts = self.indel_only[p1] + self.others[p1]
            else:
                counts = map(lambda (x, y): x + y, zip(self.indel_only[p1:p2], self.others[p1:p2]))
        else:
            if p1 == p2:
                counts = self.others[p1]
            else:
                counts = self.others[p1:p2]
        return counts

    def get_total_reads(self):
        """Return the total read count supporting a contig sequence."""
        return max(self.indel_only) + max(self.others)

    def set_superseq(self, read, nreads, start, end):
        """The read sequence is a super sequence to the current contig sequence.
        The count vectors need to be adjusted accordingly based on the read.
        Temporary count vectors are created for the read sequence and the number
        of reads with this sequence. The current count vectors are then added into
        the temporary vectors and then set as the new count vectors.
        Args:
            read: fq_read object.
            nreads: Integer for number of reads with read sequence.
            start: Integer for start alignment position of the current contig in
                   the read sequence.
            end: Integer for the end alignment position of the current contig in
                 the read sequence.
        Return: None
        """
        tmp_indel_only = [0] * len(read.seq)
        tmp_others = [0] * len(read.seq)
        if read.indel_only:
            tmp_indel_only = [nreads] * len(read.seq)
        else:
            tmp_others = [nreads] * len(read.seq)
        tmp_indel_only[start:end] = map(lambda (x, y): x + y, zip(tmp_indel_only[start:end], self.indel_only))
        tmp_others[start:end] = map(lambda (x, y): x + y, zip(tmp_others[start:end], self.others))
        self.indel_only = tmp_indel_only
        self.others = tmp_others

    def set_counts(self, start, end, nreads, indel_only):
        """Add the read count to the stored contig sequence count vectors.
        With paired end reads, there are reads that can contribute to indels only
        or to all types of variation. The counts are added according to how the
        read has been defined.
        Args:
            start: Integer of the start of the sequence to add count.
            end: Integer of the end of the sequence to add count.
            nreads: Integer of total number of reads that should be added.
            indel_only: Boolean to indicate if read should only support indels.
        Return: None
        """
        if indel_only:
            self.indel_only[start:end] = map(lambda x: x + nreads, self.indel_only[start:end])
        else:
            self.others[start:end] = map(lambda x: x + nreads, self.others[start:end])

    def extend_counts(self, extend_size, nreads, indel_only, direction):
        """Increase the size of the count vectors when the contig sequence is grown.
        If the direction is 'post', the count vectors must be increased at the end.
        If the direction is 'pre', the count vectors must be increased at the beginning.
        Args:
            extend_size: Integer for number of positions to increase the vectors.
            nreads: Integer for count to add to the count vectors.
            indel_only: Boolean to indicate if count should only be added to indel only vector.
            direction: String to indicate which side the count vector is extended.
        Return: None
        """
        fill_counts = [0] * extend_size
        ecounts = [nreads] * extend_size
        if indel_only:
            if direction == 'post':
                self.indel_only.extend(ecounts)
                self.others.extend(fill_counts)
            else:
                ecounts.extend(self.indel_only)
                self.indel_only = ecounts
                fill_counts.extend(self.others)
                self.others = fill_counts
        else:
            if direction == 'post':
                self.indel_only.extend(fill_counts)
                self.others.extend(ecounts)
            else:
                ecounts.extend(self.others)
                self.others = ecounts
                fill_counts.extend(self.indel_only)
                self.indel_only = fill_counts


class Builder:
    """A class to perform all the contig building functions and store temporary data structures.
    Attributes:
        read_batch:     ReadBatch object
        seq:            String of the consensus sequence.
        counts:         ContigCounts object to manage all the read counts supporting the consensus sequence.
        checked_kmers:  List of kmer sequences that had previously been checked while building the contig.
        kmerLen:        Integer of the kmer length.
        kmers:          List of kmer sequences that have contributed to building the contig.
        kmer_locs:      List of integers representing the positions of the kmers in the contig seq.
    """

    def __init__(self, kmerObj, readAlignValues):
        """
        Args:
            kmerObj:            A Kmer object containing seq and count information for a kmer.
            readAlignValues:    A dictionary containing information about where a kmer exists in a read.
        """
        self.read_batch = ReadBatch(readAlignValues['read'], readAlignValues['align_pos'])
        self.seq = readAlignValues['read'].seq
        self.counts = ContigCounts(readAlignValues['read'], readAlignValues['nreads'])
        self.checked_kmers = [kmerObj.seq]
        self.kmerLen = kmerObj.kmerLen
        self.kmers = []
        self.kmer_locs = []

    def check_read(self, kmerObj, readAlignValues, alignType):
        """Determine if the read should be added to the assembly or not.
        If the read aligns to the contig, set the fq_read status to used and indicate
        the AssemblyRead has been aligned. If the kmer is in more than 1 read and the
        current read has not been used in any other contigs then store for later
        analysis to build another contig possibly. Otherwise, discard the read for
        further analysis.
        Args:
            kmerObj:       Kmer object containing kmer seq specific values.
            readAlignValues: Dictionary containing:
                         - 'read': fq_read object that contains kmer sequence.
                         - 'align_pos': Integer position of kmer in read sequence
                         - 'nreads': Integer of number of reads with the same sequence.
            type: String indicating the state of this function.
        Return:
            hit: String value 'remove' or ''.
        """
        hit = ''
        self.read_batch.check_kmer_read(readAlignValues['align_pos'], readAlignValues['read'])
        if self.check_align(kmerObj, readAlignValues, alignType):
            hit = 'remove'
            readAlignValues['read'].used = True
            self.read_batch.set_last_read_aligned()
        elif kmerObj.counts > 2 and not readAlignValues['read'].used:
            self.read_batch.alt.append((readAlignValues['read'], readAlignValues['nreads']))
        else:
            self.read_batch.delete.add(readAlignValues['read'])
        return hit

    def check_align(self, kmerObj, readAlignValues, alignType='setup'):
        """Check the alignment of the read sequence to the contig sequence.
        The read sequence must match at least 25% of the shortest sequence between
        the contig and the read and an identity at least 90%. If there is clear
        alignment between the read sequence and the contig sequence, then the
        assembly consensus sequence is appropriately changed.
        Args:
            kmerObj: Dictionary containing:
                         - 'seq': String kmer sequence value.
                         - 'counts': Integer of reads containing kmer sequence.
                         - 'kmer_set': Set with all kmer sequences.
                         - 'len': Integer for kmer length.
            readAlignValues: Dictionary containing:
                         - 'read': fq_read object that contains kmer sequence.
                         - 'align_pos': Integer position of kmer in read sequence
                         - 'nreads': Integer of number of reads with the same sequence.
            type: String indicating the state of this function.
        Return:
            match: Boolean indicating if the read aligns sufficiently with the
                   contig sequence and will be added.
        """
        match = False
        queryRead = readAlignValues['read']

        minScore = float(min(len(self.seq), len(queryRead.seq))) / 4.0
        alignManager = olcAssembly.AlignManager(self.seq, queryRead.seq, minScore, 0.90)

        if alignManager.check_align_thresholds():
            # Alignment fails thresholds.
            return False
        if alignManager.same_seqs():
            # Read and contigs sequences are the same.
            return True
        if alignManager.same_max_scores():
            # Alignments both ways had equal scores.
            match = True
            if alignManager.read_is_superseq():
                # Check if the read sequence fully contains the contig sequence.
                self.set_superseq(queryRead, readAlignValues['nreads'], alignManager.get_alignment_values(0, 'i'), alignManager.get_alignment_values(0, 'prei'))
                if alignType == 'grow':
                    # Contig sequence has changed, set the kmers.
                    self.set_kmers(kmerObj.kmerSeqSet)
            # Check if the contig sequence full contains the read sequence.
            elif alignManager.read_is_subseq():
                self.add_subseq(alignManager.get_alignment_values(1, 'i'), alignManager.get_alignment_values(1, 'prei'), readAlignValues['nreads'], queryRead.indel_only)
            # There appears to be overlap, figure out how to assemble.
            else:
                match = False
                indx1 = alignManager.get_kmer_align_index(0, kmerObj.seq)
                indx2 = alignManager.get_kmer_align_index(1, kmerObj.seq)
                if indx1[0] > -1 and indx1[1] > -1:
                    # Read overlaps off front of contig sequence.
                    if (indx2[0] == -1 and indx2[1] == -1) or (abs(indx2[0] - indx2[1]) > abs(indx1[0] - indx1[1])):
                        match = True
                        self.contig_overlap_read(alignManager.get_alignment(0), queryRead, readAlignValues['nreads'], kmerObj.kmerSeqSet, alignType)
                elif indx2[0] > -1 and indx2[1] > -1:
                    # Read overlaps off end of contig sequence.
                    if (indx1[0] == -1 and indx1[1] == -1) or (abs(indx2[0] - indx2[1]) < abs(indx1[0] - indx1[1])):
                        match = True
                        self.read_overlap_contig(alignManager.get_alignment(1), queryRead, readAlignValues['nreads'], kmerObj.kmerSeqSet, alignType)
        # Read sequence overlaps off the front of the contig sequence.
        elif alignManager.better_align():
            match = True
            self.contig_overlap_read(alignManager.get_alignment(0), queryRead, readAlignValues['nreads'], kmerObj.kmerSeqSet, alignType)
        # Read sequence overlaps off the end of the contig sequence.
        else:
            match = True
            self.read_overlap_contig(alignManager.get_alignment(1), queryRead, readAlignValues['nreads'], kmerObj.kmerSeqSet, alignType)
        return match

    def set_superseq(self, read, nreads, start, end):
        """The read sequence contains the current contig sequence.
        Args:
            read: fq_read object.
            nreads: Integer for the number of reads with the same sequence as the read passed in.
            start: Integer for the start position the contig sequence aligns to the read sequence.
            end: Integer for the end position the contig sequence aligns to the read sequence.
        Return: None
        """
        self.seq = read.seq
        self.counts.set_superseq(read, nreads, start, end)

    def add_subseq(self, start, end, nreads, indel_only):
        """The read checked against the contig was found to be a subsequence of the
        contig. The nreads with the checked sequence are added to the count vectors.
        Args:
            start: Integer for start of sequence to add counts.
            end: Integer for end of sequence to add counts.
            nreads: Integer of number of reads to add.
            indel_only: Boolean to indicate whether to add to indel only count vector.
        Return: None
        """
        self.counts.set_counts(start, end, nreads, indel_only)

    def add_postseq(self, post_seq, start, end, nreads, indel_only):
        """Sequence is appended to the end of the current contig sequence. The
        read support vectors are appropriately incremented.
        Args:
            post_seq: String of sequence to add to the end of the assembled contig.
            start: Integer for start position of contig to add counts.
            end: Integer for end position of the contig to add counts.
            nreads: Integer for number of reads to add to count vectors.
            indel_only: Boolean to indicate whether read only supports indel events.
        Return: None
        """
        self.seq += post_seq
        self.counts.set_counts(start, end, nreads, indel_only)
        self.counts.extend_counts(len(post_seq), nreads, indel_only, 'post')

    def add_preseq(self, pre_seq, start, end, nreads, indel_only):
        """Sequence is append to the front of the current contig sequence. The read
        support vectors are appropriately changed.
        Args:
            pre_seq: String of sequence to add to the front of the assembly contig.
            start: Integer for start position of contig to add counts.
            end: Integer for end position of contig to add counts.
            nreads: Integer for number of reads to add to count vectors.
            indel_only: Boolean to indicate whether read only supports indel events.
        Return: None
        """
        self.seq = pre_seq + self.seq
        self.counts.set_counts(start, end, nreads, indel_only)
        self.counts.extend_counts(len(pre_seq), nreads, indel_only, 'pre')

    def finalize_reads(self, contig_reads, fq_recs, contigBuffer):
        """Sort out the reads to keep for reporting and remove the others.
        Aligned and non-redundant reads are removed from the contig read set. The
        variables in read_batch are cleared.
        Args:
            contig_reads: Set of fq_read objects.
            fq_recs: Dictionary of fq_read objects.
            contigBuffer: ContigBuffer class object.
        Return:
            contig_reads: Set of fq_reads objects
        """
        rm_reads = map(lambda y: y.read, filter(lambda x: x.redundant, self.read_batch.reads))
        keep_reads = filter(lambda x: x.aligned and not x.redundant, self.read_batch.reads)
        add_reads = map(lambda y: y.read, keep_reads)
        # Merge add_reads into contig_reads
        contig_reads = contig_reads | set(add_reads)
        # Remove rm_reads
        contig_reads = contig_reads - set(rm_reads)
        self.read_batch.clean(fq_recs, contigBuffer, keep_reads[-1])
        return contig_reads

    def contig_overlap_read(self, alignment, query_read, nreads, kmer_seqs, assemblyType):
        """Assembled consensus and read sequences, where the consensus end overlaps
        with the read sequence beginning.
        Args:
            alignment: olc.Align object
            query_read: fq_read object
            nreads: Integer for number of reads to add to count vectors.
            kmer_seqs: Set of kmer sequence values.
            type: String for source of call to function.
        Return: None
        """
        if alignment.prej == len(self.seq) and alignment.j == 0:
            self.set_superseq(query_read, nreads, alignment.i, alignment.prei)
            if assemblyType == 'grow':
                self.set_kmers(kmer_seqs)
        else:
            post_seq = query_read.seq[alignment.prei:]
            nseq = self.seq[(len(self.seq) - (self.kmerLen - 1)):] + post_seq
            self.add_postseq(post_seq, alignment.j, alignment.prej, nreads, query_read.indel_only)
            if assemblyType == 'grow':
                nkmers = get_read_kmers(nseq, self.kmerLen, kmer_seqs, 'for')
                self.kmers.extend(nkmers)

    def read_overlap_contig(self, alignment, query_read, nreads, kmer_seqs, type):
        """Assemble consensus and read sequences togheter, where the consensus
        beginning overlaps the read sequence end.
        Args:
            alignment: olc.Align object
            query_read: fq_read object
            nreads: Integer for number of reads to add to count vectors.
            kmer_seqs: Set of kmer sequence values.
            type: String for source of call to function.
        Return: None
        """

        if alignment.prej == len(query_read.seq) and alignment.j == 0:
            self.add_subseq(alignment.i, alignment.prei, nreads, query_read.indel_only)
        else:
            pre_seq = query_read.seq[0:alignment.j]
            nseq = pre_seq + self.seq[0:(self.kmerLen - 1)]
            self.add_preseq(query_read.seq[0:alignment.j], alignment.i, alignment.prei, nreads, query_read.indel_only)
            if type == 'grow':
                nkmers = get_read_kmers(nseq, self.kmerLen, kmer_seqs, 'rev')
                self.kmers.extend(nkmers)

    def check_alternate_reads(self, kmerTracker, contigBuffer, contigKmers):
        """Iterate through the buffered reads that were not aligned to the contig
        and determine if a new contig should be created.
        Args:
            kmerTracker: KmerTracker object contains all the kmer sequence values.
            contigBuffer: ContigBuffer object
            contigKmers: List of kmer sequence used in the contig assembly.
        """
        newContigs = []
        kmerSet = set()
        for read, nreads in self.read_batch.alt:
            altKmers = get_read_kmers(read.seq, self.kmerLen, kmerTracker.kmerSeqs, '')
            newKmers = set(altKmers) - set(contigKmers) - contigBuffer.used_kmers - kmerSet
            if len(newKmers) > 0:
                for kmerSeq in list(newKmers):
                    readCount = kmerTracker.get_count(kmerSeq)
                    if readCount > 1:
                        kmerPos = read.seq.find(kmerSeq)
                        kmerObj = assemblyUtils.Kmer(kmerSeq, kmerTracker.get_count(kmerSeq), kmerTracker.kmerSeqs, self.kmerLen)
                        read_align_values = {'read': read,
                                             'align_pos': kmerPos,
                                             'nreads': nreads}
                        newContigs.append((read, Contig(kmerObj, read_align_values)))
                        kmerSet = kmerSet | newKmers
                        break
        return newContigs

    def set_kmers(self, kmer_seqs):
        """Wrapper function to get_read_kmers function to parse a sequence string
        and generate all relevant kmers from the sequence.
        Args:
            kmer_seqs: Set of kmer sequences.
        Return: None
        """

        self.kmers = get_read_kmers(str(self.seq), self.kmerLen, kmer_seqs, 'mid')

    def set_kmer_locs(self):
        """Add the start alignment positions of each kmer sequence in the kmers list
        to the kmer_locs list.
        Args: None
        Return: None
        """

        self.kmer_locs = [0] * len(self.seq)
        for kmer in self.kmers:
            kmerPos = self.seq.find(kmer[0])
            self.kmer_locs[kmerPos:(kmerPos + self.kmerLen)] = map(lambda x: x + 1, self.kmer_locs[kmerPos:(kmerPos + self.kmerLen)])

    def refresh_kmers(self):
        """Return a list of kmer_sequences that have not been checked already.
        Args: None
        Return:
            List of kmer sequences.
        """
        print self.kmers
        return filter(lambda x: x[0] not in set(self.checked_kmers), self.kmers)

    def get_seq(self):
        """Return the final consensus sequence."""
        return self.seq

    def get_kmers(self):
        """Return the final kmer sequence list."""
        return self.kmers

    def get_kmer_locs(self):
        """Return the final kmer locations list."""
        return self.kmer_locs

    def get_total_reads(self):
        """Return total number of reads supporting contig."""
        return self.counts.get_total_reads()


class Meta:
    """A class to track the contig information for downstream calling and writing
    to file.
    Attributes:
        params: Param object with all BreaKmer parameters.
        path: String of path to write all files.
        id: String for contig ID.
        target_region: Tuple containing target information:
                       1. String chromosome ID
                       2. Integer of target region start position.
                       3. Integer of target region end position.
                       4. String target region name.
                       5. List of target intervals tuples.
        fq_fn: String of the fastq file containing sequence reads used to build contig.
        fa_fn: String of the fasta file containing the contig sequence.
    """

    def __init__(self):
        self.loggingName = 'breakmer.assembly.contig'
        self.params = None
        self.path = None
        self.id = None
        self.chr = None
        self.start = None
        self.end = None
        self.targetName = None
        self.regionBuffer = 0
        self.fq_fn = None
        self.fa_fn = None
        self.readVariation = None

    def set_values(self, contigId, params, queryRegionValues, contigPath, readVariation):
        """Sets the contig values after contig has been compeleted and ready for
        realignment.
        Args:
            contigId:           String containing contid ID.
            params:             Param object.
            queryRegionValues:  Tuple containing the target region information
            contigPath:         String of the path to the contig directory to store files.
        Return: None
        """
        self.params = params
        self.id = contigId
        self.readVariation = readVariation
        self.chr = queryRegionValues[0]
        self.start = int(queryRegionValues[1])
        self.end = int(queryRegionValues[2])
        self.targetName = queryRegionValues[3]
        self.regionBuffer = queryRegionValues[5]
        self.path = os.path.join(contigPath, self.id)
        logger = logging.getLogger('breakmer.assembly.contig')
        utils.log(self.loggingName, 'info', 'Setting up contig path %s' % self.path)

        if not os.path.exists(self.path):
            os.makedirs(self.path)
        self.fq_fn = os.path.join(contigPath, self.id, self.id + '.fq')
        self.fa_fn = os.path.join(contigPath, self.id, self.id + '.fa')

    def get_target_region_coordinates(self):
        """ """
        return (self.chr, self.start, self.end, self.targetName, self.regionBuffer)

    def write_files(self, cluster_fn, kmers, reads, seq):
        """Write cluster, read fastq, and contig fasta files.
        Args:
            cluster_fn: String of the file to write kmer clusters to.
            kmers: List of kmers used in the building of the contig.
            reads: List of reads used in the building of the contig.
            seq: String of contig sequence.
        Return: None
        """
        logger = logging.getLogger('breakmer.assembly.contig')
        cluster_f = open(cluster_fn, 'w')
        cluster_f.write(self.id + ' ' + str(len(kmers)) + '\n')
        cluster_f.write(','.join([x[0] for x in kmers]) + '\n')
        cluster_f.write(','.join([x.id for x in reads]) + '\n\n')
        cluster_f.close()
        assembly_fq = open(self.fq_fn, 'w')
        logger.info('Writing reads containing kmers to fastq %s' % self.fq_fn)
        for read in reads:
            assembly_fq.write(read.id + '\n' + read.seq + '\n+\n' + read.qual + '\n')
        assembly_fq.close()
        logger.info('Writing contig fasta file for blatting %s' % self.fa_fn)
        blat_f = open(self.fa_fn, 'w')
        blat_f.write('>' + self.id + '\n' + seq)
        blat_f.close()

    def write_result(self, svEventResult, outputPath):
        resultFn = os.path.join(self.path, self.id + "_svs.out")
        utils.log(self.loggingName, 'info', 'Writing %s result file %s' % (self.id, resultFn))
        resultFile = open(resultFn, 'w')

        # A string of output values for writing to file.
        headerStr, formattedResultValuesStr = svEventResult.get_formatted_output_values()
        resultFile.write(headerStr + '\n' + formattedResultValuesStr + '\n')
        resultFile.close()
        shutil.copyfile(resultFn, os.path.join(outputPath, self.id + "_svs.out"))

    def write_bam(self, outputPath, svBamReadsFn, reads):
        bamOutFn = os.path.join(outputPath, self.id + "_reads.bam")
        utils.log(self.loggingName, 'info', 'Writing contig reads bam file %s' % bamOutFn)
        bam_out_sorted_fn = os.path.join(outputPath, self.id + "_reads.sorted.bam")
        bamFile = pysam.Samfile(svBamReadsFn, 'rb')
        bam_out_f = pysam.Samfile(bamOutFn, 'wb', template=bamFile)
        for bam_read in bamFile.fetch():
            for read in reads:
                rid, idx = read.id.lstrip("@").split("/")
                ridx, indel_only_read = idx.split("_")
                if (bam_read.qname == rid) and ((ridx == '2' and bam_read.is_read2) or (ridx == '1' and bam_read.is_read1)):
                    bam_out_f.write(bam_read)
        bamFile.close()
        bam_out_f.close()
        utils.log(self.loggingName, 'info', 'Sorting bam file %s to %s' % (bamOutFn, bam_out_sorted_fn))
        pysam.sort(bamOutFn, bam_out_sorted_fn.replace('.bam', ''))
        utils.log(self.loggingName, 'info', 'Indexing bam file %s' % bam_out_sorted_fn)
        pysam.index(bam_out_sorted_fn)
        return bam_out_sorted_fn


class Contig:
    """Interface class to assemble a contig and store data all the relevant data
    for the assembly.
    Attributes:
        meta:       Meta class object to store all the interface related data.
        kmer_locs:  List of integers indicating the start alignment position of the kmers
                    in the contig sequence.
        setup:      Boolean to indicate whether a contig has gone through the setup process.
        build:      Builder class object that handles all the assembly functions.
        seq:        String for assembled sequence.
        kmers:      List of kmer sequences used to build contig.
        reads:      Set of read IDs that have been used to build contig.
        buffer:     Set of read IDs used in a batch of processing for building a contig. This is flushed.
    """

    def __init__(self, kmerObj, readAlignValues):
        self.meta = Meta()
        self.setup = False
        self.builder = Builder(kmerObj, readAlignValues)
        self.seq = None
        self.kmers = []
        self.kmer_locs = []
        self.reads = set()
        self.buffer = set([readAlignValues['read'].id])
        self.svEventResult = None
        self.realignment = None

    def check_read(self, kmerObj, readAlignValues, fncType='setup'):
        """Check if the read passed in can be added to the current contig.
        Wrapper function to Builder class check_read function.
        Args:
            kmerObj:         Instance of Kmer object with attributes for kmer sequence.
            readAlignValues: Dictionary containing:
                         - 'read': fq_read object that contains kmer sequence.
                         - 'align_pos': Integer position of kmer in read sequence
                         - 'nreads': Integer of number of reads with the same sequence.
            fncType: String indicating the state of this function.
        Return:
            String containing 'hit' or '' indicating that read matched contig seq
            or did not, respectively.
        """
        self.buffer.add(readAlignValues['read'].id)
        return self.builder.check_read(kmerObj, readAlignValues, fncType)

    def check_valid(self, read_count_thresh, read_len):
        """Determine if the finished contig sequence meets minimum requirements for
        length and read count support.
        Args:
            read_count_thresh: Integer for minimum reads that must support a contig.
            read_len: Integer for read length.
        Return:
            Boolean indicating whether it meets thresholds or not.
        """
        if self.get_total_read_support() < int(read_count_thresh) or len(self.seq) <= int(read_len):
            return True
        else:
            return False

    def finalize(self, fq_recs, kmerTracker, contigBuffer, source='setup'):
        """Finish an assembly and add the buffered contigs that were created from
        non-aligned reads to the contigBuffer.
        Args:
            fq_recs: Dicionary of fq_read objects.
            kmerTracker: KmerTracker object with all kmer sequence values.
            contigBuffer: ContigBuffer object.
            source: String for the source of function call.
        Return: None
        """
        if source == 'setup':
            self.set_kmers(kmerTracker.kmerSeqs)
        # Get alternate read kmers and see if any are different from contig kmers.
        new_contigs = self.builder.check_alternate_reads(kmerTracker, contigBuffer, self.kmers)
        for new_contig in new_contigs:
            contigBuffer.add_contig(new_contig[0], new_contig[1])
        self.reads = self.builder.finalize_reads(self.reads, fq_recs, contigBuffer)

    def set_kmers(self, kmer_seqs):
        """Wrapper function to Builder class set_kmers function.
        Args:
            kmer_seqs: Set of all kmer sequences.
        Return: None
        """
        self.setup = True
        self.builder.set_kmers(kmer_seqs)

    def set_kmer_locs(self):
        """Wrapper function to Builder class set_kmer_locs function.
        Args: None
        Return: None
        """
        self.builder.set_kmer_locs()

    def refresh_kmers(self):
        """A wrapper function to Builder class refresh kmers.
        Args: None
        Return:
            List of kmers that have not been previously checked.
        """
        return self.builder.refresh_kmers()

    def get_kmer_reads(self, kmer_values, read_items):
        """
        Args:
            kmer_values: Tuple containing the alignment information of a kmer sequence
                         in a read sequence.
                         1. String kmer sequence.
                         2. Integer kmer alignment position in sequence.
                         3. Boolean whether kmer align position is below midpoint of sequence.
                         4. Integer of absolute difference between align position and midpoint.
                         5. String of the order for tuples in a list.
            read_items: List of tuples for sequence reads:
                        1. String read sequence.
                        2. List of fq_read objects with read sequence.
        Return:
            reads: List of tuples containing:
                   1. read object,
                   2. start position of kmer match in read seq
                   3. Boolean that a match was found.
                   4. Length of the read sequence.
                   5. Number of reads with this sequence.
        """
        kmer, kmerPos, lessThanHalf, dist_half, order = kmer_values
        read_order = 'for'
        if order == 'mid':
            if lessThanHalf == 0:
                read_order = 'rev'
        elif order == 'for':
            read_order = 'rev'
        reads = assemblyUtils.find_reads(kmer, read_items, self.buffer, read_order)
        return reads

    def grow(self, fqRecs, kmerTracker, kmerLen, contigBuffer):
        """Iterates through new sample only kmers in a contig assembly and tries to
        add more relevant reads to extend the contig assembly sequence.
        For each 'new' kmer, assess the reads that have the kmer. When this function
        is complete, the contig is done assemblying.
        Args:
            fqRecs:         Dictionary of fq_read objects key = sequence, value = list of fq_reads
            kmerTracker:    KmerTracker object containing all the kmer sequences.
            kmerLen:        Integer of kmer size.
            contigBuffer:   ContigBuffer object.
        Return: None
        """
        logger = logging.getLogger('breakmer.assembly.contig')
        if not self.setup:
            self.set_kmers(kmerTracker.kmerSeqs)
        newKmers = self.refresh_kmers()
        while len(newKmers) > 0:
            iter = 0
            for kmer_lst in newKmers:
                print kmer_lst
                kmerSeq, kmerPos, lessThanHalf, dist_half, order = kmer_lst
                reads = self.get_kmer_reads(kmer_lst, fqRecs.items())
                contigBuffer.add_used_mer(kmerSeq)
                kmerObj = assemblyUtils.Kmer(kmerSeq, kmerTracker.get_count(kmerSeq), kmerTracker.kmerSeqs, kmerLen)
                for read_lst in reads:
                    read, kmerPos, bool, rlen, nreads = read_lst
                    contigBuffer.add_used_read(read.id)
                    readAlignValues = {'read': read,
                                       'align_pos': kmerPos,
                                       'nreads': nreads}
                    hit = self.check_read(kmerObj, readAlignValues, 'grow')
                    if hit == 'remove':
                        contigBuffer.remove_contig(read.id)
                self.finalize(fqRecs, kmerTracker, contigBuffer, 'grow')
                self.builder.checked_kmers.append(kmerSeq)
                iter += 1
            newKmers = self.refresh_kmers()
            logger.debug("%d kmers left to check" % len(newKmers))
        self.set_kmer_locs()
        self.set_final_values()
        logger.info('Contig done with contig seq %s. Supported by %d read(s).' % (self.seq, len(self.reads)))
        logger.info('Read IDs: %s' % (",".join([x.id for x in list(self.reads)])))

    def set_meta_information(self, contig_id, params, query_region_values, contig_path, kmer_cluster_fn, readVariation):
        """Sets the contig ID, params, target region values and contig path variables for later use.
        Output files are also written with the contig assembly information.
        Args:
            contig_id: String containing contid ID.
            params: Param object.
            query_region_values: Tuple containing the target region information
            contig_path: String of the path to the contig directory to store files.
            kmer_cluster_fn: String of the path to write the kmer clustering information to.
        Return: None
        """
        self.meta.set_values(contig_id, params, query_region_values, contig_path, readVariation)
        self.meta.write_files(kmer_cluster_fn, self.kmers, self.reads, self.seq)

    def set_final_values(self):
        """Set the seq, kmers, kmer_locs variables when the contig is done assemblying.
        Args: None
        Return: None
        """
        self.seq = self.builder.get_seq()
        self.kmers = self.builder.get_kmers()
        self.kmer_locs = self.builder.get_kmer_locs()

    def query_ref(self, targetRefFns):
        """
        Args:
        Return:
        """
        self.realignment = realigner.RealignManager(self.meta.params, targetRefFns)
        self.realignment.realign(self)

    def make_calls(self):
        """
        """
        contigCaller = sv_caller.ContigCaller(self.realignment, self, self.meta.params)
        self.svEventResult = contigCaller.call_svs()

    def filter_calls(self):
        """
        """
        if self.svEventResult is not None:
            svFilter = self.meta.params.filter
            # print 'contig.py check_filters'
            svFilter.check_filters(self.svEventResult)

    def annotate_calls(self):
        """ """
        if self.svEventResult and self.meta.params.get_param('gene_annotation_file') and self.meta.params.get_param('bedtools'):
            annotator.annotate_event(self.svEventResult, self.meta)

    def output_calls(self, outputPath, svReadsBamFn):
        """ """
        if self.svEventResult:
            self.meta.write_result(self.svEventResult, outputPath)
            readBamFn = self.meta.write_bam(outputPath, svReadsBamFn, self.reads)
            if self.meta.params.get_param('generate_image'):
                svplotter.generate_pileup_img(self.svEventResult, readBamFn, outputPath, self.get_id())

    def get_total_read_support(self):
        """Return the total read count supporting assembly."""
        return self.builder.get_total_reads()

    def get_contig_len(self):
        """Return length of contig sequence."""
        return len(self.seq)

    def get_kmer_locs(self):
        """Return the kmer locations in the contig sequence."""
        return self.kmer_locs

    def has_fa_fn(self):
        """Check if fasta file has been written for contig."""
        return self.meta.fa_fn

    def get_path(self):
        """Return file path to contig results"""
        return self.meta.path

    def get_id(self):
        """Return contig id"""
        return self.meta.id

    def get_target_name(self):
        """ """
        return self.meta.targetName

    def get_contig_count_tracker(self):
        """ """
        return self.builder.counts

    def get_disc_reads(self):
        """ """
        return self.meta.readVariation.get_disc_reads()

    def get_read_variation(self):
        """ """
        return self.meta.readVariation

    def get_var_reads(self, sampleType):
        """ """
        return self.meta.readVariation.get_var_reads(sampleType)

    def get_sample_bam_fn(self):
        """ """
        return self.meta.params.get_param('sample_bam_file')

    def get_target_region_coordinates(self):
        """ """
        return self.meta.get_target_region_coordinates()

    def get_chr(self):
        """ """
        return self.meta.chr

    def get_target_start(self):
        """ """
        return self.meta.start

    def get_target_buffer(self):
        """ """
        return self.meta.regionBuffer
