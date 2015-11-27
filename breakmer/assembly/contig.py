#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import logging
import shutil
import pysam
import breakmer.assembly.olc as olcAssembly
import breakmer.assembly.utils as assemblyUtils
import breakmer.assembly.assembler as assembler
import breakmer.realignment.realigner as realigner
import breakmer.caller.sv_caller as sv_caller
import breakmer.utils as utils
import breakmer.annotation.sv_annotation as annotator
import breakmer.plotting.sv_viz as svplotter

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"

'''
contig.py module contains classes to track assembly and post-assembly information for a
particular contig sequence.

It handles the actually assembly of a contig and all the information from realignment to the
reference sequence and the variant calling.

Classes:
    - AssemblyRead
    - ReadBatch
    - ContigCounts
    - Builder
    - Meta
    - Contig

Functions:
    - get_read_kmers
'''


def get_read_kmers(newSeq, kmerLen, kmerSeqs, order='for'):
    """Return new sample kmers from the existing contig sequence to help extend
    the contig sequence.

    All the k-length mers are determined from the newSeq. These kmer sequences are put
    into a set and intersected with the kmer sequences in the sample, ordered according
    to the position of the kmer in the newSeq string and returned.

    Args:
        newSeq (str):   Contig sequence to split into kmers.
        kmerLen (int):  Kmer length.
        kmerSeqs (set): The set of kmer sequences from the pool of extracted reads.
        order (str):    Direction to order the new set of kmer sequences. A None
                        value indicates no ordering.

    Returns:
        kmers: List of tuples containing:
               1. String kmer seq
               2. Integer kmer position
               3. Boolean if kmer seq is in the first half of the sequence
               4. Integer of position distance to middle of sequence
               5. String of how to order tuples in the list
    """

    m = len(newSeq) / 2
    # Generate a list of tuples.
    kmers = map(lambda x: (newSeq[x:(x + kmerLen)], x, int(x < m), abs(x - m), order), range(0, (len(newSeq) - kmerLen)))
    # Get unique kmer sequences.
    ks = set(map(lambda x: x[0], kmers))
    # Intersection of contig kmers and unused kmers from extracted reads.
    ss = ks & kmerSeqs
    # Filter out used kmers.
    kmers = filter(lambda x: x[0] in ss, kmers)
    if order == 'rev':
        kmers.reverse()
    elif order == 'mid':
        kmers = sorted(kmers, key=lambda x: (x[2], x[3]))
    return kmers


class AssemblyRead:
    """Wrapper class for a sequence read used in a contig assembly. This will
    track meta information about the sequence read.

    Attributes:
        read (fq_read object):  A previously generate fastq read object in utils.py
        redundant (boolean):    Indicator whether the read is duplicated.
        alignChecked (boolean): Indicator if the read has been checked against
                                the contig sequence.
        aligned (boolean):      Indicator if the read is aligned to the contig sequence.
    """

    def __init__(self, read, redundant, checked, aligned):
        self.read = read
        self.redundant = redundant
        self.alignChecked = checked
        self.aligned = aligned


class ReadBatch:
    """A class to track the reads that are being considered for building a contig
    sequence.

    This class stores the sequence reads that are passed in as fq_read objects and
    re-stores them as AssemblyRead objects.

    Attributes:
        delete (set): Set of fq_read objects to remove from further analysis.
        alt (list):   List of tuples containing (fq_read object, integer of nreads with the same sequence)
        reads (list): List of AssemblyRead objects containing the reads used to build a contig.
    """

    def __init__(self, read, merPos):
        """Initiate ReadBatch object with a fq_read object and the position of the kmer sequence
        within the read sequence.

        Args:
            read (fq_read): The read containing the kmer sequence.
            merPos (int):   The position of the kmer of interest in the read sequence.

        """

        self.delete = set()
        self.alt = []
        self.reads = [AssemblyRead(read, False, True, True)]

    def check_kmer_read(self, read):
        """Adds AssemblyRead to reads list.

        Args:
            read (fq_read): fq_read object.

        Returns:
            None
        """

        check = True
        # redund_read = False
        # add_read = True
        # add_to_pos_d = False
        self.reads.append(AssemblyRead(read, False, check, False))

    def set_last_read_aligned(self):
        """Sets the last read added to the reads list as aligned.

        Args:
            None

        Returns:
            None
        """

        self.reads[-1].aligned = True

    def clean(self, fq_reads, contigBuffer, last_keep_read):
        """Remove all data from data structures.

        Iterate through reads in delete set and delete them from the fq dictionary.
        Check if the delete reads are in the contigBuffer contig dictionary. If the
        contig associated with the read is not setup then delete the read from the dictionary.

        Args:
            fq_reads (dict):             Dictionary containing the extracted reads.
            contigBuffer (ContigBuffer): ContigBuffer object.
            last_keep_read (fq_read):    fq_read object kept for further use.

        Returns:
            None
        """

        map(fq_reads.__delitem__, map(lambda x: x.seq, list(self.delete)))
        for read_id in filter(lambda x: x in contigBuffer.contigs, list(self.delete)):
            if not contigBuffer.contigs[read_id].setup:
                del contigBuffer.contigs[read_id]
        self.delete = set()
        self.alt = []
        self.reads = [last_keep_read]


class ContigCounts:
    """A class to track the number of read sequences that support a consensus sequence.

    Initially set counts for the first read in the contig.

    Attributes:
        indel_only (list): List of integers, providing the count for number of indel only reads support
                           the given position in the consensus sequence.
        others (list): List of integers, providing the count of non indel only reads that are assembled
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
            p1 (int):      Integer indicating the first position.
            p2 (int):      Integer indicating the second position.
            sv_type (str): String indicating what kind of event the count is intended to support.

        Returns:
            counts (list): List of integers for counts of reads assembled at the provided range.
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
        """Return the total read count supporting a contig sequence.
        """

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

        Returns:
            None
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

    def set_counts(self, start, end, nReads, indelOnly):
        """Add the read count to the stored contig sequence count vectors.

        With paired end reads, there are reads that can contribute to indels only
        or to all types of variation. The counts are added according to how the
        read has been defined.

        Args:
            start (int):         Start of the sequence to add count.
            end (int):           End of the sequence to add count.
            nReads (int):        Total number of reads that should be added.
            indelOnly (boolean): Indicates if read should only support indels.

        Returns:
            None
        """

        if indelOnly:
            self.indelOnly[start:end] = map(lambda x: x + nReads, self.indelOnly[start:end])
        else:
            self.others[start:end] = map(lambda x: x + nReads, self.others[start:end])

    def extend_counts(self, extendSize, nReads, indelOnly, direction):
        """Increase the size of the count vectors when the contig sequence is grown.
        If the direction is 'post', the count vectors must be increased at the end.
        If the direction is 'pre', the count vectors must be increased at the beginning.

        Args:
            extendSize (int):    Number of positions to increase the vectors.
            nReads (int):        Count to add to the count vectors.
            indelOnly (boolean): Indicates if count should only be added to indel only vector.
            direction (str):     Indicates which side the count vector is extended.

        Returns:
            None
        """

        fill_counts = [0] * extendSize
        ecounts = [nReads] * extendSize
        if indelOnly:
            if direction == 'post':
                self.indelOnly.extend(ecounts)
                self.others.extend(fill_counts)
            else:
                ecounts.extend(self.indelOnly)
                self.indelOnly = ecounts
                fill_counts.extend(self.others)
                self.others = fill_counts
        else:
            if direction == 'post':
                self.indelOnly.extend(fill_counts)
                self.others.extend(ecounts)
            else:
                ecounts.extend(self.others)
                self.others = ecounts
                fill_counts.extend(self.indelOnly)
                self.indelOnly = fill_counts


class Builder:
    """A class to perform all the contig building functions and store temporary data structures.

    Attributes:
        readBatch (ReadBatch): ReadBatch object
        seq (str):             String of the consensus sequence.
        counts (ContigCounts): ContigCounts object to manage all the read counts supporting the consensus sequence.
        checkedKmers (list):   List of kmer sequences that had previously been checked while building the contig.
        kmerLen (int):         Integer of the kmer length.
        kmers (list):          List of kmer sequences that have contributed to building the contig.
        kmer_locs (list):      List of integers representing the positions of the kmers in the contig seq.
    """

    def __init__(self, assemblyKmer, readAlignValues):
        """Initiate the Builder object. This will create a ReadBatch object and ContigCounts object.

        Args:
            assemblyKmer (AssemblyKmer): Instance of AssemblyKmer object with attributes for kmer sequence.
            readAlignValues (dict):      Contains keys:
                                          - 'read': fq_read object that contains kmer sequence.
                                          - 'align_pos': Integer position of kmer in read sequence
                                          - 'nreads': Integer of number of reads with the same sequence.
        """

        self.readBatch = ReadBatch(readAlignValues['read'], readAlignValues['align_pos'])
        self.seq = readAlignValues['read'].seq
        self.counts = ContigCounts(readAlignValues['read'], readAlignValues['nreads'])
        self.checkedKmers = [assemblyKmer.seq]
        self.kmerLen = assemblyKmer.kmerLen
        self.kmers = []
        self.kmer_locs = []

    def check_read(self, assemblyKmer, readAlignValues, alignType):
        """Determine if the read should be added to the assembly or not.

        If the read aligns to the contig, set the fq_read status to used and indicate
        the AssemblyRead has been aligned. If the kmer is in more than 1 read and the
        current read has not been used in any other contigs then store for later
        analysis to build another contig possibly. Otherwise, discard the read for
        further analysis.

        Args:
            assemblyKmer (AssemblyKmer): Kmer object containing kmer seq specific values.
            readAlignValues:             Dictionary containing:
                                          - 'read': fq_read object that contains kmer sequence.
                                          - 'align_pos': Integer position of kmer in read sequence
                                          - 'nreads': Integer of number of reads with the same sequence.
            type (str):                  Indicates the state of this function.

        Returns:
            hit (str): String value 'remove' or ''.
        """

        hit = ''
        self.readBatch.check_kmer_read(readAlignValues['read'])  # Add sequence read with kmer to ReadBatch list.
        if self.check_align(assemblyKmer, readAlignValues, alignType):
            hit = 'remove'
            readAlignValues['read'].used = True
            self.readBatch.set_last_read_aligned()
        elif assemblyKmer.counts > 2 and not readAlignValues['read'].used:  # Read did not align to contig, save for later.
            self.readBatch.alt.append((readAlignValues['read'], readAlignValues['nreads']))
        else:  # Read is the only one with kmer sequence or it is used. Remove it.
            self.readBatch.delete.add(readAlignValues['read'])
        return hit

    def check_align(self, assemblyKmer, readAlignValues, alignType='setup'):
        """Check the alignment of the read sequence to the contig sequence.

        The read sequence must match at least 25% of the shortest sequence between
        the contig and the read and an identity at least 90%. If there is clear
        alignment between the read sequence and the contig sequence, then the
        assembly consensus sequence is appropriately changed.

        There are a couple of possibilities when trying to align a read sequence to
        an initiate contig sequence.

        Return a True to indicate the read sequence aligned to the contig sequence.
        In the situation that the contig is already established, the read ID will be
        removed from the contigBuffer dictionary to indicate that it is no longer in the
        buffer for considering an alternate contig sequence.

        1. Alignment thresholds are not met - read sequence does not align (or well enough)
           to the contig sequence.
            - Return False
        2. The read sequence and contig sequence are the same exact sequence.
            - Return True
        3. The alignment meets thresholds and there is overlap between them
            - Return True
            1. Alignment scores are the same when doing it both ways (i.e., read seq vs. contig seq, contig seq vs. read seq)
                1. The read sequence is a super sequence of the contig sequence - replace contig sequence with read sequence.
                2. The read sequence is a subset sequence of the contig sequence - update the counts vector.
            2. Read sequence overlaps off the front of the contig sequence
            3. Read sequence overlaps off the end of the contig sequence.

        Args:
            assemblyKmer (AssemblyKmer): Kmer object containing kmer seq specific values.
            readAlignValues:             Dictionary containing:
                                          - 'read': fq_read object that contains kmer sequence.
                                          - 'align_pos': Integer position of kmer in read sequence
                                          - 'nreads': Integer of number of reads with the same sequence.
            type (str):                  Indicates the state of this function - 'setup' or 'grow'

        Returns:
            match (boolean): Indicates if the read aligns sufficiently with the
                             contig sequence and will be added.
        """

        match = False
        queryRead = readAlignValues['read']

        minScore = float(min(len(self.seq), len(queryRead.seq))) / 4.0
        alignManager = olcAssembly.AlignManager(self.seq, queryRead.seq, minScore, 0.90)

        if alignManager.check_align_thresholds():
            return False  # Alignment failed minimum alignment thresholds.
        if alignManager.same_seqs():
            return True   # Read and contigs sequences are the same.
        if alignManager.same_max_scores():  # Alignments both ways had equal scores.
            match = True
            if alignManager.read_is_superseq():  # Check if the read sequence fully contains the contig sequence.
                self.set_superseq(queryRead, readAlignValues['nreads'], alignManager.get_alignment_values(0, 'i'), alignManager.get_alignment_values(0, 'prei'))
                if alignType == 'grow':
                    # Contig sequence has changed, set the kmers.
                    self.set_kmers(assemblyKmer.kmerSeqSet)
            elif alignManager.read_is_subseq():  # Check if the contig sequence full contains the read sequence.
                self.add_subseq(alignManager.get_alignment_values(1, 'i'), alignManager.get_alignment_values(1, 'prei'), readAlignValues['nreads'], queryRead.indel_only)
            else:  # There appears to be overlap, figure out how to assemble.
                match = False
                indx1 = alignManager.get_kmer_align_indices(0, assemblyKmer.seq)
                indx2 = alignManager.get_kmer_align_indices(1, assemblyKmer.seq)
                if indx1[0] > -1 and indx1[1] > -1:
                    # Read overlaps off front of contig sequence.
                    if (indx2[0] == -1 and indx2[1] == -1) or (abs(indx2[0] - indx2[1]) > abs(indx1[0] - indx1[1])):
                        match = True
                        self.contig_overlap_read(alignManager.get_alignment(0), queryRead, readAlignValues['nreads'], assemblyKmer.kmerSeqSet, alignType)
                elif indx2[0] > -1 and indx2[1] > -1:
                    # Read overlaps off end of contig sequence.
                    if (indx1[0] == -1 and indx1[1] == -1) or (abs(indx2[0] - indx2[1]) < abs(indx1[0] - indx1[1])):
                        match = True
                        self.read_overlap_contig(alignManager.get_alignment(1), queryRead, readAlignValues['nreads'], assemblyKmer.kmerSeqSet, alignType)
        elif alignManager.better_align():  # Read sequence overlaps off the front of the contig sequence.
            match = True
            self.contig_overlap_read(alignManager.get_alignment(0), queryRead, readAlignValues['nreads'], assemblyKmer.kmerSeqSet, alignType)
        else:  # Read sequence overlaps off the end of the contig sequence.
            match = True
            self.read_overlap_contig(alignManager.get_alignment(1), queryRead, readAlignValues['nreads'], assemblyKmer.kmerSeqSet, alignType)
        return match

    def set_superseq(self, read, nreads, start, end):
        """The read sequence contains the current contig sequence. Change the contig sequence to the
        read seqeunce and modify the counts list.

        Args:
            read (fq_read): fq_read object.
            nreads (int):   The number of reads with the same sequence as the read passed in.
            start (int):    The start position the contig sequence aligns to the read sequence.
            end (int):      The end position the contig sequence aligns to the read sequence.

        Returns:
            None
        """

        self.seq = read.seq
        self.counts.set_superseq(read, nreads, start, end)

    def add_subseq(self, start, end, nreads, indelOnly):
        """The read checked against the contig was found to be a subsequence of the
        contig. The nreads with the checked sequence are added to the count vectors.

        Args:
            start (int):         Start of sequence to add counts.
            end (int):           End of sequence to add counts.
            nreads (int):        Number of reads to add.
            indelOnly (boolean): Indicates whether to add to indel only count vector.

        Returns:
            None
        """

        self.counts.set_counts(start, end, nreads, indelOnly)

    def add_postseq(self, postSeq, start, end, nreads, indelOnly):
        """Sequence is appended to the end of the current contig sequence. The
        read support vectors are appropriately incremented.

        Args:
            postSeq (str):       Sequence to add to the end of the assembled contig.
            start (int):         Start of sequence to add counts.
            end (int):           End of sequence to add counts.
            nreads (int):        Number of reads to add.
            indelOnly (boolean): Indicates whether to add to indel only count vector.

        Returns:
            None
        """

        self.seq += postSeq
        self.counts.set_counts(start, end, nreads, indelOnly)
        self.counts.extend_counts(len(postSeq), nreads, indelOnly, 'post')

    def add_preseq(self, preSeq, start, end, nreads, indelOnly):
        """Sequence is append to the front of the current contig sequence. The read
        support vectors are appropriately changed.

        Args:
            preSeq:             Sequence to add to the front of the assembly contig.
            start (int):         Start of sequence to add counts.
            end (int):           End of sequence to add counts.
            nreads (int):        Number of reads to add.
            indelOnly (boolean): Indicates whether to add to indel only count vector.

        Returns:
            None
        """

        self.seq = preSeq + self.seq
        self.counts.set_counts(start, end, nreads, indelOnly)
        self.counts.extend_counts(len(preSeq), nreads, indelOnly, 'pre')

    def finalize_reads(self, contigReads, fqRecs, contigBuffer):
        """Sort out the reads to keep for reporting and remove the others.
        Aligned and non-redundant reads are removed from the contig read set. The
        variables in read_batch are cleared.

        Args:
            contigReads (set):           Set of fq_read objects.
            fqRecs (dict):               Dictionary of fq_read objects.
            contigBuffer (ContigBuffer): ContigBuffer class object.

        Returns:
            contigReads: Set of fq_reads objects
        """

        rm_reads = map(lambda y: y.read, filter(lambda x: x.redundant, self.readBatch.reads))
        keep_reads = filter(lambda x: x.aligned and not x.redundant, self.readBatch.reads)
        add_reads = map(lambda y: y.read, keep_reads)
        # Merge add_reads into contigReads
        contigReads = contigReads | set(add_reads)
        # Remove rm_reads
        contigReads = contigReads - set(rm_reads)
        self.readBatch.clean(fqRecs, contigBuffer, keep_reads[-1])
        return contigReads

    def contig_overlap_read(self, alignment, queryRead, nreads, kmerSeqs, assemblyType):
        """Assembled consensus and read sequences, where the consensus end overlaps
        with the read sequence beginning.

        Args:
            alignment (olc.Align): olc.Align object
            queryRead (fq_read):   fq_read object
            nreads (int):          The number of reads to add to count vectors.
            kmerSeqs (set):        Set of kmer sequence values.
            type (str):            Indicates the source of call to function.

        Returns:
            None
        """

        if alignment.prej == len(self.seq) and alignment.j == 0:
            self.set_superseq(queryRead, nreads, alignment.i, alignment.prei)
            if assemblyType == 'grow':
                self.set_kmers(kmerSeqs)
        else:
            postSeq = queryRead.seq[alignment.prei:]
            nseq = self.seq[(len(self.seq) - (self.kmerLen - 1)):] + postSeq
            self.add_postseq(postSeq, alignment.j, alignment.prej, nreads, queryRead.indel_only)
            if assemblyType == 'grow':
                nkmers = get_read_kmers(nseq, self.kmerLen, kmerSeqs, 'for')
                self.kmers.extend(nkmers)

    def read_overlap_contig(self, alignment, queryRead, nreads, kmerSeqs, type):
        """Assemble consensus and read sequences togheter, where the consensus
        beginning overlaps the read sequence end.

        Args:
            alignment (olc.Align): olc.Align object
            queryRead (fq_read):   fq_read object
            nreads (int):          The number of reads to add to count vectors.
            kmerSeqs (set):        Set of kmer sequence values.
            type (str):            Indicates the source of call to function.

        Return: None
        """

        if alignment.prej == len(queryRead.seq) and alignment.j == 0:
            self.add_subseq(alignment.i, alignment.prei, nreads, queryRead.indel_only)
        else:
            preSeq = queryRead.seq[0:alignment.j]
            nseq = preSeq + self.seq[0:(self.kmerLen - 1)]
            self.add_preseq(queryRead.seq[0:alignment.j], alignment.i, alignment.prei, nreads, queryRead.indel_only)
            if type == 'grow':
                nkmers = get_read_kmers(nseq, self.kmerLen, kmerSeqs, 'rev')
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
        for read, nreads in self.readBatch.alt:
            altKmers = get_read_kmers(read.seq, self.kmerLen, kmerTracker.kmerSeqs, '')
            altKmerSeqs = set(map(lambda x: x[0], altKmers))
            newKmers = set(altKmerSeqs) - set(contigKmers) - contigBuffer.used_kmers - kmerSet
            if len(newKmers) > 0:
                for kmerSeq in list(newKmers):
                    readCount = kmerTracker.get_count(kmerSeq)
                    if readCount > 1:
                        kmerPos = read.seq.find(kmerSeq)
                        assemblyKmer = assembler.AssemblyKmer(kmerSeq, kmerTracker.get_count(kmerSeq), kmerTracker.kmerSeqs, self.kmerLen)
                        read_align_values = {'read': read,
                                             'align_pos': kmerPos,
                                             'nreads': nreads}
                        newContigs.append((read, Contig(assemblyKmer, read_align_values)))
                        kmerSet = kmerSet | newKmers
                        break
        return newContigs

    def set_kmers(self, kmerSeqSet):
        """Wrapper function to get_read_kmers function to parse a sequence string
        and generate all relevant kmers from the sequence.

        Args:
            kmerSeqSet (set): Set of kmer sequences.

        Returns:
            None
        """

        self.kmers = get_read_kmers(str(self.seq), self.kmerLen, kmerSeqSet, 'mid')

    def set_kmer_locs(self):
        """Add the start alignment positions of each kmer sequence in the kmers list
        to the kmer_locs list.

        Args:
            None

        Returns:
            None
        """

        self.kmer_locs = [0] * len(self.seq)
        for kmer in self.kmers:
            kmerPos = self.seq.find(kmer[0])
            self.kmer_locs[kmerPos:(kmerPos + self.kmerLen)] = map(lambda x: x + 1, self.kmer_locs[kmerPos:(kmerPos + self.kmerLen)])

    def refresh_kmers(self):
        """Return a list of kmer_sequences that have not been checked already.

        Args:
            None

        Returns:
            List of kmer sequences.
        """

        return filter(lambda x: x[0] not in set(self.checkedKmers), self.kmers)

    def get_seq(self):
        """Return the final consensus sequence.
        """

        return self.seq

    def get_kmers(self):
        """Return the final kmer sequence list.
        """

        return self.kmers

    def get_kmer_locs(self):
        """Return the final kmer locations list.
        """

        return self.kmer_locs

    def get_total_reads(self):
        """Return total number of reads supporting contig.
        """

        return self.counts.get_total_reads()


class Meta:
    """A class to track the contig information for downstream calling and writing
    to file.

    Attributes:
        params:        Param object with all BreaKmer parameters.
        path:          String of path to write all files.
        id:            String for contig ID.
        target_region: Tuple containing target information:
                       1. String chromosome ID
                       2. Integer of target region start position.
                       3. Integer of target region end position.
                       4. String target region name.
                       5. List of target intervals tuples.
        fq_fn:         String of the fastq file containing sequence reads used to build contig.
        fa_fn:         String of the fasta file containing the contig sequence.
    """

    def __init__(self):
        """
        """

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

        Returns:
            None
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
        """
        """

        return (self.chr, self.start, self.end, self.targetName, self.regionBuffer)

    def write_files(self, cluster_fn, kmers, reads, seq):
        """Write cluster, read fastq, and contig fasta files.

        Args:
            cluster_fn: String of the file to write kmer clusters to.
            kmers: List of kmers used in the building of the contig.
            reads: List of reads used in the building of the contig.
            seq: String of contig sequence.

        Returns:
            None
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
        """
        """

        resultFn = os.path.join(self.path, self.id + "_svs.out")
        utils.log(self.loggingName, 'info', 'Writing %s result file %s' % (self.id, resultFn))
        resultFile = open(resultFn, 'w')

        # A string of output values for writing to file.
        headerStr, formattedResultValuesStr = svEventResult.get_formatted_output_values()
        resultFile.write(headerStr + '\n' + formattedResultValuesStr + '\n')
        resultFile.close()
        shutil.copyfile(resultFn, os.path.join(outputPath, self.id + "_svs.out"))

    def write_bam(self, outputPath, svBamReadsFn, reads):
        """
        """

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
        meta (Meta):       Meta class object to store all the interface related data.
        kmer_locs (list):  List of integers indicating the start alignment position of the kmers
                           in the contig sequence.
        setup (boolean):   Boolean to indicate whether a contig has gone through the setup process.
        builder (Builder): Builder class object that handles all the assembly functions.
        seq (str):         String for assembled sequence.
        kmers (list):      List of kmer sequences used to build contig.
        reads (set):       Set of read IDs that have been used to build contig.
        buffer (set):      Set of read IDs used in a batch of processing for building a contig. This is flushed.
    """

    def __init__(self, assemblyKmer, readAlignValues):
        """Initiate a Contig object with a kmer sequence and the reads that contain the kmer.

        Args:
            assemblyKmer (AssemblyKmer): Instance of AssemblyKmer object with attributes for kmer sequence.
            readAlignValues (dict):      Contains keys:
                                          - 'read': fq_read object that contains kmer sequence.
                                          - 'align_pos': Integer position of kmer in read sequence
                                          - 'nreads': Integer of number of reads with the same sequence.
        """

        self.meta = Meta()
        self.setup = False
        self.builder = Builder(assemblyKmer, readAlignValues)
        self.seq = None
        self.kmers = []
        self.kmer_locs = []
        self.reads = set()
        self.buffer = set([readAlignValues['read'].id])
        self.svEventResult = None
        self.realignment = None

    def check_read(self, assemblyKmer, readAlignValues, fncType='setup'):
        """Wrapper function to Builder class check_read function.

        Check if the read passed in can be added to the current contig.

        Args:
            assemblyKmer (AssemblyKmer): Instance of AssemblyKmer object with attributes for kmer sequence.
            readAlignValues (dict):      Contains keys:
                                          - 'read': fq_read object that contains kmer sequence.
                                          - 'align_pos': Integer position of kmer in read sequence
                                          - 'nreads': Integer of number of reads with the same sequence.
            fncType (str):               Indicates the state of this function.

        Returns:
            String containing 'hit' or '' indicating that read matched contig seq
            or did not, respectively.
        """

        self.buffer.add(readAlignValues['read'].id)
        return self.builder.check_read(assemblyKmer, readAlignValues, fncType)

    def finalize(self, fqRecs, kmerTracker, contigBuffer, source='setup'):
        """Finish an assembly and add the buffered contigs that were created from
        non-aligned reads to the contigBuffer.

        Args:
            fqRecs (dict):               Dicionary of fq_read objects.
            kmerTracker (KmerTracker):   KmerTracker object with all kmer sequence values.
            contigBuffer (ContigBuffer): ContigBuffer object.
            source (str):                Source of function call - 'setup' or 'grow'.

        Returns:
            None
        """

        if source == 'setup':
            self.set_kmers(kmerTracker.kmerSeqs)

        # Get alternate read kmers and see if any are different from contig kmers.
        newContigs = self.builder.check_alternate_reads(kmerTracker, contigBuffer, self.kmers)
        for newContig in newContigs:
            contigBuffer.add_contig(newContig[0], newContig[1])
        self.reads = self.builder.finalize_reads(self.reads, fqRecs, contigBuffer)

    def check_invalid(self, readCountThreshold, readLen):
        """Determine if the finished contig sequence meets minimum requirements for
        length and read count support.

        Args:
            readCountThreshold (int): Minimum number of reads that must support a contig.
            readLen (int):            Read length.

        Returns:
            Boolean indicating whether it meets thresholds or not.
        """

        if (self.get_total_read_support() < int(readCountThreshold)) or (len(self.seq) <= int(readLen)):
            return True
        else:
            return False

    def set_kmers(self, kmerSeqSet):
        """Wrapper function to Builder class set_kmers function.

        The Builder object will use these to determine the set of kmers that still needs to
        be analyzed that exist in the contig.

        Args:
            kmerSeqSet (set): Set of all unused kmer sequences.

        Returns:
            None
        """

        self.setup = True
        self.builder.set_kmers(kmerSeqSet)

    def set_kmer_locs(self):
        """Wrapper function to Builder class set_kmer_locs function.

        Args:
            None
        Returns:
            None
        """

        self.builder.set_kmer_locs()

    def refresh_kmers(self):
        """A wrapper function to Builder class refresh kmers.

        Args:
            None

        Returns:
            List of kmers that have not been previously checked.
        """

        return self.builder.refresh_kmers()

    def get_kmer_reads(self, kmerValues, readItems):
        """
        Args:
            kmerValues: Tuple containing the alignment information of a kmer sequence
                         in a read sequence.
                         1. String kmer sequence.
                         2. Integer kmer alignment position in sequence.
                         3. Boolean whether kmer align position is below midpoint of sequence.
                         4. Integer of absolute difference between align position and midpoint.
                         5. String of the order for tuples in a list.
            readItems: List of tuples for sequence reads:
                        1. String read sequence.
                        2. List of fq_read objects with read sequence.

        Returns:
            reads: List of tuples containing:
                   1. read object,
                   2. start position of kmer match in read seq
                   3. Boolean that a match was found.
                   4. Length of the read sequence.
                   5. Number of reads with this sequence.
        """

        kmer, kmerPos, lessThanHalf, dist_half, order = kmerValues
        read_order = 'for'
        if order == 'mid':
            if lessThanHalf == 0:
                read_order = 'rev'
        elif order == 'for':
            read_order = 'rev'
        reads = assemblyUtils.find_reads(kmer, readItems, self.buffer, read_order)
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

        Returns:
            None
        """

        logger = logging.getLogger('breakmer.assembly.contig')
        if not self.setup:
            self.set_kmers(kmerTracker.kmerSeqs)
        newKmers = self.refresh_kmers()
        while len(newKmers) > 0:
            iter = 0
            for kmer_lst in newKmers:
                kmerSeq, kmerPos, lessThanHalf, dist_half, order = kmer_lst
                reads = self.get_kmer_reads(kmer_lst, fqRecs.items())
                contigBuffer.add_used_mer(kmerSeq)
                assemblyKmer = assembler.AssemblyKmer(kmerSeq, kmerTracker.get_count(kmerSeq), kmerTracker.kmerSeqs, kmerLen)
                for readLst in reads:
                    read, kmerPos, bool, rlen, nreads = readLst
                    contigBuffer.add_used_read(read.id)
                    readAlignValues = {'read': read,
                                       'align_pos': kmerPos,
                                       'nreads': nreads}
                    hit = self.check_read(assemblyKmer, readAlignValues, 'grow')
                    if hit == 'remove':
                        contigBuffer.remove_contig(read.id)
                self.finalize(fqRecs, kmerTracker, contigBuffer, 'grow')
                self.builder.checkedKmers.append(kmerSeq)
                iter += 1
            newKmers = self.refresh_kmers()
            logger.debug("%d kmers left to check" % len(newKmers))
        self.set_kmer_locs()
        self.set_final_values()
        logger.info('Contig done with contig seq %s. Supported by %d read(s).' % (self.seq, len(self.reads)))
        logger.info('Read IDs: %s' % (",".join([x.id for x in list(self.reads)])))

    def set_meta_information(self, contigId, params, queryRegionValues, contigPath, kmerClusterFn, readVariation):
        """Sets the contig ID, params, target region values and contig path variables for later use.
        Output files are also written with the contig assembly information.

        Args:
            contigId (str):            Contid ID.
            params (ParamManger):      ParamManager object.
            queryRegionValues (tuple): Tuple containing the target region information
            contigPath (str):          Path to the contig directory to store files.
            kmerClusterFn (str):       Path to write the kmer clustering information to.

        Returns:
            None
        """

        self.meta.set_values(contigId, params, queryRegionValues, contigPath, readVariation)
        self.meta.write_files(kmerClusterFn, self.kmers, self.reads, self.seq)

    def set_final_values(self):
        """Set the seq, kmers, kmer_locs variables when the contig is done assemblying.

        Args:
            None

        Returns:
            None
        """

        self.seq = self.builder.get_seq()
        self.kmers = self.builder.get_kmers()
        self.kmer_locs = self.builder.get_kmer_locs()

    def query_ref(self, targetRefFns):
        """

        Args:
            targetRefFns (list):

        Returns:
            None
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
            svFilter.check_filters(self.svEventResult)

    def annotate_calls(self):
        """
        """

        if self.svEventResult and self.meta.params.get_param('gene_annotation_file') and self.meta.params.get_param('bedtools'):
            annotator.annotate_event(self.svEventResult, self.meta)

    def output_calls(self, outputPath, svReadsBamFn):
        """
        """

        if self.svEventResult:
            self.meta.write_result(self.svEventResult, outputPath)
            readBamFn = self.meta.write_bam(outputPath, svReadsBamFn, self.reads)
            if self.meta.params.get_param('generate_image') and not self.svEventResult.is_filtered():
                # Generate image if option is set and the result is not being filtered out.
                svplotter.generate_pileup_img(self.svEventResult, readBamFn, outputPath, self.get_id())

    def get_total_read_support(self):
        """Return the total read count supporting assembly.
        """

        return self.builder.get_total_reads()

    def get_contig_len(self):
        """Return length of contig sequence.
        """

        return len(self.seq)

    def get_kmer_locs(self):
        """Return the kmer locations in the contig sequence.
        """

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
        """
        """

        return self.meta.get_target_region_coordinates()

    def get_chr(self):
        """
        """

        return self.meta.chr

    def get_target_start(self):
        """
        """

        return self.meta.start

    def get_target_buffer(self):
        """
        """

        return self.meta.regionBuffer
