#! /usr/bin/local/python
# -*- coding: utf-8 -*-

import re
import logging
from collections import OrderedDict
import breakmer.assembly.contig as contig_assembler

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


def init_assembly(kmers, fqRecs, kmerLen, rcThresh, readLen):
    """Entry function for assemblying a contiguous sequence from
    a pool of sample only kmers and the reads that contain them.
    A kmer tracker object is instantiated containing all the kmer seqs and
    their associated counts. These are sorted by
    Args:
        kmers: Dictionary of kmers only in the sample key = kmer, value = count in reads
        fqRecs: Dictionary with sequence values as keys and a list of fq_read objects.
        kmerLen: Integer of kmer size.
        rcThresh: Integer representing the minimum readcount threshold for keeping a contig.
        readLen: Integer of the read length.
    Return:
        contigs: List of contig objects.
    """
    logger = logging.getLogger('breakmer.assembly.assembler')
    contigs = []

    # Return if no kmers to analyze.
    if len(kmers) == 0:
        logger.info('No kmers to built contigs, returning.')
        return contigs

    # Store kmers in KmerTracker object.
    kmerTracker = KmerTracker()
    for kmer in kmers:
        kmerTracker.add_kmer(kmer, kmers[kmer])

    # While there are kmers to analyze continue to build contigs.
    contigBuffer = ContigBuffer()
    # Sort all the kmers by count and store in order.
    kmerTracker.set_all_kmer_values()
    # Check if there are any kmers left to seed the build process.
    while kmerTracker.has_mers():
        # Update the set of kmers to consider for building.
        kmerTracker.update_kmer_set()
        # Get kmer seed for new contig.
        kmer, kmer_count = kmerTracker.get_kmer()
        # Only analyze contigs that exist in 2 or more reads.
        if kmer_count < 2:
            continue
        logger.info('Initiating kmer %s, found in %d reads' % (kmer, kmer_count))
        setup_contigs(kmer, fqRecs, kmerLen, kmerTracker, contigBuffer)

        # Deal with buffered contig objects that need to be grown or completed.
        while len(contigBuffer.contigs) > 0:
            contig = contigBuffer.get_contig()
            contig.grow(fqRecs, kmerTracker, kmerLen, contigBuffer)
            if contig.check_valid(rcThresh, readLen):
                logger.info('Contig did not meet the read count threshold %d, with %d or contig length (%d) < readLen (%d)' % (rcThresh, len(contig.reads), len(contig.seq), readLen))
            else:
                logger.info('Adding contig to buffer')
                contigs.append(contig)

        # Clean up the data to free up memory.
        contigBuffer.remove_kmers(kmerTracker)
        contigBuffer.remove_reads(fqRecs)
    return contigs


def setup_contigs(kmerSeq, fqRecs, kmerLen, kmerTracker, contigBuffer):
    """Create a contig instance starting with a seed kmer and associated reads.
    First find the reads containing the kmerSeq value, iterate through reads and
    either create a new contig or add to existing contig.
    Args:
        kmerSeq:        String of kmer sequence.
        fqRecs:         Dictionary with sequence values as keys and a list of fq_read objects.
        kmerLen:        Integer of kmer size.
        kmerTracker:    KmerTracker object that contains all the kmer values.
        contigBuffer:   ContigBuffer object to track the buffered contig objects.
    Return: None
    """
    logger = logging.getLogger('breakmer.assembly.assembler')
    contig = None
    # Find all reads with kmer sequence passed in.
    # kmerReads contains a list of tuples.
    #   1. fq_read object defined in breakmer.utils.py
    #   2. Starting position of the kmer match in the read sequence
    #   3. Boolean that a match was found.
    #   4. Length of the read sequence.
    #   5. Number of reads with this sequence.
    kmerReads = find_reads(kmerSeq, fqRecs.items(), set())
    contigBuffer.add_used_mer(kmerSeq)

    kmerObj = Kmer(kmerSeq, kmerTracker.get_count(kmerSeq), kmerTracker.kmerSeqs, kmerLen)
    for readVals in kmerReads:
        read, kmerPos, matchFound, seqLen, nReadsWithSeq = readVals
        readAlignValues = {'read': read,
                           'align_pos': kmerPos,
                           'nreads': nReadsWithSeq}
        contigBuffer.add_used_read(read.id)
        # If no contig, build one.
        if not contig:
            contig = contig_assembler.Contig(kmerObj, readAlignValues)
            contigBuffer.add_contig(read, contig)
        # Check if read should be added to the existing contig.
        else:
            contig.check_read(kmerObj, readAlignValues, 'setup')
    if contig:
        contig.finalize(fqRecs, kmerTracker, contigBuffer, 'setup')


def find_reads(kmerSeq, readItems, usedReads, order='for'):
    """Return a list of tuples containing information from reads with the kmer sequence.
    First search all the read sequences for the given kmer sequence. Then,
    filter out used reads and order them according to position of the kmer
    sequence in the read sequence.
    Args:
        kmerSeq: String of kmer sequence.
        readItems: List of fq_recs (key, value) tuples.
        usedReads: Set of read IDs that have been previously used.
        order: String indicating how the list of the identified reads
               should be ordered.
    Return:
        kmerReads: List of tuples containing:
                    1. read object,
                    2. start position of kmer match in read seq
                    3. Boolean that a match was found.
                    4. Length of the read sequence.
                    5. Number of reads with this sequence.
    """
    kmerReads = []
    # Filter all the reads not containing the kmerSeq
    mappedReads = filter(lambda x: x[2], map(read_search, [kmerSeq] * len(readItems), readItems))
    # Extract the read ids of the reads containing kmerSeq
    mappedReadIds = map(lambda x: x[0].id, mappedReads)
    filterIds = set(mappedReadIds) - set(usedReads)
    matchedReads = filter(lambda x: (x[0].id in filterIds), mappedReads)
    if order == 'rev':
        kmerReads = sorted(matchedReads, key=lambda z: (-z[1], -z[3]))
    else:
        kmerReads = sorted(matchedReads, key=lambda z: (z[1], -z[3]))
    return kmerReads


def read_search(kmerSeq, readItems):
    """Return a tuple containing information regarding the alignment of the kmerSeq
    in a sequence read.
    This uses regex searching function re.search to determine if the kmerSeq
    is contained in the read sequence. If so, then it returns a 5 element
    tuple about information regarding this alignment. If no match, then return
    a 3 element tuple with None values.
    Args:
        kmerSeq: String of kmer sequence.
        readItems: List of fq_recs (key, value) tuples.
    Return:
        searchResult: Tuple of result information.
                       1. read object,
                       2. start position of kmer match in read seq
                       3. Boolean that a match was found.
                       4. Length of the read sequence.
                       5. Number of reads with this sequence.
    """
    searchResult = (None, None, None)
    seq, reads = readItems
    x = re.search(kmerSeq, seq)
    if x:
        searchResult = (reads[0], x.start(), True, len(reads[0].seq), len(reads))
    return searchResult


class Kmer:
    """Class to track value associated with a particular kmer sequence.
    Attributes:

    """
    def __init__(self, seq, counts, kmerSeqSet, kmerLen):
        self.seq = seq
        self.counts = counts
        self.kmerSeqSet = kmerSeqSet
        self.kmerLen = kmerLen


class ContigBuffer:
    """A class to track the used kmers and reads and their relation to contigs.
    Attributes:
        used_kmers: Set of kmer sequences that have been used to build contigs.
        used_reads: Set of read IDs that have been used to build contigs.
        contigs: OrderedDict to track reads and the contigs they contribute to.
    """
    def __init__(self):
        self.used_kmers = set()
        self.used_reads = set()
        self.contigs = OrderedDict()

    def add_contig(self, read, contig):
        """Add read to contigs dict with contig object it is connected to.
        Set key to read ID and value to the contig object. Set the read used to True.
        Args:
            read:   fq_read object
            contig: Contig object.
        Return:
            None
        """
        # Tie a contig to the seed read ID and store in dictionary.
        if read.id not in self.contigs and not read.used:
            self.contigs[read.id] = contig
            read.used = True

    def remove_contig(self, read_id):
        """Remove read ID from contigs dictionary.
        Args:
            read_id: String of read ID.
        Return: None
        """
        if read_id in self.contigs:
            del self.contigs[read_id]

    def get_contig(self):
        """Return the contig associated with the first record the contigs dictionary.
        Delete the entry.
        Args: None
        Return:
            Contig object to grow or complete.
        """
        read_id = self.contigs.keys()[0]
        contig = self.contigs[read_id]
        del self.contigs[read_id]
        return contig

    def add_used_read(self, read_id):
        """Add read ID to used set.
        Args:
            read_id: String for read ID.
        Return: None
        """
        self.used_reads.add(read_id)

    def add_used_mer(self, kmer_seq):
        """Add kmer sequence to used set.
        Args:
            kmer_seq: String for kmer sequence.
        Return: None
        """
        self.used_kmers.add(kmer_seq)

    def remove_kmers(self, kmer_tracker):
        """Remove used kmer sequences from the kmer tracking object and reset used
        kmer set.
        Args:
            kmer_tracker: KmerTracker object.
        Return: None
        """
        map(kmer_tracker.remove_kmer, list(self.used_kmers))
        self.used_kmers = set()

    def remove_reads(self, fqReads):
        """Remove the used reads from the fq_reads dictionary.
        Args:
            fqReads: Dictionary of fq_reads.
        Return: None
        """
        del_used = filter(lambda x: x in fqReads, list(self.used_reads))
        map(fqReads.__delitem__, del_used)
        self.used_reads = set()


class KmerTracker:
    """Wrapper class for storing the kmer objects. Useful for adding
    and extracting kmers.
    Attributes:
        kmers:          List of tuples containing kmer count, kmer, kmer object.
        orderedKmers:   OrderedDict object containing kmer seq as key and kmer count as value.
                        The top values are the most frequence kmer values.
        kmerSeqs:       Set of kmer seq values that exist in orderedKmers.
    """
    def __init__(self):
        self.kmers = []
        self.orderedKmers = OrderedDict()
        self.kmerSeqs = set()

    def add_kmer(self, mer, count):
        """Add a kmer object to the list. Stores a tuple with kmer count and kmer sequence string.
        This allows easy sorting.
        Args:
            mer:    String kmer sequence value.
            count:  Integer of number of reads kmer is within.
        Return:
            None
        """
        if len(set(mer)) > 1:
            self.kmers.append((int(count), mer))

    def set_all_kmer_values(self):
        """Sort the kmer list by number of reads (descending) first and then
        by sequence value and store them in an ordered dictionary.
        Args:
            None
        Return:
            None
        """
        kmersSorted = sorted(self.kmers, key=lambda x: (int(x[0]), x[1]), reverse=True)
        for kmer in kmersSorted:
            self.orderedKmers[kmer[1]] = kmer[0]

    def has_mers(self):
        """Check if there are any kmers left in the dictionary.
        Args:
            None
        Return:
            True if there are items in the dictionary and the counts of those items are > 1.
            False if there are no items in the dictionary or the counts of those items are <= 1.
        """
        if len(self.orderedKmers) > 0 and max(self.orderedKmers.values()) > 1:
            return True
        else:
            return False

    def update_kmer_set(self):
        """Update the set of kmer values. The orderedKmers dictionary
        dynamically changes as kmers are taken out.
        Args:
            None
        Return:
            None
        """
        self.kmerSeqs = set(self.orderedKmers.keys())

    def get_kmer(self):
        """Return the first kmer in the ordered dictionary"""
        return self.orderedKmers.items()[0]

    def get_count(self, kmerSeq):
        """Return the number of reads the kmer_seq is within."""
        return self.orderedKmers[kmerSeq]

    def remove_kmer(self, kmerSeq):
        """Delete the record associated with kmer sequence."""
        del self.orderedKmers[kmerSeq]
