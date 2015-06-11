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
    kmerTracker.set_all_kmer_values()
    while kmerTracker.has_mers():
        kmerTracker.update_kmer_set()
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


def setup_contigs(kmer_seq, fq_recs, kmer_len, kmer_tracker, contig_buffer):
    """Create a contig instance starting with a seed kmer and associated reads.
    First find the reads containing the kmer_seq value, iterate through reads and
    either create a new contig or add to existing contig.
    Args:
        kmer_seq: String of kmer sequence.
        fq_recs: Dictionary with sequence values as keys and a list of fq_read objects.
        kmer_len: Integer of kmer size.
        kmer_tracker: KmerTracker object that contains all the kmer values.
        contig_buffer: ContigBuffer object to track the buffered contig objects.
    Return: None
    """

    logger = logging.getLogger('breakmer.assembly.assembler')
    contig = None

    # Find all reads with kmer sequence passed in.
    kmer_reads = find_reads(kmer_seq, fq_recs.items(), set())
    contig_buffer.add_used_mer(kmer_seq)

    kmer_values = {'seq': kmer_seq,
                   'counts': kmer_tracker.get_count(kmer_seq),
                   'kmer_set': kmer_tracker.kmer_set,
                   'len': kmer_len}
    for read_vals in kmer_reads:
        read, kmer_pos, bool, rlen, nreads = read_vals
        read_align_values = {'read': read,
                             'align_pos': kmer_pos,
                             'nreads': nreads}
        contig_buffer.add_used_read(read.id)
        # If no contig, build one.
        if not contig:
            contig = contig_assembler.Contig(kmer_values, read_align_values)
            contig_buffer.add_contig(read, contig)
        # Check if read should be added to the existing contig.
        else:
            contig.check_read(kmer_values, read_align_values, 'setup')
    if contig:
        contig.finalize(fq_recs, kmer_tracker, contig_buffer, 'setup')


def find_reads(kmer_seq, read_items, used_reads, order='for'):
    """Return a list of tuples containing information from reads with the kmer sequence.
    First search all the read sequences for the given kmer sequence. Then,
    filter out used reads and order them according to position of the kmer
    sequence in the read sequence.
    Args:
        kmer_seq: String of kmer sequence.
        read_items: List of fq_recs (key, value) tuples.
        used_reads: Set of read IDs that have been previously used.
        order: String indicating how the list of the identified reads
               should be ordered.
    Return:
        kmer_reads: List of tuples containing:
                    1. read object,
                    2. start position of kmer match in read seq
                    3. Boolean that a match was found.
                    4. Length of the read sequence.
                    5. Number of reads with this sequence.
    """

    kmer_reads = []
    mapped_reads = filter(lambda x: x[2], map(read_search, [kmer_seq]*len(read_items), read_items))
    mapped_read_ids = map(lambda x: x[0].id, mapped_reads)
    filter_ids = set(mapped_read_ids) - set(used_reads)
    matched_reads = filter(lambda x: (x[0].id in filter_ids), mapped_read_ids)
    if order == 'rev':
        kmer_reads = sorted(matched_reads, key=lambda z: (-z[1], -z[3]))
    else:
        kmer_reads = sorted(matched_reads, key=lambda z: (z[1], -z[3]))
    return kmer_reads


def read_search(kmer_seq, read_items):
    """Return a tuple containing information regarding the alignment of the kmer_seq
    in a sequence read.
    This uses regex searching function re.search to determine if the kmer_seq
    is contained in the read sequence. If so, then it returns a 5 element
    tuple about information regarding this alignment. If no match, then return
    a 3 element tuple with None values.
    Args:
        kmer_seq: String of kmer sequence.
        read_items: List of fq_recs (key, value) tuples.
    Return:
        search_result: Tuple of result information.
                       1. read object,
                       2. start position of kmer match in read seq
                       3. Boolean that a match was found.
                       4. Length of the read sequence.
                       5. Number of reads with this sequence.
    """

    search_result = (None, None, None)
    seq, reads = read_values
    x = re.search(kmer_seq, seq)
    if x:
        search_result = (reads[0], x.start(), True, len(reads[0].seq), len(reads))
    return search_result


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
            read: fq_read object
            contig: Contig object.
        Return: None
        """

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

        self.used_mers.add(kmer_seq)

    def remove_kmers(self, kmer_tracker):
        """Remove used kmer sequences from the kmer tracking object and reset used
        kmer set.
        Args:
            kmer_tracker: KmerTracker object.
        Return: None
        """

        map(kmer_tracker.remove_kmer, list(self.used_kmers))
        self.used_kmers = set()

    def remove_reads(self, fq_reads):
        """Remove the used reads from the fq_reads dictionary.
        Args:
            fq_reads: Dictionary of fq_reads.
        Return: None
        """

        del_used = filter(lambda x: x in fq_reads, list(self.used_reads))
        map(fq_reads.__delitem__, del_used)
        self.used_reads = set()


class KmerTracker:
    """Wrapper class for storing the kmer objects. Useful for adding
    and extracting kmers.
    Attributes:
        kmers: List of tuples containing kmer count, kmer, kmer object
        ordered_kmers: OrderedDict object containing kmer seq as key and kmer count as value.
                       The top values are the most frequence kmer values.
        kmer_seqs: Set of kmer seq values that exist in ordered_kmers.
    """

    def __init__(self):
        self.kmers = []
        self.ordered_kmers = OrderedDict()
        self.kmer_seqs = set()

    def add_kmer(self, mer, count):
        """Add a kmer object to the list.
        Stores a tuple with kmer count and kmer sequence string.
        This allows easy sorting.
        Args:
            mer: String kmer sequence value.
            count: Integer of number of reads kmer is within.
        Return:
            None
        """
        if len(set(mer)) > 1:
            self.kmers.append((int(count), mer))

    def set_all_kmer_values(self):
        """Sort the kmer list by number of reads (descending) they are in first and then
        by sequence value and store them in an ordered dictionary.
        Args: None
        Return: None
        """

        kmers_sorted = sorted(self.kmers, key=lambda x: (int(x[0]), x[1]), reverse=True)
        for kmer in kmers_sorted:
            self.ordered_kmers[kmer[1]] = kmer[0]

    def has_mers(self):
        """Check if there are any kmers left in the dictionary.
        Args:
            None
        Return:
            True if there are items in the dictionary and the counts of those items are > 1.
            False if there are no items in the dictionary or the counts of those items are <= 1.
        """

        if len(self.ordered_kmers) > 0 and max(self.ordered_kmers.values()) > 1:
            return True
        else:
            return False

    def update_kmer_set(self):
        """Update the set of kmer values. The ordered_kmers dictionary
        dynamically changes as kmers are taken out.
        Args: None
        Return: None
        """

        self.kmer_set = set(self.ordered_kmers.keys())

    def get_kmer(self):
        """Return the first kmer in the ordered dictionary"""
        return self.ordered_kmer.items()[0]

    def get_count(self, kmer_seq):
        """Return the number of reads the kmer_seq is within."""
        return self.ordered_kmers[kmer_seq]

    def remove_kmer(self, kmer_seq):
        """Delete the record associated with kmer sequence."""
        del self.ordered_kmers[kmer_seq]