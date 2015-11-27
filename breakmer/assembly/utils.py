#! /usr/bin/local/python
# -*- coding: utf-8 -*-

import re
import logging
from collections import OrderedDict

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


# class Kmer:
#     """Class to track values associated with a particular kmer sequence.
#     Attributes:
#         seq (str):
#         counts ():
#         kmerSeqSet ():
#         kmerLen ():
#     """

#     def __init__(self, seq, counts, kmerSeqSet, kmerLen):
#         self.seq = seq
#         self.counts = counts
#         self.kmerSeqSet = kmerSeqSet
#         self.kmerLen = kmerLen


def find_reads(kmerSeq, readItems, usedReads, order='for'):
    """Return a list of tuples containing information from reads with the kmer sequence.
    First search all the read sequences for the given kmer sequence. Then,
    filter out used reads and order them according to position of the kmer
    sequence in the read sequence.

    Args:
        kmerSeq (str):    Kmer sequence.
        readItems (list): List of dictionary items. These items will containing
                          be a tuple of (read sequence, list of fq_read objects with read sequence).
        usedReads (set):  Previously used read IDs.
        order (str):      Indicates how the list of the identified reads
                          should be ordered.
    Returns:
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
    # Remove the used read IDs.
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
        kmerSeq (str):    Kmer sequence.
        readItems (list): List of dictionary items. These items will containing
                          be a tuple of (read sequence, list of fq_read objects with read sequence).

    Returns:
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
