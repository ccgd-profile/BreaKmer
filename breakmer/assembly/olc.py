#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"

match_award = 1
mismatch_penalty = -2
gap_penalty = -2


# Creates empty matrix with all zeros
def zeros(shape):
    retval = []
    for x in range(shape[0]):
        retval.append([])
        for y in range(shape[1]):
            retval[-1].append(0)
    return retval


# No substituition matrix, just simple linear gap penalty model
def match_score(alpha, beta):
    if alpha == beta:
        return match_award
    elif alpha == '-' or beta == '-':
        return gap_penalty
    else:
        return mismatch_penalty


def nw(seq1, seq2):
    global max_score
    # lengths of two sequences
    m = len(seq1)
    n = len(seq2)

    # Generate DP (Dynamic Programmed) table and traceback path pointer matrix
    # the DP table
    score = zeros((n + 1, m + 1))
    for i in range(n + 1):
        score[i][0] = 0
    for j in range(m + 1):
        score[0][j] = 0

    # Traceback matrix
    # to store the traceback path
    pointer = zeros((n + 1, m + 1))
    for i in range(n + 1):
        pointer[i][0] = 1
    for j in range(m + 1):
        pointer[0][j] = 2

    # Calculate DP table and mark pointers
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score_diagonal = score[i - 1][j - 1] + match_score(seq1[j - 1], seq2[i - 1])
            score_up = score[i][j - 1] + gap_penalty
            score_left = score[i - 1][j] + gap_penalty
            score[i][j] = max(score_left, score_up, score_diagonal)

            if score[i][j] == score_diagonal:
                # 3 means trace diagonal
                pointer[i][j] = 3
            elif score[i][j] == score_up:
                #  2 means trace left
                pointer[i][j] = 2
            elif score[i][j] == score_left:
                # 1 means trace up
                pointer[i][j] = 1

    # Finding the right-most match which represents a longest overlap
    # Note that .index() will find the index of the first item in the list that matches,
    # so if you had several identical "max" values, the index returned would be the one for the first.
    max_i = -200
    for ii in range(n + 1):
        if score[ii][-1] >= max_i:
            max_i = score[ii][-1]
            i = ii

    prei = i
    prej = j
    # Traceback, follow pointers in the traceback matrix
    # initial sequences
    align1, align2 = '', ''
    while 1:
        if pointer[i][j] == 3:
            align1 = seq1[j - 1] + align1
            align2 = seq2[i - 1] + align2
            i -= 1
            j -= 1
        elif pointer[i][j] == 2:
            # 2 means trace left
            align2 = '-' + align2
            align1 = seq1[j - 1] + align1
            j -= 1
        elif pointer[i][j] == 1:
            # 1 means trace up
            align2 = seq2[i - 1] + align2
            align1 = '-' + align1
            i -= 1
        if (i == 0 or j == 0):
            break

    return (align1, align2, prej, j, prei, i, max_i)


class Align:
    """
    """
    def __init__(self, seq1, seq2):
        self.seq1 = se1
        self.seq2 = seq2
        self.align1 = None
        self.align2 = None
        self.prej = None
        self.j = None
        self.prei = None
        self.i = None
        self.max = None
        self.ident = None
        self.align()

    def align(self):
        self.align1, self.align2, self.prej, self.j, self.prei, self.i, self.max = nw(self.seq1, self.seq2)
        self.ident = round(float(self.max) / float(self.prej - self.j), 2)


class AlignManager:
    """
    """
    def __init__(self, seq1, seq2, scoreThresh, identThresh):
        self.seq1 = seq1
        self.seq2 = seq2
        self.scoreThresh = scoreThresh
        self.identThresh = identThresh
        self.aligns = None
        self.align_seqs()

    def align_seqs(self):
        """
        """
        self.aligns = (Align(self.seq1, self.seq2), Align(self.seq2, self.seq1))

    def check_align_thresholds(self):
        align1Check = self.aligns[0].max < self.scoreThresh or self.aligns[0].ident < self.identThresh
        align2Check = self.aligns[1].max < self.scoreThresh or self.aligns[1].ident < self.identThresh
        return align1Check and align2Check

    def same_seqs(self):
        hitEnds = self.aligns[0].j == 0 and self.aligns[0].i == 0
        equalLens = len(self.seq1) == len(self.seq2)
        return self.same_max_scores() and hitEnds and equalLens

    def same_max_scores(self):
        return self.aligns[0].max == self.aligns[1].max

    def read_is_superseq(self):
        """Check if seq2 is a superseq of seq1
        This function should be called based on the assumption that
        seq1 is checked against seq2.
        Args: None
        Return: Boolean of check
        """
        return len(self.seq1) < len(self.seq2) or (self.aligns[0].prej == len(self.seq1) and self.aligns[0].j == 0)

    def read_is_subseq(self):
        """Check if seq2 is a subseq of seq1
        Args: None
        Return: Boolean of check
        """
        return len(self.seq1) > len(self.seq2) or (self.aligns[1].prej == len(self.seq2) and self.aligns[1].j == 0)

    def better_align(self):
        """Check if align1 or align2 has a better score.
        """
        return self.aligns[0].max > self.aligns[1].max

    def get_alignment(self, index):
        """Return align object"""
        return self.aligns[index]

    def get_alignment_values(self, index, value):
        """Return a specific value from the alignment results
        Args:
            index: Integer index of the alignment - 0,1
            value: String value for the alignment value - i, prei
        """
        returnVal = None
        alignment = self.aligns[index]
        if value == 'i':
            returnVal = alignment.i
        elif value == 'prei':
            returnVal = alignment.prei
        return returnVal


    def get_kmer_align_indices(self, align_index, kmer_seq):
        """Return the alignment index of the kmer sequence with the
        sequences.
        """
        return (self.aligns[align_index].align1.replace('-', '').find(kmer_seq), self.aligns[align_index].align2.replace('-', '').find(kmer_seq))
