#! /usr/bin/python
# -*- coding: utf-8 -*-

"""bam_handler.py module

This module contains the classes and functions to handle the
"""

import pysam

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


def trim_qual(read, min_qual, min_len):
    qual_str = read.qual
    q = []
    coords = [0, len(qual_str)]
    start = seq_trim(qual_str, min_qual)
    if start == len(qual_str):
        return None
    else:
        end = len(qual_str) - seq_trim(qual_str[::-1], min_qual)
        lngth = end - start
        if lngth < min_len:
            return None
        nseq = read.seq[start:end]
        nqual = qual_str[start:end]
        read.seq = nseq
        read.qual = nqual
    return read


def fq_line(read, indel_only, min_len, trim=True):
    add_val = '0'
    if indel_only:
        add_val = '1'
    lineout = None
    if trim:
        read = trim_qual(read, 5, min_len)
    if read:
        lineout = "@" + get_seq_readname(read) + "_" + add_val + "\n" + read.seq + "\n+\n" + read.qual + "\n"
    return lineout


def get_seq_readname(read):
    """ """
    end = '1'
    if read.is_read2:
        end = '2'
    return read.qname + "/" + end


def check_pair_overlap(mate_seq, read, coords, trim_dir):
    """ """
    nmisses = 0
    add_clip = True
    clip_seq = read.seq[coords[0]:coords[1]]
    clip_len = coords[1] - coords[0]

    if abs(read.isize) < len(read.seq):
        if abs(len(read.seq) - (abs(read.isize) + 1)) >= clip_len:
            add_clip = False
    else:
        while check_overlap(trim_dir, mate_seq, clip_seq) and nmisses < 5 and len(clip_seq) > 0:
            if trim_dir == 'back':
                clip_seq = clip_seq[0:(len(clip_seq) - 1)]
            else:
                clip_seq = clip_seq[1:len(clip_seq)]
            nmisses += 1
        if len(clip_seq) == 0 or nmisses == 5:
            add_clip = True
        else:
            add_clip = False
    return add_clip


def check_overlap(dir, mate_seq, clip_seq):
    """ """
    if dir == 'back':
        return mate_seq.find(clip_seq) != (len(mate_seq) - len(clip_seq))
    else:
        return mate_seq.find(clip_seq) != 0


def get_clip_coords(read):
    """This will parse a cigar string for a read and determine the coordinates
    of the read that are not softclipped by the aligner.

    Read cigar is a list of tuples [(4,5),(0,80),(4,15)] 5 bp clipped in the start, 80 bp matching, 15 bp clipped at the end
    Start: coords = [0,0]
    Iter 1: coords = [5,5]
    Iter 2: coords = [5,85]
    Iter 3: coords = [5,85]

    Args:
        read:        pysam read object.
    Return:
        clip_coords: List with two integer values indicating the coordinates of
                     the sequence read that are not clipped.
    """

    clip_coords = [0, len(read.qual)]
    # First value is start index, second value is end index.
    coords = [0, 0]
    for i in range(len(read.cigar)):
        code, clen = read.cigar[i]
        # Inc coords if not deletion or softclip
        if not code == 2 and not code == 4:
            coords[1] += clen
        # First value is softclip, increment both by clip amount.
        if code == 4: 
            if i == 0:
                coords[0] = clen
                coords[1] += clen
        clip_coords = coords
    return clip_coords


def seq_trim(qualStr, minQual):
    """Find the first position in a list of quality values that is above the minimum
    quality value input.
    Iterate over the list of quality values, starting at the first position, and
    return the position where the quality if greater than minQual.
    Args:
        qualStr: List of quality values from pysam read object (i.e., read.qual).
                  These are Phred-based and assumed to be offset by 33.
        minQual: Integer value of the minimum acceptable quality
    Return:
        counter: Integer representing the position in the list.
    """

    counter = 0
    while (ord(qualStr[counter]) - 33) < minQual:
        counter += 1
        if counter == len(qualStr):
            break
    return counter


def trim_coords(qualStr, minQual):
    """Searches quality values of a sequence read start->end and end->start to
    determine if there is a string of low quality sequences.

    Scan along the qualStr and continue while the quality is < minQual and
    return the index of the last low quality score in the string.

    qualStr = [1-1-1-2-2-2-2-20-20-20-30-30-30]
    seq_trim(qualStr, 3) will return 6 for the start and len(qualStr) for the end.

    Args:
        qualStr (list): List of quality values from pysam read object (i.e., read.qual).
                        These are Phred-based and assumed to be offset by 33.
        minQual (int):  Value of the minimum acceptable Phred quality score.
    Return:
        three element tuple:
            1. Position start where sequence quality is good (> minQual)
            2. Position end where sequence quality is good (> minQual)
            3. Length of the sequence that has good quality.
    """

    # Scan from the start of the qualStr and stop when the base qual > minQual
    start = seq_trim(qualStr, minQual)
    if start == len(qualStr):
        return (0, 0, 0)
    else:
        # Reverse qualStr and scan from the start and stop when the base qual > minQual
        end = len(qualStr) - seq_trim(qualStr[::-1], minQual)
        trimLength = end - start
        return (start, end, trimLength)


def pe_meta(read):
    """Checks if the read is from a proper paired-end mapping, assuming an Illumina
    library.

    If the read is mapped in a proper pair, check if it overlaps with its paired read.

       Args:
            read: pysam read object
        Return:
            proper_map: Boolean to indicate that the read-pair is properly mapped
            overlap_read: Boolean to indicate that the read-pair overlap (i.e.,
                          insert size < 2*read_len
    """

    properMap = False
    overlapReads = False
    if (((read.flag == 83 or read.flag == 147) and read.tlen < 0) or ((read.flag == 99 or read.flag == 163) and read.tlen > 0)):
        properMap = True
        if abs(read.tlen) < (2 * len(read.seq)):
            overlapReads = True
    return properMap, overlapReads


def get_region_reads(bamFile, chrom, start, end):
    """Open Bam file using pysam and fetch aligned reads in the
    specified region.

    Args:
        bamFile (str): Bam file full path, index must be in the same location
        chrom (str):    Chromosome name for region
        start (int):    Region's start position.
        end (int):      Region's end position.
    Return:
        reads (list): List containing pysam read objects
        bamF (pysam bam object): Open pysam bam file object.
    """

    bamF = pysam.Samfile(bamFile, 'rb')
    reads = bamF.fetch(chrom, start, end)
    return (reads, bamF)


def get_variant_reads(bamFile, chrom, start, end, insertSizeThresh):
    """Get the softclipped, discordant read pairs, and unmapped reads.
    These reads are stored in the VarReadTracker object.

    Iterate through all the reads in a region. Skip the duplicates and 
    qc failed reads. Store all the unmapped reads. All other reads pass
    to the check_read function.

    Args:
        bamFile (str):  Path to the bam file to open, must be indexed!
        chrom (str):    Chromosome of the region to extract
        start (int):    Region start location to extract.
        end (int):      Region end location to extract.

    Returns:
        varReadTracker (VariantReadTracker): VarReadTracker object
    """

    reads, bamF = get_region_reads(bamFile, chrom, start, end)
    varReadTracker = VariantReadTracker(bamF, insertSizeThresh)
    for read in reads:
        skip = False
        if read.mate_is_unmapped or read.rnext == -1:
            read.mate_is_unmapped = True
        if read.is_duplicate or read.is_qcfail:
            skip = True
        if read.is_unmapped:
            varReadTracker.add_unmapped_read(read)
            skip = True
        if skip:  # Skip a read if it is marked as a duplicate, qc failed, or unmapped
            continue
        varReadTracker.check_read(read)
    return varReadTracker


def get_strand_str(isReverseBoolean):
    strand = '+'
    if isReverseBoolean:
        strand = '-'
    return strand


def get_strand_key(read, ordered=False):
    strands = []
    readStrand = '+'
    if read.is_reverse:
        readStrand = '-'
    mateStrand = '+'
    if read.mate_is_reverse:
        mateStrand = '-'
    strands = [readStrand, mateStrand]

    if ordered:
            strands.reverse()
    return ':'.join(strands)


def cluster_regions(dReadLst, idx, clusterType):
    distBuffer = None
    clusterLst = []
    for dRead in dReadLst:
        if distBuffer is None:
            distBuffer = dRead.readLen
        # trgtStart = dRead.pos[0]
        # mateStart = dRead.pos[1]
        # print 'cluster_regions dRead', dRead.pos, dRead.readLen, clusterType
        if len(clusterLst) == 0:
            clusterLst.append([dRead.pos[idx], dRead.pos[idx] + dRead.readLen, [dRead.readInfoStr]])
            # print 'Initial cluster list', clusterLst
        else:
            # Check overlap
            add = False
            for i, c in enumerate(clusterLst):
                # print 'Checking read pos against cluster region', c, dRead.pos
                startWithin = dRead.pos[idx] >= c[0] and dRead.pos[idx] <= c[1]
                withinBuffer = dRead.pos[idx] > c[1] and dRead.pos[idx] - c[1] <= distBuffer
                # print 'in check', startWithin, withinBuffer
                if startWithin or withinBuffer:
                    readInfoLst = clusterLst[i][2]
                    readInfoLst.append(dRead.readInfoStr)
                    # print 'Add read to cluster region', clusterLst[i]
                    clusterLst[i] = [c[0], dRead.pos[idx] + dRead.readLen, readInfoLst]
                    add = True
            if not add:
                # print 'No add, creating new cluster region'
                clusterLst.append([dRead.pos[idx], dRead.pos[idx] + dRead.readLen, [dRead.readInfoStr]])
    return clusterLst


def get_cluster_membership(item, clusters, idx):
    for i, cluster in enumerate(clusters):
        # print cluster
        # print item.pos
        if item.pos[idx] >= cluster[0] and item.pos[idx] <= cluster[1]:
            return i


class discReadPair:
    def __init__(self, read, orderType):
        self.pos = []
        self.strands = []
        self.readName = read.qname
        self.readLen = read.rlen
        self.readInfoStr = ''
        # self.read = read
        self.set_values(read, orderType)

    def set_values(self, read, orderType):
        # print 'bam_handler.py set_values() for discReadPair', read.pos, read.mpos
        self.pos = [read.pos, read.mpos]
        self.strands = [get_strand_str(read.is_reverse), get_strand_str(read.mate_is_reverse)]
        if (orderType == 'ordered') and (read.mpos < read.pos):
            # Store the read and mate ordered by chrom alignment position
            self.pos.reverse()
            self.strands.reverse()
        self.readInfoStr = '|'.join([str(x) for x in [read.qname, self.strands[0], self.strands[1], read.tlen, read.mpos]])
        # print 'bam_hanlder.py set_values() readInfoStr', self.readInfoStr


class discReads:
    """
    """
    def __init__(self, insertSizeThresh):
        self.reads = {'inter': {}, 'intra': {}}
        self.insertSizeThresh = insertSizeThresh
        self.checkedIds = {}
        self.clusters = {}
        self.disc = {}

    def add_inter_discread(self, bam, read):
        # print 'bam_handler.py add_inter_discread()', read
        dRead = discReadPair(read, 'unordered')
        mateRefId = bam.getrname(read.rnext)
        if mateRefId not in self.reads['inter']:
            self.reads['inter'][mateRefId] = {}
        strandKey = get_strand_key(read)
        if strandKey not in self.reads['inter'][mateRefId]:
            self.reads['inter'][mateRefId][strandKey] = []
        self.reads['inter'][mateRefId][strandKey].append(dRead)
        # print 'bam_handler.py add_inter_discread() self.reads inter', mateRefId, strandKey, '\n'
        # for dRead in self.reads['inter'][mateRefId][strandKey]:
            # print '\t', dRead.readInfoStr

        if mateRefId not in self.disc:
            self.disc[mateRefId] = []
        self.disc[mateRefId].append((read.pos, read.mpos))
        # print 'bam_handler.py add_inter_discread() disc dictionary', mateRefId, self.disc[mateRefId]

    def add_intra_discread(self, read, overlapping_reads):
        discType = 'other'
        dRead = discReadPair(read, True)
        disc_ins_size = abs(read.tlen) >= self.insertSizeThresh
        strandKey = ''
        if (read.is_reverse and read.mate_is_reverse) or (not read.is_reverse and not read.mate_is_reverse):
            discType = 'inv'
            strandKey = get_strand_key(read)
        elif (read.is_reverse and not read.mate_is_reverse and read.pos < read.mpos) or (not read.is_reverse and read.mate_is_reverse and read.pos > read.mpos):
            discType = 'td'
            strandKey = '-:+'
        elif disc_ins_size:
            discType = 'dist'
            strandKey = get_strand_key(read, True)
        elif (read.is_reverse and not read.mate_is_reverse and read.pos < read.mpos) or (not read.is_reverse and read.mate_is_reverse and read.mpos < read.pos):
            discType = 'other'
            strandKey = get_strand_key(read, True)
        else:
            dRead = None

        if dRead is None:
            return

        if discType not in self.reads['intra']:
            self.reads['intra'][discType] = {}
        if strandKey not in self.reads['intra'][discType]:
            self.reads['intra'][discType][strandKey] = []
        self.reads['intra'][discType][strandKey].append(dRead)

        if read.tid not in self.disc:
            self.disc[read.tid] = []
        self.disc[read.tid].append((read.pos, read.mpos))

    def add_read_pair(self, bam, read, overlapping_reads):
        """
        Args:
            read:
        Return:
            None
        """
        if read.qname not in self.checkedIds:
            self.checkedIds[read.qname] = read.qname
        else:
            return

        if read.mapq == 0 or read.mate_is_unmapped:
            return

        # Extract read-pairs that are mapped to different chromosomes or fair apart.
        diff_chroms = read.rnext != -1 and read.tid != read.rnext
        if read.tid == read.rnext and not overlapping_reads:
            self.add_intra_discread(read, overlapping_reads)
        elif diff_chroms:
            # print 'bam_handler.py add_read_pair(), diff_chroms', diff_chroms, read.rnext, read.tid, read.rnext
            self.add_inter_discread(bam, read)

    def cluster_discreads(self):
        """self.reads is a dictionary with 3 levels
        1. Inter / intra
        2. Chrom (inter) / inv, td, dist, other (intra)
        3. -:+, -:-, +:+, +:-
        4. List of discRead objects
        """
        # print 'cluster_discreads()', '*'*25
        for key1 in self.reads:
            # print 'key1', key1
            d1 = self.reads[key1]
            for key2 in d1:
                # print 'key2', key2
                d2 = d1[key2]
                interClusterClusters = {}
                for key3 in d2:
                    # print 'key3', key3
                    dReadsLst = d2[key3]
                    # print 'read list', dReadsLst
                    srt1 = sorted(dReadsLst, key=lambda x: x.pos[0])
                    srt2 = sorted(dReadsLst, key=lambda x: x.pos[1])
                    c1 = cluster_regions(srt1, 0, 'target')
                    c2 = cluster_regions(srt2, 1, 'mate')
                    for item in dReadsLst:
                        # print 'Disc read pair obj', item.readInfoStr
                        cIdx1 = get_cluster_membership(item, c1, 0)
                        cIdx2 = get_cluster_membership(item, c2, 1)
                        regionPairKey = '|'.join([key1, key2, key3, str(cIdx1), str(cIdx2)])
                        # print 'regionPairKey', regionPairKey
                        leftBrkpt = c1[cIdx1][0]
                        rightBrkpt = c2[cIdx2][0]
                        leftStrand, rightStrand = key3.split(':')
                        if leftStrand == '+':
                            leftBrkpt = c1[cIdx1][1]
                        if rightStrand == '+':
                            rightBrkpt = c2[cIdx2][1]
                        if regionPairKey not in self.clusters:
                            self.clusters[regionPairKey] = {'readCount': 0,  # Read count for sub cluster, based on strand
                                                            'interClusterCount': 0,  # Read count for cluster ignoring strands, this will be used for interchrom clustering.
                                                            'leftBounds': c1[cIdx1][0:2],
                                                            'rightBounds': c2[cIdx2][0:2],
                                                            'leftBrkpt': leftBrkpt,
                                                            'rightBrkpt': rightBrkpt,
                                                            'clusterId': len(self.clusters) + 1}
                            if key1 == 'inter':
                                # print 'Inter check clustering', interClusterClusters
                                matchFound = False
                                for clusterKey in interClusterClusters:
                                    # print 'Checking clustering of inter clusters', clusterKey, self.clusters[clusterKey]['leftBrkpt'], regionPairKey, leftBrkpt
                                    # print 'Checking clustering of inter clusters', clusterKey, self.clusters[clusterKey]['rightBrkpt'], regionPairKey, rightBrkpt
                                    if (abs(self.clusters[clusterKey]['leftBrkpt'] - leftBrkpt) < 1000) and (abs(self.clusters[clusterKey]['rightBrkpt'] - rightBrkpt) < 1000):
                                        # Merge the clusters
                                        interClusterClusters[clusterKey].append(regionPairKey)
                                        matchFound = True
                                        break
                                if not matchFound:
                                    # print 'No match', regionPairKey
                                    interClusterClusters[regionPairKey] = [regionPairKey]
                        self.clusters[regionPairKey]['readCount'] += 1
                        self.clusters[regionPairKey]['interClusterCount'] += 1
                if len(interClusterClusters) > 0:
                    for clusterKey in interClusterClusters:
                        totalCounts = 0
                        for cKey in interClusterClusters[clusterKey]:
                            totalCounts += self.clusters[cKey]['readCount']
                        for cKey in interClusterClusters[clusterKey]:
                            self.clusters[cKey]['interClusterCount'] = totalCounts
                            self.clusters[cKey]['clusterId'] = self.clusters[clusterKey]['clusterId']
        # print 'Complete clusters', self.clusters
        return self.clusters

    def check_inv_readcounts(self, brkpts):
        """ """
        brkpt1 = min(brkpts)
        brkpt2 = max(brkpts)
        counts = 0
        bpBuffer = 50
        # print 'Inversion reads', self.reads['intra']['inv']
        # print 'Brkpts', brkpts
        if 'inv' not in self.reads['intra']:
            return counts
        for strand in self.reads['intra']['inv']:
            lStrand, rStrand = strand.split(':')
            strandReads = self.reads['intra']['inv'][strand]
            for dRead in strandReads:
                # print strand, dRead.pos
                if lStrand == '+' and rStrand == '+':
                    if (dRead.pos[0] <= (brkpt1 + bpBuffer)) and (dRead.pos[1] <= (brkpt2 + bpBuffer) and dRead.pos[1] >= (brkpt1 - bpBuffer)):
                        counts += 1
                else:
                    # print dRead.pos, brkpt1, brkpt2
                    if (dRead.pos[0] <= (brkpt2 + bpBuffer) and dRead.pos[0] >= (brkpt1 - bpBuffer)) and dRead.pos[1] >= (brkpt2 - bpBuffer):
                        counts += 1
        # print 'Counts', counts
        return counts

    def check_td_readcounts(self, brkpts):
        """ """
        brkpt1 = min(brkpts)
        brkpt2 = max(brkpts)
        counts = 0
        bpBuffer = 50
        if 'td' not in self.reads['intra']:
            return counts
        for dRead in self.reads['intra']['td']['-:+']:
            if (dRead.pos[0] >= (brkpt1 - bpBuffer) and dRead.pos[0] <= (brkpt2 + bpBuffer)) and (dRead.pos[1] <= (brkpt2 + bpBuffer) and dRead.pos[1] >= (brkpt1 - bpBuffer)):
                counts += 1
        return counts

    def check_other_readcounts(self, brkpts):
        """ """
        counts = [0] * len(brkpts)
        for i in range(len(brkpts)):
            b = brkpts[i]
            if 'other' not in self.reads['intra']:
                return max(counts)
            for strand in self.reads['intra']['other']:
                lStrand, rStrand = strand.split(':')
                strandReads = self.reads['intra']['other'][strand]
                for dRead in strandReads:
                    if abs(dRead.pos[0] - b) <= 300 or abs(dRead.pos[1] - b) <= 300:
                        counts[i] += 1
        return max(counts)

    def check_inter_readcounts(self, targetBrkptChr, targetBrkptBp, nonTargetBrkpts):
        """ """
        # counts = [0] * len(brkpts)
        # for i in range(len(brkpts)):
        #     b = brkpts[i]
        #     if 'other' not in self.reads['intra']:
        #         break
        #     for strand in self.reads['intra']['other']:
        #         lStrand, rStrand = strand.split(':')
        #         strandReads = self.reads['intra']['other'][strand]
        #         for dRead in strandReads:
        #             if abs(dRead.pos[0] - b) <= 300 or abs(dRead.pos[1] - b) <= 300:
        #                 counts[i] += 1
        # return max(counts)
        discReadCount = 0
        # print 'sv_caller.py get_disc_read_count', targetBrkptChr, targetBrkptBp
        # print 'Read storage dict', self.reads['inter']
        for otherBrkpts in nonTargetBrkpts:
            nonTargetBrkptChr = otherBrkpts[0].replace('chr', '')
            nonTargetBrkptBps = otherBrkpts[1:]
            # print 'Non-target brkpts', nonTargetBrkptChr, nonTargetBrkptBps
            for nonTargetBrkptBp in nonTargetBrkptBps:
                # print 'non-target brkpt', nonTargetBrkptBp
                if nonTargetBrkptChr in self.reads['inter']:
                    for strand in self.reads['inter'][nonTargetBrkptChr]:
                        for discReadPair in self.reads['inter'][nonTargetBrkptChr][strand]:
                            d1 = abs(targetBrkptBp - discReadPair.pos[0])
                            d2 = abs(nonTargetBrkptBp - discReadPair.pos[1])
                            # print 'distances', d1, d2
                            if d1 <= 1000 and d2 <= 1000:
                                discReadCount += 1
        return discReadCount


class VariantReadTracker:
    """A class to track the reads that are identified to be 'misaligned' to
    the reference sequence.

    Attributes:
        pair_indices (dict):  Dictionary of a dictionary tracking the index of paired
                              reads in the valid list.
        valid (list):         List of read objects that are valid to consider for extraction.
        disc (dict):          Dictionary of read IDs for read-pairs that are discordantly mapped.
        unmapped (dict):      Dictionary of unmapped reads with mapped mate in the region.
        unmapped_keep (list): List containing names of reads that are mapped but their mate is unmapped and wasn't
                              kept on the first pass.
        inv (list):           List of tuples, each containing read-pair information that have alignments
                              suggestive of an inversion event.
        td (list):            List of tuples, each containing read-pair information that have alignments
                              suggestive of a tandem dup event.
        other (list):         List of tuples, each containing read-pair information that have alignments
                              suggestive of some uncategorized event.
        sv (dict):            Dictionary
        bam (str):            Bam file source the reads came from.
    """

    def __init__(self, bamFile, insertSizeThresh):
        """
        """

        self.pair_indices = {}
        self.valid = []
        self.discReadTracker = discReads(insertSizeThresh)
        self.unmapped = {}
        self.unmapped_keep = []
        self.sv = {}
        self.bam = bamFile

    def check_read(self, read):
        """Stores all reads in the self.pair_indices dictionary if it is
        mapped.

        Check if the read is part of a discordantly mapped read pair.

        Check if the read is properly mapped, as indicated by bam encoding, and
        whether the read overlaps with its pair.

        self.valid = [(read, proper_map, overlapping_reads), (read, proper_map, overlapping_reads), ...]
        self.pair_indices[read.qname][1 (read1)/0 (read2)] = index of read in self.valid

        Args:
            read (pysam read obj): An aligned sequence read.

        Returns:
            None
        """

        proper_map, overlapping_reads = pe_meta(read)
        if read.qname not in self.pair_indices and not read.mate_is_unmapped:
            self.discReadTracker.add_read_pair(self.bam, read, overlapping_reads)

        self.valid.append((read, proper_map, overlapping_reads))
        if read.qname not in self.pair_indices and not read.mate_is_unmapped:
            self.pair_indices[read.qname] = {}
        if read.qname in self.pair_indices:
            self.pair_indices[read.qname][int(read.is_read1)] = len(self.valid) - 1

    def add_unmapped_read(self, read):
        """Add read to unmapped dictionary with name as the key, object as the value.

        Args:
            read (pysam read obj): pysam read object.

        Returns:
            None
        """

        self.unmapped[read.qname] = read

    def check_clippings(self, kmer_size, region_start_pos, region_end_pos):
        """
        """

        for read_vals in self.valid:
            read, proper_map, overlap_reads = read_vals
            if read.cigar or len(read.cigar) > 1:
                good_qual_coords = trim_coords(read.qual, 3)  # Get the (start, end, length) of the high-quality sequence bases.
                clip_coords = get_clip_coords(read)  # Get the [start, end] of the non-clipped sequence bases.
                self.extract_clippings(read_vals, clip_coords, good_qual_coords, kmer_size)

            if (read.pos >= region_start_pos and read.pos <= region_end_pos) and read.mapq > 0 and read.mate_is_unmapped:
                self.unmapped_keep.append(read.qname)

    def extract_clippings(self, read_vals, clip_coords, good_qual_coords, kmer_size):
        """
        """

        read, proper_map, overlap_reads = read_vals
        clip_seqs = {'clipped': [], 'buffered': []}

        if clip_coords[0] <= good_qual_coords[0] and clip_coords[1] >= good_qual_coords[1]:
            return

        new_clip_coords = [0, 0]
        add_clip = [False, False]
        indel_only = False
        start_clip = clip_coords[0] > 0
        end_clip = clip_coords[1] < len(read.qual)
        if start_clip and end_clip:
            add_clip = [True, True]
        else:
            if start_clip:
                add_clip[0] = True
                new_clip_coords = [0, clip_coords[0]]
                if overlap_reads and read.is_reverse:
                    mate_seq = self.valid[self.pair_indices[read.qname][int(read.is_read1)]][0].seq
                    add_clip[0] = check_pair_overlap(mate_seq, read, [0, clip_coords[0]], 'back')
                if proper_map:
                    if read.is_reverse:
                        indel_only = True
                    else:
                        indel_only = False
            elif end_clip:
                new_clip_coords = [clip_coords[1], len(read.seq)]
                add_clip[1] = True
                if overlap_reads and not read.is_reverse:
                    mate_seq = self.valid[self.pair_indices[read.qname][int(read.is_read1)]][0].seq
                    add_clip[1] = check_pair_overlap(mate_seq, read, [clip_coords[1], len(read.seq)], 'front')
                if proper_map:
                    if read.is_reverse:
                        indel_only = indel_only and False
                    else:
                        indel_only = indel_only and True
        final_add = add_clip[0] or add_clip[1]
        if add_clip[0]:
            clip_seqs['buffered'].append(read.seq[0:(clip_coords[0] + kmer_size)])
            clip_seqs['clipped'].append(read.seq[0:clip_coords[0]])
        if add_clip[1]:
            clip_seqs['buffered'].append(read.seq[(clip_coords[1] - kmer_size):len(read.seq)])
            clip_seqs['clipped'].append(read.seq[clip_coords[1]:len(read.seq)])
        if final_add:
            self.sv[get_seq_readname(read)] = (read, clip_seqs, new_clip_coords, indel_only)

    def write_seqs(self, clipped_fa, reads_fq, sv_bam, kmer_size):
        """
        """

        for name in self.unmapped_keep:
            if name in self.unmapped:
                read = self.unmapped[name]
                self.sv[get_seq_readname(read)] = (read, None, None, False)
                lout = ">" + read.qname + "\n" + str(read.seq)
                clipped_fa.write(lout + "\n")

        for name in self.sv:
            read, clip_seqs, clip_coords, indel_only = self.sv[name]
            if sv_bam:
                sv_bam.write(read)
            lout = fq_line(read, indel_only, kmer_size, True)
            if lout:
                reads_fq.write(lout)
            if clip_seqs:
                for clip in clip_seqs['buffered']:
                    clipped_fa.write(">" + name + "\n" + clip + "\n")
        self.bam.close()

    def clear_sv_reads(self):
        """
        """

        self.sv = None

    def get_disc_reads(self):
        """This function needs to be updated to handle the new disc read storage.
        """

        return self.discReadTracker.disc

    def cluster_discreads(self):
        """
        """

        dReadClusters = self.discReadTracker.cluster_discreads()
        return dReadClusters

    def check_inv_readcounts(self, brkpts):
        """
        """

        return self.discReadTracker.check_inv_readcounts(brkpts)

    def check_td_readcounts(self, brkpts):
        """ """
        return self.discReadTracker.check_td_readcounts(brkpts)

    def check_other_readcounts(self, brkpts):
        """ """
        return self.discReadTracker.check_other_readcounts(brkpts)

    def check_inter_readcounts(self, targetChr, targetBps, nonTargetBrkpts):
        """ """
        return self.discReadTracker.check_inter_readcounts(targetChr, targetBps, nonTargetBrkpts)
