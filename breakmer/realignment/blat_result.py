#! /usr/bin/python
# -*- coding: utf-8 -*-

import math
import breakmer.utils as utils

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


class Matches:
    """
    """
    def __init__(self, values):
        self.match = int(values[0])
        self.mismatch = int(values[1])
        self.repeat = int(values[2])

    def get_total_matches(self):
        """Sum all match values"""
        return self.match + self.mismatch + self.repeat

    def get_total_matching(self):
        """Sum all values where the is match"""
        return self.match + self.repeat

    def get_mismatches(self):
        return self.mismatch

    def get_nmatches(self, matchType):
        nmatches = self.match
        if matchType == "mismatch":
            nmatches = self.mismatch
        elif matchType == 'repeat':
            nmatches = self.repeat
        return int(nmatches)


class Gaps:
    """Track the number and the size of gaps in the realignment
    Attributes:
        ref:    List of integers for gaps in reference alignment [number of gaps, total bp of gaps]
        query:  List of integers for gaps in query alignment [number of gaps in query sequence, total bp of gaps]
    """
    def __init__(self, values):
        """
        """
        self.ref = [int(values[6]), int(values[7])]
        self.query = [int(values[4]), int(values[5])]

    def get_ngaps(self, alignType):
        """Return a gaps"""
        ngaps = self.ref[0]
        if alignType == 'query':
            ngaps = self.query[0]
        return int(ngaps)

    def get_total_size(self):
        """Return the total bp of gap in the realignment"""
        return self.ref[1] + self.query[1]

    def get_gap_sizes(self):
        """ """
        return (self.ref[1], self.query[1])

    def get_total_num_gaps(self):
        """Return the total number of gaps in the realignment"""
        return self.ref[0] + self.query[0]


class AlignFragments:
    """
    """
    def __init__(self, values):
        self.blockSizes = [int(x) for x in values[18].rstrip(",").split(",")]
        self.ref = []
        self.query = []
        self.count = len(self.blockSizes)
        self.set_values(values)

    def set_values(self, values):
        print 'Set fragment values', values
        refStarts = [int(x) for x in values[20].rstrip(",").split(",")]
        queryStarts = [int(x) for x in values[19].rstrip(",").split(",")]
        for qstart, tstart, blocksize in zip(queryStarts, refStarts, self.blockSizes):
            self.ref.append((tstart, tstart + blocksize))
            self.query.append((qstart, qstart + blocksize))


class SVBreakpoint:
    """
    """
    def __init__(self, coords, svType, targetKey):
        self.chrom = coords[0]
        self.genomicCoords = coords[1:]
        self.svType = svType
        self.targetKey = targetKey
        self.annotated_trxs = {}

    def store_annotation(self, trxList, distList, coordIdx):
        """ """
        self.annotated_trxs[coordIdx] = [trxList, distList]


class Breakpoints:
    """
    """
    def __init__(self):
        self.contigBreakpoints = []
        self.genomicBreakpoints = []
        self.svBreakpoints = []

    def add_brkpts(self, which, bp):
        """ """
        if which == 'contig':
            self.contigBreakpoints.append(bp)
        else:
            # Tuple of coordinates for indel breakpoints with chromosome number as 'chr#'
            # - Deletion = (chr, bp1, bp2)
            # - Insertion = (chr, bp1)
            self.genomicBreakpoints.append(bp)

    def reverse_breakpts(self, querySeqSize):
        for i in range(len(self.contigBreakpoints)):
            contigBrkpts = self.contigBreakpoints[i]
            contigBrkpts[0] = querySeqSize - contigBrkpts[0]
            if len(contigBrkpts) > 1:
                contigBrkpts[1] = querySeqSize - contigBrkpts[1]
            self.contigBreakpoints[i] = contigBrkpts

    def set_sv_brkpt(self, coords, svType, targetKey):
        """ """
        self.svBreakpoints.append(SVBreakpoint(coords, svType, targetKey))
        # print self.svBreakpoints


class AlignValues:
    """
    """
    def __init__(self, values):
        self.ref = {'seqName': values[13],
                    'seqSize': int(values[14]),
                    'alignCoords': [int(values[15]), int(values[16])]}
        self.query = {'seqName': values[9],
                      'seqSize': int(values[10]),
                      'alignCoords': [int(values[11]), int(values[12])]}

    def get_coords(self, alignType, index=None):
        """Return the coordinate of the alignment for either the reference
        or the query sequence.
        Args:
            alignType:  String indicating the reference or query
            index:      Integer indicating the start(0) or end(1)
        """

        coords = self.ref['alignCoords']
        if alignType == 'query':
            coords = self.query['alignCoords']
        if index is not None:
            coords = coords[index]
        return coords

    def get_seq_name(self, alignType):
        """
        """
        name = str(self.ref['seqName'])
        if alignType == 'query':
            name = str(self.query['seqName'])
        return name

    def get_seq_size(self, alignType):
        """
        """
        size = str(self.ref['seqSize'])
        if alignType == 'query':
            size = str(self.query['seqSize'])
        return size


class BlatResult:
    """
    """
    def __init__(self, blatResultValues, refName, offset):
        self.loggingName = 'breakmer.realignment.blat_result'
        self.values = self.set_values(blatResultValues, refName, offset)
        self.matches = Matches(self.values)
        self.gaps = Gaps(self.values)
        self.alignVals = AlignValues(self.values)
        self.fragments = AlignFragments(self.values)
        self.strand = blatResultValues[8]
        self.breakpts = Breakpoints()
        # Sort results based on alignScore, percentIdent, number of gaps
        self.perc_ident = 100.0 - self.calcMilliBad()
        self.alignScore = self.get_nmatch_total() + (float(self.get_nmatch_total()) / float(self.get_seq_size('query')))
        self.ngaps = self.get_total_num_gaps()

        self.meanCov = 0.0
        self.seg_overlap = [0, 0]
        self.cigar = ''

        self.genes = ''
        self.in_target = False
        self.valid = True
        self.rep_man = None
        self.in_repeat = False
        self.repeat_overlap = 0.0
        self.repeat_coords = None
        self.filter_reps_edges = [False, False]
        self.interval = None

        self.indel_sizes = []
        self.indel_maxevent_size = [0, '']
        self.indel_flank_match = [0, 0]
        self.set_indel_locs()

    def set_values(self, blatResultValues, refName, offset):
        """Modify the blat values if refName and offset are not None
        Args:
            blatResultValues: List of values from a blat BlatResult
            refName:          String of chromosome AlignFragments
            offset:           Integer of genomic position for target alignment
        """
        rName = blatResultValues[13].replace('chr', '')
        if refName is not None:
            rName = refName
        blatResultValues[13] = rName

        coordOffset = 0
        if offset is not None:
            coordOffset = offset
        blatResultValues[15] = coordOffset + int(blatResultValues[15])
        blatResultValues[16] = coordOffset + int(blatResultValues[16])

        tstarts = [coordOffset + int(x) for x in blatResultValues[20].rstrip(",").split(",")]
        blatResultValues[20] = ",".join([str(x) for x in tstarts]) + ","
        return blatResultValues

    def set_sv_brkpt(self, coords, svType, targetKey):
        """ """
        self.breakpts.set_sv_brkpt(coords, svType, targetKey)

    def get_sv_brkpts(self):
        """ """
        return self.breakpts.svBreakpoints

    def calcMilliBad(self):
        badAlign = 0.0
        queryAlignSize = self.qend() - self.qstart()
        refAlignSize = self.tend() - self.tstart()
        minAlignSize = min(queryAlignSize, refAlignSize)
        if minAlignSize <= 0:
            return 0.0
        sizeDif = queryAlignSize - refAlignSize
        if sizeDif < 0:
            sizeDif = 0
        insertFactor = self.gaps.get_ngaps('query')
        totalMatches = self.matches.get_total_matches()
        if totalMatches != 0:
            badAlign = (1000 * (self.matches.get_mismatches() + insertFactor + round(3 * math.log(1 + sizeDif)))) / totalMatches
        return badAlign * 0.1

    def set_mean_cov(self, meanCov):
        self.meanCov = meanCov

    def set_segment_overlap(self, right, left):
        self.seg_overlap = [left, right]

    def qstart(self):
        """Query coordinate start"""
        return int(self.alignVals.get_coords('query', 0))

    def qend(self):
        """Query coordinate end"""
        return int(self.alignVals.get_coords('query', 1))

    def tstart(self):
        return self.alignVals.get_coords('reference', 0)

    def tend(self):
        return self.alignVals.get_coords('reference', 1)

    def get_seq_name(self, alignType):
        """Return the seq name of the query or reference. Typically the chromosome number"""
        return self.alignVals.get_seq_name(alignType)

    def get_seq_size(self, alignType):
        return int(self.alignVals.get_seq_size(alignType))

    def get_query_span(self):
        """Length of query sequence alignment to reference.
        """
        return self.qend() - self.qstart()

    def get_query_coverage(self):
        """Return percentage of query sequence realigned to reference"""
        return round((float(self.get_query_span()) / float(self.get_seq_size('query'))) * 100, 2)

    def spans_query(self):
        """Return boolean whether full query sequence is aligned"""
        return self.get_seq_size('query') == (self.qend() - self.qstart())

    def get_total_gap_size(self):
        return self.gaps.get_total_size()

    def get_total_num_gaps(self):
        return self.gaps.get_total_num_gaps()

    def get_nmatch_total(self):
        return self.matches.get_total_matching()

    def get_nmatches(self, matchType):
        return int(self.matches.get_nmatches(matchType))

    def sum_indel_flank_matches(self, flank_str):
        """ """
        m_indxs = []
        match_sum = 0
        for i in range(len(flank_str)):
            if flank_str[i] == "M":
                m_indxs.append(i)
        for indx in m_indxs:
            nmatch = ''
            windx = indx - 1
            while windx > -1 and utils.is_number(flank_str[windx]):
                nmatch = flank_str[windx] + nmatch
                windx = windx - 1
            match_sum += int(nmatch)
        return match_sum

    def set_indel_flank_matches(self):
        """ """
        if self.indel_maxevent_size[0] > 0:
            csplit = self.cigar.split(str(self.indel_maxevent_size[0]) + self.indel_maxevent_size[1])
            lflank = csplit[0]
            self.indel_flank_match[0] += self.sum_indel_flank_matches(lflank)
            rflank = csplit[-1]
            self.indel_flank_match[1] += self.sum_indel_flank_matches(rflank)

    def set_indel_locs(self):
        """ """
        chrom = 'chr' + self.get_seq_name('reference')
        for i in range(self.fragments.count - 1):
            if i == 0 and self.fragments.query[i][0] > 0:
                self.cigar = str(self.fragments.query[i][0]) + "S"
            qend1 = int(self.fragments.query[i][1])
            qstart2 = int(self.fragments.query[i + 1][0])
            tend1 = int(self.fragments.ref[i][1])
            tstart2 = int(self.fragments.ref[i + 1][0])
            ins_bp = qstart2 - qend1
            del_bp = tstart2 - tend1
            bp1 = tend1
            bp2 = tstart2
            self.cigar += str(self.fragments.blockSizes[i]) + "M"
            if ins_bp > 0:
                self.breakpts.add_brkpts('genomic', (chrom, bp1))
                self.indel_sizes.append("I" + str(ins_bp))
                self.breakpts.add_brkpts('contig', [qend1, qstart2])
                # self.breakpts.add_brkpts('contig', qstart2)
                self.cigar += str(ins_bp) + "I"
                if ins_bp > self.indel_maxevent_size[0]:
                    self.indel_maxevent_size = [ins_bp, "I"]
            if del_bp > 0:
                self.breakpts.add_brkpts('genomic', (chrom, bp1, bp2))
                self.indel_sizes.append("D" + str(del_bp))
                self.breakpts.add_brkpts('contig', [qend1])
                self.cigar += str(del_bp) + "D"
                if del_bp > self.indel_maxevent_size[0]:
                    self.indel_maxevent_size = [del_bp, "D"]

        self.cigar += str(self.fragments.blockSizes[-1]) + "M"
        # endClipped = self.get_seq_size('query') - self.qend()
        # if endClipped > 0:
        #     self.cigar += str(endClipped) + "S"

        self.set_indel_flank_matches()
        if self.strand == "-":
            self.breakpts.reverse_breakpts(self.get_seq_size('query'))

    # def add_query_brkpt(self, brkpt):
    #     """ """
    #     if brkpt not in self.query_brkpts:
    #         self.query_brkpts.append(brkpt)

    def get_genomic_brkpts(self):
        """ """
        return self.breakpts.genomicBreakpoints

    def get_brkpt_str(self, with_sizes=False):
        """ """
        brkpt_out = []
        bp_str = []
        chrm = 'chr' + str(self.get_name('hit'))
        if len(self.breakpts) > 0:
            for b, s in zip(self.breakpts, self.indel_sizes):
                if len(b) > 1:
                    bb = "-".join([str(x) for x in b])
                else:
                    bb = str(b[0])
                bstr = chrm + ":" + bb
                if with_sizes:
                    bstr += " " + "(" + s + ")"
                bp_str.append(bstr)
            brkpt_out.append(",".join(bp_str))
        return ",".join(brkpt_out)

    def get_brkpt_locs(self):
        brkpt_locs = []
        for b in self.breakpts:
            brkpt_locs.extend(b)
        return brkpt_locs

    def get_gene_anno(self):
        return self.genes

    def get_blat_output(self):
        return "\t".join([str(x) for x in self.values])

    def get_len(self):
        return self.qend - self.qstart

    def get_coords(self, alignType):
        """ """
        return self.alignVals.get_coords(alignType)

    # def set_repeats(self, target_rep_mask, all_rep_mask):
    #     self.rep_man = blat_repeat_manager()
    #     if self.matches['rep'] > 0:
    #         self.in_repeat = True
    #     if target_rep_mask and all_rep_mask:
    #         # Check rep_mask if it exists.
    #         rmask = target_rep_mask
    #         if not self.in_target:
    #             rmask = None
    #             if self.vals['hit']['name'] in all_rep_mask:
    #                 rmask = all_rep_mask[self.vals['hit']['name']]
    #         if rmask:
    #             self.rep_man.setup(self.get_coords('hit'), rmask)
    #             self.in_repeat, self.repeat_overlap, self.repeat_coords, self.filter_reps_edges = self.rep_man.other_values

    def in_target_region(self, targetRegionCoordinates):
        """ """
        refCoordStart, refCoordEnd = self.get_coords('ref')
        regionStart = targetRegionCoordinates[1] - targetRegionCoordinates[4]
        regionEnd = targetRegionCoordinates[2] + targetRegionCoordinates[4]
        start_in = refCoordStart >= regionStart and refCoordStart <= regionEnd
        end_in = refCoordEnd <= regionEnd and refCoordEnd >= regionStart
        if targetRegionCoordinates[0] == self.get_seq_name('reference') and (start_in or end_in):
            self.in_target = True
            self.genes = targetRegionCoordinates[3]

    def set_gene_annotations(self, targetRegionCoordinates, annotations):
        """
        Args:
            targetRegionCoordinates: Tuple containing (chr, start, end, name, intervals, regionBufferSize)
            annotations:             Object containing gene annotations. If none, then no annotation is done.
        """
        br_start = self.get_coords('hit')[0]
        br_end = self.get_coords('hit')[1]
        regionStart = targetRegionCoordinates[1] - targetRegionCoordinates[5]
        regionEnd = targetRegionCoordinates[2] + targetRegionCoordinates[5]
        start_in = br_start >= regionStart and br_start <= regionEnd
        end_in = br_end <= regionEnd and br_end >= regionStart
        if targetRegionCoordinates[0] == self.get_seq_name('reference') and (start_in or end_in):
            self.in_target = True
            self.genes = targetRegionCoordinates[3]
        else:
            ann_genes = []
            chrom = self.get_name('hit')
            pos = self.get_coords('hit')
            if chrom.find('chr') == -1:
                chrom = 'chr' + str(chrom)
            for g in annotations.genes:
                gs = annotations.genes[g][1]
                ge = annotations.genes[g][2]
                if chrom == annotations.genes[g][0]:
                    if int(pos[0]) >= gs and int(pos[0]) <= ge:
                        ann_genes.append(g)
                        break
            if len(ann_genes) == 0:
                ann_genes = ['intergenic']
                self.valid = False
            self.genes = ",".join(ann_genes)

    def check_indel(self, nBlatResults):
        indel = False
        utils.log(self.loggingName, 'info', 'Checking if blat result contains an indel variant')
        if self.spans_query() or (nBlatResults == 1 and self.in_target):
            utils.log(self.loggingName, 'info', 'Blat result spans query (%r) or only one blat result (%r) and blat result in target (%r)' % (self.spans_query(), (nBlatResults == 1), self.in_target))
            indel = True
        return indel

        # indel = False
        # indel_size_thresh = int(self.meta_dict['params'].opts['indel_size'])
        # self.logger.info('Checking if blat result contains an indel variant')
        # nhits = 0
        # for i in self.hit_freq:
        #     if i > 0:
        #         nhits += 1
        # if br.spans_query() or (len(self.blat_results) == 1 and br.in_target):
        #     self.logger.info('Blat result spans query (%r) or only one blat result (%r) and blat result in target (%r)' % (br.spans_query(), (len(self.blat_results) == 1), br.in_target))
        #     indel = True
        #     keep_br = br.valid and br.mean_cov < 2 and br.in_target and (br.indel_maxevent_size[0] >= indel_size_thresh) and (not br.rep_man.breakpoint_in_rep[0] and not br.rep_man.breakpoint_in_rep[1])
        #     self.logger.debug('Keep blat result %r' % keep_br)
        #     if keep_br:
        #         brkpt_cov = [self.meta_dict['contig_vals'][1].get_counts(x, x, 'indel') for x in br.query_brkpts]
        #         low_cov = min(brkpt_cov) < self.meta_dict['params'].get_sr_thresh('indel')
        #         flank_match_thresh = True
        #         for fm in br.indel_flank_match:
        #             fm_perc = round((float(fm) / float(br.get_size('query'))) * 100, 2)
        #             if fm_perc < 10.0:
        #                 flank_match_thresh = False
        #             self.logger.info('Indel result has matching flanking sequence of largest indel event of %d (%d of query)' % (fm, fm_perc))
        #         self.logger.info('Indel result has matching flanking sequence of largest indel event (10 perc of query) on both sides (%r)' % flank_match_thresh)
        #         in_ff, span_ff = filter_by_feature(br.get_brkpt_locs(), self.meta_dict['query_region'], self.meta_dict['params'].opts['keep_intron_vars'])
        #         if not in_ff and not low_cov and flank_match_thresh:
        #             self.se = sv_event(br, self.meta_dict['query_region'], self.meta_dict['contig_vals'], self.meta_dict['sbam'])
        #             self.logger.debug('Top hit contains whole query sequence, indel variant')
        #         else:
        #             self.logger.debug('Indel in intron (%r) or low coverage at breakpoints (%r) or minimum segment size < 20 (%r), filtering out.' % (in_ff, low_cov, min(br.query_blocksizes)))
        #     else:
        #         self.logger.debug('Indel failed checking criteria: in annotated gene: %r, mean query coverage < 2: %r, in target: %r, in repeat: %r, indel size < %d: %r' % (br.valid, br.mean_cov, br.in_target, ",".join([str(x) for x in br.rep_man.breakpoint_in_rep]), indel_size_thresh, br.indel_maxevent_size[0] < indel_size_thresh))
        # return indel


# class blat_repeat_manager:
#     def __init__(self):
#         # Booleans for both breakpoints and whether they land in simple repeats
#         self.breakpoint_in_rep = [False, False]
#         self.total_rep_overlap = 0.0
#         self.simple_rep_overlap = 0.0
#         self.other_values = [False, 0.0, [], [False, False]]

#     def setup(self, coords, repeat_locs):
#         self.check_repeat_regions(coords, repeat_locs)

#     def check_repeat_regions(self, coords, repeat_locs):
#         start, end = coords
#         seg_len = float(end - start)
#         in_repeat = False
#         rep_overlap = 0.0
#         simple_overlap = 0.0
#         rep_coords = []
#         filter_reps_edges = [False, False]
#         for rloc in repeat_locs:
#             rchr, rbp1, rbp2, rname = rloc
#             if (rbp1 >= start and rbp1 <= end) or (rbp2 >= start and rbp2 <= end) or (rbp1 <= start and rbp2 >= end):
#                 in_repeat = True
#                 rep_overlap += float(min(rbp2, end) - max(rbp1, start))
#                 rep_coords.append((rbp1, rbp2))
#                 # Simple or low complexity seq repeat for filtering
#                 if rname.find(")n") > -1 or rname.find("_rich") > -1:
#                     simple_overlap += float(min(rbp2, end) - max(rbp1, start))
#                     if (rbp1 <= start and rbp2 >= start):
#                         filter_reps_edges[0] = True
#                     elif (rbp1 <= end and rbp2 >= end):
#                         filter_reps_edges[1] = True
# #        if rep_overlap >= seg_len :
# #          break
#         roverlap = round((float(min(rep_overlap, seg_len)) / float(seg_len)) * 100, 2)
#         self.total_rep_overlap = roverlap
#         self.simple_rep_overlap = round((float(min(simple_overlap, seg_len)) / float(seg_len)) * 100, 2)
#         self.breakpoint_in_rep = filter_reps_edges
#         self.other_values = [in_repeat, roverlap, rep_coords, filter_reps_edges]
