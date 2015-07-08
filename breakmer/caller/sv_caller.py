#! /usr/bin/local/python
# -*- coding: utf-8 -*-

import sys
import os
import math
import pysam
import breakmer.utils as utils

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


class FilterValues:
    """
    """
    def __init__(self):
        self.maxEventSize = None
        self.resultMeanHitFreq = None
        self.brkptCoverages = None
        self.flankMatchPercents = None
        self.minSegmentLen = None
        self.minBrkptKmers = None
        self.seqComplexity = None
        self.startEndMissingQueryCoverage = None
        self.missingQueryCoverage = None
        self.maxSegmentOverlap = None
        self.maxMeanCoverage = None
        self.nReadStrands = None
        self.maxRealignmentGap = None

    def set_indel_values(self, blatResult, brkptCoverages):
        """ """
        self.resultMeanHitFreq = blatResult.meanCov
        self.maxEventSize = blatResult.indel_maxevent_size[0]
        self.brkptCoverages = [min(brkptCoverages), max(brkptCoverages)]
        self.flankMatchPercents = []
        for flankMatch in blatResult.indel_flank_match:
            self.flankMatchPercents.append(round((float(flankMatch) / float(blatResult.get_seq_size('query'))) * 100, 2))

    def set_trl_values(self, svEvent):
        """ """
        blatResult = svEvent.blatResultsSorted[0][0]
        breakpoints = svEvent.brkpts
        self.minSegmentLen = blatResult.get_nmatch_total()
        # Set the min to be the surrounding area of breakpoints, and max to be the direct breakpoints
        self.brkptCoverages = [min(breakpoints.counts['n']), max(breakpoints.counts['d'])]
        self.minBrkptKmers = min(breakpoints.kmers)
        # Sequence complexity of the shortest blat aligned sequence
        self.seqComplexity = svEvent.get_seq_complexity()
        self.startEndMissingQueryCoverage = svEvent.get_startend_missing_query_coverage()
        self.missingQueryCoverage = svEvent.get_missing_query_coverage()
        self.maxSegmentOverlap = max(blatResult.seg_overlap)
        self.nReadStrands = svEvent.check_read_strands()
        self.maxRealignmentGap = max(blatResult.gaps.get_gap_sizes())
        # Use this value to determine the uniqueness of the realignment
        self.maxMeanCoverage = svEvent.get_max_meanCoverage()

    def set_rearr_values(self, svEvent):
        """ """
        breakpoints = svEvent.brkpts
        blatResult = svEvent.blatResultsSorted[0][0]
        self.minBrkptKmers = min(breakpoints.kmers)
        self.minSegmentLen = blatResult.get_nmatch_total()
        self.missingQueryCoverage = svEvent.get_missing_query_coverage()
        self.maxSegmentOverlap = max(blatResult.seg_overlap)
        self.maxMeanCoverage = svEvent.get_max_meanCoverage()

    def get_formatted_output_values(self, svType, svSubtype):
        """ """
        outputValues = {}
        if svType == 'indel':
            outputValues['maxeventSize'] = self.maxEventSize
            outputValues['meanHitFreq'] = self.resultMeanHitFreq
            # Store the minimum value.
            outputValues['breakpointCoverages'] = self.brkptCoverages[0]
            outputValues['minSeqEdgeRealignmentPercent'] = min(self.flankMatchPercents)
        elif svType == 'rearrangement':
            outputValues['minBrkptKmers'] = self.minBrkptKmers
            outputValues['minSegmentLen'] = self.minSegmentLen
            outputValues['missingQueryCoverage'] = self.missingQueryCoverage
            outputValues['maxSegmentOverlap'] = self.maxSegmentOverlap
            outputValues['maxSegmentMeanHitFreq'] = self.maxMeanCoverage
            if svSubtype == 'trl':
                outputValues['breakpointCoverages'] = self.brkptCoverages
                outputValues['sequenceComplexity'] = self.seqComplexity
                outputValues['startEndMissingQueryCoverage'] = self.startEndMissingQueryCoverage
                outputValues['nReadStrands'] = self.nReadStrands
                outputValues['maxRealignmentGapSize'] = self.maxRealignmentGap

        outputList = []
        for key, value in outputValues.items():
            outputList.append(key + '=' + str(value))
        return ';'.join(outputList)


class SVResult:
    """
    """
    def __init__(self):
        self.loggingName = 'breakmer.caller.sv_caller'
        self.fullBreakpointStr = None
        self.targetBreakpointStr = None
        self.alignCigar = None
        self.totalMismatches = None
        self.strands = None
        self.totalMatching = None
        self.svType = ''
        self.svSubtype = None
        self.splitReadCount = None
        self.nKmers = None
        self.discReadCount = None
        self.contigId = None
        self.contigSeq = None
        self.targetName = None
        self.breakpointCoverageDepth = None
        self.description = None
        self.genes = None
        self.repeatOverlapPercent = None
        self.realignmentUniqueness = None
        self.filtered = {'status': False, 'reason': ''}
        self.filterValues = FilterValues()

    def format_indel_values(self, svEvent):
        self.targetName = svEvent.contig.get_target_name()
        self.contigSeq = svEvent.get_contig_seq()
        self.contigId = svEvent.get_contig_id()
        blatResult = svEvent.blatResults[0][1]
        self.genes = blatResult.get_gene_anno()
        self.repeatOverlapPercent = 0.0
        self.totalMatching = blatResult.get_nmatch_total()
        self.realignmentUniqueness = blatResult.meanCov
        self.totalMismatches = blatResult.get_nmatches('mismatch')
        self.strands = blatResult.strand
        self.targetBreakpointStr = svEvent.get_brkpt_str('target')
        self.breakpointCoverageDepth = svEvent.get_brkpt_depths()
        # List of insertion or deletion sizes that coorespond with the breakpoints
        self.description = blatResult.indel_sizes
        self.alignCigar = blatResult.cigar
        self.svType = 'indel'
        contigCountTracker = svEvent.contig.get_contig_count_tracker()
        contigBrkpts = []
        for x in blatResult.breakpts.contigBreakpoints:
            for bp in x:
                contigBrkpts.append(bp)
        self.splitReadCount = [contigCountTracker.get_counts(x, x, 'indel') for x in contigBrkpts]
        self.filterValues.set_indel_values(blatResult, self.splitReadCount)
        print 'Formatting indel values', self.contigSeq
        print 'contig id', self.contigId
        print 'target', self.targetName

    def format_rearrangement_values(self, svEvent):
        """ """
        utils.log(self.loggingName, 'info', 'Resolving SVs call from blat results')
        # Sort the stored blat results by the number of matches to the reference sequence.
        blatResSorted = sorted(svEvent.blatResults, key=lambda x: x[0])
        # brkpts = SVBreakpoints() # {'t':{'in_target':None, 'other':None }, 'formatted':[], 'r':[], 'q': [[0,0],[]], 'chrs':[], 'brkpt_str':[], 'tcoords':[], 'f': []}
        # res_values = {'target_breakpoints':[], 'align_cigar':[], 'sv_type':'', 'strands':[], 'mismatches':[], 'repeat_matching':[], 'anno_genes': [], 'disc_read_count': 0}
        resultValid = {'valid': True, 'repeatValid': True}
        maxRepeat = 0.0

        self.totalMatching = []
        self.repeatOverlapPercent = []
        self.realignmentUniqueness = []
        self.genes = []
        self.alignCigar = []
        self.strands = []
        self.totalMismatches = []

        for i, blatResultTuple in enumerate(blatResSorted):
            blatResult = blatResultTuple[1]
            resultValid['valid'] = resultValid['valid'] and blatResult.valid
            # resultValid['repeatValid'] = resultValid['repeatValid'] and (blatResult.rep_man.simple_rep_overlap > 75.0)
            maxRepeat = max(maxRepeat, blatResult.repeat_overlap)
            # self.repeatMatching.append(":".join([str(blatResult.repeat_overlap), str(blatResult.get_nmatch_total()), str(round(blatResult.mean_cov, 3))]))
            self.repeatOverlapPercent.append(blatResult.repeat_overlap)
            self.realignmentUniqueness.append(blatResult.meanCov)
            self.totalMatching.append(blatResult.get_nmatch_total())
            self.genes.append(blatResult.get_gene_anno())
            self.alignCigar.append(blatResult.cigar)
            self.strands.append(blatResult.strand)
            self.totalMismatches.append(blatResult.get_nmatches('mismatch'))
            svEvent.brkpts.update_brkpt_info(blatResult, i, i == (len(blatResSorted) - 1))

        # Sort the blatResultsSorted list by the lowest matching result to the highest matching result
        svEvent.blatResultsSorted = sorted(svEvent.blatResultsSorted, key=lambda x: x[1])
        if svEvent.brkpts.diff_chr():
            # translocation event
            svEvent.set_brkpt_counts('trl')
            self.discReadCount = svEvent.get_disc_read_count()
            self.svType = 'rearrangement'
            self.svSubtype = 'trl'
            self.filterValues.set_trl_values(svEvent)
        else:
            svEvent.set_brkpt_counts('rearr')
            rearrType, discReadCount = svEvent.define_rearr()
            self.svType = 'rearrangement'
            if rearrType != 'rearrangement':
                self.svSubtype = rearrType
            self.discReadCount = discReadCount
            self.genes = list(set(self.genes))
            self.filterValues.set_rearr_values(svEvent)
        self.targetName = svEvent.contig.get_target_name()
        self.fullBreakpointStr = svEvent.get_brkpt_str()
        self.targetBreakpointStr = svEvent.get_brkpt_str('target')
        self.breakpointCoverageDepth = svEvent.get_brkpt_depths()
        self.splitReadCount = svEvent.get_splitread_count()
        self.contigSeq = svEvent.get_contig_seq()
        self.contigId = svEvent.get_contig_id()

    # def get_output_type(self):
    #     """ """
    #     outputType = self.svType
    #     if self.svSubtype != '':
    #         outputType += '_' + self.svSubtype
    #     return outputType

    def set_filtered(self, filterReason):
        """ """
        self.filtered['status'] = True
        self.filtered['reason'] = filterReason

    def get_old_formatted_output_values(self):
        """ """
        headerStr = ['genes',
                     'target_breakpoints',
                     'align_cigar',
                     'mismatches',
                     'strands',
                     'rep_overlap_segment_len',
                     'sv_type',
                     'split_read_count',
                     'nkmers',
                     'disc_read_count',
                     'breakpoint_coverages',
                     'contig_id',
                     'contig_seq'
                     ]

        brkptStr = ','.join([str(x) for x in item])
        if self.svType == 'indel':
            brkptStr += ' (' + ','.join([str(x) for x in self.descript]) + ')'

        repOverlap_segLen_hitFreq = []
        for i in self.totalMatching:
            repOverlap_segLen_hitFreq.append('0.0:' + str(matchLen) + ':0.0')

        nkmers = '0'

        outList = [self.targetName,
                   self.brkptStr,
                   self.alignCigar,
                   self.totalMismatches,
                   self.strands,
                   repOverlap_segLen_hitFreq,
                   self.svType,
                   self.splitReadCount,
                   nkmers,
                   self.discReadCount,
                   self.breakpointCoverageDepth,
                   self.contigId,
                   self.contigSeq,
                   ]

        outListStr = []
        for item in outList:
            if not isinstance(item, list):
                outListStr.append(str(item))
            else:
                outListStr.append(','.join([str(x) for x in item]))

        formattedFilterValsStr = self.filterValues.get_formatted_output_values(self.svType, self.svSubtype)
        outListStr.append(formattedFilterValsStr)
        return ('\t'.join(headerStr), '\t'.join(outListStr))

    def get_formatted_output_values(self):
        """ """
        headerStr = ['Target_Name',
                     'SV_type',
                     'SV_subtype',
                     'Description',
                     'All_genomic_breakpoints',
                     'Target_genomic_breakpoints',
                     'Split_read_counts',
                     'Discordant_read_counts',
                     'Read_depth_at_genomic_breakpoints',
                     'Align_cigar',
                     'Strands',
                     'Total_mismatches',
                     'Total_matching',
                     'Contig_ID',
                     'Contig_length',
                     'Contig_sequence',
                     'Filtered',
                     'Filtered_reason',
                     'Filter_values'
                     ]

        outList = [self.targetName,
                   self.svType,
                   self.svSubtype,
                   self.description,
                   self.fullBreakpointStr,
                   self.targetBreakpointStr,
                   self.splitReadCount,
                   self.discReadCount,
                   self.breakpointCoverageDepth,
                   self.alignCigar,
                   self.strands,
                   self.totalMismatches,
                   self.totalMatching,
                   self.contigId,
                   len(self.contigSeq),
                   self.contigSeq,
                   self.filtered['status'],
                   self.filtered['reason']
                   ]

        outListStr = []
        for item in outList:
            if not isinstance(item, list):
                outListStr.append(str(item))
            else:
                outListStr.append(','.join([str(x) for x in item]))

        formattedFilterValsStr = self.filterValues.get_formatted_output_values(self.svType, self.svSubtype)
        outListStr.append(formattedFilterValsStr)
        return ('\t'.join(headerStr), '\t'.join(outListStr))

    # def format_result(self, values):
    #     res_lst = []
    #     if values:
    #         for v in values:
    #             if not isinstance(values[v], list):
    #                 values[v] = [values[v]]
    #             self.resultValues[v] = ",".join([str(x) for x in values[v]])
    #     if self.resultValues['sv_subtype']:
    #         self.resultValues['sv_type'] += '_' + self.resultValues['sv_subtype']
    #     return self.get_values()

    # def get_values(self):
    #     lst = ['anno_genes', 'target_breakpoints', 'align_cigar', 'mismatches', 'strands', 'repeat_matching', 'sv_type', 'split_read_count', 'nkmers', 'disc_read_count', 'breakpoint_coverages', 'contig_id', 'contig_seq']
    #     out_lst = []
    #     for l in lst:
    #         if self.resultValues[l] == 'trl':
    #             self.resultValues[l] = 'rearrangement'
    #         out_lst.append(str(self.resultValues[l]))
    #         if l == 'target_breakpoints':
    #             self.resultValues['breakpoint_coverages'] = self.get_brkpt_coverages()
    #     return out_lst


class SVBreakpoints:
    def __init__(self):
        self.loggingName = 'breakmer.caller.sv_caller'
        self.t = {'target': None, 'other': None}
        self.formatted = []
        self.r = []
        self.q = [[0, 0], []]
        self.chrs = []
        self.brkptStr = []
        self.tcoords = []
        self.f = []
        self.counts = {'n': [], 'd': [], 'b': []}
        self.kmers = []
        # Standard format for storing genomic breakpoints for outputtting rsults
        # List of tuples containing ('chr#', bp1, bp2), there will be multiple bp for deletions and
        # only one bp for insertions or rearrangment breakpoints.
        self.genomicBrkpts = {'target': [], 'other': []}

    def update_brkpt_info(self, br, i, last_iter):
        """Infer the breakpoint information from the blat result for rearrangments.
        """
        chrom = 'chr' + br.get_seq_name('ref')
        ts, te = br.get_coords('ref')
        qs, qe = br.get_coords('query')
        targetKey = 'target' if br.in_target else 'other'
        self.chrs.append(br.get_seq_name('ref'))
        self.tcoords.append((ts, te))
        tbrkpt = []
        filt_rep_start = None
        if i == 0:
            self.q[0] = [max(0, qs - 1), qe]
            self.q[1].append([qe, qe - self.q[0][0], None])
            tbrkpt = [te]
            filt_rep_start = br.filter_reps_edges[0]
            if br.strand == '-':
                tbrkpt = [ts]
                filt_rep_start = br.filter_reps_edges[0]
            self.genomicBrkpts[targetKey].append((chrom, tbrkpt[0]))
            br.set_sv_brkpt((chrom, tbrkpt[0]), 'rearrangement', targetKey)
        elif last_iter:
            self.q[1][-1][2] = qe - self.q[1][-1][0]
            self.q[1].append([qs, qs - self.q[0][0], qe - qs])
            tbrkpt = [ts]
            filt_rep_start = br.filter_reps_edges[0]
            if br.strand == '-':
                tbrkpt = [te]
                filt_rep_start = br.filter_reps_edges[1]
            self.genomicBrkpts[targetKey].append((chrom, tbrkpt[0]))
            br.set_sv_brkpt((chrom, tbrkpt[0]), 'rearrangement', targetKey)
        else:
            self.q[1][-1][2] = qe - self.q[1][-1][1]
            self.q[1].append([qs, qs - self.q[0][0], qe - qs])
            self.q[1].append([qe, qe - qs, None])
            self.q[0] = [qs, qe]
            tbrkpt = [ts, te]
            self.genomicBrkpts[targetKey].append((chrom, ts, te))
            if br.strand == '+':
                br.set_sv_brkpt((chrom, ts, te), 'rearrangement', targetKey)
            if br.strand == '-':
                filt_rep_start = br.filter_reps_edges[1]
                tbrkpt = [te, ts]
                self.genomicBrkpts[targetKey].append((chrom, te, ts))
                br.set_sv_brkpt((chrom, te, ts), 'rearrangement', targetKey)

        self.brkptStr.append('chr' + str(br.get_seq_name('ref')) + ":" + "-".join([str(x) for x in tbrkpt]))
        self.r.extend(tbrkpt)
        self.f.append(filt_rep_start)
        self.t[targetKey] = (br.get_seq_name('ref'), tbrkpt[0])
        self.formatted.append('chr' + str(br.get_seq_name('ref')) + ":" + "-".join([str(x) for x in tbrkpt]))

    def set_indel_brkpts(self, blatResult):
        """ """
        # List of tuples for indel breakpoints parsed from the blat result ('chr#', bp1, bp2)
        self.genomicBrkpts['target'] = blatResult.get_genomic_brkpts()
        for brkpt in self.genomicBrkpts['target']:
            blatResult.set_sv_brkpt(brkpt, 'indel', 'target')
        # brkpt_out = []
        # bp_str = []
        # chrm = 'chr' + str(blatResult.get_seq_name('ref'))
        # if len(self.breakpts) > 0:
        #     for b, s in zip(self.breakpts, self.indel_sizes):
        #         if len(b) > 1:
        #             bb = "-".join([str(x) for x in b])
        #         else:
        #             bb = str(b[0])
        #         bstr = chrm + ":" + bb
        #         if with_sizes:
        #             bstr += " " + "(" + s + ")"
        #         bp_str.append(bstr)
        #     brkpt_out.append(",".join(bp_str))
        # return ",".join(brkpt_out)

    def diff_chr(self):
        """Determine if the stored realignment results are on multiple chromosomes - indicating a
        translocation event.
        """
        if len(set(self.chrs)) == 1:
            return False
        else:
            return True

    def get_target_brkpt(self, key):
        """ """
        return self.target[key]

    def get_brkpt_str(self, targetKey):
        """ """
        if targetKey is None:
            brkptStr = ''
            for key in self.genomicBrkpts:
                outStr = self.get_brkpt_str(key)
                if brkptStr == '':
                    brkptStr = outStr
                elif outStr != '':
                    brkptStr += ',' + outStr
            return brkptStr
        else:
            brkptStr = []
            for genomicBrkpts in self.genomicBrkpts[targetKey]:
                chrom = genomicBrkpts[0]
                bps = genomicBrkpts[1:]
                brkptStr.append(chrom + ':' + '-'.join([str(x) for x in bps]))
            return ','.join(brkptStr)

    def get_brkpt_depths(self, sampleBamFn):
        """ """
        depths = []
        bamfile = pysam.Samfile(sampleBamFn, 'rb')
        for key in self.genomicBrkpts:
            for genomicBrkpt in self.genomicBrkpts[key]:
                chrom = genomicBrkpt[0].strip('chr')
                bps = genomicBrkpt[1:]
                for bp in bps:
                    alignedDepth = 0
                    alignedReads = bamfile.fetch(str(chrom), int(bp), int(bp) + 1)
                    for alignedRead in alignedReads:
                        if alignedRead.is_duplicate or alignedRead.is_qcfail or alignedRead.is_unmapped or alignedRead.mapq < 10:
                            continue
                        alignedDepth += 1
                    depths.append(alignedDepth)
        return depths

    def get_splitread_count(self):
        """ """
        return self.counts['b']

    def set_counts(self, svType, contig):
        """ """
        # avg_comp, comp_vec = calc_contig_complexity(self.contig_seq)
        # brkpt_rep_filt = False
        # brkpt_counts = {'n': [], 'd': [], 'b': []}
        # brkpt_kmers = []
        contigCountTracker = contig.get_contig_count_tracker()
        for qb in self.q[1]:
            left_idx = qb[0] - min(qb[1], 5)
            right_idx = qb[0] + min(qb[2], 5)
            bc = contigCountTracker.get_counts(left_idx, right_idx, svType)
            self.counts['n'].append(min(bc))
            self.counts['d'].append(min(contigCountTracker.get_counts((qb[0] - 1), (qb[0] + 1), svType)))
            self.counts['b'].append(contigCountTracker.get_counts(qb[0], qb[0], svType))
            self.kmers.append(contig.get_kmer_locs()[qb[0]])
            # brkpt_rep_filt = brkpt_rep_filt or (comp_vec[qb[0]] < (avg_comp / 2))
            utils.log(self.loggingName, 'debug', 'Read count around breakpoint %d : %s' % (qb[0], ",".join([str(x) for x in bc])))
        utils.log(self.loggingName, 'debug', 'Kmer count around breakpoints %s' % (",".join([str(x) for x in self.kmers])))
        # brkpt_rep_filt = brkpt_rep_filt or (len(filter(lambda x: x, brkpts['f'])) > 0)
        # return brkpt_counts, brkpt_kmers, brkpt_rep_filt


class SVEvent:
    def __init__(self, blatResult, contig, svType):
        self.loggingName = 'breakmer.caller.sv_caller'
        self.svType = svType
        self.svSubtype = ''
        self.events = []
        self.blatResults = []
        self.blatResultsSorted = []
        self.annotated = False
        # self.sample_bam = sample_bam
        self.qlen = 0
        self.nmatch = 0
        self.in_target = False
        self.contig = contig
        # self.query_region = query_region
        # self.contig_seq = contig_vals[0]
        # self.contig_rcounts = contig_vals[1]
        # self.contig_id = contig_vals[2]
        # self.contig_reads = contig_vals[3]
        # self.nkmers = contig_vals[4]
        # self.contig_kmer_locs = contig_vals[5]
        self.valid = True
        self.in_rep = True
        self.querySize = None
        self.queryCoverage = [0] * len(contig.seq)
        self.homology = {'non': [False, None], 'in': [False, None]}
        self.brkpts = SVBreakpoints()
        self.resultValues = SVResult()
        self.add(blatResult)

    def add(self, blatResult):
        queryStartCoord = blatResult.alignVals.get_coords('query', 0)
        queryEndCoord = blatResult.alignVals.get_coords('query', 1)
        self.blatResults.append((queryStartCoord, blatResult))

        # Add the number of hits to the query region
        for i in range(queryStartCoord, queryEndCoord):
            # print 'incrementing query cov', i, self.queryCoverage, len(self.queryCoverage), blatResult.get_coords('query')
            self.queryCoverage[i] += 1
        if not self.querySize:
            self.querySize = blatResult.get_seq_size('query')
        self.qlen += blatResult.get_query_span()
        self.nmatch += blatResult.get_nmatch_total()
        self.in_target = self.in_target or blatResult.in_target
        self.in_rep = self.in_rep and (blatResult.repeat_overlap > 75.0)
        self.valid = self.valid and blatResult.valid
        self.blatResultsSorted.append((blatResult, blatResult.get_nmatch_total()))

    def result_valid(self):
        valid = False
        if (len(self.blatResults) > 1) and self.in_target:
            valid = True
            # nMissingQueryCoverage = len(filter(lambda y: y, map(lambda x: x == 0, self.queryCoverage)))
            # if nMissingQueryCoverage < self.meta_dict['params'].get_min_segment_length('trl'):
            #     valid = True
        return valid

    def has_annotations(self):
        """ """
        return True

    def get_genomic_brkpts(self):
        """ """
        return self.brkpts.genomicBrkpts

    def check_previous_add(self, br):
        ncoords = br.get_coords('query')
        prev_br, prev_nmatch = self.blatResultsSorted[-1]
        prev_coords = prev_br.get_coords('query')
        if ncoords[0] == prev_coords[0] and ncoords[1] == prev_coords[1]:
            n_nmatch = br.get_nmatch_total()
            if abs(prev_nmatch - n_nmatch) < 10:
                if not prev_br.in_target and br.in_target:
                    self.blatResultsSorted[-1] = (br, n_nmatch)
                    self.blatResults[-1] = (ncoords[0], br)
                    self.in_target = True

    def format_indel_values(self):
        """
        """
        self.brkpts.set_indel_brkpts(self.blatResults[0][1])
        self.resultValues.format_indel_values(self)

    def format_rearr_values(self):
        """
        """
        self.resultValues.format_rearrangement_values(self)

    def get_disc_read_count(self):
        """
        """
        discReads = self.contig.get_disc_reads()
        discReadCount = 0
        nonTargetBrkptChr, nonTargetBrkptBp = self.brkpts.get_target_brkpt('other')
        targetBrkptChr, targetBrkptBp = self.brkpts.get_target_brkpt('in_target')
        if nonTargetBrkptChr in discReads:
            for p1, p2 in discReads[nonTargetBrkptChr]:
                d1 = abs(p1 - targetBrkptBp)
                d2 = abs(p2 - nonTargetBrkptBp)
                if d1 <= 1000 and d2 <= 1000:
                    discReadCount += 1
        return discReadCount

    def get_brkpt_str(self, targetKey=None):
        """ """
        return self.brkpts.get_brkpt_str(targetKey)

    def get_brkpt_depths(self):
        """
        """
        return self.brkpts.get_brkpt_depths(self.contig.get_sample_bam_fn())

    def get_splitread_count(self):
        """ """
        return self.brkpts.get_splitread_count()

    def set_filtered(self, filterReason):
        """ """
        self.resultValues.set_filtered(filterReason)

    def get_missing_query_coverage(self):
        """ """
        return len(filter(lambda y: y, map(lambda x: x == 0, self.queryCoverage)))

    def get_formatted_output_values(self):
        """ """
        return self.resultValues.get_formatted_output_values()

    def get_contig_seq(self):
        """ """
        return self.contig.seq

    def get_contig_id(self):
        """ """
        return self.contig.get_id()

    def set_brkpt_counts(self, svType):
        """ """
        self.brkpts.set_counts(svType, self.contig)

    def check_overlap(self, coord1, coord2):
        contained = False
        if coord1[0] >= coord2[0] and coord1[1] <= coord2[1]:
            contained = True
        elif coord2[0] >= coord1[0] and coord2[1] <= coord1[1]:
            contained = True
        return contained

    def define_rearr(self):
        # vrt = self.contig.get_read_variation()
        varReads = self.contig.get_var_reads('sv')
        strands = self.resultValues.strands
        brkpts = self.brkpts.r
        tcoords = self.brkpts.tcoords
        svType = 'rearrangement'
        rs = 0
        hit = False
        if len(strands) < 3:
            if not self.check_overlap(tcoords[0], tcoords[1]):
                utils.log(self.loggingName, 'debug', 'Checking rearrangement svType, strand1 %s, strand2 %s, breakpt1 %d, breakpt %d' % (strands[0], strands[1], brkpts[0], brkpts[1]))
                if (strands[0] != strands[1]) and (brkpts[0] < brkpts[1]):
                    # Inversion
                    # Get discordantly mapped read-pairs
                    utils.log(self.loggingName, 'debug', 'Inversion event identified.')
                    hit = True
                    svType = 'inversion'
                    rs = varReads.check_inv_readcounts(brkpts)
                    # for readPair in varReads.inv:
                    #     r1p, r2p, r1s, r2s, qname = readPair
                    #     if r1s == 1 and r2s == 1:
                    #         if (r1p <= brkpts[0]) and (r2p <= brkpts[1] and r2p >= brkpts[0]):
                    #             rs += 1
                    #     else:
                    #         if (r1p <= brkpts[1] and r1p >= brkpts[0]) and r2p >= brkpts[1]:
                    #             rs += 1
                elif (strands[0] == '+' and strands[1] == '+') and (brkpts[0] > brkpts[1]):
                    utils.log(self.loggingName, 'debug', 'Tandem duplication event identified.')
                    hit = True
                    svType = 'tandem_dup'
                    # Tandem dup
                    rs = varReads.check_td_readcounts(brkpts)
                    # for readPair in varReads.td:
                    #     r1p, r2p, r1s, r2s, qname = readPair
                    #     if (r1p <= brkpts[0] and r1p >= brkpts[1]) and (r2p <= brkpts[1] and r2p >= brkpts[0]):
                    #         rs += 1
        if not hit:
            utils.log(self.loggingName, 'debug', 'Not inversion or tandem dup, checking for odd read pairs around breakpoints')
            rs = varReads.check_other_readcounts(brkpts)
            # rs = [0] * len(brkpts)
            # for i in range(len(brkpts)):
            #     b = brkpts[i]
                # rs = varReads.check_other_readcounts(b)
                # for readPair in varReads.reads.other:
                #     r1p, r2p, r1s, r2s, qname = readPair
                #     if abs(r1p - b) <= 300 or abs(r2p - b) <= 300:
                #         utils.log(self.loggingName, 'debug', 'Adding read support from read %s, with strands %s, %s and positions %d, %d for breakpoint at %d' % (qname, r1s, r2s, r1p, r2p, b))
                #         rs[i] += 1
            # rs = max(rs)
        return svType, rs

    def get_max_meanCoverage(self):
        """Return the highest mean hit frequency among all blat results stored.
        """
        maxMeanCov = 0
        for blatResult, nBasesAligned in self.blatResultsSorted:
            if int(blatResult.meanCov) > int(maxMeanCov):
                maxMeanCov = int(blatResult.meanCov)

    def check_read_strands(self):
        """
        """
        same_strand = False
        strands = []
        for read in self.contig.reads:
            strand = read.id.split("/")[1]
            strands.append(strand)
        if len(set(strands)) == 1:
            same_strand = True
        utils.log(self.loggingName, 'debug', 'Checking read strands for contig reads %s' % (",".join([read.id for read in self.contig_reads])))
        utils.log(self.loggingName, 'debug', 'Reads are on same strand: %r' % same_strand)
        return len(set(strands))

    def get_seq_complexity(self):
        """Get the 3-mer complexity of the shortest aligned blat sequence.
        """
        blatResult, nBasesAligned = self.blatResultsSorted[0]
        alignedSeq = self.contig.seq[blatResult.qstart():blatResult.qend()]
        merSize = 3
        utils.log(self.loggingName, 'debug', 'Checking sequence complexity of blat result segment %s using %d-mers' % (alignedSeq, merSize))
        nmers = {}
        totalMersPossible = len(alignedSeq) - 2
        for i in range(len(alignedSeq) - (merSize - 1)):
            nmers[str(alignedSeq[i:i + merSize]).upper()] = True
        complexity = round((float(len(nmers)) / float(totalMersPossible)) * 100, 4)
        utils.log(self.loggingName, 'debug', 'Complexity measure %f, based on %d unique %d-mers observed out of a total of %d %d-mers possible' % (complexity, len(nmers), merSize, totalMersPossible, merSize))
        return complexity

    def get_startend_missing_query_coverage(self):
        """Calculate the percentage of the contig sequence that is not realigned to the reference, only examining the
        beginning and end of the contig sequence.
        """
        missingCov = 0
        for i in self.queryCoverage:
            if i == 0:
                missingCov += 1
            else:
                break
        for i in reversed(self.queryCoverage):
            if i == 0:
                missingCov += 1
            else:
                break
        percentMissing = round((float(missingCov) / float(len(self.contig.seq))) * 100, 4)
        utils.log(self.loggingName, 'debug', 'Calculated %f missing coverage of blat query sequence at beginning and end' % percentMissing)
        return percentMissing


class ContigCaller:
    """
    """
    def __init__(self, realignment, contig, params):
        self.realignment = realignment
        self.contig = contig
        self.params = params
        self.clippedQs = []
        self.svEvent = None
        self.loggingName = 'breakmer.caller.sv_caller'

    def call_svs(self):
        """ """
        print 'realigment results', self.realignment.has_results()
        if not self.realignment.has_results():
            utils.log(self.loggingName, 'info', 'No blat results file exists, no calls for %s.' % self.contig.get_id())
        else:
            print 'Results exist'
            utils.log(self.loggingName, 'info', 'Making variant calls from blat results %s' % self.realignment.get_result_fn())
            if self.check_indels():
                print 'Check indels passed'
                # self.result = self.realignment.get_indel_result()
                self.svEvent.format_indel_values()
            elif self.check_svs():
                print 'Check svs passed'
                # self.result = self.realignment.get_svs_result()
                self.svEvent.format_rearr_values()
        # Format the result into a
        print 'result', self.svEvent
        return self.svEvent

    def check_indels(self):
        """ """
        hasIndel = False
        blatResults = self.realignment.get_blat_results()
        for i, blatResult in enumerate(blatResults):
            if i == 0 and blatResult.check_indel(len(blatResults)):
                hasIndel = True
                utils.log(self.loggingName, 'info', 'Contig has indel, returning %r' % hasIndel)
                self.svEvent = SVEvent(blatResult, self.contig, 'indel')
                return hasIndel
            else:
                utils.log(self.loggingName, 'debug', 'Storing clipped blat result start %d, end %d' % (blatResult.qstart(), blatResult.qend()))
                self.clippedQs.append((blatResult.qstart(), blatResult.qend(), blatResult, i))
                # self.realignment.store_clipped_queryseq((blatResult.qstart(), blatResult.qend(), blatResult, i))
        utils.log(self.loggingName, 'info', 'Contig does not have indel, return %r' % hasIndel)
        return hasIndel

    def check_svs(self):
        """ """
        utils.log(self.loggingName, 'info', 'Checking for SVs')
        gaps = [(0, self.realignment.get_qsize())]
        if len(self.clippedQs) > 1:
            utils.log(self.loggingName, 'debug', 'Iterating through %d clipped blat results.' % len(self.clippedQs))
            mergedClip = [0, None]
            for i, clippedQs in enumerate(self.clippedQs):
                qs, qe, blatResult, idx = clippedQs
                utils.log(self.loggingName, 'debug', 'Blat result with start %d, end %d, chrom %s' % (qs, qe, blatResult.get_seq_name('ref')))
                gaps = self.iter_gaps(gaps, self.clippedQs[i], i)
                if self.svEvent.qlen > mergedClip[0]:
                    mergedClip = [self.svEvent.qlen, self.svEvent]
            self.svEvent = mergedClip[1]
        else:
            utils.log(self.loggingName, 'info', 'There are no more than 1 clipped blat results, not continuing with SVs calling.')
        if self.svEvent and self.svEvent.result_valid():
            return True
        else:
            print 'svs result invalid'
            self.svEvent = None
            return False

    def iter_gaps(self, gaps, clippedQuerySeqVals, iterIdx):
        """ """
        new_gaps = []
        qs, qe, blatResult, idx = clippedQuerySeqVals
        hit = False
        for gap in gaps:
            gs, ge = gap
            utils.log(self.loggingName, 'debug', 'Gap coords %d, %d' % (gs, ge))
            startWithinGap = (qs >= gs and qs <= ge)
            endWithinGap = (qe <= ge and qe >= gs)
            gapEdgeDistStart = (qs <= gs) and ((gs - qs) < 15)
            gapEdgeDistEnd = (qe >= ge) and ((qe - ge) < 15)
            if startWithinGap or endWithinGap or (gapEdgeDistStart and (endWithinGap or gapEdgeDistEnd)) or (gapEdgeDistEnd and (startWithinGap or gapEdgeDistStart)):
                ngap = []
                if qs > gs:
                    if (qs - 1 - gs) > 10:
                        ngap.append((gs, qs - 1))
                if qe < ge:
                    if (ge - qe + 1) > 10:
                        ngap.append((qe + 1, ge))
                if iterIdx == 0:
                    utils.log(self.loggingName, 'debug', 'Creating SV event from blat result with start %d, end %d' % (qs, qe))
                    self.svEvent = SVEvent(blatResult, self.contig, 'rearrangement')
                    new_gaps.extend(ngap)
                    hit = True
                elif self.check_add_br(qs, qe, gs, ge, blatResult):
                    utils.log(self.loggingName, 'debug', 'Adding blat result to event')
                    new_gaps.extend(ngap)
                    self.svEvent.add(blatResult)
                    hit = True
                else:
                    new_gaps.append(gap)
            else:
                new_gaps.append(gap)
            utils.log(self.loggingName, 'debug', 'New gap coords %s' % (",".join([str(x) for x in new_gaps])))
        if not hit:
            self.svEvent.check_previous_add(blatResult)
        return new_gaps

    def check_add_br(self, qs, qe, gs, ge, blatResult):
        """ """
        utils.log(self.loggingName, 'info', 'Checking to add blat result with start %d, end %d' % (qs, qe))
        add = False
        # Calc % of segment overlaps with gap
        over_perc = round((float(min(qe, ge) - max(qs, gs)) / float(qe - qs)) * 100)
        # Check overlap with other aligned segments
        ov_right = 0
        if qe > ge:
            ov_right = abs(qe - ge)
        ov_left = 0
        if qs < gs:
            ov_left = abs(qs - gs)
        blatResult.set_segment_overlap(ov_left, ov_right)
        max_seg_overlap = max(ov_right, ov_left)
        utils.log(self.loggingName, 'debug', 'Blat query segment overlaps gap by %f' % over_perc)
        utils.log(self.loggingName, 'debug', 'Max segment overlap %f' % max_seg_overlap)
        utils.log(self.loggingName, 'debug', 'Event in target %r and blat result in target %r' % (self.svEvent.in_target, blatResult.in_target))
        if over_perc >= 50 and (max_seg_overlap < 15 or (blatResult.in_target and self.svEvent.in_target)):
            add = True
        utils.log(self.loggingName, 'debug', 'Add blat result to SV event %r' % add)
        return add
