#! /usr/bin/local/python
# -*- coding: utf-8 -*-

import sys
import os
import math
import breakmer.utils as utils

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


class FilterValues:
    # Filter indel
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
    # Filter TRL 
    #         filter = br_valid[1] or (max(brkpt_counts['d']) < params.get_sr_thresh('trl'))
#         self.logger.debug('Check translocation filter')
#         self.logger.debug('All blat result segments are within annotated or pre-specified regions %r' % br_valid[0])
#         self.logger.debug('All blat result segments are within simple repeat regions that cover > 75.0 percent of the segment %r'%br_valid[1])
#         self.logger.debug('The maximum read count support around breakpoints %d meets split read threshold %d'%(max(brkpt_counts['d']),params.get_sr_thresh('trl')))
#         self.logger.debug('The minimum number of kmers at breakpoints %d' % min(brkpt_kmers))
#         self.logger.debug('The maximum repeat overlap by a blat result: %f' % max_repeat)
#         if not filter:
#             self.logger.debug('Filter %r, checking discordant read counts %d' % (filter, disc_read_count))
#             if disc_read_count < 2:
# #        print 'Filter due to repeat', rep_filt
#                 if (self.blatResultsSorted[0][1] < params.get_min_segment_length('trl')) or (min(brkpt_counts['n']) < params.get_sr_thresh('trl')) or (min(brkpt_kmers) == 0) or rep_filt:
#                     self.logger.debug('Shortest segment is < %d bp with %d discordant reads. Filtering.'%(params.get_min_segment_length('trl'), disc_read_count))
#                     self.logger.debug('The minimum read count support for breakpoints %d meets split read threshold %d'%(min(brkpt_counts['n']),params.get_sr_thresh('trl')))
#                     self.logger.debug('The minimum number of kmers at breakpoints %d'%min(brkpt_kmers))
#                     filter = True
#                 elif disc_read_count == 0:
#                     # Check a number of metrics for shortest blat segment
#                     br_qs = self.blatResultsSorted[0][0].qstart()
#                     br_qe = self.blatResultsSorted[0][0].qend()
#                     low_complexity = self.minseq_complexity(self.contig_seq[br_qs:br_qe],3) < 25.0 # Complexity of blat segment
#                     missing_qcov = self.missing_query_coverage() > 5.0
#                     short = self.blatResultsSorted[0][1] <= round(float(len(self.contig_seq))/float(4.0))
#                     self.logger.debug('Checking length of shortest sequence, considered too short %r, %d, %f'%(short, self.blatResultsSorted[0][1], round(float(len(self.contig_seq))/float(4.0))) )
#                     overlap = max(self.blatResultsSorted[0][0].seg_overlap) > 5
#                     gaps_exist = max(self.blatResultsSorted[0][0].gaps['query'][0], self.blatResultsSorted[0][0].gaps['hit'][0]) > 0
#                     low_uniqueness = self.check_uniqueness()
#                     intergenic_regions = 'intergenic' in anno_genes
#                     read_strand_bias = self.check_read_strands()
#                     check_values = [low_complexity, missing_qcov, short, overlap, gaps_exist, low_uniqueness, read_strand_bias, intergenic_regions]
#                     self.logger.debug('Discordant read count of 0 checks %s'%(",".join([str(x) for x in check_values])))
#                     num_checks = 0
#                     for check in check_values:
#                         if check:
#                             num_checks += 1
#                     if num_checks > 1:
#                         self.logger.info('Two or more filter checks, setting filtering to true for contig')
#                         filter = True
    def __init__(self):
        self.maxEventSize = None
        self.resultMeanHitFreq = None
        self.brkptCoverages = None
        self.flankMatchPercents = None
        self.minSegmentLen = None
        self.minBrkptKmers = None
        self.seqComplexity = None
        self.missingQueryCoverage = None
        self.maxSegmentOverlap = None
        self.maxMeanCoverage = None
        self.nReadStrands = None
        self.realignmentGaps = None

    def set_indel_values(self, blatResult, brkptCoverages):
        self.resultMeanHitFreq = blatResult.meanCov
        self.maxEventSize = blatResult.indel_maxevent_size[0]
        self.brkptCoverages = [min(brkptCoverages), max(brkptCoverages)]
        self.flankMatchPercents = []
        for flankMatch in blatResult.indel_flank_match:
            self.flankMatchPercents.append(round((float(flankMatch) / float(blatResult.get_seq_size('query'))) * 100, 2))

    def set_trl_values(self, blatResult, breakpoints):
        self.minSegmentLen = blatResult.get_nmatch_total()
        # Set the min to be the surrounding area of breakpoints, and max to be the direct breakpoints
        self.brkptCoverages = [min(brkptCoverages.counts['n']), max(brkptCoverages.counts['d'])]
        self.minBrkptKmers = min(brkptCoverages.kmers)
        nmers = {}
        total_possible = len(seq) - 2
        for i in range(len(seq) - (N - 1)):
            nmers[str(seq[i:i+N]).upper()] = True
        self.seqComplexity = round((float(len(nmers))/float(total_possible))*100, 4)
        self.missingQueryCoverage = 

    # def check_uniqueness(self):
    #     low_unique = False
    #     for br_vals in self.blatResultsSorted :
    #         if not br_vals[0].in_target :
    #             if br_vals[0].mean_cov > 4 : low_unique = True
    #         else :
    #             if br_vals[0].mean_cov > 10 : low_unique = True
    #     return low_unique

    # def check_read_strands(self):
    #     same_strand = False
    #     strands = []
    #     for read in self.contig_reads :
    #         strand = read.id.split("/")[1] 
    #         strands.append(strand)
    #     if len(set(strands)) == 1 :
    #         same_strand = True
    #     self.logger.debug('Checking read strands for contig reads %s'%(",".join([read.id for read in self.contig_reads])))
    #     self.logger.debug('Reads are on same strand: %r'%same_strand)
    #     return same_strand

    # def minseq_complexity(self, seq, N) :
    #     self.logger.debug('Checking sequence complexity of blat result segment %s using %d-mers'%(seq,N))
    #     nmers = {}
    #     total_possible = len(seq) - 2
    #     for i in range(len(seq) - (N - 1)) :
    #         nmers[str(seq[i:i+N]).upper()] = True
    #     complexity = round((float(len(nmers))/float(total_possible))*100,4)
    #     self.logger.debug('Complexity measure %f, based on %d unique %d-mers observed out of a total of %d %d-mers possible'%(complexity,len(nmers),N,total_possible,N))
    #     return complexity

    # def missing_query_coverage(self) :
    #     missing_cov = 0
    #     for i in self.queryCoverage :
    #         if i == 0 :
    #             missing_cov += 1
    #         else :
    #             break

    #     for i in reversed(self.queryCoverage) :
    #         if i == 0 :
    #             missing_cov += 1
    #         else : 
    #             break

    #     perc_missing = round((float(missing_cov)/float(len(self.contig_seq)))*100, 4)
    #     self.logger.debug('Calculated %f missing coverage of blat query sequence at beginning and end'%perc_missing)
    #     return perc_missing

class SVResult:
    def __init__(self):
        self.loggingName = 'breakmer.caller.sv_caller'
        self.genes = None
        self.targetBreakpoints = None
        self.alignCigar = None
        self.mismatches = None
        self.strands = None
        self.repeatMatching = None
        self.svType = {'type': None, 'subtype': None}
        self.splitReadCount = None
        self.nKmers = None
        self.discReadCount = None
        self.contigId = None
        self.breakpointCoverages = None
        self.filterValues = FilterValues()

    def format_event_values(self, svEvent):
        """ """
        if svEvent.svType == 'indel':
            self.format_indel_values(svEvent)
        else:
            self.format_rearrangement_values(svEvent)

    def format_indel_values(self, svEvent):
        blatResult = svEvent.blatResults[0][1]
        self.genes = blatResult.get_gene_anno()
        self.repeatMatching = '0.0:' + str(blatResult.get_nmatch_total())
        self.mismatches = blatResult.get_nmatches('mismatch')
        self.strands = blatResult.strand
        self.targetBreakpoints = blatResult.get_breakpoint_str(True)
        self.alignCigar = blatResult.cigar
        self.svType['type'] = 'indel'
        contigCountTracker = svEvent.contig.get_contig_count_tracker()
        self.splitReadCount = [contigCountTracker.get_counts(x, x, 'indel') for x in blatResult.contigBreakpoints]
        self.filterValues.set_indel_values(blatResult, self.splitReadCount)

    def format_rearrangement_values(self, svEvent):
        """ """
        utils.log(self.loggingName, 'info', 'Resolving SVs call from blat results')
        # Sort the stored blat results by the number of matches to the reference sequence.
        blatResSorted = sorted(self.blatResults, key=lambda x: x[1])
        # brkpts = SVBreakpoints() # {'t':{'in_target':None, 'other':None }, 'formatted':[], 'r':[], 'q': [[0,0],[]], 'chrs':[], 'brkpt_str':[], 'tcoords':[], 'f': []}
        # res_values = {'target_breakpoints':[], 'align_cigar':[], 'sv_type':'', 'strands':[], 'mismatches':[], 'repeat_matching':[], 'anno_genes': [], 'disc_read_count': 0}
        resultValid = {'valid': True, 'repeatValid': True}
        maxRepeat = 0.0

        self.repeatMatching = []
        self.genes = []
        self.alignCigar = []
        self.strands = []
        self.mismatches = []
        for i, blatResultTuple in enumerate(blatResSorted):
            blatResult = blatResultTuple[1]
            resultValid['valid'] = resultValid['valid'] and blatResult.valid
            resultValid['repeatValid'] = resultValid['repeatValid'] and (blatResult.rep_man.simple_rep_overlap > 75.0)
            maxRepeat = max(maxRepeat, blatResult.repeat_overlap)
            self.repeatMatching.append(":".join([str(blatResult.repeat_overlap), str(blatResult.get_nmatch_total()), str(round(blatResult.mean_cov, 3))]))
            self.genes.append(blatResult.get_gene_anno())
            self.alignCigar.append(blatResult.cigar)
            self.strands.append(blatResult.strand)
            self.mismatches.append(blatResult.get_nmatches('mismatch'))
            svEvent.update_brkpt_info(blatResult, i, i == (len(blatResSorted) - 1))

        if svEvent.diff_chr():
            # translocation event
            svEvent.set_brkpt_counts('trl')
            self.discReadCount = svEvent.get_disc_read_count()
            self.svType['type'] = ['trl']
            self.filterValues.set_trl_values()
            # self.targetBreakpoints = svEvent.get_brkpt_str() # brkpts['brkpt_str']
            # self.splitReadCount = svEvent.get_splitread_count() # brkpt_counts['b']
        else:
            svEvent.set_brkpt_counts('rearr')
            rearrType, discReadCount = svEvent.define_rearr()
            self.svType = 'rearrangement'
            if rearrType != 'rearrangement':
                self.svType['subtype'] = rearrType
            self.discReadCount = discReadCount
            self.genes = list(set(self.genes))
        self.targetBreakpoints = svEvent.get_brkpt_str()
        self.splitReadCount = svEvent.get_splitread_count()

        # result = None
        # self.blatResultsSorted = sorted(self.blatResultsSorted, key=lambda br: br[1])
        # if not self.multiple_genes(brkpts['chrs'], brkpts['r'], res_values['anno_genes']):
        #     brkpt_counts, brkpt_kmers, brkpt_rep_filt = self.get_brkpt_counts_filt(brkpts, 'rearr')
        #     rearr_type, disc_read_support = self.define_rearr(brkpts['r'], res_values['strands'], brkpts['tcoords'], disc_reads)
        #     if not self.filter_rearr(query_region, params, brkpts['r'], brkpt_counts, brkpt_kmers, rearr_type, disc_read_support):
        #         res_values['sv_type'] = 'rearrangement'
        #         if rearr_type != 'rearrangement':
        #             res_values['sv_subtype'] = rearr_type
        #         res_values['disc_read_count'] = disc_read_support
        #         res_values['anno_genes'] = list(set(res_values['anno_genes']))
        #         res_values['target_breakpoints'] = brkpts['brkpt_str']
        #         res_values['split_read_count'] = brkpt_counts['b']
        #         if 'rearrangement' in params.opts['var_filter']:
        #             result = self.format_result(res_values)
        # elif max(self.contig_rcounts.others) >= params.get_sr_thresh('trl'):
        #     brkpt_counts, brkpt_kmers, brkpt_rep_filt = self.get_brkpt_counts_filt(brkpts, 'trl')
        #     disc_read_count = self.check_disc_reads(brkpts['t'], query_region, disc_reads['disc'])
        #     if not self.filter_trl(br_valid, query_region, params, brkpt_counts, brkpt_kmers, disc_read_count, res_values['anno_genes'], max_repeat, brkpt_rep_filt):
        #         res_values['disc_read_count'] = disc_read_count
        #         res_values['sv_type'] = ['trl']
        #         res_values['target_breakpoints'] = brkpts['brkpt_str']
        #         res_values['split_read_count'] = brkpt_counts['b']
        #         if 'trl' in params.opts['var_filter']:
        #             result = self.format_result(res_values)
        # return result

    def format_result(self, values):
        res_lst = []
        if values:
            for v in values:
                if not isinstance(values[v], list):
                    values[v] = [values[v]]
                self.resultValues[v] = ",".join([str(x) for x in values[v]])
        if self.resultValues['sv_subtype']:
            self.resultValues['sv_type'] += '_' + self.resultValues['sv_subtype']
        return self.get_values()

    def get_brkpt_coverages(self):
        brkpts = []
        tbp = self.resultValues['target_breakpoints']
        if self.resultValues['target_breakpoints'].find("(") > -1:
            tbp = self.resultValues['target_breakpoints'].split()[0]
        tbp = tbp.split(',')
        for bp in tbp:
            chrom, locs = bp.split(':')
            chrom = chrom.replace('chr', '')
            ll = locs.split('-')
            if len(ll) > 1:
                brkpts.append((chrom, int(ll[0]), int(ll[0]) + 1))
                brkpts.append((chrom, int(ll[1]), int(ll[1]) + 1))
            else:
                brkpts.append((chrom, int(ll[0]), int(ll[0]) + 1))

        bamfile = Samfile(self.sample_bam, 'rb')

        covs = [0] * len(brkpts)
        bp_index = 0
        for bp in brkpts:
            cov = 0
            c, s, e = bp
            areads = bamfile.fetch(str(c), s, e)

            for aread in areads:
                if aread.is_duplicate or aread.is_qcfail or aread.is_unmapped or aread.mapq < 10:
                    continue
                cov += 1
            covs[bp_index] = cov
            bp_index += 1
        return ",".join([str(x) for x in covs])

    def get_values(self):
        lst = ['anno_genes', 'target_breakpoints', 'align_cigar', 'mismatches', 'strands', 'repeat_matching', 'sv_type', 'split_read_count', 'nkmers', 'disc_read_count', 'breakpoint_coverages', 'contig_id', 'contig_seq']
        out_lst = []
        for l in lst:
            if self.resultValues[l] == 'trl':
                self.resultValues[l] = 'rearrangement'
            out_lst.append(str(self.resultValues[l]))
            if l == 'target_breakpoints':
                self.resultValues['breakpoint_coverages'] = self.get_brkpt_coverages()
        return out_lst


class SVBreakpoints:
    def __init__(self):
        self.t = {'in_target': None, 'other': None}
        self.formatted = []
        self.r = []
        self.q = [[0, 0], []]
        self.chrs = []
        self.brkptStr = []
        self.tcoords = []
        self.f = []
        self.counts = {'n': [], 'd': [], 'b': []}
        self.kmers = []

    def update_brkpt_info(self, br, i, last_iter):
        ts, te = br.get_coords('hit')
        qs, qe = br.get_coords('query')
        target_key = 'in_target' if br.in_target else 'other'
        brkpt_d['chrs'].append(br.get_name('hit'))
        brkpt_d['tcoords'].append((ts, te))
        tbrkpt = []
        filt_rep_start = None
        if i == 0:
            brkpt_d['q'][0] = [max(0, qs - 1), qe]
            brkpt_d['q'][1].append([qe, qe - brkpt_d['q'][0][0], None])
            tbrkpt = [te]
            filt_rep_start = br.filter_reps_edges[0]
            if br.strand == '-':
                tbrkpt = [ts]
                filt_rep_start = br.filter_reps_edges[0]
        elif last_iter:
            brkpt_d['q'][1][-1][2] = qe - brkpt_d['q'][1][-1][0]
            brkpt_d['q'][1].append([qs, qs - brkpt_d['q'][0][0], qe - qs])
            tbrkpt = [ts]
            filt_rep_start = br.filter_reps_edges[0]
            if br.strand == '-':
                tbrkpt = [te]
                filt_rep_start = br.filter_reps_edges[1]
        else:
            brkpt_d['q'][1][-1][2] = qe - brkpt_d['q'][1][-1][1]
            brkpt_d['q'][1].append([qs, qs-brkpt_d['q'][0][0], qe - qs])
            brkpt_d['q'][1].append([qe, qe-qs, None])
            brkpt_d['q'][0] = [qs, qe]
            tbrkpt = [ts, te]
            if br.strand == '-':
                filt_rep_start = br.filter_reps_edges[1]
                tbrkpt = [te, ts]

        brkpt_d['brkpt_str'].append('chr'+str(br.get_name('hit')) + ":" + "-".join([str(x) for x in tbrkpt]))
        brkpt_d['r'].extend(tbrkpt)
        brkpt_d['f'].append(filt_rep_start)
        brkpt_d['t'][target_key] = (br.get_name('hit'), tbrkpt[0])
        brkpt_d['formatted'].append('chr'+str(br.get_name('hit')) + ":" + "-".join([str(x) for x in tbrkpt]))
        return brkpt_d

    def diff_chr(self):
        if len(set(self.chrs)) == 1:
            return False
        else:
            return True

    def get_target_brkpt(self, key):
        """ """
        return self.target[key]

    def get_brkpt_str(self):
        """ """
        return self.brkptStr

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
            self.kmers.append(self.contig.get_kmer_locs()[qb[0]])
            # brkpt_rep_filt = brkpt_rep_filt or (comp_vec[qb[0]] < (avg_comp / 2))
            utils.log(self.loggingName, 'debug', 'Read count around breakpoint %d : %s' % (qb[0], ",".join([str(x) for x in bc])))
        utils.log(self.loggingName, 'debug', 'Kmer count around breakpoints %s' % (",".join([str(x) for x in self.kmers])))
        # brkpt_rep_filt = brkpt_rep_filt or (len(filter(lambda x: x, brkpts['f'])) > 0)
        # return brkpt_counts, brkpt_kmers, brkpt_rep_filt


class SVEvent:
    def __init__(self, blatResult, contig, svType):
        self.loggingName = 'breakmer.caller.sv_caller'
        self.svType = svType
        self.events = []
        self.blatResults = []
        self.blatResultsSorted = []
        self.sample_bam = sample_bam
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
        self.blatResultSorted.append((blatResult, blatResult.get_nmatch_total()))

    def result_valid(self):
        valid = False
        if (len(self.blatResults) > 1) and self.in_target:
            valid = True
            # nMissingQueryCoverage = len(filter(lambda y: y, map(lambda x: x == 0, self.queryCoverage)))
            # if nMissingQueryCoverage < self.meta_dict['params'].get_min_segment_length('trl'):
            #     valid = True
        return valid

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

    def get_brkpt_str(self):
        """ """
        return self.brkpts.get_brkpt_str()

    def get_splitread_count(self):
        """ """
        return self.brkpts.get_splitread_count()

    # def get_brkpt_info(self, br, brkpt_d, i, last_iter):
    #     ts, te = br.get_coords('hit')
    #     qs, qe = br.get_coords('query')
    #     target_key = 'in_target' if br.in_target else 'other'
    #     brkpt_d['chrs'].append(br.get_name('hit'))
    #     brkpt_d['tcoords'].append((ts, te))
    #     tbrkpt = []
    #     filt_rep_start = None
    #     if i == 0:
    #         brkpt_d['q'][0] = [max(0, qs - 1), qe]
    #         brkpt_d['q'][1].append([qe, qe - brkpt_d['q'][0][0], None])
    #         tbrkpt = [te]
    #         filt_rep_start = br.filter_reps_edges[0]
    #         if br.strand == '-':
    #             tbrkpt = [ts]
    #             filt_rep_start = br.filter_reps_edges[0]
    #     elif last_iter:
    #         brkpt_d['q'][1][-1][2] = qe - brkpt_d['q'][1][-1][0]
    #         brkpt_d['q'][1].append([qs, qs - brkpt_d['q'][0][0], qe - qs])
    #         tbrkpt = [ts]
    #         filt_rep_start = br.filter_reps_edges[0]
    #         if br.strand == '-' :
    #             tbrkpt = [te]
    #             filt_rep_start = br.filter_reps_edges[1]
    #     else :
    #         brkpt_d['q'][1][-1][2] = qe - brkpt_d['q'][1][-1][1]
    #         brkpt_d['q'][1].append([qs,qs-brkpt_d['q'][0][0],qe-qs])
    #         brkpt_d['q'][1].append([qe, qe-qs, None])
    #         brkpt_d['q'][0] = [qs, qe]
    #         tbrkpt = [ts, te]
    #         if br.strand == "-" :
    #             filt_rep_start = br.filter_reps_edges[1]
    #             tbrkpt = [te, ts]

    #     brkpt_d['brkpt_str'].append('chr'+str(br.get_name('hit')) + ":" + "-".join([str(x) for x in tbrkpt]))
    #     brkpt_d['r'].extend(tbrkpt)
    #     brkpt_d['f'].append(filt_rep_start)
    #     brkpt_d['t'][target_key] = (br.get_name('hit'),tbrkpt[0])
    #     brkpt_d['formatted'].append( 'chr'+str(br.get_name('hit')) + ":" + "-".join([str(x) for x in tbrkpt]))
    #     return brkpt_d

    # def get_svs_result(self, query_region, params, disc_reads):
        # self.logger.info('Resolving SVs call from blat results')
        # blat_res = self.blatResults
        # blat_res_sorted = sorted(blat_res, key=lambda blat_res: blat_res[0])
        # brkpts = {'t':{'in_target':None, 'other':None }, 'formatted':[], 'r':[], 'q': [[0,0],[]], 'chrs':[], 'brkpt_str':[], 'tcoords':[], 'f': []}
        # res_values = {'target_breakpoints':[], 'align_cigar':[], 'sv_type':'', 'strands':[], 'mismatches':[], 'repeat_matching':[], 'anno_genes': [], 'disc_read_count': 0}
        # br_valid = [True, True]
        # max_repeat = 0.0

        # for i in range(len(blat_res_sorted)):
        #     br = blat_res_sorted[i][1]
        #     br_valid[0] = br_valid[0] and br.valid
        #     br_valid[1] = br_valid[1] and (br.rep_man.simple_rep_overlap > 75.0)
        #     max_repeat = max(max_repeat, br.repeat_overlap)
        #     res_values['repeat_matching'].append(":".join([str(br.repeat_overlap), str(br.get_nmatch_total()), str(round(br.mean_cov, 3))]))
        #     res_values['anno_genes'].append(br.get_gene_anno())
        #     res_values['align_cigar'].append(br.cigar)
        #     res_values['strands'].append(br.strand)
        #     res_values['mismatches'].append(br.get_nmatches('mis'))
        #     brkpts = self.get_brkpt_info(br, brkpts, i, i == (len(blat_res_sorted) - 1))

        # result = None
        # self.blatResultsSorted = sorted(self.blatResultsSorted, key=lambda br: br[1])
        # if not self.multiple_genes(brkpts['chrs'], brkpts['r'], res_values['anno_genes']):
        #     brkpt_counts, brkpt_kmers, brkpt_rep_filt = self.get_brkpt_counts_filt(brkpts, 'rearr')
        #     rearr_type, disc_read_support = self.define_rearr(brkpts['r'], res_values['strands'], brkpts['tcoords'], disc_reads)
        #     if not self.filter_rearr(query_region, params, brkpts['r'], brkpt_counts, brkpt_kmers, rearr_type, disc_read_support):
        #         res_values['sv_type'] = 'rearrangement'
        #         if rearr_type != 'rearrangement':
        #             res_values['sv_subtype'] = rearr_type
        #         res_values['disc_read_count'] = disc_read_support
        #         res_values['anno_genes'] = list(set(res_values['anno_genes']))
        #         res_values['target_breakpoints'] = brkpts['brkpt_str']
        #         res_values['split_read_count'] = brkpt_counts['b']
        #         if 'rearrangement' in params.opts['var_filter']:
        #             result = self.format_result(res_values)
        # elif max(self.contig_rcounts.others) >= params.get_sr_thresh('trl'):
        #     brkpt_counts, brkpt_kmers, brkpt_rep_filt = self.get_brkpt_counts_filt(brkpts, 'trl')
        #     disc_read_count = self.check_disc_reads(brkpts['t'], query_region, disc_reads['disc'])
        #     if not self.filter_trl(br_valid, query_region, params, brkpt_counts, brkpt_kmers, disc_read_count, res_values['anno_genes'], max_repeat, brkpt_rep_filt):
        #         res_values['disc_read_count'] = disc_read_count
        #         res_values['sv_type'] = ['trl']
        #         res_values['target_breakpoints'] = brkpts['brkpt_str']
        #         res_values['split_read_count'] = brkpt_counts['b']
        #         if 'trl' in params.opts['var_filter']:
        #             result = self.format_result(res_values)
        # return result

    def set_brkpt_counts(self, svType):
        """ """
        self.brkpts.set_counts(svType)

    # def call_rearr(self):
    #         rearr_type, disc_read_support = self.define_rearr(brkpts['r'], res_values['strands'], brkpts['tcoords'], disc_reads)
    #         if not self.filter_rearr(query_region, params, brkpts['r'], brkpt_counts, brkpt_kmers, rearr_type, disc_read_support):
    #             res_values['sv_type'] = 'rearrangement'
    #             if rearr_type != 'rearrangement' : res_values['sv_subtype'] = rearr_type
    #             res_values['disc_read_count'] = disc_read_support
    #             res_values['anno_genes'] = list(set(res_values['anno_genes']))
    #             res_values['target_breakpoints'] = brkpts['brkpt_str']
    #             res_values['split_read_count'] = brkpt_counts['b']
    #             if 'rearrangement' in params.opts['var_filter']:
    #                 result = self.format_result(res_values)

    # def call_trl(self):
    #         disc_read_count = self.check_disc_reads(brkpts['t'], query_region, disc_reads['disc'])
    #         if not self.filter_trl(br_valid, query_region, params, brkpt_counts, brkpt_kmers, disc_read_count, res_values['anno_genes'], max_repeat, brkpt_rep_filt):
    #             res_values['disc_read_count'] = disc_read_count
    #             res_values['sv_type'] = ['trl']
    #             res_values['target_breakpoints'] = brkpts['brkpt_str']
    #             res_values['split_read_count'] = brkpt_counts['b']
    #             if 'trl' in params.opts['var_filter']:
    #                 result = self.format_result(res_values)

    def check_overlap(self, coord1, coord2):
        contained = False
        if coord1[0] >= coord2[0] and coord1[1] <= coord2[1]:
            contained = True
        elif coord2[0] >= coord1[0] and coord2[1] <= coord1[1]:
            contained = True
        return contained

    def define_rearr(self):
        vrt = self.contig.get_variant_read_tracker()
        strands = self.svResultValues.strands
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
                    for readPair in vrt.inv:
                        r1p, r2p, r1s, r2s, qname = readPair
                        if r1s == 1 and r2s == 1:
                            if (r1p <= brkpts[0]) and (r2p <= brkpts[1] and r2p >= brkpts[0]):
                                rs += 1
                        else:
                            if (r1p <= brkpts[1] and r1p >= brkpts[0]) and r2p >= brkpts[1]:
                                rs += 1
                elif (strands[0] == '+' and strands[1] == '+') and (brkpts[0] > brkpts[1]):
                    utils.log(self.loggingName, 'debug', 'Tandem duplication event identified.')
                    hit = True
                    svType = 'tandem_dup'
                    # Tandem dup
                    for readPair in vrt.td:
                        r1p, r2p, r1s, r2s, qname = readPair
                        if (r1p <= brkpts[0] and r1p >= brkpts[1]) and ():
                            rs += 1
        if not hit:
            utils.log(self.loggingName, 'debug', 'Not inversion or tandem dup, checking for odd read pairs around breakpoints')
            rs = [0] * len(brkpts)
            for i in range(len(brkpts)):
                b = brkpts[i]
                for readPair in vrt.other:
                    r1p, r2p, r1s, r2s, qname = readPair
                    if abs(r1p - b) <= 300 or abs(r2p - b) <= 300:
                        utils.log(self.loggingName, 'debug', 'Adding read support from read %s, with strands %d, %d and positions %d, %d for breakpoint at %d' % (qname, r1s, s, r1p, r2p, b))
                        rs[i] += 1
            rs = max(rs)
        return svType, rs

    # def filter_rearr(self, query_region, params, brkpts, brkpt_counts, brkpt_kmers, rearr_type, disc_read_count):
    #     in_ff, span_ff = filter_by_feature(brkpts, query_region, params.opts['keep_intron_vars'])
    #     filter = (min(brkpt_counts['n']) < params.get_sr_thresh('rearrangement')) or self.blatResultsSorted[0][1] < params.get_min_segment_length('rearr') or (in_ff and span_ff) or (disc_read_count < 1) or (rearr_type == 'rearrangement') or (min(brkpt_kmers) == 0)
    #     self.logger.info('Check filter for rearrangement')
    #     self.logger.info('Filter by feature for being in exon (%r) or spanning exon (%r)'%(in_ff, span_ff))
    #     self.logger.info('Split read threshold %d, breakpoint read counts %d'%(min(brkpt_counts['n']),params.get_sr_thresh('rearrangement')))
    #     self.logger.info('Minimum segment length observed (%d) meets threshold (%d)'%(self.blatResultsSorted[0][1], params.get_min_segment_length('rearr')))
    #     self.logger.info('Minimum discordant read pairs for rearrangement (%d)'%(disc_read_count))
    #     return filter

#     def filter_trl(self, br_valid, query_region, params, brkpt_counts, brkpt_kmers, disc_read_count, anno_genes, max_repeat, rep_filt) :
#         filter = br_valid[1] or (max(brkpt_counts['d']) < params.get_sr_thresh('trl')) #or not br_valid[0]
#         self.logger.debug('Check translocation filter')
#         self.logger.debug('All blat result segments are within annotated or pre-specified regions %r'%br_valid[0])
#         self.logger.debug('All blat result segments are within simple repeat regions that cover > 75.0 percent of the segment %r'%br_valid[1])
#         self.logger.debug('The maximum read count support around breakpoints %d meets split read threshold %d'%(max(brkpt_counts['d']),params.get_sr_thresh('trl')))
#         self.logger.debug('The minimum number of kmers at breakpoints %d'%min(brkpt_kmers))
#         self.logger.debug('The maximum repeat overlap by a blat result: %f'%max_repeat)
#         if not filter :
#             self.logger.debug('Filter %r, checking discordant read counts %d'%(filter, disc_read_count)) 
#             if disc_read_count < 2 :
# #        print 'Filter due to repeat', rep_filt
#                 if (self.blatResultsSorted[0][1] < params.get_min_segment_length('trl')) or (min(brkpt_counts['n']) < params.get_sr_thresh('trl')) or (min(brkpt_kmers)==0) or rep_filt :
#                     self.logger.debug('Shortest segment is < %d bp with %d discordant reads. Filtering.'%(params.get_min_segment_length('trl'), disc_read_count))
#                     self.logger.debug('The minimum read count support for breakpoints %d meets split read threshold %d'%(min(brkpt_counts['n']),params.get_sr_thresh('trl')))
#                     self.logger.debug('The minimum number of kmers at breakpoints %d'%min(brkpt_kmers))
#                     filter = True
#                 elif disc_read_count == 0 : 
#                     # Check a number of metrics for shortest blat segment
#                     br_qs = self.blatResultsSorted[0][0].qstart()
#                     br_qe = self.blatResultsSorted[0][0].qend()
#                     low_complexity = self.minseq_complexity(self.contig_seq[br_qs:br_qe],3) < 25.0 # Complexity of blat segment
#                     missing_qcov = self.missing_query_coverage() > 5.0
#                     short = self.blatResultsSorted[0][1] <= round(float(len(self.contig_seq))/float(4.0))
#                     self.logger.debug('Checking length of shortest sequence, considered too short %r, %d, %f'%(short, self.blatResultsSorted[0][1], round(float(len(self.contig_seq))/float(4.0))) )
#                     overlap = max(self.blatResultsSorted[0][0].seg_overlap) > 5
#                     gaps_exist = max(self.blatResultsSorted[0][0].gaps['query'][0], self.blatResultsSorted[0][0].gaps['hit'][0]) > 0
#                     low_uniqueness = self.check_uniqueness()
#                     intergenic_regions = 'intergenic' in anno_genes
#                     read_strand_bias = self.check_read_strands()
#                     check_values = [low_complexity, missing_qcov, short, overlap, gaps_exist, low_uniqueness, read_strand_bias, intergenic_regions]
#                     self.logger.debug('Discordant read count of 0 checks %s'%(",".join([str(x) for x in check_values])))
#                     num_checks = 0
#                     for check in check_values :
#                         if check : 
#                             num_checks += 1
#                     if num_checks > 1 : 
#                         self.logger.info('Two or more filter checks, setting filtering to true for contig')
#                         filter = True
#         return filter

    # def check_uniqueness(self):
    #     low_unique = False
    #     for br_vals in self.blatResultsSorted :
    #         if not br_vals[0].in_target :
    #             if br_vals[0].mean_cov > 4 : low_unique = True
    #         else :
    #             if br_vals[0].mean_cov > 10 : low_unique = True
    #     return low_unique

    # def check_read_strands(self):
    #     same_strand = False
    #     strands = []
    #     for read in self.contig_reads :
    #         strand = read.id.split("/")[1] 
    #         strands.append(strand)
    #     if len(set(strands)) == 1 :
    #         same_strand = True
    #     self.logger.debug('Checking read strands for contig reads %s'%(",".join([read.id for read in self.contig_reads])))
    #     self.logger.debug('Reads are on same strand: %r'%same_strand)
    #     return same_strand

    # def minseq_complexity(self, seq, N) :
    #     self.logger.debug('Checking sequence complexity of blat result segment %s using %d-mers'%(seq,N))
    #     nmers = {}
    #     total_possible = len(seq) - 2
    #     for i in range(len(seq) - (N - 1)) :
    #         nmers[str(seq[i:i+N]).upper()] = True
    #     complexity = round((float(len(nmers))/float(total_possible))*100,4)
    #     self.logger.debug('Complexity measure %f, based on %d unique %d-mers observed out of a total of %d %d-mers possible'%(complexity,len(nmers),N,total_possible,N))
    #     return complexity

    # def missing_query_coverage(self) :
    #     missing_cov = 0
    #     for i in self.queryCoverage :
    #         if i == 0 :
    #             missing_cov += 1
    #         else :
    #             break

    #     for i in reversed(self.queryCoverage) :
    #         if i == 0 :
    #             missing_cov += 1
    #         else : 
    #             break

    #     perc_missing = round((float(missing_cov)/float(len(self.contig_seq)))*100, 4)
    #     self.logger.debug('Calculated %f missing coverage of blat query sequence at beginning and end'%perc_missing)
    #     return perc_missing

    # def multiple_genes(self, chrs, brkpts, anno_genes):
    #     mult_genes = True
    #     if len(set(anno_genes)) == 1:
    #         self.logger.debug('One annotated gene among SVs breakpoints: %s' % (",".join(anno_genes)))
    #         mult_genes = False
    #     elif self.dup_gene_names(anno_genes) and len(set(chrs)) == 1 and ((max(brkpts) - min(brkpts)) < 10000):
    #         self.logger.debug('Anno genes are not the same, but similar and breakpoints are less than 10Kb from each other %s' % (",".join(anno_genes)))
    #         mult_genes = False
    #     self.logger.debug('Test whether SVs breakpoints are in multiple genes %r' % mult_genes)
    #     return mult_genes

    # def dup_gene_names(self, anno_genes):
    #     found_dup = False
    #     for i in range(len(anno_genes)-1) :
    #         g1 = anno_genes[i]
    #         for g2 in anno_genes[(i+1):] :
    #             if (g1.find(g2) > -1) or (g2.find(g1) > -1) :
    #                 found_dup = True
    #     return found_dup 

    # def check_disc_reads(self, brkpts, query_region, disc_reads):
    #     disc_read_count = 0
    #     if brkpts['other'][0] in disc_reads:
    #         for p1, p2 in disc_reads[brkpts['other'][0]]:
    #             d1 = abs(p1 - brkpts['in_target'][1])
    #             d2 = abs(p2 - brkpts['other'][1])
    #             if d1 <= 1000 and d2 <= 1000:
    #                 disc_read_count += 1
    #     return disc_read_count


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
        if not self.realignment.has_results():
            self.logger.info('No blat results file exists, no calls for %s.' % self.contig.get_id())
        else:
            self.logger.info('Making variant calls from blat results %s' % self.realignment.get_result_fn())
            if self.check_indels():
                # self.result = self.realignment.get_indel_result()
                self.svEvent.format_indel_values()
            elif self.check_svs():
                # self.result = self.realignment.get_svs_result()
                self.svEvent.format_rearr_values()
        # Format the result into a
        return self.svEvent

    def check_indels(self):
        """ """
        hasIndel = False
        blatResults = self.realignment.get_blat_results()
        for i, blatResult in enumerate(blatResults):
            if i == 0 and blatResult.check_indel(params.get_param('indel_size'), len(blatResults)):
                hasIndel = True
                utils.log(self.loggingName, 'info', 'Contig has indel, returning %r' % has_indel)
                self.svEvent = SVEvent(blatResult, self.contig, 'indel')
                return has_indel
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
                utils.log(self.loggingName, 'debug', 'Blat result with start %d, end %d, chrom %s' % (qs, qe, blatResult.get_name('hit')))
                gaps = self.iter_gaps(gaps, self.clippedQs[i], i)
                if self.svEvent.qlen > mergedClip[0]:
                    mergedClip = [self.svEvent.qlen, self.svEvent]
            self.svEvent = mergedClip[1]
        else:
            utils.log(self.loggingName, 'info', 'There are no more than 1 clipped blat results, not continuing with SVs calling.')
        if self.svEvent and self.svEvent.result_valid():
            return True
        else:
            return False

    def iter_gaps(self, gaps, clippedQuerySeqVals, iterIdx):
        """ """
        new_gaps = []
        qs, qe, blatResult, idx = clippedQuerySeqVals
        hit = False
        for gap in gaps:
            gs, ge = gap
            utils.log(self.loggingName, 'debug', 'Gap coords %d, %d' % (gs, ge))
            if (qs >= gs and qs <= ge) or (qe <= ge and qe >= gs):
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

    def check_add_br(self, qs, qe, gs, ge, br):
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
