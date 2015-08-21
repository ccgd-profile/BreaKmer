#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
import subprocess
import breakmer.realignment.blat_result as blat_result
import breakmer.utils as utils

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


class AlignParams:
    """
    """

    def __init__(self, params, targetRefFns):
        self.program = {'target': 'blat', 'genome': 'blat'}
        self.extension = {'target': 'psl', 'genome': 'psl'}
        self.binary = {'target': None, 'genome': None}
        self.binaryParams = {'target': None, 'genome': None}
        self.ref = {'target': None, 'genome': None}
        self.set_values(params, targetRefFns)

    def set_values(self, params, targetRefFns):
        """
        """
        self.binary['target'] = params.get_param('blat')
        blast = params.get_param('blast')
        if blast is not None:
            self.program['target'] = 'blast'
            self.binary['target'] = blast
            self.extension['target'] = 'txt'

        self.binary['genome'] = params.get_param('gfclient')
        self.binaryParams['genome'] = {'hostname': params.get_param('blat_hostname'),
                                       'port': int(params.get_param('blat_port'))}
        # Use the forward sequence for blatting targeted sequences
        self.ref['target'] = targetRefFns[0]
        self.ref['genome'] = params.get_param('reference_fasta_dir')

    def get_values(self, type):
        return (self.program[type], self.extension[type], self.binary[type], self.binaryParams[type], self.ref[type])


class RealignManager:
    """
    """

    def __init__(self, params, targetRefFns):
        self.realignment = None
        self.alignParams = AlignParams(params, targetRefFns)

    def realign(self, contig):
        """
        """
        if not contig.has_fa_fn():
            return

        self.realignment = Realignment(contig)
        if not self.realignment.align(self.alignParams.get_values('target'), 'target'):
            # Log
            return
        # print 'realigner.py realign() check target_aligned()', self.realignment.target_aligned()
        if not self.realignment.target_aligned():
            self.realignment.align(self.alignParams.get_values('genome'), 'genome')
        else:
            if self.realignment.targetHit and self.alignParams.get_values('target')[0] == 'blast':
                self.realignment.check_record_merge()

    def get_result_fn(self):
        resultFn = None
        if self.realignment.has_results():
            resultFn = self.realignment.get_result_fn()
        return resultFn

    def has_results(self):
        """ """
        return self.realignment.has_results()

    def get_blat_results(self):
        """
        """
        return self.realignment.get_blat_results()

    def store_clipped_queryseq(self, blatResultValues):
        """
        """
        self.realignment.store_clipped_queryseq(blatResultValues)

    def get_qsize(self):
        """ """
        return self.realignment.results.querySize


class Realignment:
    """
    """
    def __init__(self, contig):
        self.loggingName = 'breakmer.realignment.realigner'
        self.scope = None
        self.results = None
        self.targetHit = False
        self.resultFn = None
        self.alignParams = None
        self.contig = contig

    def align(self, alignParams, scope):
        """
        """
        self.alignParams = alignParams
        alignProgram, alignExt, alignBinary, binaryParams, alignRef = self.alignParams
        self.scope = scope

        self.resultFn = os.path.join(self.contig.get_path(), '%s_res.%s.%s' % (alignProgram, scope, alignExt))
        utils.log(self.loggingName, 'info', 'Running realignment with %s, storing results in %s' % (alignProgram, self.resultFn))
#        self.results = AlignResults(alignProgram, scope, resultFn)

        cmd = ''
        if alignProgram == 'blast':
            cmd = "%s -task 'blastn-short' -db %s -query %s -evalue 0.01 -out %s -outfmt '7 qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore gaps sstrand qseq sseq'" % (alignBinary, alignRef, self.contig.meta.fa_fn, self.resultFn)
        elif alignProgram == 'blat':
            if scope == 'genome':
                # all blat server
                cmd = '%s -t=dna -q=dna -out=psl -minScore=20 -nohead %s %d %s %s %s' % (alignBinary, binaryParams['hostname'], binaryParams['port'], alignRef, self.contig.meta.fa_fn, self.resultFn)
            elif scope == 'target':
                # target
                cmd = '%s -t=dna -q=dna -out=psl -minScore=20 -stepSize=10 -minMatch=2 -repeats=lower -noHead %s %s %s' % (alignBinary, alignRef, self.contig.meta.fa_fn, self.resultFn)

        utils.log(self.loggingName, 'info', 'Realignment system command %s' % cmd)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, errors = p.communicate()
        utils.log(self.loggingName, 'info', 'Realignment output file %s' % self.resultFn)
        if errors != '':
            utils.log(self.loggingName, 'info', 'Realignment errors %s' % errors)

        if not os.path.isfile(self.resultFn):
            return False
        else:
            self.results = AlignResults(alignProgram, scope, self.resultFn, self.contig, alignRef)
            # print self.results.resultFn
            return True

    def get_result_fn(self):
        """ """
        if self.results is not None:
            return self.resultFn

    def has_results(self):
        """ """
        if self.results is not None:
            return True
        else:
            return False

    def target_aligned(self):
        """
        """
        noAlignmentResults = False
        targetHit = False
        if self.results is None:
            noAlignmentResults = True
        else:
            self.results.modify_blat_result_file()
            if self.results.target_hit():
                self.targetHit = True
                utils.log(self.loggingName, 'debug', 'Top hit contains whole query sequence, indicating an indel variant within the target region.')

        # If there was a sufficient target hit or no alignment at all then return True
        # this effectively prevents a genome alignment.
        return self.targetHit or noAlignmentResults

    def get_blat_results(self):
        """
        """
        return self.results.sortedResults

    def store_clipped_queryseq(self, blatResultValues):
        self.clippedQs.append(blatResultValues)

    def check_record_merge(self):
        """ """
        # Check if need to merge indels from blast results
        if len(self.results.mergedRecords) > 0:
            # Re-write a blast-style file and re-procress with merged records.
            newResultFn = open(self.resultFn + '.merged_recs', 'w')
            linesOut = self.results.merge_records()
            for line in linesOut:
                newResultFn.write(line)
            newResultFn.close()
            self.resultFn = self.resultFn + '.merged_recs'
            alignProgram, alignExt, alignBinary, binaryParams, alignRef = self.alignParams
            self.results = AlignResults('blat', 'blast_target', self.resultFn, self.contig, alignRef)


class AlignResults:
    def __init__(self, program, scope, alignResultFn, contig, alignRefFn):
        self.loggingName = 'breakmer.realignment.realigner'
        self.resultFn = alignResultFn
        self.program = program
        self.scope = scope
        self.querySize = 0
        self.alignmentFreq = []
        self.nmismatches = 0
        self.ngaps = 0
        self.hasResults = True
        self.results = []
        self.sortedResults = []
        self.clippedQs = []
        self.contig = contig
        self.alignRefFn = alignRefFn
        self.mergedRecords = [] # List of tuples containing indices of realignment records that need to be merged.
        self.targetSegmentsSorted = None
        self.targetHit = False
        self.set_values()

    def set_values(self):
        """ """
        if not self.resultFn:
            self.hasResults = False
        elif len(open(self.resultFn, 'rU').readlines()) == 0:
            self.hasResults = False
        else:
            self.parse_result_file()

    def modify_blat_result_file(self):
        """ """
        blatFile = open(self.resultFn + '.mod', 'w')
        for result in self.results:
            blatFile.write(result.get_blat_output() + "\n")
        blatFile.close()
        self.resultFn = self.resultFn + '.mod'

    def target_hit(self):
        """ """
        if self.hasResults:
            # print 'realigner.py target_hit(), results', self.results
            cond1 = self.results[0].spans_query() and (self.ngaps > 0)
            cond2 = (len(self.results) == 1) and self.get_query_coverage() >= 90.0 and (self.ngaps > 0)
            cond3 = (len(self.results) > 1) and (self.get_query_coverage() >= 90.0)
            self.targetHit = cond1 or cond2 or cond3
        # print 'realigner.py target_hit() indelHit', indelHit, self.results[0].spans_query(), len(self.results), self.get_query_coverage()
        utils.log(self.loggingName, 'debug', 'Checking if query is a target hit or not %r' % self.targetHit)

        if self.targetHit:
            print len(self.results), self.get_query_coverage(), self.program
            if ((len(self.results) > 1) and (self.get_query_coverage() >= 90.0)) and (self.program == 'blast'):
                # Check for a gapped Blast result.
                segments = []
                for i, result in enumerate(self.results):
                    resultOverlap = 0
                    addSegment = True
                    for segment in segments:
                        overlapSeg = (result.qstart() >= segment[0] and result.qstart() <= segment[1]) or (result.qend() >= segment[0] and result.qend() <= segment[1])
                        containSeg = result.qstart() >= segment[0] and result.qend() >= segment[1]
                        withinSeg = result.qstart() >= segment[0] and result.qend() <= segment[1]
                        print i, result, containSeg, withinSeg, overlapSeg
                        if containSeg or withinSeg:
                            addSegment = False
                        elif overlapSeg:
                            if (result.qstart() >= segment[0] and result.qstart() <= segment[1]):
                                overlapBp = segment[1] - result.qstart()
                                print 'Overlapbp', overlapBp
                                if overlapBp < 20:
                                    resultOverlap += overlapBp
                                    addSegment = True
                            elif (result.qend() >= segment[0] and result.qend() <= segment[1]):
                                overlapBp = result.end() - segment[0]
                                print 'Overlapbp', overlapBp
                                if overlapBp < 20:
                                    resultOverlap += overlapBp
                                    addSegment = True
                    if addSegment and (result.get_query_span() - resultOverlap) > 20:
                        segments.append((result.qstart(), result.qend(), result))
                print segments
                self.targetSegmentsSorted = sorted(segments, key=lambda x: x[0])
                for i in range(1, len(self.targetSegmentsSorted)):
                    lResult = self.targetSegmentsSorted[i-1][2]
                    rResult = self.targetSegmentsSorted[i][2]
                    qgap = rResult.qstart() - lResult.qend()
                    tgap = rResult.tstart() - lResult.tend()
                    print rResult.tstart(), lResult.tend()
                    print rResult.qstart(), lResult.qend()
                    print 'tgap', tgap
                    print 'qgap', qgap
                    print (tgap < 0 and (abs(tgap) > abs(qgap)))
                    print (lResult.strand != rResult.strand)
                    if (tgap < 0 and (abs(tgap) > abs(qgap))) or (lResult.strand != rResult.strand):
                        # Tandem dup or inversion
                        break
                    else:
                        self.mergedRecords.append((i-1, i))
        print 'Merged record indices', self.mergedRecords
        return self.targetHit

    def parse_result_file(self):
        """ """
        refName = None
        offset = None
        if self.scope == 'target':
            # Need to reset the chrom name and coordinates for blat results.
            refName = self.contig.get_chr()
            offset = self.contig.get_target_start() - self.contig.get_target_buffer()

        for line in open(self.resultFn, 'r'):
            if line.find('#') > -1:
                continue
            line = line.strip()
            parsedResult = blat_result.BlatResult(line.split('\t'), refName, offset, self.program, self.alignRefFn, self.contig.seq, self.scope)
            parsedResult.in_target_region(self.contig.get_target_region_coordinates())
            # parsedBlatResult.set_gene_annotations(self.contig.get_target_region_coordinates(), self.contig.get_gene_annotations())
            # parsedBlatResult.set_repeats(self.contig.get_repeat_annotations())
            self.process_blat_result(parsedResult)
            self.results.append(parsedResult)
        # Update to use class attributes as sorting categories
        self.sortedResults = sorted(self.results, key=lambda x: (-x.alignScore, -x.perc_ident, x.get_total_num_gaps()))

        for i, blatResult in enumerate(self.sortedResults):
            blatResult.set_mean_cov(self.get_mean_cov(blatResult.qstart(), blatResult.qend()))

        if len(self.results) == 0:
            self.hasResults = False

    def merge_records(self):
        """ """
        mergedResults = []
        mapResults = {}
        for mergeIdx in self.mergedRecords:
            lResult = self.targetSegmentsSorted[mergeIdx[0]][2].resultValues
            rResult = self.targetSegmentsSorted[mergeIdx[1]][2].resultValues
            print 'left', lResult
            print 'right', rResult
            newMergedIdx = len(mergedResults)
            if mergeIdx[0] in mapResults:
                # Merge a result and a previously merged result.
                lResult = mergedResults[mapResults[mergedIdx[0]]]
                newMergedIdx = mapResults[mergedIdx[0]]
                mergedResults[newMergedIdx] = self.merge_record_fields(lResult, rResult)
            else:
                # Merge right and left values
                mergedResults.append(self.merge_record_fields(lResult, rResult))
            mapResults[mergeIdx[0]] = newMergedIdx
            mapResults[mergeIdx[1]] = newMergedIdx
        outputIdx = []
        lOuts = []
        for i, result in enumerate(self.targetSegmentsSorted):
            if i in mapResults:
                print 'New result', i, mergedResults[mapResults[i]]
                if mapResults[i] not in outputIdx:
                    outputIdx.append(mapResults[i])
                    lOuts.append(self.format_to_blat_output(mergedResults[mapResults[i]]))
            else:
                print 'Result', self.targetSegmentsSorted[i][2].resultValues
                lOuts.append(self.format_to_blat_output(self.targetSegmentsSorted[i][2].resultValues))
        print lOuts
        return lOuts

    def format_to_blat_output(self, resultValues):
        """ """
#         matches': int(values[0]),
# 'mismatches': int(values[1]),
# 'repmatches': int(values[2]),
# 'ncount': int(values[3]),
# 'qNumInsert': int(values[4]),
# 'qBaseInsert': int(values[5]),
# 'tNumInsert': int(values[6]),
# 'tBaseInsert': int(values[7]),
# 'strand': values[8],
# 'qName': values[9],
# 'qSize': int(values[10]),
# 'qStart': int(values[11]),
# 'qEnd': int(values[12]),
# 'tName': values[13].replace('chr', ''),
# 'tSize': int(values[14]),
# 'tStart': int(values[15]),
# 'tEnd': int(values[16]),
# 'blockCount': int(values[17]),
# 'blockSizes': values[18],
# 'qStarts': values[19],
# 'tStarts': values[20],
        outStr = []
        keys = ['matches', 'mismatches', 'repmatches', 'ncount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand', 'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts']
        for key in keys:
            outStr.append(resultValues[key])
        return '\t'.join([str(x) for x in outStr]) + '\n'

    def merge_record_fields(self, lResult, rResult):
        """ """
        tGap = rResult['tStart'] - lResult['tEnd']
        qGap = rResult['qStart'] - lResult['qEnd']
        lqEnd = lResult['qEnd']

        lResult['qStart'] = int(lResult['qStart']) - 1
        rResult['qStart'] = int(rResult['qStart']) - 1

        if tGap > 1:
            # Del
            newqStart = str(lResult['qEnd'] - 1)
            if qGap > 0:
                newqStart = str(rResult['qStart'])
            lResult['qStarts'] += newqStart + ','
            lResult['tStarts'] += str(rResult['tStart']) + ','
            lResult['tNumInsert'] += 1
            lResult['tBaseInsert'] += tGap
        if qGap > 1:
            # Ins
            newtStart = str(lResult['tEnd'] - 1)
            if tGap > 0:
                newtStart = str(rResult['tStart'])
            lResult['qStarts'] += str(rResult['qStart']) + ','
            lResult['tStarts'] += str(newtStart) + ','
            lResult['qNumInsert'] += 1
            lResult['qBaseInsert'] += qGap

        lResult['blockSizes'] = str(int(lResult['blockSizes'].rstrip(',')) - 1) + ','
        if qGap < 0:
            rResult['blockSizes'] = str(int(rResult['blockSizes'].rstrip(',')) + qGap) + ','
            rResult['matches'] += qGap - 1
        keys = ['matches', 'mismatches', 'repmatches', 'ncount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert']
        for key in keys:
            lResult[key] += rResult[key]
        lResult['qEnd'] = rResult['qEnd']
        lResult['tEnd'] = rResult['tEnd']
        lResult['blockCount'] += 1
        lResult['blockSizes'] += rResult['blockSizes']
        # lResult['qStarts'] += rResult['qStarts']
        # lResult['tStarts'] += rResult['tStarts']

        return lResult

    # def parse_blast_results(self):
    #     """ """
    #     refName = None
    #     offset = None
    #     if self.scope == 'target':
    #         # Need to reset the chrom name and coordinates for blat results.
    #         refName = self.contig.get_chr()
    #         offset = self.contig.get_target_start() - self.contig.get_target_buffer()

    #     for line in open(self.resultFn, 'r'):
    #         if line.find('#') > -1:
    #             continue
    #         line = line.strip()
    #         parsedBlastResult = blat_result.BlatResult(line.split('\t'), refName, offset)
    #         parsedBlastResult.in_target_region(self.contig.get_target_region_coordinates())
    #         # parsedBlatResult.set_gene_annotations(self.contig.get_target_region_coordinates(), self.contig.get_gene_annotations())
    #         # parsedBlatResult.set_repeats(self.contig.get_repeat_annotations())
    #         self.process_blat_result(parsedResult)
    #         self.results.append(parsedResult)
    #     # Update to use class attributes as sorting categories
    #     self.sortedResults = sorted(self.results, key=lambda x: (-x.alignScore, -x.perc_ident, x.get_total_num_gaps()))

    #     for i, blatResult in enumerate(self.sortedResults):
    #         blatResult.set_mean_cov(self.get_mean_cov(blatResult.qstart(), blatResult.qend()))

    # def parse_blat_results(self):
    #     """ """
    #     refName = None
    #     offset = None
    #     if self.scope == 'target':
    #         # Need to reset the chrom name and coordinates for blat results.
    #         refName = self.contig.get_chr()
    #         offset = self.contig.get_target_start() - self.contig.get_target_buffer()

    #     for line in open(self.resultFn, 'r'):
    #         line = line.strip()
    #         parsedBlatResult = blat_result.BlatResult(line.split('\t'), refName, offset)
    #         parsedBlatResult.in_target_region(self.contig.get_target_region_coordinates())
    #         # parsedBlatResult.set_gene_annotations(self.contig.get_target_region_coordinates(), self.contig.get_gene_annotations())
    #         # parsedBlatResult.set_repeats(self.contig.get_repeat_annotations())
    #         self.process_blat_result(parsedBlatResult)
    #         self.results.append(parsedBlatResult)
    #     # Update to use class attributes as sorting categories
    #     self.sortedResults = sorted(self.results, key=lambda x: (-x.alignScore, -x.perc_ident, x.get_total_num_gaps()))

    #     for i, blatResult in enumerate(self.sortedResults):
    #         blatResult.set_mean_cov(self.get_mean_cov(blatResult.qstart(), blatResult.qend()))

    def process_blat_result(self, blatResultObj):
        """Summarize metrics from all alignments.
        """
        self.nmismatches += blatResultObj.get_nmatches('mismatch')
        self.ngaps += blatResultObj.get_total_num_gaps()
        if not self.querySize:
            self.querySize = blatResultObj.get_seq_size('query')
            self.alignmentFreq = [0] * self.querySize
        for i in range(blatResultObj.qstart(), blatResultObj.qend()):
            self.alignmentFreq[i] += 1

    def get_query_coverage(self):
        nhits = 0
        for i in self.alignmentFreq:
            if i > 0:
                nhits += 1
        return round((float(nhits) / float(self.querySize)) * 100, 2)

    def get_mean_cov(self, s, e):
        return float(sum(self.alignmentFreq[s:e])) / float(len(self.alignmentFreq[s:e]))

    # def check_indels(self):
    #     has_indel = False
    #     for i, blatResult in enumerate(self.sortedResults):
    #         # nmatch, ngaps, in_target, br, perc_ident = self.blat_results[i]
    #         blatResult.set_mean_cov(self.get_mean_cov(blatResult.qstart(), blatResult.qend()))
    #         # keep_clipped = (mean_cov<4 and ((br.get_nmatch_total()<30 and not br.in_repeat) or br.get_nmatch_total()>=30))
    #         # keep_clipped = keep_clipped or (br.in_target and mean_cov<10)
    #         # print nmatch, ngaps, br.mean_cov
    #         if i == 0 and self.check_blat_indel(br):
    #             has_indel = True
    #             self.logger.info('Contig has indel, returning %r' % has_indel)
    #             return has_indel
    #         else: #if keep_clipped :
    #             self.logger.debug('Storing clipped blat result start %d, end %d' % (br.qstart(), br.qend()))
    #             self.clipped_qs.append((br.qstart(), br.qend(), br, i))
    #     self.logger.info('Contig does not have indel, return %r' % has_indel)
    #     return has_indel

    # def check_blat_indel(self, br):
    #     indel = False
    #     indel_size_thresh = int(self.meta_dict['params'].opts['indel_size'])
    #     self.logger.info('Checking if blat result contains an indel variant')
    #     nhits = 0
    #     for i in self.hit_freq:
    #         if i > 0:
    #             nhits += 1
    #     if br.spans_query() or (len(self.blat_results) == 1 and br.in_target):
    #         self.logger.info('Blat result spans query (%r) or only one blat result (%r) and blat result in target (%r)' % (br.spans_query(), (len(self.blat_results) == 1), br.in_target))
    #         indel = True
    #         keep_br = br.valid and br.mean_cov < 2 and br.in_target and (br.indel_maxevent_size[0] >= indel_size_thresh) and (not br.rep_man.breakpoint_in_rep[0] and not br.rep_man.breakpoint_in_rep[1])
    #         self.logger.debug('Keep blat result %r' % keep_br)
    #         if keep_br:
    #             brkpt_cov = [self.meta_dict['contig_vals'][1].get_counts(x, x, 'indel') for x in br.query_brkpts]
    #             low_cov = min(brkpt_cov) < self.meta_dict['params'].get_sr_thresh('indel')
    #             flank_match_thresh = True
    #             for fm in br.indel_flank_match:
    #                 fm_perc = round((float(fm) / float(br.get_size('query'))) * 100, 2)
    #                 if fm_perc < 10.0:
    #                     flank_match_thresh = False
    #                 self.logger.info('Indel result has matching flanking sequence of largest indel event of %d (%d of query)' % (fm, fm_perc))
    #             self.logger.info('Indel result has matching flanking sequence of largest indel event (10 perc of query) on both sides (%r)' % flank_match_thresh)
    #             in_ff, span_ff = filter_by_feature(br.get_brkpt_locs(), self.meta_dict['query_region'], self.meta_dict['params'].opts['keep_intron_vars'])
    #             if not in_ff and not low_cov and flank_match_thresh:
    #                 self.se = sv_event(br, self.meta_dict['query_region'], self.meta_dict['contig_vals'], self.meta_dict['sbam'])
    #                 self.logger.debug('Top hit contains whole query sequence, indel variant')
    #             else:
    #                 self.logger.debug('Indel in intron (%r) or low coverage at breakpoints (%r) or minimum segment size < 20 (%r), filtering out.' % (in_ff, low_cov, min(br.query_blocksizes)))
    #         else:
    #             self.logger.debug('Indel failed checking criteria: in annotated gene: %r, mean query coverage < 2: %r, in target: %r, in repeat: %r, indel size < %d: %r' % (br.valid, br.mean_cov, br.in_target, ",".join([str(x) for x in br.rep_man.breakpoint_in_rep]), indel_size_thresh, br.indel_maxevent_size[0] < indel_size_thresh))
    #     return indel

    # def get_indel_result(self):
    #     if self.se:
    #         return self.se.get_indel_result()
    #     else:
    #         return None

    # def get_svs_result(self):
    #     if self.se:
    #         return self.se.get_svs_result(self.meta_dict['query_region'], self.meta_dict['params'], self.meta_dict['disc_reads'])
    #     else:
    #         return None

    # def check_svs(self):
    #     self.logger.info('Checking for SVs')
    #     gaps = [(0, self.qsize)]
    #     if len(self.clipped_qs) > 1:
    #         self.logger.debug('Iterating through %d clipped blat results.' % len(self.clipped_qs))
    #         merged_clip = [0, None]
    #         for i in range(len(self.clipped_qs)):
    #             qs, qe, br, idx = self.clipped_qs[i]
    #             self.logger.debug('Blat result with start %d, end %d, chrom %s' % (qs, qe, br.get_name('hit')))
    #             gaps = self.iter_gaps(gaps, self.clipped_qs[i], i)
    #             if self.se.qlen > merged_clip[0]: # and self.se.in_target :
    #                 merged_clip = [self.se.qlen, self.se]
    #         self.se = merged_clip[1]
    #     else:
    #         self.logger.info('There are no more than 1 clipped blat results, not continuing with SVs calling.')

    #     if self.se_valid():
    #         return True
    #     else:
    #         return False

    # def se_valid(self):
    #     valid = False
    #     if self.se and len(self.se.blat_res) > 1 and self.se.in_target:
    #         nmissing_query_cov = len(filter(lambda y: y, map(lambda x: x == 0, self.se.query_cov)))
    #         if nmissing_query_cov < self.meta_dict['params'].get_min_segment_length('trl'):
    #             valid = True
    #     return valid

    # def check_add_br(self, qs, qe, gs, ge, br) :
    #     self.logger.info('Checking to add blat result with start %d, end %d'%(qs, qe))
    #     add = False
    #     over_perc = round((float(min(qe,ge)-max(qs,gs)) / float(qe-qs)) * 100) # Calc % of segment overlaps with gap
    #     ov_right = 0 # Check overlap with other aligned segments
    #     if qe > ge : ov_right = abs(qe-ge)
    #     ov_left = 0
    #     if qs < gs : ov_left = abs(qs-gs)
    #     br.set_segment_overlap(ov_left, ov_right)
    #     max_seg_overlap = max(ov_right,ov_left)
    #     self.logger.debug('Blat query segment overlaps gap by %f'%over_perc)
    #     self.logger.debug('Max segment overlap %f'%max_seg_overlap)
    #     self.logger.debug('Event in target %r and blat result in target %r'%(self.se.in_target, br.in_target))
    #     if over_perc >= 50 and (max_seg_overlap < 15 or (br.in_target and self.se.in_target) ) : # and (self.se.in_target or br.in_target) : 
    #         add = True
    #     self.logger.debug('Add blat result to SV event %r'%add)
    #     return add

    # def iter_gaps(self, gaps, cq, iter) :
    #     new_gaps = []
    #     qs, qe, br, idx = cq
    #     hit = False
    #     for gap in gaps :
    #         gs, ge = gap
    #         self.logger.debug('Gap coords %d, %d'%(gs, ge))
    #         if (qs >= gs and qs <= ge) or (qe <= ge and qe >= gs) :
    #             ngap = []
    #             if qs > gs : 
    #                 if (qs-1-gs) > 10 : 
    #                     ngap.append((gs,qs-1))
    #             if qe < ge : 
    #                 if (ge-qe+1) > 10 :
    #                     ngap.append((qe+1,ge))
    #             if iter == 0 : 
    #                 self.logger.debug('Creating SV event from blat result with start %d, end %d'%(qs, qe))
    #                 self.se = sv_event(br, self.meta_dict['query_region'], self.meta_dict['contig_vals'], self.meta_dict['sbam'])
    #                 new_gaps.extend(ngap)
    #                 hit = True
    #             elif self.check_add_br(qs, qe, gs, ge, br) :
    #                 self.logger.debug('Adding blat result to event')
    #                 new_gaps.extend(ngap)
    #                 self.se.add(br)
    #                 hit = True
    #             else :
    #                 new_gaps.append(gap)
    #         else :
    #             new_gaps.append(gap)
    #         self.logger.debug('New gap coords %s'%(",".join([str(x) for x in new_gaps])))
    #     if not hit :
    #         self.se.check_previous_add(br)
    #     return new_gaps
