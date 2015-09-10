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
            return
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
        self.mergedRecords = []  # List of tuples containing indices of realignment records that need to be merged.
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
            cond1 = self.results[0].spans_query() and (self.ngaps > 0)
            cond2 = (len(self.results) == 1) and self.get_query_coverage() >= 90.0 and (self.ngaps > 0)
            cond3 = (len(self.results) > 1) and (self.get_query_coverage() >= 90.0)
            self.targetHit = cond1 or cond2 or cond3
        utils.log(self.loggingName, 'debug', 'Checking if query is a target hit or not %r' % self.targetHit)

        if self.targetHit:
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
                        # print i, result, containSeg, withinSeg, overlapSeg
                        if containSeg or withinSeg:
                            addSegment = False
                        elif overlapSeg:
                            if (result.qstart() >= segment[0] and result.qstart() <= segment[1]):
                                overlapBp = segment[1] - result.qstart()
                                # print 'Overlapbp', overlapBp
                                if overlapBp < 20:
                                    resultOverlap += overlapBp
                                    addSegment = True
                            elif (result.qend() >= segment[0] and result.qend() <= segment[1]):
                                overlapBp = result.qend() - segment[0]
                                # print 'Overlapbp', overlapBp
                                if overlapBp < 20:
                                    resultOverlap += overlapBp
                                    addSegment = True
                    if addSegment and (result.get_query_span() - resultOverlap) > 20:
                        segments.append((result.qstart(), result.qend(), result))
                self.targetSegmentsSorted = sorted(segments, key=lambda x: x[0])
                for i in range(1, len(self.targetSegmentsSorted)):
                    lResult = self.targetSegmentsSorted[i - 1][2]
                    rResult = self.targetSegmentsSorted[i][2]
                    qgap = rResult.qstart() - lResult.qend()
                    tgap = rResult.tstart() - lResult.tend()
                    if (tgap < 0 and (abs(tgap) > abs(qgap))) or (lResult.strand != rResult.strand):
                        # Tandem dup or inversion
                        break
                    else:
                        self.mergedRecords.append((i - 1, i))
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
            blatResult.set_realign_freq(self.get_align_freq(blatResult.qstart(), blatResult.qend()))

        if len(self.results) == 0:
            self.hasResults = False

    def merge_records(self):
        """ """
        mergedResults = []
        mapResults = {}
        for mergeIdx in self.mergedRecords:
            lResult = self.targetSegmentsSorted[mergeIdx[0]][2].resultValues
            rResult = self.targetSegmentsSorted[mergeIdx[1]][2].resultValues
            # print 'left', lResult
            # print 'right', rResult
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
                # print 'New result', i, mergedResults[mapResults[i]]
                if mapResults[i] not in outputIdx:
                    outputIdx.append(mapResults[i])
                    lOuts.append(self.format_to_blat_output(mergedResults[mapResults[i]]))
            else:
                # print 'Result', self.targetSegmentsSorted[i][2].resultValues
                lOuts.append(self.format_to_blat_output(self.targetSegmentsSorted[i][2].resultValues))
        # print lOuts
        return lOuts

    def format_to_blat_output(self, resultValues):
        """ """
        outStr = []
        keys = ['matches',
                'mismatches',
                'repmatches',
                'ncount',
                'qNumInsert',
                'qBaseInsert',
                'tNumInsert',
                'tBaseInsert',
                'strand',
                'qName',
                'qSize',
                'qStart',
                'qEnd',
                'tName',
                'tSize',
                'tStart',
                'tEnd',
                'blockCount',
                'blockSizes',
                'qStarts',
                'tStarts']
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
        return lResult

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

    def get_align_freq(self, s, e):
        return float(sum(self.alignmentFreq[s:e])) / float(len(self.alignmentFreq[s:e]))
