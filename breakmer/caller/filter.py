#! /usr/bin/local/python
# -*- coding: utf-8 -*-

import breakmer.utils as utils

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


class Filter:
    """
    """

    def __init(self):
        self.name = None
        self.svType = None
        self.breakpoints = None
        self.description = None

    def set_values(self, values):
        name, svType, breakpoints, description = values
        self.name = name.lower()
        self.svType = svType
        self.breakpoints = self.parse_breakpoints(breakpoints)
        self.description = description

    def parse_breakpoints(self, breakpoints):
        """ """
        breakpointList = breakpoints.split(',')
        brkpts = []
        for breakpoint in breakpointList:
            chrom, bps = breakpoint.split(':')
            bpSplit = bps.split('-')
            brkptList = [chrom]
            for bp in bpSplit:
                brkptList.append(bp)
            brkpts.append(brkptList)
        return brkpts


class ResultFilter:
    """
    """

    def __init__(self, filterFn, params):
        self.loggingName = 'breakmer.caller.filter'
        self.filterFn = filterFn
        self.filters = []
        self.params = params
        self.setup()

    def setup(self):
        if self.filterFn:
            for line in open(self.filterFn, 'rU'):
                line = line.strip()
                resultFilter = Filter()
                resultFilter.set_values(line.split('\t'))
                self.filters.append(resultFilter)

# Filters for events
# 1. Contig complexity
        # avg_comp, comp_vec = calc_contig_complexity(self.contig_seq)
                # brkpt_rep_filt = False
                    # brkpt_rep_filt = brkpt_rep_filt or (comp_vec[qb[0]] < (avg_comp / 2))
                # brkpt_rep_filt = brkpt_rep_filt or (len(filter(lambda x: x, brkpts['f'])) > 0)

    def check_filters(self, svEvent):
        """
        """
        if len(self.filters) > 0:
            # Check if event is in the pre-defined filters
            self.check_defined_filters(svEvent)

        print 'check filters', svEvent.resultValues.svType, svEvent.resultValues.svSubtype
        if svEvent.svType == 'indel':
            self.filter_indel(svEvent)
        elif svEvent.resultValues.svType == 'rearrangement':
            if svEvent.resultValues.svSubtype == 'trl':
                print 'Checking trl filters'
                self.filter_trl(svEvent)
            else:
                self.filter_rearr(svEvent)

    def check_defined_filters(self, svEvent):
        """The event must match
            1. target name
            2. SV type (indel, trl, rearrangement_inversion, rearrangement_tandem_dup)
            3. breakpoints
            4. description
        """
        for SVFilter in self.filters:
            nameMatch = svEvent.contig.get_target_name().lower() == SVFilter.name
            typeMatch = svEvent.svType == SVFilter.svType
            eventBrkpts = svEvent.get_genomic_brkpts()
            ebps = []
            for key in eventBrkpts:
                ebps.extend(eventBrkpts[key])
            bpMatches = True
            for eventBrkpt in ebps:
                # Should be a tuple with chr, bp1, bp2 or chr, bp1
                match = False
                for filterBrkpt in SVFilter.breakpoints:
                    if len(eventBrkpt) == len(filterBrkpt):
                        bpMatch = True
                        # Check chroms
                        bpMatch = eventBrkpt[0].replace('chr', '') == filterBrkpt[0].replace('chr', '')
                        if bpMatch:
                            for v1, v2 in zip(eventBrkpt[1:], filterBrkpt[1:]):
                                if int(v1) != int(v2):
                                    bpMatch = False
                                    break
                            if bpMatch:
                                match = True
                                break
                boMatches = bpMatches and match
            if nameMatch and typeMatch and bpMatches:
                svEvent.set_filtered('Matched input filter variant')

    def filter_indel(self, svEvent):
        """ """
        indelSizeThresh = int(self.params.get_param('indel_size'))
        utils.log(self.loggingName, 'info', 'Checking if blat result contains an indel variant')
        blatResult = svEvent.blatResults[0][1]
        keep_br = blatResult.valid and blatResult.meanCov < 2 and blatResult.in_target and (blatResult.indel_maxevent_size[0] >= indelSizeThresh)
        utils.log(self.loggingName, 'debug', 'Keep blat result %r' % keep_br)

        # Determine the uniqueness of the realignment.
        svFilterValues = svEvent.resultValues.filterValues
        uniqRealignment = svFilterValues.resultMeanHitFreq < 2
        indelSize = svFilterValues.maxEventSize >= indelSizeThresh
        brkptCoverages = svFilterValues.brkptCoverages[0] >= self.params.get_sr_thresh('indel')
        minFlankMatches = min(svFilterValues.flankMatchPercents) >= 10.0

        if uniqRealignment and indelSize and brkptCoverages and minFlankMatches:
            utils.log(self.loggingName, 'debug', 'Indel meets basic filtering requirements.')
        else:
            utils.log(self.loggingName, 'debug', 'Indel filtered due to non-unique realignment (%r), less than input size threshold (%r), low coverage at breakpoints (%r), or contig edge realignment not long enough (%r), filter status set to True.' % (uniqRealignment, indelSize, brkptCoverages, minFlankMatches))
            filterReasons = []
            if not uniqRealignment:
                filterReasons.append('Non-unique realignment (%d) > 2' % svFilterValues.resultMeanHitFreq)
            if not indelSize:
                filterReasons.append('Max indel size (%d) is less than %d' % (svFilterValues.maxEventSize, indelSizeThresh))
            if not brkptCoverages:
                filterReasons.append('Minimum coverage at breakpoints (%d) less than input threshold %d' % (svFilterValues.brkptCoverages[0], self.params.get_sr_thresh('indel')))
            if not minFlankMatches:
                filterReasons.append('Minimum percentage of contig sequence that realigns to the reference to the left or right of the indel event less than 10.0 percent (%d)' % min(svFilterValues.flankMatchPercents))
            svEvent.set_filtered(','.join(filterReasons))

    def filter_rearr(self, svEvent):
        # in_ff, span_ff = filter_by_feature(brkpts, query_region, params.opts['keep_intron_vars'])
        # filter = (min(brkpt_counts['n']) < params.get_sr_thresh('rearrangement')) or self.blatResultsSorted[0][1] < params.get_min_segment_length('rearr') or (in_ff and span_ff) or (disc_read_count < 1) or (rearr_type == 'rearrangement') or (min(brkpt_kmers) == 0)
        svFilterValues = svEvent.resultValues.filterValues
        missingQueryCoverage = svFilterValues.missingQueryCoverage >= self.params.get_min_segment_length('rearr')
        brkptCoverages = svFilterValues.brkptCoverages[0] >= self.params.get_sr_thresh('rearrangement')
        minSegmentLen = svFilterValues.minSegmentLen >= self.params.get_min_segment_length('rearr')
        discReadCount = svEvent.discReadCount >= 1
        minBrkptKmers = svEvent.minBrkptKmers > 0

        if brkptCoverages and minSegmentLen and discReadCount and minBrkptKmers:
            utils.log(self.loggingName, 'info', 'Rearrangement meets basic filtering requirements.')
        else:
            filteredReasons = []
            if not missingQueryCoverage:
                logMsg = 'No realignment for %d bases in the contig sequence, more than threshold %d.' % (svFilterValues.missingQueryCoverage, self.params.get_min_segment_length('rearr'))
                utils.log(self.loggingName, 'info', logMsg)
                filteredReasons.append(logMsg)
            if not brkptCoverages:
                logMsg = 'Minimum coverage at breakpoints (%d) less than input threshold %d.' % (svFilterValues.brkptCoverages[0], self.params.get_sr_thresh('rearrangement'))
                utils.log(self.loggingName, 'info', logMsg)
                filteredReasons.append(logMsg)
            if not minSegmentLen:
                logMsg = 'The minimum realigned segment length (%d) is less than the input threshold %d.' % (svFilterValues.minSegmentLen, self.params.get_min_segment_length('rearr'))
                utils.log(self.loggingName, 'info', logMsg)
                filteredReasons.append(logMsg)
            if not discReadCount:
                logMsg = 'The number of discordant read pairs supporting the event (%d) is less than 1.' % svEvent.discReadCount
                utils.log(self.loggingName, 'info', logMsg)
                filteredReasons.append(logMsg)
            if not minBrkptKmers:
                logMsg = 'There were no variant kmers at the one or more of the breakpoint locations.'
                utils.log(self.loggingName, 'info', logMsg)
                filteredReasons.append(logMsg)
            svEvent.set_filtered(','.join(filteredReasons))

    def filter_trl(self, svEvent):
        svFilterValues = svEvent.resultValues.filterValues
        maxBrkptCoverages = svFilterValues.brkptCoverages[1] >= self.params.get_sr_thresh('trl')
        minBrkptCoverages = svFilterValues.brkptCoverages[0] >= self.params.get_sr_thresh('trl')
        discReadCount = svEvent.resultValues.discReadCount
        minSegmentLen = svFilterValues.minSegmentLen >= self.params.get_min_segment_length('trl')
        minBrkptKmers = svFilterValues.minBrkptKmers > 0
        minSeqComplexity = svFilterValues.seqComplexity >= 25.0
        startEndMissingQueryCoverage = svFilterValues.startEndMissingQueryCoverage <= 5.0
        maxSegmentOverlap = svFilterValues.maxSegmentOverlap < 5
        maxMeanHitFreq = svFilterValues.maxMeanCoverage < 10
        nReadStrands = svFilterValues.nReadStrands > 1
        maxRealignmentGaps = svFilterValues.maxRealignmentGaps

        strictFilter = [minSeqComplexity, startEndMissingQueryCoverage, minSegmentLen, maxRealignmentGaps, maxMeanHitFreq, nReadStrands]
        nStrictFiltersFail = 0
        for f in strictFilter:
            if not f:
                nStrictFiltersFail += 1

        print 'max breakpoint coverages', maxBrkptCoverages
        if not maxBrkptCoverages:
            logMsg = 'Maximum breakpoint coverages (%d) did not meet input threshold %d.' % (svFilterValues.brkptCoverages[1], self.params.get_sr_thresh('trl'))
            utils.log(self.loggingName, 'info', logMsg)
            svEvent.set_filtered(logMsg)
        else:
            if discReadCount >= 2 and minSegmentLen and minBrkptCoverages and minBrkptKmers:
                utils.log(self.loggingName, 'info', 'Translocation event passed all basic requirements. Filter set 1.')
            elif discReadCount == 0 and nStrictFiltersFail <= 1:
                utils.log(self.loggingName, 'info', 'Translocation event passed all basic requirements. Filter set 2.')
            else:
                logMsg = 'Translocation failed to pass filters.'
                filteredReasons = [logMsg]
                utils.log(self.loggingName, 'info', logMsg)
                svEvent.set_filtered(logMsg)
