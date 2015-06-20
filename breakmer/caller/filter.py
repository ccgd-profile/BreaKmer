#! /usr/bin/local/python
# -*- coding: utf-8 -*-

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


class Filter:
    """
    """

    def __init(self, name, type, breakpoints, description):
        self.name = name
        self.type = type
        self.breakpoints = breakpoints
        self.description = description


class ResultFilter:
    """
    """

    def __init__(self, filterFn):
        self.filterFn = filterFn
        self.filters = []
        self.setup()

    def setup(self):
        if self.filterFn:
            for line in open(self.filterFn, 'rU'):
                line = line.strip()
                resultFilter = Filter(line.split('\t'))
                self.filters.append(resultFilter)

# Filters for events
# 1. Check this for non-indels:
            # nMissingQueryCoverage = len(filter(lambda y: y, map(lambda x: x == 0, self.queryCoverage)))
            # if nMissingQueryCoverage < self.meta_dict['params'].get_min_segment_length('trl'):
            #     valid = True
# 2. Contig complexity 
        # avg_comp, comp_vec = calc_contig_complexity(self.contig_seq)
                # brkpt_rep_filt = False
                    # brkpt_rep_filt = brkpt_rep_filt or (comp_vec[qb[0]] < (avg_comp / 2))
                # brkpt_rep_filt = brkpt_rep_filt or (len(filter(lambda x: x, brkpts['f'])) > 0)


    def check_filters(self, svEvent):
        """
        """

    def filter_rearr(self, query_region, params, brkpts, brkpt_counts, brkpt_kmers, rearr_type, disc_read_count):
        # in_ff, span_ff = filter_by_feature(brkpts, query_region, params.opts['keep_intron_vars'])
        # filter = (min(brkpt_counts['n']) < params.get_sr_thresh('rearrangement')) or self.blatResultsSorted[0][1] < params.get_min_segment_length('rearr') or (in_ff and span_ff) or (disc_read_count < 1) or (rearr_type == 'rearrangement') or (min(brkpt_kmers) == 0)
        filter = (min(brkpt_counts['n']) < params.get_sr_thresh('rearrangement')) or self.blatResultsSorted[0][1] < params.get_min_segment_length('rearr') or (disc_read_count < 1) or (rearr_type == 'rearrangement') or (min(brkpt_kmers) == 0)
        self.logger.info('Check filter for rearrangement')
        self.logger.info('Filter by feature for being in exon (%r) or spanning exon (%r)'%(in_ff, span_ff))
        self.logger.info('Split read threshold %d, breakpoint read counts %d'%(min(brkpt_counts['n']),params.get_sr_thresh('rearrangement')))
        self.logger.info('Minimum segment length observed (%d) meets threshold (%d)'%(self.blatResultsSorted[0][1], params.get_min_segment_length('rearr')))
        self.logger.info('Minimum discordant read pairs for rearrangement (%d)'%(disc_read_count))
        return filter

    def filter_trl(self, br_valid, query_region, params, brkpt_counts, brkpt_kmers, disc_read_count, anno_genes, max_repeat, rep_filt) :
        filter = br_valid[1] or (max(brkpt_counts['d']) < params.get_sr_thresh('trl'))
        self.logger.debug('Check translocation filter')
        self.logger.debug('All blat result segments are within annotated or pre-specified regions %r' % br_valid[0])
        self.logger.debug('All blat result segments are within simple repeat regions that cover > 75.0 percent of the segment %r'%br_valid[1])
        self.logger.debug('The maximum read count support around breakpoints %d meets split read threshold %d'%(max(brkpt_counts['d']),params.get_sr_thresh('trl')))
        self.logger.debug('The minimum number of kmers at breakpoints %d' % min(brkpt_kmers))
        self.logger.debug('The maximum repeat overlap by a blat result: %f' % max_repeat)
        if not filter :
            self.logger.debug('Filter %r, checking discordant read counts %d' % (filter, disc_read_count))
            if disc_read_count < 2 :
#        print 'Filter due to repeat', rep_filt
                if (self.blatResultsSorted[0][1] < params.get_min_segment_length('trl')) or (min(brkpt_counts['n']) < params.get_sr_thresh('trl')) or (min(brkpt_kmers)==0) or rep_filt :
                    self.logger.debug('Shortest segment is < %d bp with %d discordant reads. Filtering.'%(params.get_min_segment_length('trl'), disc_read_count))
                    self.logger.debug('The minimum read count support for breakpoints %d meets split read threshold %d'%(min(brkpt_counts['n']),params.get_sr_thresh('trl')))
                    self.logger.debug('The minimum number of kmers at breakpoints %d'%min(brkpt_kmers))
                    filter = True
                elif disc_read_count == 0:
                    # Check a number of metrics for shortest blat segment
                    br_qs = self.blatResultsSorted[0][0].qstart()
                    br_qe = self.blatResultsSorted[0][0].qend()
                    low_complexity = self.minseq_complexity(self.contig_seq[br_qs:br_qe],3) < 25.0 # Complexity of blat segment
                    missing_qcov = self.missing_query_coverage() > 5.0
                    short = self.blatResultsSorted[0][1] <= round(float(len(self.contig_seq))/float(4.0))
                    self.logger.debug('Checking length of shortest sequence, considered too short %r, %d, %f'%(short, self.blatResultsSorted[0][1], round(float(len(self.contig_seq))/float(4.0))) )
                    overlap = max(self.blatResultsSorted[0][0].seg_overlap) > 5
                    gaps_exist = max(self.blatResultsSorted[0][0].gaps['query'][0], self.blatResultsSorted[0][0].gaps['hit'][0]) > 0
                    low_uniqueness = self.check_uniqueness()
                    intergenic_regions = 'intergenic' in anno_genes
                    read_strand_bias = self.check_read_strands()
                    check_values = [low_complexity, missing_qcov, short, overlap, gaps_exist, low_uniqueness, read_strand_bias, intergenic_regions]
                    self.logger.debug('Discordant read count of 0 checks %s'%(",".join([str(x) for x in check_values])))
                    num_checks = 0
                    for check in check_values :
                        if check : 
                            num_checks += 1
                    if num_checks > 1 : 
                        self.logger.info('Two or more filter checks, setting filtering to true for contig')
                        filter = True
        return filter

    def check_uniqueness(self):
        low_unique = False
        for br_vals in self.blatResultsSorted:
            if not br_vals[0].in_target:
                if br_vals[0].mean_cov > 4:
                    low_unique = True
            else:
                if br_vals[0].mean_cov > 10:
                    low_unique = True
        return low_unique

    def check_read_strands(self):
        same_strand = False
        strands = []
        for read in self.contig_reads :
            strand = read.id.split("/")[1] 
            strands.append(strand)
        if len(set(strands)) == 1 :
            same_strand = True
        self.logger.debug('Checking read strands for contig reads %s'%(",".join([read.id for read in self.contig_reads])))
        self.logger.debug('Reads are on same strand: %r'%same_strand)
        return same_strand

    def minseq_complexity(self, seq, N):
        self.logger.debug('Checking sequence complexity of blat result segment %s using %d-mers' % (seq, N))
        nmers = {}
        total_possible = len(seq) - 2
        for i in range(len(seq) - (N - 1)):
            nmers[str(seq[i:i+N]).upper()] = True
        complexity = round((float(len(nmers))/float(total_possible))*100,4)
        self.logger.debug('Complexity measure %f, based on %d unique %d-mers observed out of a total of %d %d-mers possible'%(complexity,len(nmers),N,total_possible,N))
        return complexity

    def missing_query_coverage(self):
        missing_cov = 0
        for i in self.queryCoverage:
            if i == 0:
                missing_cov += 1
            else:
                break

        for i in reversed(self.queryCoverage):
            if i == 0:
                missing_cov += 1
            else:
                break

        perc_missing = round((float(missing_cov)/float(len(self.contig_seq)))*100, 4)
        self.logger.debug('Calculated %f missing coverage of blat query sequence at beginning and end' % perc_missing)
        return perc_missing

    def multiple_genes(self, chrs, brkpts, anno_genes):
        mult_genes = True
        if len(set(anno_genes)) == 1:
            self.logger.debug('One annotated gene among SVs breakpoints: %s' % (",".join(anno_genes)))
            mult_genes = False
        elif self.dup_gene_names(anno_genes) and len(set(chrs)) == 1 and ((max(brkpts) - min(brkpts)) < 10000):
            self.logger.debug('Anno genes are not the same, but similar and breakpoints are less than 10Kb from each other %s' % (",".join(anno_genes)))
            mult_genes = False
        self.logger.debug('Test whether SVs breakpoints are in multiple genes %r' % mult_genes)
        return mult_genes

    def dup_gene_names(self, anno_genes):
        found_dup = False
        for i in range(len(anno_genes)-1) :
            g1 = anno_genes[i]
            for g2 in anno_genes[(i+1):] :
                if (g1.find(g2) > -1) or (g2.find(g1) > -1) :
                    found_dup = True
        return found_dup 