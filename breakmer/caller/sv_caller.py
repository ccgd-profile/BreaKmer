#! /usr/bin/local/python
# -*- coding: utf-8 -*-

import sys
import os
import math

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


class sv_event :
    def __init__(self, blat_res, query_region, contig_vals, sample_bam) :
        self.events = []
        self.blat_res = []
        self.br_sorted = []
        self.sample_bam = sample_bam
        self.qlen = 0
        self.nmatch = 0
        self.in_target = False
        self.query_region = query_region
        self.contig_seq = contig_vals[0]
        self.contig_rcounts = contig_vals[1]
        self.contig_id = contig_vals[2]
        self.contig_reads = contig_vals[3]
        self.nkmers = contig_vals[4]
        self.contig_kmer_locs = contig_vals[5]
        self.logger = logging.getLogger('root')
        self.valid = True
        self.in_rep = False
        self.query_size = None
        self.query_cov = [0]*len(self.contig_seq)
        self.homology = {'non':[False,None], 'in':[False,None]}
        self.result_values = {'anno_genes':None, 'target_breakpoints':None, 'align_cigar':'', 'mismatches':0, 'strands':None, 'repeat_matching':None, 'sv_type':None, 'split_read_count':None, 'nkmers':self.nkmers, 'disc_read_count':0, 'contig_id':query_region[3] + "_" + self.contig_id,'contig_seq':self.contig_seq, 'sv_subtype':None, 'breakpoint_coverages':0}
        self.add(blat_res)
    #*********************************************************

    #*********************************************************
    def add(self, br) :
        self.blat_res.append((br.get_coords('query')[0],br))
        for i in range(br.get_coords('query')[0], br.get_coords('query')[1]) : 
            #print 'incrementing query cov', i, self.query_cov, len(self.query_cov), br.get_coords('query')
            self.query_cov[i] += 1
        if not self.query_size :
            self.query_size = br.get_size('query')
        self.qlen += br.get_query_span()
        self.nmatch += br.get_nmatch_total()
        self.in_target = self.in_target or br.in_target
        self.in_rep = self.in_rep and (br.repeat_overlap > 75.0)
        self.valid = self.valid and br.valid
        self.br_sorted.append((br,br.get_nmatch_total()))
    #*********************************************************

    def check_previous_add(self, br) :
#    print 'Checking previous br add'
        ncoords = br.get_coords('query')
#    print 'br new coords', ncoords
        prev_br, prev_nmatch = self.br_sorted[-1]
        prev_coords = prev_br.get_coords('query')
#    print 'Previous br coords', prev_coords
        if ncoords[0] == prev_coords[0] and ncoords[1] == prev_coords[1] :
#      print 'Coords match'
            n_nmatch = br.get_nmatch_total()
            if abs(prev_nmatch-n_nmatch) < 10 : 
#        print 'score different is less than 10'
#        print 'prev br', prev_br, prev_br.in_target
#        print 'new br', br, br.in_target
                if not prev_br.in_target and br.in_target :
                    self.br_sorted[-1] = (br, n_nmatch)
                    self.blat_res[-1] = (ncoords[0], br)
                    self.in_target = True
#          self.homology['non'] = [True, br, prev_br.get_name('hit'), len(self.br_sorted)-1]
#        elif not br.in_target and prev_br.in_target :
#          self.homology['in'] = [True, br.get_name('hit')]
#    print 'homology check', self.homology
#    if self.homology['non'][0] and self.homology['in'][0] :
#      if self.homology['non'][2] == self.homology['in'][1] :
#          print 'Swap brs', br, self.br_sorted[-1], self.blat_res[-1] #self.homology['non'][3], self.homology['non'][1]
#        self.br_sorted[self.homology['non'][3]] = self.homology['non'][1]  

    #*********************************************************
    def set_result_value(self, key, value) :
        self.result_values[key] = value
    #*********************************************************

    #*********************************************************
    def format_result(self, values) :
        res_lst = []
        if values : 
            for v in values :
                if not isinstance(values[v], list) : values[v] = [values[v]]
                self.result_values[v] = ",".join([str(x) for x in values[v]])
        if self.result_values['sv_subtype'] : self.result_values['sv_type'] += '_' + self.result_values['sv_subtype']
        return self.get_values()
    #*********************************************************

    #*********************************************************
    def get_brkpt_coverages(self) :
        brkpts = []
        tbp = self.result_values['target_breakpoints']
        if self.result_values['target_breakpoints'].find("(") > -1 :
            tbp = self.result_values['target_breakpoints'].split()[0]
        tbp = tbp.split(',')
        for bp in tbp :
            chrom,locs = bp.split(':')
            chrom = chrom.replace('chr','')
            ll = locs.split('-')
            if len(ll) > 1 : 
                brkpts.append((chrom,int(ll[0]),int(ll[0])+1))
                brkpts.append((chrom,int(ll[1]),int(ll[1])+1))
            else :
                brkpts.append((chrom, int(ll[0]), int(ll[0])+1))
            
        bamfile = Samfile(self.sample_bam,'rb')

        covs = [0]*len(brkpts)
        bp_index = 0
        for bp in brkpts :
            cov = 0
            c,s,e = bp
            areads = bamfile.fetch(str(c), s, e)

            for aread in areads :
                if aread.is_duplicate or aread.is_qcfail or aread.is_unmapped or aread.mapq < 10 : 
                    continue
                cov += 1
            covs[bp_index] = cov
            bp_index += 1
        return ",".join([str(x) for x in covs])
    #*********************************************************

    #*********************************************************
    def get_values(self) :
        lst = ['anno_genes', 'target_breakpoints', 'align_cigar', 'mismatches', 'strands', 'repeat_matching', 'sv_type', 'split_read_count', 'nkmers', 'disc_read_count', 'breakpoint_coverages', 'contig_id', 'contig_seq']
        out_lst = []
        for l in lst :
            if self.result_values[l] == 'trl' :
                self.result_values[l] = 'rearrangement'  
            out_lst.append(str(self.result_values[l]))
            if l == 'target_breakpoints' :
                self.result_values['breakpoint_coverages'] = self.get_brkpt_coverages()
        return out_lst
    #*********************************************************

    #*********************************************************
    def get_indel_result(self) :
        br = self.blat_res[0][1]
        self.set_result_value('anno_genes', br.get_gene_anno()) 
        self.set_result_value('repeat_matching', '0.0:'+str(br.get_nmatch_total()))
        self.set_result_value('mismatches', br.get_nmatches('mis'))
        self.set_result_value('strands', br.strand)
        self.set_result_value('target_breakpoints', br.get_brkpt_str(True))
        self.set_result_value('align_cigar', br.cigar)
        self.set_result_value('sv_type', 'indel')
        self.set_result_value('split_read_count',",".join([str(self.contig_rcounts.get_counts(x,x,'indel')) for x in br.query_brkpts]))
        return self.format_result(None)
    #*********************************************************

    #*********************************************************
    def get_brkpt_info(self, br, brkpt_d, i, last_iter) :
        ts, te = br.get_coords('hit')
        qs, qe = br.get_coords('query')
        target_key = 'in_target' if br.in_target else 'other'
        brkpt_d['chrs'].append(br.get_name('hit'))
        brkpt_d['tcoords'].append((ts,te))
        tbrkpt = []
        filt_rep_start = None
        if i == 0 :
            brkpt_d['q'][0] = [max(0,qs-1),qe]
            brkpt_d['q'][1].append([qe, qe-brkpt_d['q'][0][0], None])
            tbrkpt = [te]
            filt_rep_start = br.filter_reps_edges[0]
            if br.strand == '-' :
                tbrkpt = [ts]
                filt_rep_start = br.filter_reps_edges[0]
        elif last_iter : 
            brkpt_d['q'][1][-1][2] = qe - brkpt_d['q'][1][-1][0]
            brkpt_d['q'][1].append([qs,qs-brkpt_d['q'][0][0],qe-qs])
            tbrkpt = [ts]
            filt_rep_start = br.filter_reps_edges[0]
            if br.strand == '-' :
                tbrkpt = [te] 
                filt_rep_start = br.filter_reps_edges[1]
        else :
            brkpt_d['q'][1][-1][2] = qe - brkpt_d['q'][1][-1][1]
            brkpt_d['q'][1].append([qs,qs-brkpt_d['q'][0][0],qe-qs])
            brkpt_d['q'][1].append([qe, qe-qs, None])
            brkpt_d['q'][0] = [qs, qe]
            tbrkpt = [ts, te]
            if br.strand == "-" :
                filt_rep_start = br.filter_reps_edges[1]
                tbrkpt = [te, ts]

        brkpt_d['brkpt_str'].append('chr'+str(br.get_name('hit')) + ":" + "-".join([str(x) for x in tbrkpt]))
        brkpt_d['r'].extend(tbrkpt)
        brkpt_d['f'].append(filt_rep_start)
        brkpt_d['t'][target_key] = (br.get_name('hit'),tbrkpt[0])
        brkpt_d['formatted'].append( 'chr'+str(br.get_name('hit')) + ":" + "-".join([str(x) for x in tbrkpt]))
        return brkpt_d
    #*********************************************************

    #*********************************************************
    def get_svs_result(self, query_region, params, disc_reads) :
        self.logger.info('Resolving SVs call from blat results')
        blat_res = self.blat_res
        blat_res_sorted = sorted(blat_res, key=lambda blat_res: blat_res[0]) 
        brkpts = {'t':{'in_target':None, 'other':None }, 'formatted':[], 'r':[], 'q': [[0,0],[]], 'chrs':[], 'brkpt_str':[], 'tcoords':[], 'f':[]}
        res_values = {'target_breakpoints':[], 'align_cigar':[], 'sv_type':'', 'strands':[], 'mismatches':[], 'repeat_matching':[], 'anno_genes':[], 'disc_read_count':0 }
        br_valid = [True, True]
        max_repeat = 0.0

        for i in range(len(blat_res_sorted)) :
            br = blat_res_sorted[i][1]
            br_valid[0] = br_valid[0] and br.valid
            br_valid[1] = br_valid[1] and (br.rep_man.simple_rep_overlap > 75.0)
            max_repeat = max(max_repeat, br.repeat_overlap)
            res_values['repeat_matching'].append(":".join([str(br.repeat_overlap), str(br.get_nmatch_total()), str(round(br.mean_cov,3))]))
            res_values['anno_genes'].append(br.get_gene_anno())
            res_values['align_cigar'].append(br.cigar)
            res_values['strands'].append(br.strand)
            res_values['mismatches'].append(br.get_nmatches('mis'))
            brkpts = self.get_brkpt_info(br, brkpts, i, i==(len(blat_res_sorted)-1))

        result = None
        self.br_sorted = sorted(self.br_sorted, key=lambda br: br[1])
        if not self.multiple_genes(brkpts['chrs'], brkpts['r'], res_values['anno_genes']) : 
            brkpt_counts, brkpt_kmers, brkpt_rep_filt = self.get_brkpt_counts_filt(brkpts, 'rearr') 
            rearr_type, disc_read_support = self.define_rearr(brkpts['r'], res_values['strands'], brkpts['tcoords'], disc_reads) 
            if not self.filter_rearr(query_region, params, brkpts['r'], brkpt_counts, brkpt_kmers, rearr_type, disc_read_support) :
                res_values['sv_type'] = 'rearrangement'
                if rearr_type != 'rearrangement' : res_values['sv_subtype'] = rearr_type
                res_values['disc_read_count'] = disc_read_support
                res_values['anno_genes'] = list(set(res_values['anno_genes']))
                res_values['target_breakpoints'] = brkpts['brkpt_str']
                res_values['split_read_count'] = brkpt_counts['b']
                if 'rearrangement' in params.opts['var_filter'] : 
                    result = self.format_result(res_values)
        elif max(self.contig_rcounts.others) >= params.get_sr_thresh('trl') :
            brkpt_counts, brkpt_kmers, brkpt_rep_filt = self.get_brkpt_counts_filt(brkpts, 'trl') 
            disc_read_count = self.check_disc_reads(brkpts['t'], query_region, disc_reads['disc'])
            if not self.filter_trl(br_valid, query_region, params, brkpt_counts, brkpt_kmers, disc_read_count, res_values['anno_genes'], max_repeat, brkpt_rep_filt) : 
                res_values['disc_read_count'] = disc_read_count
                res_values['sv_type'] = ['trl']
                res_values['target_breakpoints'] = brkpts['brkpt_str']
                res_values['split_read_count'] = brkpt_counts['b']
                if 'trl' in params.opts['var_filter'] : 
                    result = self.format_result(res_values)
        return result
    #*********************************************************

    #*********************************************************
    def get_brkpt_counts_filt(self, brkpts, sv_type) :
#    print 'Contig seq', self.contig_seq
#    print 'Breakpoint simple repeat filter', brkpts['f'] 
        avg_comp, comp_vec = calc_contig_complexity(self.contig_seq) 
#    print 'Contig avg complexity', avg_comp
#    print 'Contig complexity vec', comp_vec
        brkpt_rep_filt = False
        brkpt_counts = {'n':[],'d':[],'b':[]}
        brkpt_kmers = []
        for qb in brkpts['q'][1] :
            left_idx = qb[0] - min(qb[1],5)
            right_idx = qb[0] + min(qb[2],5)
#      print self.contig_rcounts.others
#      print self.contig_rcounts.indel_only
#      print qb[0], left_idx, right_idx
            bc = self.contig_rcounts.get_counts(left_idx, right_idx, sv_type)
            brkpt_counts['n'].append(min(bc))
            brkpt_counts['d'].append(min(self.contig_rcounts.get_counts((qb[0]-1), (qb[0]+1), sv_type)))
#      print 'Others counts', self.contig_rcounts.others, qb[0]
#      print 'Indel only counts', self.contig_rcounts.indel_only, qb[0]
            brkpt_counts['b'].append(self.contig_rcounts.get_counts(qb[0], qb[0], sv_type))
            brkpt_kmers.append(self.contig_kmer_locs[qb[0]])
#      print 'Breakpoint in contig', qb[0]
            brkpt_rep_filt = brkpt_rep_filt or (comp_vec[qb[0]] < (avg_comp/2))
#      print 'Breakpoint rep filter', brkpt_rep_filt, comp_vec[qb[0]]
            self.logger.debug('Read count around breakpoint %d : %s'%(qb[0],",".join([str(x) for x in bc])))
        self.logger.debug('Kmer count around breakpoints %s'%(",".join([str(x) for x in brkpt_kmers])))
        brkpt_rep_filt = brkpt_rep_filt or (len(filter(lambda x: x, brkpts['f'])) > 0)
        return brkpt_counts, brkpt_kmers, brkpt_rep_filt
    #*********************************************************

    #*********************************************************
    def call_rearr(self) :
            rearr_type, disc_read_support = self.define_rearr(brkpts['r'], res_values['strands'], brkpts['tcoords'], disc_reads) 
            if not self.filter_rearr(query_region, params, brkpts['r'], brkpt_counts, brkpt_kmers, rearr_type, disc_read_support) :
                res_values['sv_type'] = 'rearrangement'
                if rearr_type != 'rearrangement' : res_values['sv_subtype'] = rearr_type
                res_values['disc_read_count'] = disc_read_support
                res_values['anno_genes'] = list(set(res_values['anno_genes']))
                res_values['target_breakpoints'] = brkpts['brkpt_str']
                res_values['split_read_count'] = brkpt_counts['b']
                if 'rearrangement' in params.opts['var_filter'] : 
                    result = self.format_result(res_values)

    #*********************************************************

    #*********************************************************
    def call_trl(self) :
            disc_read_count = self.check_disc_reads(brkpts['t'], query_region, disc_reads['disc'])
            if not self.filter_trl(br_valid, query_region, params, brkpt_counts, brkpt_kmers, disc_read_count, res_values['anno_genes'], max_repeat, brkpt_rep_filt) : 
                res_values['disc_read_count'] = disc_read_count
                res_values['sv_type'] = ['trl']
                res_values['target_breakpoints'] = brkpts['brkpt_str']
                res_values['split_read_count'] = brkpt_counts['b']
                if 'trl' in params.opts['var_filter'] : 
                    result = self.format_result(res_values)

    #*********************************************************


    #*********************************************************
    def check_overlap(self, coord1, coord2 ) :
        contained = False
        if coord1[0] >= coord2[0] and coord1[1] <= coord2[1] :
            contained = True
        elif coord2[0] >= coord1[0] and coord2[1] <= coord1[1] :
            contained = True
        return contained 
    #*********************************************************
 
    #*********************************************************
    def define_rearr(self, brkpts, strands, tcoords, disc_reads) :
        type = 'rearrangement'
        rs = 0
        hit = False
        if len(strands) < 3 :
            if not self.check_overlap(tcoords[0], tcoords[1]) :
                self.logger.debug('Checking rearrangement type, strand1 %s, strand2 %s, breakpt1 %d, breakpt %d'%(strands[0], strands[1], brkpts[0], brkpts[1]))
                if (strands[0] != strands[1]) and (brkpts[0] < brkpts[1]) :
                    # Inversion
                    # Get discordantly mapped read-pairs
                    self.logger.debug('HIT INVERSION')
                    hit = True
                    type = 'inversion'
                    for read_pair in disc_reads['inv'] :
                        r1p, r2p, r1s, r2s, qname = read_pair
                        if r1s == 1 and r2s == 1 :
                            if (r1p <= brkpts[0]) and (r2p <= brkpts[1] and r2p >= brkpts[0]) :
                                rs += 1
                        else :
                            if (r1p <= brkpts[1] and r1p >= brkpts[0]) and r2p >= brkpts[1] :
                                rs += 1 
                elif (strands[0] == "+" and strands[1] == "+") and (brkpts[0] > brkpts[1]) :
                    self.logger.debug('HIT TANDEM DUP')
                    hit = True
                    type = 'tandem_dup'
                    # Tandem dup
                    for read_pair in disc_reads['td'] :
                        r1p, r2p, r1s, r2s, qname = read_pair 
                        if (r1p <= brkpts[0] and r1p >= brkpts[1]) and () : rs += 1
        if not hit :
            self.logger.debug('Not inversion or tandem dup, checking for odd read pairs around breakpoints')
            rs = [0] * len(brkpts)
            for i in range(len(brkpts)) : 
                b = brkpts[i]
                for read_pair in disc_reads['other'] :
                    r1p, r2p, r1s, r2s, qname = read_pair
                    if abs(r1p-b) <= 300 or abs(r2p-b) <= 300 :
                        self.logger.debug('Adding read support from read %s, with strands %d, %d and positions %d, %d for breakpoint at %d'%(qname,r1s,r2s,r1p,r2p,b))
                        rs[i] += 1
            rs = max(rs)
        return type, rs   
    #*********************************************************

    #*********************************************************
    def filter_rearr(self, query_region, params, brkpts, brkpt_counts, brkpt_kmers, rearr_type, disc_read_count ) :
        in_ff, span_ff = filter_by_feature(brkpts, query_region, params.opts['keep_intron_vars'])
        filter = (min(brkpt_counts['n']) < params.get_sr_thresh('rearrangement')) or self.br_sorted[0][1] < params.get_min_segment_length('rearr') or (in_ff and span_ff) or (disc_read_count < 1) or (rearr_type == 'rearrangement') or (min(brkpt_kmers) == 0)
        self.logger.info('Check filter for rearrangement')
        self.logger.info('Filter by feature for being in exon (%r) or spanning exon (%r)'%(in_ff, span_ff))
        self.logger.info('Split read threshold %d, breakpoint read counts %d'%(min(brkpt_counts['n']),params.get_sr_thresh('rearrangement')))
        self.logger.info('Minimum segment length observed (%d) meets threshold (%d)'%(self.br_sorted[0][1], params.get_min_segment_length('rearr')))
        self.logger.info('Minimum discordant read pairs for rearrangement (%d)'%(disc_read_count))
        return filter
    #*********************************************************

    #*********************************************************
    def filter_trl(self, br_valid, query_region, params, brkpt_counts, brkpt_kmers, disc_read_count, anno_genes, max_repeat, rep_filt) :
        filter = br_valid[1] or (max(brkpt_counts['d']) < params.get_sr_thresh('trl')) #or not br_valid[0]
        self.logger.debug('Check translocation filter')
        self.logger.debug('All blat result segments are within annotated or pre-specified regions %r'%br_valid[0])
        self.logger.debug('All blat result segments are within simple repeat regions that cover > 75.0 percent of the segment %r'%br_valid[1])
        self.logger.debug('The maximum read count support around breakpoints %d meets split read threshold %d'%(max(brkpt_counts['d']),params.get_sr_thresh('trl')))
        self.logger.debug('The minimum number of kmers at breakpoints %d'%min(brkpt_kmers))
        self.logger.debug('The maximum repeat overlap by a blat result: %f'%max_repeat)
        if not filter :
            self.logger.debug('Filter %r, checking discordant read counts %d'%(filter, disc_read_count)) 
            if disc_read_count < 2 :
#        print 'Filter due to repeat', rep_filt
                if (self.br_sorted[0][1] < params.get_min_segment_length('trl')) or (min(brkpt_counts['n']) < params.get_sr_thresh('trl')) or (min(brkpt_kmers)==0) or rep_filt :
                    self.logger.debug('Shortest segment is < %d bp with %d discordant reads. Filtering.'%(params.get_min_segment_length('trl'), disc_read_count))
                    self.logger.debug('The minimum read count support for breakpoints %d meets split read threshold %d'%(min(brkpt_counts['n']),params.get_sr_thresh('trl')))
                    self.logger.debug('The minimum number of kmers at breakpoints %d'%min(brkpt_kmers))
                    filter = True
                elif disc_read_count == 0 : 
                    # Check a number of metrics for shortest blat segment
                    br_qs = self.br_sorted[0][0].qstart()
                    br_qe = self.br_sorted[0][0].qend()
                    low_complexity = self.minseq_complexity(self.contig_seq[br_qs:br_qe],3) < 25.0 # Complexity of blat segment
                    missing_qcov = self.missing_query_coverage() > 5.0
                    short = self.br_sorted[0][1] <= round(float(len(self.contig_seq))/float(4.0))
                    self.logger.debug('Checking length of shortest sequence, considered too short %r, %d, %f'%(short, self.br_sorted[0][1], round(float(len(self.contig_seq))/float(4.0))) )
                    overlap = max(self.br_sorted[0][0].seg_overlap) > 5
                    gaps_exist = max(self.br_sorted[0][0].gaps['query'][0], self.br_sorted[0][0].gaps['hit'][0]) > 0
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
    #*********************************************************

    #*********************************************************
    def check_uniqueness(self) :
        low_unique = False
        for br_vals in self.br_sorted :
            if not br_vals[0].in_target :
                if br_vals[0].mean_cov > 4 : low_unique = True
            else :
                if br_vals[0].mean_cov > 10 : low_unique = True
        return low_unique
    #*********************************************************

    #*********************************************************
    def check_read_strands(self) :
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
    #*********************************************************

    #*********************************************************
    def minseq_complexity(self, seq, N) :
        self.logger.debug('Checking sequence complexity of blat result segment %s using %d-mers'%(seq,N))
        nmers = {}
        total_possible = len(seq) - 2
        for i in range(len(seq) - (N - 1)) :
            nmers[str(seq[i:i+N]).upper()] = True
        complexity = round((float(len(nmers))/float(total_possible))*100,4)
        self.logger.debug('Complexity measure %f, based on %d unique %d-mers observed out of a total of %d %d-mers possible'%(complexity,len(nmers),N,total_possible,N))
        return complexity
    #*********************************************************

    #*********************************************************
    def missing_query_coverage(self) :
        missing_cov = 0
        for i in self.query_cov :
            if i == 0 :
                missing_cov += 1
            else :
                break

        for i in reversed(self.query_cov) :
            if i == 0 :
                missing_cov += 1
            else : 
                break

        perc_missing = round((float(missing_cov)/float(len(self.contig_seq)))*100, 4)
        self.logger.debug('Calculated %f missing coverage of blat query sequence at beginning and end'%perc_missing)
        return perc_missing
    #*********************************************************

    #*********************************************************
    def multiple_genes(self, chrs, brkpts, anno_genes) :
        mult_genes = True
        if len(set(anno_genes)) == 1 : 
            self.logger.debug('One annotated gene among SVs breakpoints: %s'%(",".join(anno_genes)))
            mult_genes = False
        elif self.dup_gene_names(anno_genes) and len(set(chrs)) == 1 and ((max(brkpts) - min(brkpts)) < 10000 ) :
            self.logger.debug('Anno genes are not the same, but similar and breakpoints are less than 10Kb from each other %s'%(",".join(anno_genes)))
            mult_genes = False
        self.logger.debug('Test whether SVs breakpoints are in multiple genes %r'%mult_genes)
        return mult_genes
    #*********************************************************

    #*********************************************************
    def dup_gene_names(self, anno_genes) :
        found_dup = False
        for i in range(len(anno_genes)-1) :
            g1 = anno_genes[i]
            for g2 in anno_genes[(i+1):] :
                if (g1.find(g2) > -1) or (g2.find(g1) > -1) :
                    found_dup = True
        return found_dup 
    #*********************************************************

    #*********************************************************
    def check_disc_reads(self, brkpts, query_region, disc_reads ) :
        disc_read_count = 0
        if brkpts['other'][0] in disc_reads :
            for p1, p2 in disc_reads[brkpts['other'][0]] :
                d1 = abs(p1 - brkpts['in_target'][1])
                d2 = abs(p2 - brkpts['other'][1])
                if d1 <= 1000 and d2 <= 1000 : disc_read_count += 1
        return disc_read_count 
    #*********************************************************


class blat_manager :
    def __init__(self, meta_dict) :
        self.meta_dict = meta_dict
        self.fn = meta_dict['query_res_fn']
        self.qsize = 0
        self.hit_freq = []
        self.nmismatches = 0
        self.ngaps = 0
        self.has_blat_results = True
        self.blat_results = []
        self.clipped_qs = []
        self.se = None
        self.logger = logging.getLogger('root')
        self.set_values()

#     #*********************************************************
#     def target_hit(self) :
#         #single_res = (len(self.blat_results)==1) and (self.get_query_coverage() >= 90.0)
#         #mult_res1 = (len(self.blat_results)>1 and (self.get_query_coverage()>=95.0) and (self.nmismatches<5) and (self.ngaps<3))
#         #mult_res2 = self.blat_results[0][3].spans_query() and (self.blat_results[0][3].get_nmatches('mis')<5) and (self.blat_results[0][3].get_num_gaps()<3)
#         # Single hit with 90% identity.
#         indel_hit = self.blat_results[0][3].spans_query() or (len(self.blat_results) == 1 and self.get_query_coverage() >= 90.0)
# #    single_res = (len(self.blat_results)==1) and (self.get_query_coverage() >= 90.0)
# #    mult_res1 = False
# #    mult_res2 = False
#         # Or multiple hits
# #    if len(self.blat_results) > 1 :
# #      mult_res1 = self.blat_results[0][3].get_query_coverage() >= 90.0
# #      print 'blat percent ident', self.blat_results[0][3].perc_ident
# #      qcov = [0]*self.qsize
# #      for br in self.blat_results :
#                 # Calculate the total query coverage for blat results with 95% identity
# #        print 'Multiple hit', br[3].perc_ident
# #        if br[3].perc_ident >= 95.0 :
# #          qcov[br[3].qstart():br[3].qend()] = map(lambda x: x+1, qcov[br[3].qstart():br[3].qend()])
# #      nhits = 0
# #      for i in qcov :
# #        if i>0 : nhits += 1
# #      total_query_cov = round((float(nhits)/float(self.qsize))*100,2) 
# #      mult_res2 = total_query_cov >= 95.0
# #    sys.exit()
#         self.logger.debug('Checking if query is a target hit or not %r'%(indel_hit))
#         return indel_hit #single_res or mult_res1 or mult_res2
#     #*********************************************************

#     #*********************************************************
#     def set_values(self) :
#         if not self.fn :
#             self.has_blat_results = False
#         elif len(open(self.fn, 'rU').readlines()) == 0 : 
#             self.has_blat_results = False
#         else :
#             qres_f = open(self.fn, 'rU')
#             qres_lines = qres_f.readlines()
#             for line in qres_lines :
#                 line = line.strip()
#                 self.meta_dict['blat_values'] = line.split("\t")
#                 br = blat_res(self.meta_dict)
#                 score_raw = br.get_nmatch_total()
#                 ngaps = br.get_ngap_total()
#                 in_target = 1 if br.in_target else 0
#                 score_frac = float(score_raw)/float(br.get_size('query'))
#                 score = score_raw + score_frac
#                 perc_ident = br.perc_ident
#                 self.blat_results.append((score, ngaps, in_target, br, perc_ident))
#                 self.nmismatches += br.get_nmatches('mis')
#                 self.ngaps += br.get_num_gaps()
#                 if not self.qsize:
#                     self.qsize = br.get_size('query')
#                     self.hit_freq = [0]*self.qsize
#                 for i in range(br.qstart(), br.qend()) : 
#                     self.hit_freq[i] += 1
#             qres_f.close()
#         self.blat_results = sorted(self.blat_results, key=lambda blat_results: (-blat_results[0], -blat_results[4], blat_results[1]) ) 

    # def write_mod_result_file(self, fn) :
    #     blat_ff = open(fn,'w')
    #     for br in self.blat_results :
    #         blat_ff.write(br[3].get_blat_output() + "\n")
    #     blat_ff.close()

    def get_query_coverage(self) :
        nhits = 0
        for i in self.hit_freq :
            if i>0 : nhits += 1
        return round((float(nhits)/float(self.qsize))*100,2) 
    #*********************************************************

    #*********************************************************
    def get_mean_cov(self,s,e) :
        return float(sum(self.hit_freq[s:e])) / float(len(self.hit_freq[s:e]))
    #*********************************************************

    #*********************************************************
    def check_blat_indel(self, br) :
        indel = False
        indel_size_thresh = int(self.meta_dict['params'].opts['indel_size']) 
        self.logger.info('Checking if blat result contains an indel variant')
        nhits = 0
        for i in self.hit_freq :
            if i > 0 : nhits += 1  
        if br.spans_query() or (len(self.blat_results) == 1 and br.in_target) :
            self.logger.info('Blat result spans query (%r) or only one blat result (%r) and blat result in target (%r)'%(br.spans_query(), (len(self.blat_results) == 1), br.in_target))
            indel = True
            keep_br = br.valid and br.mean_cov < 2 and br.in_target and (br.indel_maxevent_size[0] >= indel_size_thresh) and (not br.rep_man.breakpoint_in_rep[0] and not br.rep_man.breakpoint_in_rep[1])
            self.logger.debug('Keep blat result %r'%keep_br)
            if keep_br :
                brkpt_cov = [self.meta_dict['contig_vals'][1].get_counts(x, x, 'indel') for x in br.query_brkpts]  
                low_cov = min(brkpt_cov) < self.meta_dict['params'].get_sr_thresh('indel')  
                flank_match_thresh = True
                for fm in br.indel_flank_match :
                    fm_perc = round((float(fm)/float(br.get_size('query')))*100,2)
                    if fm_perc < 10.0 : 
                        flank_match_thresh = False
                    self.logger.info('Indel result has matching flanking sequence of largest indel event of %d (%d of query)'%(fm,fm_perc))
                self.logger.info('Indel result has matching flanking sequence of largest indel event (10 perc of query) on both sides (%r)'%flank_match_thresh)
                in_ff, span_ff = filter_by_feature(br.get_brkpt_locs(), self.meta_dict['query_region'], self.meta_dict['params'].opts['keep_intron_vars'])
                if not in_ff and not low_cov and flank_match_thresh :
                    self.se = sv_event(br, self.meta_dict['query_region'], self.meta_dict['contig_vals'], self.meta_dict['sbam'])
                    self.logger.debug('Top hit contains whole query sequence, indel variant')
                else : 
                    self.logger.debug('Indel in intron (%r) or low coverage at breakpoints (%r) or minimum segment size < 20 (%r), filtering out.'%(in_ff, low_cov, min(br.query_blocksizes)) )
            else : 
                self.logger.debug('Indel failed checking criteria: in annotated gene: %r, mean query coverage < 2: %r, in target: %r, in repeat: %r, indel size < %d: %r'%(br.valid, br.mean_cov, br.in_target, ",".join([str(x) for x in br.rep_man.breakpoint_in_rep]), indel_size_thresh, br.indel_maxevent_size[0] < indel_size_thresh))
        return indel

    def get_indel_result(self):
        if self.se:
            return self.se.get_indel_result()
        else:
            return None

    def get_svs_result(self):
        if self.se :
            return self.se.get_svs_result(self.meta_dict['query_region'], self.meta_dict['params'], self.meta_dict['disc_reads']) 
        else :
            return None

    def check_indels(self):
        has_indel = False
#    blat_results_sorted = sorted(self.blat_results, key=lambda blat_results: (-blat_results[0], blat_results[1]) )
        for i in range(len(self.blat_results)): #blat_results_sorted)) :
            nmatch, ngaps, in_target, br, perc_ident = self.blat_results[i] #blat_results_sorted[i]
            br.mean_cov = self.get_mean_cov(br.qstart(), br.qend())
            #keep_clipped = (mean_cov<4 and ((br.get_nmatch_total()<30 and not br.in_repeat) or br.get_nmatch_total()>=30))
            #keep_clipped = keep_clipped or (br.in_target and mean_cov<10)
            #print nmatch, ngaps, br.mean_cov
            if i == 0 and self.check_blat_indel(br):
                has_indel = True
                self.logger.info('Contig has indel, returning %r'%has_indel)
                return has_indel
            else: #if keep_clipped :
                self.logger.debug('Storing clipped blat result start %d, end %d'%(br.qstart(), br.qend()))
                self.clipped_qs.append((br.qstart(),br.qend(),br,i))
        self.logger.info('Contig does not have indel, return %r'%has_indel)
        return has_indel

    def check_svs(self):
        self.logger.info('Checking for SVs')
        gaps = [(0, self.qsize)]
        if len(self.clipped_qs) > 1:
            self.logger.debug('Iterating through %d clipped blat results.'%len(self.clipped_qs))
            merged_clip = [0,None]
            for i in range(len(self.clipped_qs)) :
                qs, qe, br, idx = self.clipped_qs[i]
                self.logger.debug('Blat result with start %d, end %d, chrom %s'%(qs,qe,br.get_name('hit')))
                gaps = self.iter_gaps(gaps, self.clipped_qs[i], i) 
                if self.se.qlen > merged_clip[0] : # and self.se.in_target : 
                    merged_clip = [self.se.qlen,self.se]
            self.se = merged_clip[1]
        else :
            self.logger.info('There are no more than 1 clipped blat results, not continuing with SVs calling.')

        if self.se_valid() : 
            return True
        else:
            return False

    def se_valid(self) :
        valid = False
        if self.se and len(self.se.blat_res) > 1 and self.se.in_target :
            nmissing_query_cov = len(filter(lambda y: y, map(lambda x: x==0, self.se.query_cov)))
            if nmissing_query_cov < self.meta_dict['params'].get_min_segment_length('trl') : 
                valid = True 
        return valid

    def check_add_br(self, qs, qe, gs, ge, br) :
        self.logger.info('Checking to add blat result with start %d, end %d'%(qs, qe))
        add = False
        over_perc = round((float(min(qe,ge)-max(qs,gs)) / float(qe-qs)) * 100) # Calc % of segment overlaps with gap
        ov_right = 0 # Check overlap with other aligned segments
        if qe > ge : ov_right = abs(qe-ge)
        ov_left = 0
        if qs < gs : ov_left = abs(qs-gs)
        br.set_segment_overlap(ov_left, ov_right)
        max_seg_overlap = max(ov_right,ov_left)
        self.logger.debug('Blat query segment overlaps gap by %f'%over_perc)
        self.logger.debug('Max segment overlap %f'%max_seg_overlap)
        self.logger.debug('Event in target %r and blat result in target %r'%(self.se.in_target, br.in_target))
        if over_perc >= 50 and (max_seg_overlap < 15 or (br.in_target and self.se.in_target) ) : # and (self.se.in_target or br.in_target) : 
            add = True
        self.logger.debug('Add blat result to SV event %r'%add)
        return add

    def iter_gaps(self, gaps, cq, iter) :
        new_gaps = []
        qs, qe, br, idx = cq
        hit = False
        for gap in gaps :
            gs, ge = gap
            self.logger.debug('Gap coords %d, %d'%(gs, ge))
            if (qs >= gs and qs <= ge) or (qe <= ge and qe >= gs) :
                ngap = []
                if qs > gs : 
                    if (qs-1-gs) > 10 : 
                        ngap.append((gs,qs-1))
                if qe < ge : 
                    if (ge-qe+1) > 10 :
                        ngap.append((qe+1,ge))
                if iter == 0 : 
                    self.logger.debug('Creating SV event from blat result with start %d, end %d'%(qs, qe))
                    self.se = sv_event(br, self.meta_dict['query_region'], self.meta_dict['contig_vals'], self.meta_dict['sbam'])
                    new_gaps.extend(ngap)
                    hit = True
                elif self.check_add_br(qs, qe, gs, ge, br) :
                    self.logger.debug('Adding blat result to event')
                    new_gaps.extend(ngap)
                    self.se.add(br)
                    hit = True
                else :
                    new_gaps.append(gap)
            else :
                new_gaps.append(gap)
            self.logger.debug('New gap coords %s'%(",".join([str(x) for x in new_gaps])))
        if not hit :
            self.se.check_previous_add(br)
        return new_gaps


class ContigCaller:
    """
    """
    def __init__(self):

    def call_svs(self):
        contig_id = self.meta_dict['contig_vals'][2]
        if not self.bm.has_blat_results:
            self.logger.info('No blat results file %s, no calls for %s.' % (self.query_res_fn, contig_id))
        else:
            self.logger.info('Making variant calls from blat results %s'%self.query_res_fn)
            if self.bm.check_indels():
                self.result = self.bm.get_indel_result()
            elif self.bm.check_svs():
                self.result = self.bm.get_svs_result()
        return self.result


class align_manager:
    def __init__(self, meta_dict):
        self.meta_dict = meta_dict
        self.query_res_fn = meta_dict['query_res_fn']
        self.logger = logging.getLogger('root')
        self.result = None
        self.bm = blat_manager(meta_dict)

    def check_target_results(self):
        hit = False
        self.logger.info('Checking if target blat contains most of query or if whole genome needs to queried.')
        if not self.bm.has_blat_results:
            self.query_res_fn = None
            hit = True
        else:
            self.bm.write_mod_result_file(self.query_res_fn + '.mod')
            if self.bm.target_hit():
                hit = True
                self.logger.debug('Top hit contains whole query sequence, indel variant')
                self.query_res_fn += '.mod'
#      else :
#        if max(contig_counts.others) < params.get_sr_thresh('trl') :
#          hit = True
#          self.query_res_fn = None
        return hit, self.query_res_fn

    def get_result(self):
        contig_id = None
        if 'contig_vals' in self.meta_dict:
            contig_id = self.meta_dict['contig_vals'][2]
        if not self.bm.has_blat_results:
            self.logger.info('No blat results file %s, no calls for %s.' % (self.query_res_fn, contig_id))
        else:
            self.logger.info('Making variant calls from blat results %s'%self.query_res_fn)
            if self.bm.check_indels():
                self.result = self.bm.get_indel_result()
            elif self.bm.check_svs():
                self.result = self.bm.get_svs_result()
        return self.result

# def parse_config_f(config_fn, opts) :
#     param_opts = {}
#     config_f = open(config_fn, 'rU')
#     flines = config_f.readlines()
#     for line in flines :
#         line = line.strip()
#         linesplit = line.split("=")
#         if len(linesplit) == 1 : 
#             sys.exit()
#         else :
#             k,v = linesplit
#             param_opts[k] = v
#     config_f.close()
#     for opt in vars(opts) : 
#         param_opts[opt] = vars(opts)[opt]
#     return param_opts