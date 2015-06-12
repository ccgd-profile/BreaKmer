#! /usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


class BlatResult:
    def __init__(self, blatResultValues):
        self.values = blatResultValues
        self.matches = {}
        self.gaps = {}
        self.vals = {}
        self.fragments = {}
        self.strand = ''
        self.genes = ''
        self.in_target = False
        self.rep_man = None
        self.in_repeat = False
        self.repeat_overlap = 0.0
        self.repeat_coords = None
        self.filter_reps_edges = [False, False]
        self.valid = True
        self.interval = None
        self.breakpts = []
        self.query_brkpts = []
        self.query_blocksizes = []
        self.indel_sizes = []
        self.indel_maxevent_size = [0, '']
        self.mean_cov = 0.0
        self.perc_ident = 0.0
        self.seg_overlap = [0, 0]
        self.cigar = ''
        self.indel_flank_match = [0, 0]
        self.set_values()

    def set_values(self):
        self.matches['match'] = int(self.values[0])
        self.matches['mis'] = int(self.values[1])
        self.matches['rep'] = int(self.values[2])
        self.gaps['hit'] = [int(self.values[6]), int(self.values[7])]
        self.gaps['query'] = [int(self.values[4]), int(self.values[5])]

        tname = self.values[13].replace('chr', '')
        if 'tname' in d:
            tname = d['tname']
            self.values[13] = tname

        toffset = 0
        if 'offset' in d:
            toffset = d['offset']
        tcoords = [toffset + int(self.values[15]), toffset + int(self.values[16])]
        self.values[15] = tcoords[0]
        self.values[16] = tcoords[1]

        self.vals['hit'] = {'name': tname,
                            'size': int(self.values[14]),
                            'coords': tcoords}
        self.vals['query'] = {'name': self.values[9],
                              'size': int(self.values[10]),
                              'coords': [int(self.values[11]), int(self.values[12])]}
        self.strand = self.values[8]
        self.query_blocksizes = [int(x) for x in self.values[18].rstrip(",").split(",")]
        tstarts = [toffset + int(x) for x in self.values[20].rstrip(",").split(",")]
        self.values[20] = ",".join([str(x) for x in tstarts]) + ","
        qstarts = [int(x) for x in self.values[19].rstrip(",").split(",")]
        self.fragments['hit'] = []
        self.fragments['query'] = []
        for qstart, tstart, blocksize in zip(qstarts, tstarts, self.query_blocksizes):
            self.fragments['hit'].append((tstart, tstart + blocksize))
            self.fragments['query'].append((qstart, qstart + blocksize))
        self.fragments['count'] = len(self.query_blocksizes)

        if 'params' in d:
            self.set_gene_anno(d['params'].gene_annotations, d['query_region'])
            if 'repeat_mask' in d:
                self.set_repeat(d['repeat_mask'], d['params'].repeat_mask)

#    if 'gene_anno' in d : self.set_gene_anno(d['gene_anno'], d['query_region'])
#    if 'repeat_mask' in d : self.set_repeat(d['repeat_mask'],d['repeat_mask_all'])
        self.set_indel_locs()
        self.perc_ident = 100.0 - self.calcMilliBad()

    def calcMilliBad(self):
        bad = 0.0
        qAliSize = self.qend() - self.qstart()
        tAliSize = self.tend() - self.tstart()
        aliSize = min(qAliSize, tAliSize)
        if aliSize <= 0:
            return 0.0
        sizeDif = qAliSize - tAliSize
        if sizeDif < 0:
            sizeDif = 0
        insertFactor = self.gaps['query'][0]
        total = self.matches['match'] + self.matches['rep'] + self.matches['mis']
        if total != 0:
            bad = (1000 * (self.matches['mis'] + insertFactor + round(3 * math.log(1 + sizeDif)))) / total
        return bad * 0.1

    def set_segment_overlap(self, right, left):
        self.seg_overlap = [left, right]
            
    def set_repeat(self, target_rep_mask, all_rep_mask):
        self.rep_man = blat_repeat_manager()
        if self.matches['rep'] > 0 : 
            self.in_repeat = True
        if target_rep_mask and all_rep_mask :
            # Check rep_mask if it exists.
            rmask = target_rep_mask
            if not self.in_target : 
                rmask = None
                if self.vals['hit']['name'] in all_rep_mask : 
                    rmask = all_rep_mask[self.vals['hit']['name']]
            if rmask :
                self.rep_man.setup(self.get_coords('hit'), rmask) 
                self.in_repeat, self.repeat_overlap, self.repeat_coords, self.filter_reps_edges = self.rep_man.other_values 

    def get_coords(self,type) : return self.vals[type]['coords']
    
    def qstart(self) : return self.get_coords('query')[0]

    def qend(self) : return self.get_coords('query')[1]

    def tstart(self) : return self.get_coords('hit')[0]

    def tend(self) : return self.get_coords('hit')[1]
    
    def get_name(self,type) : return self.vals[type]['name']
    
    def get_size(self,type) : return self.vals[type]['size']

    def get_query_span(self) : 
        return self.get_coords('query')[1] - self.get_coords('query')[0]

    def get_query_coverage(self) : 
        return round((float(self.get_query_span())/float(self.get_size('query')))*100,2) 

    def spans_query(self) :
        return self.get_size('query') == (self.get_coords('query')[1]- self.get_coords('query')[0]) 

    def get_ngap_total(self) :
        return self.gaps['hit'][1] + self.gaps['query'][1]

    def get_num_gaps(self) :
        return self.gaps['hit'][0] + self.gaps['query'][0]

    def get_nmatch_total(self) :
        return self.matches['match'] + self.matches['rep']

    def get_nmatches(self,type) : return self.matches[type]

    def sum_indel_flank_matches(self,flank_str) :
        m_indxs = []
        match_sum = 0
        for i in range(len(flank_str)) :  
            if flank_str[i] == "M" : m_indxs.append(i)
        for indx in m_indxs :
            nmatch = ''
            windx = indx-1
            while windx > -1 and is_number(flank_str[windx]) :
                nmatch = flank_str[windx] + nmatch
                windx = windx - 1
            match_sum += int(nmatch)
        return match_sum

    def set_indel_flank_matches(self) :
        if self.indel_maxevent_size[0] > 0 :
            csplit = self.cigar.split(str(self.indel_maxevent_size[0]) + self.indel_maxevent_size[1])
            lflank = csplit[0]
            self.indel_flank_match[0] += self.sum_indel_flank_matches(lflank)
            rflank = csplit[-1]
            self.indel_flank_match[1] += self.sum_indel_flank_matches(rflank)

    def set_indel_locs(self) :
        chrm = 'chr'+str(self.get_name('hit'))
        for i in range(self.fragments['count']-1) :
            if i==0 and self.fragments['query'][i][0] > 0 : 
                self.cigar = str(self.fragments['query'][i][0]) + "S"
            qend1 = int(self.fragments['query'][i][1])
            qstart2 = int(self.fragments['query'][i+1][0])
            tend1 = int(self.fragments['hit'][i][1])
            tstart2 = int(self.fragments['hit'][i+1][0])
            ins_bp = qstart2 - qend1
            del_bp = tstart2 - tend1
            bp1 = tend1
            bp2 = tstart2
            self.cigar += str(self.query_blocksizes[i]) + "M"
            if ins_bp > 0 :
                self.breakpts.append([bp1])
                self.indel_sizes.append("I"+str(ins_bp))
                self.add_query_brkpt(qend1)
                self.add_query_brkpt(qstart2)
                self.cigar += str(ins_bp) + "I"
                if ins_bp > self.indel_maxevent_size[0] : self.indel_maxevent_size = [ins_bp,"I"]
            if del_bp > 0 :
                self.breakpts.append([bp1,bp2])
                self.indel_sizes.append("D"+str(del_bp))
                self.add_query_brkpt(qend1)
                self.cigar += str(del_bp) + "D"
                if del_bp > self.indel_maxevent_size[0] : self.indel_maxevent_size = [del_bp,"D"]
            
        self.cigar += str(self.query_blocksizes[-1]) + "M"
        end_clipped = self.get_size('query') - self.get_coords('query')[1] 
        if end_clipped > 0 : 
            self.cigar += str(end_clipped) + "S"

        self.set_indel_flank_matches()

        if self.strand == "-" :
            for i in range(len(self.query_brkpts)) :
                self.query_brkpts[i] = self.get_size('query') - self.query_brkpts[i]  

    def add_query_brkpt(self, brkpt) :
        if brkpt not in self.query_brkpts : 
            self.query_brkpts.append(brkpt)

    def get_brkpt_str(self, with_sizes=False) :
        brkpt_out = []
        bp_str = []
        chrm = 'chr' + str(self.get_name('hit'))
        if len(self.breakpts) > 0 :
            for b,s in zip(self.breakpts, self.indel_sizes) : 
                if len(b) > 1 : bb = "-".join([str(x) for x in b])
                else : bb = str(b[0])
                bstr = chrm + ":" + bb
                if with_sizes : 
                    bstr += " " + "(" + s + ")"
                bp_str.append(bstr)
            brkpt_out.append(",".join(bp_str))
        return ",".join(brkpt_out)
    
    def get_brkpt_locs(self) :
        brkpt_locs = []
        for b in self.breakpts :
            brkpt_locs.extend(b)
        return brkpt_locs

    def get_gene_anno(self) : return self.genes

    def get_blat_output(self) :
        return "\t".join([str(x) for x in self.values])

    def get_len(self) :
        return self.qend - self.qstart

    def set_gene_anno(self, annotations, query_region) :
        br_start = self.get_coords('hit')[0]
        br_end = self.get_coords('hit')[1]
        start_in = br_start >= (query_region[1]-200) and br_start <= (query_region[2]+200)
        end_in = br_end <= (query_region[2]+200) and br_end >= (query_region[1]-200)
        if query_region[0] == self.get_name('hit') and (start_in or end_in) :
            self.in_target = True
            self.genes = query_region[3]
        else :
            ann_genes = []
            chrom = self.get_name('hit')
            pos = self.get_coords('hit')
            if chrom.find('chr') == -1 : chrom = 'chr'+str(chrom)
            for g in annotations.genes :
                gs = annotations.genes[g][1]
                ge = annotations.genes[g][2]
                if chrom == annotations.genes[g][0] :
                    if int(pos[0]) >= gs and int(pos[0]) <= ge :
                        ann_genes.append(g)
                        break
            if len(ann_genes) == 0 :
                ann_genes = ['intergenic']
                self.valid = False
            self.genes = ",".join(ann_genes)