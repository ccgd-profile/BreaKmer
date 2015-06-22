"""

"""

from glob import glob
import os
import sys
from os.path import join as jp


class filter():
    def __init__(self, targetName, svType, descriptions) :
        self.name = targetName
        self.type = svType
        self.descriptions = descriptions.split('|')

    def check_rec(self, rec) :
        filterCheck = (rec.type == self.type) and (rec.targetName.find(self.name) > -1)
        descMatch = False
        if self.descriptions[0] != '*' :
            for desc in self.descriptions :
                if rec.svDescription :
                    descMatch = descMatch or rec.svDescription.find(desc) > -1
        else :
            descMatch = True

        return (filterCheck and descMatch)


class filters() :
    def __init__(self, filterFn) :
        self.filterFn = filterFn
        self.filters = []
        self.setup()

    def setup(self) :
        filterF = open(filterFn, 'r') 
        for line in filterF.readlines() :
            line = line.strip()
            targetName, svType, descriptions = line.split(" ")
            self.filters.append(filter(targetName, svType, descriptions))

    def check_rec(self, rec) :
        filterCheck = False
        for filter in self.filters :
            filterCheck = filterCheck or filter.check_rec(rec)

        keepIndel = True
        if rec.type == "indel" :
            #print rec.values
            #print rec.segments[0].annoVals, abs(int(rec.segments[0].annoVals['exon'][1])) 
            #print rec.svDescription, rec.svDescription.find(',')
            if abs(int(rec.segments[0].annoVals['exon'][1])) > 0 and rec.svDescription.find(',') == -1 :
                keepIndel = False    
            #print keepIndel

        return filterCheck or (not keepIndel)


class segment() :
    def __init__(self, recId, segId, svType, geneName, brkpt, svDescription, cigar, mm, strand, alignMetrics, srCount, covCount) :
        self.recId = recId
        self.segId = segId
        self.svType = svType
        self.geneName = geneName
        self.chr = brkpt.split(":")[0]
        self.bp = brkpt.split(":")[1]
        self.indelDesc = svDescription
        self.cigar = cigar
        self.mm = mm
        self.strand = strand
        self.repOverlap = alignMetrics.split(":")[0]
        self.segLen = alignMetrics.split(":")[1]
        self.uniqAlign = None
        self.srCount = srCount
        self.covCount = covCount
        self.annoVals = None
        if self.svType.find("rearr") > -1 :
            self.uniqAlign = alignMetrics.split(":")[2]

    def get_meta_vals(self) :
        # Segment meta values encode the annotation and realigment values
        # s1=chr+cytoband:bp:gene_name:distance:strand:repeat_overlap:uniq_align_metric, s2=...
        segId = 'seg' + str(self.segId+1)
        segLoc = self.chr + self.annoVals['cyto'][0]
        segVals = [segId+'='+segLoc, self.bp]
        geneVals = ''
        geneGeneStatus = 'gene-gene'
        if 'cluster' in self.annoVals :
            geneVals += ":".join(self.annoVals['cluster'])
            segGeneNames = self.annoVals['cluster'][0]
            geneGeneStatus = 'gene-gene'
        else :
            if abs(self.annoVals['gene_up'][1]) > 0 :
                # Segment lands in another gene
                geneVals += ":".join([str(x) for x in self.annoVals['gene_up']])
                geneVals += '|' + ":".join([str(x) for x in self.annoVals['gene_down']])
                segGeneNames = "|".join([self.annoVals['gene_up'][0], self.annoVals['gene_down'][0]])
                geneGeneStatus = 'gene-intergenic'
            else :
                geneVals += ":".join([str(x) for x in self.annoVals['gene_up']])
                segGeneNames = self.annoVals['gene_up'][0]
                geneGeneStatus = 'gene-gene'
        segVals.append(geneVals)
        segVals.append(self.repOverlap)
        segVals.append(self.uniqAlign)
        return segVals, segGeneNames, geneGeneStatus

    def get_gene_name(self) :
        geneName = ''
        if 'cluster' in self.annoVals :
            geneName = self.annoVals['cluster'][0]
        else :
            if abs(self.annoVals['gene_up'][1]) > 0 :
                # Segment lands in another gene
                geneName = "|".join([self.annoVals['gene_up'][0], self.annoVals['gene_down'][0]])
            else :
                geneName = self.annoVals['gene_up'][0]
        return geneName

class breakmerRec() :
    def __init__(self, subProjId, sampleId, recIter, values) :
        '''
        1. genes - comma-delimited list
        2. target_breakpoints - comma delimited list
        3. align_cigar  - comma delimited list
        4. mismatches - comma delimited list
        4. strands - comma delimited list
        5. rep_overlap_segment_len - comma delimited list repeat overlap percentage:segment length:uniqueness score
        6. sv_type - string value
        7. split_read_count - comma delimited list entry for each breakpoint
        8. nkmers - integer
        9. disc_read_count - integer value
        10. breakpoint_coverages - comma delimited list entry for each breakpoint
        11. contig_id - string <target name>_<contig id>
        12. contig_seq - contig sequence
        '''
        self.subProjId = subProjId
        self.sampleId = sampleId
        self.recId = 'B' + str(recIter)
        self.values = values
        self.genes = values[0].split(',')
        self.brks = [values[1]]
        self.cigars = values[2].split(',')
        self.mismatches = [int(x) for x in values[3].split(',')]
        self.strands = values[4].split(',')
        self.repOverlapSegLens = values[5].split(',')
        self.type = values[6]
        self.splitReadCounts = [int(x) for x in values[7].split(',')]
        self.kmerCount = int(values[8])
        self.discReadCounts = [int(x) for x in values[9].split(',')]
        self.coverages = [int(x) for x in values[10].split(',')]
        self.contigId = values[11].split('_')[1]
        self.targetName = values[11].split('_')[0]
        self.contigSeq = values[12]
        self.segments = []
        self.svDescription = None
        self.format_results()

    def get_target_name(self) :
        # Suffers from annotation issue, if breakpoint lands outside of any
        # target capture regions it won't hit any of the target regions and
        # therefore will not have a target name associated.
        # For now this function is deprecated
        tName = ''
        for seg in self.segments :
            if 'target' in seg.annoVals :
                tName = seg.annoVals['target'][0]
        return tName

    def get_cytobands(self) :
        cytoNames = []
        for seg in self.segments :
            cytoNames.append(seg.annoVals['cyto'][0])
        return ','.join(cytoNames)

    def format_trl_str(self) :
        outStr = ''
        chrs = []
        bands = []
        for seg in self.segments :
            chrs.append(seg.chr)
            bands.append(seg.annoVals['cyto'][0])
        outStr = 't(' + ';'.join(chrs) + ')(' + ';'.join(bands) + ')'
        return outStr

    def get_segment_meta(self) :
        segMeta = []
        segGeneNames = []
        geneGeneStatus = ''
        for seg in self.segments :
            segVals, geneNames, geneGeneStatus = seg.get_meta_vals()
            segGeneNames.append(geneNames)
            '''
            # Segment meta values encode the annotation and realigment values
            # s1=chr+cytoband:bp:gene_name:distance:strand:repeat_overlap:uniq_align_metric, s2=...
            segVals = ['seg' + str(seg.segId+1) + '=' + seg.chr + seg.annoVals['cyto'][0], seg.bp]
            geneVals = ''
            geneGeneStatus = 'gene-gene'
            if 'cluster' in seg.annoVals :
                print seg.annoVals['cluster']
                geneVals += ":".join(seg.annoVals['cluster'])
                segGeneNames.append(seg.annoVals['cluster'][0])
                geneGeneStatus = 'gene-gene'
            else :
                if abs(seg.annoVals['gene_up'][1]) > 0 :
                    # Segment lands in another gene
                    geneVals += ":".join([str(x) for x in seg.annoVals['gene_up']])
                    geneVals += '|' + ":".join([str(x) for x in seg.annoVals['gene_down']])
                    segGeneNames.append("|".join([seg.annoVals['gene_up'][0], seg.annoVals['gene_down'][0]]))
                    geneGeneStatus = 'gene-intergenic'
                else :
                    geneVals += ":".join([str(x) for x in seg.annoVals['gene_up']])
                    segGeneNames.append(seg.annoVals['gene_up'][0])
                    geneGeneStatus = 'gene-gene'
            segVals.append(geneVals)
            segVals.append(seg.repOverlap)
            segVals.append(seg.uniqAlign)
            '''
            segMeta.append(":".join([str(x) for x in segVals]))
        segGeneNames.sort()
        return ",".join(segMeta), ','.join(segGeneNames), geneGeneStatus

    def get_circos_vals(self):
        targetSeg = [None,0]
        iter = 0
        for seg in self.segments:
            if 'cluster' in seg.annoVals and self.targetName.find(seg.annoVals['cluster'][0]) > -1: 
                targetSeg = [seg, iter]
                break 
            elif (seg.annoVals['gene_up'][1]==0) and (self.targetName.find(seg.annoVals['gene_up'][0]) > -1):
                targetSeg = [seg, iter]
                break
            elif self.targetName.find(seg.annoVals['gene_up'][0]) > -1 or self.targetName.find(seg.annoVals['gene_down'][0]) > -1:
                targetSeg = [seg, iter]
                break
            iter += 1

        if not targetSeg[0]:
            return None, None
        targetBp = targetSeg[0].bp.split('-')[0]
        targetLink = "\t".join(['hs' + targetSeg[0].chr.replace('chr',''), targetBp, targetBp])
        links = []
        anno = ["\t".join(['hs' + targetSeg[0].chr.replace('chr',''), targetBp, targetBp, self.targetName])]
        for seg in self.segments :
            if targetSeg[0] == seg:
                continue
            bp = seg.bp
            if bp.find('-') > -1:
                bp = bp.split('-')[0]
            anno.append("\t".join(['hs' + seg.chr.replace('chr',''), bp, bp, seg.get_gene_name()]))
            links.append(targetLink + "\t" + "\t".join(['hs' + seg.chr.replace('chr',''), bp, bp, self.sampleId]))
        return links, anno 

    '''
    def get_trl_partner_genes(self) :
        partnerGeneVals = []
        partnerGeneNames = []
        for seg in self.segments :
            segVals = []
            # Segment is outside of target region
            segVals.append(seg.chr + seg.annoVals['cyto'][0])
            geneVals = ''
            if 'cluster' in seg.annoVals :
                geneVals += ":".join(seg.annoVals['cluster']) 
                partnerGeneNames.append(seg.annoVals['cluster'][0])
            else:
                if abs(seg.annoVals['gene_up'][1]) > 0 :
                    # Segment lands in another gene
                    geneVals += ":".join([str(x) for x in seg.annoVals['gene_up']])
                    geneVals += '|' + ":".join([str(x) for x in seg.annoVals['gene_down']])
                    partnerGeneNames.append("|".join([seg.annoVals['gene_up'][0], seg.annoVals['gene_down'][0]]))
                else :
                    geneVals += ":".join([str(x) for x in seg.annoVals['gene_up']])
                    partnerGeneNames.append(seg.annoVals['gene_up'][0])
            segVals.append(geneVals)
            partnerGeneVals.append('s' + str(seg.segId) + "=" + ":".join(segVals))
        partnerGeneNames.sort()
        return ",".join(partnerGeneVals), ','.join(partnerGeneNames)
    '''

    def get_segment_strands(self) :
        strands = []
        for seg in self.segments :
            strands.append(seg.strand)
        return ','.join(strands)

    def get_segment_brkpts(self) :
        brkpts = []
        for seg in self.segments :
            brkpts.append(seg.chr + ':' + seg.bp)
        return ','.join(brkpts)

    def filter_result(self, filterList) :
        filter = False
        primaryChrs = ['chr' + str(x) for x in range(1,23)] + ['chrX', 'chrY']
        # check if non-primary chromosomes
        for seg in self.segments :
            if seg.chr not in primaryChrs :
                filter = True
                break
        filter = filter or filterList.check_rec(self)
        return filter

        '''
        if self.type == 'indel' :
            aurka_indel = self.targetName == 'AURKA' and self.svDescription == '(I61)'
            aurkc_indel = self.targetName == 'AURKC' and self.svDescription == '(D14,I1)'
            notch2_indel = self.targetName == 'NOTCH2' and self.svDescription == '(I20)'
            flt4_indel = self.targetName == 'FLT4' and self.svDescription == '(I22,D2)'
            ptprd_indel = self.targetName == 'PTPRD' and self.svDescription == '(I21)'
            hla_indel = self.targetName.find('HLA') > - 1

            if aurka_indel or aurkc_indel or notch2_indel or flt4_indel or ptprd_indel or hla_indel :
                filter = True
        return filter
        '''

    def get_rep_overlaps(self) :
        repOverlaps = []
        for seg in self.segments :
            repOverlaps.append(float(seg.repOverlap))
        return repOverlaps

    def get_uniq_align(self) :
        uniqAligns = []
        for seg in self.segments :
            uniqAligns.append(float(seg.uniqAlign))
        return uniqAligns

    def get_report(self, filterList) :
#        targetGeneName = self.get_target_name()
        cytoBands = self.get_cytobands()
        descStr = ''
        segMeta = '' 
        segGenes = ''
        if self.type.find('rearr') > -1 :
            descStr = self.format_trl_str()
            self.svDescription = descStr
            # geneMeta, geneList = self.get_trl_partner_genes()
            segMeta, segGenes, geneGeneStatus = self.get_segment_meta()
        elif self.type == 'indel' :
            descStr = self.svDescription
            segMeta = self.segments[0].chr + self.segments[0].annoVals['cyto'][0] + ':' + self.segments[0].bp
            geneGeneStatus = 'NA'
            segGenes = self.targetName + ":" + self.segments[0].annoVals['exon'][0] + ":" + self.segments[0].annoVals['exon'][1]

        recKey = None
        lout = [self.subProjId, self.sampleId, self.type, descStr, self.get_segment_strands(), self.targetName, segGenes, geneGeneStatus, ','.join([str(x) for x in self.splitReadCounts]), ','.join([str(x) for x in self.discReadCounts]), segMeta, self.contigSeq]
        brkpts = self.get_segment_brkpts()
        recKey = self.targetName + '|' + descStr + '|' + brkpts

        if self.filter_result(filterList) :
            lout.append('filtered')

        circosVals = {'links':[], 'annos':[]}
        if self.type == 'rearrangement' :
            links, anno = self.get_circos_vals()
            if links and anno:
                circosVals['links'] = links
                circosVals['annos'] = anno
        lout = "\t".join(lout)
        
        return lout, recKey, circosVals

    def get_all_dataframe(self, filterList) :
        cytoBands = self.get_cytobands()
        descStr = ''
        segMeta = '' 
        segGenes = ''
        maxUniqAlign = 'NA'
        minUniqAlign = 'NA'
        maxRepOverlap = 'NA'

        if self.type.find('rearr') > -1 :
            descStr = self.format_trl_str()
            self.svDescription = descStr
            segMeta, segGenes, geneGeneStatus = self.get_segment_meta()
            maxUniqAlign = max(self.get_uniq_align())
            minUniqAlign = min(self.get_uniq_align())
            maxRepOverlap = max(self.get_rep_overlaps())
        elif self.type == 'indel' :
            descStr = self.svDescription
            descStr = descStr.lstrip('(')
            descStr = descStr.rstrip(')')

        lout = None
        recKey = None

        if not self.filter_result(filterList) :
            lout = [self.subProjId, self.sampleId, self.type, descStr, self.get_segment_strands(), self.targetName, segGenes, max(self.splitReadCounts), max(self.discReadCounts), max(self.coverages), maxRepOverlap, maxUniqAlign, minUniqAlign, len(self.contigSeq)]
            brkpts = self.get_segment_brkpts()
            recKey = self.targetName + '|' + descStr + '|' + brkpts
            lout = "\t".join([str(x) for x in lout])
        return lout, recKey

    def format_results(self) :
        if self.targetName == 'MLL':
            self.targetName = 'KMT2A'
        if self.targetName == 'MLL3':
            self.targetName = 'KMT2C'
        genes = self.genes
        if self.type.find("rearr") > -1 :
            self.brks = self.brks[0].split(",")
            if len(self.genes) != len(self.brks) :
                genes = [self.genes] * len(self.brks)
        else :
            self.brks[0], self.svDescription = self.brks[0].split(' ')
            self.svDescription = self.svDescription.lstrip('(').rstrip(')')
        for i in range(len(self.brks)) :
            self.segments.append(segment(self.recId, i, self.type, genes[i], self.brks[i], self.svDescription, self.cigars[i], self.mismatches[i], self.strands[i], self.repOverlapSegLens[i], self.splitReadCounts[i], self.coverages[i]))

    def write_brkpt_bed(self, fbed) :
        for segment in self.segments :
            bedline = [segment.chr]
            bp1 = ''
            bp2 = ''
            if segment.bp.find('-') > -1 : 
                bp1,bp2 = segment.bp.split('-')
            else :
                bp1 = int(segment.bp) - 1
                bp2 = int(segment.bp) + 1
            if bp1 > bp2 : 
                tmp = bp1
                bp1 = bp2
                bp2 = tmp
            bedline.extend([bp1, bp2, self.sampleId + ":" + segment.recId+":s"+str(segment.segId), '0', segment.strand])
            fbed.write("\t".join([str(x) for x in bedline]) + "\n")

#------------------------------------------------------------------------------
def process_res(subProjId, sampleId, fin, resD, recIter) :
    res = open(fin, 'rU')
    resLines = res.readlines()[1:]
    for line in resLines :
        line = line.strip()
        linesplit = line.split("\t")
        bRec = breakmerRec(subProjId, sampleId, recIter, linesplit)
        resD.append(bRec)
        recIter += 1
    return recIter

#------------------------------------------------------------------------------
def process_bed_anno(f, d, tag) :
    '''
    chr1    2493182 2493215 B1448_s0    0   +   chr1    2493111 2493254 TNFRSF14    exon    33
    chr1    11300368    11300370    B5258_s1    0   +   chr1    11300359    11300604    MTOR    exon    2
    chr1    11850950    11850951    B912_s0 0   +   chr1    11850736    11850955    MTHFR   exon    1
    '''
    for line in open(f, 'r').readlines() :
        line = line.strip()
        linesplit = line.split('\t')
        sampleId, brkptId, segId = linesplit[3].split(":")
        overlapId = linesplit[9]
        metaVals = [None, None]
        if tag == 'target' :
            metaVals[0] = linesplit[10]
        elif tag == 'cluster' :
            metaVals[0] = '0'
            metaVals[1] = linesplit[11]
        elif tag.find('gene') > -1 :
            metaVals[0] = int(linesplit[-1])
            metaVals[1] = linesplit[11]
            # print 'gene', linesplit[11], linesplit[-1]
        elif tag == 'exon' :
            metaVals[0] = linesplit[-1]
            metaVals[1] = linesplit[11]

        key = ":".join([sampleId, brkptId, segId])
        if key not in d :
            d[key] = {}
        if tag not in d[key]:
            v = [overlapId] + metaVals 
            if tag.find('gene') > -1:
                d[key][tag] = [v]
            else:
                d[key][tag] = v
        else:
            if tag.find('gene') > -1 :
                d[key][tag].append([overlapId] + metaVals)
    return d

#------------------------------------------------------------------------------
def annotate_recs(resD, annoDictFs) :
    annoIndx = {}
    for item in annoDictFs :
        annoIndx = process_bed_anno(annoDictFs[item], annoIndx, item)

    for key in annoIndx :
        sampleId, brkptId, segId = key.split(":")
        annoVals = annoIndx[key]
        rec = resD[sampleId]['recs'][int(brkptId.lstrip('B'))-1]
        seg = rec.segments[int(segId.lstrip('s'))]
        seg.annoVals = annoVals
        # Select only one annotated gene if there are two or more.
        # This is for instances where the target region is not 
        # designated.

        if len(seg.annoVals['gene_up']) > 1 :
            for vals in seg.annoVals['gene_up'] :
                if vals[0] == rec.targetName :
                    seg.annoVals['gene_up'] = vals
                    break
        if seg.annoVals['gene_up'][0].__class__ == list:
            seg.annoVals['gene_up'] = seg.annoVals['gene_up'][0]

        if len(seg.annoVals['gene_down']) > 1 :
            for vals in seg.annoVals['gene_down'] :
                if vals[0] == rec.targetName :
                    seg.annoVals['gene_down'] = vals
                    break
        if seg.annoVals['gene_down'][0].__class__ == list:
            seg.annoVals['gene_down'] = seg.annoVals['gene_down'][0]


#------------------------------------------------------------------------------
def report_df(resD, fout, filterList) :
    '''
    '''
    header = ['subprojectid','sampleid', 'rearr_type', 'event_desc', 'realign_strands', 'target_name', 'gene_list', 'max_splitread_counts', 'max_discread_counts', 'max_read_coverage', 'max_repeat_overlap', 'max_uniq_align', 'min_uniq_align', 'len_contig']
    fout.write("\t".join(header) + "\n")
    for sid in resD :
        recs = resD[sid]['recs']
        recStor = {}
        nsampleRecs = 0
        for rec in recs :
            recData, recKey = rec.get_all_dataframe(filterList)
            if recData :
                if recKey not in recStor :
                    # Store the key to eliminate duplicate calls
                    recStor[recKey] = recData
                    fout.write(recData + "\n")
                    nsampleRecs += 1
        if nsampleRecs == 0 :
            lout = resD[sid]['values'] + ['NA']*(len(header)-2)
            fout.write("\t".join(lout) + "\n")
    fout.close()
#------------------------------------------------------------------------------

def add_circos_values(circosValues, circosVals) :
    circosValues['links'].extend(circosVals['links'])
    for anno in circosVals['annos'] :
        if anno not in circosValues['annos'] :
            circosValues['annos'].append(anno)

#------------------------------------------------------------------------------
def report(resD, fileDict, filterList) :
    '''
    '''

    header = ['subprojectid', 'sampleid', 'rearr_type', 'event_description', 'realign_strands', 'target_name', 'gene_list', 'gene-gene_type', 'splitread_counts', 'discordantread_counts', 'segment_values(realign_chr_cytoband:bp:refseq_annotated_gene:distance to annotated_gene:realign_strand)', 'contig_sequence', 'filtered']
    legend = {}
    legend['segment_values'] = 'Comma delimited string with values relating to the realignment of the contig segments - realign_chr_cytoband:bp:refseq_annotated_gene:distance to annotated_gene (0=within, negative=upstream gene, positive=downstream gene):realign_strand:repeat_overlap_percentage (overlap with annotated repeats):uniq_alignment_metric (lower = more unique alignment)'
    legend['gene_list'] = 'Comma delimited string with RefSeq gene names that the segments either map directly to or are inbetween. A breakpoint that is intergenic will have an upstream gene and downstream gene separated by "|"'
    legend['splitread_counts'] = 'Comma delimited string with each value indicating the number of splitreads supporting a breakpoint.'
    legend['discordantread_counts'] = 'Number of read pairs where one read maps upstream of the breakpoint(s) and the other read maps downstream.'
    legend['Software'] = 'BreaKmer v0.0.6 (Abo, RP, Ducar, M, Garcia, EP, Thorner, AR, Rojas-Rudilla, V, Lin, L, Sholl, LM, Hahn, WC, Meyerson, M, Lindeman, NI, Van Hummelen, P, MacConaill, LE (2015). BreaKmer: detection of structural variation in targeted massively parallel sequencing data using kmers. Nucleic Acids Res., 43, 3:e19.)'

    reportOut = fileDict['report']
    filteredOut = fileDict['filtered']
    sampleOut = fileDict['sample']
    targetOut = fileDict['target']
    circosLinkOut = fileDict['link']
    circosAnnoOut = fileDict['anno']

    reportOut.write("\t".join(header) + "\n")
    filteredOut.write("\t".join(header) + "\n")

    circosValues = {'links':[], 'annos':[]}
    summary_values = {'nsamples':0, 'type_counts':{}, 'samples_with_event':[]}
    for sid in resD :
        summary_values['nsamples'] += 1
        recs = resD[sid]['recs']
        recStor = {}
        nsample_recs = 0
        for rec in recs :
            recReport, recKey, circosVals = rec.get_report(filterList)
            if recKey not in recStor :
                # Store the key to eliminate duplicate calls
                recStor[recKey] = recReport
                if recReport.find('filtered') == -1 :
                    reportOut.write(recReport + "\n")
                    add_circos_values(circosValues, circosVals)
                    if rec.type not in summary_values['type_counts']:
                        summary_values['type_counts'][rec.type] = 0
                    summary_values['type_counts'][rec.type] += 1
                    if sid not in summary_values['samples_with_event'] :
                        summary_values['samples_with_event'].append(sid)
                    nsample_recs += 1
                else :
                    filteredOut.write(recReport + "\n")
        if nsample_recs == 0 :
            reportOut.write('\t'.join(resD[sid]['values']) + "\n")
        # print sid, len(recs), nsample_recs
    for item in legend :
        reportOut.write(item + ': ' + legend[item] + "\n")

    for item in summary_values :
        if item == 'nsamples' :
            print 'Number of samples:', summary_values[item]
        elif item == 'type_counts' :
            for type in summary_values[item] :
                print 'Number of', type, 'events:', summary_values[item][type]
        elif item == 'samples_with_event' :
            print 'Number of samples with event:', len(summary_values['samples_with_event'])

    for link in circosValues['links']:
        circosLinkOut.write(link + '\n')
    for annos in circosValues['annos']:
        circosAnnoOut.write(annos + '\n')

    reportOut.close()
    filteredOut.close()
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def collapse_bed(targetBedFn) :
    collapsedBedPath = os.path.dirname(targetBedFn)
    fn = os.path.basename(targetBedFn).split('.bed')[0]
    collapsedBedFn = os.path.join(collapsedBedPath, fn + '.collapsed.bed')
    newBed = open(collapsedBedFn, 'w')

    geneRegion = [None, None, None, None]
    for line in open(targetBedFn, 'r').readlines() :
        line = line.strip()
        linesplit = line.split("\t")
        chr, start, stop, name, feature = linesplit
        if not geneRegion[3] :
            geneRegion = [chr, int(start), int(stop), name]
        elif name != geneRegion[3] :
            newBed.write("\t".join([str(x) for x in geneRegion]) + "\n")
            geneRegion = [chr, int(start), int(stop), name]
        else :
            if int(start) < geneRegion[1] :
                geneRegion[1] = int(start)
            if int(stop) > geneRegion[2] :
                geneRegion[2] = int(stop)
    newBed.write("\t".join([str(x) for x in geneRegion]) + "\n")
    newBed.close()
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def run_bedtools(brkptBedFn) :
    cytoBandBedFn = '/data/ccgd/reference/human/refseq/GRCh37/annotation/bed/cytoBandIdeo.primary_chrs.bed'
    geneAnnoBedFn = '/data/ccgd/reference/human/refseq/GRCh37/annotation/bed/refseq_refFlat.bed'
    canonTrxAnnoBedFn = '/data/ccgd/reference/human/ensembl/GRCh37-75/annotation/bed/ensembl_GRCh37_canon_trx_anno.v75.primary_chrs.sorted.bed'
    targetBedFn = '/data/ccgd/breakmer/ccgd_projects/panel_beds/P50.bed'
    clusterRegionBedFn = '/data/ccgd/breakmer/ccgd_projects/panel_beds/cluster_regions.bed'
    targetCollapsedBedFn = collapse_bed(targetBedFn)

    brkptTrgtAnno = 'brkpt_target_anno.bed'
    brkptGeneAnnoUp = 'brkpt_gene_anno_upstream.bed'
    brkptGeneAnnoDown = 'brkpt_gene_anno_downstream.bed'
    brkptCytoBand = 'brkpt_cytoband_anno.bed'
    brkptGeneAnnoIntersect = 'brkpt_gene_anno_intersect.bed'
    brkptExonAnnoIntersect = 'brkpt_exon_anno_intersect.bed'
    brkptClusterAnnoIntersect = 'brkpt_cluster_anno_intersect.bed'

    annoDict = {'target': brkptTrgtAnno, 'gene_int': brkptGeneAnnoIntersect, 'gene_up': brkptGeneAnnoUp, 'gene_down': brkptGeneAnnoDown, 'cyto': brkptCytoBand, 'exon': brkptExonAnnoIntersect, 'cluster': brkptClusterAnnoIntersect}

    # Identify segments that overlap the targeted regions
    cmd = 'sort -k1,1 -k2,2n %s | bedtools intersect -wo -a - -b %s > %s' % (brkptBedFn, targetBedFn, 'brkpt_target_anno.bed')
    os.system(cmd)
    # Intersect with gene annotation bed file
    cmd = 'sort -k1,1 -k2,2n %s | bedtools intersect -wo -a - -b %s > %s' % (brkptBedFn, geneAnnoBedFn, 'brkpt_gene_anno_intersect.bed')
    os.system(cmd)
    # Bedtools upstream genes
    cmd = 'sort -k1,1 -k2,2n %s | bedtools closest -D a -id -a - -b %s > %s' % (brkptBedFn, geneAnnoBedFn, 'brkpt_gene_anno_upstream.bed') 
    os.system(cmd)
    # Bedtools downstream genes
    cmd = 'sort -k1,1 -k2,2n %s | bedtools closest -D a -iu -a - -b %s > %s' % (brkptBedFn, geneAnnoBedFn, 'brkpt_gene_anno_downstream.bed') 
    os.system(cmd)
    # Cytoband 
    cmd = 'sort -k1,1 -k2,2n %s | bedtools closest -a - -b %s > %s' % (brkptBedFn, cytoBandBedFn, 'brkpt_cytoband_anno.bed') 
    os.system(cmd)

    # Cluster regions 
    cmd = 'sort -k1,1 -k2,2n %s | bedtools intersect -wo -a - -b %s > %s' % (brkptBedFn, clusterRegionBedFn, 'brkpt_cluster_anno_intersect.bed') 
    os.system(cmd)

    # Intersect with Ensembl canonical trxs
    cmd = 'sort -k1,1 -k2,2n %s | bedtools closest -D a -a - -b %s > %s' % (brkptBedFn, canonTrxAnnoBedFn, 'brkpt_exon_anno_intersect.bed')
    os.system(cmd)
    return annoDict
#------------------------------------------------------------------------------

def format_results(sampleListFn, sampleType, filterList) :
    resD = {}
    for line in open(sampleListFn, 'rU') :
        line = line.strip()
        subProjId, sampleId = line.split("\t")
        outFs = glob(jp(resultsDir, sampleId + '*svs.out'))
        resD[sampleId] = {'recs':[], 'values':[]}
        recIter = 1
        resD[sampleId]['values'] = [subProjId, sampleId]
        for outF in outFs :
            recIter = process_res(subProjId, sampleId, outF, resD[sampleId]['recs'], recIter)

    fbed = open('brkpts.%s.bed' % sampleType, 'w')
    for sid in resD :
        for rec in resD[sid]['recs'] :
            rec.write_brkpt_bed(fbed)
    fbed.close()

    annoDictFs = run_bedtools('brkpts.%s.bed' % sampleType)
    annotate_recs(resD, annoDictFs)

    if sampleType == "normal" :
        for sid in resD :
            for rec in resD[sid]['recs'] :
                descriptions = rec.svDescription
                if rec.type.find("rearr") > -1:
                    descriptions = rec.format_trl_str()
                filterList.filters.append(filter(rec.targetName, rec.type, descriptions))
    return resD

#------------------------------------------------------------------------------
if __name__ == "__main__" :
    args = sys.argv[1:]
    projectId = args[0]
    tumorSamplesFn = args[1]
    normalSamplesFn = args[2]

    resultsDir = os.path.join('/data/ccgd/breakmer/ccgd_projects/', projectId + '_breakmer')

#    sampleListFn = os.path.join(resultsDir, 'sample_list.txt')
    filterFn = os.path.join(resultsDir, 'filters.txt')
    filterList = filters(filterFn)

    resD = format_results(os.path.join(resultsDir, tumorSamplesFn), 'tumor', filterList)
    if normalSamplesFn != '' :
        normals = format_results(os.path.join(resultsDir, normalSamplesFn), 'normal', filterList)

#    sampleListF = open(sampleListFn, 'rU')
    '''
    resD = {}
    for line in sampleListF.readlines() :
        line = line.strip()
        subProjId, sampleId = line.split("\t")
        outFs = glob(jp(resultsDir, sampleId + '*svs.out'))
        resD[sampleId] = {'recs':[], 'values':[]}
        recIter = 1
        resD[sampleId]['values'] = [subProjId, sampleId]
        for outF in outFs :
            recIter = process_res(subProjId, sampleId, outF, resD[sampleId]['recs'], recIter)

    fbed = open('brkpts.bed', 'w')
    for sid in resD :
        for rec in resD[sid]['recs'] :
            rec.write_brkpt_bed(fbed)
    fbed.close()

    filterList = filters(filterFn)
    annoDictFs = run_bedtools('brkpts.bed')
    annotate_recs(resD, annoDictFs)
    '''

    dataFrameFn = '%s_breakmer.df.txt' % projectId
    reportFn = '%s_breakmer.txt' % projectId
    filteredFn = '%s_breakmer.filtered.txt' % projectId
    sampleSummaryFn = '%s_breakmer.sample_summary.txt' % projectId
    targetSummaryFn = '%s_breakmer.gene_summary.txt' % projectId
    circosLinkFn = '%s_circos_links.txt' % projectId
    circosAnnoFn = '%s_circos_annos.txt' % projectId

    outFileDict = {'report':reportFn, 'df':dataFrameFn, 'filtered':filteredFn, 'sample':sampleSummaryFn, 'target':targetSummaryFn, 'link':circosLinkFn, 'anno':circosAnnoFn}
    for f in outFileDict :
        outFileDict[f] = open(outFileDict[f], 'w')

    report(resD, outFileDict, filterList)
    report_df(resD, outFileDict['df'], filterList)

    for f in outFileDict :
        outFileDict[f].close()
