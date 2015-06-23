#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import re
from math import log
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import patches
import breakmer.assembly.olc as olcAssembly

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


class Segment:
    def __init__(self, alignResult, segmentColor):
        """ """
        self.alignResult = alignResult
        self.color = segmentColor
        self.queryCoordinates = [alignResult.qstart(), alignResult.qend()]
        self.indelCoordinates = alignResult.breakpts.contigBreakpoints
        self.strand = alignResult.strand
        self.alignLen = alignResult.get_query_span()
        self.genomicCoords = alignResult.alignVals.get_coords('ref')


class AlignSegments:
    def __init__(self, svEventResult):
        """ """
        self.svEventResult = svEventResult
        self.segments = None
        self.colors = ['green', 'orange', 'blue', 'orange', 'purple']
        self.orderedSeqs = None
        self.setup()

    def setup(self):
        """ """
        for i, blatResult in enumerate(self.svEventResult.blatResults):
            segment = self.Segment(blatResult, self.colors[i])

    def get_contig_seq(self):
        """ """
        return self.svEventResult.contig.seq

    def set_orderedseqs(self, orderedSeqs):
        """ """
        self.orderedSeqs = orderedSeqs

    def get_segment_color(self, nucIter):
        """ """
        colors = []
        for segment in self.segments:
            if nucIter >= segment.queryCoordinates[0] and nucIter <= segment.queryCoordinates[1]:
                colors.append(segment.color)
        if len(colors) > 1:
            colors = 'black'
        else:
            colors = colors[0]
        return colors


def generate_pileup_img(svEventResult, bamReadsFn, outPath, contigId):
    """ """
    segmentManager = AlignSegments(svEventResult)
    bamFile = pysam.Samfile(bamReadsFn, "rb")
    segmentManager.set_orderedseqs(pile_reads(bamFile.fetch(), svEventResult.contig.seq))
    plot_pileup(segmentManager, os.path.join(outPath, contigId))


def pile_reads(reads, contigSeq):
    """ """
    orderedSeqs = []
    for read in reads:
        idx = contigSeq.find(read.seq)
        seq = read.seq
        add = True
        if idx == -1:
            aln1 = olcAssembly.nw(contigSeq, read.seq)
            aln2 = olcAssembly.nw(read.seq, contigSeq)
            idx = aln1[3]
            seq = aln1[1]
            if aln1[-1] < aln2[-1]:
                idx = aln2[5]
                seq = aln2[0]
        if add:
            orderedSeqs.append((idx, ' ' * idx + seq))
    os = sorted(orderedSeqs, key=lambda x: x[0])
    return os


def plot_pileup(segmentManager, outBaseFn):
    """ """
    # Determine coordinate constants
    seqPlotSize = (len(segmentManager.orderedSeqs) + 1) * 0.75
    plotHeight = 5 # seqPlotSize*1.5
    if len(orderedSeqs) > 10:
        plotHeight = 15

    # Setup figure
    fig = plt.figure(figsize=(30, plotHeight), frameon=False)
    ax = fig.add_subplot(111)
    ax.axis('off')

    # Set the y-index for the sequence plotting
    yCoord = 0
    # Start plotting at unit 2 on the x-axis
    xOffset = 2
    # Increment text by 1 unit
    xInc = 1
    plot_realignment_strands(ax, yCoord + 0.5, xOffset, segmentManager)
    plot_contig_seq(ax, yCoord, xOffset, segmentManager)
    plot_pileup_seq(ax, yCoord, xOffset, segmentManager)
    plot_segments(ax, yCoord + 1, xOffset, segmentManager)

#     annoYidx = seqYidx + len(cSeq.segments) + 1
#     # Vertical breakpojnt lines, colors match the segments.
#     brkptLines = []
#     contigGenomicCoords = []
#     ycoord = -seqPlotSize-1
#     for i, seg in enumerate(cSeq.segments):
#         coords = seg.coords
#         gCoords = seg.genomicBrkpts
#         yoffset = float(i)/float(5)
#         seg.yidx = annoYidx - yoffset
#         plt.plot((xoffset+coords[0], xoffset+coords[0]+seg.plotLen), (seg.yidx, seg.yidx), color=seg.color, linewidth=2)
#         lr = ['left', 'right']
#         for coord, gCoord, leftRight in zip(coords, gCoords, lr):
#             print cSeq.contigBrkpts
#             if coord in cSeq.contigBrkpts:
#                 brkptLines.append((coord, seg.color, seg.yidx))
#             contigGenomicCoords.append((coord, gCoord, leftRight, ycoord))
#         ycoord -= 0.8

#     brkptIter = 0
#     bLines = sorted(brkptLines, key=lambda x: x[0])
#     print 'Blines', bLines
#     bIdxs = []
#     for coord, scolor, syidx in bLines:
#         xidx = coord + xoffset + brkptIter
#         if coord in bIdxs :
#             xidx = coord + xoffset
#         bIdxs.append(coord)
#         plt.plot((xidx, xidx), (-seqPlotSize-3, syidx), color=scolor, linewidth=2)
#         brkptIter += 1

#     for coord, gCoord, leftRight, coordYidx in contigGenomicCoords:
#         xidx = coord + xoffset
#         if leftRight == 'right' :
#             xidx += 0.5
#         ax.text(xidx, coordYidx, gCoord, ha=leftRight, va='bottom', size=12)

#     pa = plotAnnot(svRes, annot, cSeq)
#     pa.get_coding_features()

#     exonXoffset = xoffset
#     iter = 0
#     for segAnnot in pa.segAnnots:
#         for exCoords in segAnnot.mappedExons:
#             ycoord = annoYidx - (float(iter)/float(5))
#             rect = patches.Rectangle((segAnnot.segment.coords[0]+exonXoffset+exCoords[0], ycoord), exCoords[1]-exCoords[0], 1, color=segAnnot.segment.color)
#             ax.add_patch(rect)
#             ax.text(segAnnot.segment.coords[0]+exonXoffset+exCoords[0], annoYidx+1, exCoords[2], ha='left', va='bottom', size=7, rotation=60)
#         if segAnnot.geneName:
#             offset = segAnnot.segment.coords[0]+exonXoffset
#             midDist = float(abs(segAnnot.segment.coords[1] - segAnnot.segment.coords[0])) / float(2)
#             ax.text(offset + midDist, annoYidx+3, segAnnot.geneName +'('+ segAnnot.geneStrand+')', ha='center', va='bottom', size=12, style='italic')
# #        exonXoffset += segAnnot.segment.plotLen
#         if segAnnot.genomicEnd:
#             offsets = [exonXoffset+segAnnot.segment.coords[0], exonXoffset+segAnnot.segment.coords[1]]
#             coordTexts = [segAnnot.genomicEnd, segAnnot.bp]
#             leftRight = ['left', 'right']
#             if segAnnot.pos == 'last':
#                 coordTexts.reverse()
#             for offset, coordText, lr in zip(offsets, coordTexts, leftRight):
#                 if lr == 'right':
#                     offset += 0.5
#                 ax.text(offset, annoYidx-(iter*0.8)-0.5, coordText, ha=lr, va='top', size=12)
#         iter += 1

    ySize = (len(segmentManager.orderedSeqs) + 1) * 0.75 + 1.5 + 10
    ax.axis([0, len(svRes.contig) + 10, -seqPlotSize - 5, 8])
    plt.savefig(outName + '.pdf', bbox_inches='tight', dpi=300)
    plt.savefig(outName + '.png', bbox_inches='tight', dpi=300)


def plot_realignment_strands(ax, seqYidx, xOffset, segmentManager):
    """Plot the strand the segments were realigned to the reference - +/-"""
    for segment in segmentManager.segments:
        queryStartCoord, queryEndCoord = segment.queryCoordinates
        midDist = float(abs(queryEndCoord - queryStartCoord)) / float(2)
        xCoord = xOffset + queryStartCoord + midDist
        ax.text(xCoord, seqYidx + 0.5, segment.strand, ha='center', va='bottom', size=14, family='monospace')


def add_seq_text(ax, x, y, char, color='black'):
    """ """
    ax.text(x, y, char, ha='center', va='center', size=11, family='monospace', color=color)


def plot_contig_seq(ax, seqYidx, xOffset, segmentManager):
    """ """
    xInc = 1
    # Add 5' designation to contig sequence.
    add_seq_text(ax, 1, seqYidx, "5'")
    # Iterate over the nucleotides of the contig sequence.
    for nucIter, nuc in enumerate(segmentManager.get_contig_seq()):
        add_seq_text(ax, xOffset, seqYidx, nuc)
        xOffset += xInc
        # Insert a pipe character for the breakpoint in the contig seq.
        # if nucIter in cSeq.contigBrkpts:
        #     add_seq_text(ax, xOffset, seqYidx, ' ')
        #     xOffset += xinc
        nucIter += 1
    # Add 3' designation to contig sequence.
    add_seq_text(ax, xOffset, seqYidx, "3'")


def plot_segments(ax, yCoord, xOffset, segmentManager):
    """ """
    for segment in segmentManager.segments:
        # Plot rectangles for each realignment result
        rect = patches.Rectangle((segAnnot.segment.coords[0]+exonXoffset+exCoords[0], ycoord), exCoords[1]-exCoords[0], 1, color=segAnnot.segment.color)
        ax.add_patch(rect)


def plot_pileup_seq(ax, seqYidx, xOffset, segmentManager):
    """ """
    yInc = 0.75
    xInc = 1
    # Iterate through sequences.
    for idx, seq in segmentManager.orderedSeqs:
        seqTextOff = xOffset
        seqYidx = seqYidx - yInc
        segIdx = 0
        brkIdx = 0
        for nuc in seq:
            nucColor = segmentManager.get_segment_color(nucIter)
            add_seq_text(ax, seqTextOff, seqYidx, nuc, nucColor)
            seqTextOff += xInc


# class annotation():
#     def __init__(self, annoFn):
#         self.trxs = {}
#         self.setup_anno(annoFn)

#     def setup_anno(self, annoFn):
#         for line in open(annoFn, 'r'):
#             if line.find('#') > -1:
#                 continue
#             line = line.strip()
#             chr, src, type, start, stop, fill, strand, fill2, meta = line.split('\t')
#             meta = meta.split(';')
#             trxId = meta[1].split(' ')[2].lstrip('"').rstrip('"')
#             gId = meta[0].split(' ')[1].lstrip('"').rstrip('"')
#             gName = meta[4].split(' ')[2].lstrip('"').rstrip('"')
#             gStatus = meta[3].split(' ')[2].lstrip('"').rstrip('"')
#             if gStatus != "KNOWN":
#                 continue
#             if chr not in self.trxs:
#                 self.trxs[chr] = {}
#             if gId not in self.trxs[chr]:
#                 self.trxs[chr][gId] = {'trx':{}, 'canon':'', 'trxSizes':[], 'name':gName, 'strand':strand}
#             if trxId not in self.trxs[chr][gId]['trx']:
#                 self.trxs[chr][gId]['trx'][trxId] = {'bounds':[int(start), int(stop)], 'features':{}, 'strand':strand}

#             if type == 'transcript':
#                 self.trxs[chr][gId]['trx'][trxId]['bounds'] = [int(start), int(stop)]
#                 self.trxs[chr][gId]['trxSizes'].append((abs(int(stop)-int(start)), trxId))
#             elif type != 'gene' and type != 'transcript' and type != 'CDS':
#                 if type not in self.trxs[chr][gId]['trx'][trxId]['features']:
#                     self.trxs[chr][gId]['trx'][trxId]['features'][type] = []
#                 self.trxs[chr][gId]['trx'][trxId]['features'][type].append((int(start), int(stop)))

#     def get_trx_exons(self, bp, strand, segLoc):
#         chr, bp = bp.split(":")
#         if bp.find('-') > -1:
#             return None
#         nearestGene = [None, None]
#         for achr in self.trxs:
#             if achr == chr:
#                 for gene in self.trxs[achr]:
#                     retVals = {}
#                     trxSizes = sorted(self.trxs[achr][gene]['trxSizes'], key=lambda x: x[0])
#                     canonTrxId = trxSizes[-1][1]
#                     trxStart = self.trxs[achr][gene]['trx'][canonTrxId]['bounds'][0]
#                     trxStop = self.trxs[achr][gene]['trx'][canonTrxId]['bounds'][1]
#                     #print bp, trxStart, trxStop
#                     #print trx, trxStart, trxStop, self.trxs[achr][trx]['strand']
#                     retVals['features'] = self.trxs[achr][gene]['trx'][canonTrxId]['features']
#                     retVals['name'] = self.trxs[achr][gene]['name']
#                     retVals['strand'] = self.trxs[achr][gene]['strand']
#                     retVals['trxId'] = canonTrxId
#                     if int(bp) >= int(trxStart) and int(bp) <= int(trxStop):
#                         return retVals
# #                        return self.trxs[achr][gene]['trx'][canonTrxId]['features'], self.trxs[achr][gene]['name'], self.trxs[achr][gene]['strand'], canonTrxId
#                     else:
#                         if segLoc == "first":
#                             if strand == '-':
#                                 # Get downstream gene
#                                 dist = int(trxStart) - int(bp)
#                                 if dist > 0 and not nearestGene[0]:
#                                     nearestGene = [dist, retVals]
# #                                    print 'first', nearestGene, dist, bp, trxStart
#                                 elif dist > 0 and dist < nearestGene[0]:
#                                     nearestGene = [dist, retVals]
# #                                    print 'change', nearestGene, dist, bp, trxStop
#                             elif strand == '+':
#                                 # Get upstream genes
#                                 dist = int(bp) - int(trxStop)
#                                 if dist > 0 and not nearestGene[0]:
#                                     nearestGene = [dist, retVals]
# #                                    print 'first', nearestGene, dist, bp, trxStop
#                                 elif dist > 0 and dist < nearestGene[0]:
#                                     nearestGene = [dist, retVals]
# #                                    print 'change', nearestGene, dist, bp, trxStop
#                         elif segLoc == "last":
#                             if strand == '-':
#                                 dist = int(bp) - int(trxStop)
#                                 if dist > 0 and not nearestGene[0]:
#                                     nearestGene = [dist, retVals]
#                                 elif dist > 0 and dist < nearestGene[0]:
#                                     nearestGene = [dist, retVals]
#                             elif strand == '+':
#                                 dist = int(trxStart) - int(bp)
#                                 if dist > 0 and not nearestGene[0]:
#                                     nearestGene = [dist, retVals]
#                                 elif dist > 0 and dist < nearestGene[0]:
#                                     nearestGene = [dist, retVals]
# #        print bp
# #        print nearestGene
# #        sys.exit()
#         return nearestGene[1]

# class segAnnot():
#     def __init__(self, annoRes, bp,  realignStrand, pos, segLen, segment):
#         self.exons = None
#         self.geneName = None
#         self.trxId = None
#         self.geneStrand = None
#         self.bp = bp
#         self.realignStrand = realignStrand
#         self.pos = pos
#         self.segLen = segLen
#         self.mappedExons = []
#         self.genomicEnd = None
#         self.segment = segment
#         self.features = []
#         self.setup(annoRes)

#     def setup(self, annoRes):
#         if annoRes:
#             self.features = annoRes['features']
#             self.geneName = annoRes['name']
#             self.trxId = annoRes['trxId']
#             self.geneStrand = annoRes['strand']

#     def get_brkpt(self):
#         return self.bp.split(':')[1]

#     def get_exons(self):
#         features = []
#         bp = self.get_brkpt()
#         for feature in self.features:
#             print 'Feature', feature
#             if feature != 'exon' :
#                 continue
#             fElements = self.features[feature]
#             reverse = False
#             if self.realignStrand == '-':
#                 reverse = True
#             fElements_sorted = sorted(fElements, key=lambda x: x[0], reverse=reverse)
#             eIter = 1
#             for element in fElements:
#                 add = False
#                 if self.pos == "first":
#                     if int(element[1]) >= int(bp) and self.realignStrand=='-':
#                         add = True
#                     elif int(element[0]) <= int(bp) and self.realignStrand=='+':
#                         add = True
#                     if add:
#                         features.append((element[0], element[1], feature+str(eIter)))
#                 elif self.pos == "last":
#                     if int(element[0]) <= int(bp) and self.realignStrand=='-':
#                         add = True
#                     elif int(element[1]) >= int(bp) and self.realignStrand=='+':
#                         add = True
#                     if add:
#                         features.append((element[0], element[1], feature+str(eIter)))
#                 eIter += 1
#         print 'Brkpt', bp
#         print 'Features', features
#         return features

#     def map_coding_features(self, exons):
#         bp = self.get_brkpt()
#         reverse = False
#         print self.realignStrand
#         if self.realignStrand == '-':
#             reverse = True
#         l = sorted(exons, key=lambda x: x[1], reverse=reverse)
#         genomicStart = self.get_brkpt()

#         print self.pos
#         print l
#         if self.pos == "first":
#             end = l[0][0]
#             genomicStart = end
#         elif self.pos == "last":
#             end = l[-1][1]

#         print exons
#         print end, bp
#         genomicLen = log(abs(end - int(bp)), 2)
#         bpUnits = float(self.segment.plotLen)/float(genomicLen)
#         print self.segment.plotLen
#         print genomicLen
#         print bpUnits
#         pos = []

#         previousCoord = 0
#         if self.pos == "last":
#             previousCoord = log(int(bp), 2)

#         intronLen = 0
#         for coords in l:
#             if coords[2].find('exon') == -1:
#                 continue
#             ll = [log(coords[0], 2), log(coords[1], 2)]
#             e1 = log(max(abs(int(genomicStart) - int(coords[0])), 1), 2)*bpUnits
#             e2 = log(max(abs(int(genomicStart) - int(coords[1])), 1), 2)*bpUnits
#             eCoords = [e1, e2]
#             eCoords.sort()
#             eCoords.append(coords[2])
#             pos.append(eCoords)
#         self.mappedExons = pos
#         self.genomicEnd = self.bp.split(':')[0] + ':' + str(end)

# class plotAnnot():
#     def __init__(self, res, annot, contigSeq):
#         self.segAnnots = []
#         self.res = res
#         self.run(res, annot, contigSeq)

#     def run(self, res, annot, contigSeq):
#         iter = 0
#         for strand, brkpt, segLen, segment in zip(res.strands, res.brkpts, res.segLens, contigSeq.segments):
#             if iter == 0:
#                 pos = "first"
#             elif iter == (len(res.strands)-1):
#                 pos = "last"
#             r = annot.get_trx_exons(brkpt, strand, pos)
#             self.segAnnots.append(segAnnot(r, brkpt, strand, pos, segLen, segment))
#             iter += 1

#     def get_coding_features(self):
#         for iter, segAnnot in enumerate(self.segAnnots):
#             e = segAnnot.get_exons()
#             segAnnot.map_coding_features(e)

# class svResult():
#     def __init__(self, fn, contigId):
#         self.contig = None
#         self.segMetrics = []
#         self.brkpts = []
#         self.strands = []
#         self.segLens = []
#         self.setup_results(fn, contigId)

#     def setup_results(self, fn, contigId):
#         i=0
#         for line in open(fn, 'r'):
#             if i==0:
#                 i+=1
#                 continue
#             if line.find(contigId) == -1:
#                 continue
#             line = line.strip()
#             linesplit = line.split('\t')
#             print 'Result values', linesplit
#             self.contig = linesplit[-1]
#             self.segMetrics = linesplit[5].split(',')
#             self.brkpts = linesplit[1].split(',')
#             self.cigar = linesplit[2].split(',')
#             self.strands = linesplit[4].split(',')
#         self.set_segment_lengths()

#     def set_segment_lengths(self):
#         for item in self.segMetrics:
#             self.segLens.append(int(item.split(":")[1]))

#     def get_segment_lengths(self):
#         lens = []
#         for item in self.segMetrics:
#             lens.append(int(item.split(":")[1]))
#         return lens

#     def get_contig_brkpts(self):
#         # Return the start and end positions relative to the contig sequence
#         # where the segmenets map to the reference genome.
#         # Note that there will usually be overlap between two segments.
#         l = len(self.contig)
#         contigBps = []
#         overlaps = []
#         overlapIndices = []
#         breaks = []
#         genomicBps = []
#         segLens = self.get_segment_lengths()
#         iter = 0
#         for seqCigar in self.cigar:
#             split = re.split('(\D+)', seqCigar)
#             matchLen = 0
#             for i, ssplit in enumerate(split):
#                 if ssplit == 'M' or ssplit == 'I':
#                     matchLen += int(split[i-1])
# #            matchLen = int(split[split.index('M')-1])
#             scLen = 0
#             if seqCigar.find('S') > -1 and split[1] != 'S':
#                 scLen = int(split[split.index('S')-1])
#             start = l - (matchLen + scLen)
#             end = start + matchLen
#             if len(contigBps) > 0:
#                 if contigBps[-1][1] < end:
#                     overlaps.append(end - contigBps[-1][1])
#                     overlapIndices.extend(range(contigBps[-1][1], end))
#             contigBps.append((start, end))
#             if iter > 0 :
#                 breaks.append(start)
#             elif iter != (len(self.cigar)-1):
#                 breaks.append(end)

#             gBps = []
#             if self.brkpts[iter].find('-') > -1:
#                 gBps = self.brkpts[iter].split('-')
#                 if self.strands[iter] == '-':
#                     gBps.reverse()
#             else:
#                 chr, bp = self.brkpts[iter].split(':')
#                 sl = segLens[iter]
#                 if self.strands[iter] == '-':
#                     if iter == 0:
#                         gBps = [chr+':'+str(int(bp)+int(sl)), chr+':'+bp]
#                     else:
#                         gBps = [chr+':'+bp, chr+':'+str(int(bp)-int(sl))]
#                 elif self.strands[iter] == '+':
#                     if iter == 0:
#                         gBps = [chr+':'+str(int(bp)-int(sl)), chr+':'+bp]
#                     else:
#                         gBps = [chr+':'+bp, chr+':'+str(int(bp)+int(sl))]
#             genomicBps.append(gBps)
#             iter += 1
#         breakIdx = list(set(breaks))
#         breakIdx.sort()
#         print 'Breaks', breakIdx
#         return contigBps, overlaps, overlapIndices, breakIdx, genomicBps

# class segment:
#     def __init__(self, len, coords, pos, color, plotLen, gBp, alignStrand):
#         self.len = len
#         self.coords = coords
#         self.pos = pos
#         self.color = color
#         self.plotLen = plotLen
#         self.yidx = 0
#         self.genomicBrkpts = gBp
#         self.alignStrand = alignStrand

# class contigSeq:
#     def __init__(self, svRes):
#         self.contigSeqLen = len(svRes.contig)
#         self.totalPlotLen = 0
#         self.segments = []
#         self.overlaps = []
#         self.overlapIndices = []
#         self.contigBrkpts = []
#         self.genomicBrkpts = []
#         self.setup(svRes)

#     def setup(self, svRes):
#         nSegs = len(svRes.segMetrics)
#         segLens = svRes.get_segment_lengths()
#         segStartEnd, self.overlaps, self.overlapIndices, self.contigBrkpts, self.genomicBrkpts = svRes.get_contig_brkpts()
#         self.totalPlotLen = self.contigSeqLen + len(self.contigBrkpts)

#         segColors = ['green', 'orange', 'blue']
#         for i in range(nSegs):
#             segPlotLen = segStartEnd[i][1] - segStartEnd[i][0] #segLens[i]
#             segAlignStrand = svRes.strands[i]
#             for brkpt in self.contigBrkpts:
#                 if brkpt < segStartEnd[i][1]:
#                     segPlotLen += 1
#             self.segments.append(segment(segLens[i], segStartEnd[i], i, segColors[i], segPlotLen, self.genomicBrkpts[i], segAlignStrand))

#     def get_segment(self, nucIter):
#         # Determine which segment the iterator is within.
#         segColor = ''
#         brkptColor = ''
#         segs = []
#         for segment in self.segments:
#             if nucIter <= segment.coords[1] and nucIter > segment.coords[0]:
#                 segs.append(segment)
#             if (nucIter == segment.coords[1]) or (nucIter == segment.coords[0]):
#                 brkptColor = segment.color
#                 if brkptColor != '':
#                     brkptColor = 'black'

#         if len(segs) == 0:
#             segColor = 'grey'
#         elif len(segs) > 1:
#             # Overlapping region
#             segColor = 'black'
#         else:
#             segColor = segs[0].color
#         return segColor, brkptColor, segs

# def map_coding_features(f, s, iter, bp, len):
#     reverse = False
#     if s == '-':
#         reverse = True
#     l = sorted(f, key=lambda x: x[1], reverse=reverse)
#     genomicStart = bp
#     if iter == 1:
#         end = l[0][1]
#         genomicStart = end
#     elif iter == 2:
#         end = l[-1][1]

#     print end, bp, iter
#     genomicLen = log(abs(end - int(bp)), 2)
#     print genomicLen
#     bpUnits = float(len)/float(genomicLen)
#     print bpUnits
#     intronUnits = bpUnits
#     exonUnits = bpUnits
#     pos = []

#     previousCoord = 0
#     if iter ==2:
#         previousCoord = log(int(bp), 2)

#     intronLen = 0
#     for coords in l:
#         if coords[2].find('exon') == -1:
#             continue
#         ll = [log(coords[0], 2), log(coords[1], 2)]
#         e1 = log(max(abs(int(genomicStart) - int(coords[0])), 1), 2)*bpUnits
#         e2 = log(max(abs(int(genomicStart) - int(coords[1])), 1), 2)*bpUnits
#         eCoords = [e1, e2]
#         eCoords.sort()
#         eCoords.append(coords[2])
#         pos.append(eCoords)
#     return pos, end

# def get_coding_features(f, g, s, iter, bp):
#     features = []
#     for feature in f:
#         if feature != 'exon' :
#             continue

#         fElements = f[feature]
#         reverse = False
#         if s == '-':
#             reverse = True
#         fElements_sorted = sorted(fElements, key=lambda x: x[0], reverse=reverse)
#         eIter = 1
#         for element in fElements:
#             add = False
#             if iter == 1 :
#                 if int(element[1]) >= int(bp) and s=='-':
#                     add = True
#                 elif int(element[0]) <= int(bp) and s=='+':
#                     add = True
#                 if add:
#                     features.append((element[0], element[1], feature+str(eIter)))
#             elif iter == 2:
#                 if int(element[0]) <= int(bp) and s=='-':
#                     add = True
#                 elif int(element[1]) >= int(bp) and s=='+':
#                     add = True
#                 if add:
#                     features.append((element[0], element[1], feature+str(eIter)))
#             eIter += 1
#     return features
