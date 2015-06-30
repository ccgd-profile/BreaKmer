#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import re
import pysam
from math import log
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import patches
import breakmer.assembly.olc as olcAssembly

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


class TrxBrkpt:
    def __init__(self, distToTrx, svBrkpt, brkptIdx, svType):
        self.dist = distToTrx
        self.svBrkpt = svBrkpt
        self.brkptIdx = brkptIdx
        self.svType = svType

    def get_genomic_coord(self):
        return int(self.svBrkpt.genomicCoords[self.brkptIdx])


class AnnoTrx:
    def __init__(self, trx, trxDist, svBreakpoint, brkptIdx, svType):
        self.trx = trx
        self.brkpts = [TrxBrkpt(trxDist, svBreakpoint, brkptIdx, svType)]

    def add_brkpt(self, trxDist, svBreakpoint, brkptIdx, svType):
        self.brkpts.append(TrxBrkpt(trxDist, svBreakpoint, brkptIdx, svType))


def check_add_trx(trx, trxItems, trxIds, trxDist, svBreakpoint, brkptIdx, svType):
    if trx.id in trxIds:
        idx = trxIds.index(trx.id)
        trxItems[idx].add_brkpt(trxDist, svBreakpoint, brkptIdx, svType)
    else:
        trxItems.append(AnnoTrx(trx, trxDist, svBreakpoint, brkptIdx, svType))
        trxIds.append(trx.id)
    return trxItems, trxIds


class Segment:
    def __init__(self, alignResult, segmentColor, segmentIdx, nSegments):
        """ """
        self.alignResult = alignResult
        self.color = segmentColor
        self.queryCoordinates = [alignResult.qstart(), alignResult.qend()]
        self.genomicCoordinates = [alignResult.tstart(), alignResult.tend()]
        self.chromName = alignResult.get_seq_name('ref')
        self.indelCoordinates = alignResult.breakpts.contigBreakpoints
        self.indelSizes = alignResult.indel_sizes
        self.strand = alignResult.strand
        self.alignLen = alignResult.get_query_span()
        self.idx = segmentIdx
        self.nSegments = nSegments
        self.genomicCoords = alignResult.alignVals.get_coords('ref')

    def get_len(self):
        """ """
        return self.queryCoordinates[1] - self.queryCoordinates[0]

    def get_segment_trxs(self):
        svBreakpoints = self.alignResult.get_sv_brkpts()
        # Determine the number of transcripts for this segment based on the sv breakpoints
        trxItems = []
        trxIds = []
        for svBreakpoint in svBreakpoints:
            annotatedTrxsDict = svBreakpoint.annotated_trxs
            dKeys = annotatedTrxsDict.keys()
            dKeys.sort()
            # If svBreakpoint type is 'indel' then there can be two trxs associated with the svBreakpoint
            # If the type is 'rearrangement' then there should only be one - this needs to be inferred by the realignment
            # strand and the index of the segment.
            brkptTrxs = []

            if svBreakpoint.svType == 'rearrangement':
                if len(dKeys) > 1:
                    # This indicates that the segment is in the middle - key innner transcripts if the breakpoints are intergenic.
                    leftBpTrxList, leftBpDistList = annotatedTrxsDict[0]
                    rightBpTrxList, rightBpDistList = annotatedTrxsDict[1]
                    keepIdx = 0
                    if len(leftBpTrxList) > 1:
                        # Take the inner trxs
                        keepIdx = 1
                    # Left
                    trxItems, trxIds = check_add_trx(leftBpTrxList[keepIdx], trxItems, trxIds, leftBpDistList[keepIdx], svBreakpoint, 0, 'rearr')
                    keepIdx = 0
                    if len(rightBpTrxList) > 1:
                        # Take the inner trxs
                        keepIdx = 0
                    # Right
                    trxItems, trxIds = check_add_trx(rightBpTrxList[keepIdx], trxItems, trxIds, rightBpDistList[keepIdx], svBreakpoint, 1, 'rearr')
                else:
                    # Single breakpoint
                    trxList, distList = annotatedTrxsDict[0]
                    if len(trxList) > 1:
                        # Pick which transcript to keep based on strands, breakpoint is outside of a transcript
                        if self.idx == 0:
                            # First
                            if self.alignResult.strand == '-':
                                # Get downstream gene
                                trx = trxList[1]
                            else:
                                trx = trxList[0]
                        elif self.idx == (self.nSegments - 1):
                            if self.alignResult.strand == '-':
                                # Get upstream gene
                                trx = trxList[0]
                            else:
                                trx = trxList[1]
                        trxItems, trxIds = check_add_trx(trx, trxItems, trxIds, trxDist, svBreakpoint, 0, 'rearr')
                    else:
                        # lands in a single trancript
                        trxItems, trxIds = check_add_trx(trxList[0], trxItems, trxIds, distList[0], svBreakpoint, 0, 'rearr')
            elif svBreakpoint.svType == 'indel':
                if len(dKeys) == 1:
                    # Insertion with one genomic breakpoint
                    trxList, distList = annotatedTrxsDict[0]
                    trxItems, trxIds = check_add_trx(trxList[0], trxItems, trxIds, distList[0], svBreakpoint, 0, 'ins')
                else:
                    # Deletion with two genomic breakpoints, if intergenic then keep the outer transcripts
                    leftBpTrxList, leftBpDistList = annotatedTrxsDict[0]
                    rightBpTrxList, rightBpDistList = annotatedTrxsDict[1]
                    # Take the first trx no matter what
                    trx = leftBpTrxList[0]
                    trxDist = leftBpDistList[0]
                    trxItems, trxIds = check_add_trx(trx, trxItems, trxIds, trxDist, svBreakpoint, 0, 'del')
                    keepIdx = 0
                    if len(rightBpTrxList) > 1:
                        # Take the outer trx
                        keepIdx = 1
                    trx = rightBpTrxList[keepIdx]
                    trxDist = rightBpDistList[keepIdx]
                    trxItems, trxIds = check_add_trx(trx, trxItems, trxIds, trxDist, svBreakpoint, 1, 'del')
        return trxItems, trxIds


class AlignSegments:
    def __init__(self, svEventResult):
        """ """
        self.svEventResult = svEventResult
        self.segments = []
        self.colors = ['green', 'orange', 'blue', 'orange', 'purple']
        self.orderedSeqs = None
        self.setup()

    def has_annotations(self):
        """ """
        return self.svEventResult.annotated

    def setup(self):
        """ """
        for i, blatResult in enumerate(self.svEventResult.blatResults):
            self.segments.append(Segment(blatResult[1], self.colors[i], i, len(self.svEventResult.blatResults)))

    def get_contig_seq(self):
        """ """
        return self.svEventResult.contig.seq

    def get_contig_id(self):
        """ """
        return self.svEventResult.contig.get_id()

    def set_orderedseqs(self, orderedSeqs):
        """ """
        self.orderedSeqs = orderedSeqs

    def get_segment_color(self, nucIter):
        """ """
        colors = ['grey']
        for segment in self.segments:
            # print nucIter
            # print segment.queryCoordinates
            if (nucIter >= segment.queryCoordinates[0]) and (nucIter < segment.queryCoordinates[1]):
                colors.append(segment.color)
                # print colors
        if len(colors) > 2:
            colors = 'black'
        else:
            colors = colors[1]
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


def plot_pileup(segmentManager, outName):
    """ """
    # Determine coordinate constants
    seqPlotSize = (len(segmentManager.orderedSeqs) + 1) * 0.75
    plotHeight = round(seqPlotSize) + 5 # seqPlotSize*1.5
    # if len(segmentManager.orderedSeqs) > 10:
    #     plotHeight = 20

    # Setup figure
    fig = plt.figure(figsize=(25, plotHeight), frameon=False)
    ax = fig.add_subplot(111)
    ax.axis('off')

    # Set the y-index for the sequence plotting
    yCoord = 0
    # Start plotting at unit 2 on the x-axis
    xOffset = 2
    # Increment text by 1 unit
    xInc = 1
    # plot_realignment_strands(ax, yCoord + 0.5, xOffset, segmentManager)
    plot_contig_seq(ax, yCoord, xOffset, segmentManager)
    plot_pileup_seq(ax, yCoord, xOffset, segmentManager)
    plot_segments(ax, yCoord + 1, xOffset, segmentManager)
    plot_indel_track(ax, yCoord + 1, xOffset, segmentManager)
    plot_annotation_track(ax, yCoord + 5, xOffset, segmentManager)
    plot_global_trx_track(ax, yCoord + 8, xOffset, segmentManager)
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
    ax.axis([0, len(segmentManager.get_contig_seq()) + 8, -seqPlotSize - 1, 10])
    plt.savefig(outName + '.pdf', bbox_inches='tight', dpi=300)
    plt.savefig(outName + '.png', bbox_inches='tight', dpi=300)
    plt.savefig(outName + '.svg')


# def plot_realignment_strands(ax, seqYidx, xOffset, segmentManager):
#     """Plot the strand the segments were realigned to the reference - +/-"""
#     for segment in segmentManager.segments:
#         queryStartCoord, queryEndCoord = segment.queryCoordinates
#         midDist = float(abs(queryEndCoord - queryStartCoord)) / float(2)
#         xCoord = xOffset + queryStartCoord + midDist
#         ax.text(xCoord, seqYidx + 0.5, segment.strand, ha='center', va='bottom', size=14, family='monospace')


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
    for i, segment in enumerate(segmentManager.segments):
        # Plot rectangles for each realignment result
        xCoord = xOffset + segment.queryCoordinates[0]
        yCoord = yCoord + ((i + 0.75) * 0.50)
        rectLen = segment.queryCoordinates[1] - segment.queryCoordinates[0]
        rectHeight = 0.25
        lenText = str(rectLen) + 'bp'
        xCoordLabel = xCoord + (float(rectLen) / float(2))
        rect = patches.Rectangle((xCoord, yCoord), rectLen, rectHeight, color=segment.color)
        ax.add_patch(rect)
        ax.text(xCoordLabel, yCoord - 0.15, lenText + ' (' + segment.strand + ')', ha='center', va='top', size=10)
        # Plot genomic coordinates of the segment
        gCoordOrder = [xCoord, xCoord + rectLen]
        if segment.strand == '-':
            gCoordOrder = [xCoord + rectLen, xCoord]
        horizAlign = ['left', 'right']
        if segment.strand == '-':
            horizAlign.reverse()
        segCoordStart = segment.chromName + ':' + str(segment.genomicCoordinates[0])
        segCoordEnd = segment.chromName + ':' + str(segment.genomicCoordinates[1])
        ax.text(gCoordOrder[0], yCoord - 0.15, segCoordStart, ha=horizAlign[0], va='top', size=10)
        ax.text(gCoordOrder[1], yCoord - 0.15, segCoordEnd, ha=horizAlign[1], va='top', size=10)


def plot_indel_track(ax, yCoord, xOffset, segmentManager):
    """ """
    for i, segment in enumerate(segmentManager.segments):
        indelCoordinates = segment.indelCoordinates
        for coord in indelCoordinates:
            xCoord = xOffset + coord[0]
            yCoord = yCoord + ((i + 1) * 0.5)
            rectLen = 1
            indelType = 'D'
            if len(coord) == 2:
                rectLen = coord[1] - coord[0]
                indelType = 'I'
            rectHeight = 0.5
            rect = patches.Rectangle((xCoord, yCoord), rectLen, rectHeight, color='red')
            ax.add_patch(rect)
            xCoordLabel = xCoord + (float(rectLen) / float(2))
            ax.text(xCoordLabel, yCoord + 1, self.indelSizes[i], ha='center', va='top', size=10)


def plot_pileup_seq(ax, seqYidx, xOffset, segmentManager):
    """ """
    yInc = 0.75
    xInc = 1
    # Iterate through sequences.
    print segmentManager.get_contig_seq(), len(segmentManager.get_contig_seq())
    for idx, seq in segmentManager.orderedSeqs:
        seqTextOff = xOffset
        seqYidx = seqYidx - yInc
        segIdx = 0
        brkIdx = 0
        nucIter = 0
        for nuc in seq:
            nucColor = segmentManager.get_segment_color(nucIter)
            add_seq_text(ax, seqTextOff, seqYidx, nuc, nucColor)
            seqTextOff += xInc
            nucIter += 1
        print seq, idx, nucIter


class AnnotationBrkpt:
    def __init__(self, trxBrkpts, segPos, segStrand):
        self.segPos = segPos
        self.segStrand = segStrand
        self.trxBrkpts = trxBrkpts
        self.other_brkpts = None
        self.bps = []
        self.bounds = []
        self.setup()

    def setup(self):
        """ """
        for bp in self.trxBrkpts:
            exonCode = 'right'
            if self.segPos == 'only':
                if bp.svType == 'del':
                    if bp.brkptIdx == 0:
                        exonCode = 'left'
                elif bp.svType == 'ins':
                    exonCode = 'all'
            elif self.segPos == 'first':
                if bp.svType == 'rearr':
                    if self.segStrand == '+':
                        exonCode = 'left'
            elif self.segPos == 'middle':
                if bp.svType == 'rearr':
                    self.bounds.append(bp.get_genomic_coord())
                    if bp.brkptIdx == 1:
                        if self.segStrand == '+':
                            exonCode = 'left'
                    elif bp.brkptIdx == 0:
                        if self.segStrand == '-':
                            exonCode = 'left'
            elif self.segPos == 'last':
                if bp.svType == 'rearr':
                    if self.segStrand == '-':
                        exonCode = 'left'
            print 'sv_viz.py setup() adding breakpoint', (bp, bp.get_genomic_coord(), exonCode)
            self.bps.append((bp, bp.get_genomic_coord(), exonCode))

    def add_brkpt(self, trxBrkpts):
        """ """
        self.other_brkpts = trxBrkpts

    def select_exons(self, exons):
        selectedExons = {}
        if len(self.bounds) > 0:
            print 'Bounds', self.bounds
            self.bounds.sort()
            bpCoordKey = '-'.join([str(x) for x in self.bounds])
            selectedExons[bpCoordKey] = {'coords': []}
            selectedExons[bpCoordKey]['coords'].append((self.bounds[0] - 1, self.bounds[0], 'breakpoint', None))
            selectedExons[bpCoordKey]['coords'].append((self.bounds[1] - 1, self.bounds[1], 'breakpoint', None))
            eIter = 1
            maxminCoords = [self.bounds[0], self.bounds[1], self.bounds[0], self.bps[0][2]]
            bpOverlap = [False, None]
            for exon in exons:
                # print 'Check exon', exon.start, exon.stop, exon.featureType
                add = False
                estart = int(exon.start)
                estop = int(exon.stop)
                exonCoords = [estart, estop]
                if (estart >= self.bounds[0] and estart <= self.bounds[1]):
                    add = True
                    if estop > self.bounds[1]:
                        bpOverlap = [True, self.bounds[1] - 1]
                        exonCoords[1] = self.bounds[1]
                elif (estop >= self.bounds[0] and estop <= self.bounds[1]):
                    add = True
                    if estart < self.bounds[0]:
                        bpOverlap = [True, self.bounds[0] - 1]
                        exonCoords[0] = self.bounds[0]
                if add:
                    # print 'sv_viz.py keep exon', bp, estart, estop, exonCode, exon.featureType
                    # absDist = abs(bpCoord - int(exonCoords[0]))
                    # if len(firstLastExons['nearest_exon']) == 0:
                    #     firstLastExons['nearest_exon'] = [absDist, len(selectedExons), 'exon' + str(eIter)]
                    # elif absDist < firstLastExons['nearest_exon'][0]:
                    #     firstLastExons['nearest_exon'] = [absDist, len(selectedExons), 'exon' + str(eIter)]
                    # if len(firstLastExons['furthest_exon']) == 0:
                    #     firstLastExons['furthest_exon'] = [absDist, len(selectedExons), 'exon' + str(eIter)]
                    # elif absDist > firstLastExons['furthest_exon'][0]:
                    #     firstLastExons['furthest_exon'] = [absDist, len(selectedExons), 'exon' + str(eIter)]

                    selectedExons[bpCoordKey]['coords'].append([int(exonCoords[0]), int(exonCoords[1]), 'exon' + str(eIter)], bpOverlap[1])
                    if maxminCoords[0] > int(exonCoords[0]):
                        maxminCoords[0] = int(exonCoords[0])
                    if maxminCoords[1] < int(exonCoords[1]):
                        maxminCoords[1] = int(exonCoords[1])
                eIter += 1
            selectedExons[bpCoordKey]['maxmincoords'] = maxminCoords
        else:
            for bp in self.bps:
                maxminCoords = []
                bpObj, bpCoord, exonCode = bp
                selectedExons[bpCoord] = {'coords': []}
                selectedExons[bpCoord]['coords'].append((bpCoord - 1, bpCoord, 'breakpoint', None))
                if len(maxminCoords) == 0:
                    maxminCoords = [bpCoord - 1, bpCoord, bpCoord, exonCode]
                eIter = 1
                firstLastExons = {'nearest_exon': [], 'furthest_exon': []}
                bpOverlap = [False, None]
                for exon in exons:
                    # print 'Check exon', exon.start, exon.stop, exon.featureType
                    add = False
                    estart = int(exon.start)
                    estop = int(exon.stop)
                    exonCoords = [estart, estop]
                    if (exonCode == 'left') and (estart <= bpCoord):
                        # Get all exons with start < bp
                        if bpCoord < estop:
                            # Breakpoint intersects with exon, reduce feature count to 2
                            bpOverlap = [True, bpCoord - 1]
                            exonCoords[1] = bpCoord
                        add = True
                    elif (exonCode == 'right') and (estop >= bpCoord):
                        if bpCoord > estart:
                            bpOverlap = [True, bpCoord - 1]
                            exonCoords[0] = bpCoord
                        add = True
                    elif exonCode == 'all':
                        # Single insertion in a gene
                        if bpCoord >= estart and bpCoord <= estop:
                            bpOverlap = [True, bpCoord - 1]
                        add = True
                    if add:
                        # print 'sv_viz.py keep exon', bp, estart, estop, exonCode, exon.featureType
                        # absDist = abs(bpCoord - int(exonCoords[0]))
                        # if len(firstLastExons['nearest_exon']) == 0:
                        #     firstLastExons['nearest_exon'] = [absDist, len(selectedExons), 'exon' + str(eIter)]
                        # elif absDist < firstLastExons['nearest_exon'][0]:
                        #     firstLastExons['nearest_exon'] = [absDist, len(selectedExons), 'exon' + str(eIter)]
                        # if len(firstLastExons['furthest_exon']) == 0:
                        #     firstLastExons['furthest_exon'] = [absDist, len(selectedExons), 'exon' + str(eIter)]
                        # elif absDist > firstLastExons['furthest_exon'][0]:
                        #     firstLastExons['furthest_exon'] = [absDist, len(selectedExons), 'exon' + str(eIter)]
                        selectedExons[bpCoord]['coords'].append([int(exonCoords[0]), int(exonCoords[1]), 'exon' + str(eIter), bpOverlap[1]])
                        if maxminCoords[0] > int(exonCoords[0]):
                            maxminCoords[0] = int(exonCoords[0])
                        if maxminCoords[1] < int(exonCoords[1]):
                            maxminCoords[1] = int(exonCoords[1])
                    eIter += 1
                # selectedExons[bpCoord]['coords'][firstLastExons['nearest_exon'][1]][2] = firstLastExons['nearest_exon'][2]
                # selectedExons[bpCoord]['coords'][firstLastExons['furthest_exon'][1]][2] = firstLastExons['furthest_exon'][2]
                selectedExons[bpCoord]['maxmincoords'] = maxminCoords
        return selectedExons


def determine_annotation_brkpts(trxBrkpts, segPos, segStrand):
    """ """
    abrkpt = None
    brkptTypes = {}
    for brkpt in trxBrkpts:
        if brkpt.svType not in brkptTypes:
            brkptTypes[brkpt.svType] = []
        brkptTypes[brkpt.svType].append(brkpt)

    if 'rearr' in brkptTypes:
        abrkpt = AnnotationBrkpt(brkptTypes['rearr'], segPos, segStrand)
        if 'del' in brkptTypes:
            abrkpt.add_brkpt(brkptTypes['del'])
        if 'ins' in brkptTypes:
            abrkpt.add_brkpt(brkptTypes['ins'])
    elif 'del' in brkptTypes:
        abrkpt = AnnotationBrkpt(brkptTypes['del'], segPos, segStrand)
    else:
        abrkpt = AnnotationBrkpt(brkptTypes['ins'], segPos, segStrand)
    return abrkpt


def get_neighbor_exons(exons):
    """ """
    leftExonBuffer = []
    rightExonBuffer = []
    bpExonBuffer = {}
    currentBp = None
    bpOverlaps = []
    for exon in exons:
        print exon
        start, end, name, bpOverlapCoord = exon
        if name == 'breakpoint':
            bpExonBuffer[start] = {'left': leftExonBuffer, 'right': rightExonBuffer, 'add_to_list': False}
            currentBp = start
            leftExonBuffer = []
        elif currentBp is None:
            leftExonBuffer.append(exon)
            if bpOverlapCoord is not None:
                bpOverlaps.append(bpOverlapCoord)
        else:
            bpExonBuffer[currentBp]['right'].append(exon)
            if bpOverlapCoord is not None:
                bpOverlaps.append(bpOverlapCoord)
    finalList = []
    for item in bpExonBuffer:
        left = bpExonBuffer[item]['left']
        right = bpExonBuffer[item]['right']
        if len(left) > 0:
            finalList.extend(left[len(left) - 2: len(left)])
        if item not in bpOverlaps:
            finalList.append((item, item + 1, 'breakpoint', None))
        if len(right) > 0:
            finalList.extend(right[0:2])
    return finalList


def plot_global_trx_track(ax, yCoord, xOffset, segmentManager):
    """ """
    if not segmentManager.has_annotations():
        return

    segStarts = []
    for i, segment in enumerate(segmentManager.segments):
        segStarts.append((segment.queryCoordinates[0], segment))

    sortedSegs = sorted(segStarts, key=lambda x: x[0])

    for i, segmentTuple in enumerate(sortedSegs):
        print 'sv_viz.py plot_annotation_track segment', i
        segment = segmentTuple[1]
        segmentPos = 'only'
        if len(sortedSegs) > 1:
            if i == 0:
                segmentPos = 'first'
            elif i > 0 and i < (len(sortedSegs) - 1):
                segmentPos = 'middle'
            elif i == (len(sortedSegs) - 1):
                segmentPos = 'last'
        print 'segment position', segmentPos, 'segmentStrand', segment.strand

        segTrxs, segTrxIds = segment.get_segment_trxs()
        print 'Segment transcript ids', segTrxIds
        segLen = segment.get_len()
        segStart, segEnd = segment.queryCoordinates
        reverse = False

        if segment.strand == '-':
            reverse = True

        trxOffset = segStart + xOffset
        segTrxIter = 0
        for segTrx in segTrxs:
            yCoord = yCoord + ((i + 0.5) * 0.25)
            trxLen = float(segLen) / float(len(segTrxs))
            print 'TRX len', trxLen
            trxOffset += segTrxIter * (trxLen)
            # rect = patches.Rectangle((trxOffset, yCoord + 0.15), trxLen, 0.05, color=segment.color)
            # ax.add_patch(rect)
            print 'TRX offset', trxOffset
            trx = segTrx.trx
            brkpts = segTrx.brkpts
            exons = sorted(trx.exons, key=lambda x: x.start)

            parsedExons = []
            for exon in exons:
                parsedExons.append((int(exon.start), int(exon.stop), 'exon'))

            bpPlotBins = []
            for brkpt in brkpts:
                print 'SV breakpoints for segTrx', brkpt.dist, brkpt.svBrkpt.chrom, brkpt.svBrkpt.svType, brkpt.svBrkpt.genomicCoords[brkpt.brkptIdx], brkpt.brkptIdx, segment.strand
                gCoord = brkpt.get_genomic_coord()
                if gCoord < trx.start or gCoord > trx.stop:
                    parsedExons.append((int(gCoord) - 1, int(gCoord), 'breakpoint'))
                else:
                    for i, exon in enumerate(exons):
                        if gCoord >= exon.start and gCoord <= exon.stop:
                            # within exon
                            print 'Gcoord, exon', gCoord, exon.start, exon.stop, i
                            bpPlotBins.append(('exon', i))
                            break
                        elif gCoord < exon.start:
                            print 'Gcoord, exon.start', gCoord, exon.start, i
                            bpPlotBins.append(('intron', i - 1))
                            break

            newExons = sorted(parsedExons, key=lambda x: x[0])

            binSize = trxLen / (2 * len(exons) - 1)
            offset = trxOffset
            ycoord = int(yCoord) - (float(segTrxIter) / float(5))
            # labelStr = trx.geneName + ':' + trx.id + ' (' + trx.strand + ')'
            # ax.text(trxOffset + (float(trxLen) / float(2)), yCoord + 2, labelStr, ha='center', va='center', size=12)
            ax.text(trxOffset, yCoord - 0.25, str(trx.start), ha='left', va='center', size=10)
            ax.text(trxOffset + trxLen, yCoord - 0.25, str(trx.stop), ha='right', va='center', size=10)
            exonLabel = 'exon1'
            if trx.strand == '-':
                exonLabel = 'exon' + str(len(exons))
            ax.text(trxOffset, yCoord + 0.5, exonLabel, ha='left', va='center', size=8)
            trxElements = []
            for i, exon in enumerate(newExons):
                rectLen = binSize
                start = offset
                color = segment.color
                height = 0.35
                if exon[2] == 'breakpoint':
                    rectLen = 0.5
                    height = 5
                    exonStr = ''
                    if i == (len(plotExons) - 1):
                        start += binSize
                    ax.vlines(x=start, ymin=yCoord - 0.35, ymax=yCoord + 0.35, color='grey', linewidth=1.5, zorder=2)
                    if int(exon[0]) >= int(trx.start) and int(exon[1]) <= int(trx.stop):
                        trxElements.append(start)
                offset += binSize + rectLen + (binSize - rectLen)
                if exon[2] != 'breakpoint':
                    rect = patches.Rectangle((start, yCoord), rectLen, height, color=color)
                    ax.add_patch(rect)
                    trxElements.append(start)
                    trxElements.append(start + binSize)
            segTrxIter += 1
            # This guarantees that intergenic breakpoints don't appear to be in the transcript.
            print 'TRX elements', trxElements, trxOffset, trxLen
            trxMin = max(min(trxElements), trxOffset)
            trxMax = min(max(trxElements), trxOffset + trxLen)
            print 'TRX max, min', trxMin, trxMax
            rect = patches.Rectangle((trxMin, yCoord + 0.125), trxMax - trxMin, 0.125, color=segment.color)
            ax.add_patch(rect)

            for bp in bpPlotBins:
                print 'BP', bp
                mult = 1
                if bp[0] == 'exon':
                    multiplier = 2
                start = trxOffset + (binSize * mult * bp[1]) + (float(binSize) / float(2))
                print 'Start coord', start
                ax.vlines(x=start, ymin=yCoord - 0.35, ymax=yCoord + 0.35, color='grey', linewidth=1.5, zorder=2)


def plot_annotation_track(ax, yCoord, xOffset, segmentManager):
    """ """
    if not segmentManager.has_annotations():
        return

    segStarts = []
    for i, segment in enumerate(segmentManager.segments):
        segStarts.append((segment.queryCoordinates[0], segment))

    sortedSegs = sorted(segStarts, key=lambda x: x[0])

    for i, segmentTuple in enumerate(sortedSegs):
        print 'sv_viz.py plot_annotation_track segment', i
        segment = segmentTuple[1]
        segmentPos = 'only'
        if len(sortedSegs) > 1:
            if i == 0:
                segmentPos = 'first'
            elif i > 0 and i < (len(sortedSegs) - 1):
                segmentPos = 'middle'
            elif i == (len(sortedSegs) - 1):
                segmentPos = 'last'
        print 'segment position', segmentPos, 'segmentStrand', segment.strand

        segTrxs, segTrxIds = segment.get_segment_trxs()
        print 'Segment transcript ids', segTrxIds
        segLen = segment.get_len()
        segStart, segEnd = segment.queryCoordinates
        reverse = False

        if segment.strand == '-':
            reverse = True

        trxOffset = segStart + xOffset
        if (segmentPos == 'first' or segmentPos == 'only'):
            trxOffset += 3
        segTrxIter = 0
        for segTrx in segTrxs:
            if (segmentPos == 'first' or segmentPos == 'only') and segTrxIter == 0:
                segLen = segLen - 3
                for i in range(3):
                    rect = patches.Rectangle((trxOffset - 3 + i, yCoord + 0.125), 0.25, 0.1, color=segment.color)
                    ax.add_patch(rect)
            if segTrxIter == (len(segTrxs) - 1) and (segmentPos == 'last' or segmentPos == 'only'):
                # Last segment and trx
                segLen = segLen - 3
                print 'HELLO', '@'*20
                for i in range(3):
                    rect = patches.Rectangle((trxOffset + 0.5 + segLen + i, yCoord + 0.125), 0.2, 0.1, color=segment.color)
                    ax.add_patch(rect)

            trxLen = float(segLen) / float(len(segTrxs))
            print 'TRX len', trxLen
            trxOffset += segTrxIter * (trxLen)
            print 'TRX offset', trxOffset
            trx = segTrx.trx
            brkpts = segTrx.brkpts
            trx_reverse = False
            if trx.strand == '-':
                trx_reverse = True
            exons = sorted(trx.exons, key=lambda x: x.start, reverse=trx_reverse)

            for brkpt in brkpts:
                print 'SV breakpoints for segTrx', brkpt.dist, brkpt.svBrkpt.chrom, brkpt.svBrkpt.svType, brkpt.svBrkpt.genomicCoords[brkpt.brkptIdx], brkpt.brkptIdx, segment.strand

            abrkpt = determine_annotation_brkpts(segTrx.brkpts, segmentPos, segment.strand)
            selectedExons = abrkpt.select_exons(exons)
            print selectedExons

            # genomicLen = log(abs(maxminCoords[0] - maxminCoords[1]), 2)
            # bpUnits = float(trxLen) / float(genomicLen)

            mergedExons = []
            for item in selectedExons:
                mergedExons.extend(selectedExons[item]['coords'])
            allExons = sorted(mergedExons, key=lambda x: x[0], reverse=reverse)
            print 'All Exons sorted', 10 * '#'
            print allExons
            plotExons = get_neighbor_exons(allExons)
            print plotExons

            binSize = trxLen / (2 * len(plotExons) - 1)
            offset = trxOffset
            ycoord = int(yCoord) - (float(segTrxIter) / float(5))
            labelStr = trx.geneName + ':' + trx.id + ' (' + trx.strand + ')'
            ax.text(trxOffset + (float(trxLen) / float(2)), yCoord + 2, labelStr, ha='center', va='center', size=12)
            trxElements = []
            for i, exon in enumerate(plotExons):
                rectLen = binSize
                start = offset
                color = segment.color
                height = 0.5
                exonStr = exon[2]
                if exon[2] == 'breakpoint':
                    rectLen = 0.5
                    height = 5
                    color = 'black'
                    exonStr = ''
                    if i == (len(plotExons) - 1):
                        start += binSize #- rectLen
                    ax.vlines(x=start, ymin=0.2, ymax=yCoord + 0.5, color='grey', linewidth=1.5, zorder=2)
                    if int(exon[0]) >= int(trx.start) and int(exon[1]) <= int(trx.stop):
                        trxElements.append(start)
                offset += binSize + rectLen + (binSize - rectLen)
                print 'Rect plot coords', start, yCoord, start + rectLen, binSize
                if exon[2] != 'breakpoint':
                    rect = patches.Rectangle((start, yCoord), rectLen, height, color=color)
                    ax.add_patch(rect)
                    print 'Exon', exon
                    ax.text(start + (float(binSize) / float(2)), yCoord + 0.75, exonStr, ha='center', va='center', size=8)

                if exonStr != '':
                    exstart = exon[0]
                    exend = exon[1]
                    if segment.strand == '-':
                        exstart = exon[1]
                        exend = exon[0]
                    exstart = segment.chromName + ':' + str(exstart)
                    exend = segment.chromName + ':' + str(exend)
                    ax.text(start, yCoord - 0.25, str(exstart), ha='left', va='center', size=8)
                    ax.text(start + binSize, yCoord - 0.25, str(exend), ha='right', va='center', size=8)
                    if int(exon[0]) >= int(trx.start) and int(exon[1]) <= int(trx.stop):
                        trxElements.append(start)
                        trxElements.append(start + binSize)
                if exon[3] is not None:
                    if i == (len(plotExons) - 1):
                        start += binSize
                    ax.vlines(x=start, ymin=0.2, ymax=yCoord + 0.5, color='grey', linewidth=1.5, zorder=2)
            # This guarantees that intergenic breakpoints don't appear to be in the transcript.
            trxMin = max(min(trxElements), trxOffset)
            trxMax = min(max(trxElements), trxOffset + trxLen)
            rect = patches.Rectangle((trxMin, yCoord + 0.125), trxMax - trxMin, 0.125, color=segment.color)
            ax.add_patch(rect)
                    # rect = patches.Rectangle((start, yCoord), 0.1, 5, color='black')
                    # ax.add_patch(rect)


            # for exon in selectedExons:
            #     genomicStart = maxminCoords[2]
            #     if maxminCoords[3] == 'all':
            #         genomicStart = int(trx.start)
            #     # if reverse:
            #     #     genomicStart = maxminCoords[1]
            #     # ll = [log(int(exon[0]), 2), log(int(exon[1]), 2)]
            #     e1 = log(max(abs(int(genomicStart) - int(exon[0])), 1), 2) * bpUnits
            #     e2 = log(max(abs(int(genomicStart) - int(exon[1])), 1), 2) * bpUnits
            #     if segmentPos == 'first':
            #         e1 = trxLen - e1
            #         e2 = trxLen - e2
            #     print 'genomic start', genomicStart, bpUnits
            #     print 'Exon', exon[0], exon[1], trxOffset + e1, exon[0] - exon[1]
            #     eCoords = [e1, e2]
            #     print 'Mapped coords', e1, e2
            #     eCoords.sort()
            #     if reverse:
            #         print 'Reversing exon coords'
            #         tmp = e2
            #         e2 = e1
            #         e1 = tmp
            #         print 'Mapped exon coords', e1, e2
            #     ycoord = int(yCoord) - (float(segTrxIter) / float(5))
            #     color = segment.color
            #     rectLen = e2 - e1
            #     if exon[2] == 'breakpoint':
            #         color = 'black'
            #         # rectLen = 0.25

            #     rect = patches.Rectangle((trxOffset + e1, ycoord), rectLen, 1, color=color)
            #     ax.add_patch(rect)
            #     if exon[2] != '' and exon[2] != 'breakpoint':
            #         ax.text(trxOffset + e1, ycoord, exon[2], ha='center', va='top', size=8)
            # ax.text(trxOffset, yCoord + 1, trx.strand, ha='center', va='top', size=10)
            segTrxIter += 1


# def map_coding_features(exons):
#     bp = self.get_brkpt()
#     reverse = False
#     print self.realignStrand
#     if self.realignStrand == '-':
#         reverse = True
#     l = sorted(exons, key=lambda x: x[1], reverse=reverse)
#     genomicStart = self.get_brkpt()

#     print self.pos
#     print l
#     if self.pos == "first":
#         end = l[0][0]
#         genomicStart = end
#     elif self.pos == "last":
#         end = l[-1][1]

#     print exons
#     print end, bp
#     genomicLen = log(abs(end - int(bp)), 2)
#     bpUnits = float(self.segment.plotLen)/float(genomicLen)
#     print self.segment.plotLen
#     print genomicLen
#     print bpUnits
#     pos = []

#     previousCoord = 0
#     if self.pos == "last":
#         previousCoord = log(int(bp), 2)

#     intronLen = 0
#     for coords in l:
#         if coords[2].find('exon') == -1:
#             continue
#         ll = [log(coords[0], 2), log(coords[1], 2)]
#         e1 = log(max(abs(int(genomicStart) - int(coords[0])), 1), 2) * bpUnits
#         e2 = log(max(abs(int(genomicStart) - int(coords[1])), 1), 2) * bpUnits
#         eCoords = [e1, e2]
#         eCoords.sort()
#         eCoords.append(coords[2])
#         pos.append(eCoords)
#     self.mappedExons = pos
#     self.genomicEnd = self.bp.split(':')[0] + ':' + str(end)

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
