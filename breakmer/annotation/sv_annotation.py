#! /usr/bin/local/python
# -*- coding: utf-8 -*-

import os
import shutil
import subprocess
import breakmer.utils as utils

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


class Exon:
    def __init__(self, values):
        self.chr = ''
        self.start = ''
        self.stop = ''
        self.featureType = ''
        self.set_values(values)

    def set_values(self, values):
        """ """
        self.chr, self.src, featureType, self.start, self.stop, fill, self.strand, fill2, meta = values
        self.start = int(self.start)
        self.stop = int(self.stop)
        self.featureType = featureType


class Transcript:
    def __init__(self, values):
        self.chr = ''
        self.src = ''
        self.start = ''
        self.stop = ''
        self.strand = ''
        self.id = ''
        self.geneName = ''
        self.geneId = ''
        self.geneStatus = ''
        self.len = 0
        self.exons = []
        self.set_values(values)

    def set_values(self, values):
        """ """
        self.chr, self.src, featureType, self.start, self.stop, fill, self.strand, fill2, meta, dist = values
        self.start = int(self.start)
        self.stop = int(self.stop)
        meta = meta.split(';')
        self.id = meta[1].split(' ')[2].lstrip('"').rstrip('"')
        self.geneId = meta[0].split(' ')[1].lstrip('"').rstrip('"')
        self.geneName = meta[4].split(' ')[2].lstrip('"').rstrip('"')
        self.geneStatus = meta[3].split(' ')[2].lstrip('"').rstrip('"')
        self.len = int(self.stop) - int(self.start)

    def get_exons(self, annotationFn, tmpFilePath):
        """ """
        # Grep the exons and UTRs from the annotationFn
        exonSelect = '$3 == "exon"' #' || $3 == "UTR")'
        outFn = os.path.join(tmpFilePath, self.id + '.exons')
        cmd = 'cat ' + annotationFn + " | awk '" + exonSelect + "' | grep '" + self.id + "' > " + os.path.join(tmpFilePath, self.id + '.exons')
        os.system(cmd)
        # ' 'bedtools multicov -bams ' + args.bam + ' -bed ' + args.intervals
        for line in open(outFn, 'r'):
            self.exons.append(Exon(line.strip().split('\t')))
        os.remove(outFn)


def annotate_event(svEventResult, contigMeta):
    """ """
    if svEventResult.is_filtered():
        svEventResult.annotated = False
    else:
        svEventResult.annotated = True
        # Make sure annotation file is sorted for bedtools use.
        bedtools = contigMeta.params.get_param('bedtools')
        annotationFn = contigMeta.params.get_param('gene_annotation_file')
        brkptBedFn = os.path.join(contigMeta.path, contigMeta.id + '_breakpoints.bed')

        # Dictionary with 'targets' and 'other' breakpoint lists
        # Deletions have two breakpoints in reference.
        # Insertions have one breakpoint in reference.
        # Rearrangements have breakpoints for each segment that is rearranged.
        #  genomicBrkpts = svEventResult.get_genomic_brkpts()
        bpMap = write_brkpt_bed_file(brkptBedFn, svEventResult.blatResults)
        # print 'sv_annotation.py bpMap', bpMap
        outputFiles = run_bedtools(bedtools, annotationFn, brkptBedFn, contigMeta.path)
        trxMap = parse_bedtools_output(outputFiles)
        store_annotations(svEventResult, bpMap, trxMap, annotationFn, contigMeta.params, contigMeta.path)
        # Remove temporary bedtools output files.
        # print 'annotate_event, sv_annotation.py', svEventResult
        # svEventResult.set_annotations()
        # print 'svEvent annotated', svEventResult.annotated


def store_annotations(svEventResult, bpMap, trxMap, annotationFn, params, tmpFilePath):
    for bpKey in bpMap:
        blatResult, svBrkptIdx, coordIdx = bpMap[bpKey]
        # print 'sv_annotation store_annotations', bpKey, bpMap[bpKey]
        if bpKey not in trxMap:
            print 'Missing a breakpoint annotation', bpKey
            svEventResult.set_failed_annotation()
            svEventResult.set_filtered('Breakpoints are not fully annotated. Typically due to non-primary chromosome.')
        else:
            svBreakpoint = blatResult.get_sv_brkpts()[svBrkptIdx]
            trxMappings = trxMap[bpKey]
            intersect = trxMap[bpKey]['intersect']
            upstream = trxMap[bpKey]['upstream']
            downstream = trxMap[bpKey]['downstream']
            # print 'Intersect', intersect
            # print 'Downstream', downstream
            # print 'Upstream', upstream
            if intersect is not None:
                trx, dist = intersect
                if params.get_param('generate_image') or True:
                    trx.get_exons(annotationFn, tmpFilePath)
                # print blatResult, blatResult.get_sv_brkpts()
                blatResult.get_sv_brkpts()[svBrkptIdx].store_annotation([trx], [dist], coordIdx)
            else:
                upTrx, upDist = upstream
                downTrx, downDist = downstream
                # print 'Up', upTrx.id, upDist
                # print 'Down', downTrx.id, downDist
                if params.get_param('generate_image') or True:
                    upTrx.get_exons(annotationFn, tmpFilePath)
                    downTrx.get_exons(annotationFn, tmpFilePath)
                blatResult.get_sv_brkpts()[svBrkptIdx].store_annotation([upTrx, downTrx], [upDist, downDist], coordIdx)


def write_brkpt_bed_file(bpBedFn, blatResults):
    """ """
    bpMap = {}
    bpBedFile = open(bpBedFn, 'w')
    bpIter = 1
    for queryStartCoord, blatResult in blatResults:
        svBreakpoints = blatResult.get_sv_brkpts()
        svBrkptIdx = 0
        for svBreakpoint in svBreakpoints:
            chrom = svBreakpoint.chrom
            brkptCoords = svBreakpoint.genomicCoords
            # print 'write_brkpt_bed_file', chrom, brkptCoords, svBreakpoint.svType
            # brkptKey = 'BP' + str(bpIter) + '|' + chrom + ':' + '-'.join([str(x) for x in brkptCoords])
            coordIdx = 0
            for coord in brkptCoords:
                bpKey = chrom + ':' + str(coord) + '_BP' + str(bpIter) + '_' + str(svBrkptIdx)
                # print 'write_brkpt_bed_file', bpKey
                bpStr = [chrom, coord, int(coord) + 1, bpKey]
                bpBedFile.write('\t'.join([str(x) for x in bpStr]) + '\n')
                bpMap[bpKey] = (blatResult, svBrkptIdx, coordIdx)
                coordIdx += 1
            svBrkptIdx += 1
        bpIter += 1
    bpBedFile.close()
    cmd = 'sort -k1,1 -k2,2n %s > %s' % (bpBedFn, bpBedFn + '.sorted')
    os.system(cmd)
    shutil.move(bpBedFn + '.sorted', bpBedFn)
    return bpMap


def run_bedtools(bedtools, annotationFn, brkptBedFn, tmpFilePath):
    """ """

    # Identify the transcripts first
    trxSelect = '$3 == "transcript"'
    knownGeneSelect = 'gene_status "KNOWN"'

    outputFiles = {'intersect': os.path.join(tmpFilePath, 'bedtools.intersect.out'),
                   'upstream': os.path.join(tmpFilePath, 'bedtools.upstream.out'),
                   'downstream': os.path.join(tmpFilePath, 'bedtools.downstream.out')}
    # Intersecting transcripts
    cmd = 'cat ' + annotationFn + " | awk '" + trxSelect + "' | grep '" + knownGeneSelect + "' | " + bedtools + ' intersect -wo -a %s -b - > %s' % (brkptBedFn, outputFiles['intersect'])
    os.system(cmd)
    # Upstream transcripts
    cmd = 'cat ' + annotationFn + " | awk '" + trxSelect + "' | grep '" + knownGeneSelect + "' | " + bedtools + ' closest -D a -id -a %s -b - > %s' % (brkptBedFn, outputFiles['upstream'])
    os.system(cmd)
    # Downstream transcripts
    cmd = 'cat ' + annotationFn + " | awk '" + trxSelect + "' | grep '" + knownGeneSelect + "' | " + bedtools + ' closest -D a -iu -a %s -b - > %s' % (brkptBedFn, outputFiles['downstream'])
    os.system(cmd)
    return outputFiles


def parse_bedtools_file(fn, fileKey, trxMap):
    for line in open(fn, 'r'):
        line = line.strip()
        linesplit = line.split('\t')
        bpChrom, bpStart, bpEnd, bpKey = linesplit[0:4]

        # No value found for this breakpoint. This could be due to the chromsome not existing in the
        # annotation file.
        if linesplit[4] == '.':
            return

        if bpKey not in trxMap:
            trxMap[bpKey] = {'intersect': None, 'upstream': None, 'downstream': None}
        trx = Transcript(linesplit[4:])
        dist = int(linesplit[-1])
        checkStorage = ((fileKey != 'intersect') and (trxMap[bpKey]['intersect'] is None)) or (fileKey == 'intersect')
        if checkStorage:
            if trxMap[bpKey][fileKey] is None:
                trxMap[bpKey][fileKey] = [trx, dist]
            else:
                # Check if trx is longer (i.e. canonical) vs. current stored
                if trx.len > trxMap[bpKey][fileKey][0].len:
                    trxMap[bpKey][fileKey] = [trx, dist]


def parse_bedtools_output(outputFileDict):
    """ """
    trxMap = {}
    # Map each bp to a transcript (or two) if it is intergenic.
    parse_bedtools_file(outputFileDict['intersect'], 'intersect', trxMap)
    parse_bedtools_file(outputFileDict['upstream'], 'upstream', trxMap)
    parse_bedtools_file(outputFileDict['downstream'], 'downstream', trxMap)
    return trxMap
