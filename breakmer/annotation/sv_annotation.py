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


# class segAnnot():
#     def __init__(self, annoRes, bp, realignStrand, pos, segLen, segment):
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


def annotate_event(svEventResult, contigMeta):
    """ """
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
    print 'sv_annotation.py bpMap', bpMap
    outputFiles = run_bedtools(bedtools, annotationFn, brkptBedFn, contigMeta.path)
    trxMap = parse_bedtools_output(outputFiles)
    store_annotations(bpMap, trxMap, annotationFn, contigMeta.params, contigMeta.path)
    # Remove temporary bedtools output files.


def store_annotations(bpMap, trxMap, annotationFn, params, tmpFilePath):
    for bpKey in bpMap:
        blatResult, svBrkptIdx, coordIdx = bpMap[bpKey]
        if bpKey not in trxMap:
            print 'Missing a breakpoint annotation', bpKey
        else:
            svBreakpoint = blatResult.get_sv_brkpts()[svBrkptIdx]
            trxMappings = trxMap[bpKey]
            intersect = trxMap[bpKey]['intersect']
            upstream = trxMap[bpKey]['upstream']
            downstream = trxMap[bpKey]['downstream']
            if intersect is not None:
                trx, dist = intersect
                if params.get_param('generate_image') or True:
                    trx.get_exons(annotationFn, tmpFilePath)
                blatResult.get_sv_brkpts()[svBrkptIdx].store_annotation([trx], [dist], coordIdx)
            else:
                upTrx, upDist = upstream
                downTrx, downDist = downstream
                if params.get_param('generate_image') or True:
                    upTrx.get_exons(annotationFn, tmpFilePath)
                    downTrx.get_exons(annotationFn, tmpFilePath)
                blatResult.get_sv_brkpts()[svBrkptIdx].store_annotation([downTrx, upTrx], [upDist, downDist], coordIdx)


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
            # brkptKey = 'BP' + str(bpIter) + '|' + chrom + ':' + '-'.join([str(x) for x in brkptCoords])
            coordIdx = 0
            for coord in brkptCoords:
                bpKey = chrom + ':' + str(coord) + '_BP' + str(bpIter)
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
        trx = Transcript(linesplit[4:])
        dist = int(linesplit[-1])
        if bpKey not in trxMap:
            trxMap[bpKey] = {'intersect': None, 'upstream': None, 'downstream': None}

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


# def get_exons(annotationFn, )

# #         chr, bp = bp.split(":")
# #         if bp.find('-') > -1:
# #             return None
# #         nearestGene = [None, None]
# #         for achr in self.trxs:
# #             if achr == chr:
# #                 for gene in self.trxs[achr]:
# #                     retVals = {}
# #                     trxSizes = sorted(self.trxs[achr][gene]['trxSizes'], key=lambda x: x[0])
# #                     canonTrxId = trxSizes[-1][1]
# #                     trxStart = self.trxs[achr][gene]['trx'][canonTrxId]['bounds'][0]
# #                     trxStop = self.trxs[achr][gene]['trx'][canonTrxId]['bounds'][1]
# #                     #print bp, trxStart, trxStop
# #                     #print trx, trxStart, trxStop, self.trxs[achr][trx]['strand']
# #                     retVals['features'] = self.trxs[achr][gene]['trx'][canonTrxId]['features']
# #                     retVals['name'] = self.trxs[achr][gene]['name']
# #                     retVals['strand'] = self.trxs[achr][gene]['strand']
# #                     retVals['trxId'] = canonTrxId
# #                     if int(bp) >= int(trxStart) and int(bp) <= int(trxStop):
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