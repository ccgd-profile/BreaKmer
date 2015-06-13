#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
import breakmer.realignment.blat_result as blat_result
import breakmer.utils as utils

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


class AlignParams:
    """
    """
    def __init__(self, params, targetrefFns):
        self.program = {'target': 'blat', 'genome': 'blat'}
        self.extension = {'target': 'psl', 'genome': 'psl'}
        self.binary = {'target': None, 'genome': None}
        self.ref = {'target': None, 'genome': None}
        self.set_values(params, targetRefFns)

    def set_values(self, params, targetRefFns):
        """
        """
        blast = params.get_param('blast')
        if blast:
            self.program['target'] = 'blast'
            self.binary['target'] = blast
            self.extension['target'] = 'xml'

        blat = params.get_param('blat')
        self.binary['genome'] = {'gfclient': params.get_param('gfclient'),
                                 'hostname': params.get_param('blat_hostname'),
                                 'port': params.get_param('blat_port')
                                 }
        # Use the forward sequence for blatting targeted sequences
        self.ref['target'] = targetRefFns[0]
        self.ref['genome'] = params.get_param('reference_fasta_dir')

    def get_values(self, type):
        return (self.program[type], self.extension[type], self.binary[type])


class RealignManager:
    """
    """
    def __init__(self, params, targetRefFns):
        self.realignment = None
        self.alignParams = AlignParams(params, targetRefFns)

    def realign(contig):
        """
        """
        if not contig.has_fa_fn():
            return

        self.realignment = Realignment(contig)
        if not self.realignment.align(self.alignParams.get_values('target'), 'target'):
            # Log
            return
        if not self.realignment.target_aligned():
            self.realignment.align(self.alignParams.get_values('genome'), 'genome')


class Realignment:
    """
    """
    def __init__(self, contig):
        self.scope = None
        self.results = None
        self.contig = contig

    def align(self, alignParams, scope):
        """
        """
        alignProgram, alignExt, alignBinary = alignParams
        self.scope = scope
        # update
        self.logger.info('Running blat %s, storing results in %s' % (self.params.opts['gfclient'], self.query_res_fn))

        resultFn = os.path.join(contig.get_path(), '%s_res.%s.%s' % (alignProgram, scope, alignExt))
        self.results = AlignResults(alignProgram, scope, resultFn)

        cmd = ''
        if alignprogram == 'blast':
            cmd = ''
        elif self.alignprogram == 'blat':
            if scope == 'target':
                # all blat server
                cmd = '%s -t=dna -q=dna -out=psl -minScore=20 -nohead localhost %d %s %s %s' % (self.params.opts['gfclient'], self.params.opts['blat_port'], self.params.opts['reference_fasta_dir'], self.contig_fa_fn, self.query_res_fn)
            elif scope == 'genome':
                # target
                cmd = '%s -t=dna -q=dna -out=psl -minScore=20 -stepSize=10 -minMatch=2 -repeats=lower -noHead %s %s %s' % (self.params.opts['blat'], db, self.contig_fa_fn, self.query_res_fn)

        # update
        utils.log(self.loggingName, 'info', 'Realignment system command %s' % cmd)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, errors = p.communicate()
        utils.log(self.loggingName, 'info', 'Realignment output file %s' % output)
        if errors != '':
            utils.log(self.loggingName, 'info', 'Realignment errors %s' % errors)

        if not os.path.isfile(resultFn):
            return False
        else:
            self.results = AlignResults(alignProgram, scope, resultFn)


class AlignResults:
    def __init__(self, program, scope, alignResultFn):
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
        self.clippedQs = []
        self.set_values()

    def set_values(self):
        if not self.resultFn:
            self.hasResults = False
        elif len(open(self.resultFn, 'rU').readlines()) == 0:
            self.hasResults = False
        else:
            self.parse_result_file()

    def parse_result_file(self):
        if self.program == 'blat':
            self.parse_blat_results()
        elif self.program == 'blast':
            self.parse_blast_results()

    def parse_blat_results(self):
        """
        """
        refName = None
        offset = None
        if self.scope == 'target':
            # Need to reset the chrom name and coordinates for blat results.
            refName = self.contig.get_chr()
            offset = self.contig.get_target_start() - self.contig.get_target_buffer()

        for line in open(self.resultFn, 'r'):
            line = line.strip()
            parsedBlatResult = blat_result.BlatResult(line.split('\t'), refName, offset)
            parsedBlatResult.set_gene_annotations(self.contig.get_target_region_coordinates(), self.contig.get_gene_annotations())
            parsedBlatResult.set_repeats(self.contig.get_repeat_annotations())
            self.process_blat_result(parsedBlatResult)

            # Move this to blatResult class.
            # score_raw = br.get_nmatch_total()
            # ngaps = br.get_ngap_total()
            # in_target = 1 if br.in_target else 0
            # score_frac = float(score_raw) / float(br.get_size('query'))
            # score = score_raw + score_frac
            # perc_ident = br.perc_ident
            self.results.append(parsedBlatResults) #(score, ngaps, in_target, br, perc_ident))
        # Update to use class attributes as sorting categories
        self.blat_results = sorted(self.blat_results, key=lambda blat_results: (-blat_results[0], -blat_results[4], blat_results[1]))

    def process_blat_result(self, blatResultObj):
        """Summarize metrics from all alignments.
        """
        self.nmismatches += blatResultObj.get_nmatches('mismatch')
        self.ngaps += blatResultObj.get_num_gaps()
        if not self.querySize:
            self.querySize = blatResultObj.get_seq_size('query')
            self.alignmentFreq = [0] * self.querySize
        for i in range(blatResultObj.qstart(), blatResultObj.qend()):
            self.alignmentFreq[i] += 1
