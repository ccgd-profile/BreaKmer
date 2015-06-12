#! /usr/bin/python
# -*- coding: utf-8 -*-

import os


__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


class AlignParams:
    """
    """
    def __init__(self, params, target):
        self.program = {'target': 'blat', 'genome': 'blat'}
        self.extension = {'target': 'psl', 'genome': 'psl'}
        self.binary = {'target': None, 'genome': None}
        self.ref = {'target': None, 'genome': None}
        self.set_values(params, target)

    def set_values(self, params, target):
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
        self.ref['target'] = target.files['target_ref_fn']
        self.ref['genome'] = params.get_param('reference_fasta_dir')

    def get_values(self, type):
        return (self.program[type], self.extension[type], self.binary[type])


class Realignment:
    """
    """
    def __init__(self, params, target, contig):
        self.scope = None
        self.results = None
        self.contig = contig
        self.alignParams = AlignParams(params, target)

    def align(self, scope):
        """
        """
        self.scope = scope
        # update
        self.logger.info('Running blat %s, storing results in %s' % (self.params.opts['gfclient'], self.query_res_fn))

        resultFn = os.path.join(contig.get_path(), '%s_res.%s.%s' % (self.alignParams.program, scope, self.alignParams.extension))
        self.results = AlignResults(self.alignParams.program, scope, resultFn)

        if self.alignParams.program == 'blast':
            cmd = ''
        elif self.alignParams.program == 'blat':
            # all blat server
            cmd = '%s -t=dna -q=dna -out=psl -minScore=20 -nohead localhost %d %s %s %s' % (self.params.opts['gfclient'], self.params.opts['blat_port'], self.params.opts['reference_fasta_dir'], self.contig_fa_fn, self.query_res_fn)
            # target
            cmd = '%s -t=dna -q=dna -out=psl -minScore=20 -stepSize=10 -minMatch=2 -repeats=lower -noHead %s %s %s' % (self.params.opts['blat'], db, self.contig_fa_fn, self.query_res_fn)

        # update
        self.logger.info('Blat system command %s' % cmd)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, errors = p.communicate()
        self.logger.info('Blat output %s' % output)
        if errors != '':
            self.logger.info('Blat errors %s' % errors)

        if not os.path.isfile(self.result_fn):
            return False
        else:
            self.result_manager = AlignResults(self.result_fn, self.alignParams.program, scope)

"""
    if self.contig_fa_fn :
      self.run_blat(target_ref_fn, 'target') # Run blat against target reference sequence first for speed.
      if not self.query_res_fn :
        self.logger.info('No blat results file %s, no calls for %s.'%(self.query_res_fn, self.id))
        return
      if not self.check_target_blat(query_region) :
        # Blat against whole genome reference fasta
        self.run_blat(self.params.opts['reference_fasta'], 'all')
  #*********************************************************

  #*********************************************************
  def run_blat(self, db, name) :
    self.query_res_fn = os.path.join(self.path,'blat_res.'+name+'.psl')
    if not os.path.isfile(self.query_res_fn) :
      if name == 'all' : 
        self.logger.info('Running blat %s, storing results in %s'%(self.params.opts['gfclient'],self.query_res_fn))
        cmd = '%s -t=dna -q=dna -out=psl -minScore=20 -nohead localhost %d %s %s %s'%(self.params.opts['gfclient'],self.params.opts['blat_port'], self.params.opts['reference_fasta_dir'], self.contig_fa_fn, self.query_res_fn)
      else :
        self.logger.info('Running blat %s, storing results in %s'%(self.params.opts['blat'],self.query_res_fn))
        cmd = '%s -t=dna -q=dna -out=psl -minScore=20 -stepSize=10 -minMatch=2 -repeats=lower -noHead %s %s %s'%(self.params.opts['blat'], db, self.contig_fa_fn, self.query_res_fn)

      self.logger.info('Blat system command %s'%cmd)
      p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
      output, errors = p.communicate()
      self.logger.info('Blat output %s'%output)
      if errors != '' : self.logger.info('Blat errors %s'%errors)
    else : self.logger.info('Blat already run, results file %s exists, continuing'%self.query_res_fn)
"""


class AlignResults:
    def __init__(self, program, scope, alignResultFn):
        self.resultFn = alignResultFn
        self.program = program
        self.scope = scope
        self.query_size = 0
        self.alignment_freq = []
        self.nmismatches = 0
        self.ngaps = 0
        self.hasResults = True
        self.results = []
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

        elif self.program == 'blast':
            self.parse_blast_results()

    def parse_blat_results(self):
        for line in open(self.resultFn, 'rU').readlines():
            line = 


class RealignManager:
    """
    """

    def __init__(self, params, target):
        self.realignment_scope = None
        self.align_params = AlignParams(params, target)

    def realign(contig, params):
        """
        """

        if not contig.has_fa_fn():
            return

        realignment = Realignment(contig)
        if not realignment.align(self.align_params.get_values('target'), 'target'):
            # Log
            return
        if not realignment.target_aligned():
            realignment.align(self.align_params.get_values('genome'), 'genome')
