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
        blast = params.get_param('blast')
        if blast :
            self.program['target'] = 'blast'
            self.binary['target'] = blast
            self.extension['target'] = 'xml'

        blat = params.get_param('blat')
        self.binary['genome'] = {'gfclient':params.get_param('gfclient'), 'hostname':params.get_param('blat_hostname'), 'port':params.get_param('blat_port')}
        self.ref['target'] = target.files['target_ref_fn']
        self.ref['genome'] = params.get_param('reference_fasta_dir')

    def get_values(self, type) :
        return (self.program[type], self.extension[type], self.binary[type])

class Realignment :
    def __init__(self, contig) :
        self.scope = None
        self.result_manager = None

    def align(self, align_params, scope) :
        """
        """
        self.scope = scope

        alignment_program, extension, bin = align_params

        self.result_fn = os.path.join(path, '%s_res.%s.%s'%(alignment_program, type, extension))
        if alignment_program == 'blast' :
            cmd = 
        elif alignment_program == 'blat' :
            cmd = 

        if not os.path.isfile(self.result_fn) :
            return False
        else :
            self.result_manager = AlignResults(self.result_fn, alignment_program, scope)


class AlignResults :
    def __init__(self, align_result_fn, program, scope) :
        self.result_fn = align_result_fn
        self.program = program
        self.scope = scope
        self.query_size = 0
        self.alignment_freq = []
        self.nmismatches = 0
        self.ngaps = 0
        self.has_results = True
        self.results = []
        self.set_values()

    def set_values(self) :
        if not self.result_fn :
            self.has_results = False
        elif len(open(self.result_fn, 'rU').readlines()) == 0 :
            self.has_results = False
        else :
            self.parse_result_file()

    def parse_result_file(self) :
        if self.program == 'blat' :
            
        elif self.program == 'blast' :
            self.parse_blast_results()

    def parse_blat_results(self) :
        for line in open(self.result_fn, 'rU').readlines() :
            line = 


class RealignManager :
    """
    """

    def __init__(self, params, target) :
        self.realignment_scope = None
        self.align_params = AlignParams(params, target)

    def realign(contig, params) :
        """
        """

        if not contig.has_fa_fn() : 
            return

        realignment = Realignment(contig)
        if not realignment.align(self.align_params.get_values('target'), 'target') :
            # Log
            return
        if not realignment.target_aligned() :
            realignment.align(self.align_params.get_values('genome'), 'genome')
