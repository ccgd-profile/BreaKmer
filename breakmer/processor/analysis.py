#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
import logging
import time
import math
import multiprocessing
import breakmer.processor.target as target
import breakmer.utils as utils

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


def analyze_targets(targetList):
    """Analyze a list of targets.
    A list of TargetManager objects are passed in to be analyzed independently.
    Each target ref data is set, if necessary, then the reads are extracted,
    contigs built, and calls made.
    Args:
        targetList: A list of TargetManager objects, representing target regions.
    Returns:
        None
    """

    for targetRegion in targetList:
        utils.log('breakmer.processor.analysis', 'info', 'Analyzing %s' % targetRegion.name)
        targetRegion.set_ref_data()

        # If only presetting data stop here.
        if targetRegion.params.get_param('preset_ref_data'):
            continue
        # targetRegion fails to find any interesting reads to use. Exiting.
        if not targetRegion.find_sv_reads():
            continue
        targetRegion.compare_kmers()
        targetRegion.resolve_sv()
        # targetRegion.complete_analysis()


class RunTracker:
    """Class to manage the running of all the target region analyses.
    The params object is passed in with all the input information.
    The run() function creates the target region objects from the
    param inputs and then starts the analysis for each target.
    Args:
        params: ParamManager object.
    Returns:
        None
    """
    def __init__(self, params):
        self.params = params
        self.results = []
        self.targets = {}
        self.summary = {}
        self.summary_header = ''
        self.logger = logging.getLogger('breakmer.processor.analysis')

    def run(self):
        """Create and analyze the target regions.
        The target objects are made and grouped for multiprocessing (if set)
        and these are all analyzed independently. This is where the analysis
        starts and ends.
        Args:
            None
        Returns
            None
        """
        startTime = time.clock()
        targetAnalysisList = self.create_targets()

        if not self.params.get_param('preset_ref_data'):
            self.params.start_blat_server()

        nprocs = int(self.params.get_param('nprocs'))
        if nprocs > 1:
            self.logger.info('Creating all reference data.')
            p = multiprocessing.Pool(nprocs)
            p.map(analyze_targets, targetAnalysisList)
        else:
            analyze_targets(targetAnalysisList)

#        self.write_output()

        # Perform any post-primary analysis scripts here.


        self.logger.info('Analysis complete in %s' % str(time.clock() - startTime))

        if not self.params.get_param('keep_blat_server'):
            cmd = '%s stop localhost %d' % (self.params.opts['gfserver'], int(self.params.get_param('blat_port')))
            os.system(cmd)
        print 'Analysis complete!'

    def create_targets(self):
        """Create target objects and group them by the number of
        multiprocs that are specified (i.e. n=1 for 1 processor.)
        Args:
            None
        Returns:
            trgtGroups: List of lists of target objects. Each list is
                        analyzed by a processor.
        """
        nprocs = int(self.params.get_param('nprocs'))
        multiprocs = nprocs > 1
        ngroups = nprocs
        ntargets = len(self.params.targets)
        ntargetsPerGroup = ntargets / nprocs
        modval = math.fmod(ntargets, nprocs)
        if modval > 0:
            ngroups += 1
        trgtGroups = []
        trgtGroup = []

        targetNames = self.params.targets.keys()
        targetNames.sort()
        for targetName in targetNames:
            trgt = target.TargetManager(targetName, self.params)
            self.targets[targetName] = trgt
            if multiprocs:
                if len(trgtGroup) == ntargetsPerGroup:
                    trgtGroups.append(trgtGroup)
                    trgtGroup = []
                trgtGroup.append(trgt)
            else:
                trgtGroups.append(trgt)

        if multiprocs:
            if len(trgtGroup) < ntargetsPerGroup:
                trgtGroups[-1].extend(trgtGroup)
            else:
                trgtGroups.append(trgtGroup)
        return trgtGroups

    def write_output(self):
        """Write the information for all results to file.
        Header is written at the top of the file if option to remove is not
        specified.
        Args:
            None
        Returns:
            None
        """

        resultFiles = {}
        for res in self.results:
            tag = res[6]
            if tag not in resultFiles:
                header = "\t".join(['genes', 'target_breakpoints', 'align_cigar', 'mismatches', 'strands', 'rep_overlap_segment_len', 'sv_type', 'split_read_count', 'nkmers', 'disc_read_count', 'breakpoint_coverages', 'contig_id', 'contig_seq']) + "\n"
                res_fn = os.path.join(self.params.paths['output'], self.params.opts['analysis_name'] + "_" + tag + "_svs.out")
                self.logger.info('Writing %s output file %s' % (tag, res_fn))
                resultFiles[tag] = open(res_fn, 'w')
                if not self.params.opts['no_output_header']:
                    resultFiles[tag].write(header)
            resultFiles[tag].write("\t".join([str(x) for x in res]) + "\n")

        for f in resultFiles:
            resultFiles[f].close()

#        summary_fn = os.path.join(self.params.paths['output'], self.params.opts['analysis_name']+"_summary.out")
#        summary_f = open(summary_fn,'w')
#        self.logger.info('Writing summary file to %s'%summary_fn)
#        summary_f.write(self.summary_header+"\n")
#        keys = self.summary.keys()
#        keys.sort()
#        for gene in keys :
#            summary_f.write(self.summary[gene]+"\n")
#        summary_f.close()
