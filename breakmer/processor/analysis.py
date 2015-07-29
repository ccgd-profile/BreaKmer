#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
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


def check_status(results):
    """Check the status of the multiple processors running analysis

    Args:
        results (list): A list of results from the multiprocessing analysis.
    Returns:
        notReady (int): An integer value indicating the number of processors not complete.
    Raises:
        None
    """

    notReady = 0
    for r in results:
        if not r.ready():
            notReady += 1
    return (notReady)


def wait(results):
    """Check if the nprocs are complete.

    Args:
        results (list): A list of results from the multiprocessing analysis.
    Returns:
        None
    Raises:
        None
    """

    njobs = check_status(results)
    while njobs > 0:
        time.sleep(10)
        jobs = check_status(results)
        if jobs < njobs:
            njobs = jobs
        # else:
        #     sys.stdout.write('.')
        #     sys.stdout.flush()


def analyze_targets(targetList):
    """Analyze a list of targets.
    A list of TargetManager objects are passed in to be analyzed independently.
    Each target ref data is set, if necessary, then the reads are extracted,
    contigs built, and calls made.

    Args:
        targetList (list):          A list of TargetManager objects, representing target regions.
    Returns:
        aggregateResults (dict):    A dictoinary containing lists of formatted output strings for the
                                    contig-based calls and the discordant-only read clusters.
    Raises:
        None
    """

    aggregateResults = {'contigs': [], 'discreads': []}
    for targetRegion in targetList:
        print 'Analyzing', targetRegion.name
        utils.log('breakmer.processor.analysis', 'info', 'Analyzing %s' % targetRegion.name)
        targetRegion.set_ref_data()
        if targetRegion.params.fncCmd == 'prepare_reference_data':  # Stop here if only preparing ref data.
            continue
        if not targetRegion.find_sv_reads():  # No SV reads extracted. Exiting.
            continue
        targetRegion.compare_kmers()
        targetRegion.resolve_sv()
        if targetRegion.has_results():
            outputs = targetRegion.get_formatted_output()
            for key in outputs:
                aggregateResults[key].extend(outputs[key])
        targetRegion.complete_analysis()
    return aggregateResults


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
        self.loggingName = 'breakmer.processor.analysis'

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

        self.params.start_blat_server()
        if self.params.fncCmd == 'start_blat_server':
            print 'Server started!'
            return

        targetAnalysisList = self.create_targets()

        aggResults = {'contigs': [], 'discreads': []}
        nprocs = int(self.params.get_param('nprocs'))
        if nprocs > 1:
            utils.log(self.loggingName, 'info', 'Creating all reference data.')
            p = multiprocessing.Pool(nprocs)
            # p.map(analyze_targets, targetAnalysisList)
            multiprocResults = []
            for targetList in targetAnalysisList:
                multiprocResults.append(p.apply_async(analyze_targets, (targetList, )))
            wait(multiprocResults)
            for multiprocResult in multiprocResults:
                a = multiprocResult.get()
                aggResults['contigs'].extend(a['contigs'])
                aggResults['discreads'].extend(a['discreads'])
        else:
            aggResults = analyze_targets(targetAnalysisList)

        if self.params.fncCmd == 'prepare_reference_data':
            print 'Reference data setup!'
            return

        self.write_aggregated_output(aggResults)

        utils.log(self.loggingName, 'info', 'Analysis complete in %s' % str(time.clock() - startTime))

        if not self.params.get_param('keep_blat_server'):
            cmd = '%s stop %s %d' % (self.params.opts['gfserver'], self.params.get_param('blat_hostname'), int(self.params.get_param('blat_port')))
            os.system(cmd)

        print 'Analysis complete!'

    def create_targets(self):
        """Create target objects and group them by the number of
        multiprocs that are specified (i.e. n=1 for 1 processor.)

        Args:
            None
        Returns:
            trgtGroups (list): A list of lists containing target objects. Each toplevel list is
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

        targetNames = self.params.get_target_names()
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

    def write_aggregated_output(self, aggregateResults):
        """Write the information for all results to file.
        Header is written at the top of the file if option to remove is not
        specified.

        Args:
            aggregateResults (dict):    A dictionary containing the formatted output string values.
        Returns:
            None
        """

        if len(aggregateResults['contigs']) > 0:
            resultFn = os.path.join(self.params.paths['output'], self.params.opts['analysis_name'] + "_svs.out")
            utils.log(self.loggingName, 'info', 'Writing %s aggregated results file %s' % (self.params.opts['analysis_name'], resultFn))
            resultFile = open(resultFn, 'w')
            for i, formattedResultStr in enumerate(aggregateResults['contigs']):
                headerStr, formattedResultValuesStr = formattedResultStr
                if not self.params.get_param('no_output_header') and i == 0:
                    resultFile.write(headerStr + '\n')
                resultFile.write(formattedResultValuesStr + '\n')
            resultFile.close()

        if len(aggregateResults['discreads']) > 0:
            resultFn = os.path.join(self.params.paths['output'], self.params.opts['analysis_name'] + "_discreads.out")
            utils.log(self.loggingName, 'info', 'Writing %s aggregated results file %s' % (self.params.opts['analysis_name'], resultFn))
            resultFile = open(resultFn, 'w')
            for i, formattedResultStr in enumerate(aggregateResults['discreads']):
                headerStr, formattedResultValuesStr = formattedResultStr
                if not self.params.get_param('no_output_header') and i == 0:
                    resultFile.write(headerStr + '\n')
                resultFile.write(formattedResultValuesStr + '\n')
            resultFile.close()
