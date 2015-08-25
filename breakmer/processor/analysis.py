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


def analyze_targets(targetList):
    """Analyze a list of targets.

    A list of TargetManager objects are passed in to be analyzed independently.
    Each target ref data is set, if necessary, then the reads are extracted,
    contigs built, and calls made.

    This function performs all the top level functions on the target regions being analyzed.

    Args:
        targetList (list):          A list of TargetManager objects, representing target regions.
    Returns:
        aggregateResults (dict):    A dictoinary containing lists of formatted output strings for the
                                    contig-based calls and the discordant-only read clusters.
    Raises:
        None
    """

    aggregateResults = {'contigs': [], 'discreads': []}  # Formatted output strings for contig based calls and discordant read calls are different.
    for targetRegion in targetList:
        print 'Analyzing', targetRegion.name
        utils.log('breakmer.processor.analysis', 'info', 'Analyzing %s' % targetRegion.name)
        targetRegion.set_ref_data()
        if targetRegion.fnc == 'prepare_reference_data':  # Stop here if only preparing ref data.
            continue
        if not targetRegion.find_sv_reads():  # No SV reads extracted. Exiting.
            continue
        targetRegion.compare_kmers()  # Perform kmer subtraction.
        targetRegion.resolve_sv()  # Assemble extracted reads and make calls.
        if targetRegion.has_results():
            outputs = targetRegion.get_formatted_output()
            for key in outputs:
                aggregateResults[key].extend(outputs[key])
        targetRegion.complete_analysis()  # Write results out to file.
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
        self.__params = params
        self.__results = []
        # self.__targets = {}
        self.__loggingName = 'breakmer.processor.analysis'

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

        startTime = time.clock()  # Track the run time.

        self.__params.start_blat_server()
        if self.__params.fncCmd == 'start_blat_server':
            print 'Server started!'
            return

        targetAnalysisList = self.__create_targets()

        aggResults = {'contigs': [], 'discreads': []}  # Buffer the formatted output strings for each target to write out in batch.
        nprocs = int(self.__params.get_param('nprocs'))
        if nprocs > 1:  # Make use of multiprocessing by mapping targets to n jobs.
            utils.log(self.loggingName, 'info', 'Creating all reference data.')
            p = multiprocessing.Pool(nprocs)
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

        if self.__params.fncCmd == 'prepare_reference_data':
            print 'Reference data setup!'
            return

        self.__write_aggregated_output(aggResults)
        utils.log(self.loggingName, 'info', 'Analysis complete in %s' % str(time.clock() - startTime))

        if not self.__params.get_param('keep_blat_server'):  # Keep blat server is specified.
            cmd = '%s stop %s %d' % (self.__params.get_param('gfserver'), self.__params.get_param('blat_hostname'), int(self.__params.get_param('blat_port')))
            os.system(cmd)
        print 'Analysis complete!'

    def __create_targets(self):
        """Create target objects and group them by the number of
        multiprocs that are specified (i.e. n=1 for 1 processor.)

        If multiprocs are used then split the list into n batches:
        [1,2,3,4,5,....,20] = [[targetGroup1], [targetGroup2],...]

        N target groups (+1 if there is a remainder).

        Store all the TargetManager instances in self.targets dictionary.
        self.targets[<target name>] = TargetManager instance.

        Args:
            None
        Returns:
            trgtGroups (list): A list of lists containing target objects. Each toplevel list is
                               analyzed by a processor.
        """

        nprocs = int(self.__params.get_param('nprocs'))
        multiprocs = nprocs > 1
        ngroups = nprocs
        ntargets = len(self.__params.targets)
        ntargetsPerGroup = ntargets / nprocs
        modval = math.fmod(ntargets, nprocs)
        if modval > 0:
            ngroups += 1
        trgtGroups = []
        trgtGroup = []

        # Iterate through the target name list, sorted alphabetically.
        targetNames = self.__params.get_target_names()
        targetNames.sort()
        for targetName in targetNames:
            targetManager = target.TargetManager(targetName, self.__params)
            if multiprocs:
                if len(trgtGroup) == ntargetsPerGroup:
                    trgtGroups.append(trgtGroup)
                    trgtGroup = []
                trgtGroup.append(targetManager)
            else:
                trgtGroups.append(targetManager)

        # For the last batch, check if there are less elements than each group has
        # if so, then extend the last group to add them, otherwise create a new batch.
        if multiprocs:
            if len(trgtGroup) < ntargetsPerGroup:
                trgtGroups[-1].extend(trgtGroup)
            else:
                trgtGroups.append(trgtGroup)
        return trgtGroups

    def __write_aggregated_output(self, aggregateResults):
        """Write the SV calls to a top level file in the specified output directory.
        Header is written at the top of the file if option to remove is not
        specified.

        The output files are:
            <output_dir>/<analysis_name>_svs.out
            <output_dir>/<analysis_name>_discreads.out

        Args:
            aggregateResults (dict):    A dictionary containing the formatted output string values.
        Returns:
            None
        """

        # Write assembled contig-based SV calls.
        if len(aggregateResults['contigs']) > 0:
            resultFn = os.path.join(self.__params.paths['output'], self.__params.get_param('analysis_name') + "_svs.out")
            utils.log(self.loggingName, 'info', 'Writing %s aggregated results file %s' % (self.__params.get_param('analysis_name'), resultFn))
            resultFile = open(resultFn, 'w')
            for i, formattedResultStr in enumerate(aggregateResults['contigs']):
                headerStr, formattedResultValuesStr = formattedResultStr
                if not self.__params.get_param('no_output_header') and i == 0:
                    resultFile.write(headerStr + '\n')
                resultFile.write(formattedResultValuesStr + '\n')
            resultFile.close()

        # Write discordant read pair clusters.
        if len(aggregateResults['discreads']) > 0:
            resultFn = os.path.join(self.__params.paths['output'], self.__params.get_param('analysis_name') + "_discreads.out")
            utils.log(self.loggingName, 'info', 'Writing %s aggregated results file %s' % (self.__params.get_param('analysis_name'), resultFn))
            resultFile = open(resultFn, 'w')
            for i, formattedResultStr in enumerate(aggregateResults['discreads']):
                headerStr, formattedResultValuesStr = formattedResultStr
                if not self.__params.get_param('no_output_header') and i == 0:
                    resultFile.write(headerStr + '\n')
                resultFile.write(formattedResultValuesStr + '\n')
            resultFile.close()
