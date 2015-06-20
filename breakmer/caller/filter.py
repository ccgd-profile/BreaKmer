#! /usr/bin/local/python
# -*- coding: utf-8 -*-

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


class Filter:
    """
    """

    def __init(self, name, type, breakpoints, description):
        self.name = name
        self.type = type
        self.breakpoints = breakpoints
        self.description = description


class ResultFilter:
    """
    """

    def __init__(self, filterFn):
        self.filterFn = filterFn
        self.filters = []
        self.setup()

    def setup(self):
        if self.filterFn:
            for line in open(self.filterFn, 'rU'):
                line = line.strip()
                resultFilter = Filter(line.split('\t'))
                self.filters.append(resultFilter)

# Filters for events
# 1. Check this for non-indels:
            # nMissingQueryCoverage = len(filter(lambda y: y, map(lambda x: x == 0, self.queryCoverage)))
            # if nMissingQueryCoverage < self.meta_dict['params'].get_min_segment_length('trl'):
            #     valid = True
# 2. Contig complexity 
        # avg_comp, comp_vec = calc_contig_complexity(self.contig_seq)


    def check_filters(self, svEvent):
        """
        """
        