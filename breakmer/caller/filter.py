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
        for line in open(self.filterFn, 'rU'):
            line = line.strip()
            resultFilter = Filter(line.split('\t'))
            self.filters.append(resultFilter)
