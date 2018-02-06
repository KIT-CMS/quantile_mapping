#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT
import copy
import logging

logger = logging.getLogger(__name__)

class QuantileShifter(object):
    def __init__(self, inputfile, source, target, use_bisect=False):
        quantilefile = ROOT.TFile(inputfile, "READ")
        self._source = copy.deepcopy(quantilefile.Get(source))
        self._target = copy.deepcopy(quantilefile.Get(target))
        quantilefile.Close()
        self._use_bisect = use_bisect
        if self._use_bisect:
            logger.info("Using bisect.")
    
    def _bisect(self, quantile, up , down, steps):
        if steps <= 0:
            return down + (quantile - self._target.Eval(down)) / (self._target.Eval(up) - self._target.Eval(down)) * (up - down)
        middle = (up + down) / 2.0
        if quantile > self._target.Eval(middle):
            return self._bisect(quantile, up, middle, steps - 1)
        else:
            return self._bisect(quantile, middle, down, steps - 1)
    
    def shift(self, value):
        xup = ROOT.Double()
        yup = ROOT.Double()
        xdown = ROOT.Double()
        ydown = ROOT.Double()
        
        npoints = self._target.GetNp()
        self._source.GetKnot(0, xdown, ydown)
        self._source.GetKnot(npoints-1, xup, yup)
        if value < xdown or value > xup:
            logger.warning("Input value %f out of range [%f, %f]. No correction applied."%(value, xdown, xup))
            return value
        
        quantile = max(0.0, min(self._source.Eval(value), 1.0))
        self._target.GetKnot(0, xdown, ydown)
        for index in range(npoints-1):
            self._target.GetKnot(index, xdown, ydown)
            self._target.GetKnot(index+1, xup, yup)
            if quantile <= yup and yup != 0.0:
                break            
                
        steps = 5
        if self._use_bisect:
            return self._bisect(quantile, xup, xdown, steps)
        
        result = xdown + (quantile - ydown) / (yup - ydown) * (xup - xdown)
        derivative = self._target.Derivative(result)
        if derivative == 0.0:
            logger.warning("Default inversion method fails due to zero derivative. Bisecting is used instead.")
            return self._bisect(quantile, xup, xdown, steps)
        correction = (quantile - self._target.Eval(result)) / self._target.Derivative(result)
        result += correction
        if (abs(correction) > (xup - xdown) / 2.0 or result < xdown or result > xup):
            logger.warning("Default inversion method yields too large corrections. Bisecting is used instead.")
            return self._bisect(quantile, xup, xdown, steps)
        return result
