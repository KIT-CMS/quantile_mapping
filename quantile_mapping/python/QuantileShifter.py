#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT
import copy
import logging

logger = logging.getLogger(__name__)

class QuantileShifter(object):
    def __init__(self, inputfile, source, target, use_bisect=False):
        percentagefile = ROOT.TFile(inputfile, "READ")
        self._source = copy.deepcopy(percentagefile.Get(source))
        self._target = copy.deepcopy(percentagefile.Get(target))
        percentagefile.Close()
        self._use_bisect = use_bisect
        if self._use_bisect:
            logger.info("Using bisect.")
    
    def _bisect(self, percentage, up , down, steps):
        if steps <= 0:
            return down + (percentage - self._target.Eval(down)) / (self._target.Eval(up) - self._target.Eval(down)) * (up - down)
        middle = (up + down) / 2.0
        if percentage > self._target.Eval(middle):
            return self._bisect(percentage, up, middle, steps - 1)
        else:
            return self._bisect(percentage, middle, down, steps - 1)
    
    def _FindY(self, percentage, search_in_target): #similar to FindX method of TSpline
        bindown = 0
        binup = self._target.GetNp() - 1 if search_in_target else self._source.GetNp() - 1
    
        while binup - bindown > 1:
            x = ROOT.Double()
            y = ROOT.Double()
            binhalf = int((binup + bindown) / 2)
            if search_in_target:
                self._target.GetKnot(binhalf, x, y)
            else:
                self._source.GetKnot(binhalf, x, y)
            if y < percentage:
                bindown = binhalf;
            else:
                binup = binhalf;
        
        return bindown

    
    def shift(self, value, linear_interpolation_threshold = 0.0):
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
        
        percentage = max(0.0, min(self._source.Eval(value), 1.0))
        
        if percentage < linear_interpolation_threshold or percentage > 1.0 - linear_interpolation_threshold:
            knot = self._source.FindX(value)
            if knot >= npoints - 1: #this can happen if an input exactly hits the upper boundary
                npoints = self._target.GetNp()
                self._target.GetKnot(npoints-1, xup, yup)
                return xup
            else:
                self._source.GetKnot(knot, xdown, ydown)
                self._source.GetKnot(knot + 1, xup, yup);
                percentage = ydown + (value - xdown) * (yup - ydown) / (xup - xdown)
                
                knot = self._FindY(percentage, True)
                self._target.GetKnot(knot, xdown, ydown)
                self._target.GetKnot(knot + 1, xup, yup)
                if ydown == yup:
                    return xdown
                return xdown + (percentage - ydown) / (yup - ydown) * (xup - xdown)
        
        self._target.GetKnot(0, xdown, ydown)
        for index in range(npoints-1):
            self._target.GetKnot(index, xdown, ydown)
            self._target.GetKnot(index+1, xup, yup)
            if percentage <= yup and yup != 0.0:
                break            
                
        steps = 5
        if self._use_bisect:
            return self._bisect(percentage, xup, xdown, steps)
        
        result = xdown + (percentage - ydown) / (yup - ydown) * (xup - xdown)
        derivative = self._target.Derivative(result)
        if derivative == 0.0:
            logger.warning("Default inversion method fails due to zero derivative. Bisecting is used instead.")
            return self._bisect(percentage, xup, xdown, steps)
        correction = (percentage - self._target.Eval(result)) / self._target.Derivative(result)
        result += correction
        if (abs(correction) > (xup - xdown) / 2.0 or result < xdown or result > xup):
            logger.warning("Default inversion method yields too large corrections. Bisecting is used instead.")
            return self._bisect(percentage, xup, xdown, steps)
        return result