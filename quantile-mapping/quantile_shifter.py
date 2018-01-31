#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT
import numpy
import argparse
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
        bin_index = 0
        for index in range(npoints-1):
            bin_index = index
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

def setup_logging(output_file, level=logging.DEBUG):
    logger.setLevel(level)
    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    file_handler = logging.FileHandler(output_file, "w")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Append new branch to a TTree with a quantile-mapped variable.")

    parser.add_argument(
        "-i",
        "--input-file",
        required=True,
        type=str,
        help="ROOT file containing the events to be corrected.")
    parser.add_argument(
        "-t",
        "--tree",
        required=True,
        type=str,
        help="Path to tree within the ROOT file (-i).")
    parser.add_argument(
        "-v",
        "--variable",
        required=True,
        type=str,
        help="Leaf of the tree to be corrected.")
    parser.add_argument(
        "-n",
        "--new-name",
        required=True,
        type=str,
        help="Name of the corrected leaf.")
    
    parser.add_argument(
        "-c",
        "--correction-file",
        required=True,
        type=str,
        help="ROOT file containg the CDF splines.")
    parser.add_argument(
        "-cs",
        "--source",
        required=True,
        type=str,
        help="Name of the source CDF spline.")
    parser.add_argument(
        "-ct",
        "--target",
        required=True,
        type=str,
        help="Name of the target CDF spline.")
    parser.add_argument(
        "-b",
        "--bisect-method",
        action="store_true",
        help="Always use bisect method for the inversion of the target CDF spline")

    return parser.parse_args()

def main(args):
    shifter = QuantileShifter(args.correction_file, args.source, args.target, args.bisect_method)
    rootfile = ROOT.TFile(args.input_file, "UPDATE")
    tree = rootfile.Get(args.tree)        
    var = numpy.zeros(1, dtype=float)
    branch = tree.Branch(args.new_name, var, "%s/D"%args.new_name)
    for event in tree:
        var[0] = shifter.shift(getattr(event, args.variable))
        branch.Fill()
    tree.Write("", ROOT.TObject.kOverwrite)
    rootfile.Close()

if __name__ == "__main__":
    args = parse_arguments()
    setup_logging("quantile_shifter.log", logging.WARNING)
    main(args)
