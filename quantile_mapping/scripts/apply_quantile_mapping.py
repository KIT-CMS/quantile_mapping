#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT
import numpy
import argparse
from quantile_mapping.quantile_mapping.QuantileShifter import QuantileShifter
import logging

logger = logging.getLogger(__name__)


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
    parser.add_argument(
        "-l",
        "--linear-extrapolation-threshold",
        type=float,
        default=0.0,
        help="Fraction of events in the lower and upper rim respectively for that linear interpolation is used instead of splines to approximate CDFs")

    return parser.parse_args()

def main(args):
    shifter = QuantileShifter(args.correction_file, args.source, args.target, args.bisect_method)
    rootfile = ROOT.TFile(args.input_file, "UPDATE")
    tree = rootfile.Get(args.tree)        
    var = numpy.zeros(1, dtype=float)
    branch = tree.Branch(args.new_name, var, "%s/D"%args.new_name)
    for event in tree:
        var[0] = shifter.shift(getattr(event, args.variable), args.linear_extrapolation_threshold)
        branch.Fill()
    tree.Write("", ROOT.TObject.kOverwrite)
    rootfile.Close()

if __name__ == "__main__":
    args = parse_arguments()
    setup_logging("quantile_shifter.log", logging.WARNING)
    main(args)
