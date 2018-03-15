#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT
import argparse
from quantile_mapping.quantile_mapping.QuantileShifter import QuantileShifter
from array import array

import logging
logger = logging.getLogger("")

def parse_arguments():
    parser = argparse.ArgumentParser(
        description=
        "Derive cumulative distribution function as splines from histograms."
    )
    parser.add_argument(
        "-c",
        "--correction-file",
        type=str,
        required=True,
        help="ROOT file with CDF splines")
    parser.add_argument(
        "-s",
        "--sources",
        nargs="+",
        type=str,
        required=True,
        help="Names of source CDF splines")
    parser.add_argument(
        "-t",
        "--targets",
        nargs="+",
        type=str,
        default=None,
        help="Names of target CDF splines")
    parser.add_argument(
        "--labels",
        nargs="+",
        type=str,
        default=None,
        help="Custom legend labels. Use '{default}' in your string, in case you want to include default label.")
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="qm_transformations",
        help="Output file name (.png is automatically appended). Default: %(default)s")
    parser.add_argument(
        "-b",
        "--bisect-method",
        action="store_true",
        help="Always use bisect method for the inversion of the target CDF spline")
    parser.add_argument(
        "-bi",
        "--bisect-method-individual",
        nargs="+",
        choices=['0', '1'],
        help="Use bisect method for the inversion of the target CDF spline for indiviual inputs. Use as many 0/1 flags as inputs given. Overrides -b option.")
    parser.add_argument(
        "-l",
        "--linear-interpolation-threshold",
        nargs="+",
        type=float,
        default=[0.0],
        help="Fraction of events in the lower and upper rim respectively for that linear interpolation is used instead of splines to approximate CDFs. Default: 0.0")
    parser.add_argument(
        "--x-min",
        type=float,
        default=None,
        help="lower boundary of x-range")
    parser.add_argument(
        "--x-max",
        type=float,
        default=None,
        help="upper boundary of x-range")

    return parser.parse_args()

def setup_logging(output_file, level=logging.DEBUG):
    logger.setLevel(level)
    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    file_handler = logging.FileHandler(output_file, "w")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

def SetCanvasStyle():
    # For the canvas:
    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
    ROOT.gStyle.SetCanvasDefH(600)  # Height of canvas
    ROOT.gStyle.SetCanvasDefW(600)  # Width of canvas
    ROOT.gStyle.SetCanvasDefX(0)  # POsition on screen
    ROOT.gStyle.SetCanvasDefY(0)

    # For the Pad:
    ROOT.gStyle.SetPadBorderMode(0)
    ROOT.gStyle.SetPadColor(ROOT.kWhite)
    ROOT.gStyle.SetPadGridX(False)
    ROOT.gStyle.SetPadGridY(False)
    ROOT.gStyle.SetGridColor(1)
    ROOT.gStyle.SetGridStyle(3)
    ROOT.gStyle.SetGridWidth(1)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    ROOT.gStyle.SetPadRightMargin(0.05)

    # For the frame:
    ROOT.gStyle.SetFrameBorderMode(0)
    ROOT.gStyle.SetFrameBorderSize(1)
    ROOT.gStyle.SetFrameFillColor(0)
    ROOT.gStyle.SetFrameFillStyle(0)
    ROOT.gStyle.SetFrameLineColor(1)
    ROOT.gStyle.SetFrameLineStyle(1)
    ROOT.gStyle.SetFrameLineWidth(1)

def main(args):
    # check and prepare inputs
    if len(args.sources)!=len(args.targets):
        logger.FATAL("Number of sources and targets must be equal!")
        raise Exception
    if args.labels!=None and len(args.labels)!=len(args.targets):
        logger.FATAL("Number of legend labels must be equal to number of sources/targets!")
        raise Exception
    if len(args.linear_interpolation_threshold)!=len(args.targets) and len(args.linear_interpolation_threshold)!=1:
        logger.FATAL("Number of -l arguments must be either 1 or equal to number of sources/targets!")
        raise Exception
    lit = args.linear_interpolation_threshold * len(args.targets) if len(args.linear_interpolation_threshold)==1 else args.linear_interpolation_threshold
    bmi = [args.bisect_method] * len(args.targets)
    if args.bisect_method_individual!=None:
        if len(args.bisect_method_individual)!=len(args.targets):
            logger.FATAL("Number of -bi arguments must be equal to number of sources/targets!")
            raise Exception
        else:
            bmi = [(True if entry == '1' else False) for entry in args.bisect_method_individual]
    
    # book quantile shifters, scan input range and store in TGraphs
    graphs = []
    for source, target, lit_arg, bmi_arg in zip(args.sources, args.targets, lit, bmi):
        shifter = QuantileShifter(args.correction_file, source, target, bmi_arg)
        x_min = shifter._source.GetXmin() if args.x_min==None else max(args.x_min, shifter._source.GetXmin())
        x_range = (shifter._source.GetXmax() if args.x_max==None else min(args.x_max, shifter._source.GetXmax())) - x_min
        x_vals = array('d')
        y_vals = array('d')
        for i in range(101):
            x_vals.append(ROOT.Double(x_min + i * x_range / 100.))
            y_vals.append(ROOT.Double(shifter.shift(x_vals[i], lit_arg)) - x_vals[i])
        graphs.append(ROOT.TGraph(101, x_vals, y_vals))
    
    # create histogram with event distribution
    shifter = QuantileShifter(args.correction_file, args.sources[0], args.targets[0], bmi[0])
    x_min = shifter._source.GetXmin() if args.x_min==None else max(args.x_min, shifter._source.GetXmin())
    x_max = shifter._source.GetXmax() if args.x_max==None else min(args.x_max, shifter._source.GetXmax())
    nbins = 40
    hist = ROOT.TH1F("h1", "h1", nbins, x_min, x_max)
    for i in range(nbins):
        hist.SetBinContent(i+1, shifter._source.Derivative(x_min + (i+0.5) / nbins * (x_max - x_min)))
    global_max = max(gr.GetHistogram().GetMaximum() for gr in graphs)
    global_min = min(gr.GetHistogram().GetMinimum() for gr in graphs)
    rescale = global_max / hist.GetMaximum() * 0.9 #graphs[0].GetHistogram().GetMaximum() / hist.GetMaximum() * 0.9
    hist.Scale(rescale)
    hist.SetFillColorAlpha(ROOT.kBlack, 0.35)
    hist.SetLineColorAlpha(ROOT.kBlack, 0.35)
    
    # draw graphs
    SetCanvasStyle()
    canvas = ROOT.TCanvas()
    canvas.SetGrid()
    graphs[0].GetYaxis().SetRangeUser(global_min, global_min + 4./3. * (global_max - global_min))
    graphs[0].GetXaxis().SetTitle("input (arb. units)")
    graphs[0].GetYaxis().SetTitle("output - input (arb. units)")
    graphs[0].Draw()
    graphs[0].SetTitle("")
    hist.Draw("same hist")
    for i, graph in enumerate(graphs):
        graph.SetLineColor(i+1)
        graph.SetLineWidth(2)
        graph.Draw("same")
    legend = ROOT.TLegend(0.15,0.7,0.95,0.9, '', 'NBNDC')
    legend.AddEntry(hist, "expected event distribution", "f")
    for i, graph in enumerate(graphs):
        defaultlabel = "%s#rightarrow %s"%(args.sources[i], args.targets[i])
        if bmi[i] and lit[i]>0.0 and lit[i]<=0.5:
            defaultlabel += " (#it{b, l = %3.1f})" %min(lit[i], 1.0)
        elif lit[i]>0.0:
            defaultlabel += " (#it{l=%3.1f})"%min(lit[i], 1.0)
        elif bmi[i]:
            defaultlabel += " (#it{b})"
        legend.AddEntry(graph, defaultlabel if args.labels==None else args.labels[i].format(default=defaultlabel), "l")
    legend.Draw("same")
    canvas.SaveAs("%s.png"%args.output)

if __name__ == "__main__":
    args = parse_arguments()
    setup_logging("cumdistfunc.log", logging.WARNING)
    main(args)
