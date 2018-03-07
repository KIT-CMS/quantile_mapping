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
        "-o",
        "--output",
        type=str,
        default="qm_transformations",
        help="Output file name (.png is automatically appended). Default: %{default}")
    parser.add_argument(
        "-b",
        "--bisect-method",
        action="store_true",
        help="Always use bisect method for the inversion of the target CDF spline")
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
    # check inputs
    if len(args.sources)!=len(args.targets):
        logger.FATAL("Number of sources and targets must be equal!")
        raise Exception
    
    # book quantile shifters, scan input range and store in TGraphs
    graphs = []
    for source, target in zip(args.sources, args.targets):
        shifter = QuantileShifter(args.correction_file, source, target, args.bisect_method)
        x_min = shifter._source.GetXmin() if args.x_min==None else max(args.x_min, shifter._source.GetXmin())
        x_range = (shifter._source.GetXmax() if args.x_max==None else min(args.x_max, shifter._source.GetXmax())) - x_min
        x_vals = array('d')
        y_vals = array('d')
        for i in range(101):
            x_vals.append(ROOT.Double(x_min + i * x_range / 100.))
            y_vals.append(ROOT.Double(shifter.shift(x_vals[i])) - x_vals[i])
        graphs.append(ROOT.TGraph(101, x_vals, y_vals))
    
    # create histogram with event distribution
    shifter = QuantileShifter(args.correction_file, args.sources[0], args.targets[0], args.bisect_method)
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
        legend.AddEntry(graph, "%s#rightarrow %s"%(args.sources[i], args.targets[i]), "l")
    legend.Draw("same")
    canvas.SaveAs("%s.png"%args.output)

if __name__ == "__main__":
    args = parse_arguments()
    setup_logging("cumdistfunc.log", logging.WARNING)
    main(args)
