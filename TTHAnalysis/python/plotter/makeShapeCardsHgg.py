#!/usr/bin/env python
from CMGTools.TTHAnalysis.plotter.mcAnalysis import *
from CMGTools.TTHAnalysis.plotter.histoWithNuisances import _cloneNoDir
import re, sys, os, os.path
systs = {}

from optparse import OptionParser
parser = OptionParser(usage="%prog [options] mc.txt cuts.txt var bins proc")
addMCAnalysisOptions(parser)
parser.add_option("--od", "--outdir", dest="outdir", type="string", default=None, help="output directory name") 
parser.add_option("--asimov", dest="asimov", type="string", default=None, help="Use an Asimov dataset of the specified kind: including signal ('signal','s','sig','s+b') or background-only ('background','bkg','b','b-only')")
parser.add_option("--bbb", dest="bbb", type="string", default=None, help="Options for bin-by-bin statistical uncertainties with the specified nuisance name")
parser.add_option("--amc", "--autoMCStats", dest="autoMCStats", action="store_true", default=False, help="use autoMCStats")
parser.add_option("--autoMCStatsThreshold", dest="autoMCStatsValue", type="int", default=10, help="threshold to put on autoMCStats")
parser.add_option("--infile", dest="infile", action="store_true", default=False, help="Read histograms to file")
parser.add_option("--savefile", dest="savefile", action="store_true", default=False, help="Save histos to file")
parser.add_option("--categorize", dest="categ", type="string", nargs=3, default=None, help="Split in categories. Requires 3 arguments: expression, binning, bin labels")
parser.add_option("--regularize", dest="regularize", action="store_true", default=False, help="Regularize templates")
(options, args) = parser.parse_args()
options.weight = True
options.final  = True

if "/functions_cc.so" not in ROOT.gSystem.GetLibraries(): 
    ROOT.gROOT.ProcessLine(".L %s/src/CMGTools/TTHAnalysis/python/plotter/functions.cc+" % os.environ['CMSSW_BASE']);

mca  = MCAnalysis(args[0],options)
cuts = CutsFile(args[1],options)
proc = args[4]

binname = os.path.basename(args[1]).replace(".txt","") if options.binname == 'default' else options.binname
if binname[0] in "1234567890": raise RuntimeError("Bins should start with a letter.")
outdir  = options.outdir+"/" if options.outdir else ""
if not os.path.exists(outdir): os.mkdir(outdir)

procId = "_".join(os.path.basename(mca.treename).split("_")[:-1])

report={}
if options.infile:
    infile = ROOT.TFile(outdir+binname+".bare.root","read")
    for p in mca.listSignals(True)+mca.listBackgrounds(True)+['data']:
        variations = mca.getProcessNuisances(p) if p != "data" else []
        h = readHistoWithNuisances(infile, "x_"+p, variations, mayBeMissing=True)
        if h: report[p] = h
else:
    if options.categ:
       cexpr, cbins, _ = options.categ
       report = mca.getPlotsRaw("x", cexpr+":"+args[2], makeBinningProductString(args[3],cbins), cuts.allCuts(), nodata=options.asimov) 
    else:
       report = mca.getPlotsRaw("x", args[2], args[3], cuts.allCuts(), nodata=options.asimov) 
    for p,h in report.iteritems(): h.cropNegativeBins()

if options.savefile:
    savefile = ROOT.TFile(outdir+binname+".bare.root","recreate")
    for k,h in report.iteritems(): 
        h.writeToFile(savefile, takeOwnership=False)
    savefile.Close()

if options.asimov:
    if options.asimov in ("s","sig","signal","s+b"):
        asimovprocesses = mca.listSignals() + mca.listBackgrounds()
    elif options.asimov in ("b","bkg","background", "b-only"):
        asimovprocesses = mca.listBackgrounds()
    else: raise RuntimeError("the --asimov option requires to specify signal/sig/s/s+b or background/bkg/b/b-only")
    tomerge = None
    for p in asimovprocesses:
        if p in report: 
            if tomerge is None: 
                tomerge = report[p].raw().Clone("x_data_obs"); tomerge.SetDirectory(None)
            else: tomerge.Add(report[p].raw())
    report['data_obs'] = HistoWithNuisances(tomerge)
else:
    report['data_obs'] = report['data'].Clone("x_data_obs") 

if options.categ:
    allreports = dict()
    catlabels = options.categ[2].split(",")
    if len(catlabels) != report["data_obs"].GetNbinsY(): raise RuntimeError("Mismatch between category labels and bins")
    for ic,lab in enumerate(catlabels):
        allreports["%s_%s"%(binname,lab)] = dict( (k, h.projectionX("x_"+k,ic+1,ic+1)) for (k,h) in report.iteritems() )
else:
    allreports = {binname:report}

outfile = ROOT.TFile.Open("{d}/output_{p}.root".format(d=outdir,p=proc), "RECREATE")
outfile.mkdir("tagsDumper")
outfile.cd("tagsDumper")
wsp = ROOT.RooWorkspace("cms_hgg_13TeV")

for catname, report in allreports.iteritems():
  if options.bbb:
    if options.autoMCStats: raise RuntimeError("Can't use --bbb together with --amc/--autoMCStats")
    for p,h in report.iteritems(): 
      if p not in ("data", "data_obs"):
        h.addBinByBin(namePattern="%s_%s_%s_bin{bin}" % (options.bbb, catname, p), conservativePruning = True)
  for p,h in report.iteritems():
    for b in xrange(1,h.GetNbinsX()+1):
      h.SetBinError(b,min(h.GetBinContent(b),h.GetBinError(b))) # crop all uncertainties to 100% to avoid negative variations
  nuisances = sorted(listAllNuisances(report))

  allyields = dict([(p,h.Integral()) for p,h in report.iteritems()])
  procs = []; iproc = {}
  for i,s in enumerate(mca.listSignals()):
    if s != proc: continue
    if s not in allyields: continue
    if allyields[s] == 0: continue
    procs.append(s); iproc[s] = i-len(mca.listSignals())+1
  for i,b in enumerate(mca.listBackgrounds()):
    if b != proc: continue
    if b not in allyields: continue
    if allyields[b] == 0: continue
    procs.append(b); iproc[b] = i+1
  #for p in procs: print "%-10s %10.4f" % (p, allyields[p])

  if proc=='data': towrite = report["data_obs"].raw()
  else: towrite = [ report[p].raw() for p in procs ]

  name = h.GetName()
  var = ROOT.RooRealVar(args[2],args[2],h.GetXaxis().GetBinLowEdge(1),h.GetXaxis().GetBinUpEdge(h.GetNbinsX()))
  rdsname = "{p}_{cat}".format(p=procId,cat=catname.replace(binname+"_",""))
  rds = ROOT.RooDataHist(rdsname, rdsname, ROOT.RooArgList(var), h.raw())
  getattr(wsp,'import')(rds)

wsp.Write()
outfile.Close()
ROOT.gDirectory.Add(wsp) # needed not to crash
print "Wrote to {d}/output_{p}.root".format(d=outdir,p=proc)

