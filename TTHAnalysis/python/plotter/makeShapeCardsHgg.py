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
parser.add_option("--onlyfracs", dest="onlyfracs", action="store_true", default=False, help="Only print the category fractions, do not create the workspaces")
(options, args) = parser.parse_args()
options.weight = True
options.final  = True

if "/functions_cc.so" not in ROOT.gSystem.GetLibraries(): 
    ROOT.gROOT.ProcessLine(".L %s/src/CMGTools/TTHAnalysis/python/plotter/functions.cc+" % os.environ['CMSSW_BASE']);

mca  = MCAnalysis(args[0],options)
cuts = CutsFile(args[1],options)
fitvar = args[2]
fitbins = [float(s.replace("[","").replace("]","")) for s in args[3].split(",")]
proc = args[4]

if len(options.processes)>1 or options.processes[0]!=proc:
    raise RuntimeError("There should be only one process, and the 4th argument coincide with -p <proc> argument")

binname = os.path.basename(args[1]).replace(".txt","") if options.binname == 'default' else options.binname
if binname[0] in "1234567890": raise RuntimeError("Bins should start with a letter.")
outdir  = options.outdir+"/" if options.outdir else ""
if not os.path.exists(outdir): os.mkdir(outdir)

if "H" in proc:
    plotsraw = {}
    if options.categ:
        allplotsraw = dict()
        cexpr, cbins, clabels = options.categ
        catbins = [float(s.replace("[","").replace("]","")) for s in cbins.split(",")]
        catlabels = clabels.split(",")
        if len(catlabels) != len(catbins)-1: raise RuntimeError("Mismatch between category labels and bins")
        plotsraw = mca.getPlotsRaw("x", cexpr+":"+args[2], makeBinningProductString(args[3],cbins), cuts.allCuts(), nodata=options.asimov) 
        for p,h in plotsraw.iteritems(): h.cropNegativeBins()
        for ibin,binname in enumerate(catlabels):
            allplotsraw["%s_%s"%(proc,binname)] = plotsraw[proc].projectionX("x_"+proc,ibin+1,ibin+1)

    allyields = dict([(p,h.Integral()) for p,h in allplotsraw.iteritems()])
    procyield = sum([y for k,y in allyields.iteritems()])

    lines=[]
    for k,y in allyields.iteritems():
        lines.append("globalXSBRMap['AC']['%s'] = {'mode':'%s','factor':%.4f}" % (k,proc,y/procyield))
    
    fracstxt = open(outdir+proc+"_fracs.txt","w")
    for l in sorted(lines):
        fracstxt.write(l+'\n')
    print "Wrote fractions to: "+outdir+proc+"_fracs.txt"

    if options.onlyfracs:
        print "Category fractions DONE. Exiting without making the workspaces."
        exit(0)
    

# make the RooDataSets
MH=0
if "H" in proc:
    prod,MH = proc.split("H")
    outfilename = "{d}/output_M{MH}_13TeV_{p}.root".format(d=outdir,MH=MH,p=prod)
    MH=int(MH)
    procId = "{prod}_{MH}_13TeV".format(prod=prod,MH=MH)
else:
    outfilename = "{d}/output_Data_13TeV.root".format(d=outdir)
    procId = "Data_13TeV"
outfile = ROOT.TFile.Open(outfilename, "RECREATE")
outfile.mkdir("tagsDumper")
outfile.cd("tagsDumper")
wsp = ROOT.RooWorkspace("cms_hgg_13TeV")


vars = ["{fitv}[{xmin},{xmax}]".format(fitv=fitvar,xmin=fitbins[0],xmax=fitbins[-1]), # fit variable (mgg)
        "dZ[-2000,2000]",
        "centralObjectWeight[-999999.,999999.]",
        "weight[-999999.,999999.]"
    ]

allreports = {}
if options.categ:
    cexpr, cbins, clabels = options.categ
    catbins = [float(s.replace("[","").replace("]","")) for s in cbins.split(",")]
    catlabels = clabels.split(",")
    if len(catlabels) != len(catbins)-1: raise RuntimeError("Mismatch between category labels and bins")
    reports = {}

    for ibin,binname in enumerate(catlabels):
        cexprfull = "({catexpr}=={ic})*({cut})".format(catexpr=cexpr,ic=ibin+1,cut=cuts.allCuts())
        rdsname = "{p}_{cat}".format(p=procId,cat=binname)
        reports = mca.getRooDataSet(rdsname,vars,cexprfull,proc)
        print "got RooDataSet named %-40s for: %-20s category: %-40s: %8d entries, equiv to %.2f weighted events." % (rdsname,proc,binname,reports[proc].numEntries(),reports[proc].sumEntries())
        # WARNING! This is very error prone, but it is linked to the way it HAS to be run:
        # only 1 proc / command. Then take the only report
        allreports[rdsname] = reports[proc]
else:
    print "Not yet implemented"


outfile.cd("tagsDumper")
for catname, report in allreports.iteritems():
    getattr(wsp,'import')(report)
wsp.Print()

wsp.Write()
outfile.Close()
ROOT.gDirectory.Add(wsp) # needed not to crash
print "Wrote to ",outfilename

