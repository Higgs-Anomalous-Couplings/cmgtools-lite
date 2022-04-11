#!/usr/bin/env python
# USAGE: 
# 1. [shapes comparison] ==> python vbfH-diphoton/vbfhgg_plots.py plots/2022-04-11-vbfPresel 2016 shape --sP 'vbfDNN.*,D.*,dipho_mass_narrow'
# 2. [data/MC scaled] ==> python vbfH-diphoton/vbfhgg_plots.py plots/2022-04-11-vbfPresel 2016 lumi --sP 'vbfDNN.*,D.*,dipho_mass_narrow'

import sys
import re
import os

ODIR=sys.argv[1]
YEAR=sys.argv[2]

dowhat = "plots"
#dowhat = "yields"

ORIGIN="/eos/cms/store/group/dpg_ecal/comm_ecal/localreco/vbfhgg/HiggsCouplings/Trees_21032022";
lumiMap = {'2016':36.33, '2017':41.48, '2018':59.35, 'combined':138, 'merged':138}

def base(mode):
    CORE = "-P '{ORIGIN}/{YEAR}' -F 'Friends' '{{P}}/friends/{{cname}}_Friend.root' -L vbfH-diphoton/functionsH2G.cc ".format(ORIGIN=ORIGIN,YEAR=YEAR)
    CORE += " -f -j 8 -l %s " % lumiMap[YEAR] 
    RATIO  = " --maxRatioRange 0.5 1.5 --ratioYNDiv 505 "
    RATIO2 = " --showRatio --attachRatioPanel --fixRatioRange "
    LEGEND = " --allProcInLegend --n-column-legend 2 --setLegendCoordinates '0.2,0.72,0.95,0.92' "
    LEGEND2 = " --legendFontSize 0.03 "
    SPAM = " --noCms --topSpamSize 1.1 --lspam '#scale[1.1]{#bf{CMS}} #scale[0.9]{#it{Preliminary}}' "
    if dowhat == "plots": CORE+=RATIO+RATIO2+LEGEND+LEGEND2+SPAM+"  --showMCError "
    
    GO="%s vbfH-diphoton/mca-vbfhgg-%s.txt vbfH-diphoton/vbfhgg.txt vbfH-diphoton/vbfhgg_plots.txt "%(CORE,YEAR)
    GO+=" -W weight --xp '.*H120,.*H130' "
    if dowhat=='plots':
        GO+= " vbfH-diphoton/vbfhgg_plots.txt "
        if mode=='shapes':
            GO+= " --xp 'data' --contentAxisTitle 'Arbitrary units' --forceFillColorNostackMode 'vbfH125,DiPhoton' --plotmode norm --ratioNums 'ggH125,vbfH0M,ttH125,wzH125,DiPhoton' --ratioDen 'vbfH125' --ratioYLabel 'X/(SM VBF)' "
        elif mode=='lumi':
            GO+= " --xp '.*ALT.*' --preFitData dipho_mass_narrow "
            GO+= " --ratioYLabel 'X/(SM VBF)' "
        else:
            raise RuntimeError, 'Unknown selection'
        return GO

def runIt(GO,name,plots=[],noplots=[]):
    if dowhat == "plots":  
        print 'python mcPlots.py',"--pdir %s/%s/%s"%(ODIR,YEAR,name),GO,' '.join(['--sP %s'%p for p in plots]),' '.join(['--xP %s'%p for p in noplots]),' '.join(sys.argv[4:])
    elif dowhat == "yields": print 'echo %s; python mcAnalysis.py'%name,GO,' '.join(sys.argv[4:])
def add(GO,opt):
    return '%s %s'%(GO,opt)

if __name__ == '__main__':

    torun = sys.argv[3]
    if 'shapes' in torun:
        x = base('shapes')
    if 'lumi' in torun:
        x = base('lumi')
    runIt(x,'%s'%torun)
