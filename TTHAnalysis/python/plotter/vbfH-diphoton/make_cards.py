# USAGE: python vbfH-diphoton/make_cards.py cards_fithgg 2018

import os, sys
nCores=8
submit = '{command}'

ORIGIN="/eos/cms/store/group/dpg_ecal/comm_ecal/localreco/vbfhgg/HiggsCouplings/Trees_21032022"; 

if len(sys.argv) < 3: 
    print 'Syntax is %s [outputdir] [year] [other]'%sys.argv[0]
    raise RuntimeError 
OUTNAME=sys.argv[1]
YEAR=sys.argv[2]
OTHER=sys.argv[3:] if len(sys.argv) > 3 else ''

if YEAR not in ['2016','2017','2018']:
    raise RuntimeError("Wrong year %s"%YEAR)
LUMI="1" # this is what is needed by the finalfits, which re-normalizes assuming L=1fb-1

#print "Normalizing to {LUMI}/fb".format(LUMI=LUMI);
OPTIONS=" -j {J} -l {LUMI} -f ".format(LUMI=LUMI,J=nCores)
os.system("test -d cards/{OUTNAME} || mkdir -p cards/{OUTNAME}".format(OUTNAME=OUTNAME))
OPTIONS="{OPTIONS} --od cards/{OUTNAME} ".format(OPTIONS=OPTIONS, OUTNAME=OUTNAME)
T2G="-P '{ORIGIN}/{YEAR}' -F 'Friends' '{{P}}/friends/{{cname}}_Friend.root' ".format(ORIGIN=ORIGIN,YEAR=YEAR)

MCAOPTION=""
SCRIPT= "makeShapeCardsHgg.py"

OPTIONS="{OPTIONS} -L 'ttH-multilepton/functionsTTH.cc' -L 'vbfH-diphoton/functionsH2G.cc'  ".format(OPTIONS=OPTIONS)
CATPOSTFIX=""

FITVAR="dipho_mass"
VARBINS='['+','.join([str(b) for b in xrange(100,181)])+']'
CATFUNCTION_2G="dnn_catIndex(vbfDNN_pbkg,vbfDNN_pbsm,D0minus)"
MCASUFFIX="-"+YEAR

DOFILE = ""

OPT_VBF='{T2G} {OPTIONS} -W "weight"'.format(T2G=T2G, OPTIONS=OPTIONS)
CATPOSTFIX=""

ncats = 8
CATBINS="["+",".join([str(i+0.5) for i in xrange(ncats+1)])+"]"
NAMES  = ','.join( 'RECO_%s_%s_Tag0'%(x,y) for x in 'DCP0,DCP1'.split(',') for y in 'Bsm0,Bsm1,Bsm2'.split(','))
NAMES += ','+','.join( 'RECO_%s_Tag1'%(x) for x in 'DCP0,DCP1'.split(','))

procs = {
    'ggh_120_13TeV' : 'ggH120',
    'ggh_125_13TeV' : 'ggH125',
    'ggh_130_13TeV' : 'ggH130',
    'vbfh_120_13TeV': 'vbfH120',
    'vbfh_125_13TeV': 'vbfH125',
    'vbfh_130_13TeV': 'vbfH130',
    'tth_120_13TeV' : 'ttH120',
    'tth_125_13TeV' : 'ttH125',
    'tth_130_13TeV' : 'ttH130',
    'wzh_120_13TeV' : 'wzH120',
    'wzh_125_13TeV' : 'wzH125',
    'wzh_130_13TeV' : 'wzH130',

    'vbfh_ALT_125_13TeV': 'vbfALT0MH125',

    'Data_13TeV'    : 'data',
}

for procid,proc in procs.iteritems():
    ASIMOV = ' --asimov s ' if proc!='data' else ' '
    TORUN='''python {SCRIPT} {DOFILE} vbfH-diphoton/mca-vbfhgg{MCASUFFIX}{MCAOPTION}.txt vbfH-diphoton/vbfhgg.txt "{FITVAR}" "{VARBINS}" {PROC} {OPT_VBF} --binname {PROCID} -p {PROC} {ASIMOV} --categorize "{CATFUNCTION_2G}" "{CATBINS}" "{NAMES}" '''.format(SCRIPT=SCRIPT, DOFILE=DOFILE, MCASUFFIX=MCASUFFIX, MCAOPTION=MCAOPTION, FITVAR=FITVAR, VARBINS=VARBINS, CATFUNCTION_2G=CATFUNCTION_2G, CATBINS=CATBINS, OPT_VBF=OPT_VBF,NAMES=NAMES,PROC=proc,PROCID=procid,ASIMOV=ASIMOV)
    print submit.format(command=TORUN)
    print "\n"
