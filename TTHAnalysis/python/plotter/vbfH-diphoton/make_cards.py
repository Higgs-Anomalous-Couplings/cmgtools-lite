import os, sys
nCores=8
#submit = '''sbatch -c %d -p short  --wrap '{command}' '''%nCores
submit = '{command}' 


ORIGIN="/eos/cms/store/group/dpg_ecal/comm_ecal/localreco/vbfhgg/HiggsCouplings/Trees_161221/2017ULReReco/"; 

if len(sys.argv) < 3: 
    print 'Syntax is %s [outputdir] [year] [other]'%sys.argv[0]
    raise RuntimeError 
OUTNAME=sys.argv[1]
YEAR=sys.argv[2]
OTHER=sys.argv[3:] if len(sys.argv) > 3 else ''

if   YEAR=='2016': LUMI="35.9"
elif YEAR=='2017': LUMI="41.4"
elif YEAR=='2018': LUMI="59.7"
elif YEAR=='Run2': LUMI="138.0"
#elif YEAR in '2018': LUMI="137.0" # using 2018 MC/data as proxy for tot Run 2
else:
    raise RuntimeError("Wrong year %s"%YEAR)


#print "Normalizing to {LUMI}/fb".format(LUMI=LUMI);
OPTIONS=" -j {J} -l {LUMI} -f ".format(LUMI=LUMI,J=nCores)
os.system("test -d cards/{OUTNAME} || mkdir -p cards/{OUTNAME}".format(OUTNAME=OUTNAME))
OPTIONS="{OPTIONS} --od cards/{OUTNAME} ".format(OPTIONS=OPTIONS, OUTNAME=OUTNAME)
T2G="-P '{ORIGIN}' -F 'Friends' '{{P}}/friends/{{cname}}_Friend.root' ".format(ORIGIN=ORIGIN)

MCAOPTION=""
SCRIPT= "makeShapeCardsHgg.py"

OPTIONS="{OPTIONS} -L 'ttH-multilepton/functionsTTH.cc' -L 'vbfH-diphoton/functionsH2G.cc'  ".format(OPTIONS=OPTIONS)
CATPOSTFIX=""

FITVAR="M2G"
VARBINS='['+','.join([str(b) for b in xrange(100,181)])+']'
CATFUNCTION_2G="mela_catIndex( dijet_Mjj, dipho_mva, dijet_mva_prob_VBF, dijet_mva_prob_ggH, D0minus )"
MCASUFFIX=""

DOFILE = ""

OPT_VBF='{T2G} {OPTIONS} -W "weight"'.format(T2G=T2G, OPTIONS=OPTIONS)
CATPOSTFIX=""

ncats = 18
CATBINS="["+",".join([str(i+0.5) for i in xrange(ncats+1)])+"]"
NAMES  = ','.join( 'RECO_%s_%s_%s'%(x,y,z) for x in 'MJJ_250_350,MJJ_350_700,MJJ_GE700'.split(',') for y in 'DCP_0,DCP_1,DCP_2'.split(',') for z in 'Tag0,Tag1'.split(','))

procs = {
    'vbf_125_13TeV': 'vbfH',
    'ggh_125_13TeV': 'ggH'
}

for procid,proc in procs.iteritems():
    ASIMOV = ' --asimov s ' if proc!='data' else ' '
    TORUN='''python {SCRIPT} {DOFILE} vbfH-diphoton/mca-vbfhgg{MCASUFFIX}{MCAOPTION}.txt vbfH-diphoton/vbfhgg.txt "{FITVAR}" "{VARBINS}" {PROC} {OPT_VBF} --binname {PROCID} -p {PROC} {ASIMOV} --categorize "{CATFUNCTION_2G}" "{CATBINS}" "{NAMES}" '''.format(SCRIPT=SCRIPT, DOFILE=DOFILE, MCASUFFIX=MCASUFFIX, MCAOPTION=MCAOPTION, FITVAR=FITVAR, VARBINS=VARBINS, CATFUNCTION_2G=CATFUNCTION_2G, CATBINS=CATBINS, OPT_VBF=OPT_VBF,YEAR=YEAR,NAMES=NAMES,PROC=proc,PROCID=procid,ASIMOV=ASIMOV)
    print submit.format(command=TORUN)
            
