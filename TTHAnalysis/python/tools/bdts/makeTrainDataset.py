from __future__ import print_function

import ROOT as r 
import numpy as np
import pickle,math,os,tqdm
from keras.utils import np_utils
from multiprocessing import Pool

testCut   = lambda ev: ev.event%5==0
trainCut  = lambda ev: ev.event%5!=0

commonFeatureList = [
    "leadJPt",
    "leadJEta",
    "subJPt",
    "subJEta",
    "leadGPtOverM",
    "subGPtOverM",
    "mJJ",
    "centrality",
    "dijet_dipho_dphi_trunc",
    "dijet_abs_dEta",
    "dijet_dphi",
    "dijet_minDRJetPho",
    "dipho_PToM",
 ]


features = {

    "leadJPt"          : lambda ev : ev.dijet_LeadJPt if ev.dijet_nj >= 2 else -9,
    "leadJEta"         : lambda ev : abs(ev.dijet_leadEta) if ev.dijet_nj >= 2 else -9,
    "subJPt"           : lambda ev : ev.dijet_SubJPt if ev.dijet_nj >= 2 else -9,
    "subJEta"          : lambda ev : abs(ev.dijet_subleadEta) if ev.dijet_nj >= 2 else -9,
    "leadGPtOverM"     : lambda ev : ev.leadPho_PToM, 
    "subGPtOverM"      : lambda ev : ev.sublPho_PToM,
    "mJJ"              : lambda ev : ev.dijet_Mjj,
    "centrality"       : lambda ev : ev.dijet_centrality_gg,
    "dijet_dipho_dphi_trunc" : lambda ev : ev.dijet_dipho_dphi_trunc if ev.dijet_nj >= 2 else -9,
    "dijet_abs_dEta"   : lambda ev : ev.dijet_abs_dEta if ev.dijet_nj >= 2 else -9,
    "dijet_dphi"       : lambda ev : ev.dijet_dphi if ev.dijet_nj >= 2 else -9,
    "dijet_minDRJetPho": lambda ev : ev.dijet_minDRJetPho if ev.dijet_nj >= 2 else -9,
    "dipho_PToM"       : lambda ev : ev.dipho_PToM,
    }

cuts = {
    'vbfH' : lambda ev : (abs(ev.dipho_leadEta) < 2.5 and abs(ev.dipho_subleadEta) < 2.5 and (abs(ev.dipho_leadEta) < 1.44 or abs(ev.dipho_leadEta) > 1.57) and (abs(ev.dipho_subleadEta) < 1.44 or abs(ev.dipho_subleadEta) > 1.57)) and ev.dipho_mass > 100.0 and ev.dipho_mass < 180.0 and ev.dipho_lead_ptoM > 0.333 and ev.dipho_sublead_ptoM > 0.25 and ev.dipho_leadIDMVA > -0.2 and ev.dipho_subleadIDMVA > -0.2 and ev.dijet_abs_dEta > 0.0 and abs(ev.dijet_leadEta) < 4.7 and abs(ev.dijet_subleadEta) < 4.7 and ev.dijet_minDRJetPho > 0.4 and ev.dijet_LeadJPt > 40. and ev.dijet_SubJPt > 30. and ev.dijet_Mjj > 250.0,

    'vbfL1H' : lambda ev : (abs(ev.dipho_leadEta) < 2.5 and abs(ev.dipho_subleadEta) < 2.5 and (abs(ev.dipho_leadEta) < 1.44 or abs(ev.dipho_leadEta) > 1.57) and (abs(ev.dipho_subleadEta) < 1.44 or abs(ev.dipho_subleadEta) > 1.57)) and ev.dipho_mass > 100.0 and ev.dipho_mass < 180.0 and ev.dipho_lead_ptoM > 0.333 and ev.dipho_sublead_ptoM > 0.25 and ev.dipho_leadIDMVA > -0.2 and ev.dipho_subleadIDMVA > -0.2 and ev.dijet_abs_dEta > 0.0 and abs(ev.dijet_leadEta) < 4.7 and abs(ev.dijet_subleadEta) < 4.7 and ev.dijet_minDRJetPho > 0.4 and ev.dijet_LeadJPt > 40. and ev.dijet_SubJPt > 30. and ev.dijet_Mjj > 250.0,
    
    'ggH' : lambda ev : (abs(ev.dipho_leadEta) < 2.5 and abs(ev.dipho_subleadEta) < 2.5 and (abs(ev.dipho_leadEta) < 1.44 or abs(ev.dipho_leadEta) > 1.57) and (abs(ev.dipho_subleadEta) < 1.44 or abs(ev.dipho_subleadEta) > 1.57)) and ev.dipho_mass > 100.0 and ev.dipho_mass < 180.0 and ev.dipho_lead_ptoM > 0.333 and ev.dipho_sublead_ptoM > 0.25 and ev.dipho_leadIDMVA > -0.2 and ev.dipho_subleadIDMVA > -0.2 and ev.dijet_abs_dEta > 0.0 and abs(ev.dijet_leadEta) < 4.7 and abs(ev.dijet_subleadEta) < 4.7 and ev.dijet_minDRJetPho > 0.4 and ev.dijet_LeadJPt > 40. and ev.dijet_SubJPt > 30. and ev.dijet_Mjj > 250.0,

    'other' : lambda ev: (abs(ev.dipho_leadEta) < 2.5 and abs(ev.dipho_subleadEta) < 2.5 and (abs(ev.dipho_leadEta) < 1.44 or abs(ev.dipho_leadEta) > 1.57) and (abs(ev.dipho_subleadEta) < 1.44 or abs(ev.dipho_subleadEta) > 1.57)) and ev.dipho_mass > 100.0 and ev.dipho_mass < 180.0 and ev.dipho_lead_ptoM > 0.333 and ev.dipho_sublead_ptoM > 0.25 and ev.dipho_leadIDMVA > -0.2 and ev.dipho_subleadIDMVA > -0.2 and ev.dijet_abs_dEta > 0.0 and abs(ev.dijet_leadEta) < 4.7 and abs(ev.dijet_subleadEta) < 4.7 and ev.dijet_minDRJetPho > 0.4 and ev.dijet_LeadJPt > 40. and ev.dijet_SubJPt > 30. and ev.dijet_Mjj > 250.0,
}

classes = {
    'vbfH'      : { 'cut': cuts['vbfH'],  'lst_train' : [], 'lst_test' : [] , 'lst_y_train' : [], 'lst_y_test' : [] },
    'ggH'       : { 'cut': cuts['ggH'],   'lst_train' : [], 'lst_test' : [] , 'lst_y_train' : [], 'lst_y_test' : [] },
    'other'     : { 'cut': cuts['other'], 'lst_train' : [], 'lst_test' : [] , 'lst_y_train' : [], 'lst_y_test' : [] },
}

vbfL1Hclass = {
    'vbfL1H'    : { 'cut': cuts['vbfL1H'],'lst_train' : [], 'lst_test' : [] , 'lst_y_train' : [], 'lst_y_test' : [] },
}

sampleDir='/eos/cms/store/group/phys_higgs/emanuele/vbfhgg_ac/VBFHiggs_UL2017_09July2021/'

vbfSMSamples   = []
vbfCPSamples   = []
vbfL1Samples = []
ggHSamples   = []
otherSamples = []

vbfSMSamples.extend( ['output_VBFHToGG_M-125_TuneCP5_13TeV-powheg-pythia8.root', 'output_VBFHToGG_M125_TuneCP5_13TeV-amcatnlo-pythia8.root'] )
vbfCPSamples.extend( ['output_VBFHiggs0Mf05ph0ToGG_M125_TuneCP5_13TeV-JHUGenV7011-pythia8.root', 'output_VBFHiggs0MToGG_M125_TuneCP5_13TeV-JHUGenV7011-pythia8.root', 'output_VBFHiggs0PHf05ph0ToGG_M125_TuneCP5_13TeV-JHUGenV7011-pythia8.root', 'output_VBFHiggs0PHToGG_M125_TuneCP5_13TeV-JHUGenV7011-pythia8.root'] )
vbfL1Samples.extend( ['output_VBFHiggs0L1ToGG_M125_TuneCP5_13TeV-JHUGenV7011-pythia8.root', 'output_VBFHiggs0L1f05ph0ToGG_M125_TuneCP5_13TeV-JHUGenV7011-pythia8.root', 'output_VBFHiggs0L1Zgf05ph0ToGG_M125_TuneCP5_13TeV-JHUGenV7011-pythia8.root'] )
ggHSamples.extend( ['output_GluGluHToGG_M-125_TuneCP5_13TeV-powheg-pythia8.root', 'output_GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8.root'] )
otherSamples.extend( ['output_DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa.root','output_DiPhotonJetsBox_M40_80-sherpa.root'] )


def toNumpy(featureList,maxEntries,task):
    print('starting', task)
    fil, typs = task
    path = os.path.dirname(fil)
    fname = os.path.basename(fil)
    print('List of features for', featureList + eval('featureList'))
    print("file = ",fil,"  typs = ",typs)
    tdir = 'vbfTagDumper/trees/'
    treenames = {'vbfH'  : 'vbf_125_13TeV_GeneralDipho',
                 'vbfL1H': 'vbf_125_13TeV_GeneralDipho',
                 'ggH'   : 'ggh_125_13TeV_GeneralDipho',
                 'other' : 'dipho_13TeV_GeneralDipho'}
    tfile = r.TFile(fil); ttree = tfile.Get('vbfTagDumper/trees/{tname}'.format(tname=treenames[typs[0]]))
    print ("number of entries for ",treenames[typs[0]]," is ",ttree.GetEntries())
    results = {}
    for ty in typs: 
        results[ty + '_test']  = []
        results[ty + '_train'] = []

    print("start looping on file ",fil)
    for iev,ev in enumerate(ttree):
        if iev%1000==0: print('Processing event ',iev,'...')
        if iev>maxEntries: break
        tstr = 'test' if testCut(ev) else 'train'
        for ty in typs:
            if classes[ty]['cut'](ev):
                results[ty+'_'+tstr].append([ features[s](ev) for s in (featureList) ])
    tfile.Close()
    print('finishing', task)
    return results




if __name__ == "__main__":

    from functools import partial

    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options]")
    parser.add_option("--max-entries",    dest="maxEntries", default=1000000000, type="int", help="Max entries to process in each tree")
    parser.add_option("-o",   "--outfile",  dest="outfile", type="string", default="vars.pkl", help="Output pickle file (default: vars.pkl)");
    parser.add_option("--cp", "--add-cp-samples", dest="addcpsamples", action="store_true", default=False, help="Add also the alternative CP samples to the signal");
    parser.add_option("--l1", "--add-l1-samples", dest="addl1samples", action="store_true", default=False, help="Add also the alternative Lambda1 samples to the signal");
    parser.add_option("--four-classes", dest="fourclasses", action="store_true", default=False, help="Make a custom class for L1 samples");
    (options, args) = parser.parse_args()

    tasks = []
    print('Setting up the tasks')
    vbfSamples = vbfSMSamples
    if options.addcpsamples: 
        vbfSamples = vbfSamples+vbfCPSamples
    if options.addl1samples and not options.fourclasses:
        vbfSamples = vbfSamples+vbfL1Samples
    for samp in vbfSamples:
        tasks.append( (sampleDir+'/'+samp, ['vbfH']) )
    if options.fourclasses:
        for samp in vbfL1Samples:
            tasks.append( (sampleDir+'/'+samp, ['vbfL1H']) )
    for samp in ggHSamples:
        tasks.append( (sampleDir+'/'+samp, ['ggH']) )
    for samp in otherSamples:
        tasks.append( (sampleDir+'/'+samp, ['other']) )
    
    print('Numpy conversion. It will take time...')
    print("max entries = ",options.maxEntries)

    if options.fourclasses:
        classes.update(vbfL1Hclass)

    featureList = commonFeatureList
    for cl,vals in classes.iteritems():
        vals['cut'] = cuts[cl]

    ## lxplus seems to have 10 cores/each
    p =  Pool(min(30,len(vbfSamples+vbfL1Samples+ggHSamples+otherSamples)))
    func = partial(toNumpy,featureList,options.maxEntries)
    results = list(tqdm.tqdm(p.imap(func, tasks), total=len(tasks)))

    print('Now putting everything together')

    ## for later training, is useful to have types sorted in this way
    types = ['vbfH', 'vbfL1H', 'ggH', 'other'] if options.fourclasses else ['vbfH', 'ggH', 'other']
    for result in results: 
        for ty in types:
            if ty+'_train' in result:
                classes[ty]['lst_train'].extend( result[ty+'_train'])
                classes[ty]['lst_test' ].extend( result[ty+'_test'])

            
    print('Setting the indices')
    toDump = {} 
    for i, ty in enumerate(types):
        classes[ty]['lst_train'  ] = np.asarray(classes[ty]['lst_train'])
        classes[ty]['lst_y_train'] = i*np.ones((classes[ty]['lst_train'].shape[0],1))
        classes[ty]['lst_test'   ] = np.asarray(classes[ty]['lst_test'])
        classes[ty]['lst_y_test' ] = i*np.ones((classes[ty]['lst_test'].shape[0],1))

    train_x = np.concatenate( tuple( [classes[ty]['lst_train'] for ty in types] ), axis=0)
    train_y = np_utils.to_categorical( np.concatenate( tuple( [classes[ty]['lst_y_train'] for ty in types] ), axis=0), len(classes))
    test_x = np.concatenate( tuple( [classes[ty]['lst_test'] for ty in types] ), axis=0)
    test_y = np_utils.to_categorical( np.concatenate( tuple( [classes[ty]['lst_y_test'] for ty in types] ), axis=0), len(classes))
    toDump['train_x'] = train_x
    toDump['train_y'] = train_y
    toDump['test_x' ] = test_x
    toDump['test_y' ] = test_y

    ### dump to file
    print ('dump to ',options.outfile,' now...')
    pickle_out = open(options.outfile,'wb')
    pickle.dump( toDump, pickle_out)
    pickle_out.close()
