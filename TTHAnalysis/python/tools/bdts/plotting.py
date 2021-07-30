from __future__ import print_function

import numpy as np
import pickle,matplotlib,os
matplotlib.use('ps')
import matplotlib.pyplot as plt

from keras.models import load_model
from sklearn.metrics import roc_curve,auc

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options]")
    parser.add_option("-v", "--varfile", dest="varfile", type="string", default="vars.pkl", help="Input pickle file (default: vars.pkl)");
    parser.add_option("-m", "--modelfile", dest="modelfile", type="string", default="model.h5", help="Input model file (default: model.h5)");
    parser.add_option("--pdir", "--print-dir", dest="printDir", type="string", default="plots", help="print out plots in this directory");
    (options, args) = parser.parse_args()

    outname = options.printDir
    if not os.path.exists(outname):
        os.system("mkdir -p "+outname)
        if os.path.exists("/afs/cern.ch"): os.system("cp /afs/cern.ch/user/g/gpetrucc/php/index.php "+outname)
    print ("Will save plots to ",outname)

    if not os.path.exists(options.varfile) or not os.path.exists(options.modelfile):
        raise RuntimeError, "Either input pkl file or h5 model file not accessible. Exiting. "
    data = pickle.load(open(options.varfile,'rb'))
    model = load_model(options.modelfile)

    prediction = model.predict(data['test_x'])
    x = data['test_x']
    y = np.argmax(data['test_y'], axis=1)

    classifier = np.sum( prediction[:,[0,1]], axis=1)

    fpr_keras, tpr_keras, thresholds_keras = roc_curve(y<=1, classifier)
    auc_keras = auc(fpr_keras, tpr_keras)

    plt.plot([0, 1], [0, 1], 'k--')
    plt.plot(fpr_keras, tpr_keras, label='Keras (area = {:.3f})'.format(auc_keras))
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('ROC curve')
    plt.legend(loc='best')

    plt.show()
    for ext in ['png','pdf']:
        plt.savefig(outname+'/roc.'+ext)

    nvars = 10
    for var in range(nvars):
        together = np.dstack( (x[:,var], y ) )[0]
        
        class1 =  (together[together[:,1] == 0]) [:,0]
        class2 =  (together[together[:,1] == 1]) [:,0]
        class3 =  (together[together[:,1] == 2]) [:,0]

        bins = 20
        plt.clf()
        plt.hist(class1, bins, alpha=0.5,density=True, label='vbfH')
        plt.hist(class2, bins, alpha=0.5,density=True, label='ggH')
        plt.hist(class3, bins, alpha=0.5,density=True, label='other')
        plt.legend(loc='upper right')
        plt.show()
    
        for ext in ['png','pdf']:
            plt.savefig('{o}/input_{var}.{ext}'.format(o=outname,var=var,ext=ext))

    nnodes = 3
    for node in range(nnodes):
        together = np.dstack( (prediction[:,node], y ) )[0]
        
        class1 =  (together[together[:,1] == 0]) [:,0]
        class2 =  (together[together[:,1] == 1]) [:,0]
        class3 =  (together[together[:,1] == 2]) [:,0]
        
        bins = 20
        plt.clf()
        plt.hist(class1, bins, alpha=0.5,density=True, label='vbfH')
        plt.hist(class2, bins, alpha=0.5,density=True, label='ggH')
        plt.hist(class3, bins, alpha=0.5,density=True, label='other')
        plt.legend(loc='upper right')
        plt.show()

        for ext in ['png','pdf']:
            plt.savefig('{o}/output_{var}.{ext}'.format(o=outname,var=node,ext=ext))

    together = np.dstack( (classifier,y))[0]
    class1 =  (together[together[:,1] == 0]) [:,0]
    class2 =  (together[together[:,1] == 1]) [:,0]
    class3 =  (together[together[:,1] == 2]) [:,0]

    bins = 20
    plt.clf()
    plt.hist(class1, bins, alpha=0.5,density=True, label='vbfH')
    plt.hist(class2, bins, alpha=0.5,density=True, label='ggH')
    plt.hist(class3, bins, alpha=0.5,density=True, label='other')
        
    plt.legend(loc='upper left')
    plt.show()
    
    for ext in ['png','pdf']:
        plt.savefig('{o}/output_combined.{ext}'.format(o=outname,ext=ext))

