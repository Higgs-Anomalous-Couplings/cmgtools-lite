import ROOT as r 
import numpy as np
import pickle,math,os
from keras.models import Sequential
from keras.layers import Dense, Dropout

import tensorflow as tf
from keras.backend import tensorflow_backend as K

from keras import optimizers


def getCompiledModelA(nvars,nnodes):
    # optimal so far ( 0.980, 0.966)
    model = Sequential()
    model.add(Dense(10,input_dim=nvars, activation='relu'))
    model.add(Dense(6, activation='relu'))
    model.add(Dense(nnodes, activation='softmax'))
    adam = optimizers.adam(lr=1e-4) 
    model.compile(loss='categorical_crossentropy', optimizer=adam, metrics=['accuracy','categorical_crossentropy'])
    return model

def getCompiledModelB(nvars,nnodes):
    model = Sequential()
    model.add(Dense(10,input_dim=nvars, activation='relu'))
    model.add(Dropout(0.4))
    model.add(Dense(6, activation='relu'))
    model.add(Dropout(0.4))
    model.add(Dense(6, activation='relu'))
    model.add(Dense(nnodes, activation='softmax'))
    adam = optimizers.adam(lr=1e-4) 
    model.compile(loss='categorical_crossentropy', optimizer=adam, metrics=['accuracy','categorical_crossentropy'])
    return model

def freeze_session(session, keep_var_names=None, output_names=None, clear_devices=True):
    """
    Freezes the state of a session into a pruned computation graph.

    Creates a new computation graph where variable nodes are replaced by
    constants taking their current value in the session. The new graph will be
    pruned so subgraphs that are not necessary to compute the requested
    outputs are removed.
    @param session The TensorFlow session to be frozen.
    @param keep_var_names A list of variable names that should not be frozen,
                          or None to freeze all the variables in the graph.
    @param output_names Names of the relevant graph outputs.
    @param clear_devices Remove the device directives from the graph for better portability.
    @return The frozen graph definition.
    """
    from tensorflow.python.framework.graph_util import convert_variables_to_constants
    graph = session.graph
    with graph.as_default():
        freeze_var_names = list(set(v.op.name for v in tf.global_variables()).difference(keep_var_names or []))
        output_names = output_names or []
        output_names += [v.op.name for v in tf.global_variables()]
        # Graph -> GraphDef ProtoBuf
        input_graph_def = graph.as_graph_def()
        if clear_devices:
            for node in input_graph_def.node:
                node.device = ""
        frozen_graph = convert_variables_to_constants(session, input_graph_def,
                                                      output_names, freeze_var_names)
        return frozen_graph

def saveTF1p6Model(protobuffer_name):
    frozen_graph = freeze_session(K.get_session(),
                                  output_names=[out.op.name for out in model.outputs])
    tf.train.write_graph(frozen_graph, os.getcwd(), protobuffer_name, as_text=False)
    print ("Saved frozen model in ",protobuffer_name)
    
    
if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options]")
    parser.add_option("-i", "--infile", dest="infile", type="string", default="vars.pkl", help="Input pickle file (default: vars.pkl)");
    parser.add_option("-o", "--outfile", dest="outfile", type="string", default="trained_model", help="Output pb and h5 fil;e (default: trained_model)");
    (options, args) = parser.parse_args()

    nvars = 10

    data = pickle.load( open(options.infile,'rb'))
    sums = np.sum(data['train_y'],axis=0)
    print(sums)

    sig = sums[0]
    bkg = sums[1] + sums[2]

    class_weight = { 0 : float((sig+bkg)/sig),
                     1 : float((sig+bkg)/bkg),
                     2 : float((sig+bkg)/bkg)}

    print ('weights will be', class_weight)

    with tf.Session(config=tf.ConfigProto(
            intra_op_parallelism_threads=50,
            inter_op_parallelism_threads=50)) as sess:
        K.set_session(sess)
        
        nnodes = 3
        model = getCompiledModelA(nvars,nnodes)
        #model = getCompiledModelB(nvars,nnodes)

        history = model.fit( data['train_x'], data['train_y'], epochs=5, batch_size=200, validation_data=(data['test_x'], data['test_y']), class_weight=class_weight)

        modelname = os.getcwd()+'/'+os.path.basename(options.outfile).split('.')[0]
        # keras model (H5)
        model.save(modelname+'.h5')
        # tf model (PB)
        saveTF1p6Model(modelname+'.pb')
        
        pickle_out = open(modelname+'.pkl','wb')
        pickle.dump( history.history, pickle_out)
        pickle_out.close()
        
