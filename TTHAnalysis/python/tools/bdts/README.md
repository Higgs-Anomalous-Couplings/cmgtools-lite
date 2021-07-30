# Short recipe

#### Set up CMSSW 
```
cmsrel CMSSW_10_4_0 (only the first time, not if you have already the directory)
cd CMSSW_10_4_0/src
cmsenv
```

#### To convert the trees in the pickle format needed by Tensorflow
This is needed only when you want to modify the variables (add some, remove others) or change the selection of the events going to the training.
```
python makeTrainDataset.py -o vars_onlysm.pkl
```

#### To make the training 
This is the main script, where to optimize which variables, the structure of the NN, the training optiminzation
```
python trainNet.py -i vars_onlysm.pkl -o trained_model_onlysm
```

#### Plotting
This is a simple plotting script that plots the input variables of the DNN, the multiclass DNN output and the ROC curve.
This is very rough (without labels, etc)
```
python plotting.py -v vars_onlysm.pkl -m trained_model_onlysm.h5 --pdir plotdir/
```
