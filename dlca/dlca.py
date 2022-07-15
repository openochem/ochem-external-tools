
import csv
from collections import OrderedDict
import numpy as np
import math
import sys
import re
import os
import configparser
import random
import h5py
import tarfile
import shutil
import pickle

from sklearn.preprocessing import StandardScaler

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolSurf import *
from rdkit.Chem.rdMolDescriptors import *
from rdkit.Chem.Descriptors import *
from rdkit import Avalon
from rdkit.Avalon.pyAvalonTools import *
import rdkit.Chem.AtomPairs.Pairs as Pairs

version = 1;
print("Version: ", version);

if(len (sys.argv) != 2):
    print("Usage: ", sys.argv[0], "config.cfg");
    sys.exit(0);

print("Load config file: ", sys.argv[1]);

config = configparser.ConfigParser();
config.read(sys.argv[1]);

def getConfig(section, attribute, default=""):
    try:
        return config[section][attribute];
    except:
        return default;

TRAIN = getConfig("Task","train_mode");
MODEL_FILE = getConfig("Task","model_file");
TRAIN_FILE = getConfig("Task","train_data_file");
APPLY_FILE = getConfig("Task","apply_data_file", "train.csv");
RESULT_FILE = getConfig("Task","result_file", "results.csv");
NUM_EPOCHS = int(getConfig("Details","n_epochs", "100"));
BATCH_SIZE = int(getConfig("Details","batch_size", "32"));
SEED = int(getConfig("Details","seed", "657488"));
DEVICE = getConfig("Details","gpu");

np.random.seed(SEED);
random.seed(SEED);

os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID";
os.environ["CUDA_VISIBLE_DEVICES"] = DEVICE;

import tensorflow as tf
from tensorflow.keras import layers
from tensorflow.keras import backend as K
from tensorflow.keras.callbacks import ModelCheckpoint

config = tf.compat.v1.ConfigProto(allow_soft_placement=True, log_device_placement=False);
config.gpu_options.allow_growth = True;
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR);
tf.compat.v1.keras.backend.set_session(tf.compat.v1.Session(config=config));
tf.compat.v1.disable_eager_execution();

try: tf.random.set_seed(SEED);
except:
    try: tf.random.set_random(SEED);
    except: print("Cannot set seed to TensorFlow.");

#Calculated from the whole 80081 endpoints from the dataset.csv file.
chars = "#%()+-./0123456789=ABCFGHIKLMNOPRSTUVZ[\]abcdefghilnorstu";
char_to_ix = { ch:i for i,ch in enumerate(chars) }
ix_to_char = { i:ch for i,ch in enumerate(chars) }

def calcDescriptors(m):
   #returns the descriptor set as in DLCA.

   #Parameters:
   #     m - RDKit Mol object

   #Returns:
   #     [0:119] - rdkit descriptors
   #     [119:1143] - morgan fingerprint bits
   #     [1143:2167] - avalon fingerprint bits
   #     [2167:] on keys in pairs fingerprint

   Chem.AssignStereochemistry(m, flagPossibleStereoCenters=True, force=True);
   crippen = CalcCrippenDescriptors(m);

   nhs = 0;
   for atom in m.GetAtoms():
      nhs += atom.GetTotalNumHs();

   #rdkit descriptors
   descrs = [ crippen[0],
              crippen[1],
              CalcLabuteASA(m),
              CalcTPSA(m),
              MolWt(m),
              CalcExactMolWt(m),
              CalcNumLipinskiHBA(m),
              CalcNumLipinskiHBD(m),
              CalcNumRotatableBonds(m),
              CalcNumHBD(m),
              CalcNumHBA(m),
              CalcNumAmideBonds(m),
              CalcNumHeteroatoms(m),
              HeavyAtomCount(m),
              m.GetNumAtoms() + nhs,
              CalcNumAtomStereoCenters(m),
              CalcNumUnspecifiedAtomStereoCenters(m),
              CalcNumRings(m),
              CalcNumAromaticRings(m),
              CalcNumSaturatedRings(m),
              CalcNumAliphaticRings(m),
              CalcNumAromaticHeterocycles(m),
              CalcNumSaturatedHeterocycles(m),
              CalcNumAliphaticHeterocycles(m),
              CalcNumAromaticCarbocycles(m),
              CalcNumSaturatedCarbocycles(m),
              CalcNumAliphaticCarbocycles(m),
              CalcFractionCSP3(m),
              CalcChi0v(m), CalcChi1v(m), CalcChi2v(m), CalcChi3v(m), CalcChi4v(m),
              CalcChi1n(m), CalcChi2n(m), CalcChi3n(m), CalcChi4n(m),
              CalcHallKierAlpha(m),
              CalcKappa1(m), CalcKappa2(m), CalcKappa3(m),
              SlogP_VSA1(m), SlogP_VSA2(m), SlogP_VSA3(m), SlogP_VSA4(m),
              SlogP_VSA5(m), SlogP_VSA6(m), SlogP_VSA7(m), SlogP_VSA8(m),
              SlogP_VSA9(m), SlogP_VSA10(m), SlogP_VSA11(m), SlogP_VSA12(m),
              SMR_VSA1(m), SMR_VSA2(m), SMR_VSA3(m), SMR_VSA4(m), SMR_VSA5(m),
              SMR_VSA6(m), SMR_VSA7(m), SMR_VSA8(m), SMR_VSA9(m), SMR_VSA10(m),
              PEOE_VSA1(m), PEOE_VSA2(m), PEOE_VSA3(m), PEOE_VSA4(m), PEOE_VSA5(m),
              PEOE_VSA6(m), PEOE_VSA7(m), PEOE_VSA8(m), PEOE_VSA9(m), PEOE_VSA10(m),
              PEOE_VSA11(m), PEOE_VSA12(m), PEOE_VSA13(m), PEOE_VSA14(m)];
   descrs += rdMolDescriptors.MQNs_(m);

   #Morgan
   fp = [0]*1024;
   for i in AllChem.GetMorganFingerprintAsBitVect(m, 2, 1024).GetOnBits():
      fp[i] = 1;
   descrs += fp;

   #avalon
   fp = [0]*1024;
   for i in GetAvalonFP(m, 1024).GetOnBits():
       fp[i] = 1;
   descrs += fp;

   #atom pairs
   for i in Pairs.GetAtomPairFingerprintAsBitVect(m).GetOnBits():
      descrs.append(i);
   return descrs;

def CleanSmileString(smi):
    #replace as in the original code
    smi = smi.replace("[nH]", "A");
    smi = smi.replace("Cl", "L");
    smi = smi.replace("Br", "R");
    smi = smi.replace("[C@]", "C");
    smi = smi.replace("[C@@]", "C");
    smi = smi.replace("[C@@H]", "C");
    return smi;

def readDataSet(fname):
    #returns descriptors for ML

    #returns
    #   numpy descriptor matrix (normalized)
    #   y values (not normalized)
    #   info with additional parameters on pairs keys, means, stds.

    ds = [];
    props = [];
    print("Loading data.");

    line = 0;
    for row in csv.reader(open(fname, "r")):
        line += 1;
        if line == 1:
           props = row[1:];
           continue;
        try:
           m = Chem.MolFromSmiles(row[0]);
           ds.append( [CleanSmileString(row[0]), calcDescriptors(m), row[1:]])
        except Exception as e:
           print(e);

    #shuffle
    random.shuffle(ds);

    print("Loaded ", len(ds), " records.");
    print("Number of properties: ", len(ds[0][2]));

    #Generating a dictionary for Pairs fingerprint.
    keys = set();
    for d in ds:
        keys |= set(d[1][2167:]);

    counts = { k:0 for k in keys};
    for d in ds:
        for k in d[1][2167:]:
            counts[k] += 1;

    #Select Top-1024 (the most populated pairs).
    keys = list(OrderedDict(sorted(counts.items(), key=lambda x: x[1], reverse = True)))[:1024];

    #The authors did not scale or did other preprocessing for properties.
    y = np.zeros( (len(ds), len(ds[0][2])), dtype=np.float32);
    y[:] = np.nan;  #The loss functions explicitely wants nan, where there is no value.

    #Collect all the descriptors in a one big matrix.
    descrs = np.zeros( (len(ds), 2167 + len(keys)), dtype=np.float32);
    descrs_smiles = np.zeros( (len(ds), 200), dtype=np.int8);

    for i, d in enumerate(ds):
        descrs[i,:2167] = ds[i][1][:2167];
        for k in ds[i][1][2167:]:
            try:
               pos = keys.index(k);
               descrs[i, 2167 + pos] = 1;
            except:
               pass;
        for j, p in enumerate(d[2]):
           if p: y[i, j] = float(p);
        for j, c in enumerate(ds[i][0]):
           if j < 200: # IVT fix, crashes for 200
              descrs_smiles[i, j] = char_to_ix[c];

    #The list with all parameters of the data.
    # [0] Pairs uniq_keys
    # [1] Scaler parameters for RDKit descriptors.

    info = [keys, y.shape[1]];

    #DLCA authors use StandardsScaler, so do we.
    scaler = StandardScaler();
    descrs[:, :119] = scaler.fit_transform( descrs [:, :119]);
    info.append( [scaler.mean_, scaler.scale_] );

    info.append( props );

    return descrs, descrs_smiles, y, info;

def buildModel(info):

    inputs = [];
    N = info[1];

    #RDKit
    desc = layers.Input(shape=(119,), name='desc_rdkit')
    x = layers.Dense(8000, activation='relu')(desc)
    x = layers.Dropout(0.3)(x)
    x = layers.BatchNormalization()(x);
    x = layers.Dense(2000, activation='relu')(x)
    x = layers.Dense(1000, activation='relu')(x)
    x = layers.BatchNormalization()(x)
    x = layers.Dense(700, activation='relu')(x)
    out1 = layers.Dense(N, activation='linear', name='out_rdkit')(x)
    inputs.append(desc);
    
    #Morgan
    desc = layers.Input(shape=(1024,), name='desc_morgan')
    x = layers.Dense(8000, activation='relu')(desc)
    x = layers.Dense(2000, activation='relu')(x)
    x = layers.Dense(1000, activation='relu')(x)
    x = layers.Dense(700, activation='relu')(x)
    out2 = layers.Dense(N, activation='linear', name='out_morgan')(x)
    inputs.append(desc);

    #Avalon
    desc = layers.Input(shape=(1024,), name='desc_avalon')
    x = layers.Dense(8000, activation='relu')(desc)
    x = layers.Dense(2000, activation='relu')(x)
    x = layers.Dense(1000, activation='relu')(x)
    x = layers.Dense(700, activation='relu')(x)
    out3 = layers.Dense(N, activation='linear', name='out_avalon')(x)
    inputs.append(desc);

    #Pairs
    desc = layers.Input(shape=(len(info[0]),), name='desc_pairs')
    x = layers.Dense(8000, activation='relu')(desc)
    x = layers.Dense(2000, activation='relu')(x)
    x = layers.Dense(1000, activation='relu')(x)
    x = layers.Dense(700, activation='relu')(x)
    out4 = layers.Dense(N, activation='linear', name='out_pairs')(x)
    inputs.append(desc);

    desc = layers.Input(shape=(200,), name='smiles')
    x = layers.Embedding(len(chars) + 1, 128, input_length=200)(desc)
    x = layers.Conv1D(256,16,padding='valid',activation='relu',strides=1)(x)
    x = layers.GlobalMaxPooling1D()(x)
    x = layers.Dense(200, activation='relu')(x)
    out5 = layers.Dense(N, activation='linear',name='out_smiles')(x)
    inputs.append(desc);

    #Averaging...
    out6 = layers.Average()([out1, out2, out3, out4, out5])
    model = tf.keras.Model(inputs=inputs, outputs=[out1, out2, out3, out4, out5, out6]);

    #Missing values are marked with NaN values.
    def mse(y_true, y_pred):
         y_true = tf.where(tf.math.is_nan(y_true), y_pred, y_true)
         cost = tf.abs(y_pred - y_true)
         return K.sum(cost, axis=-1)

    model.compile (optimizer = 'adam', loss = mse);
    #model.summary();

    return model;

if __name__ == "__main__":

    device_str = "GPU" + str(DEVICE);

    if TRAIN == "True":
        print("Analyze training file...");
        ds, descrs_smiles, y, info = readDataSet(TRAIN_FILE);
        model = buildModel(info);

        N = len(ds);
        n_train = int(math.floor(N * 0.85));
        n_test = N - n_train;

        train_x = ds[:n_train];
        train_smiles = descrs_smiles[:n_train];
        train_y = y[:n_train];

        test_x = ds[n_train:];
        test_smiles = descrs_smiles[n_train:];
        test_y = y[n_train:];

        def train_generator():
            while True:
               for i in range(0, n_train, BATCH_SIZE):
                  y = train_y[i:i+BATCH_SIZE,:];
                  r =  [ train_x[i:i+BATCH_SIZE, :119],
                         train_x[i:i+BATCH_SIZE, 119:1143],
                         train_x[i:i+BATCH_SIZE, 1143:2167],
                         train_x[i:i+BATCH_SIZE, 2167:2167 + len(info[0])],
                         train_smiles[i: i+BATCH_SIZE, :]], [y, y, y, y, y, y];
                  yield r;

        def test_generator():
            while True:
               for i in range(0, n_test, BATCH_SIZE):
                  y = test_y[i:i+BATCH_SIZE,:];
                  yield  [ test_x[i:i+BATCH_SIZE, :119],
                           test_x[i:i+BATCH_SIZE, 119:1143],
                           test_x[i:i+BATCH_SIZE, 1143:2167],
                           test_x[i:i+BATCH_SIZE, 2167:2167 + len(info[0])],
                           test_smiles[i: i+BATCH_SIZE, :]], [y, y, y, y, y, y];

        class MessagerCallback(tf.keras.callbacks.Callback):

           def __init__(self, **kwargs):
              self.early_count = 0;
              self.early_best = 0;

           def on_epoch_end(self, epoch, logs={}):
              print("MESSAGE: train score: {} / validation score: {} / at epoch: {} {} ".format(round(float(logs["loss"]), 7),
                    round(float(logs["val_loss"]), 7), epoch +1, device_str));
              early = float(logs["val_loss"]);
              if(epoch == 0):
                    self.early_best = early;
              else:
                 if early < self.early_best :
                    self.early_count = 0;
                    self.early_best = early;
                 else:
                    self.early_count += 1;
                    if self.early_count > 0.2 * NUM_EPOCHS:
                       self.model.stop_training = True;

              if os.path.exists("stop"):
                 self.model.stop_training = True;
        # end of the callback

        history = model.fit_generator( generator = train_generator(),
                                       steps_per_epoch = n_train // BATCH_SIZE +1 ,
                                       epochs = NUM_EPOCHS,
                                       validation_data = test_generator(),
                                       validation_steps = len(train_y) //BATCH_SIZE +1,
                                       use_multiprocessing=False,
                                       shuffle = True,
                                       verbose = 0,
                                       callbacks = [ ModelCheckpoint("model/", monitor='val_loss',
                                                     save_best_only= True, save_weights_only= True,
                                                     mode='auto', period=1),
                                                     MessagerCallback()]);
        #restoring best saved model
        model.load_weights("model/");
        model.save_weights("model.h5");

        with open('model.pkl', 'wb') as f:
           pickle.dump(info, f);

        tar = tarfile.open(MODEL_FILE, "w:gz");
        tar.add("model.pkl");
        tar.add("model.h5");
        tar.close();

        shutil.rmtree("model/");
        os.remove("model.pkl");
        os.remove("model.h5");

    else:
        print("Prediction mode...");

        tar = tarfile.open(MODEL_FILE);
        tar.extractall();
        tar.close();

        info = pickle.load( open( "model.pkl", "rb" ));

        keys = info[0];

        scaler = StandardScaler();
        scaler.mean_ = info[2][0];
        scaler.scale_ = info[2][1];
        scaler.n_features_in_ = 119;

        model = buildModel(info);
        model.load_weights("model.h5");

        fp = open(RESULT_FILE, "w");
        print(",".join(info[3]), file=fp);

        line = 0;
        for row in csv.reader(open(APPLY_FILE, "r")):
            line += 1;
            if line == 1: continue;

            smiles = row[0];
            try:
               m = Chem.MolFromSmiles(smiles);
               d = np.reshape(calcDescriptors(m), (1, -1));

               #scale for rdkit descriptors
               x1 = scaler.transform(d[:, :119]);
               x2 = d[:, 119:1143];
               x3 = d[:, 1143:2167];

               #mapping for atom pairs 
               x4 = np.zeros((1, len(keys)), dtype=np.float32);
               for k in d[:, 2167:].flatten():
                  try:                      
                     pos = keys.index(k);                     
                     x4[0, pos] = 1;
                  except: pass;

               #smiles conv1D
               x5 = np.zeros((1, 200), dtype=np.int8);
               smi = CleanSmileString(smiles);
               for i, c in enumerate(smi):
                  try: x5[0,i] = char_to_ix[c];
                  except: pass;
              
               pred = model.predict([x1, x2, x3, x4, x5])[-1].flatten();

               for p in pred[:-1]:
                  print(str(p), end=",", file=fp);
               print(str(pred[-1]), file=fp);

            except Exception as e:
               print(e);
               print(",".join(["error"]*len(info[3])), file=fp);

        os.remove("model.pkl");
        os.remove("model.h5");
        fp.close();

print("Relax!");
