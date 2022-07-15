
import sys
import os
import csv
import tarfile
import configparser
import glob
from pathlib import Path

import macau
import numpy as np
from scipy.sparse import coo_matrix

#always dense?
def loadDescriptors(fname):

    x = [];
    t_cols = [];
    t_rows = [];

    line = 0;
    total = 0;

    with open(fname, "r") as fn:
        reader = csv.reader(fn, delimiter=',', quotechar='\"');
 
        #skip the first line with captions 
        next(reader, None);

        for row in reader: 
           col = 0;
           for el in row:
               total += 1;
               try:
                   u = float(el);
                   x.append(u);

                   t_cols.append(col);
                   t_rows.append(line);                   
               except:
                   pass;                   
    
               col += 1;
           line += 1;
 
    cols = int(total / line);
    rows = line;   
    
    mat = coo_matrix((x, (t_rows, t_cols)), shape=(rows, cols));
    return mat;

def cleanDirectoryFromMacau():
    dpath = Path(os.getcwd());
    for fp in dpath.glob('m-*'):
        fp.unlink();

def packMacauModel():
 
    tar = tarfile.open("macau.tar", "w");

    for fp in glob.glob('./m-*'):
        tar.add(fp);

    tar.close();

def loadMacauModel():
   
    meanvalue = np.loadtxt("m-meanvalue.csv").tolist()

    links = [];
    mean = [];
    V = [];

    for fp in glob.glob("./m-sample*-U1-link.csv"):
        links.append(np.loadtxt(fp, delimiter = ","));

    for fp in glob.glob("./m-sample*-U1-latentmean.csv"):
        mean.append(np.loadtxt(fp, delimiter = ","));

    for fp in glob.glob("./m-sample*-U2-latents.csv"):
        V.append(np.loadtxt(fp, delimiter = ","));

    return links, mean, V, meanvalue;

def readConfig():
    
    conf = configparser.ConfigParser();
    conf.read("config.ini")

    if "macau" not in conf:
        print("Please, check the config file.");
        sys.exit(0);
     
    conf = conf["macau"];

    test_persent = 0.2;
    num_latent = 32;
    precision = 0.5;
    burnin = 5;
    nsamples = 3;

    try:

        test_percent = float(conf["test_percent"]);
        num_latent   = int(conf["num_latent"]);
        try:
            precision    = float(conf["precision"]);
        except:
            precision = conf["precision"]
        burnin       = int(conf["burnin"]);
        nsamples     = int(conf["samples"]);

    except KeyError as err:
        print("Error in the config.ini file. Check the tag: ", err);
        raise;
    except ValueError as val:
        print("Error in values in the config.ini.", val);
        raise;

    return test_percent, num_latent, precision, burnin, nsamples;

def main():

    if(len(sys.argv) != 2):
        print("Usage: {} --training|--prognosis".format(sys.argv[0]));
        sys.exit(1);

    if(sys.argv[1] == "--training"):

        try:
            test_percent, num_latent, precision, burnin, nsamples = readConfig();
        except:
            sys.exit(1);

        cleanDirectoryFromMacau();

        x = loadDescriptors("descrs.csv");
        y = loadDescriptors("prop.csv");

        #running factorization (Macau)
        result = macau.macau(Y = y,
                             Ytest = test_percent,
                             side = [x, None],
            		     num_latent = num_latent,
                             precision  = precision,
                             burnin = burnin,
                             nsamples = nsamples,
                             save_prefix = "m");
                 
        #some macau exception's handler?
        packMacauModel();
        cleanDirectoryFromMacau();

        print("Probaly we have finished. All the samples are packed into the macau.tar archive. Relax!");
        sys.exit(0);

    elif(sys.argv[1] == "--prognosis"):
        print("Prognosis...");

        try:
             tar = tarfile.open("macau.tar");
             tar.extractall();
             tar.close();
        except:
             print("Problems with unpacking the macau.tar file.");

        print("Loading samples into the memory."); 
   
        links, mean, V, meanvalue = loadMacauModel();
        if( not(len(links) == len(mean) == len(V))):
            print("Cannot use the model. Sizes for links, U and V are different.");

        fr = open("results.csv", "w"); 
        writer = csv.writer(fr, delimiter=',', quotechar='\"');

        with open("descrs.csv", "r") as fn:
            reader = csv.reader(fn, delimiter=',', quotechar='\"');
            next(reader, None);

            for row in reader: 
               x = [];

               for el in row:
                   u = float(el);
                   x.append(u);

               result = [];

               for i in range(len(V)):
                  
                   X = np.dot(x, links[i].transpose()) + mean[i]; 
                   Y = X.dot(V[i]) + meanvalue;
 
                   result.append(Y);

               for i in range(result[0].shape[0]):
                   for j in range(1, len(result)):
                       result[0][i] += result[j][i];

               y = result[0] / len(result);             
               writer.writerow(y);

        fr.close();

        cleanDirectoryFromMacau();
        print("Probably finished. Check the file results.csv.");

        sys.exit(0);
    else:
        print("Unkonwn option. Valid are: --training or --prognosis");
        sys.exit(1);
             

    return;


if __name__ == "__main__":
    main();

'''
'''
