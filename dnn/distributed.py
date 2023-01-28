import argparse
from models import *
from providers import LibSVMDataset
from os.path import join
import chainer
import h5py
import numpy as np
import pickle 
from glob import glob

parser = argparse.ArgumentParser(description='Entering point for distributed training')
parser.add_argument('--file',default="train.libsvm",type=str, help='A file with train set')
parser.add_argument('--mode',choices=['create','mean'],required=True,type=str, help='create -- creates the file with random weights, mean -- mean the files best_model_*.hdf')


models_parsers = parser.add_subparsers(title="Models",description="models avaliable",help="Use models from the list. Use -h to show the params of the model i.e. train.py dense_exp -h")
model_params_parsers = {model_name:models_parsers.add_parser(model_name, help=model.help()) for model_name, model in get_avaliable_models().items()}
for k,v in model_params_parsers.items():
        get_avaliable_models()[k].add_parameters_for_parser(v)

args = parser.parse_args()
if args.mode == 'create':
    provider = LibSVMDataset(args.file,debug=True)
    model = get_avaliable_models()[args.model](provider.input_dimension,
                                               provider.output_dimension, args)

    fname = join("best_model.hdf")
    chainer.serializers.save_hdf5(fname, model)
    f = h5py.File(fname,'r+')
    f["/"].attrs['scaler'] = np.void(pickle.dumps(provider.scaler))
    f["/"].attrs['previous_result'] = 9999.
    f["/"].attrs['epochs_without_update'] = 0
    f["/"].attrs['input_dimension'] = provider.input_dimension
    f["/"].attrs['output_dimension'] = provider.output_dimension
    f.close()
else:
    weighs_files = glob("best_model_*.hdf")
    print("MESSSAGE: found {} files".format(len(weighs_files)))
    reference = weighs_files[0]
    with h5py.File(reference, 'r') as h5file:
        scaler = h5file['/'].attrs['scaler']
        input_dimension = h5file['/'].attrs['input_dimension']
        output_dimension = h5file['/'].attrs['output_dimension']
        previous_result = h5file['/'].attrs['previous_result']
        epochs_without_update = h5file["/"].attrs['epochs_without_update']

    reference_model = get_avaliable_models()[args.model](input_dimension,
                                                   output_dimension, args)
    chainer.serializers.load_hdf5(reference, reference_model)

    for model_file in weighs_files[1:]:    
        model = get_avaliable_models()[args.model](input_dimension,
                                                   output_dimension, args)
        chainer.serializers.load_hdf5(model_file, model)
        for param in model.namedparams():
            ref_data = [p[1] for p in reference_model.namedparams() if p[0] == param[0]][0].data
            another_data = param[1].data
            ref_data = np.average([ref_data,another_data],axis=0)

    chainer.serializers.save_hdf5("best_model.hdf", reference_model)
    f = h5py.File(reference,'r+')
    f["/"].attrs['scaler'] = scaler
    f["/"].attrs['previous_result'] = previous_result
    f["/"].attrs['epochs_without_update'] = epochs_without_update
    f["/"].attrs['input_dimension'] = input_dimension
    f["/"].attrs['output_dimension'] = output_dimension
    f.close()
        