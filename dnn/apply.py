import chainer
from chainer.dataset import convert

try:
    '''
     We try to import cupy because we need it in forward_gpu routine.
     If our system has only cpu: do nothing
    '''
    import cupy
    gpu_available = True
except ImportError:
    gpu_available = False

import os
import argparse
import pickle

import h5py
import numpy as np
from chainer import cuda, Variable
from chainer import serializers
from chainer.iterators import SerialIterator

from models import get_avaliable_models,cpu
from providers import LibSVMDataset

os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"   # see issue #152

_FNAME_ = 'best_model.hdf'
flatten = lambda l: [item for sublist in l for item in sublist]

parser = argparse.ArgumentParser(description='Entering point for applying deep NN in OChem')
parser.add_argument('--batch_size', default=1024, type=int,
                    help='A batch size for train iterations. Reduce if you have insufficient memory')
parser.add_argument('--apply_to', default='apply.libsvm', type=str,
                    help='A file apply to. Use only for debug!')
parser.add_argument('--best_model', default='best_model.hdf', type=str,
                    help='A path to the model file. Use only for debug!')
parser.add_argument('--outpath', default='apply_results.txt', type=str,
                    help='A path to save calc results. Use only for debug!')
parser.add_argument('--use_cpu', action='store_true',
                    help='Use this flag if you want to use only CPU for applying model')
parser.add_argument('--gpu_card', default=0, type=str,
                    help='The default GPU card to run calculations')

models_parsers = parser.add_subparsers(title="Models",
                                       description="models avaliable",
                                       help="Select a model with the same params.")

model_params_parsers = {model_name: models_parsers.add_parser(model_name, help=model.help()) for model_name, model in
                        get_avaliable_models().items()}
for k, v in model_params_parsers.items():
    get_avaliable_models()[k].add_parameters_for_parser(v)

args = parser.parse_args()

with h5py.File(args.best_model, 'r') as h5file:
    scaler = pickle.loads(h5file['/'].attrs['scaler'].tostring())
    # scalerX = pickle.loads(h5file['/'].attrs['scalerX'].tostring()) if 'scalerX' in h5file['/'].attrs else None
    input_dimension = h5file['/'].attrs['input_dimension']
    output_dimension = h5file['/'].attrs['output_dimension']
    binary_columns = h5file['/'].attrs['binary_columns'] if 'binary_columns' in h5file['/'].attrs else np.array([False]*output_dimension,dtype=bool)

def sigmod(x):
    return 1 / (1 + np.exp(-x))

examples = LibSVMDataset(args.apply_to, max_dimension=input_dimension,apply_mode=True)

model = get_avaliable_models()[args.model](input_dimension,
                                           output_dimension, args)


serializers.load_hdf5(args.best_model, model)

if not args.use_cpu:
    os.environ["CUDA_VISIBLE_DEVICES"]=str(args.gpu_card)
    cuda.cudnn_enabled = False  # Disable cuDNN
    model.to_gpu()

iterator = SerialIterator(examples, batch_size=args.batch_size, repeat=False, shuffle=False)
with chainer.using_config('train', False):
    y_pred_stack,y_real_stack = [],[]
    for batch in (iterator):
        X, y_real, types = convert.concat_examples(batch, device=model.device)
        with chainer.no_backprop_mode():
            y_pred = model(X)
            y_pred_stack.append(cpu(y_pred))
            y_real_stack.append(cpu(y_real))

y_pred_stack = np.concatenate(y_pred_stack, axis=0)
y_real_stack = np.concatenate(y_real_stack, axis=0)
y_pred_stack[:,binary_columns] = sigmod(y_pred_stack[:,binary_columns])
y_pred_stack[:,~binary_columns]  = scaler.inverse_transform(y_pred_stack[:,~binary_columns])
np.savetxt(args.outpath, y_pred_stack, delimiter=',', fmt='%.2f')
