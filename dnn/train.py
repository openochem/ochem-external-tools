try:
    '''
     We try to import cupy because we need it in forward_gpu routine.
     If our system has only cpu: do nothing
    '''
    import cupy
    gpu_available = True
except ImportError:
    gpu_available = False

import argparse
import logging
import os
from os.path import expanduser, join

import chainer
from chainer import cuda
from chainer import training
from chainer.datasets import split_dataset_random
from chainer.iterators import SerialIterator
from chainer.optimizers import RMSprop, Adam
from chainer.training import extension
from chainer.training import extensions

from evaluators import NewEvaluator
from models import *
from providers import LibSVMDataset
import time

global run_on
global start_time

os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"   # see issue #152

try:
    from cupy.cuda import cudnn
    _cudnn_version = cudnn.getVersion()
    print("MESSAGE: Found CuDNN version: {}"
          " and we can use it (or disable by --disable_cudnn)".format(_cudnn_version))
except:
    print("MESSAGE: Can't load CuDNN. Looks like CuDNN is not installed.")

log = logging.getLogger('trainer')
home = expanduser("~")

_optimizers_= {'adam':Adam(),
               'sgd':MomentumSGD(lr=0.005,momentum=0.9),
               'rmsprop':RMSprop(),
               'smorms3':SMORMS3()}

parser = argparse.ArgumentParser(description='Entering point for train deep NN in OChem')
parser.add_argument('--file',default="train.libsvm",type=str, help='A file with train set')
parser.add_argument('--batch_size',default=1024,type=int, help='A batch size for train iterations. Reduce if you have insufficient memory')
parser.add_argument('--internal_train_test_ratio',default=0.8,type=float, help='A batch size for train iterations. Reduce, if you have insufficient memory or faced with other troubles.')
parser.add_argument('--n_epochs',default=1000,type=int, help='A maximum number of epochs')
parser.add_argument('--optimizer',choices=_optimizers_.keys(),help='An optimizer with default params')
parser.add_argument('--lossfunc', choices=list(Regressor.get_lossfuncs_avaliable()) + list(Classifier.get_lossfuncs_avaliable()),default=Regressor.get_default_lossfunc(),help='A loss function')
#parser.add_argument('--lossfuncs_for_mixed', choices=,default=,help='')

parser.add_argument('--seed',type=int, help='A seed value')

parser.add_argument('--use_gradient_clipping',action='store_true', help='Use gradient clipping to prevent NAN\'s (recommended ONLY for SGD)')
parser.add_argument('--gradient_clipping',type=float,default=1., help='A value of scaling all gradient arrays to fit to the defined L2 norm threshold')

parser.add_argument('--use_gradient_noise',action='store_true', help='Add gradient noise (sometimes it produces better results)')
parser.add_argument('--gradient_noise_eta',type=float,default=0.01, help='eta value for gradient noise')


parser.add_argument('--use_cpu',action='store_true', help='Use this flag if you want to use only CPU')
parser.add_argument('--debug',action='store_true', help='Enable chainer debug mode. Use in case of errors and NANs')
parser.add_argument('--clear_cupy_cache',action='store_true', help='Remove $HOME/.cupy directory before launch. Use in case of problems with CUDA')
parser.add_argument('--disable_cudnn',action='store_true', help='Use this flag if you want to disable cuDNN library (only for GPU mode)')

parser.add_argument('--gpu_card', default=0, type=str,
                    help='The default GPU card to run calculations')

models_parsers = parser.add_subparsers(title="Models",
                                       description="models avaliable",
                                       help="Use models from the list. Use -h to show the params of the model i.e. train.py dense_exp -h")


model_params_parsers = {model_name:models_parsers.add_parser(model_name, help=model.help()) for model_name, model in get_avaliable_models().items()}
for k,v in model_params_parsers.items():
    get_avaliable_models()[k].add_parameters_for_parser(v)

args = parser.parse_args()

os.environ["CUDA_VISIBLE_DEVICES"]=str(args.gpu_card)

if args.clear_cupy_cache:
    import shutil
    try:
        shutil.rmtree(join(home,".cupy")) # BEWARE TO CHANGE IT This code removes $HOME/.cupy directory !!! #
    except:
        pass


if args.debug:
   print("MESSAGE: enable chainer debug mode!")
   chainer.set_debug(True)

gpu = False if args.use_cpu else True

if args.seed:
    print("MESSAGE: Set seed value {} for numpy".format(args.seed))
    np.random.seed(args.seed)

if gpu and args.seed:
    print("MESSAGE: Set seed value {} for cupy".format(args.seed))
    cupy.random.seed(args.seed)


provider = LibSVMDataset(args.file,debug=True)

if args.internal_train_test_ratio == 0 or args.internal_train_test_ratio == 1:
    test_provider = provider
    train_provider = provider
else:
    train_provider,test_provider = split_dataset_random(provider,int(len(provider)*args.internal_train_test_ratio),seed=args.seed)

scaler = provider.get_scaler()
output_weights = provider.get_output_weights()

train_iterator = SerialIterator(train_provider, batch_size=args.batch_size)
test_iterator = SerialIterator(test_provider,
                               batch_size=args.batch_size,repeat=False,shuffle=False)

model = get_avaliable_models()[args.model](provider.input_dimension,
                                           provider.output_dimension, args)

if args.debug:  print("MESSAGE: binary columns: ", provider.get_binary_columns())
binary_columns = provider.get_binary_columns()

if provider.type_of_a_problem() == "regression":
    model = Regressor(model, args.lossfunc,output_weights,mask=~binary_columns)
    print("MESSAGE: all columns are continuous, so we use Regressor")
elif provider.type_of_a_problem() == "binary classification":
    model = Classifier(model, args.lossfunc,output_weights,mask=binary_columns)
    print("MESSAGE: all columns are binary, so we use Classifier")
elif provider.type_of_a_problem() == "mixed":
    model = Mixing(model, ('sce','mse'),output_weights,binary_columns)
    print("MESSAGE: some columns are binary, some are continuous, so we use Mixed")
else:
    raise ValueError("Something is clearly wrong!")

optimizer = _optimizers_[args.optimizer] if args.optimizer else model.recommended_optimizer()
print ("MESSAGE: optimizer {}".format(optimizer))
optimizer.setup(model)

if gpu:
    print("MESSAGE: running on GPU")
    model.to_gpu()
    cuda.cudnn_enabled = True
else:
    print("MESSAGE: running on CPU")

if args.disable_cudnn:
    cuda.cudnn_enabled = False
    print("MESSAGE: cuDNN is disabled by user")


if args.use_gradient_clipping:
    print ("MESSAGE: We will use gradient clipping with value {}".format(args.gradient_clipping))
    optimizer.add_hook(chainer.optimizer.GradientClipping(args.gradient_clipping))

if args.use_gradient_noise:
    print ("MESSAGE: We will use gradient noise with eta {}".format(args.gradient_noise_eta))
    optimizer.add_hook(chainer.optimizer.GradientNoise(args.gradient_noise_eta))

updater = training.StandardUpdater(train_iterator, optimizer, device=0 if gpu else -1)
trainer = training.Trainer(updater, (args.n_epochs, 'epoch'), out="")


evaluator = NewEvaluator(scaler,int(args.n_epochs*0.2),iterator=test_iterator, target=model, device=0 if gpu else -1,binary_columns=binary_columns,input_dimension=provider.input_dimension)
trainer.extend(evaluator)

@extension.make_extension(trigger=(1,'epoch'), priority=-100)
def plot_epoch(trainer):
    train_loss = trainer.observation['main/loss'].array
    train_loss = train_loss if isinstance(train_loss,np.ndarray) else train_loss.get()
    print("MESSAGE: train score: {:.4f} / validation score: {:.4f} / at epoch: {} / {} / elapsed time: {:d}s ".
          format(train_loss,evaluator.current_loss,trainer.updater.epoch,run_on,int(time.time() - start_time)))

trainer.extend(plot_epoch)
trainer.extend(extensions.LogReport(log_name='train.log'))
start_time = time.time()
run_on = ('gpu' if gpu else 'cpu')
trainer.run()
