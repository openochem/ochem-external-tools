try:
    '''
     We try to import cupy because we need it in forward_gpu routine.
     If our system has only cpu: do nothing
    '''
    import cupy

    gpu_available = True
except ImportError:
    gpu_available = False
import chainer
import chainer.functions as F
from chainer import cuda
from tqdm import tqdm
import chainer.functions as F

from line_profiler import LineProfiler

profiler = LineProfiler()  # DelMe after debuging

def cpu(x):
    x1 = x.data if isinstance(x, chainer.Variable) else x
    res = x1.get()  if isinstance(x1, cuda.ndarray) else x1
    return res

def profile(func):
    def inner(*args, **kwargs):
        profiler.add_function(func)
        profiler.enable_by_count()
        return func(*args, **kwargs)

    return inner


from pudb import set_trace;
import chainer.links as L
import numpy as np
from chainer.iterators import SerialIterator
from chainer import Chain, reporter
from chainer.optimizers import MomentumSGD, SMORMS3
from bisect import bisect
from robust import our_mean_squared_error, our_mean_absolute_error, our_sigmoid_cross_entropy
import chainer
import numpy as np
from chainer.dataset import convert
from robust import our_mean_squared_error, our_mean_absolute_error

_activations_ = ['relu', 'elu', 'prelu', 'sigmoid', 'none']


def softmax(x):
    max = np.max(x, axis=1, keepdims=True)  # returns max of each row and keeps same dims
    e_x = np.exp(x - max)  # subtracts each row with its max value
    sum = np.sum(e_x, axis=1, keepdims=True)  # returns sum of each row and keeps same dims
    f_x = e_x / sum
    return f_x


def return_xp(x):
    '''
    This function is used to return cupy or numpy depending on the type of x and the availability of gpu.
    '''
    xp = None
    if gpu_available:
        if cuda.get_device_from_array(x[0].data).id >= 0:
            xp = cuda
        else:
            xp = np
    else:
        xp = np
    return xp


class LinearBlock(Chain):
    def __init__(self, input, output, activation='relu', dropout=0.5, is_bn=True):
        self.input = input
        self.output = output
        self.dropout = dropout
        self.is_bn = is_bn
        self.activation = activation
        super(LinearBlock, self).__init__(
            layer=L.Linear(input, output),
        )
        if self.is_bn:
            self.add_link("bn", L.BatchNormalization(output))
        if self.activation == 'prelu':
            self.add_link("prelu", L.PReLU())

    def __repr__(self):
        return "Linear block {}->{} with" \
               " dropout: {}," \
               " activation: {}," \
               " bn: {}".format(self.input,
                                self.output,
                                self.dropout,
                                self.activation,
                                True if self.is_bn else False)  # str(self.layer)

    def __call__(self, x):
        h = self.layer(x)
        if self.is_bn:
            h = self.bn(h)

        if self.activation == 'relu':
            h = F.relu(h)
        elif self.activation == 'prelu':
            h = self.prelu(h)
        elif self.activation == 'sigmoid':
            h = F.sigmoid(h)
        elif self.activation == 'elu':
            h = F.elu(h)
        elif self.activation == 'none':
            pass
        else:
            raise NotImplementedError()

        if self.dropout > 0.:
            h = F.dropout(h, ratio=self.dropout)

        return h


class DenseLog2Net(Chain):
    def __init__(self, input, output, parser):
        self.activation = parser.activation
        self.n = int(np.log2(input))
        self.layers = [input] + [2 ** n for n in range(self.n, 2, -1)]
        self.layers[-1] = output
        self.n = len(list(self.layers))
        dropout_layers = len(self.layers) - 2
        self.dropout = [0.5] * int(dropout_layers * 0.8)
        self.dropout += [0.25] * int(dropout_layers * 0.2)
        self.dropout += [0.1]
        self.dropout += [0.]
        super(DenseLog2Net, self).__init__()
        for i in range(1, self.n):
            self.add_link("l{}".format(i), LinearBlock(int(self.layers[i - 1]),
                                                       int(self.layers[i]),
                                                       activation=self.activation if (self.n - i) > 1 else 'none',
                                                       dropout=self.dropout[i - 1]))

    def __call__(self, x):
        h = x
        for i in range(1, self.n):
            h = self.__getitem__("l{}".format(i))(h)
        return h

    @staticmethod
    def help():
        return """
        This model reduces neurons on half at each layer.
        The number of layers is calculated automatically.
        (Recommended to use if you don't want to fit parameters)
        """

    @staticmethod
    def add_parameters_for_parser(p):
        p.add_argument('--activation', choices=_activations_, default='elu', help='An activation function for NN')
        p.set_defaults(model='dense_log2')

    @staticmethod
    def recommended_optimizer():
        return SMORMS3()

    def __repr__(self):
        s = []
        for i in range(1, self.n):
            s.append(str(self.__getitem__("l{}".format(i))))
        return "\n".join(s)


class DenseNNet(Chain):
    def __init__(self, input, output, parser):
        self.n = parser.depth
        self.a = input
        self.alpha = parser.alpha
        self.activation = parser.activation
        self.b = self.a ** ((1. / self.alpha) / -self.n)
        # self.layers = map(lambda i: self.a * (self.b ** (i)), range(0, self.n + 1)) # last change by Sergey
        self.layers = list(map(lambda i: self.a * (self.b ** (i)), range(0, self.n + 1)))
        self.n = len(list(self.layers))
        # print self.layers
        self.layers[-1] = output
        dropout_layers = len(self.layers) - 2
        self.dropout = [0.5] * int(dropout_layers * 0.7)
        self.dropout += [0.25] * int(dropout_layers * 0.3)
        self.dropout += [0.1]
        self.dropout += [0.]

        super(DenseNNet, self).__init__()
        for i in range(1, self.n):
            self.add_link("l{}".format(i), LinearBlock(int(self.layers[i - 1]),
                                                       int(self.layers[i]),
                                                       activation=self.activation if (self.n - i) > 1 else 'none',
                                                       dropout=self.dropout[i - 1]
                                                       ))

    def __call__(self, x):
        h = x
        for i in range(1, self.n):
            h = self.__getitem__("l{}".format(i))(h)
        return h

    @staticmethod
    def help():
        return """
        This model build exponentially decreasing layers. The number of layers
        is required to be defined by user. Also an alpha value to scale the network is required
        (Recommended to use if you wish to fit parameters to achieve better results)
        """

    @staticmethod
    def add_parameters_for_parser(p):
        p.add_argument('--depth', type=int, required=True, help='A depth of the DNN')
        p.add_argument('--alpha', type=float, required=True, help='A value to scale NN')
        p.add_argument('--activation', choices=_activations_, default='elu', help='An activation function for NN')
        p.set_defaults(model='dense_exp')

    @staticmethod
    def recommended_optimizer():
        return SMORMS3()

    def __repr__(self):
        s = []
        for i in range(1, self.n):
            s.append(str(self.__getitem__("l{}".format(i))))
        return "\n".join(s)


class Dense7Net(Chain):
    def __init__(self, input, output, args):
        super(Dense7Net, self).__init__(
            l1=LinearBlock(input, 512, dropout=.5, activation='elu'),
            l2=LinearBlock(512, 256, dropout=.5, activation='elu'),
            l3=LinearBlock(256, 128, dropout=.5, activation='elu'),
            l4=LinearBlock(128, 64, dropout=.5, activation='elu'),
            l5=LinearBlock(64, 32, dropout=0.25, activation='elu'),
            l6=LinearBlock(32, 16, dropout=0.1, activation='elu'),
            l7=LinearBlock(16, output, dropout=.0, is_bn=False, activation='none')
        )

    def __call__(self, x):
        h = x
        for i in range(1, 8):
            h = self.__getitem__("l{}".format(i))(h)
        return h

    @staticmethod
    def help():
        return """
        Our honor dense7 model. Not recommended to use, because we have
        another adaptive to input dimensions models. Only for baseline.
        """

    @staticmethod
    def add_parameters_for_parser(p):
        p.set_defaults(model='dense7')

    @staticmethod
    def recommended_optimizer():
        return SMORMS3()

    def __repr__(self):
        s = []
        for i in range(1, 8):
            s.append(str(self.__getitem__("l{}".format(i))))
        return "\n".join(s)


class Dense7NetNew(Chain):
    def __init__(self, input, output, args):
        neurons = [2 ** n for n in range(4, 12)]
        pre_output = neurons[bisect(neurons, output)]
        super(Dense7NetNew, self).__init__(
            l1=LinearBlock(input, max(512, pre_output), dropout=.5, activation='elu'),
            l2=LinearBlock(max(512, pre_output), max(256, pre_output), dropout=.5, activation='elu'),
            l3=LinearBlock(max(256, pre_output), max(128, pre_output), dropout=.5, activation='elu'),
            l4=LinearBlock(max(128, pre_output), max(64, pre_output), dropout=.5, activation='elu'),
            l5=LinearBlock(max(64, pre_output), max(32, pre_output), dropout=0.25, activation='elu'),
            l6=LinearBlock(max(32, pre_output), max(16, pre_output), dropout=0.1, activation='elu'),
            l7=LinearBlock(max(16, pre_output), output, dropout=.0, is_bn=False, activation='none')
        )

    def __call__(self, x):
        h = x
        for i in range(1, 8):
            h = self.__getitem__("l{}".format(i))(h)
        return h

    @staticmethod
    def help():
        return """
        Our honor dense7 model with modification.
        """

    @staticmethod
    def add_parameters_for_parser(p):
        p.set_defaults(model='dense7new')

    @staticmethod
    def recommended_optimizer():
        return SMORMS3()

    def __repr__(self):
        s = []
        for i in range(1, 8):
            s.append(str(self.__getitem__("l{}".format(i))))
        return "\n".join(s)


class Merck(Chain):
    def __init__(self, input, output, args):
        self.bn = args.use_batch_normalisation
        self.original_scaling = args.use_original_scaling
        print("MESSAGE: model Merck, batch normalisation: {},"
              " original scaling: {}".format(self.bn, self.original_scaling))

        super(Merck, self).__init__(
            l1=LinearBlock(input, 4000, dropout=.25, activation='relu', is_bn=self.bn),
            l2=LinearBlock(4000, 2000, dropout=.25, activation='relu', is_bn=self.bn),
            l3=LinearBlock(2000, 1000, dropout=.10, activation='relu', is_bn=self.bn),
            l4=LinearBlock(1000, output, dropout=.0, activation='relu', is_bn=self.bn),
        )

    def __call__(self, x):
        h = x if not self.original_scaling else F.log1p(
            x.clip(-0.9))  # Scale by x = log(y+1) as proposed in original article
        for i in range(1, 5):
            h = self.__getitem__("l{}".format(i))(h)
        return h

    @staticmethod
    def help():
        return """
        Merk's model from https://github.com/Merck/DeepNeuralNet-QSAR
        """

    @staticmethod
    def add_parameters_for_parser(p):
        p.set_defaults(model='merck')
        p.add_argument('--use_batch_normalisation', action='store_true',
                       help='In original paper Batch normalisation is NOT used')
        p.add_argument('--use_original_scaling', action='store_true',
                       help='In original paper authors make log transformation y = log(x+1)')

    @staticmethod
    def recommended_optimizer():
        sgd = MomentumSGD(lr=0.05, momentum=0.9)
        return sgd

    def __repr__(self):
        s = []
        for i in range(1, 5):
            s.append(str(self.__getitem__("l{}".format(i))))
        return "\n".join(s)


class Dense2Net(Chain):
    def __init__(self, input, output, args):
        super(Dense2Net, self).__init__(
            fc=LinearBlock(input, args.hidden_neurons, dropout=args.dropout, is_bn=not args.disable_batch_normalisation,
                           activation=args.activation),
            out=LinearBlock(args.hidden_neurons, output, dropout=0., is_bn=not args.disable_batch_normalisation,
                            activation='none'),
        )

    def __call__(self, x):
        h = self.fc(x)
        return self.out(h)

    @staticmethod
    def help():
        return """
        The simplest model with one hidden layer.
        """

    @staticmethod
    def add_parameters_for_parser(p):
        p.set_defaults(model='dense2')
        p.add_argument('--hidden_neurons', type=int, default=12,
                       help='A number of hidden neurons')

        p.add_argument('--disable_batch_normalisation', action='store_true',
                       help='Use this flag to disable batch normalisation (enabled by default)')
        p.add_argument('--activation', choices=_activations_, default='elu',
                       help='An activation function for NN (default: elu)')
        p.add_argument('--dropout', type=float, default=0.25, help='A dropout ratio on hidden layer (default: 0.25)')

    @staticmethod
    def recommended_optimizer():
        return SMORMS3()

    def __repr__(self):
        s = []
        for i in ["fc", "out"]:
            s.append(str(self.__getitem__(i)))
        return "\n".join(s)


def get_avaliable_models():
    return {"dense_exp": DenseNNet,
            "dense7": Dense7Net,
            "dense7new": Dense7NetNew,
            "dense_log2": DenseLog2Net,
            "merck": Merck,
            "dense2": Dense2Net
            }
"""
Auxiliary functions
"""

def forward(model, X):
    """
    Just apply the model to the descriptors.
    """
    return model(X)

def calc_loss(lossfunc,y,t,mask,output_weights,types,device):
    """
    Calculate the loss function.
    """
    xp = cuda.get_array_module(y)
    mask = xp.asarray(mask, dtype=bool)
    types = xp.asarray(types, dtype=xp.int32)
    output_weights = xp.asarray(output_weights, dtype=xp.float32)
    weights = xp.expand_dims(output_weights[mask], 1).repeat(t.shape[0], 1).swapaxes(0, 1)
    loss = lossfunc(y[:, mask], t[:, mask], types[:, mask],
                    weights=weights)
    return loss



def calc_losses_and_outputs(model,iterator,output_weights,dict_of_loss_functions):
    #assert loss_functions is not None and len(loss_functions) > 0 and isinstance(loss_functions, list)
    assert dict_of_loss_functions is not None and len(dict_of_loss_functions) > 0 and isinstance(dict_of_loss_functions, dict)
    y_pred_stack,y_real_stack,loss_dict = [],[],{}
    #Fill the loss_dict with the empty values and keys from the loss_functions
    for key in dict_of_loss_functions.keys(): loss_dict[key] = []
    with chainer.using_config('train', False):
        for batch in (iterator):
            X, y_real, types = convert.concat_examples(batch, device=model.device)
            with chainer.no_backprop_mode():
                y_pred = forward(model, X)
                for mask,loss in dict_of_loss_functions.items():
                    loss_dict[tuple(mask)].append(calc_loss(loss,y_pred,y_real,mask,output_weights,types,model.device))
                y_pred_stack.append(cpu(y_pred))
                y_real_stack.append(cpu(y_real))
        #Concatenate the loss_dicts for each mask and calculate the mean
        for mask,loss_list in loss_dict.items():
            loss_dict[mask] = F.mean(F.stack(loss_list)).item()
        y_pred, y_real, losses = np.vstack(y_pred_stack), np.vstack(y_real_stack),loss_dict
        return y_pred, y_real, losses

class Loss(Chain):
    def __init__(self, predictor, output_weights,mask,supress_reporter=False):
        self.output_weights = output_weights
        self.supress_reporter = supress_reporter
        self.mask = mask
        super(Loss, self).__init__(predictor=predictor)

    def __call__(self, x, t, types):
        y = forward(self.predictor, x)
        loss = calc_loss(self.lossfunc,y, t, mask=self.mask, output_weights=self.output_weights,types=types,device=self.predictor.device)
        if not self.supress_reporter:  reporter.report({'loss': loss}, self)
        return loss

    def apply(self, iterator):
        return calc_losses_and_outputs(self.predictor,iterator,self.output_weights,{tuple(self.mask):self.lossfunc})

class Regressor(Loss):
    _lossfunc_ = {'mse': our_mean_squared_error, 'mae': our_mean_absolute_error}

    @staticmethod
    def get_lossfuncs_avaliable():
        return Regressor._lossfunc_.keys()

    @staticmethod
    def get_default_lossfunc():
        return 'mse'

    def __init__(self, predictor, lossfunc, output_weighs, mask, supress_reporter=False):
        self.lossfunc = self._lossfunc_[lossfunc]
        print("MESSAGE: regressor loss function: {}".format(lossfunc))
        super(Regressor, self).__init__(predictor=predictor, output_weights=output_weighs,
                                        supress_reporter=supress_reporter,mask=mask)


class Classifier(Loss):
    _lossfunc_ = {'sce': our_sigmoid_cross_entropy}

    @staticmethod
    def get_lossfuncs_avaliable():
        return Classifier._lossfunc_.keys()

    @staticmethod
    def get_default_lossfunc():
        return 'sce'

    def __init__(self, predictor, lossfunc, output_weighs, mask, supress_reporter=False):
        self.lossfunc = self._lossfunc_[lossfunc]
        print("MESSAGE: classifier loss function: {}".format(lossfunc))
        super(Classifier, self).__init__(predictor=predictor, output_weights=output_weighs,
                                         supress_reporter=supress_reporter,mask=mask)

    def apply(self, iterator):
        y_pred, y_real, loss = super().apply(iterator)
        return F.sigmoid(y_pred).data, y_real, loss

    # def __call__(self, x, t,types): # This function also accepts types of the outputs 1 - if more than, 0 - if equal
    #     return super().__call__(x, t,types,self.mask)


class Mixing(Chain):
    def __init__(self, predictor, loss, output_weights, binary_columns):
        class_loss, reg_loss = loss
        self.binary_columns = binary_columns

        self.output_weighs = output_weights
        self.classifier = Classifier(predictor=predictor, lossfunc=class_loss, output_weighs=output_weights,
                                     mask=binary_columns, supress_reporter=True)
        self.regressor = Regressor(predictor=predictor, lossfunc=reg_loss, output_weighs=np.copy(output_weights),
                                   mask=binary_columns, supress_reporter=True)
        self.dict_of_loss_functions = {tuple(self.binary_columns): self.classifier.lossfunc, tuple(~self.binary_columns): self.regressor.lossfunc}
        print("MESSAGE: Mixing uses loss functions: {},{}".format(class_loss, reg_loss))
        super(Mixing, self).__init__(predictor=predictor)

    def __call__(self, x, t, types):  # This function also accepts types of the outputs 1 - if more than, 0 - if equal
        class_loss = self.classifier(x, t, types)
        reg_loss = self.regressor(x, t, types)
        loss = reg_loss + class_loss

        print("MESSAGE: Classifier loss: {}".format(class_loss))
        print("MESSAGE: Regressor loss: {}".format(reg_loss))
        print("MESSAGE: Total loss: {}".format(loss))

        reporter.report({'loss': loss}, self)
        return loss


    def apply(self, iterator):
        return calc_losses_and_outputs(self.predictor,iterator,self.output_weighs,self.dict_of_loss_functions)


#if isinstance(y_real, cupy.ndarray):  # Do convert back if y_real is cupy.ndarray
#    y_real = y_real.get()
#if isinstance(types, cupy.ndarray):
#    types = types.get()
#if isinstance(self.output_weights, cupy.ndarray):
#    weights = self.output_weights.get()
#else:
#    weights = self.output_weights

