try:
    '''
     We try to import cupy because we need it in forward_gpu routine.
     If our system has only cpu: do nothing
    '''
    import cupy
except:
    pass
import numpy
import sys
from pudb import set_trace; 
import numpy as np
import numpy.ma as ma
from chainer import Variable
from chainer import cuda
from chainer import function
from chainer.functions import softmax_cross_entropy,sigmoid_cross_entropy,sigmoid
import chainer.functions as F
from chainer.utils import type_check


class OurMeanAbsoluteError(function.Function):
    """Our patch to mean absolute error function to deal with NaNs."""

    def check_type_forward(self, in_types):
        type_check.expect(in_types.size() == 4)
        type_check.expect(
            in_types[0].dtype == numpy.float32,
            in_types[1].dtype == numpy.float32,
            in_types[0].shape == in_types[1].shape,
            in_types[1].shape == in_types[2].shape,
            in_types[2].shape == in_types[3].shape
        )

    def forward_cpu(self, inputs):
        y_pred, y_real,types,weights = inputs
        y_real[np.isnan(y_real)] = y_pred[np.isnan(y_real)]
        self.diff = y_pred - y_real
        self.diff[(types == 1.) & (self.diff > 0)] = 0 # Our fix to handle more_than values
        self.diff[(types == 2.) & (self.diff < 0)] = 0  # Our fix to handle less_than values
        self.diff *= weights  # Our fix to concern about weights of outputs
        diff = self.diff.ravel()
        return numpy.array(abs(diff).sum() / diff.size, dtype=diff.dtype),

    def forward_gpu(self, inputs):
        y_pred, y_real,types,weights = inputs
        y_real[cupy.isnan(y_real)] = y_pred[cupy.isnan(y_real)]
        self.diff = y_pred - y_real
        self.diff[(types == 1.) & (self.diff > 0)] = 0 # Our fix to handle more_than values
        self.diff[(types == 2.) & (self.diff < 0)] = 0  # Our fix to handle less_than values
        self.diff *= weights  # Our fix to concern about weights of outputs
        diff = self.diff.ravel()
        return abs(diff).sum() / diff.dtype.type(diff.size),

    def backward(self, inputs, gy):
        xp = cuda.get_array_module(*inputs)
        coeff = gy[0] * gy[0].dtype.type(1. / self.diff.size)
        gx0 = coeff * xp.sign(self.diff)
        return gx0, -gx0


def our_mean_absolute_error(x0, x1, types):
    """Mean absolute error function.

    This function computes mean absolute error between two variables. The mean
    is taken over the minibatch.

    """
    return OurMeanAbsoluteError()(x0, x1, types)


class OurMeanSquaredError(function.Function):
    """Our patch to Mean squared error (a.k.a. Euclidean loss) function. to deal with NaNs"""

    def check_type_forward(self, in_types):
        type_check.expect(in_types.size() == 4)
        type_check.expect(
            in_types[0].dtype == numpy.float32,
            in_types[1].dtype == numpy.float32,
            in_types[0].shape == in_types[1].shape,
            in_types[1].shape == in_types[2].shape,
            in_types[2].shape == in_types[3].shape
        )

    def forward_cpu(self, inputs):
        y_pred, y_real,types,weights = inputs
        y_real[np.isnan(y_real)] = y_pred[np.isnan(y_real)]
        self.diff = y_pred - y_real
        self.diff[(types == 1.) & (self.diff > 0)] = 0 # Our fix to handle more_than values
        self.diff[(types == 2.) & (self.diff < 0)] = 0  # Our fix to handle less_than values
        self.diff *= weights # Our fix to concern about weights of outputs
        diff = self.diff.ravel()
        return numpy.array(diff.dot(diff) / diff.size, dtype=diff.dtype),

    def forward_gpu(self, inputs):
        y_pred, y_real,types,weights = inputs
        y_real[cupy.isnan(y_real)] = y_pred[cupy.isnan(y_real)]
        self.diff = y_pred - y_real
        self.diff[(types == 1.) & (self.diff > 0)] = 0 # Our fix to handle more_than values
        self.diff[(types == 2.) & (self.diff < 0)] = 0  # Our fix to handle less_than values
        self.diff *= weights # Our fix to concern about weights of outputs
        diff = self.diff.ravel()
        return diff.dot(diff) / diff.dtype.type(diff.size),

    def backward(self, inputs, gy):
        coeff = gy[0] * gy[0].dtype.type(2. / self.diff.size)
        gx0 = coeff * self.diff
        return gx0, -gx0


def our_mean_squared_error(x0, x1 ,types,weights):
    return OurMeanSquaredError()(x0, x1, types,weights)

def our_mean_absolute_error(x0, x1 ,types,weights):
    return OurMeanAbsoluteError()(x0, x1, types,weights)

# def our_softmax_cross_entropy(x, t ,types,weights):
#     xp = cuda.get_array_module(x)
#     t[xp.isnan(t)] = -1
#     return F.softmax_cross_entropy(x, t, ignore_label=-1)

def our_sigmoid_cross_entropy(x, t, types, weights):
    xp = cuda.get_array_module(x)
    t1 = t.copy()
    t1[xp.isnan(t)] = -1
    loss_collector = None
    for i in range(x.shape[1]):
        r = F.sigmoid_cross_entropy(x[:,i],t1[:,i].astype(xp.int32))
        if loss_collector is None: loss_collector = r*weights[0,i]
        else: loss_collector += r
    return loss_collector / x.shape[1]


    #return F.bernoulli_nll(x,t,reduce='mean') #t*F.log(x)+(1-t)*F.log(1-x)

#def our_softmax_cross_entropy(x, t ,types,weights):
#    xp = cuda.get_array_module(x)
#    mask = xp.isnan(t)
#    t = t[~mask]
#    x = x[~mask]
#    # if isinstance(t,cupy.ndarray):
#    #     t[cupy.isnan(t)] = -1
#    # else:
#    #     t[np.isnan(t)] = -1
#    x = sigmoid(x)
#    # t = Variable(xp.array(t, dtype=xp.int32))
#    a = t*F.log(x)
#    b = (1-t)*F.log(1-x)
#    return a + b#t*F.log(x)+(1-t)*F.log(1-x)
#    #print(type(x))
#    #return sigmoid_cross_entropy(x, t.astype(cupy.int32), weights)

# def our_softmax_cross_entropy(x, t ,types,weights):
#     if isinstance(t,cupy.ndarray):
#         t[cupy.isnan(t)] = -1
#         class_weight = cupy.array([0.1,0.9])
#     else:
#         t[np.isnan(t)] = -1
#         class_weight = np.array([0.1,0.9])
#     return softmax_cross_entropy(x, t, ignore_label=-1,class_weight=class_weight)


class OurRobustToNanScaler():
    """
    This class is equal to StandardScaler from sklearn but can work with NaN's (ignoring it) but
    sklearn's scaler can't do it.
    """

    def fit(self, data):
        masked = ma.masked_invalid(data)
        self.means = np.mean(masked, axis=0)
        self.stds = np.std(masked, axis=0)

    def fit_transform(self, data):
        self.fit(data)
        return self.transform(data)

    def transform(self, data):
        masked = ma.masked_invalid(data)
        masked -= self.means
        masked /= self.stds
        return ma.getdata(masked)

    def inverse_transform(self, data):
        masked = ma.masked_invalid(data)
        masked *= self.stds
        masked += self.means
        return ma.getdata(masked)
