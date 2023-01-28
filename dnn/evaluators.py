try:
    import cupy
    gpu_available = True
except ImportError:
    gpu_available = False

import copy
import pickle
import sys
import os
from os.path import join

import chainer
import h5py
import numpy as np
from chainer.training.extensions import Evaluator
from sklearn.metrics import roc_auc_score, classification_report, precision_recall_curve, auc, confusion_matrix

flatten = lambda l: [item for sublist in l for item in sublist]


class NewEvaluator(Evaluator):
    def __init__(self, scaler,max_epochs_without_update, binary_columns, input_dimension, *args, **kwargs):
        self.previous_result = 1e10
        self.scaler = scaler
        # self.scalerX = scalerX

        self.epochs_without_update = 0
        self.binary_columns = binary_columns
        self.input_dimension = input_dimension
        self.max_epochs_without_update = max_epochs_without_update
        self.device = kwargs['device']
        super(NewEvaluator, self).__init__(*args, **kwargs)

    def evaluate(self):
        iterator = self._iterators['main']
        target = self._targets['main']
        # output_weights = target.output_weights
        y_pred, y_real, losses = target.apply(copy.copy(iterator))

        y_pred_continuous = y_pred[:, ~self.binary_columns]
        y_real_continuous = y_real[:, ~self.binary_columns]
        additional_metrics = {}
        if np.any(self.binary_columns):
            y_pred_binary = y_pred[:, self.binary_columns]
            y_real_binary = y_real[:, self.binary_columns]

            mask = ~np.isnan(y_real_binary)
            auc_roc = roc_auc_score(y_real_binary[mask], y_pred_binary[mask])
            # print("Classification report: {}".format(classification_report(y_real_binary[mask],y_pred_binary[mask])))
            # Data to plot precision - recall curve
            precision, recall, thresholds = precision_recall_curve(y_real_binary[mask], y_pred_binary[mask])
            # Use AUC function to calculate the area under the curve of precision recall curve
            auc_precision_recall = auc(recall, precision)
            additional_metrics = {'auc_roc': auc_roc, 'auc_precision_recall': auc_precision_recall}
            # print("AUC PR = {}".format(auc_precision_recall))
        loss = np.mean(np.array(list(losses.values()))).item()
        self.current_loss = loss
        # Merge dicts
        result = {**{'loss': loss}, **additional_metrics}

        if loss < self.previous_result:
            self.previous_result = loss
            fname = join("best_model.hdf")
            chainer.serializers.save_hdf5(fname, target.predictor)
            with h5py.File(fname, 'r+') as f:
                f["/"].attrs['previous_result'] = self.previous_result
                f["/"].attrs['epochs_without_update'] = self.epochs_without_update
                f["/"].attrs['input_dimension'] = self.input_dimension
                f["/"].attrs['output_dimension'] = int(self.binary_columns.shape[0])
                f["/"].attrs['binary_columns'] = self.binary_columns
                f["/"].attrs['scaler'] = np.void(pickle.dumps(self.scaler))
                # f["/"].attrs['scalerX'] = np.void(pickle.dumps(self.scalerX))
            self.epochs_without_update = 0
        else:
            self.epochs_without_update += 1
        if (((self.epochs_without_update >= self.max_epochs_without_update) or os.path.exists(
                'stop')) and os.path.exists('best_model.hdf')):
            if os.path.exists('stop'):
                os.remove('stop')
            print("MESSAGE: Early stopping after {} epochs without update".format(self.epochs_without_update))
            sys.exit(0)
        return result

'''
class MyEvaluator(Evaluator):
    def __init__(self, scaler, max_epochs_without_update, restart, weigths, binary_columns, *args, **kwargs):
        self.previous_result = 1e10  # For regression
        self.epochs_without_update = 0
        if restart:
            h5file = h5py.File("best_model.hdf", 'r')
            self.previous_result = h5file['/'].attrs['previous_result']
            self.epochs_without_update = h5file['/'].attrs['epochs_without_update']
            self.binary_columns = h5file['/'].attrs['binary_columns']
            h5file.close()

        self.max_epochs_without_update = max_epochs_without_update
        self.scaler = scaler
        self.device = kwargs['device']
        sys.exit()
        self.binary_scaler = np.void(pickle.dumps(scaler))
        super(MyEvaluator, self).__init__(*args, **kwargs)

    def evaluate(self):
        iterator = self._iterators['main']
        target = self._targets['main']
        it = copy.copy(iterator)
        first_batch = True
        with chainer.using_config('train', False):
            for batch in it:
                X, y_real, types = convert.concat_examples(batch)
                X = Variable(X)
                if gpu_available:
                    if isinstance(self.device, cupy.cuda.Device) >= 0:
                        X.to_gpu()
                with chainer.no_backprop_mode():
                    res = target(X)
                # if self.device >= 0:
                res.to_cpu()

                if first_batch:
                    y_pred_scaled = res.data
                    y_real_scaled = y_real
                    y_types = types
                    y_weights = weights
                    input_dimension = X.shape[-1]
                    output_dimension = y_pred_scaled.shape[-1]
                    first_batch = False
                else:
                    y_pred_scaled = np.vstack((y_pred_scaled, res.data))
                    y_real_scaled = np.vstack((y_real_scaled, y_real))
                    y_types = np.vstack((y_types, types))
                    y_weights = np.vstack((y_weights, weights))
        y_types = y_types.reshape(1, -1)[0]  # Beware! here can be an error
        y_weights = y_weights.reshape(1, -1)[0]  # Beware! here can be an error
        y_pred_normal = self.scaler.inverse_transform(y_pred_scaled).reshape(1, -1)[0]
        y_real_normal = self.scaler.inverse_transform(y_real_scaled).reshape(1, -1)[0]
        validation_rmse = float((our_mean_squared_error(y_pred_normal, y_real_normal, y_types, y_weights) ** 0.5).data)
        result = {  # "validation_scaled_rmse": our_mean_squared_error(y_real_scaled, y_pred_scaled,y_types) ** 0.5,
            "validation_rmse": validation_rmse,
            "validation_mae": float(our_mean_absolute_error(y_pred_normal, y_real_normal, y_types, y_weights).data),
            "epochs_without_update": self.epochs_without_update}

        if validation_rmse < self.previous_result:
            self.previous_result = validation_rmse
            fname = join("best_model.hdf")
            if gpu_available:
                if isinstance(self.device, cupy.cuda.Device) >= 0:
                    target.to_cpu()
            chainer.serializers.save_hdf5(fname, target)
            if gpu_available:
                if isinstance(self.device, cupy.cuda.Device) >= 0:
                    target.to_gpu()
            f = h5py.File(fname, 'r+')
            # f["/"].attrs['scaler'] = self.binary_scaler
            # f["/"].attrs['previous_result'] = self.previous_result
            # f["/"].attrs['epochs_without_update'] = self.epochs_without_update
            # f["/"].attrs['input_dimension'] = input_dimension
            # f["/"].attrs['output_dimension'] = output_dimension
            f.close()
            self.epochs_without_update = 0
        else:
            self.epochs_without_update += 1
        if (((self.epochs_without_update >= self.max_epochs_without_update) or os.path.exists(
                'stop')) and os.path.exists('best_model.hdf')):
            if os.path.exists('stop'):
                os.remove('stop')
            print("MESSAGE: Early stopping after {} epochs without update".format(self.epochs_without_update))
            sys.exit(0)
            # chainer.serializers.save_hdf5("results/model{:03.2f}_rmse.hdf".format(normal_rmse),target)
         return result
'''
