import sys
from io import BytesIO
from os.path import isfile
import os

import chainer
import numpy as np

from robust import OurRobustToNanScaler

flatten = lambda l: [item for sublist in l for item in sublist]

class LibSVMDataset(chainer.dataset.DatasetMixin):

    def _generate_transformed_outputs(self):
        def _generate_masks_for_more_than_symbol(x):
            if x.startswith(b'>'): return 1.
            elif x.startswith(b'<'): return 2.
            else: return 0.
        def _convert_to_float_with_more_than(x):
            x = x.decode('ascii')
            if x == '': return np.NaN
            elif (x.startswith('>') or x.startswith('<')): return float(x[1:])
            else: return float(x)

        def is_binary_categorial_column(x):
            return np.all(np.isin(x[np.isnan(x) == False], [0., 1.]))


        #Prepare two auxiliary functions 
        more_than_mask_processor = np.vectorize(_generate_masks_for_more_than_symbol)
        convert_to_float = np.vectorize(_convert_to_float_with_more_than)

        str_outputs = list(map(lambda v: v[0], self.data.values()))
        self.object_outputs = np.genfromtxt(BytesIO("\n".join(str_outputs).encode()), delimiter=",", dtype=object)
        self.outputs_types = more_than_mask_processor(self.object_outputs).astype(np.ubyte) 
        self.outputs = np.zeros_like(self.outputs_types, dtype=np.float32)

        #Set correct types for outputs
        self.outputs = convert_to_float(self.object_outputs).astype(np.float32)
        if self.outputs.shape == ():
            self.outputs = self.outputs.reshape((1,1))
        self.binary_columns = np.apply_along_axis(is_binary_categorial_column, 0,self.outputs) #Mark each column as binary or continuous
        if (self.binary_columns.shape) == ():
            self.binary_columns = self.binary_columns.reshape(-1,)
        if self.outputs.shape == ():
            self.outputs = self.outputs.reshape(-1,)
        if self.outputs.ndim == 1:
            self.outputs = self.outputs.reshape(-1, 1)
        self.outputs_types = self.outputs_types.reshape(self.outputs.shape)
        self.output_dimension = self.outputs.shape[1]
        self.scaler = OurRobustToNanScaler()
        #Transform only continuous columns
        self.transformed_outputs = np.zeros_like(self.outputs, dtype=np.float32)
        self.transformed_outputs[:,~self.binary_columns] = self.scaler.fit_transform(self.outputs[:,~self.binary_columns])
        self.transformed_outputs[:,self.binary_columns] = self.outputs[:,self.binary_columns]

    def _max_dimension(self):
        maximum_index = 0
        for line_id in self.data:
            line = self.data[line_id]
            records = line[1:]
            for record in records:
                index = int(record.split(':')[0])
                maximum_index = max(index, maximum_index)
        return maximum_index

    def __init__(self, path, max_dimension=None, debug=False,weights_file='weights.txt',apply_mode=False):
        self.file = path
        self.apply_mode = apply_mode
        self.debug = debug
        self.data = open(path, 'r').read().split("\n")
        self.data = list(filter(lambda str: str != "", self.data))  # Filter possible empty lines

        self.len = len(self.data)
        self.data = list(map(lambda str: str.split(" "), self.data))
        self.data = {i: v for i, v in zip(range(len(self.data)), self.data)}
        self._generate_transformed_outputs()
        if not max_dimension:
            self.input_dimension = self._max_dimension()
        else:
            self.input_dimension = max_dimension
        self.matrix = np.zeros((self.len, self.input_dimension), dtype=np.float32)
        if isfile(weights_file):
            if debug: print("MESSAGE: weights file is found")
            output_weights = np.genfromtxt(weights_file,delimiter=',').astype(np.float32)
            if (len(output_weights) != self.outputs.shape[1]) and (not self.apply_mode):
                print("MESSAGE: error, output shape is {} but weights shape is {}".format(len(output_weights),self.outputs.shape[1]))
                sys.exit(1)
        else:
            if debug: print("MESSAGE: weights file is NOT found")
            output_weights = np.ones(self.outputs.shape[1]).astype(np.float32)
        self.output_weights = output_weights
        if debug: print("MESSAGE: start matrix formation")
        for i in self.data.keys():
            sparse_records = self.data[i][1:]
            for record in sparse_records:
                splitted = record.split(":")
                pos, v = int(splitted[0]) - 1, float(splitted[1])
                self.matrix[i, pos] = v
        if debug: print("MESSAGE: matrix is formed")
        # self.matrix = self.Xscaler.fit_transform(self.matrix)
        if os.path.exists('stop') and os.path.exists('best_model.hdf'):
            print("MESSAGE: stopped by request")
            os.remove('stop')
            sys.exit(0)


    def __len__(self):
        return self.len

    def get_scaler(self):
        return self.scaler

    def set_scaler(self,scaler):
        self.scaler = scaler

    def get_binary_columns(self):
        return self.binary_columns

    def get_output_weights(self):
        return self.output_weights

    def type_of_a_problem(self):
        if np.all(~self.binary_columns): #If all columns are continuous
            return "regression"
        elif np.all(self.binary_columns): #If all columns are binary
            return "binary classification"
        else:
            return "mixed" #If some columns are binary and some are continuous

    def get_example(self, item):
        return self.matrix[item], self.transformed_outputs[item] , self.outputs_types[item]