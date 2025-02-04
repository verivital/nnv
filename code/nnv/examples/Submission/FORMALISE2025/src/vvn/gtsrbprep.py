# python standard library
import csv
import itertools
import onnxruntime
import os
import random
from functools import reduce
from collections import defaultdict
from typing import List

# third-party packages
import numpy as np

# local modules
from vvn.config import Config

# define global variables
PATH_TO_MODELS = os.path.join(os.getcwd(), 'models')
PATH_TO_DATA = os.path.join(os.getcwd(), 'data')

# set seed
# random.seed(42)

def prepare_filetree(config: Config):
    dst = 'gtsrb' 
    for va in ['relax', 'approx']:
        for length in ['4', '8', '16']:
            for eps_filename in [f'eps={e}_255' for e in range(1, 4)]:
                fp = build_output_filepath(config, eps_filename)

                # create the parent directories if they don't already exist
                os.makedirs(os.path.join(config.output_dir, dst, va, length), exist_ok=True)

    # make the results files once we know all directories have been made
    for va in ['relax', 'approx']:
        for length in ['4', '8', '16']:
            for eps_filename in [f'eps={e}_255' for e in range(1, 4)]:
                fp = build_output_filepath(config, eps_filename)

                # if the file doesn't exist yet, create it
                if not os.path.isfile(fp):
                    with open(fp, 'w', newline='') as f:
                        # write CSV headers
                        writer = csv.writer(f)
                        writer.writerow(['Sample Number', 'Result', 'Time', 'Method'])

def build_output_filepath(config: Config, filename=None, parent_only=False):
    """
    For our purposes, we will be naming the file based on
    the epsilon value used for the verification.
    """
    # error handling
    if not filename and not parent_only:
        raise Exception(f'No filename given. Please provide a filename when parent_only={parent_only}.')

    # convert all config values to str
    str_config = config.tostr()

    # get the values we need for building the output filepath
    output_dir = str_config.output_dir
    dst = str_config.ds_type
    va = str_config.ver_algorithm
    length = str_config.sample_len

    fp = os.path.join(output_dir, dst, va, length)

    return fp if parent_only else os.path.join(fp, filename) + '.csv'

def get_correct_samples(modelpath, datapath) -> tuple[list[int], list[int]]:
    outputs = []
    for sample_len in ['4', '8', '16']:
        
        # load the data + labels
        data = np.load(os.path.join(datapath, 'GTSRB', 'test', f'gtsrbvideo_{sample_len}f_test_data_seq.npy'))
        labels = np.load(os.path.join(datapath, 'GTSRB', 'test', f'gtsrbvideo_test_labels_seq.npy'))

        # select the model
        model = os.path.join(modelpath, f'gtsrb_{sample_len}f.onnx')

        # load the model + start onnx runtime session
        session = onnxruntime.InferenceSession(model)

        # specify input name for inference
        input_name = session.get_inputs()[0].name

        # run inference
        model_outputs = []

        for i in range(data.shape[0]):
            sample = data[i:i+1]
            sample = sample.transpose(0, 2, 1, 3, 4)
            sample = sample.astype(np.float32) # expects tensor(float) not tensor(double)
            output = session.run(None, {input_name: sample})
            model_outputs.append(output[0])
            
        # convert model_outputs from logits for each class to prediction
        model_outputs = [np.argmax(model_outputs[i], axis=1) for i in range(data.shape[0])]

        # get the true labels and compare them to the outputs
        labels = labels.astype(int).tolist()
        
        # filter for only correctly classified samples
        correct_samples = [i for i in range(data.shape[0]) if model_outputs[i] == labels[i]]
        
        # add the model outputs to the corresponding list of correct samples
        outputs.append(correct_samples)

    # return only samples that are correctly classified by all models
    outputs = list(reduce(lambda out, lst: out.intersection(lst), map(set, outputs)))

    return outputs

def generate_indices(config) -> tuple[list[int], list[int]]:
    # unpack config settings
    class_size = config.class_size

    # get the indices of all correctly classified samples
    correct_outputs = get_correct_samples(PATH_TO_MODELS, PATH_TO_DATA)

    # partition the correctly classified samples by class
    indices = defaultdict(list, {value: [i for i in correct_outputs] for value in range(0, 43)})

    # check that there are atleast 10 correctly classified samples for each class
    if not all(len(lst) >= class_size for lst in indices.values()):
        raise Exception("Not enough correctly classified samples for 'zoom_in'.")

    # randomly sample 10 of the correctly classified samples per class
    indices = [random.sample(indices[class_label], class_size) for class_label in indices.keys()]

    # flatten the list before returning
    indices = list(itertools.chain(*indices))

    # add 1 to all values of list because MATLAB uses 1-indexing
    indices = [v + 1 for v in indices]

    if len(indices) < class_size * 10:
        raise Exception("Not enough correctly classified samples.")

    # write the indices for the current experiment to its path (in this case it will be in random directory)
    with open(os.path.join(os.getcwd(), 'results', 'random', 'gtsrb_indices.txt'), 'w') as f:
        f.write(f'GTSRB Indices : {indices} \n')

    return indices


















