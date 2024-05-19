import pandas as pd
import logging
import chemprop
import sys
from io import StringIO

# Redirection of the standard output and standard error to avoid printing the output of the prediction function
old_stdout = sys.stdout
old_stderr = sys.stderr
sys.stdout = StringIO()
sys.stderr = StringIO()

arguments = [
        '--test_path', '/dev/null',
        '--preds_path', '/dev/null',
        '--checkpoint_dir', '../../',
        '--no_cuda',
        '--num_workers', '0',
        '--features_generator', 'rdkit_2d_normalized',  # Include the same feature generator used during training
        '--no_features_scaling'
        ]

args = chemprop.args.PredictArgs().parse_args(arguments)

model_objects = chemprop.train.load_model(args=args)

args = chemprop.args.PredictArgs().parse_args(arguments)

def prediction(smiles) :
    smiles = [[smiles]]
    preds = chemprop.train.make_predictions(args=args, smiles=smiles, model_objects=model_objects)

    return preds[0][0]

# Reset the standard output and standard error
sys.stdout = old_stdout
sys.stderr = old_stderr