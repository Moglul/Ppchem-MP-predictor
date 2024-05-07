import chemprop
import pandas as pd


arguments = [
        '--test_path', '/dev/null',
        '--preds_path', '/dev/null',
        '--checkpoint_dir', '../',
        '--no_cuda',
        '--num_workers', '0'
        ]

args = chemprop.args.PredictArgs().parse_args(arguments)

model_objects = chemprop.train.load_model(args=args)

args = chemprop.args.PredictArgs().parse_args(arguments)

def prediction(smiles) :
    smiles = [[smiles]]
    preds = chemprop.train.make_predictions(args=args, smiles=smiles, model_objects=model_objects)

    return preds[0][0]
