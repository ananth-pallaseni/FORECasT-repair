import os
import random
from subprocess import check_call
import numpy as np
import pandas as pd 
from predictor.features import calculateFeaturesForGenIndelFile, readFeaturesData

PACKAGE_DIRECTORY = os.path.dirname(os.path.realpath(__file__))
MODEL_WEIGHT_PATH = os.path.join(PACKAGE_DIRECTORY, 'model_weights.tsv')

# Load model weights as a knockout x feature matrix
MODEL_WEIGHTS = pd.read_csv(
    MODEL_WEIGHT_PATH, sep='\t'
).rename(
    columns={'Unnamed: 0': 'Gene'}
).set_index('Gene')

KOs = MODEL_WEIGHTS.index

# Generate predictions given a matrix of features, a matrix of all weights, and the name of a column in the weight matrix.
def _predict_from_features(features, weight_matrix, ko):
    a = weight_matrix.loc[ko].values[:-1].reshape(1, -1)
    b = weight_matrix.loc[ko].values[-1]
    y = (features * a).sum(axis=1) + b
    y = np.exp(y) / np.exp(y).sum()
    inserted_seq = features[[
        'I1_A','I2_AA','I2_AT','I2_AG','I2_AC',
        'I1_T','I2_TA','I2_TT','I2_TG','I2_TC',
        'I1_G','I2_GA','I2_GT','I2_GG','I2_GC',
        'I1_C','I2_CA','I2_CT','I2_CG','I2_CC',
    ]].idxmax(axis=1).mask(~y.index.str.contains('I'), '')
    inserted_seq = [s.split('_')[-1] for s in inserted_seq]
    y = pd.DataFrame({
        'Mutation': y.index,
        'Inserted Sequence': inserted_seq,
        'Prediction': y.values,
    })
    return y

# Generate features for a set of possible mutations given a target sequence and the index of the PAM in that targete.
def generate_forecast_features(target_seq, pam_idx, is_reverse=False):
    random_id = random.randint(0, 10000)
    tmp_genindels_file = f'tmp_genindels_{random_id}.txt' 
    tmp_features_file = f'tmp_features_{random_id}.txt'
    try:
        # Check that the indelgen program exists
        indelgen_path = os.environ.get("INDELGENTARGET_EXE")

        if (indelgen_path is None) or (not os.path.exists(indelgen_path)):
            raise ValueError('INDELGENTARGET_EXE environment variable not set or path does not exist. Set the environment variable INDELGENTARGET_EXE to point to the indelgentarget program.')

        # Generate temporary file containing possible mutations
        cmd = os.environ["INDELGENTARGET_EXE"] + ' %s %d %s' % (target_seq, pam_idx, tmp_genindels_file)
        check_call(cmd.split())
        assert os.path.exists(tmp_genindels_file)

        # Generate temporary feature file 
        calculateFeaturesForGenIndelFile(tmp_genindels_file, target_seq, pam_idx-3, is_reverse=is_reverse, out_file=tmp_features_file)

        # Load features
        df, _ = readFeaturesData(tmp_features_file)
        df = df.iloc[:, :-1]
        return df
        
    finally:
        # Remove temporary files
        if os.path.exists(tmp_genindels_file):
            os.remove(tmp_genindels_file)
        if os.path.exists(tmp_features_file):
            os.remove(tmp_features_file)

# Generate predicted mutation frequencies given a target sequence, the index of the PAM in that target, and the gene which was knocked out.
def predict(target_seq, pam_idx, ko):
    # Check that the model exists
    if ko not in KOs:
        raise ValueError(f'A model for {ko} deficient cells does not exist. Available options are: {KOs}')

    # Generate features
    features = generate_forecast_features(target_seq, pam_idx)

    # Predict
    predictions = _predict_from_features(features, MODEL_WEIGHTS, ko)
    return predictions




    