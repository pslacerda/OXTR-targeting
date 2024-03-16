from sklearn.model_selection import train_test_split
from matplotlib import pyplot as plt
import pandas as pd
from sklearn.neural_network import MLPClassifier, MLPRegressor
import numpy as np
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
from imblearn.under_sampling import RandomUnderSampler
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from scipy.stats import kendalltau, spearmanr, pearsonr


SEED = 0

def print_result(x, x_pred):
    print('kendall', kendalltau(x, x_pred))
    print('spearman', spearmanr(x, x_pred))
    print('pearson', pearsonr(x, x_pred))

CHUNKSIZE = 10240

est = MLPRegressor()

Test = pd.DataFrame([])
Test_y = []
Test_dg = []

peptibase = pd.read_csv("tests/peptibase.csv", sep=";").sample(frac=1)
assert 'Pdb' in peptibase.columns

with pd.read_csv("tests/fingerprints.csv", chunksize=CHUNKSIZE) as reader:
    for idx, chunk in enumerate(reader):
        df = pd.DataFrame.merge(peptibase, chunk,
            on="Pdb", how="left").drop([
                'Pdb', 'Affinity_Log', 'RMSD', 'S', 'S0', 'MD', 'CD',
                'Sequence1', 'Sequence3', 'Smiles', 'Chain', 'Resid', 'Eq',
                'Uniprot', 'Ecnumber', 'Constant'],
                axis=1)
        df = df.sample(frac=1).fillna(0)
        df = df[np.logical_and(df.Affinity != 0, df.Model == 1)]

        train, test = train_test_split(df, test_size=0.2)
        train_y = np.array(-np.log10(train.Affinity))
        test_y = np.array(-np.log10(test.Affinity))
        del train['Affinity'], test['Affinity']

        Test = pd.concat([Test, test])
        Test_y.extend(test_y)
        Test_dg.extend(test.Dg)
        
        est.partial_fit(train, train_y)
        print_result(train_y, est.predict(train))

print_result(Test_y, est.predict(Test))
print_result(Test_y, Test_dg)