from sklearn.model_selection import train_test_split
from matplotlib import pyplot as plt
import pandas as pd
from sklearn.neural_network import MLPClassifier, MLPRegressor
import numpy as np
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
from imblearn.under_sampling import RandomUnderSampler
import pandas as pd

SEED = 0


CHUNKSIZE = 40960

est_rmsd = MLPClassifier()
Test_rmsd_x = pd.DataFrame([])
Test_rmsd_y = []

real_rmsd = {}

peptibase = pd.read_csv("tests/peptibase.csv", sep=";").sample(frac=1)
with pd.read_csv("tests/fingerprints.csv", chunksize=CHUNKSIZE) as reader:
    for idx, chunk in enumerate(reader):
        assert 'Model' in chunk.columns
        assert 'Pdb' in peptibase.columns
        df = pd.DataFrame.merge(chunk, peptibase,
            on="Pdb").drop([
                'Pdb', 'Affinity_Log', 'Affinity',
                'Sequence1', 'Sequence3', 'Smiles', 'Chain', 'Resid', 'Eq',
                'Uniprot', 'Ecnumber', 'Constant'],
                axis=1)
        
        df = df.sample(frac=1).fillna(0)
        
        train_x, test_x = train_test_split(df, test_size=0.2)

        train_x = train_x
        train_y = np.array(train_x.RMSD < 2, dtype=bool)
        
        test_x = test_x
        test_y = test_x.RMSD < 2
        
        del train_x['RMSD'], test_x['RMSD']

        Test_rmsd_x = pd.concat([Test_rmsd_x, test_x])
        Test_rmsd_y.extend(test_y)
        real_rmsd.extend(test_y[test_x.Model == 1])

        s = RandomUnderSampler(random_state=SEED)
        try:
            x_res, y_res = s.fit_resample(train_x, train_y)
        except:
            continue
        
        est_rmsd.partial_fit(x_res, y_res, classes=(True, False))
        print(confusion_matrix(y_res, est_rmsd.predict(x_res)))

        rmsd = pd.DataFrame.merge(chunk[chunk.Model == 1].drop_duplicates('Pdb'), peptibase,
                                 on="Pdb").RMSD
        
        print(sum(rmsd <= 2)/len(rmsd))


s = RandomUnderSampler(random_state=SEED)
x_res, y_res = s.fit_resample(Test_rmsd_x, Test_rmsd_y)
print(confusion_matrix(y_res, est_rmsd.predict(x_res), normalize="all"))

ConfusionMatrixDisplay.from_estimator(est_rmsd, x_res, y_res, normalize="all")
plt.show()

print(sum(real_rmsd)/len(real_rmsd))
