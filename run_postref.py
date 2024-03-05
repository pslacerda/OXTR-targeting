import pandas as pd
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, kendalltau
from matplotlib import pyplot as plt
import joblib
import seaborn as sb


def parse_pdbqt_poses(fname):
    
    with open(fname) as file:
        if 'VINA' not in '\n'.join(file.readlines()):
            return
    with open(fname) as file:
        for line in file:
            if line.startswith("MODEL"):
                model = int(line.split()[1])
            elif "VINA RESULT" in line:
                dg = -float(line.split()[3])
            elif "INTER + INTRA" in line:
                inter_intra = float(line.split()[4])
            elif "INTER:" in line:
                inter = float(line.split()[2])
            elif " INTRA:" in line:
                intra = float(line.split()[2])
            elif "UNBOUND" in line:
                unbound = float(line.split()[2])
            elif "active torsions" in line:
                torsions = int(line.split()[1])
            elif line.startswith("ROOT"):
                return {
                    "Model": model,
                    "Dg": dg,
                    "Torsions": torsions,
                    "InterIntra": inter_intra,
                    "Inter": inter,
                    "Intra": intra,
                    "Unbound": unbound,
                }


bbb = pd.read_csv('bbb_vina8.tsv', sep='\t')
bbb = bbb[bbb.Permeability == 'Permeable'].Name

est = joblib.load("regressor.joblib")
chembl = pd.read_csv("chembl.csv", na_values=["NA", "None"])

names = []
exp = []
scores1 = []
scores2 = []
klass = []
with pd.read_csv(f"fingerprints_vina8.csv", chunksize=5000) as reader:
    for index, chunk in enumerate(reader):
        chunk = chunk.reset_index().drop('index', axis=1)
        chunk = chunk.fillna(0)
        # chunk = chunk[~chunk['Molecule.ChEMBL.ID'].str.upper().isin(bbb)]
        idx = chunk.groupby('Molecule.ChEMBL.ID')['S'].transform('max') == chunk['S']
        chunk = chunk[idx]
        for i, row in chunk.iterrows():
            name = row['Molecule.ChEMBL.ID']
            num_atoms = row['NumAtoms']
            
            del row['Molecule.ChEMBL.ID']
            score = est.predict([row])
            if len(score) != 1:
                continue
            if np.all(score < 0 or score > 12):
                continue
            xp = chembl[chembl['Molecule.ChEMBL.ID'] == name]['Affinity']
            if len(xp) != 1:
                continue
            names.append(name)
            scores1.append(float(score[0]))
            scores2.append(float(row['Dg']))
            exp.append(xp.iloc[0])
            klass.append(num_atoms < 50)


print('Score:', round(pearsonr(exp, scores1).statistic, 2))
print('Vina:', round(pearsonr(exp, scores2).statistic, 2))



print('Score top 20:', round(pearsonr(exp[:20], scores1[:20]).statistic, 2))
print('Vina top 20:', round(pearsonr(exp[:20], scores2[:20]).statistic, 2))


fig, (ax1) = plt.subplots(1, 1)

sb.scatterplot(x=exp, y=scores1, marker='x', color='blue', label='Post-refined', ax=ax1)
sb.scatterplot(x=exp, y=scores2, marker='+', color='red', label='Raw Vina', ax=ax1)
ax1.set_xlabel('pKi (exerimental)')
ax1.set_ylabel('pKi (predictive)')

plt.legend()
plt.show()
