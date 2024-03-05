import pandas as pd
import joblib


VERSION = 'nubbe'

scaler = joblib.load('scaler.joblib')
reg = joblib.load("regressor.joblib")

bbb = []
bbb_fname = f'bbb_{VERSION}.tsv'
with open(bbb_fname) as fbbb:
    for line in fbbb:
        name, pearm = line.strip().split('\t')
        if pearm != 'Non-Permeable':
            bbb.append(name.strip())

names = []
scores = []
with pd.read_csv(f"fingerprints_{VERSION}.csv", chunksize=5000) as reader:
    for idx, chunk in enumerate(reader):
        chunk = chunk.fillna(0)
        idx = chunk['Molecule.ChEMBL.ID'].str.rsplit('_', n=2).apply(lambda n: n[0] in bbb)
        names.extend(chunk[idx]['Molecule.ChEMBL.ID'])

        chunk = chunk[idx].drop('Molecule.ChEMBL.ID', axis=1)
        chunk = scaler.transform(chunk)
        scores.extend(reg.predict(chunk))

ret = list(zip(names, scores))
ret = list(sorted(ret, key=lambda r: -r[1]))[:10]

for (name, score) in ret:
    name = name.upper()
    print(name, round(score, 2))
