import pandas as pd
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from matplotlib import pyplot as plt
import joblib
import seaborn as sb


chembl = pd.read_csv("data/chembl.csv", na_values=["NA", "None"])
est = joblib.load("data/regressor.joblib")
scaler = joblib.load("data/scaler.joblib")

names = []
exp = []
scores1 = []
scores2 = []
klass = []
with pd.read_csv(f"data/fingerprints_chembl.csv", chunksize=5000) as reader:
    for index, chunk in enumerate(reader):
        chunk = chunk.reset_index().drop("index", axis=1)
        chunk = chunk.fillna(0)

        # chunk = chunk[~chunk['Mol ecule.ChEMBL.ID'].str.upper().isin(bbb)]
        idx = chunk.groupby("Name")["S"].transform("max") == chunk["S"]
        chunk = chunk[idx]
        for i, row in chunk.iterrows():
            name = row["Name"]
            num_atoms = row["NumAtoms"]

            del row["Name"]

            chunk = scaler.transform([row])
            score = est.predict(chunk)

            if len(score) != 1:
                continue
            if np.all(score < 0 or score > 12):
                continue
            xp = chembl[chembl["Name"] == name]["Affinity"]
            if len(xp) != 1:
                continue

            names.append(name)
            scores1.append(float(score[0]))
            scores2.append(float(row["Dg"]))
            exp.append(xp.iloc[0])
            klass.append(num_atoms < 50)


print("Score:", round(pearsonr(exp, scores1).statistic, 2))
print("Vina:", round(pearsonr(exp, scores2).statistic, 2))


print("Score top 20:", round(pearsonr(exp[:20], scores1[:20]).statistic, 2))
print("Vina top 20:", round(pearsonr(exp[:20], scores2[:20]).statistic, 2))


fig, (ax1) = plt.subplots(1, 1)

sb.scatterplot(x=exp, y=scores1, marker="x", color="blue", label="Post-refined", ax=ax1)
sb.scatterplot(x=exp, y=scores2, marker="+", color="red", label="Raw Vina", ax=ax1)
ax1.set_xlabel("pKi (exerimental)")
ax1.set_ylabel("pKi (predictive)")

plt.legend()
plt.show()
