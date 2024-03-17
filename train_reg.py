from sklearn.model_selection import train_test_split, GridSearchCV
import pandas as pd
from scipy.stats import pearsonr
import joblib
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error
from sklearn.ensemble import HistGradientBoostingRegressor, RandomForestRegressor


chembl = pd.read_csv("data/chembl.csv", na_values=["NA", "None"]).sample(frac=1)
fpts = pd.read_csv(f"data/fingerprints_chembl.csv")
fpts = fpts[fpts.Model == 1].sample(frac=1)


df = pd.DataFrame.merge(chembl, fpts, on="Name", how="inner").drop(
    [
        "Molecule.Max.Phase",
        "Molecular.Weight",
        "X.RO5.Violations",
        "Compound.Key",
        "Standard.Type",
        "Standard.Relation",
        "Standard.Value",
        "Standard.Units",
        "Data.Validity.Comment",
        "Comment",
        "Uo.Units",
        "Ligand.Efficiency.BEI",
        "Ligand.Efficiency.LE",
        "Ligand.Efficiency.LLE",
        "Ligand.Efficiency.SEI",
        "Potential.Duplicate",
        "Assay.ChEMBL.ID",
        "Assay.Description",
        "Assay.Type",
        "BAO.Format.ID",
        "BAO.Label",
        "Assay.Organism",
        "Assay.Tissue.ChEMBL.ID",
        "Assay.Tissue.Name",
        "Assay.Cell.Type",
        "Assay.Subcellular.Fraction",
        "Assay.Parameters",
        "Assay.Variant.Accession",
        "Assay.Variant.Mutation",
        "Target.ChEMBL.ID",
        "Target.Name",
        "Target.Organism",
        "Target.Type",
        "Document.ChEMBL.ID",
        "Source.ID",
        "Source.Description",
        "Document.Journal",
        "Document.Year",
        "Cell.ChEMBL.ID",
        "Properties",
        "Action.Type",
        "Permeable",
        "Smiles",
        "AlogP",
    ],
    axis=1,
)

df = df.sample(frac=1).fillna(0)

df = df[df.groupby(by="Name")["S"].transform("max") == df.S].drop(
    ["Name", "Name_"], axis=1
)
y = df["Affinity"]
del df["Affinity"]


X_train, X_test, y_train, y_test = train_test_split(
    df, y, test_size=0.2, random_state=0
)


scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)
joblib.dump(scaler, "data/scaler.joblib")

models = []


for model_klass, params in [
    # (MLPRegressor, {
    #     # 'hidden_layer_sizes': [(100,), (50, 50), (100, 50)],
    #     'solver': ['lbfgs'],
    #     # 'alpha': [0.0001, 0.001, 0.01],
    #     # 'max_iter': [500]
    # }),
    # (RandomForestRegressor, {
    #     "n_estimators": [100, 200, 500],
    #     "max_depth": [3, 5, 10],
    #     "min_samples_split": [2, 5, 10],
    #     "min_samples_leaf": [1, 2, 4],
    # }),
    (HistGradientBoostingRegressor, {})
]:
    model = model_klass(random_state=0)
    gcv = GridSearchCV(model, params, cv=10, n_jobs=-1)
    gcv.fit(X_train_scaled, y_train)
    best_model = gcv.best_estimator_
    best_params = gcv.best_params_
    y_pred = best_model.predict(X_test_scaled)

    name = model_klass.__name__
    corr = pearsonr(y_pred, y_test)
    mse = mean_squared_error(y_test, y_pred)
    print(name, mse, corr)

    models.append((best_model, name, corr, mse))

models = sorted(models, key=lambda m: -m[2].statistic)
(reg, _, _, _) = models[0]
joblib.dump(reg, "data/regressor.joblib")
