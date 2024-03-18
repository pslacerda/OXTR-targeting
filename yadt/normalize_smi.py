with open("ligands.smi") as f:
    for line in f:
        print(line)
        names, smiles = line.strip().rsplit(" ", 1)
        print(smiles, names)
