import pandas as pd
from rdkit.Chem.rdmolfiles import SDMolSupplier

def fetch_sd(project, number):
    directory = project + '/result.sd'
    data = SDMolSupplier(directory)

    for mol in data:
        print(mol)
    smiles = []

    # print(project, 'fetch!', number , 'times')
    pass

def fetch_csv(project):
    directory = project + '/smiles.csv'
    df = pd.read_csv(directory)
    df.drop(['Unnamed: 0'], axis=1, inplace=True)
    print(len(df), 'smiles fetched')
    return df.values.tolist()
    pass

def fetch_latent_space(project):
    print('fetch latent space :', project)

def fetch_fingerprints(project):
    print('fetch latent fingerprint :', project)
