# tensorflow backend
print('tensorflow backend')
from os import environ
environ['KERAS_BACKEND'] = 'tensorflow'

import pandas as pd
from chemvae import mol_utils as mu
from chemvae.vae_utils import VAEUtils
from rdkit.Chem import AllChem as Chem

def access_bit(data, num):
    base = int(num // 8)
    shift = int(num % 8)
    return (data[base] & (1<<shift)) >> shift

def encode_latent_fp(project):
    # fetch data
    print('fetching data...')
    directory = project + '/smiles.csv'
    df = pd.read_csv(directory)
    df.drop(['Unnamed: 0'], axis=1, inplace=True)
    smiles_group = df.values.tolist()
    print('fetched smiles :', len(smiles_group))

    # load autoencoder
    vae = VAEUtils(directory='models/zinc_properties')

    # encode data
    count = 0
    Xs = []
    Zs = []
    latent_group = []

    print('canonizing smiles...')
    smiles_list = [ mu.canon_smiles(data[0]) for data in smiles_group ]
    
    print('encoding mols to hot...')
    for idx, smiles in enumerate(smiles_list):
        try:
            Xs.append(vae.smiles_to_hot(smiles, canonize_smiles=True))
        except:
            count += 1
            del smiles_group[idx]
    print('failed : ', count)
    print('encoded smiles :', len(smiles_group))

    print('encoding mols...')
    count = 0
    for idx, X in enumerate(Xs):
        try:
            Zs.append(list(vae.encode(X)[0]))
        except:
            count += 1
            del smiles_group[idx]
    print('failed:', count)
    print('encoded smiles :', len(smiles_group))
    
    latent_group.extend(list(zip([ data[0] for data in smiles_group ], Zs, [ data[1] for data in smiles_group ])))

    # save latent data
    df = pd.DataFrame(latent_group)
    directory = project + '/data/latent_group.csv'
    df.to_csv(directory)

    # encode fingerprints
    mols = [ Chem.MolFromSmiles(smiles) for smiles in smiles_list ]
    latent_group = []
    rdk_fp_group = []
    pattern_fp_group = []
    layered_fp_group = []
    MACCSKeys_fp_group = []

    # rdk fp
    print('encoding rdk fp...')
    fps = [ Chem.RDKFingerprint(mol) for mol in mols ]
    fpBits = [ [ int(char) for char in fp.ToBitString() ] for fp in fps ]
    rdk_fp_group.extend(list(zip([ data[0] for data in smiles_group ], fpBits, [ data[1] for data in smiles_group ])))
    df = pd.DataFrame(rdk_fp_group)
    directory = project + '/data/rdk_fp_group.csv'
    df.to_csv(directory)

    # pattern fp
    print('encoding pattern fp...')
    fps = [ Chem.PatternFingerprint(mol) for mol in mols ]
    fpBits = [ [ int(char) for char in fp.ToBitString() ] for fp in fps ]
    pattern_fp_group.extend(list(zip([ data[0] for data in smiles_group ], fpBits, [ data[1] for data in smiles_group ])))
    df = pd.DataFrame(pattern_fp_group)
    directory = project + '/data/pattern_fp_group.csv'
    df.to_csv(directory)

    # layered fp
    print('encoding layered fp...')
    fps = [ Chem.LayeredFingerprint(mol) for mol in mols ]
    fpBits = [ [ int(char) for char in fp.ToBitString() ] for fp in fps ]
    layered_fp_group.extend(list(zip([ data[0] for data in smiles_group ], fpBits, [ data[1] for data in smiles_group ])))
    df = pd.DataFrame(layered_fp_group)
    directory = project + '/data/layered_fp_group.csv'
    df.to_csv(directory)

    # MACCSKeys fp
    print('encoding MACCSKeys fp...')
    fps = [ Chem.GetMACCSKeysFingerprint(mol) for mol in mols ]
    fpBits = [ [ int(char) for char in fp.ToBitString() ] for fp in fps ]
    MACCSKeys_fp_group.extend(list(zip([ data[0] for data in smiles_group ], fpBits, [ data[1] for data in smiles_group ])))
    df = pd.DataFrame(MACCSKeys_fp_group)
    directory = project + '/data/MACCSKeys_fp_group.csv'
    df.to_csv(directory)

    print('finish!')