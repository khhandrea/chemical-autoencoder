import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

def learn_tsne_pca(project):
    # fetch data
    print('(1/5) fetching latent space group...', end='')
    latent_group = []
    directory = project + '/data/latent_group.csv'
    df = pd.read_csv(directory)
    for i in range(len(df)):
        points = list(map(float, df.loc[i, '1'].strip('[]').split(', ')))
        latent_group.append([df.loc[i, '0'], points, int(df.loc[i, '2'])])

    print('\r(2/5) fetching rdk fp group...', end='')
    rdk_fp_group = []
    directory = project + '/data/rdk_fp_group.csv'
    df = pd.read_csv(directory)
    for i in range(len(df)):
        points = list(map(int, df.loc[i, '1'].strip('[]').split(', ')))
        rdk_fp_group.append([df.loc[i, '0'], points, int(df.loc[i, '2'])])

    print('\r(3/5) fetching pattern fp group...', end='')
    pattern_fp_group = []
    directory = project + '/data/pattern_fp_group.csv'
    df = pd.read_csv(directory)
    for i in range(len(df)):
        points = list(map(float, df.loc[i, '1'].strip('[]').split(', ')))
        pattern_fp_group.append([df.loc[i, '0'], points, int(df.loc[i, '2'])])

    print('\r(4/5) fetching layered fp group...', end='')
    layered_fp_group = []
    directory = project + '/data/layered_fp_group.csv'
    df = pd.read_csv(directory)
    for i in range(len(df)):
        points = list(map(float, df.loc[i, '1'].strip('[]').split(', ')))
        layered_fp_group.append([df.loc[i, '0'], points, int(df.loc[i, '2'])])


    print('\r(5/5) fetching MACCSKeys fp group...')
    MACCSKeys_fp_group = []
    directory = project + '/data/MACCSKeys_fp_group.csv'
    df = pd.read_csv(directory)
    for i in range(len(df)):
        points = list(map(float, df.loc[i, '1'].strip('[]').split(', ')))
        MACCSKeys_fp_group.append([df.loc[i, '0'], points, int(df.loc[i, '2'])])

    n_components = 2
    model_pca = PCA(n_components=n_components)
    model_tsne = TSNE(n_components=n_components)

    print('\r(1/10) {:30s}'.format('training latent tsne...'), end='')
    points_latent_tsne = model_tsne.fit_transform([ points[1] for points in latent_group ])
    df = pd.DataFrame(points_latent_tsne)
    df['smiles'] = pd.DataFrame([points[0] for points in latent_group])
    df['values'] = pd.DataFrame([points[2] for points in latent_group])
    directory = project + '/data/points_latent_tsne.csv'
    df.to_csv(directory)

    print('\r(2/10) {:30s}'.format('training latent pca...'), end='')
    points_latent_pca = model_pca.fit_transform([ points[1] for points in latent_group ])
    df = pd.DataFrame(points_latent_pca)
    df['smiles'] = pd.DataFrame([points[0] for points in latent_group]) 
    df['values'] = pd.DataFrame([points[2] for points in latent_group])
    directory = project + '/data/points_latent_pca.csv'
    df.to_csv(directory)

    print('\r(3/10) {:30s}'.format('training rdk fingerprint tsne...'), end='')
    points_rdk_fp_tsne = model_tsne.fit_transform([ points[1] for points in rdk_fp_group ])
    df = pd.DataFrame(points_rdk_fp_tsne)
    df['smiles'] = pd.DataFrame([points[0] for points in rdk_fp_group]) 
    df['values'] = pd.DataFrame([points[2] for points in rdk_fp_group])
    directory = project + '/data/points_rdk_fp_tsne.csv'
    df.to_csv(directory)

    print('\r(4/10) {:30s}'.format('training rdk fingerprint pca...'), end='')
    points_rdk_fp_pca = model_pca.fit_transform([ points[1] for points in rdk_fp_group ])
    df = pd.DataFrame(points_rdk_fp_pca)
    df['smiles'] = pd.DataFrame([points[0] for points in rdk_fp_group]) 
    df['values'] = pd.DataFrame([points[2] for points in rdk_fp_group])
    directory = project + '/data/points_rdk_fp_pca.csv'
    df.to_csv(directory)

    print('\r(5/10) {:30s}'.format('training pattern fingerprint tsne...'), end='')
    points_pattern_fp_tsne = model_tsne.fit_transform([ points[1] for points in pattern_fp_group ])
    df = pd.DataFrame(points_pattern_fp_tsne)
    df['smiles'] = pd.DataFrame([points[0] for points in pattern_fp_group]) 
    df['values'] = pd.DataFrame([points[2] for points in pattern_fp_group])
    directory = project + '/data/points_pattern_fp_tsne.csv'
    df.to_csv(directory)

    print('\r(6/10) {:30s}'.format('training pattern fingerprint pca...'), end='')
    points_pattern_fp_pca = model_pca.fit_transform([ points[1] for points in pattern_fp_group ])
    df = pd.DataFrame(points_pattern_fp_pca)
    df['smiles'] = pd.DataFrame([points[0] for points in pattern_fp_group]) 
    df['values'] = pd.DataFrame([points[2] for points in pattern_fp_group])
    directory = project + '/data/points_pattern_fp_pca.csv'
    df.to_csv(directory)

    print('\r(7/10) {:30s}'.format('training layered fingerprint tsne...'), end='')
    points_layered_fp_tsne = model_tsne.fit_transform([ points[1] for points in layered_fp_group ])
    df = pd.DataFrame(points_layered_fp_tsne)
    df['smiles'] = pd.DataFrame([points[0] for points in layered_fp_group]) 
    df['values'] = pd.DataFrame([points[2] for points in layered_fp_group])
    directory = project + '/data/points_layered_fp_tsne.csv'
    df.to_csv(directory)

    print('\r(8/10) {:30s}'.format('training layered fingerprint pca...'), end='')
    points_layered_fp_pca = model_pca.fit_transform([ points[1] for points in layered_fp_group ])
    df = pd.DataFrame(points_layered_fp_pca)
    df['smiles'] = pd.DataFrame([points[0] for points in layered_fp_group]) 
    df['values'] = pd.DataFrame([points[2] for points in layered_fp_group])
    directory = project + '/data/points_layered_fp_pca.csv'
    df.to_csv(directory)

    print('\r(9/10) {:30s}'.format('training MACCSKeys fingerprint tsne...'), end='')
    points_MACCSKeys_fp_tsne = model_tsne.fit_transform([ points[1] for points in MACCSKeys_fp_group ])
    df = pd.DataFrame(points_MACCSKeys_fp_tsne)
    df['smiles'] = pd.DataFrame([points[0] for points in MACCSKeys_fp_group]) 
    df['values'] = pd.DataFrame([points[2] for points in MACCSKeys_fp_group])
    directory = project + '/data/points_MACCSKeys_fp_tsne.csv'
    df.to_csv(directory)

    print('\r(10/10) {:30s}'.format('training MACCSKeys fingerprint pca...'), end='')
    points_MACCSKeys_fp_pca = model_pca.fit_transform([ points[1] for points in MACCSKeys_fp_group ])
    df = pd.DataFrame(points_MACCSKeys_fp_pca)
    df['smiles'] = pd.DataFrame([points[0] for points in MACCSKeys_fp_group]) 
    df['values'] = pd.DataFrame([points[2] for points in MACCSKeys_fp_group])
    directory = project + '/data/points_MACCSKeys_fp_pca.csv'
    df.to_csv(directory)