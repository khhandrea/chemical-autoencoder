import pandas as pd
from matplotlib import pyplot as plt
from sklearn.cluster import KMeans
from rdkit.Chem import PandasTools, SDWriter
from rdkit.Chem import AllChem as Chem
from IPython.display import display

def plot_molecules(project,
                alpha = 0.5,
                size = 5,
                color = '#74b9ff',
                bg_color = '#ffffff', # '#CAD3C8'
                fontsize = 30):
    print('(1/10) fetching latent space tsne points...', end='')
    points_latent_tsne = []
    directory = project + '/data/points_latent_tsne.csv'
    df = pd.read_csv(directory)
    for i in range(len(df)):
        points_latent_tsne.append([df.loc[i, 'smiles'], [float(df.loc[i, '0']), float(df.loc[i, '1'])], int(df.loc[i, 'values'])])

    print('\r(2/10) fetching latent space pca points...', end='')
    points_latent_pca = []
    directory = project + '/data/points_latent_pca.csv'
    df = pd.read_csv(directory)
    for i in range(len(df)):
        points_latent_pca.append([df.loc[i, 'smiles'], [float(df.loc[i, '0']), float(df.loc[i, '1'])], int(df.loc[i, 'values'])])

    print('\r(3/10) fetching rdk fingerprint tsne points...', end='')
    points_rdk_fp_tsne = []
    directory = project + '/data/points_rdk_fp_tsne.csv'
    df = pd.read_csv(directory)
    for i in range(len(df)):
        points_rdk_fp_tsne.append([df.loc[i, 'smiles'], [float(df.loc[i, '0']), float(df.loc[i, '1'])], int(df.loc[i, 'values'])])
        
    print('\r(4/10) fetching rdk fingerprint pca points...', end='')
    points_rdk_fp_pca = []
    directory = project + '/data/points_rdk_fp_pca.csv'
    df = pd.read_csv(directory)
    for i in range(len(df)):
        points_rdk_fp_pca.append([df.loc[i, 'smiles'], [float(df.loc[i, '0']), float(df.loc[i, '1'])], int(df.loc[i, 'values'])])

    print('\r(5/10) fetching pattern fingerprint tsne points...', end='')
    points_pattern_fp_tsne = []
    directory = project + '/data/points_pattern_fp_tsne.csv'
    df = pd.read_csv(directory)
    for i in range(len(df)):
        points_pattern_fp_tsne.append([df.loc[i, 'smiles'], [float(df.loc[i, '0']), float(df.loc[i, '1'])], int(df.loc[i, 'values'])])
        
    print('\r(6/10) fetching pattern fingerprint pca points...', end='')
    points_pattern_fp_pca = []
    directory = project + '/data/points_pattern_fp_pca.csv'
    df = pd.read_csv(directory)
    for i in range(len(df)):
        points_pattern_fp_pca.append([df.loc[i, 'smiles'], [float(df.loc[i, '0']), float(df.loc[i, '1'])], int(df.loc[i, 'values'])])

    print('\r(7/10) fetching layered fingerprint tsne points...', end='')
    points_layered_fp_tsne = []
    directory = project + '/data/points_layered_fp_tsne.csv'
    df = pd.read_csv(directory)
    for i in range(len(df)):
        points_layered_fp_tsne.append([df.loc[i, 'smiles'], [float(df.loc[i, '0']), float(df.loc[i, '1'])], int(df.loc[i, 'values'])])
        
    print('\r(8/10) fetching layered fingerprint pca points...', end='')
    points_layered_fp_pca = []
    directory = project + '/data/points_layered_fp_pca.csv'
    df = pd.read_csv(directory)
    for i in range(len(df)):
        points_layered_fp_pca.append([df.loc[i, 'smiles'], [float(df.loc[i, '0']), float(df.loc[i, '1'])], int(df.loc[i, 'values'])])

    print('\r(9/10) fetching MACCSKeys fingerprint tsne points...', end='')
    points_MACCSKeys_fp_tsne = []
    directory = project + '/data/points_MACCSKeys_fp_tsne.csv'
    df = pd.read_csv(directory)
    for i in range(len(df)):
        points_MACCSKeys_fp_tsne.append([df.loc[i, 'smiles'], [float(df.loc[i, '0']), float(df.loc[i, '1'])], int(df.loc[i, 'values'])])
        
    print('\r(10/10) fetching MACCSKeys fingerprint pca points...')
    points_MACCSKeys_fp_pca = []
    directory = project + '/data/points_MACCSKeys_fp_pca.csv'
    df = pd.read_csv(directory)
    for i in range(len(df)):
        points_MACCSKeys_fp_pca.append([df.loc[i, 'smiles'], [float(df.loc[i, '0']), float(df.loc[i, '1'])], int(df.loc[i, 'values'])])

    x1 = [ point[1][0] for point in points_latent_pca ]
    y1 = [ point[1][1] for point in points_latent_pca ]

    x2 = [ point[1][0] for point in points_latent_tsne ]
    y2 = [ point[1][1] for point in points_latent_tsne ]

    x3 = [ point[1][0] for point in points_rdk_fp_pca ]
    y3 = [ point[1][1] for point in points_rdk_fp_pca ]

    x4 = [ point[1][0] for point in points_rdk_fp_tsne ]
    y4 = [ point[1][1] for point in points_rdk_fp_tsne ]

    x5 = [ point[1][0] for point in points_pattern_fp_pca ]
    y5 = [ point[1][1] for point in points_pattern_fp_pca ]

    x6 = [ point[1][0] for point in points_pattern_fp_tsne ]
    y6 = [ point[1][1] for point in points_pattern_fp_tsne ]

    x7 = [ point[1][0] for point in points_layered_fp_pca ]
    y7 = [ point[1][1] for point in points_layered_fp_pca ]

    x8 = [ point[1][0] for point in points_layered_fp_tsne ]
    y8 = [ point[1][1] for point in points_layered_fp_tsne ]

    x9 = [ point[1][0] for point in points_MACCSKeys_fp_pca ]
    y9 = [ point[1][1] for point in points_MACCSKeys_fp_pca ]

    x10 = [ point[1][0] for point in points_MACCSKeys_fp_tsne ]
    y10 = [ point[1][1] for point in points_MACCSKeys_fp_tsne ]



    cmap = plt.cm.get_cmap('cool', 2)

    fig, ax = plt.subplots(5, 2)
    fig.set_size_inches((40, 16))
    plt.rcParams['axes.facecolor'] = bg_color
    plt.subplot(251)
    color = [ points[2] for points in points_latent_tsne ]
    plt.scatter(x1, y1, c=color, s=size, alpha=alpha)
    plt.title('1. latent PCA', fontsize=fontsize)

    plt.subplot(256)
    color = [ points[2] for points in points_latent_tsne ]
    plt.scatter(x2, y2, c=color, s=size, alpha=alpha)
    plt.title('2. latent TSNE', fontsize=fontsize)

    plt.subplot(252)
    color = [ points[2] for points in points_rdk_fp_tsne ]
    plt.scatter(x3, y3, c=color, s=size, alpha=alpha)
    plt.title('3. rdk fingerprint PCA', fontsize=fontsize)

    plt.subplot(257)
    color = [ points[2] for points in points_rdk_fp_tsne ]
    plt.scatter(x4, y4, c=color, s=size, alpha=alpha)
    plt.title('4. rdk fingerprint TNSE', fontsize=fontsize)

    plt.subplot(253)
    color = [ points[2] for points in points_pattern_fp_tsne ]
    plt.scatter(x5, y5, c=color, s=size, alpha=alpha)
    plt.title('5. pattern fingerprint PCA', fontsize=fontsize)

    plt.subplot(258)
    color = [ points[2] for points in points_pattern_fp_tsne ]
    plt.scatter(x6, y6, c=color, s=size, alpha=alpha)
    plt.title('6. pattern fingerprint TNSE', fontsize=fontsize)

    plt.subplot(254)
    color = [ points[2] for points in points_layered_fp_tsne ]
    plt.scatter(x7, y7, c=color, s=size, alpha=alpha)
    plt.title('7. layered fingerprint PCA', fontsize=fontsize)

    plt.subplot(259)
    color = [ points[2] for points in points_layered_fp_tsne ]
    plt.scatter(x8, y8, c=color, s=size, alpha=alpha)
    plt.title('8. layered fingerprint TNSE', fontsize=fontsize)

    plt.subplot(255)
    color = [ points[2] for points in points_MACCSKeys_fp_tsne ]
    plt.scatter(x9, y9, c=color, s=size, alpha=alpha)
    plt.title('9. MACCSKeys fingerprint PCA', fontsize=fontsize)

    plt.subplot(2, 5,10)
    color = [ points[2] for points in points_MACCSKeys_fp_tsne ]
    plt.scatter(x10, y10, c=color, s=size, alpha=alpha)
    plt.title('10. MACCSKeys fingerprint TNSE', fontsize=fontsize)

    plt.show()
    

def cluster(project, sel, n_clusters, # cluster, plot, display
            alpha = 0.5,
            size = 5,
            bg_color = '#ffffff'): # '#CAD3C8'):

    # fetch coordinates
    selected = None
    if sel == '' or sel == 1:
        selected = 'latent_group_pca'
    elif sel == 2:
        selected = 'latent_group_tsne'
    elif sel == 3:
        selected = 'rdk_fp_pca'
    elif sel == 4:
        selected = 'rdk_fp_tsne'
    elif sel == 5:
        selected = 'pattern_fp_pca'
    elif sel == 6:
        selected = 'pattern_fp_tsne'
    elif sel == 7:
        selected = 'layered_fp_pca'
    elif sel == 8:
        selected = 'layered_fp_tsne'
    elif sel == 9:
        selected = 'MACCSKeys_fp_pca'
    elif sel == 10:
        selected = 'MACCSKeys_fp_tsne'

    print('fetching', selected, 'points...')
    points_selected = []
    directory = project + '/data/points_' + selected + '.csv'
    df = pd.read_csv(directory)
    for i in range(len(df)):
        points_selected.append([df.loc[i, 'smiles'], [float(df.loc[i, '0']), float(df.loc[i, '1'])], int(df.loc[i, 'values'])])
    df = pd.DataFrame([ points[1] for points in points_selected ]) # 여기부터
    
    # k-means clusters
    kmeans = KMeans(n_clusters=n_clusters)
    kmeans.fit(df)
    labels = kmeans.predict(df)
    x = df[0]
    y = df[1]

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches((5, 5))
    plt.rcParams['axes.facecolor'] = bg_color
    plt.scatter(x, y, c=labels, alpha=alpha, s=size)
    plt.colorbar()
    plt.show()

    # display
    df_area = [ pd.DataFrame([], columns=[i]) for i in range(n_clusters) ]

    for idx, i in enumerate(labels):
        df_area[i].loc[len(df_area[i])] = [Chem.MolFromSmiles(points_selected[idx][0])]
    for label in range(n_clusters):
        print('<< label', label, '>>')
        display(PandasTools.FrameToGridImage(df_area[label][:10], column=label, legendsCol=label, molsPerRow=10))    
    
    return df_area

def export_molecules(project, clustered):
    for idx, group in enumerate(clustered):
        directory = project + '/result' + '{0:0>2}'.format(idx) + '.sd'
        print('\rsaving', directory, '...', end='')
        print(group.values)
        w = SDWriter(directory)   
        for mol in [ mol[0] for mol in group.values ]:
            w.write(mol)
    print('files are saved!')
