import os
from implements.encoder import encode_latent_fp
from implements.learner import learn_tsne_pca
from implements.plotter import plot_molecules, cluster, export_molecules

# Setup Project
project = input('input project name ')
directory = './' + project
try:
    if not os.path.exists(directory):
        os.makedirs(directory)
        os.makedirs(directory + '/data')
        print('project created :', project)
    else:
        print('project with the same name already exists :', project)
except OSError:
    print('Error : Creating directory', directory)


## Encode molecules to latent space & fingerprints
encode_latent_fp(project)

## tSNE & PCA
learn_tsne_pca('test')

## Plot molecules
plot_molecules('test', size=200, alpha=1, fontsize=12)

## Select and cluster the molecules
sel = input('enter the graph to be clustered')
sel = 1 if sel == '' else int(sel)
n_clusters = input('enter the number of group')
n_clusters = 2 if n_clusters == '' else int(n_clusters)
print('input:', sel, n_clusters)
clustered = cluster('test', sel, n_clusters, size=100, alpha=1) # cluster, plot, display

## Export the molecules
export_molecules('test', clustered)