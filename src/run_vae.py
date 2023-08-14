#!/usr/bin/env python3

import sys, getopt
sys.path.append("/groups/hongxuding/software/Python_Lib")
import torch
from vae import VAE
import numpy as np
import umap

job = ''
mode = ''
data_path = ''
gene_path = ''
tf_path = ''
model_path = ''
out_dir = ''
depth = ''
random_seed=1
opt, args = getopt.getopt(sys.argv[1:], 'j:m:i:g:t:m:o:n:s', ['job=',
                                                              'mode=',
                                                              'data_path=',
                                                              'gene_path=',
                                                              'tf_path=',
                                                              'model_path=',
                                                              'out_dir=',
                                                              'depth=',
                                                              'random_seed='])

for o, a in opt:
    if o in ('-j', '--job'):
        job = a
    if o in ('-m', '--mode'):
        mode = a
    if o in ('-i', '--data_path'):
        data_path = a
    elif o in ('-g', '--gene_path'):
        gene_path = a
    elif o in ('-t', '--tf_path'):
        tf_path = a
    elif o in ('-m', '--model_path'):
        model_path = a
    elif o in ('-o', '--out_dir'):
        out_dir = a
    elif o in ('-n', '--depth'):
        depth = int(a)
    elif o in ('-s', '--random_seed'):
        random_seed = int(a)

print("-" * 50)
print("Job Information for VAE Running")
print("-" * 50)
print("job name            : %s" % job)
print("mode                : %s" % mode)
print("tab-delimited data  : %s" % data_path)
print("gene list           : %s" % gene_path)
print("TF list             : %s" % tf_path)
print("VAE model           : %s" % model_path)
print("output directory    : %s" % out_dir)
print("number of layers    : %s" % depth)
print("seed                : %s" % random_seed)
print("\n")

# Set seed
torch.manual_seed(random_seed)
np.random.seed(random_seed)

# load model
with open(gene_path) as f:
    genes = f.read().splitlines()
with open(tf_path) as f:
    tfs = f.read().splitlines()
mask = np.isin(np.array(genes), np.array(tfs))
vae = VAE(mask_tf=mask, depth=depth)
vae.load_state_dict(torch.load(model_path, map_location=torch.device('cpu')))
vae.eval()

# Run VAE
run_data = np.genfromtxt(data_path, delimiter = " ")
if mode == 'generate':
    run_data_rec = vae.decode(torch.Tensor(run_data))
    run_data_rec = run_data_rec.detach().numpy()
    
    reducer = umap.UMAP(random_state=random_seed, min_dist=0.5, n_neighbors=15)
    umap_x_rec = reducer.fit_transform(run_data_rec)
    
    np.savetxt(out_dir+'/'+job+'_'+mode+'_x_rec_data.csv.gz', run_data_rec, delimiter=" ")
    np.savetxt(out_dir+'/'+job+'_'+mode+'_x_rec_umap.csv.gz', umap_x_rec, delimiter=" ")

elif mode == 'predict':
    run_data_rec, z, mu, logvar = vae.forward(torch.Tensor(run_data))
    run_data_rec = run_data_rec.detach().numpy()
    z = z.detach().numpy()
    
    reducer = umap.UMAP(random_state=random_seed, min_dist=0.5, n_neighbors=15)
    umap_x = reducer.fit_transform(run_data)
    umap_x_rec = reducer.fit_transform(run_data_rec)
    
    np.savetxt(out_dir+'/'+job+'_'+mode+'_x_rec_data.csv.gz', run_data_rec, delimiter=" ")
    np.savetxt(out_dir+'/'+job+'_'+mode+'_z_reparametrization.csv.gz', z, delimiter=" ")
    np.savetxt(out_dir+'/'+job+'_'+mode+'_x_umap.csv.gz', umap_x, delimiter=" ")
    np.savetxt(out_dir+'/'+job+'_'+mode+'_x_rec_umap.csv.gz', umap_x_rec, delimiter=" ")

else:
    print("mode should be generate or predict")
