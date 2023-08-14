#!/usr/bin/env python3

import sys, getopt
sys.path.append("/groups/hongxuding/software/Python_Lib")
import torch
from vae import VAE
from utils import UnsupervisedDataset, KFoldTorch
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.cm as cm

job = ''
data_path = ''
gene_path = ''
tf_path = ''
out_dir = ''
depth = 3
beta = 0.00005
alpha = 0.5
cv = 10
lr = 1e-3
random_seed=1
opt, args = getopt.getopt(sys.argv[1:], 'j:i:g:t:o:n:b:a:f:r:s', ['job=',
                                                                  'data_path=',
                                                                  'gene_path=',
                                                                  'tf_path=',
                                                                  'out_dir=',
                                                                  'depth=',
                                                                  'beta=',
                                                                  'alpha=',
                                                                  'cv=',
                                                                  'lr=',
                                                                  'random_seed='])

for o, a in opt:
    if o in ('-j', '--job'):
        job = a
    if o in ('-i', '--data_path'):
        data_path = a
    elif o in ('-g', '--gene_path'):
        gene_path = a
    elif o in ('-t', '--tf_path'):
        tf_path = a
    elif o in ('-o', '--out_dir'):
        out_dir = a
    elif o in ('-n', '--depth'):
        depth = int(a)
    elif o in ('-b', '--beta'):
        beta = float(a)
    elif o in ('-a', '--alpha'):
        alpha = float(a)
    elif o in ('-f', '--cv'):
        cv = int(a)
    elif o in ('-r', '--lr'):
        lr = float(a)
    elif o in ('-s', '--random_seed'):
        random_seed = int(a)

print("-" * 50)
print("Job Information for VAE Training")
print("-" * 50)
print("job name            : %s" % job)
print("tab-delimited data  : %s" % data_path)
print("gene list           : %s" % gene_path)
print("TF list             : %s" % tf_path)
print("output directory    : %s" % out_dir)
print("number of layers    : %s" % depth)
print("beta                : %s" % beta)
print("alpha               : %s" % alpha)
print("cv folds            : %s" % cv)
print("learning rate       : %s" % lr)
print("seed                : %s" % random_seed)
print("\n")

# Set model
torch.backends.cudnn.enabled = True
torch.manual_seed(random_seed)
dev = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print('Using device:', dev, flush=True)
    
# Read tf and data
with open(gene_path) as f:
    genes = f.read().splitlines()
with open(tf_path) as f:
    tfs = f.read().splitlines()
mask = np.isin(np.array(genes), np.array(tfs))
train_data = np.genfromtxt(data_path, delimiter = " ")
train_ds = torch.Tensor(train_data)
train_ds = UnsupervisedDataset(train_ds, targets=torch.zeros(train_data.shape[0]))
    
# Initialize CV
kfold = KFoldTorch(cv=cv, n_epochs=2000, lr=lr, train_p=25, test_p=25, num_workers=0, save_all=True, save_best=False, path_dir=out_dir, model_prefix=job)
dict_params = {'mask_tf':mask, 'depth':depth, 'init_w':True, 'beta':beta, 'alpha':alpha, 'dropout':0.3, 'path_model':None, 'device':dev}
kfold.train_kfold(VAE, dict_params, train_ds, batch_size=64)
np.save(out_dir+'/'+job+'_'+str(cv)+'CV_vae.npy', kfold.cv_res_dict)

# Train history

res = np.load(out_dir+'/'+job+'_'+str(cv)+'CV_vae.npy', allow_pickle=True).item()
colors = cm.jet(np.linspace(0, 1, cv))
for f in range(0, cv):
    loss = res[f].get('history').get('train_loss')
    epoch = [i for i in range(0, len(loss))]
    plt.plot(epoch, loss, 'o', color=colors[f])
    loss = res[f].get('history').get('valid_loss')
    epoch = [i for i in range(0, len(loss))]
    plt.plot(epoch, loss, 'x', color=colors[f])
plt.xlabel("epoch")
plt.ylabel("loss")
plt.legend(handles=[mlines.Line2D([], [], marker='o', linestyle='None', color='black', label='train'),
                    mlines.Line2D([], [], marker='x', linestyle='None', color='black', label='valid')])
plt.savefig(out_dir+'/history.pdf') 
