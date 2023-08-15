# Manatee

Manatee variational autoencoder (VAE) model predicts transcription factor (TF) perturbation-induced transcriptomes.

We use the [GSE72857](https://pubmed.ncbi.nlm.nih.gov/26627738/) dataset as an example to show the Manatee functions.

**Directory Tree**

./GSE72857/
|-- benchmark
|-- model
|-- perturb
|-- processed
|-- random
`-- screen

**Model Training:**
```
vae=./src/train_vae.py
job=GSE72857
data_path=./GSE72857/processed/data_x.csv.gz
gene_path=./GSE72857/processed/genes.txt
tf_path=./GSE72857/processed/tfs.txt
out_dir=./GSE72857/model/

python3 $vae --job=$job --data_path=$data_path --gene_path=$gene_path --tf_path=$tf_path --out_dir=$out_dir --cv=5 --lr=1e-4 --depth=3 --alpha=0.8
```
**Model Benchmarking:**

The "predict" mode reconstruct transcriptomes (X')
