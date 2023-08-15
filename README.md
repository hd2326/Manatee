# Manatee

**Overview**

Manatee variational autoencoder (VAE) model predicts transcription factor (TF) perturbation-induced transcriptomes. The workflow of Manatee is shown in the following schematics. Specifically:

During training, the Manatee latent space is constrained to approximate TF expression with an additional TF loss term.

The Manatee _in silico_ perturbation prediction start with adjusting the latent space: for TFs to be up- and down-regulated, corresponding latent values are sampled from the top and bottom Q quantile of the given reference, respectively. Adjusted profiles are subsequently flowed through the trained decoder for final results.

To assess if a perturbation yields the desired transcriptomic phenotype, the relative strength among _r(Perturb, End)_, _r(Start, End)_ and _r(Perturb, Start)_ Pearson correlations, quantified by t-statistics, is calculated. The assessment plot is used to distinguish ideal, compromised (under) and non-significant (N.S.) groups. _Start_, _End_ and _Perturb_ represent the original, target and perturbed groups, respectively

![Manatee](https://github.com/hd2326/Manatee/blob/main/images/manatee.png)

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
