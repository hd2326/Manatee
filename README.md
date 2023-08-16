# Manatee: a Variational Autoencoder Model for Predicting Transcription Factor Perturbation-Induced Transcriptomes

## Overview

Manatee variational autoencoder (VAE) model predicts transcription factor (TF) perturbation-induced transcriptomes. The workflow of Manatee is shown in the following schematics. Specifically:

During training, the Manatee latent space is constrained to approximate TF expression with an additional TF loss term.

The Manatee _in silico_ perturbation prediction start with adjusting the latent space. For TFs to be up- and down-regulated, corresponding latent values are sampled from the top and bottom Q quantile of the given reference, respectively. Adjusted profiles are subsequently flowed through the trained decoder for final results.

To assess if a perturbation yields the desired transcriptomic phenotype, the relative strength among _r(Perturb, End)_, _r(Start, End)_ and _r(Perturb, Start)_ Pearson correlations, quantified by t-statistics, is calculated. The assessment plot is used to distinguish ideal, compromised (under) and non-significant (N.S.) groups. _Start_, _End_ and _Perturb_ represent the original, target and perturbed groups, respectively

![Manatee](https://github.com/hd2326/Manatee/blob/main/images/manatee.png)

## Data Processing

We use the [GSE72857](https://pubmed.ncbi.nlm.nih.gov/26627738/) dataset as an example to show the Manatee functions.

Assuming we are at the root directory of the Manatee repository. The rawdata is in ```./GSE72857/processed/```. Further data processing for sections **_In Silico_ Perturbation** and **_In Silico_ Screening** can be done by running ```processed.R```. The yielded data files are provided in ```./GSE72857/perturb/``` and ```./GSE72857/screen/``` respectively for reproducing our results.

## Model Training:

```
vae=./src/train_vae.py
job=GSE72857
data_path=./GSE72857/processed/data_x.csv.gz
gene_path=./GSE72857/processed/genes.txt
tf_path=./GSE72857/processed/tfs.txt
out_dir=./GSE72857/model/

python3 $vae --job=$job --data_path=$data_path --gene_path=$gene_path --tf_path=$tf_path --out_dir=$out_dir --cv=5 --lr=1e-4 --depth=3 --alpha=0.8
```

A pre-trained model can be found [here](https://visualify.pharmacy.arizona.edu/Manatee/GSE72857best_fold.pt). You will need to download it into ```./GSE72857/model/``` to run the following code chunks.

## Model Benchmarking:

We examine whether Manatee is able to capture biological information, by benchmarking its two modes:

The “predict” mode encodes original transcriptomes (X), reparametrizes the latent space as decoder inputs Z*, and reconstructs X’.

The “generate” mode, on the other hand, directly decodes X’ from Z*.

```
vae=./src/train_vae.py
gene_path=./GSE72857/processed/genes.txt
tf_path=./GSE72857/processed/tfs.txt
model_path=./GSE72857/model/GSE72857best_fold.pt
x_path=./GSE72857/processed/data_x.csv.gz
z_path=./GSE72857/processed/data_z.csv.gz
out_dir=./GSE72857/benchmark/

python3 $vae --job=basic --mode=predict --data_path=$x_path --gene_path=$gene_path --tf_path=$tf_path --model_path=$model_path --out_dir=$out_dir --depth=3
python3 $vae --job=basic --mode=generate --data_path=$z_path --gene_path=$gene_path --tf_path=$tf_path --model_path=$model_path --out_dir=$out_dir --depth=3
```

The result shows that both modes could recapitulate original transcriptomes (X'=X):

![benchmark](https://github.com/hd2326/Manatee/blob/main/images/benchmark.png)

## _In Silico_ Perturbation

We use Manatee to model the hematopoietic CMP to GMP VS MEP development, which is driven by the _Gata1-Spi1_ TF module:

```
vae=./src/train_vae.py
gene_path=./GSE72857/processed/genes.txt
tf_path=./GSE72857/processed/tfs.txt
model_path=./GSE72857/model/GSE72857best_fold.pt
z_path=./GSE72857/perturb/z_perturb.csv.gz
out_dir=./GSE72857/perturb/

python3 $vae --job=perturb --mode=generate --data_path=$z_path --gene_path=$gene_path --tf_path=$tf_path --model_path=$model_path --out_dir=$out_dir --depth=3
```
The result shows the recapitulation of the lineage bifurcation:

![perturb](https://github.com/hd2326/Manatee/blob/main/images/perturb.png)

## _In Silico_ Screening

We leverage Manatee for the _in silico_ screening of perturbations that could yield the target transcriptomic phenotype. As a proof-of-concept, we screened alternative TF duos driving the above hematopoiesis process. Without losing generality, we screened all 25 combinations of the top 5 highly expressed TFs within MEP and GMP:

```
vae=./src/train_vae.py
gene_path=./GSE72857/processed/genes.txt
tf_path=./GSE72857/processed/tfs.txt
model_path=./GSE72857/model/GSE72857best_fold.pt
z_path=./GSE72857/screen/z_screen.csv.gz
out_dir=./GSE72857/screen/

python3 $vae --job=screen --mode=generate --data_path=$z_path --gene_path=$gene_path --tf_path=$tf_path --model_path=$model_path --out_dir=$out_dir --depth=3
```

The result shows the _Gata1-Cebpa_ duo to be an alternative lineage controlling module. Noticeably, _Cebpa_ has been experimentally validated as a master regulator controlling the identity of GMP:

![screen](https://github.com/hd2326/Manatee/blob/main/images/screen.png)

## A Comprehensive Model for Mouse Development

We speculate the great value of “comprehensive models”, which represent all major biological processes related to a certain species. We thus provide [a Manatee model](https://visualify.pharmacy.arizona.edu/Manatee/TOMEbest_fold.pt) trained from the [TOME](http://tome.gs.washington.edu/) mouse embryogenesis single cell dataset, which contains 468 cell populations and 113 lineages. Such a model can be used for mouse _in silico_ perturbation and _in silico_ screening analyses, as well as for further algorithmic developments.
