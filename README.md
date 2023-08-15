# Manatee

**Overview**

Manatee variational autoencoder (VAE) model predicts transcription factor (TF) perturbation-induced transcriptomes. The workflow of Manatee is shown in the following schematics. Specifically:

During training, the Manatee latent space is constrained to approximate TF expression with an additional TF loss term.

The Manatee _in silico_ perturbation prediction start with adjusting the latent space: for TFs to be up- and down-regulated, corresponding latent values are sampled from the top and bottom Q quantile of the given reference, respectively. Adjusted profiles are subsequently flowed through the trained decoder for final results.

To assess if a perturbation yields the desired transcriptomic phenotype, the relative strength among _r(Perturb, End)_, _r(Start, End)_ and _r(Perturb, Start)_ Pearson correlations, quantified by t-statistics, is calculated. The assessment plot is used to distinguish ideal, compromised (under) and non-significant (N.S.) groups. _Start_, _End_ and _Perturb_ represent the original, target and perturbed groups, respectively

![Manatee](https://github.com/hd2326/Manatee/blob/main/images/manatee.png)

We use the [GSE72857](https://pubmed.ncbi.nlm.nih.gov/26627738/) dataset as an example to show the Manatee functions.

**Model Training:**

Assuming we are at the root directory of the Manatee repository. We could train Manatee with the following code: 

```
vae=./src/train_vae.py
job=GSE72857
data_path=./GSE72857/processed/data_x.csv.gz
gene_path=./GSE72857/processed/genes.txt
tf_path=./GSE72857/processed/tfs.txt
out_dir=./GSE72857/model/

python3 $vae --job=$job --data_path=$data_path --gene_path=$gene_path --tf_path=$tf_path --out_dir=$out_dir --depth=3
```

A pre-trained model can be found at ```./GSE72857/model/```

**Model Benchmarking:**

We examined whether Manatee is able to capture biological information, by benchmarking its two modes.

The “predict” mode encodes original transcriptomes (X), reparametrizes the latent space as decoder inputs Z*, and reconstructs X’.

The “generate” mode, on the other hand, directly decodes X’ from Z*.

We could benchmark Manatee with the following code:

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

And here are the results, which show that both modes could recapitulate original transcriptomes (X'=X):

![benchmark](https://github.com/hd2326/Manatee/blob/main/images/benchmark.png)


