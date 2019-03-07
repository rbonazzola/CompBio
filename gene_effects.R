# Estimation of gene-level effects, using GWAS and eQTL data

library(rstan)
options(mc.cores = parallel::detectCores())
library(dplyr)

# Estimating the variance of gene-level effects on height, using a GWAS on height and eQTL data on Whole Blood
df <- readRDS("Genomica/misc/CompBio/DAPG_with_ldscore_by_phenotype_Whole_Blood__UKB_50_Standing_height.rds")

which_genes <- sample(1:length(unique(df$gene_id)), 2000)
df <- df %>% filter(as.integer(gene_id) %in% which_genes)

# preparing the data for Stan
stan_data <- list(beta_gwas=df$beta_gwas,
                  beta_eqtl=df$beta_eqtl,
                  N = nrow(df),
                  N_genes= length(unique(df$gene_id)),
                  gene=as.integer(factor(df$gene_id)),
                  ldscore=df$ldscore)


stan_model <- "data{
  int<lower=1> N;
  int<lower=1> N_genes;
  int<lower=1> gene[N];
  real beta_gwas[N];
  real beta_eqtl[N];
  real ldscore[N];
}
parameters{
  real<lower=0> beta0; // beta_gwas is forced to be positive
  vector[N_genes] beta_gene;
  real<lower=0> sigma_g;
  real<lower=0> sigma_e;
  real ldscore_slope;
}
model{
  real mu;
  beta_gene ~ normal(0, sigma_g); // prior on beta_gene
  for (i in 1:N) {
    mu = beta0 + beta_gene[gene[i]] * beta_eqtl[i] + ldscore_slope * ldscore[i];
    beta_gwas[i] ~ normal(mu, sigma_e);
  }
}
"

stan_fit <- stan(model_code = stan_model, data = stan_data, iter=2000, chains=4)
pp <- traceplot(stan_fit_2, pars=c("beta0", "sigma_g", "sigma_e", "ldscore_slope"), inc_warmup=F)
