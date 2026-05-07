import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import nbinom, poisson
from sklearn.decomposition import PCA


#%%
#import data
data = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19', var_names='gene_symbols')

#%%visualize highest expressed genes
sc.pl.highest_expr_genes(data, n_top=20)
sc.pp.calculate_qc_metrics(data, inplace=True)

#%%visualize counts and genes distributions
plt.hist(data.obs['total_counts'], bins=100)
plt.title("Total counts per cell")
plt.show()

plt.hist(data.obs['n_genes_by_counts'], bins=100)
plt.title("Genes per cell")
plt.show()

plt.hist(data.var['n_cells_by_counts'], bins=100)
plt.title("Cells per gene")
plt.show()

#%% filter data and revisualize

sc.pp.filter_genes(data, min_cells= 5)
sc.pp.filter_cells(data, min_genes=200)

sc.pp.calculate_qc_metrics(data, inplace=True)

plt.hist(data.obs['n_genes_by_counts'], bins=100)
plt.title("Genes per cell after filtering")
plt.show()

plt.hist(data.var['n_cells_by_counts'], bins=100)
plt.title("Cells per gene after filtering")
plt.show()

#%% perform total count normalization

df_raw = data.to_df()
sc.pp.normalize_total(data, target_sum=1e4)
data.layers['normalized'] = data.X.copy()


#%% visualize mean-variance relationship

df = data.to_df('normalized')

#use unfiltered mean and variance to plot relationship
gene_means = df.mean(axis=0)
gene_vars = df.var(axis=0)


plt.scatter(gene_means,gene_vars, alpha=0.3, s=5)
plt.xlabel("Mean expression")
plt.ylabel("Variance")
plt.title("Mean-Variance Relationship")
plt.show()
## indicates that as mean increases variance increase, poission assumption violated

plt.scatter(gene_means,gene_vars, alpha=0.3, s=5)
plt.xscale('log')
plt.yscale('log')
plt.xlabel("Log(Mean expression)")
plt.ylabel("Log(Variance)")
plt.title("Mean-Variance Relationship")
plt.show()
# increased variance at lower mean, meaning lowly expressed genes have noisy data

#%% estimate and visualize dispersion

#only consider over dispersed genes
valid = gene_vars > gene_means
df_valid = df.loc[:,valid]
gene_means_v = gene_means[valid]
gene_vars_v  = gene_vars[valid]

#compute gene-wise θ using MOM
theta_genewise = gene_means_v**2 / (gene_vars_v - gene_means_v)

plt.scatter(gene_means_v,(1/theta_genewise), alpha=0.3, s=5)
plt.xscale('log')
plt.yscale('log')
plt.xlabel("Log(Mean expression)")
plt.ylabel("Log(1/θ = overdispersion)")
plt.title("Overdispersion vs Mean Expression")
plt.show()

# decreased over dispersion at larger expression, lowly expressed genes once
# show noise and are over dispersed

#%% fit θ estimates

# --- global θ ---
theta_global = np.median(theta_genewise)

# --- binned θ  ---

# take log so bins are more balanced
log_means_v = np.log1p(gene_means_v)
n_bins = 50
bins = pd.qcut(log_means_v, q=n_bins, duplicates='drop')
theta_binned = bins.map(theta_genewise.groupby(bins).median()).astype(float)

# --- LOESS regularization ---
from statsmodels.nonparametric.smoothers_lowess import lowess
log_theta = np.log(theta_genewise)
smoothed = lowess(log_theta, log_means_v, frac=0.3)
# interpolate back to all genes
theta_loess = np.exp(np.interp(log_means_v, smoothed[:, 0], smoothed[:, 1]))


#%% plot estimates

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

x = gene_means_v

for ax, y, title in zip(
    axes,
    [theta_genewise, 1/theta_genewise],
    ["θ vs Mean", "1/θ (overdispersion) vs Mean"]
):
    ax.scatter(x, y, alpha=0.2, s=3, label='Gene-wise θ', color='grey')
    ax.axhline(theta_global if 'θ vs' in title else 1/theta_global, 
               color='red', lw=2, label='Global θ')
    ax.scatter(x, theta_binned if 'θ vs' in title else 1/theta_binned, 
               alpha=0.4, s=3, color='orange', label='Binned θ')
    ax.scatter(x, theta_loess if 'θ vs' in title else 1/theta_loess, 
               alpha=0.4, s=3, color='green', label='LOESS θ')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel("Mean expression"); ax.set_ylabel(title.split(' vs')[0])
    ax.set_title(title); ax.legend(markerscale=3)

plt.tight_layout()
plt.show()

#%% compute pearson residuals
def pearson_residuals(counts_df, mu_df, theta):
    mu = mu_df
    var = mu + mu**2 / theta
    return (counts_df - mu) / np.sqrt(var)

mu_df = gene_means_v

resid_poisson   = (df_valid - mu_df) / np.sqrt(mu_df)
resid_global    = pearson_residuals(df_valid, mu_df, theta_global)
resid_genewise  = pearson_residuals(df_valid, mu_df, theta_genewise)
resid_binned    = pearson_residuals(df_valid, mu_df, theta_binned)
resid_loess     = pearson_residuals(df_valid, mu_df, theta_loess)

models = {
    'Poisson':      {'resid': resid_poisson,  'theta': None},
    'Global θ':     {'resid': resid_global,   'theta': pd.Series(theta_global, index=df_valid.columns)},
    'Gene-wise θ':  {'resid': resid_genewise, 'theta': theta_genewise},
    'Binned θ':     {'resid': resid_binned,   'theta': theta_binned},
    'LOESS θ':      {'resid': resid_loess,    'theta': pd.Series(theta_loess, index=df_valid.columns)},
}
n_models     = len(models)


#%% goodness of fit — pearson chi-squared

def pearson_chisq_pergene(counts_df, mu, theta):
    """chi-squared statistic summed over cells, per gene"""
    var = mu + mu**2 / theta
    return ((counts_df - mu)**2 / var).sum(axis=0)

n_cells = df_valid.shape[0]

# compute per-gene chi-squared for each model
chisq = {}
for name, m in models.items():
    if name == 'Poisson':
        var = mu_df  # poisson: var = mean
        chisq[name] = ((df_valid - mu_df)**2 / var).sum(axis=0)
    else:
        chisq[name] = pearson_chisq_pergene(df_valid, mu_df, m['theta'])

# normalize by degrees of freedom
chisq_df = pd.DataFrame({name: vals / (n_cells - 1) 
                          for name, vals in chisq.items()})

# print summary statistics
summary = pd.DataFrame({
    'Median': chisq_df.median(),
    'Mean':   chisq_df.mean(),
    'Std':    chisq_df.std(),
    '90th percentile': chisq_df.quantile(0.90),
    '% genes > 2': (chisq_df > 2).mean() * 100
}).round(3)
print(summary)

# plot
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

colors = ['steelblue', 'red', 'grey', 'orange', 'green']

# left panel includes all models
for (name, vals), color in zip(chisq.items(), colors):
    axes[0].hist((vals / (n_cells - 1)).clip(0, 20), 
                 bins=80, alpha=0.5, density=True, 
                 label=name, color=color)
axes[0].axvline(1, color='black', linestyle='--', label='ideal (=1)')
axes[0].set_yscale('log')
axes[0].set_xlabel("Pearson χ²/df")
axes[0].set_ylabel("Density (log scale)")
axes[0].set_title("All Models")
axes[0].legend()

# include NB models only for right panels
for (name, vals), color in zip(list(chisq.items())[1:], colors[1:]):
    axes[1].hist((vals / (n_cells - 1)).clip(0, 5), 
                 bins=80, alpha=0.5, density=True, 
                 label=name, color=color)
axes[1].axvline(1, color='black', linestyle='--', label='ideal (=1)')
axes[1].set_yscale('log')
axes[1].set_xlabel("Pearson χ²/df")
axes[1].set_title("NB Models Only")
axes[1].legend()

plt.suptitle("Goodness of Fit: Pearson χ² per Gene")
plt.tight_layout()
plt.show()
#%% plot variance of resiuals vs mean
fig, axes = plt.subplots(1, n_models, figsize=(4 * n_models, 4), sharey=True)
for ax, (name, m) in zip(axes, models.items()):
    resid_var = m['resid'].var(axis=0)
    ax.scatter(gene_means_v, resid_var, alpha=0.3, s=3)
    ax.axhline(1, color='red', linestyle='--', label='ideal (=1)')
    ax.set_xscale('log')
    ax.set_xlabel("Mean expression")
    ax.set_title(name)
    ax.legend()
axes[0].set_ylabel("Variance of Pearson residuals")
plt.suptitle("Variance Stabilization by Model", y=1.02)
plt.tight_layout()
plt.show()

#%% fit vs observed distribution

# plot one highly expressed and one lowly expressed gene
genes_to_plot = ['MALAT1', 'RPS23']

fig, axes = plt.subplots(len(genes_to_plot), n_models,
                          figsize=(4 * n_models, 3.5 * len(genes_to_plot)))

for i, gene in enumerate(genes_to_plot):
    #use raw counts
    counts = df_raw[gene].values.astype(int)
    mu = counts.mean()
    x = np.arange(0, np.percentile(counts, 99).astype(int))  # plot up to 99th percentile

    for j, (model_name, m) in enumerate(models.items()):
        ax = axes[i, j]

        # plot histogram
        ax.hist(counts, bins=40, density=True, alpha=0.6,
                color='lightgrey', edgecolor='white')

        # plot distribution
        if model_name == 'Poisson':
            ax.plot(x, poisson.pmf(x, mu), 'r-', lw=2)
        else:
            theta = float(m['theta'][gene])
            p = theta / (theta + mu)
            ax.plot(x, nbinom.pmf(x, theta, p), 'r-', lw=2)

        ax.set_title(f"{gene}\n{model_name}", fontsize=9)
        ax.set_xlabel("Counts")
        ax.set_ylabel("Density")
        ax.set_xlim(0, x.max())  # same x-axis per gene

plt.suptitle("NB Model Fit vs Observed Count Distribution", y=1.02)
plt.tight_layout()
plt.show()

#%% PCA

#get cell type labels
sc.pp.neighbors(data)
sc.tl.leiden(data, resolution=0.5, flavor='igraph', n_iterations=2, directed=False)
labels = data.obs['leiden']

fig, axes = plt.subplots(1, n_models, figsize=(4.5 * n_models, 4))

for ax, (name, m) in zip(axes, models.items()):
    X = m['resid'].clip(-10, 10).fillna(0).values
    coords = PCA(n_components=2).fit_transform(X)

    for cluster in sorted(labels.unique(), key=int):
        mask = labels == cluster
        ax.scatter(coords[mask, 0], coords[mask, 1], s=3, alpha=0.5, label=cluster)

    var1, var2 = PCA(n_components=2).fit(X).explained_variance_ratio_[:2] * 100
    ax.set_title(f"{name}\nPC1={var1:.1f}% PC2={var2:.1f}%")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")

axes[-1].legend(title="Cluster", markerscale=3,
                bbox_to_anchor=(1.05, 1), loc='upper left')
plt.suptitle("PCA by Model", y=1.02)
plt.tight_layout()
plt.show()







