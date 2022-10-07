import numpy as np
from scipy.special import stdtr
import pandas as pd
from sklearn.decomposition import PCA

def load_screens():
    # Load screens
    screens = pd.read_csv('data/gene_effect.csv', index_col=0).T
    screens.index = screens.index.str.split(' ').str.get(0)
    # Map Broad ID to CCLE name
    cell_lines = pd.read_csv('data/sample_info.csv', index_col='DepMap_ID',
                             usecols=['DepMap_ID', 'CCLE_name'], squeeze=True)
    screens.columns = screens.columns.reindex(cell_lines)[0]
    # Bias-correct using "molecular function:olfactory receptor activity" genes
    olfactory_genes = pd.read_csv(
        'data/olfactory_genes.txt', header=None, squeeze=True)
    olfactory_data = screens.reindex(olfactory_genes).dropna()
    transformation = PCA(n_components=4)
    transformation.fit(olfactory_data)
    top_PC_effects = transformation.inverse_transform(
        transformation.transform(screens))
    screens -= top_PC_effects
    screens = screens.iloc[:, :-4]
    return screens


# Load batch-corrected screens

screens = load_screens()

# Remove cell lines with any missing genes
# (not required for DepMap 18Q3, but is for more recent releases)
# You can use other strategies to remove NaNs instead, like imputing,
# removing genes with any missing cell lines

screens.dropna(axis=1, inplace=True)

# Warp screen data and intercept based on covariance of screens

cholsigmainv = np.linalg.cholesky(np.linalg.inv(np.cov(screens.T)))
warped_screens = screens.values @ cholsigmainv
warped_intercept = cholsigmainv.sum(axis=0)

# Then just run linear regression; this implementation is based on 
# https://pingouin-stats.org/generated/pingouin.linear_regression.html

def linear_regression(warped_screens, warped_intercept):
    GLS_coef = np.empty((len(warped_screens), len(warped_screens)))
    GLS_se = np.empty((len(warped_screens), len(warped_screens)))
    ys = warped_screens.T
    for gene_index in range(len(warped_screens)):
        X = np.stack((warped_intercept, warped_screens[gene_index]), axis=1)
        coef, residues = np.linalg.lstsq(X, ys, rcond=None)[:2]
        df = warped_screens.shape[1] - 2
        GLS_coef[gene_index] = coef[1]
        GLS_se[gene_index] = \
            np.sqrt(np.linalg.pinv(X.T @ X)[1, 1] * residues / df)
    return GLS_coef, GLS_se

GLS_coef, GLS_se = linear_regression(warped_screens, warped_intercept)
df = warped_screens.shape[1] - 2
GLS_p = 2 * stdtr(df, -np.abs(GLS_coef / GLS_se))
np.fill_diagonal(GLS_p, 1)

# Save everything

np.save('data/GLS_p.npy', GLS_p)
np.save('data/GLS_sign.npy', np.sign(GLS_coef))
screens.index.to_series().to_csv('data/genes.txt', index=False, header=False)
