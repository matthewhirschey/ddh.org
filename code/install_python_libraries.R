# create a new environment
reticulate::virtualenv_create("gls_regression")

# install requirements for by code/generate_GLS.R
reticulate::virtualenv_install("gls_regression", c("numpy","scipy", "pandas", "scikit-learn"))
