# 1. median-center
meds   <- apply(lipid_mat, 2, median, na.rm=TRUE)
med_mat <- sweep(lipid_mat, 2, meds, `-`)

# 2. scale by MAD
mads    <- apply(lipid_mat, 2, mad,    na.rm=TRUE)
robust_z <- sweep(med_mat, 2, mads, `/`)

# 3. UMAP on this
set.seed(42)
umap_robust <- umap(robust_z, n_neighbors=15, min_dist=0.1)
quartz()
plot(umap_robust, col = as.factor(batch_vec), pch=19,
     main="UMAP on median/MAD-scaled data")
lipid_mat <- robust_z

