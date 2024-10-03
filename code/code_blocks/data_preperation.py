# read and save subject data as CSV
file_1 = pd.read_csv('GSM3746212_Muscle_1_Counts.csv', index_col=0)
file_2 = pd.read_csv('GSM3746213_Muscle_2_Counts.csv', index_col=0)
file_3 = pd.read_csv('GSM3746214_Muscle_3_Counts.csv', index_col=0)
file_4 = pd.read_csv('GSM3746215_Muscle_4_Counts.csv', index_col=0)

# convert CSV files into standalone AnnData objects
muscle_1 = sc.AnnData(file_1)
muscle_2 = sc.AnnData(file_2)
muscle_3 = sc.AnnData(file_3)
muscle_4 = sc.AnnData(file_4)

# initialize list of AnnData objects, storing ea/ subjects data
adatas = [muscle_1, muscle_2, muscle_3, muscle_4]

# create composite AnnData object
adata_combined = sc.concat(adatas, axis=1, label='sample', keys=['muscle_1', 'muscle_2', 'muscle_3', 'muscle_4'])
adata_combined.var_names_make_unique()

# transpose AnnData object
adata_transposed = sc.AnnData(adata_combined.T)

# filter out cells that have fewer than 200 detected genes
sc.pp.filter_cells(adata_transposed, min_genes=200)

# filter out genes that appear in fewer than 20 cells
sc.pp.filter_genes(adata_transposed, min_cells=20) 

# print resulting number of cells and genes
num_genes = adata_transposed.n_vars
print(f"Number of Genes: {num_genes}")

num_cells = adata_transposed.n_obs
print(f"Number of cells: {num_cells}")

# apply global-scaling normalization
sc.pp.normalize_total(adata_transposed, target_sum=1e4)
sc.pp.log1p(adata_transposed)

# find the 2000 most highly variable genes 
sc.pp.highly_variable_genes(adata_transposed, n_top_genes=2000, subset=True)
print(adata_transposed)

# apply z-transformation
sc.pp.scale(adata_transposed, zero_center=True)

# perform dimensionality reduction via PCA
sc.tl.pca(adata_transposed, svd_solver='arpack')
