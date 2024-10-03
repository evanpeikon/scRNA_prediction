# 🧬 Project Overview
## 🧬 Background
Physiological responses to exercise are generally predictable. For example, resistance training leads to increased strength and muscle hypertrophy, while endurance training improves aerobic capacity. However, modern exercise science still struggles to determine the optimal combination of volume, intensity, duration, frequency, and timing of exercise stimuli to maximize an individual’s health, fitness, and performance outcomes.

As George Brooks, a pioneer in the field, aptly put it:
> “It is wise to note that we are all individuals and that whereas physiological responses to particular stimuli are largely predictable, the precise responses and adaptations to those stimuli will vary among individuals. Therefore, the same training regimen may not equally benefit all those who follow it.”

Despite ongoing research, the field of exercise science has made limited progress since the early 2000s in refining these principles to account for individual variability. Achieving a deeper understanding of personalized, dose-response optimized exercise will require rigorous and systematic research that:
  1. Defines the optimal dose (i.e., volume and intensity), duration, frequency, and timing of exercise stimuli to produce specific adaptations.
  2. Intentionally investigates how combinations of potentially synergistic or antagonistic stimuli influence outcomes.
  3. Examines the effect of individual factors (e.g., age, training status, genetics, comorbidities, medication) on training adaptations.

Unfortunately, many exercise science studies are underpowered, lack proper controls, or suffer from limited funding, making it challenging to address these issues comprehensively. Another significant limitation is that most studies only measure broad end-responses to training (e.g., increases in strength or endurance) rather than examining the molecular changes driving these adaptations.

This is beginning to change with the application of single-cell RNA sequencing (scRNA-seq) techniques in skeletal muscle and exercise research. scRNA-seq allows for the analysis of gene expression at a single-cell level, providing a granular view of how different cell types within muscle tissue respond to exercise. The technique is valuable because it can capture the heterogeneity in cell populations, offering insights into specific gene expression patterns and molecular pathways involved in exercise adaptations.

However, scRNA-seq studies are still limited in scope (the specific problems they address) and breadth (the number of subjects and samples they can analyze). A key constraint is funding, as sequencing costs rise significantly with the number of subjects and the depth of sequencing required to capture more gene expression data. The more genes researchers aim to measure, the more library preparation and sequencing depth is needed—driving up costs exponentially.

## 🧬 Project Goals
The goal of this project is to explore whether we can predict the expression of a subset of genes in skeletal muscle cells using the expression levels of other genes as input to a predictive model. If successful, this approach could allow researchers to capture more samples before and after exercise with lower sequencing depth, reducing costs. By imputing or predicting gene expression levels, exercise scientists could gain deeper insights into the body’s adaptive responses to training.

While this method may not achieve the level of precision required for clinical or high-level bioinformatics research, it could provide a valuable tool for exercise scientists seeking to understand how different exercise stimuli affect gene expression. Additionally, there may be potential applications in high-performance sports, where athletes and coaches aim to optimize training regimens to maximize adaptation rates through variations in volume, intensity, and frequency.

# 🧬 Project Walkthrough

<img width="550" alt="Screenshot 2024-10-03 at 2 04 17 PM" src="https://github.com/user-attachments/assets/0b25e35f-0d93-4a9d-acc3-ca1d5a724552">

## 🧬 Data Availability 
The data for this project comes from a study called "Single-cell transcriptional profiles in human skeletal muscle," which is accessible via the Gene Expression Omnibus (GEO) via the following accession number: [GSE130646](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi). 

Additionally, we can access the Single-cell RNA sequencing data, resulting from each subject's muscle biopsy, using the following accession numbers: [GSM3746212](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3746212), [GSM3746213](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3746213), [GSM3746214](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3746214), and [GSM3746215](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3746215).

## 🧬 Import libraries 
Before loading the data above, I'll first import the following Python libraries, which will be used in downstream analyses:
```python
import pandas as pd
import scanpy as sc
import numpy as np
import anndata as an
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.pipeline import Pipeline
from sklearn.svm import SVR
from sklearn.multioutput import MultiOutputRegressor
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet
from sklearn.tree import DecisionTreeRegressor
from sklearn.model_selection import cross_val_score, KFold
import joblib
```

## 🧬 Load, Prepare, and Inspect Data
Next, we'll use Bash's wget command to retrieve each subject's Single-cell RNA-seq data, and following that, we'll use the Bash gunzip command to decompress the files:

```bash
# retrieve subect data from GEO
!wget -O GSM3746212_Muscle_1_Counts.csv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3746212&format=file&file=GSM3746212%5FMuscle%5F1%5FCounts%2Ecsv%2Egz'
!wget -O GSM3746213_Muscle_2_Counts.csv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3746213&format=file&file=GSM3746213%5FMuscle%5F2%5FCounts%2Ecsv%2Egz'
!wget -O GSM3746214_Muscle_3_Counts.csv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3746214&format=file&file=GSM3746214%5FMuscle%5F3%5FCounts%2Ecsv%2Egz'
!wget -O GSM3746215_Muscle_4_Counts.csv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3746215&format=file&file=GSM3746215%5FMuscle%5F4%5FCounts%2Ecsv%2Egz'

# decompress files
!gunzip GSM3746212_Muscle_1_Counts.csv.gz
!gunzip GSM3746213_Muscle_2_Counts.csv.gz
!gunzip GSM3746214_Muscle_3_Counts.csv.gz
!gunzip GSM3746215_Muscle_4_Counts.csv.gz
```

Now that we've decompressed the CSV files, we'll read them using the Pandas read_csv() function. Following that, we'll convert each CSV file into an individual AnnData object before combining them into a single composite object named adata_combined:

```python
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
```

We now have a composite AnnData object storing all of our subject's data. However, this composite object is formatted such that the rows are genes and the columns are cell IDs. As a result, we need to transpose our AnnData object so the rows are cell IDs and the columns are genes since this is the conventional format used for analyzing scRNA-seq data:

```python
# transpose AnnData object
adata_transposed = sc.AnnData(adata_combined.T)
```

Notably, I've saved the newly transposed AnnData object as adata_transposed, and as a result, all downstream analysis will use this newly created AnnData object instead of our original variable, adata_combined. Now, I'll print some basic summary information about our new AnnData object:

```python
num_genes = adata_transposed.n_vars # variables = columns (genes)
print(f"Number of Genes: {num_genes}")

num_cells = adata_transposed.n_obs # observations = rows (cells IDs)
print(f"Number of Cells: {num_cells}")
```
Which produces the following output:
  - Number of Genes: 15406
  - Number of Cells: 2876

## 🧬 Data Preperation
Before analyzing out single-cell data or building a model, we'll need to perform data preperation, which includes quality control, filtering, and normalization. First, we'll start with  quality control and inspect the integrity of our data and check if there are any missing values:

```python
# find indices of rows (cell IDs) with NaN values
nan_rows = np.isnan(adata_transposed.X).any(axis=1)
print(f"Number of rows with NaN values: {np.sum(nan_rows)}")

# find indices of columns (genes) with NaN values
nan_cols = np.isnan(adata_transposed.X).any(axis=0)
print(f"Number of columns with NaN values: {np.sum(nan_cols)}")
```
Which produces the following outout:
  - Number of rows with NaN values: 0
  - Number of columns with NaN values: 0

As you can see, there are no missing values in the data set. Next, we'll look at the distribution of genes per cel, cells per gene, and percent mitochondrial content per cell to determine our filtering criteria:

```python
# calculate the number of genes expressed per cell
adata_transposed.obs['n_genes'] = (adata_transposed.X > 0).sum(axis=1)

# calculate the number of cells in which each gene is expressed
adata_transposed.var['n_cells'] = (adata_transposed.X > 0).sum(axis=0)

# identify mitochondrial genes
mt_gene_mask = adata_transposed.var_names.str.startswith('MT-')

# calculate % mitochondrial genes per cell
if isinstance(adata_transposed.X, np.ndarray):
    adata_transposed.obs['percent_mito'] = np.sum(adata_transposed[:, mt_gene_mask].X, axis=1) / np.sum(adata_transposed.X, axis=1) * 100
else:
    adata_transposed.obs['percent_mito'] = np.sum(adata_transposed[:, mt_gene_mask].X.toarray(), axis=1) / np.sum(adata_transposed.X.toarray(), axis=1) * 100

# create subplots
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# plot histogram of the number of genes per cell
sns.histplot(adata_transposed.obs['n_genes'], bins=50, kde=True, ax=axes[0])
axes[0].set_xlabel('Number of Genes per Cell')
axes[0].set_ylabel('Number of Cells')
axes[0].set_title('Distribution of Number of Genes per Cell')

# plot histogram of the number of cells per gene
sns.histplot(adata_transposed.var['n_cells'], bins=50, kde=True, ax=axes[1])
axes[1].set_xlabel('Number of Cells per Gene')
axes[1].set_ylabel('Number of Genes')
axes[1].set_title('Distribution of Number of Cells per Gene')

# plot the distribution of mito gene %
sns.histplot(adata_transposed.obs['percent_mito'], bins=50, kde=True, ax=axes[2])
axes[2].set_xlabel('Percentage of Mitochondrial Genes')
axes[2].set_ylabel('Number of Cells')
axes[2].set_title('Distribution of Mitochondrial Gene % per Cell')

plt.tight_layout()
plt.show()
```

<img width="629" alt="Screenshot 2024-10-03 at 11 46 12 AM" src="https://github.com/user-attachments/assets/5fd50442-64be-47c2-872b-00ee9f314a71">

Based on the output above, it's clear that some filtering was already applied to this data before it was stored in the Gene Expression Omnibus (GEO). For example, there aren't any cells in the dataset with a mitochondrial gene percentage over ~5% (generally, we want to filter out cells with a mitochondrial gene percentage of >10% or >15%).

However, before moving on, I'll do some basic filtering, removing cells with fewer than 200 detected genes and genes that appear in fewer than 20 cells (I do not expect genes appearing in fewer than 20 cells to exist in this dataset based on the visualization above, but It's worth hedging our bets in case they were washed out of the visualization due to too few counts).

```python
# filter out cells that have fewer than 200 detected genes
sc.pp.filter_cells(adata_transposed, min_genes=200)

# filter out genes that appear in fewer than 20 cells
sc.pp.filter_genes(adata_transposed, min_cells=20) 

# print resulting number of cells and genes
num_genes = adata_transposed.n_vars
print(f"Number of Genes: {num_genes}")

num_cells = adata_transposed.n_obs
print(f"Number of cells: {num_cells}")
```
Which produces the following output:
  - Number of Genes: 10485
  - Number of Cells: 2876

As you can see in the output above, filtering our cells with fewer than 200 detected genes removed 4921 genes from the dataset, a 30% reduction in the number of genes present. However, our filtering criteria did not remove any cells from the dataset.

Following filtering, I'll apply global-scaling normalization, which consists of dividing each gene's expression level by the total expression in that cell, multiplying the result by a scaling factor to standardize the values, and then applying a log transformation to stabilize the variance. Then, following that, I'll identify the 2000 most highly variable genes in our dataset, which are useful for predicting unknown genes expression levels, and remove all other genes from the AnnData object:

```python
# apply global-scaling normalization
sc.pp.normalize_total(adata_transposed, target_sum=1e4)
sc.pp.log1p(adata_transposed)

# find the 2000 most highly variable genes 
sc.pp.highly_variable_genes(adata_transposed, n_top_genes=2000, subset=True)
print(adata_transposed)
```
Which produces the following output:
  - AnnData object with n_obs × n_vars = 2876 × 2000

As you can see in the output above, our AnnData object now contains 2187 cells and 2000 genes that exhibit the highest variability between cells.

## 🧬 Transformation and Dimensionality Reduction
As demonstrated above, single-cell RNA sequencing data is high-dimensional, consisting of thousands of genes measured across thousands of cells. High-dimensional data is challenging to analyze and interpret, so dimensionality reduction techniques are used to simplify the data by reducing the number of dimensions while simultaneously retaining the most important information in the dataset.

In the code block below, I'll tranform the dataset using Z-transforamtoon which standardizes the data so that each gene has a mean of zero and a standard deviation of one. This reduces the influence of genes with extremely high expression levels and ensures that all genes contribute equally to downstream analyses. Following that, I'll  use a dimensionality reduction technique called principal component analysis (PCA), which reduces the dimensionality of the data while capturing most of the variance.

```python
# apply z-transformation
sc.pp.scale(adata_transposed, zero_center=True)

# perform dimensionality reduction via PCA
sc.tl.pca(adata_transposed, svd_solver='arpack')
```

## 🧬 Building A Machine Learning Model To Predict Gene Expression
### Data Preperation 
First, we'll want to convert out AnnData object named ```adata_transposed``` into a DataFrame, which is the preffered data structure for many of the downstream analyses and processes we'll be performing. 

```python
# Convert expression matrix (adata_transposed.X) to a DataFrame using adata_transposed.obs_names as row labels (cell names) and adata_transposed.var_names as column labels (gene names)
df = pd.DataFrame(adata_transposed.X, index=adata_transposed.obs_names, columns=adata_transposed.var_names)
```

After the above process, we now have a DataFrame with 2876 rows, representing individual cells, and 2,000 columns with one for each gene in our dataset. Now, before evaluating different machine learning models on our dataset, we need to split our dataset in inputs and outputs. In this case, our inputs will be the first 1,950 genes in our dataset's expression levels/ 

```python
# first 1,950 genes as inputs to model
X = df.iloc[:, :1950].values
```

### Algorithm Spot-Checking
Now, any machine learning problem, we must select an algorithm to make predictions, an evaluation method to estimate a model's performance on unseen data, and an evaluation metric(s) to quantify how well the model works.

Unfortunately, we can't always know which algorithm will work best on our dataset beforehand. As a result, we have to try several algorithms, then focus our attention on those that seem most promising. Thus, it's important to have quick and easy ways to assess and compare different algorithms' performance before we select one to tune and optimize - this is where spot-checking comes in.

Spot-checking is a way to quickly discover which algorithms perform well on our machine-learning problem before selecting one to commit to. In the code block we'll we'll spot check six regression algoritms using the same [evaluation method](https://github.com/evanpeikon/machine-learning/tree/main/resampling) and [evaluation metric](https://github.com/evanpeikon/Machine-Learning/tree/main/Regression-metrics) to compare the model's performance. 

```python
# select model for spot-checking
models = {
    'LR': Pipeline([('scaler', StandardScaler()), ('model', LinearRegression())]),
    'Ridge': Pipeline([('scaler', StandardScaler()), ('model', Ridge())]),
    'Lasso': Pipeline([('scaler', StandardScaler()), ('model', Lasso())]),
    'ENR': Pipeline([('scaler', StandardScaler()), ('model', ElasticNet())]),
    'CART': Pipeline([('model', DecisionTreeRegressor())]),
    'SVM': Pipeline([('scaler', StandardScaler()), ('model', SVR())]),}

# list to store results for each model
results = []

# iterate over the last 50 output genes and calculate average MSE for each model
kfold = KFold(n_splits=5, random_state=5, shuffle=True)
for name, model in models.items():
    all_mse = []
    for gene_index in range(1950, 2000):
        y = df.iloc[:, gene_index].values
        cv_results = cross_val_score(model, X, y, cv=kfold, scoring='neg_mean_squared_error')
        all_mse.append(-cv_results.mean())
    avg_mse = np.mean(all_mse)
    std_mse = np.std(all_mse)
    results.append({'Model': name, 'Mean_MSE': avg_mse, 'Std_MSE': std_mse})

# convert results list to DF and display DF
results_df = pd.DataFrame(results)
print(results_df)

# visualize the average performance of each model across all genes
plt.figure(figsize=(10, 6))
plt.barh(results_df['Model'], results_df['Mean_MSE'], xerr=results_df['Std_MSE'], capsize=5)
plt.xlabel('Mean MSE')
plt.ylabel('Model')
plt.grid(axis='x')
plt.show()
```
Which produces the following output:

<img width="611" alt="Screenshot 2024-10-03 at 12 19 14 PM" src="https://github.com/user-attachments/assets/7c89e084-e2dd-4001-ab19-009358c4918c">

As you can see in the image above, support vector regresson (SVR) outperformed all other algorithms with the lowest average mean squarred error (MSE) across the 50 output gene's expression levels we are aiming to predict. A potential reason for this is that SVR is well-suited for handling handling high-dimensional, and noisy, data which are common characteristics of gene expression datasets. Additioanlly, SVR may capture complex relationships between genes that linear models like linear regression (which was the worst performing mode) may miss. SVR's regularization techniques also help prevent overfitting, ensuring that the model generalizes well to new data, which is critical for the use case described in the project overview section above. 

### Hyperparameter Tuning 
You can think of machine learning algorithms as systems with various knobs and dials, which you can adjust in any number of ways to change how output data (predictions) are generated from input data. The knobs and dials in these systems can be subdivided into parameters and hyperparameters.

Parameters are model settings that are learned, adjusted, and optimized automatically. Conversely, hyperparameters need to be manually set manually by whoever is programming the machine learning algorithm. Generally, tuning hyperparameters has known effects on machine learning algorithms. However, it’s not always clear how to best set a hyperparameter to optimize model performance for a specific dataset. As a result, search strategies are often used to find optimal hyperparameter configurations. In the code block below, we'll use [random search](https://github.com/evanpeikon/Machine-Learning/tree/main/hyperparameter_optimization), which is a tuning technique that randomly samples a specified number of uniformly distributed algorithm parameters.

> Note: the code below uses a subset of the dataset to tune hyperparameters to reduce computational cost and runtime. 

```python
# define the SVR pipeline
svr_pipeline = Pipeline([('scaler', StandardScaler()), ('model', SVR())])

# define the hyperparameter space for random search w/ reduced parameter range
param_distributions = {'model__C': [0.1, 1, 10], 'model__epsilon': [0.01, 0.1], 'model__kernel': ['linear', 'rbf'], 'model__gamma': ['scale', 'auto'],}

random_search = RandomizedSearchCV(estimator=svr_pipeline, param_distributions=param_distributions, n_iter=10, scoring='neg_mean_squared_error', cv=3, random_state=5, verbose=1, n_jobs=-1)

# input values
X = df.iloc[:, :1950].values

# randomly sample 20% of the rows
X_sub, _, y_sub, _ = train_test_split(X, df.iloc[:, 1950].values, test_size=0.8, random_state=5)

# perform random search on the smaller subset
random_search.fit(X_sub, y_sub)

# print the best parameters
print(f"Best parameters on subset: {random_search.best_params_}")
print(f"Best MSE on subset: {-random_search.best_score_}")
```
Which produces the following outputs:
- Best parameters on subset: {'model__kernel': 'rbf', 'model__gamma': 'auto', 'model__epsilon': 0.01, 'model__C': 10}

Now that we've have the optimal hyperparameters to optimize our models peformance, we can create pipeline. 

### Creating A Pipeline 

```python
# separate the target genes (last 50 genes)
target_genes_indices = list(range(1950, 2000))
y = df.iloc[:, target_genes_indices].values
X = df.drop(columns=df.columns[target_genes_indices]).values

# split into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# ceate a pipeline with preprocessing and model using optimal hyperparameters
pipeline = Pipeline(steps=[
    ('scaler', StandardScaler()),
    ('model', MultiOutputRegressor(SVR(kernel='rbf', C=10, gamma='auto', epsilon=0.01)))])

# fit the pipeline on the training data
pipeline.fit(X_train, y_train)

# save the pipeline to a file for later use
joblib.dump(pipeline, 'gene_expression_pipeline.pkl')
```

To use the saved pipeline for predicting gene expression counts on a new dataset, you'll need to follow a few straightforward steps. Here’s how you can do it:

Step-by-Step Instructions
Load the Saved Pipeline: Use joblib to load the pipeline you previously saved.
Prepare the New Dataset: Ensure that your new dataset is preprocessed in the same way as your training data (e.g., same feature columns, handling missing values if applicable).
Make Predictions: Use the loaded pipeline to make predictions on the new data.


# 🧬 Conclusion
This project demonstrates the potential of our predictive model for predicting gene expression and achieving the broader objectives discussed in the project overview above. However, to enhance its utility in practical applications, it is crucial to gather more training data and develop a more modular framework that can be tailored to specific study or field-based requirements. Additionally, there are several improvements I envision for future iterations that could significantly bolster the model’s robustness and broaden its applicability. Unfortunately, due to computational and resource constraints, these enhancements were not feasible within the scope of this project. One notable direction for future work would be to cluster skeletal muscle cells based on known marker genes, then to create individual models to predict gene expression for each cell type. This approach could potentially enhance the model's accuracy and overall utility.

I welcome any suggestions for improvements, alternative approaches, or ideas for extending this project. Please feel free to reach out to me at evanpeikon@gmail.com. 
