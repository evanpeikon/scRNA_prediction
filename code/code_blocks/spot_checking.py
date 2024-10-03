# Convert expression matrix (adata_transposed.X) to a DataFrame using adata_transposed.obs_names as row labels (cell names) and adata_transposed.var_names as column labels (gene names)
df = pd.DataFrame(adata_transposed.X, index=adata_transposed.obs_names, columns=adata_transposed.var_names)

# first 1,950 genes as inputs to model
X = df.iloc[:, :1950].values

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
