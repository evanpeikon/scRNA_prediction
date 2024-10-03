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
