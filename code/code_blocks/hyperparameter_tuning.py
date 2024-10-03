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
