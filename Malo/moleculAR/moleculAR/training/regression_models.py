import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.linear_model import LinearRegression
from xgboost import XGBRegressor
from lightgbm import LGBMRegressor
import matplotlib.pyplot as plt

def compare_regression_models(X, y, test_size=0.2, random_state=12, models=None):
    """
    Trains various regression models and evaluates them on the provided dataset.
    
    Parameters:
    - X (pd.DataFrame): Feature data.
    - y (pd.Series): Target data.
    - test_size (float or int): The proportion of the dataset to include in the test split.
    - random_state (int): Controls the shuffling applied to the data before applying the split.
    - models (list or None): List of model names to try. If None, all models will be tried.
    
    Returns:
    - pd.DataFrame: DataFrame containing the R2 and RMSE scores for each model.
    """
    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=random_state)
    
    # Initialize models
    all_models = {
        'Linear Regression': LinearRegression(),
        'Random Forest': RandomForestRegressor(n_estimators=100),
        'Gradient Boosting': GradientBoostingRegressor(n_estimators=100),
        'XGBoost': XGBRegressor(objective='reg:squarederror'),
        'LightGBM': LGBMRegressor()
    }
    
    if models is None:
        models = list(all_models.keys())  # If models are not specified, try all models
    
    # Filter models based on user input
    models_to_try = {name: model for name, model in all_models.items() if name in models}
    
    # Dictionary to store scores
    scores = {}

    # Train and evaluate each model
    for name, model in models_to_try.items():
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        rmse = float(np.sqrt(mean_squared_error(y_test, y_pred)))
        r2 = float(r2_score(y_test, y_pred))
        scores[name] = {'R2 Score': r2, 'RMSE': rmse}

    # Convert scores dictionary to DataFrame for nicer formatting
    scores_df = pd.DataFrame(scores).transpose()
    scores_df.sort_values(by='R2 Score', ascending=False, inplace=True)
    
    return scores_df

def plot_predictions(model, X, y, test_size=0.2, random_state = 42):
    """
    Trains the given model and plots predicted vs. actual values along with displaying R² and RMSE.

    Parameters:
    - model: A regression model instance that follows the scikit-learn model interface.
    - X (pd.DataFrame or np.array): Feature data.
    - y (pd.Series or np.array): Actual target values.
    - Test size (float or int) : The test size for the data train-test spliting
    - Random state (int) : The random state for the data train-test spliting

    Outputs:
    - R2 score and RMSE with a plot showing actual vs. predicted values.
    """
    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=random_state)
    
    # Train the model
    model.fit(X_train, y_train)
    
    # Make predictions
    y_pred = model.predict(X_test)
    
    # Calculate R2 and RMSE
    r2 = r2_score(y_test, y_pred)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))
    
    # Plotting
    plt.figure(figsize=(10, 6))
    plt.scatter(y_test, y_pred, alpha=0.6, edgecolors="k", color="blue", label="Predictions")  
    plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--', lw=2, label="Perfect Predictions") 
    plt.xlabel("Actual Values", fontsize=12, color="black")  
    plt.ylabel("Predicted Values", fontsize=12, color="black")  
    plt.title("Actual vs. Predicted Values", fontsize=14, color="black")  
    plt.grid(True)

    # Print R2 and RMSE
    print(f'R² Score: {r2:.4f}')
    print(f'RMSE: {rmse:.4f}')

    # Set legend and legend fontsize
    plt.legend(fontsize=10)

    plt.show()


