import pytest
import pandas as pd
import numpy as np
from sklearn.datasets import make_regression
from scripts.regression_models import compare_regression_models, plot_predictions
from sklearn.linear_model import LinearRegression

# Create a sample regression dataset
X, y = make_regression(n_samples=100, n_features=1, noise=0.1, random_state=42)
X = pd.DataFrame(X, columns=['Feature'])
y = pd.Series(y, name='Target')

def test_compare_regression_models():
    # Test with all models
    result = compare_regression_models(X, y)
    assert isinstance(result, pd.DataFrame)
    assert 'R2 Score' in result.columns
    assert 'RMSE' in result.columns

    # Test with specific models
    result = compare_regression_models(X, y, models=['Linear Regression', 'Random Forest'])
    assert result.shape[0] == 2

    # Test with invalid model name
    with pytest.raises(KeyError):
        compare_regression_models(X, y, models=['Invalid Model'])

def test_plot_predictions(capsys):
    # Test with Linear Regression model
    plot_predictions(LinearRegression(), X, y)
    captured = capsys.readouterr()
    assert 'R² Score:' in captured.out
    assert 'RMSE:' in captured.out

    # Test with invalid model
    with pytest.raises(AttributeError):
        plot_predictions('Invalid Model', X, y)