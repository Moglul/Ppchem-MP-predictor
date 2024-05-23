import pytest
import pandas as pd
import numpy as np
from sklearn.datasets import make_regression
import sys
from regression_models import compare_regression_models

# Create a sample regression dataset
X, y = make_regression(n_samples=100, n_features=3, noise=0.3, random_state=142)
X = pd.DataFrame(X, columns=['Feature 1', 'Feature 2', 'Feature 3'])
y = pd.Series(y, name='Target')

def test_compare_regression_models():
    # Test with all models
    result = compare_regression_models(X, y)
    assert isinstance(result, pd.DataFrame)
    assert 'R2 Score' in result.columns
    assert 'RMSE' in result.columns

    # Test with specific models
    result = compare_regression_models(X, y, models=['Linear Regression', 'XGBoost'])
    assert result.shape[0] == 2

    # Test with invalid model name
    with pytest.raises(KeyError):
        compare_regression_models(X, y, models=['Invalid Model'])
