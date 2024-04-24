# MoleculAR

MoleculAR is a Python package designed for cheminformatics applications, particularly for machine learning tasks involving molecular data in a form of a dataframe with a "smiles" and "target" column to create a prediction of model. It provides functionalities for datacleaning, descriptor generation, feature preparation, regression modeling, and SMILES cleaning.

## Installation

You can install MoleculAR using pip:

```
pip install moleculAR
```

## Usage

### 1. SMILES Cleaning

The `smiles_cleaning` module provides functions for cleaning SMILES strings, including canonicalization and removal of invalid SMILES and unwanted substructures from a "smiles" or "SMILES" column in a dataframe.

```python
from moleculAR.preprocessing import smiles_cleaning

# Example usage:

# canonical_smiles = smiles_cleaning.canonical_smiles(smiles)

# cleaned_df = smiles_cleaning.remove_invalid_smiles(df, smiles_column='SMILES')

# substructure_removed_df = smiles_cleaning.remove_unwanted_substructures(df, smiles_column='SMILES', substructure_smiles="C(=O)N")
```

### 2. Descriptor Generation

The `descriptor_generation` module provides functions for generating molecular descriptors from SMILES strings.

```python
from moleculAR.preprocessing import descriptor_generation

# Example usage:

# descriptors_df = descriptor_generation.generate_extended_descriptors(df, smiles_column='SMILES')

# fingerprint_df = descriptor_generation.generate_fingerprint_descriptors(df, smiles_column='SMILES', fingerprint_type='Morgan', n_bits=2048)
```

### 3. Feature Preparation

The `feature_preparation` module offers functions for preparing features for machine learning tasks, including descriptor selection and handling missing data.

```python
from moleculAR.preprocessing import feature_preparation

# Example usage:

# selected_features_df = feature_preparation.descriptor_selection(df, exclude_columns=['Target', "smiles], min_variance=0.05, max_correlation=0.9)

# cleaned_df = feature_preparation.handle_missing_data(df, threshold=0.8, fill_method='mean')
```

### 4. Regression Models

The `regression_models` module contains functions for training and evaluating regression models.

```python
from moleculAR.models import regression_models

# Example usage:

# regression_scores_df = regression_models.compare_regression_models(X, y)

# regression_plot = regression_models.plot_predictions(model, X, y)
```