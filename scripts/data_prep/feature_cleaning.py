import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, MACCSkeys, RDKFingerprint
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler

def remove_low_variance_features(df, excluded_cols = [], threshold=0.05):
    """
    Removes low variance columns from a DataFrame, ignoring excluded_cols columns.

    Parameters:
    df (pd.DataFrame): The input DataFrame.
    excluded_cols (list): List of column names to exclude from variance analysis.
    threshold (float): The variance threshold.

    Returns:
    pd.DataFrame: The DataFrame with low variance features removed.
    """
    # Separate excluded columns from the DataFrame
    df_excluded = df[excluded_cols]
    df_selected = df.drop(columns=excluded_cols)

    # Remove low variance features
    selector = VarianceThreshold(threshold=threshold)
    df_reduced = selector.fit_transform(df_selected)
    selected_columns = df_selected.columns[selector.get_support()]

    # Convert the transformed array back to a DataFrame
    df_reduced = pd.DataFrame(df_reduced, columns=selected_columns)

    # Combine the reduced DataFrame with the excluded columns
    final_df = pd.concat([df_excluded, df_reduced], axis=1)

    return final_df

def remove_highly_correlated_features(df, excluded_cols = None, threshold=0.95):
    """
    Removes highly correlated columns from a DataFrame, keeping only one column from each group of highly correlated columns.
    Excludes string columns and certain specified columns from the correlation analysis.

    Parameters:
    df (pd.DataFrame): The input DataFrame.
    excluded_cols (list): List of column names to exclude from correlation analysis.
    threshold (float): The correlation coefficient threshold to consider for removing columns.

    Returns:
    pd.DataFrame: The DataFrame with highly correlated columns removed.
    """
    # Exclude string columns and any additional specified columns from the correlation analysis
    numerical_cols = df.select_dtypes(include=[np.number]).columns.difference(excluded_cols)
    
    # Calculate the correlation matrix for the filtered DataFrame
    df_filtered = df[numerical_cols]
    corr_matrix = df_filtered.corr().abs()

    # Select upper triangle of correlation matrix
    upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))

    # Find columns with correlation greater than the threshold
    to_drop = [column for column in upper.columns if any(upper[column] > threshold)]

    # Drop highly correlated columns but keep the excluded columns intact
    df_reduced = df.drop(columns=to_drop)

    return df_reduced

def descriptor_selection(df, exclude_cols = None, min_variance=0.05, max_correlation=0.9):
    """
    Selects descriptors based on variance and inter-correlation, excluding specified columns.

    Parameters:
    - df (pd.DataFrame): DataFrame containing descriptors and any other columns.
    - exclude_columns (list of str): List of column names to exclude from feature selection.
    - min_variance (float): Minimum variance threshold for feature selection.
    - max_correlation (float): Maximum allowable correlation between features to keep them.

    Returns:
    - pd.DataFrame: DataFrame with selected descriptors and the excluded columns.

    Note : All features will be Standardize in the output.
    """
    # Ensure all excluded columns are in the DataFrame
    if not all(col in df.columns for col in exclude_cols):
        raise ValueError("Some excluded columns are not in the DataFrame.")

    # Filter out the excluded columns when considering feature selection
    features = df.drop(columns=exclude_cols)

    # Apply variance threshold to reduce the feature set
    selector = VarianceThreshold(threshold=min_variance)
    features_reduced = selector.fit_transform(features)
    features_reduced = pd.DataFrame(features_reduced, columns=features.columns[selector.get_support(indices=True)])

    # Normalize features before calculating correlation to avoid scale impacts
    scaler = StandardScaler()
    features_scaled = pd.DataFrame(scaler.fit_transform(features_reduced), columns=features_reduced.columns)

    # Calculate correlation matrix and remove highly correlated features
    corr_matrix = features_scaled.corr().abs()
    upper_triangle = corr_matrix.where(pd.np.triu(pd.np.ones(corr_matrix.shape), k=1).astype(bool))
    to_drop = [column for column in upper_triangle.columns if any(upper_triangle[column] > max_correlation)]
    features_selected = features_reduced.drop(columns=to_drop)

    # Combine the selected features with the excluded columns
    final_df = pd.concat([df[exclude_cols], features_selected], axis=1)

    return final_df

def handle_missing_data(df, excluded_cols = [], threshold=0.8, fill_method="mean"):
    """
    Handles missing data in a DataFrame based on the percentage of non-missing data per column.
    Columns with non-missing data above a specified threshold are filled with either the mean or median,
    while columns below the threshold are removed from the DataFrame.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the data with potential missing values.
    - excluded_cols (list): List of column names to exclude from missing data handling.
    - threshold (float): Proportion of non-missing values required to keep and fill the column. 
                         Values range from 0 to 1, where 1 means no missing values are allowed for the column to be retained.
    - fill_method (str): Method for filling missing values. Options: 'mean', 'median'

    Returns:
    - pd.DataFrame: DataFrame with handled missing data, where columns with too many missing values have been removed and others filled.
    """
    df_selected = df.copy()
    df_excluded = pd.DataFrame()
    # Separate excluded columns from the DataFrame if needed
    if len(excluded_cols) != 0 :
        df_excluded = df[excluded_cols]
        df_selected = df.drop(columns=excluded_cols)

    # Iterate over each column and decide to fill or remove
    for column in df_selected.columns:
        # Calculate the proportion of non-missing values
        non_missing_ratio = df_selected[column].notna().mean()
        
        if non_missing_ratio >= threshold:
            # If the non-missing ratio is above the threshold, fill missing values
            if fill_method == "mean":
                df_selected[column].fillna(df_selected[column].mean(), inplace=True)
            elif fill_method == "median":
                df_selected[column].fillna(df_selected[column].median(), inplace=True)
            else:
                raise ValueError("Invalid fill method. Choose 'mean' or 'median'.")
        else:
            # If the non-missing ratio is below the threshold, drop the column
            df_selected.drop(column, axis=1, inplace=True)
    
    return pd.concat([df_excluded, df_selected], axis=1)