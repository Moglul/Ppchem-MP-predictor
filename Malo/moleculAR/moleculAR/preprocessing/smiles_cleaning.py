from rdkit import Chem
from rdkit import Chem
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem

def canonical_smiles(smiles):
    """
    Converts a given SMILES string to its canonical form using RDKit.

    Parameters:
    - smiles (str): A SMILES string representing a chemical molecule.

    Returns:
    - str: The canonical SMILES string if the input is valid; otherwise, None.

    Raises:
    - ValueError: If the input is not a valid SMILES string.
    """
  
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Chem.MolToSmiles(mol, canonical=True)
    else:
        return None

def clean_smiles_dataframe(df, column = "SMILES"):
    """
    Remove rows with invalid and duplicate canonical SMILES from the DataFrame.

    Parameters:
    - df (pd.DataFrame): DataFrame containing a smiles column.
    - column (str) : The name of the column containing all smiles value. Default : "SMILES"

    Returns:
    - pd.DataFrame: DataFrame with invalid and duplicate canonical SMILES removed.
    """
    # Save the initail number of rows
    initial_rows = df.shape[0]

    # Check if "smiles" column exists
    if column not in df.columns:
        raise ValueError(f"DataFrame does not contain a {column} column.")

    # Remove rows with invalid SMILES
    valid_smiles_mask = df[column].apply(lambda x: Chem.MolFromSmiles(x) is not None)
    df_valid_smiles = df.loc[valid_smiles_mask]

    # Remove rows with duplicate canonical SMILES
    df_valid_smiles.loc[:, 'canonical_smiles'] = df_valid_smiles[column].apply(canonical_smiles)
    df_unique = df_valid_smiles.drop_duplicates(subset='canonical_smiles')

    # Print the number of rows removed
    removed_rows = initial_rows - df_unique.shape[0]
    print(f"Removed {removed_rows} rows with invalid or duplicate SMILES.")

    return df_unique.drop(columns=['canonical_smiles'])


def remove_unwanted_substructure(df, substructure_smiles, smiles_column='SMILES'):
    """
    Removes molecules containing specific substructures.

    Parameters:
    - df (pd.DataFrame): DataFrame containing SMILES.
    - substructure_smiles (str): SMILES of the substructure to be removed.
    - smiles_column (str): The column containing SMILES strings.
    

    Returns:
    - pd.DataFrame: DataFrame without the molecules containing the substructure.
    """
    substructure = Chem.MolFromSmiles(substructure_smiles)
    def contains_substructure(smiles):
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None and not mol.HasSubstructMatch(substructure)

    return df[df[smiles_column].apply(contains_substructure)]

def keep_only_substructure(df, substructure_smiles, smiles_column='SMILES'):
    """
    Keeps molecules containing a specific substructure and removes others.

    Parameters:
    - df (pd.DataFrame): DataFrame containing SMILES.
    - substructure_smiles (str): SMILES of the substructure to keep.
    - smiles_column (str): The column containing SMILES strings.

    Returns:
    - pd.DataFrame: DataFrame with only the molecules containing the specified substructure.
    """
    substructure = Chem.MolFromSmiles(substructure_smiles)
    
    def contains_substructure(smiles):
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None and mol.HasSubstructMatch(substructure)

    return df[df[smiles_column].apply(contains_substructure)]

def remove_unwanted_atoms(df, unwanted_atoms, smiles_column='SMILES'):
    """
    Removes molecules containing specific atoms.

    Parameters:
    - df (pd.DataFrame): DataFrame containing SMILES.
    - unwanted_atoms (list): List of atomic symbols to be removed.
    - smiles_column (str): The column containing SMILES strings.

    Returns:
    - pd.DataFrame: DataFrame without the molecules containing the unwanted atoms.
    """
    def contains_unwanted_atoms(smiles):
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None and not any(atom.GetSymbol() in unwanted_atoms for atom in mol.GetAtoms())

    return df[df[smiles_column].apply(contains_unwanted_atoms)]

def keep_only_atoms(df, wanted_atoms, smiles_column='SMILES'):
    """
    Keeps molecules containing specific atoms.

    Parameters:
    - df (pd.DataFrame): DataFrame containing SMILES.
    - wanted_atoms (list): List of atomic symbols to be kept.
    - smiles_column (str): The column containing SMILES strings.

    Returns:
    - pd.DataFrame: DataFrame with only the molecules containing the wanted atoms.
    """
    def contains_wanted_atoms(smiles):
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None and all(atom.GetSymbol() in wanted_atoms for atom in mol.GetAtoms())

    return df[df[smiles_column].apply(contains_wanted_atoms)]
