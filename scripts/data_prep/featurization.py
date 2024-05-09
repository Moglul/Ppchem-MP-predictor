import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, MACCSkeys, RDKFingerprint, Descriptors
import pandas as pd
from rdkit.ML.Descriptors import MoleculeDescriptors

def generate_mordred_columns(df, smiles_column = "SMILES", ignore_3D = True):
    """
    Generates Mordred features for a DataFrame with a 'smiles' or 'SMILES' column.

    Parameters:
    df (pd.DataFrame): The DataFrame with a 'smiles' or 'SMILES' column.
    smiles_column (str): The name of the column containing SMILES strings.
    ignore_3D (bool): Whether to ignore 3D descriptors.

    Returns:
    pd.DataFrame: The DataFrame augmented with Mordred features.
    """
    def calculate_all_mordred_descriptors(smiles):
        mol = Chem.MolFromSmiles(smiles)
        calc = Calculator(descriptors, ignore_3D=ignore_3D)
        mordred_descriptors = calc(mol)
        return mordred_descriptors.fill_missing().asdict()

    # List to store dictionaries of descriptors for each molecule
    descriptor_dicts = []

    # Iterate over each SMILES string and calculate descriptors
    for smiles in df[smiles_column]:
        descriptor_dict = calculate_all_mordred_descriptors(smiles)
        descriptor_dicts.append(descriptor_dict)

    # Convert the list of dictionaries into a DataFrame
    df_descriptors = pd.DataFrame(descriptor_dicts)

    # Reset the indices of df and df_descriptors
    df.reset_index(drop=True, inplace=True)
    df_descriptors.reset_index(drop=True, inplace=True)

    # Concatenate the original DataFrame with the descriptor DataFrame
    df_combined = pd.concat([df.reset_index(drop=True), df_descriptors.reset_index(drop=True)], axis=1)

    return df_combined

def generate_fingerprint_descriptors(df, smiles_column='SMILES', fingerprint_type='Morgan', n_bits=2048):
    """
    Generates fingerprint descriptors for SMILES strings in a DataFrame.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the SMILES strings.
    - smiles_column (str): Column name for SMILES strings.
    - fingerprint_type (str): Type of fingerprint to generate. Supported types: 'Morgan', 'MACCS', 'Topological'.
    - n_bits (int): Number of bits in the fingerprint. Default is 2048 for Morgan fingerprints.

    Returns:
    - pd.DataFrame: DataFrame with fingerprint descriptors as lists of bits.
    """
    def compute_fingerprint(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        if fingerprint_type == "Morgan":
            return list(rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=n_bits))
        elif fingerprint_type == "MACCS":
            return list(MACCSkeys.GenMACCSKeys(mol))
        elif fingerprint_type == "Topological":
            return list(RDKFingerprint(mol))
        else:
            raise ValueError("Unsupported fingerprint type")

    # Apply the compute_fingerprint function to the SMILES column
    fingerprints = df[smiles_column].apply(compute_fingerprint)

    # Create a new DataFrame for the fingerprints
    fingerprint_df = pd.DataFrame(fingerprints.tolist(), columns=[f'{fingerprint_type}_bit_{i}' for i in range(len(fingerprints.iloc[0]))])

    # Set the index of the fingerprints DataFrame to match the original DataFrame
    fingerprint_df.index = df.index

    # Add fingerprint columns to the original DataFrame
    df = pd.concat([df, fingerprint_df], axis=1)

    return df

def RDkit_descriptors(df, smiles_column='SMILES'):
    """
    Calculates RDKit descriptors for a DataFrame containing SMILES.

    Parameters:
    - df (pd.DataFrame): DataFrame containing SMILES.
    - smiles_column (str): The column containing SMILES strings.

    Returns:
    - pd.DataFrame: DataFrame with the calculated descriptors.
    """
    # Get the SMILES strings
    smiles = df[smiles_column].tolist()
    # Get the mol objects
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    # Get the descriptor names
    descriptor_names = [x[0] for x in Descriptors._descList]
    # Remove Ipc descriptor from the list because it can cause problems with the newest numpy version
    descriptor_names.remove('Ipc')
    # Calculate the descriptors
    calculator = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)
    descriptors = [calculator.CalcDescriptors(mol) for mol in mols]
    # Create a DataFrame
    descriptors_df = pd.DataFrame(descriptors, columns=descriptor_names)
    # Add the descriptor to the original DataFrame
    df = df.reset_index(drop=True)
    df = pd.concat([df, descriptors_df], axis=1)
    return df