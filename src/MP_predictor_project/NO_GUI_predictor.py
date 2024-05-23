"""
This script is used to start the prediction of the melting point of a molecule without using any GUI.
"""
# Important for making paths generalizable
from pathlib import Path

# Needed to import the predictor function
import sys

# Used to check the SMILES
import rdkit
from rdkit import Chem

# Import of the function
from predictor import prediction

# Variable to control the status of the program
status = True

# Loop to make multiple predictions, or getting out of the program
i = True

# Print the welcome message
print("Welcome to the Melting Point Predictor !\nIf you want to exit the programm, type 'exit'.")

# Loop to make predictions
while status:
    # Input of the SMILES, and check if it is valid
    smiles = input("Enter the SMILES of the molecule: ")
    if smiles == 'exit':
        status = False
        break
    elif Chem.MolFromSmiles(smiles) is not None:
        print(f"The melting of this molecule is: {round(prediction(smiles),2)}°C")
    else:
        print("The SMILES is not valid.")
        i = False
    # Loop to check if the user wants to make another prediction
    while i:
        choice = input("Do you want to make another prediction? (y/n): ")
        if choice.lower() == 'y':
            i = False
        elif choice.lower() == 'n':
            i = False
            status = False
        # Case if the input is not valid
        else:
            print("Invalid input. Please enter 'y' or 'n'.")
    # Reset the loop variable
    i = True
    
