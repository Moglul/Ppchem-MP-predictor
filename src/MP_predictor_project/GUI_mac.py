# Libraries used for creating the GUI
import tkinter as tk
import customtkinter
from tkinter import *

# Important for making paths generalizable
from pathlib import Path

# Needed to import the predictor function
import sys

# Allows you to import images as well as logos
from PIL import ImageTk, Image

# This allows you to create binds
import keyboard

# Allows you to check that the SMILES is acceptable, but also to display the molecule
from rdkit import Chem
from rdkit.Chem import Draw

# Useful for making a link clickable on the interface (github)
import webbrowser

# Allows you to perform two processes simultaneously, such as interface and prediction
import threading

# Get the directory of the current script
current_dir = Path(__file__).resolve().parent
# Construct the relative path to the scripts directory
scripts_dir = current_dir.parent.parent / 'scripts'
# Add the scripts directory to the Python path
sys.path.append(str(scripts_dir))
# Import of the function
from predictor import prediction


# start_gui function to create the GUI

def start_gui():
    
    ######## General settings ########

    # Path for the different assets to load
    script_dir = Path(__file__).resolve().parent
    assets_dir = script_dir.parent.parent / 'assets'


    # Creation of the window, using customtkinter
    window = customtkinter.CTk()
    customtkinter.set_appearance_mode("#1A1A1A")
    window.title("Melting point predictor")
    window.geometry("800x600")
    window.resizable(True, True)


    # Definition of a new icon
    path_icon = assets_dir / 'icon.ico'
    window.iconbitmap(path_icon)


    ########## Functions ##########


    # Main function, takes the smiles, controls and displays the melting point
    def getEntry_and_Test(entry):

        # Makes the submit button disappear
        submit_button.place_forget()
        submit_button.config(state=tk.DISABLED)

        # Makes the previous prediction, or error message disappear
        error_label.config(text="", fg="red", font=("Gill Sans MT", 13 * -1))
        answer_label.config(text="")

        # Makes appear the progressbar
        make_appear_progressbar()

        # User input recovery
        smiles = entry.get()

        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                # Molecule exists, draw and predicts...
                print_molecule(mol)
                print_prediciton(smiles)
            else:
                # Molecule doesn't exist, prints error and no molecule...
                error_label.config(text="Invalid SMILES", fg="red", font=("Gill Sans MT", 13 * -1))
                print_empty_molecule()
        finally:
            # Reappears the submit button
            submit_button.config(state=tk.NORMAL)
            submit_button.place(x=275, y=217)


    # Function used to display the molecule on the GUI, taking mol as argument of getEntry_and_Test function
    def print_molecule(mol):

        # Draws the molecule using the argument, and saves it as an image
        img = Draw.MolToImage(mol, size=(200, 140))
        img = ImageTk.PhotoImage(img)

        # Takes the image and displays it
        canvas.create_image(120, 330, anchor=NW, image=img)
        canvas.image = img


    # Function used to erase the space dedicated to drawing the molecule
    def print_empty_molecule():

        # Takes an empty SMILES and converts it to a molecule (null)
        smiles = ""
        mol = Chem.MolFromSmiles(smiles)

        # Draws the molecule (null) using mol defined before, and saves it as an image
        img = Draw.MolToImage(mol, size=(200, 140))
        img = ImageTk.PhotoImage(img)

        # Takes the image and displays it
        canvas.create_image(120, 330, anchor=NW, image=img)
        canvas.image = img

        # Makes disappear the progressbar
        make_disappear_progressbar()


    # Function that launches the prediction and displays the result
    def print_prediciton(smiles):

        # Launching of the prediction, rounds the result to 2 decimals
        answer = prediction(smiles)
        answer = round(answer, 2)

        # Displays the result
        answer_label.config(text=str(answer) + "°C", fg="white", font=("Gill Sans MT", 13 * -1))

        # Makes disappear the progressbar
        make_disappear_progressbar()


    # Function which defines and displays the progressbar
    def make_appear_progressbar():

        # Settings of the progressbar
        global progress_bar
        progress_bar = customtkinter.CTkProgressBar(window, 
            orientation="HORIZONTAL", 
            mode="indeterminate", 
            height=3, 
            width=640, 
            fg_color="#1A1A1A", 
            bg_color="#1A1A1A", 
            border_color="#1A1A1A", 
            progress_color="#305EA2", 
            indeterminate_speed=(1)
        )

        # Defines the position and displays it
        progress_bar.pack(pady=73)
        progress_bar.start()


    # Function that stops and deletes the progressbar
    def make_disappear_progressbar():
        progress_bar.stop()
        progress_bar.pack_forget()


    # Defines a callback function for the github link
    def callback(url):
        webbrowser.open_new_tab(url)
    

    # Simulate button press animation by changing the relief style temporarily
    def button_clicked_animation(button_name):
        button_name.config(relief=tk.SUNKEN, activebackground="#2A2A2A", borderwidth=0, highlightthickness=0)
        window.after(70, lambda: button_name.config(relief="flat", activebackground="#2A2A2A", borderwidth=0, highlightthickness=0))


    ########## Objects ##########


    canvas = Canvas(
    window,
    bg = "#1A1A1A",
    height = 600,
    width = 800,
    bd = 0,
    highlightthickness = 0,
    relief = "ridge"
)
    canvas.place(x = 0, y = 0)

    image_left_container = PhotoImage(
        file=assets_dir / 'left_container.png')
    left_container = canvas.create_image(
        220.0,
        337.0,
        image=image_left_container
    )

    image_upper_container = PhotoImage(
        file=assets_dir / 'upper_container.png')
    upper_container = canvas.create_image(
        400.0,
        46.0,
        image=image_upper_container
    )

    image_moldraw_box = PhotoImage(
        file=assets_dir / 'moldraw_box.png')
    moldraw_box = canvas.create_image(
        219.0,
        401.0,
        image=image_moldraw_box
    )

    image_epfl_logo = PhotoImage(
        file=assets_dir / 'epfl_logo.png')
    epfl_logo = canvas.create_image(
        735.0,
        48.0,
        image=image_epfl_logo
    )

    image_up_right_container = PhotoImage(
        file=assets_dir / 'up_right_container.png')
    up_right_container = canvas.create_image(
        603.0,
        249.0,
        image=image_up_right_container
    )

    image_down_right_container = PhotoImage(
        file=assets_dir / 'down_right_container.png')
    down_right_container = canvas.create_image(
        603.0,
        423.0,
        image=image_down_right_container
    )

    image_prediction_textbox = PhotoImage(
        file=assets_dir / 'prediction_textbox.png')
    prediction_textbox = canvas.create_image(
        603.0,
        248.0,
        image=image_prediction_textbox
    )

    canvas.create_text(
        18.0,
        24.0,
        anchor="nw",
        text="Melting point predictor\n\n\n",
        fill="#FFFFFF",
        font=("Gill Sans MT", 32 * -1)
    )

    canvas.create_text(
        86.0,
        195.0,
        anchor="nw",
        text="SMILES of the molecule",
        fill="#FFFFFF",
        font=("Gill Sans MT", 13 * -1)
    )

    canvas.create_text(
        86.0,
        308.0,
        anchor="nw",
        text="Overview of the molecule",
        fill="#FFFFFF",
        font=("Gill Sans MT", 13 * -1)
    )

    canvas.create_text(
        507.0,
        213.0,
        anchor="nw",
        text="Melting point [°C]",
        fill="#FFFFFF",
        font=("Gill Sans MT", 13 * -1)
    )

    canvas.create_text(
        465.0,
        372.0,
        anchor="nw",
        text="More about the program",
        fill="#FFFFFF",
        font=("Gill Sans MT", 15 * -1)
    )

    canvas.create_text(
        465.0,
        401.0,
        anchor="nw",
        text="This program has been developed to predict the \nmelting point of certain organic molecules. \nOur program was trained by Mordred descriptors,\nmore details in",
        fill="#FFFFFF",
        font=("Gill Sans MT", 11 * -1)
    )

    entry_left_container = PhotoImage(
        file=assets_dir / 'input_textbox.png')
    entry_bg_1 = canvas.create_image(
        169.08280181884766,
        229.7234992980957,
        image=entry_left_container
    )
    input_textbox = Entry(
        bd=0,
        bg="#616161",
        fg="#FFFFFF",
        highlightthickness=0
    )
    input_textbox.place(
        x=93.41810750961304,
        y=216.0,
        width=151.32938861846924,
        height=25.446998596191406
    )

    button_left_container = PhotoImage(
        file=assets_dir / 'submit_button.png')
    submit_button = Button(
        image=button_left_container,
        borderwidth=0,
        highlightthickness=0,
        activebackground="#2A2A2A",
        disabledforeground="#FFFFFF",
        
        # Allows you to make the prediction while running the interface
        command= lambda: threading.Thread(target= lambda: getEntry_and_Test(input_textbox)).start(),
        relief="flat"
    )
    submit_button.place(
        x=275.0,
        y=217.0,
        width=53.0,
        height=27.0
    )


    # Defining the labels for the functions presented before, without text
    error_label = tk.Label(window, fg="red", bg="#2A2A2A")
    error_label.place(x=90, y=250)
    answer_label = tk.Label(window, fg="white", bg="#616161")
    answer_label.place(x=573, y=237)

    #Create a Label to display the link
    link = Label(window, text="https://github.com/Moglul/Ppchem-MP-predictor",font=('Gill Sans MT', 12 * -1, UNDERLINE), fg="#FFFFFF", bg ="#2A2A2A", cursor="hand2")
    link.place(x=462, y=452)
    link.bind("<Button-1>", lambda e:
    callback("https://github.com/Moglul/Ppchem-MP-predictor"))


    # If the enter button is pressed, the button submit is pressed, also the animation is launched
    input_textbox.bind('<Return>', lambda x: submit_button.invoke() and button_clicked_animation(submit_button))


    # Closes the loop of the GUI
    window.mainloop()


# Launching of the GUI
if __name__ == "__main__":
    start_gui()