import tkinter as tk

def on_button_click():
    label.config(text="Hello, " + entry.get())

# Creation of the window
window = tk.Tk()
window.title("Interface for melting point")

# Ajouter des composants
label = tk.Label(window, text="Entrez votre nom :")
label.pack()

entry = tk.Entry(window)
entry.pack()

button = tk.Button(window, text="Valider", command=on_button_click)
button.pack()

# Lancer la boucle principale
window.mainloop()
