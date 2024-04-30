
# This file was generated by the Tkinter Designer by Parth Jadhav
# https://github.com/ParthJadhav/Tkinter-Designer


from pathlib import Path

# from tkinter import *
# Explicit imports to satisfy Flake8
from tkinter import Tk, Canvas, Entry, Text, Button, PhotoImage


OUTPUT_PATH = Path(__file__).parent
ASSETS_PATH = OUTPUT_PATH / Path(r"C:\Users\rikim\git\Ppchem-MP-predictor\Ricardo\Graphic_interface\build\assets\frame0")


def relative_to_assets(path: str) -> Path:
    return ASSETS_PATH / Path(path)


window = Tk()

window.geometry("800x600")
window.configure(bg = "#1A1A1A")


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
image_image_1 = PhotoImage(
    file=relative_to_assets("image_1.png"))
image_1 = canvas.create_image(
    220.0,
    337.0,
    image=image_image_1
)

image_image_2 = PhotoImage(
    file=relative_to_assets("image_2.png"))
image_2 = canvas.create_image(
    603.0,
    249.0,
    image=image_image_2
)

image_image_3 = PhotoImage(
    file=relative_to_assets("image_3.png"))
image_3 = canvas.create_image(
    603.0,
    423.0,
    image=image_image_3
)

entry_image_1 = PhotoImage(
    file=relative_to_assets("entry_1.png"))
entry_bg_1 = canvas.create_image(
    169.08280181884766,
    229.7234992980957,
    image=entry_image_1
)
entry_1 = Entry(
    bd=0,
    bg="#616161",
    fg="#000716",
    highlightthickness=0
)
entry_1.place(
    x=93.41810750961304,
    y=216.0,
    width=151.32938861846924,
    height=25.446998596191406
)

button_image_1 = PhotoImage(
    file=relative_to_assets("button_1.png"))
button_1 = Button(
    image=button_image_1,
    borderwidth=0,
    highlightthickness=0,
    command=lambda: print("button_1 clicked"),
    relief="flat"
)
button_1.place(
    x=275.0,
    y=216.0,
    width=53.0,
    height=27.0
)

canvas.create_text(
    86.0,
    195.0,
    anchor="nw",
    text="SMILE of the molecule",
    fill="#FFFFFF",
    font=("Quicksand SemiBold", 13 * -1)
)

image_image_4 = PhotoImage(
    file=relative_to_assets("image_4.png"))
image_4 = canvas.create_image(
    400.0,
    46.0,
    image=image_image_4
)

canvas.create_text(
    86.0,
    308.0,
    anchor="nw",
    text="Overview of the molecule",
    fill="#FFFFFF",
    font=("Quicksand SemiBold", 13 * -1)
)

image_image_5 = PhotoImage(
    file=relative_to_assets("image_5.png"))
image_5 = canvas.create_image(
    219.0,
    401.0,
    image=image_image_5
)

image_image_6 = PhotoImage(
    file=relative_to_assets("image_6.png"))
image_6 = canvas.create_image(
    735.0,
    45.0,
    image=image_image_6
)

canvas.create_text(
    18.0,
    25.0,
    anchor="nw",
    text="Melting point predictor\n\n\n",
    fill="#FFFFFF",
    font=("Quicksand SemiBold", 32 * -1)
)

canvas.create_text(
    507.0,
    213.0,
    anchor="nw",
    text="Melting point [K]",
    fill="#FFFFFF",
    font=("Quicksand SemiBold", 13 * -1)
)

image_image_7 = PhotoImage(
    file=relative_to_assets("image_7.png"))
image_7 = canvas.create_image(
    603.0,
    248.0,
    image=image_image_7
)

canvas.create_text(
    465.0,
    401.0,
    anchor="nw",
    text="This program has been developed to predict the melting point of certain organic molecules. Our program was trained by Mordred descriptors, more details in https://github.com/Moglul/Ppchem-MP-predictor.\n",
    fill="#FFFFFF",
    font=("Quicksand SemiBold", 11 * -1)
)

canvas.create_text(
    465.0,
    377.0,
    anchor="nw",
    text="More about the program",
    fill="#FFFFFF",
    font=("Quicksand SemiBold", 13 * -1)
)
window.resizable(False, False)
window.mainloop()
