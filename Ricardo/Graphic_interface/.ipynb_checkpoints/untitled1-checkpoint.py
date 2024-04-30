import subprocess
import sys

def create_executable(script_path):
    try:
        subprocess.check_call(["pyinstaller", "--onefile", script_path])
    except Exception as e:
        print("Une erreur s'est produite lors de la création de l'exécutable :", e)

if __name__ == "__main__":
    script_path = "graphic_interface.py"  # Remplacez "graphic_interface.py" par le chemin de votre script Python
    create_executable(script_path)