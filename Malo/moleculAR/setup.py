from setuptools import setup, find_packages

setup(
    name='moleculAR',
    version='0.1',
    author='Malo Gfeller',
    author_email='malo.gfeller@epfl.ch',
    description='A cheminformatics machine learning package',
    packages=find_packages(),
    install_requires=[
        'rdkit', 
        'pandas',
        'scikit-learn',
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)
