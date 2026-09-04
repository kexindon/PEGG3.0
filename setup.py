from setuptools import setup
import setuptools

setup(
    include_package_data=True,
    name = 'pegg',
    author = 'Samuel Gould',
    author_email = 'samgould@mit.edu',
    url = 'https://github.com/kexindon/PEGG3.0',
    version = '3.0.1',
    description = 'Prime Editing Guide Generator, with silent bystander mutations',
    long_description = open('README.rst', encoding='utf-8').read(),
    long_description_content_type = 'text/x-rst',
    #package_dir = {'pegg': ''},
    #packages = setuptools.find_packages(), #['pegg'],
    #package_dir = {'pegg': '', 'azimuth': 'pegg/bin/Azimuth-2.0/azimuth'},
    packages = ['pegg'],
    #data loaded at runtime: the scoring models (.pkl/.pickle), the canonical
    #transcript tables (.json), and the safe-target / non-targeting sets
    package_data = {'pegg': ['*.json', '*.pkl', '*.pickle', '*.csv', '*.txt']},
    #py_modules= ["pegg.prime", "pegg.base", "pegg.library", "pegg.crisporEffScores"],

    install_requires = ["Bio>=1.4.0",
        "cyvcf2==0.30.18",
        "matplotlib>=3.5.1",
        "mock>=4.0.3",
        "numpy>=1.21.5,<2",
        "pandas>=1.4.2",
        "seaborn>=0.11.2",
        "setuptools>=61.2.0",
        "Sphinx>=4.4.0",
        "scikit-learn==1.1.1",
        "regex>=2023.8.8",
        "gffutils>=0.11",
    ],

    classifiers=[
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.9",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent"
        ]
)