# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 11:02:56 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""

import setuptools 

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(name="SPTAnalyser",
                 version="230714",
                 author="Johanna Rahm and Sebastian Malkusch",
                 author_email="johanna-rahm@gmx.de",
                 description="a package for analysing SPT data",
                 long_description=long_description,
                 long_description_content_type="text/markdown",
                 licence="GNU GENERAL PUBLIC LICENSE",
                 url="https:/github.com/JohannaRahm/SPTAnalyser",
                 packages=setuptools.find_packages(),
                 install_requires=["numpy>=1.20.0",
                                   "matplotlib",
                                   "seaborn",
                                   "scipy",
                                   "h5py",
                                   "tqdm",
                                   "ipywidgets",
                                   "IPython",
                                   "hide_code",
                                   "tornado",
                                   "pandas",
                                   "graphviz"],
                   classifiers=["Programming Language :: Python :: 3.9",
                                  "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
                                  "Operating System :: OS Independent",],)
