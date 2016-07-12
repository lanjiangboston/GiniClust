#!/usr/bin/env python


# GiniClust/setup.py


# Author: Gregory Giecold
# Affiliation: Harvard University
# Contact: ggiecold@jimmy.harvard.edu


from codecs import open
from os import path

from setuptools import setup


here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README'), encoding = 'utf-8') as f:
    long_description = f.read()
    

setup(name = 'GiniClust',
      version = '0.1.7',
      
      description = "Detecting rare cell-types from single-cell "
          "gene expression data",
      long_description = long_description,
                    
      url = 'https://github.com/lanjiangboston/GiniClust',
      download_url = 'https://github.com/lanjiangboston/GiniClust',
      
      author = 'Lan Jiang, Qian Zhu, Gregory Giecold',
      author_email = "lan_jiang@hms.harvard.edu, qzhu@princeton.edu, "
          "ggiecold@jimmy.harvard.edu",
      maintainer = 'Lan Jiang, Qian Zhu',
      maintainer_email = "lan_jiang@hms.harvard.edu, qzhu@princeton.edu",
      
      license = 'MIT License',
      
      platforms = ('Any',),
      install_requires = ['Gooey>=0.9.2.3', 'setuptools>=24.0.2'],
                          
      classifiers = ['Development Status :: 4 - Beta',
                     'Environment :: Console',
                     'Intended Audience :: End Users/Desktop',
                     'Intended Audience :: Developers',
                     'Intended Audience :: Healthcare Industry',
                     'Intended Audience :: Science/Research',          
                     'License :: OSI Approved :: MIT License',
                     'Natural Language :: English',
                     'Operating System :: MacOS :: MacOS X',
                     'Operating System :: Microsoft',
                     'Operating System :: POSIX',
                     'Programming Language :: Python :: 2.7',
                     'Topic :: Scientific/Engineering',
                     'Topic :: Scientific/Engineering :: Bio-Informatics',
                     'Topic :: Scientific/Engineering :: Visualization',
                     'Topic :: Software Development :: User Interfaces', ],
                   
      packages = ['GiniClust'],
      package_dir = {'GiniClust': 'src/GiniClust'},
      include_package_data = True,

      keywords = "bioinformatics biology clustering genomics "
          "machine-learning outliers PCR qPCR RNASeq single-cell",
                 
      entry_points = {
          'console_scripts': ['GiniClust = GiniClust.__main__:main'],
          }    
)



