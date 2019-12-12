# ProXAS-GUI
Python codes for QXAS processing at Super-XAS

Installation Instructions
The easiest method is to install anaconda python 3 distribution making sure to tick the box to add to system path 
https://repo.anaconda.com/archive/Anaconda3-5.3.0-Windows-x86_64.exe

Once complete there are some additional packages are required. These can be installed easily from the command prompt using the install_script.py contained within the download.

Alternatively these can in installed from command prompt:

conda install -c GSECARS xraylarch

conda install -c conda-forge numpy-indexed

To run the main program from command prompt first navigate to the directory contraining the program then use: python ProXAS-GUI.py
Alternatively if the python interpretter has been set as the default for opening .py files, simply double click or run from command line the ProXAS_v2.34.py file

To run this program the subroutine files are also required (in the same folder as the main script): batch_split_subroutine_v2.3.py, batch_extract_subroutine_v2.9.1.py

Package Requirements:
larch >= 0.9.40, matplotlib >= 2.2.2, numpy_indexed >= 0.3.5, numpy >= 1.15.4, pandas >= 0.22.0, peakutils >= 1.1.0, matplotlib >= 2.2.2, psutil >= 5.4.7, re >= 2.2.1, scipy >= 1.1.0, xraydb >= 1.3

A user manual is supplied but at this stage is under constant revision and update.
