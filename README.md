# ProXAS-GUI
Python codes for QXAS processing at Super-XAS

Installation Instructions
The easiest method is to install anaconda python 3 distribution making sure to tick the box to add to system path 
https://repo.anaconda.com/archive/Anaconda3-5.3.0-Windows-x86_64.exe

Once complete there are some additional packages are required. These can be installed easily from the command prompt using the following commands:

conda install –c conda-forge numpy-indexed

conda update pandas

conda update scipy

conda update -yc GSECARS xraylarch

Finally, the program internally uses some data that is stored in a Xray Database for elements. Download the following:
https://github.com/scikit-beam/XrayDB
Extract the ZIP file. From command prompt change to the directory where the XrayDB master file is. Within there is a folder called python. Execute the script from command prompt with the following command:
Python setup.py install

To run the main program from command prompt use: python ProXAS-GUI.py
Alternatively if the python interpretter has been set as the default for opening .py files, simply double click or run from command line the ProXAS_v2.15.py file

To run this program the subroutine files are also required (in the same folder as the main script): batch_split_subroutine_v2.1.py, batch_extract_subroutine_v2.8.py

Package Requirements:
larch >= 0.9.40, matplotlib >= 2.2.2, numpy_indexed >= 0.3.5, numpy >= 1.15.4, pandas >= 0.22.0, peakutils >= 1.1.0, matplotlib >= 2.2.2, psutil >= 5.4.7, re >= 2.2.1, scipy >= 1.1.0, xraydb >= 1.3
