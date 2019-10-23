import subprocess, shutil, os

subprocess.call(['conda', 'install', '-c', 'GSECARS', 'xraylarch'])
subprocess.call(['conda', 'install', '-c', 'conda-forge', 'numpy-indexed'])
subprocess.call(['conda', 'update', 'scipy'])
subprocess.call(['conda', 'update', 'pandas'])
