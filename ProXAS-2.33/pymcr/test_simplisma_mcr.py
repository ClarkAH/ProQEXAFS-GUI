import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random

from pymcr.simplisma import svd, simplisma
from pymcr.mcr import McrAls
from pymcr.regressors import OLS, NNLS
from pymcr.constraints import ConstraintNonneg, ConstraintNorm

#Number of Spectral Components
nSVD = 10
#Allowed Noise Percentage
noise = 5	
#manual
manual = False

D = np.asarray((pd.read_csv('Total_MCR_CuSSZ13.dat', sep='\t', header=None)).values)

if manual == False:
	#Run SVD
	eigens, explained_variance_ratio = svd.svd(D, nSVD)
	nPure = np.int(input('Number of Principle Components :'))
	#Run Simplisma
	S, C_u, C_c = simplisma.pure(D.T, nPure, noise, True)
else:
	S = np.asarray((pd.read_csv('sopt_5c2.dat', sep='\t', header=None)).values).T
	

#Run MCR
mcrals = McrAls(max_iter=50, st_regr='NNLS', c_regr='NNLS', 
                c_constraints=[ConstraintNonneg(), ConstraintNorm()])

mcrals.fit(D, ST=S.T, verbose=True)
print('\nFinal MSE: {:.7e}'.format(mcrals.err[-1]))

plt.subplot(2, 1, 1)
plt.plot(mcrals.ST_opt_.T)

plt.subplot(2, 1, 2)
plt.plot(mcrals.C_opt_)

plt.show()