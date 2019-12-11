import numpy as np
import matplotlib.pyplot as plt
import random, sys
import scipy.optimize as optimize
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

#Main Algorithm
class lcf:
	def __init__(self,D, S, constrain):
		"""Perform SIMPLISMA
		
		D = data matrix
		
		S = spectra
		
		constrained LCF
		
		"""
		print('Performing Simplisma')
		
	def lcf(D,S,constrain):
		print(constrain)

		def model(S, *wt0):
			s = np.zeros(S.shape[0])
			for i in range(len(wt0)):
				s = s + wt0[i]*S[:,i]
			if constrain == True:
				return s/sum(wt0)
			else:
				return s
			
		CS = np.zeros(D.shape)
		if constrain == True:
			lower = np.zeros(S.shape[1])
			upper = np.ones(S.shape[1])*1
		else:
			lower = -np.ones(S.shape[1])*1
			upper = np.ones(S.shape[1])*1
		
		LOF = np.zeros(np.shape(D)[1])
		
		C_u = np.dot(D.T, np.linalg.pinv(S.T))
		
		if constrain == True:
			for i in range(np.shape(D)[1]):
				wt0 = np.abs(C_u[i,:])/sum(np.abs(C_u[i,:]))			
				popt,pcov = optimize.curve_fit(model,S,D[:,i],p0=wt0, method='trf', bounds=(lower, upper))
				CS[:,i] = model(S, *popt)
				LOF[i] = 100*np.sqrt(np.sum(D[:,i]-CS[:,i])**2/np.sum(D[:,i])**2)
				sys.stdout.write("\r" + 'Percent complete (%): '+str(np.around(100*(i+1)/np.shape(D)[1], decimals=0))+' : '+str(np.around(LOF[i], decimals=2)))
				sys.stdout.flush()
				wt0 = popt
			
			C_r = np.dot(CS.T, np.linalg.pinv(S.T))
				
		if constrain == False:
			CS = np.dot(C_u, S.T).T
			for i in range(np.shape(D)[1]):
				#LOF[i] = 100*np.sqrt(np.sum(D[:,i]-CS[:,i])**2/np.sum(D[:,i])**2)
				sys.stdout.write("\r" + 'Percent complete (%): '+str(np.around(100*(i+1)/np.shape(D)[1], decimals=0))+' : '+str(np.around(LOF[i], decimals=2)))
				sys.stdout.flush()
				
			LOF = 100*np.sqrt(((D-CS)**2)/((D)**2))
			
			C_r = C_u
		
		return C_r, LOF