import numpy as np
from scipy.interpolate import Rbf
import numpy_indexed as npi

#Main Algorithm
class interpolate:
	def __init__(self,x,y,xnew):
		"""Perform Interpolation
		
		x = xdata
		
		y = ydata
		
		xnew = new xgrid
		
		"""
		print('Performing Interpolation')
		
	def localised(x,y,xnew):
		windowloc=200
		x = np.asarray(x)
		y = np.asarray(y)
		
		if x[0] > x[-1]:
			x=x[::-1]
			y=y[::-1]
		
		ynew = np.zeros(len(xnew))
		
		#Localised Radial Basis functions by looping through subset of data
		localisation = int(np.ceil(len(x)/windowloc))
		factor = np.ceil(len(x)/localisation)
		for idloc in range(localisation):
			#localisation allows for some overlapping to improve the consistency upon reforming total data
			#determines if the sub-data is the last sub-data group
			if idloc == localisation-1:
				start = int((idloc*factor)-5)
				end = int(len(x))
				#determines if the sub-data is the first sub-data group
			elif idloc == 0:
				start = 0
				end = (int(1*factor)+5)
			else:
				start = int((idloc*factor)-5)
				end = int(((idloc+1)*factor)+5)
		
			ydata_loop = y[start:end]
			xdata_loop = x[start:end]
		
			#selects the region of x data to extract the new y data
			xnew_cut_idx = npi.indices(xnew, xnew[(xnew >= xdata_loop[0]) & (xnew < xdata_loop[-1])])
			xnew_cut=xnew[xnew_cut_idx]

			#caries out the radial basis function interpolation on the localised sub-data
			rbf = Rbf(xdata_loop, ydata_loop, function='linear', smooth=0)
			ynew[xnew_cut_idx] = rbf(xnew_cut)
		
		return [xnew, ynew]	