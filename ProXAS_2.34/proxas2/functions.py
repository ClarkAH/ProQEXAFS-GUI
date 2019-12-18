import numpy as np
from scipy.interpolate import Rbf
import numpy_indexed as npi

#Main Algorithm
class polynomials:
	def __init__(self,x,*terms):
		"""
		Defines polynomials
		
		"""
	def constant(x, *terms):
		return terms[0]*np.ones(len(x))
	
	def linear(x, *terms):
		return (terms[0]*x) + terms[1]
	
	def quadratic(x, *terms):
		return (terms[0]*x) + (terms[1]*x**2) + terms[2]
	
	def cubic(x, *terms):
		return (terms[0]*x) + (terms[1]*x**2) + (terms[2]*x**3) + terms[3]
		
	def quartic(x, *terms):
		return (terms[0]*x) + (terms[1]*x**2) + (terms[2]*x**3) + (terms[3]*x**4) + terms[4]
	
	def victoreen(x, *terms):
		f = 1.23986*10**4
		return ((terms[0]*f**3)/(x**3)) - ((terms[1]*f**4)/(x**4)) + terms[2]
		
