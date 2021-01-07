import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
import os

folder = r'C:\Users\clark_a\Documents\Analysis\XPDF\PDF Newton\IR'
file = '210816_1.0-reformat.dat'

import_data = pd.read_csv(folder+'/'+file, sep='\t', header=0)

period = 480
nphase = 181                    

start_period = 1
end_period = 999

phase_delay = 30

w = 2*np.pi/period
max_n = 21

t1 = 140
t2 = 480
d = t1/t2

alpha = 1.05

headers = list(import_data.columns.values)

del headers[0]

headers = [int(i) for i in headers]
headers = np.asarray(headers)

period_average = pd.DataFrame()
try:
	period_average['wn'] = import_data['wn']
except:
	period_average['wn'] = import_data['Wavenumber']
	
def intergrand_function(chi_data,nw,t,phiPSD):
	chi_ft = np.asarray(chi_data*np.sin(nw*t + phiPSD))
	return chi_ft

for j in range(int(max_n)): 
	n = j+1
	print(n, np.sin(n*d*np.pi))
	
	#plt.figure(n)

	if j == 0:
		for i in range(period+1):
			columns = headers[np.int(i+(phase_delay))::(period)]
			columns=list(columns[(columns >= ((period)*start_period)+phase_delay)&(columns < ((period)*(end_period+1))+phase_delay)])
			columns = [str(i) for i in columns]
			data_to_average = import_data[columns]
			
			period_average[str(i)] = data_to_average.mean(axis = 1)
			#if i%10 == 0:
				#plt.plot(period_average['E'], period_average[str(i)])	
					
	x =  np.asarray(range(period_average.shape[1] - 1)) 
	energy = np.asarray(period_average['wn'].values)
	
	chi_data = np.asarray(period_average.drop(['wn'], axis=1).values)
		
	phiPSD_grid = np.linspace(0, 2*np.pi, nphase, endpoint=True)
	phiPSD_grid = phiPSD_grid[0:-1]
	
	ft_out_1D = np.zeros(len(energy))
			
	#ax=plt.subplot(111)
	#ax.text(0.9, 0.9, 'n='+str(n),
    #    horizontalalignment='right',
    #    verticalalignment='top',
    #    transform=ax.transAxes)
		
	output_psd = pd.DataFrame()
	output_psd['wn'] = period_average['wn']
		
	for k in range(len(phiPSD_grid)):
		for i in range(len(energy)):	
			y = (np.sin(n*d*np.pi))*(4*0.5/(n*np.pi))*(2/(period))*intergrand_function(chi_data[i], n*w, x, phiPSD_grid[k])
			ft_out_1D[i] = (2*alpha)*simps(y,x)
			
		output_psd[str(k)] = ft_out_1D.real
			
		#if phiPSD_grid[k] <= np.pi:
		#	ax.plot(output_psd['wn'], output_psd[str(k)], label='\u03C6'+' = '+str(int(np.around(180*(phiPSD_grid[k].real)/(np.pi), decimals=1)))+'\u00B0')
		#else:
		#	ax.plot(output_psd['wn'], output_psd[str(k)], ls='--', label='\u03C6'+' = '+str(int(np.around(180*(phiPSD_grid[k].real)/(np.pi), decimals=1)))+'\u00B0')
	
	#ax.legend(loc='upper left', ncol=2)
	#		
	#plt.ion()
	#
	#if n == max_n:
	#	plt.ioff()
	#	
	#plt.show()
	
	if not os.path.exists(folder+'/output_psd'):
		os.makedirs(folder+'/output_psd')
		
	if j == 0:
		output_psd_tot = output_psd
	else:
		for k in range(len(phiPSD_grid)):
			output_psd_tot[str(k)] = output_psd[str(k)] + output_psd_tot[str(k)]
	
	output_psd.to_csv(folder+'/output_psd/'+str(os.path.splitext(file)[0])+'_psd_n='+str(n)+'_s'+str(start_period)+'_e'+str(end_period)+'_pd'+str(phase_delay)+'.dat', sep='\t', index=False)
	output_psd_tot.to_csv(folder+'/output_psd/'+str(os.path.splitext(file)[0])+'_psd_tot_s'+str(start_period)+'_e'+str(end_period)+'_pd'+str(phase_delay)+'.dat', sep='\t', index=False)
	period_average.to_csv(folder+'/output_psd/'+str(os.path.splitext(file)[0])+'_s'+str(start_period)+'_e'+str(end_period)+'_period_average_pd'+str(phase_delay)+'.dat', sep='\t', index=False)
	
plt.figure(n+1)
ax=plt.subplot(111)
ax.text(0.9, 0.9, 'n='+str(n),
       horizontalalignment='right',
       verticalalignment='top',
       transform=ax.transAxes)
	
for k in range(len(phiPSD_grid)):
	if phiPSD_grid[k] <= np.pi:
		ax.plot(output_psd_tot['wn'], output_psd_tot[str(k)], label='\u03C6'+' = '+str(int(np.around(180*(phiPSD_grid[k].real)/(np.pi), decimals=1)))+'\u00B0')
	else:
		ax.plot(output_psd_tot['wn'], output_psd_tot[str(k)], ls='--', label='\u03C6'+' = '+str(int(np.around(180*(phiPSD_grid[k].real)/(np.pi), decimals=1)))+'\u00B0')
		
ax.legend(loc='upper left', ncol=2)
plt.show()