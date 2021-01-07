import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
import os

folder = r'E:\ProXAS-2\output_Ir_L3_IrO2_WEpr1_CEpr1_RE13_1V_1pt6V_500cycles\Export'
file = 'IrO2_WEpr1_CEpr1_RE13_1V_1pt6V_500cycles_sam_matrix_Up.dat'
import_data = pd.read_csv(folder+'/'+file, sep='\t', header=0)

period = 100
nphase = 25

start_period = 140
end_period = 150

phase_delay = 50

w = 2*np.pi/period
max_n = 1

alpha = 1.05

headers = list(import_data.columns.values)

del headers[0]

headers = [int(i) for i in headers]
headers = np.asarray(headers)

period_average = pd.DataFrame()
try:
	period_average['E'] = import_data['E']
except:
	period_average['E'] = import_data['Energy']
	
def intergrand_function(chi_data,nw,t,phiPSD):
	chi_ft = np.asarray(chi_data*np.sin(nw*t + phiPSD))
	return chi_ft

for j in range(int((max_n-1)/2)+1): 
	n = (2*j)+1
	
	plt.figure(n)

	if j == 0:
		for i in range(period+1):
			columns = headers[(i+phase_delay)::(period+1)]
			columns=list(columns[(columns >= (period+1)*start_period)&(columns <= (period+1)*end_period)])
			columns = [str(i) for i in columns]
			data_to_average = import_data[columns]
			
			period_average[str(i)] = data_to_average.mean(axis = 1)
			#if i%10 == 0:
				#plt.plot(period_average['E'], period_average[str(i)])	
					
	x =  np.asarray(range(period_average.shape[1] - 1)) 
	energy = np.asarray(period_average['E'].values)
	
	chi_data = np.asarray(period_average.drop(['E'], axis=1).values)
		
	phiPSD_grid = np.linspace(0, 2*np.pi, nphase, endpoint=True)
	ft_out_1D = np.zeros(len(energy))
			
	ax=plt.subplot(111)
	ax.text(0.9, 0.9, 'n='+str(n),
        horizontalalignment='right',
        verticalalignment='top',
        transform=ax.transAxes)
		
	output_psd = pd.DataFrame()
	output_psd['E'] = period_average['E']
		
	for k in range(len(phiPSD_grid)):
		for i in range(len(energy)):	
			y = (2/(period))*intergrand_function(chi_data[i], n*w, x, phiPSD_grid[k])
			ft_out_1D[i] = (2*alpha)*(4/np.pi)*simps(y,x)/((2*n)-1)
			
		output_psd[str(k)] = ft_out_1D.real
			
		if phiPSD_grid[k] <= np.pi:
			ax.plot(output_psd['E'], output_psd[str(k)], label='\u03C6'+' = '+str(int(np.around(180*(phiPSD_grid[k].real)/(np.pi), decimals=1)))+'\u00B0')
		else:
			ax.plot(output_psd['E'], output_psd[str(k)], ls='--', label='\u03C6'+' = '+str(int(np.around(180*(phiPSD_grid[k].real)/(np.pi), decimals=1)))+'\u00B0')
	
	ax.legend(loc='upper left', ncol=2)
			
	plt.ion()
	
	if n == max_n:
		plt.ioff()
		
	plt.show()
	
	if not os.path.exists(folder+'/output_psd'):
		os.makedirs(folder+'/output_psd')
	
	output_psd.to_csv(folder+'/output_psd/'+str(os.path.splitext(file)[0])+'_psd_n='+str(n)+'_s'+str(start_period)+'_e'+str(end_period)+'.dat', sep='\t', index=False)
	period_average.to_csv(folder+'/output_psd/'+str(os.path.splitext(file)[0])+'_s'+str(start_period)+'_e'+str(end_period)+'_period_average.dat', sep='\t', index=False)