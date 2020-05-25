import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps

folder = r'F:\XAS\4th_beamtime_Fe_PtFe\output_Fe_K_PtFe_4_Tuesday_monday_catalyst_90_co_cut_cycles\Export'
file = 'PtFe_4_Tuesday_monday_catalyst_90_co_cut_cycles_sam_matrix_Up_normalised_5_average_2_average_1_average_1.dat'
import_data = pd.read_csv(folder+'/'+file, sep='\t', header=0)

period = 60
nphase = 25
start_period = 0
phase_delay = 0

w = 2*np.pi/(period)
n = 1

headers = list(import_data.columns.values)

del headers[0]

headers = [int(i) for i in headers]
headers = np.asarray(headers)

period_average = pd.DataFrame()
period_average['E'] = import_data['E']

plt.subplot(211)

for i in range(period+1):
	columns = headers[(i+phase_delay)::(period+1)]
	columns=list(columns[(columns >= (period+1)*start_period)])
	columns = [str(i) for i in columns]
	data_to_average = import_data[columns]
	
	period_average[str(i)] = data_to_average.mean(axis = 1)
	if i%10 == 0:
		plt.plot(period_average['E'], period_average[str(i)])

def intergrand_function(chi_data,nw,t,phiPSD):
	chi_ft = np.asarray(chi_data*np.sin(nw*t + phiPSD))
	return chi_ft		
				
x =  np.asarray(range(period_average.shape[1] - 1)) 
energy = np.asarray(period_average['E'].values)

chi_data = np.asarray(period_average.drop(['E'], axis=1).values)
	
phiPSD_grid = np.linspace(0, 2*np.pi, nphase, endpoint=True)
ft_out_1D = np.zeros(len(energy))
		
ax=plt.subplot(212)
output_psd = pd.DataFrame()
output_psd['E'] = period_average['E']
	
for k in range(len(phiPSD_grid)):
	for i in range(len(energy)):	
		y = (2/(period))*intergrand_function(chi_data[i], n*w, x, phiPSD_grid[k])
		ft_out_1D[i] = simps(y,x)
		
	output_psd[str(k)] = ft_out_1D.real
		
	if phiPSD_grid[k] <= np.pi:
		ax.plot(output_psd['E'], output_psd[str(k)], label='\u03C6'+' = '+str(int(np.around(180*(phiPSD_grid[k].real)/(np.pi), decimals=1)))+'\u00B0')
	else:
		ax.plot(output_psd['E'], output_psd[str(k)], ls='--', label='\u03C6'+' = '+str(int(np.around(180*(phiPSD_grid[k].real)/(np.pi), decimals=1)))+'\u00B0')

ax.legend(loc='upper left', ncol=2)
		
#plt.subplot(313)
#plt.plot(period*phiPSD_grid/(2*np.pi))
plt.show()

output_psd.to_csv(folder+'/output_psd.dat', sep='\t', index=False)
period_average.to_csv(folder+'/period_average_psd.dat', sep='\t', index=False)

				

