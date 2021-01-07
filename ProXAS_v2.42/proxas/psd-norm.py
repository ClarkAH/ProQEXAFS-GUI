import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

folder = r'C:\Users\clark_a\Documents\Analysis\XPDF\PDF Newton\IR\output_psd'
file = '210816_1.0-reformat_psd_tot_s1_e999_pd0.dat'

import_data = pd.read_csv(folder+'/'+file, sep='\t', header=0)

dtheta = 2

wn = import_data['wn']
import_data = import_data.drop(['wn'], axis=1)
max_array = []

for i in range(import_data.shape[0]):
	import_data.iloc[i,:] = import_data.iloc[i,:]/np.max(import_data.iloc[i,:])

	inphase = (dtheta*np.argmax(import_data.iloc[i,:].values)-180) % (360) + 180
	
	max_array.append(inphase)
	
for i in range(import_data.shape[1]):
	plt.plot(wn, import_data.iloc[:,i])

plt.xlim(np.max(wn), np.min(wn))	
plt.show()

max_array = np.asarray(max_array)

plt.plot(wn, max_array, 'b')
plt.xlim(np.max(wn), np.min(wn))	
plt.show()

inphase_pd = pd.DataFrame()
inphase_pd['wn'] = wn
inphase_pd['inphase'] = max_array
inphase_pd.to_csv(folder+'/'+file.split('.dat')[0]+'_inphase_angles.dat', sep='\t', header=False)