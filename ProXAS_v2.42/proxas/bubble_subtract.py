import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

export = 'n'

while export != 'y':
	folder = r'C:\Users\clark_a\Desktop\output_Pd_K_Cell5_N2sat_KHCO3_-300mV\Export'
	file_avg = r'Cell5_N2sat_KHCO3_-300mV_sam_matrix_Both_normalised_all.dat'
	file_raw = r'Cell5_N2sat_KHCO3_-300mV_sam_matrix_Both_normalised_2.dat'
	
	import_data_avg = pd.read_csv(folder+'/'+file_avg, sep='\t', header=0)
	import_data_raw = pd.read_csv(folder+'/'+file_raw, sep='\t', header=0)
	import_data_raw.drop(['E'], axis=1, inplace = True)
	
	diff_data = import_data_raw.sub(import_data_avg['0'], axis=0)
	diff_sum = np.abs(diff_data.sum(axis = 0))
	
	plt.plot(diff_sum)
	plt.show()
	diff_std = np.std(diff_sum)
	print('absolute difference Std =', diff_std)
	
	msk = np.float(input('input masking value: '))
	
	bubble_list = []
	diff_sum_k = []
	
	for i in range(len(diff_sum)):
		if diff_sum[i] > msk:
			continue
		else:
			bubble_list.append(str(i))
			diff_sum_k.append(diff_sum[i])
			
	print(len(bubble_list))
	print(bubble_list)
	
	ex_bubble = import_data_raw[bubble_list]
	ex_bubble.insert(0, 'E', import_data_avg['E'])
	
	plt.plot(ex_bubble['E'], ex_bubble.iloc[:,1:])
	plt.show()
	
	
	export = input('Export data? y/n: ')
	if export == 'y':
		ex_bubble.to_csv(folder+'/'+file_raw.split('.dat')[0]+'_bubble_subtract.dat', sep='\t', index=False)
	

