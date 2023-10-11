import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

export = 'n'

while export != 'y':
	folder = r'Z:\11_2020\Dino\BSCFW622\Try4_WE16D_CE16_REN1_W\output_W_L3_BSCFW622_W_Try4_WE16D_CE16_REN1_916mV_17__1000mV__1550mV_MEXAS\Export'
	file_avg = r'BSCFW622_W_Try4_WE16D_CE16_REN1_916mV_17__1000mV__1550mV_MEXAS_ref_matrix_Up_normalised_all.dat'
	file_raw = r'BSCFW622_W_Try4_WE16D_CE16_REN1_916mV_17__1000mV__1550mV_MEXAS_ref_matrix_Up_normalised_4.dat'
	
	import_data_avg = pd.read_csv(folder+'/'+file_avg, sep='\t', header=0)
	import_data_raw = pd.read_csv(folder+'/'+file_raw, sep='\t', header=0)
	import_data_raw.drop(['E'], axis=1, inplace = True)
	
	diff_data = import_data_raw.sub(import_data_avg['0'], axis=0)
	diff_max = np.abs(diff_data.max(axis = 0))
	
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(diff_max)
	ax.set_xticklabels([])
	plt.show()
	diff_std = np.std(diff_max)
	print('absolute difference Std =', diff_std)
	
	msk = np.float(input('input masking value: '))
	
	bubble_list = []
	diff_max_k = []
	
	for i in range(len(diff_max)):
		if diff_max[i] > msk:
			if i > 100:
				test = import_data_raw[str(i-100)]
				bubble_list.append(str(i-100))
				diff_max_k.append(diff_max[i-100])
			else:
				bubble_list.append(str(i+100))
				diff_max_k.append(diff_max[i+100])
		else:
			bubble_list.append(str(i))
			diff_max_k.append(diff_max[i])
			
	print(len(bubble_list))
	
	ex_bubble = import_data_raw[bubble_list]
	shape = ex_bubble.shape[1]
	col_list = []
	for i in range(shape):
		col_list.append(str(i))
		
	ex_bubble.columns = col_list
	ex_bubble.insert(0, 'E', import_data_avg['E'])
	
	plt.plot(ex_bubble['E'], ex_bubble.iloc[:,1:])
	plt.show()
	
	
	export = input('Export data? y/n: ')
	if export == 'y':
		ex_bubble.to_csv(folder+'/'+file_raw.split('.dat')[0]+'_bubble_subtract.dat', sep='\t', index=False)
	

