import pandas as pd
import os

folder = 'E:\ProXAS-2\output_Fe_K_D123_activation_35_to_450'
file = 'Fe_K_D123_activation_35_to_450'
job_length = 9997

ProcessSamUsevar = 1

df_sam = pd.DataFrame()

updownstring = 'Both'
	
for i in range(job_length):
	try:
		data_read = pd.read_csv(folder+'/Export/Individual_Both/'+str(int(i))+'.dat', sep = '\t', header=0)
		df_sam[str(i)] = data_read['mu_sam']
	except:
		print('missing')
	print(i)
	
try:
	if ProcessSamUsevar == 1:
		df_sam.insert(0, 'E', data_read['E'])
		df_sam.to_csv(folder+'/Export/'+file+'_sam_matrix_'+updownstring+'.dat', sep='\t', index=False)
		print('Sample matrix saved')
except:
	print('failed to save sample matrix')