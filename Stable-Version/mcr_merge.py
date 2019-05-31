import pandas as pd
import os, shutil, sys

#Folder name
folder='Ti'

#Extension of output files
extension='mcr.dat'

#Creates list of files in folder
cfn = os.walk(folder, topdown=True, onerror=None, followlinks=False)

mcr_files = []
mcr_files_full = []

for root_o, dir_o, files_o in cfn:
	for file in files_o:
		if file.endswith(extension):
			mcr_files.append(file)
			mcr_files_full.append(root_o + '/' + file)
						
print(mcr_files)

total_mcr = pd.DataFrame()
index_list = pd.DataFrame()
index_list['file'] = mcr_files
index_start = []
index_end = []

for i in range(len(mcr_files)):
	file = mcr_files_full[i]
	readdata = pd.read_csv(file, sep='\t', header=None)
	if i == 0:
		total_mcr = readdata
		index_start.append(0)
	else:
		index_start.append(len(total_mcr.index))
		total_mcr = pd.concat([total_mcr, readdata])
	index_end.append(len(total_mcr.index)-1)
	print('percent complete:', 100*(i+1)/len(mcr_files),' %')

print('saving mcr file')
index_list['index start'] = index_start
index_list['index end'] = index_end
total_mcr.to_csv('Total_MCR.dat', header=None, index=False, sep='\t')
index_list.to_csv('MCR_Exp_Index.dat', header=None, index=False, sep='\t')

print(index_list)
