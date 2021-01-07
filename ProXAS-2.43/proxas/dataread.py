import numpy as np
import pandas as pd
from proxas.headerread import headerread 
		
class dataread:
	def __init__(self,root_i, minr):
		"""
		Read Data
		
		"""
		
	def data_read_qex(self, root_i, minr):
		headerSize, line_bytes, dt, nData = headerread.header_read_qex(self, self.data_file)
	
		root = int(self.calibroots[root_i+3])
		root_end = int(self.calibroots[root_i+4])
		minr = int(root_end - root)

		qex_file = codecs.open(self.data_file+'.qex', 'rb', encoding='cp1252') 

		qex_file.seek(int(headerSize)+(int(line_bytes)*root))
	
		data = np.fromfile(qex_file, dtype=dt, count = minr)
		
		ang = data['encoder']
		
		toggle = int(self.togglevar.get())
		mu = np.zeros((minr, 1))
		if toggle == 0:
			numervar = str(self.sam_numer.get())
			denomvar = str(self.sam_denom.get())
			logvar = int(self.sam_logvar.get())
		else:
			numervar = str(self.ref_numer.get())
			denomvar = str(self.ref_denom.get())
			logvar = int(self.ref_logvar.get())
		
		for j in range(int(self.nChannels)+1):
			if numervar == self.choices[j]:
				numercol = int(j)
			if denomvar == self.choices[j]:
				denomcol = int(j)
		
		if logvar == 1:
			if denomcol == int(len(self.choices)-1):
				mu[:,0] = np.log(data['CH '+str(numercol)])
			else:
				mu[:,0] = np.log(data['CH '+str(numercol)]/data['CH '+str(denomcol)])
		else:
			if denomcol == int(len(self.choices)-1):
				mu[:,0] = data['CH '+str(numercol)]
			else:
				mu[:,0] = data['CH '+str(numercol)]/data['CH '+str(denomcol)]
			
		RawData = pd.DataFrame()
		RawData['ang'] = np.around(ang, decimals=5)
		RawData['mu'] = mu
		
		RawData.dropna(how='any', inplace = True)
		
		RawData = RawData.groupby('ang', as_index=False).mean()
		#RawData.columns = ['ang','mu','sem']
		RawData = RawData[(RawData['ang'] >= np.float(self.min_ang)) & (RawData['ang'] <= np.float(self.max_ang))]
		
		RawData.sort_values(by='ang', ascending=False)
		RawData.reset_index(drop=True)
		
		#print(RawData.shape, np.min(RawData['ang']), np.max(RawData['ang']))
		return np.around(RawData['ang'].values, decimals=5), RawData['mu'].values
		 
	def data_read_bin(self, root_i, minr):
		headerread.header_read_bin(self, self.data_file)
		headerread.header_read_bin_ch(self, self.data_file)
		
		root = self.calibroots[root_i+3]
		root_end = int(self.calibroots[root_i+4])
		minr = int(root_end - root)
		
		g = open(self.data_file+'.bin', 'rb')
		g.seek(int(self.headerSize+(4*self.nChannels*root)))
		data_D = (np.fromfile(g, dtype='f4', count = int(self.nChannels*minr)))
		
		encoder_bin = self.data_file+'_Encoder'
		f = open(encoder_bin+'.bin', 'rb')
		f.seek(int(self.headerSize+(4*root)))
		ang = np.around(np.fromfile(f, dtype='f4', count = int(minr)), decimals=5)
		
		toggle = int(self.togglevar.get())
		mu = np.zeros((minr, 1))
		if toggle == 0:
			numervar = str(self.sam_numer.get())
			denomvar = str(self.sam_denom.get())
			logvar = int(self.sam_logvar.get())
		else:
			numervar = str(self.ref_numer.get())
			denomvar = str(self.ref_denom.get())
			logvar = int(self.ref_logvar.get())
		
		for j in range(int(self.nChannels)+1):
			if numervar == self.choices[j]:
				numercol = int(j)
			if denomvar == self.choices[j]:
				denomcol = int(j)
		
		if logvar == 1:
			if denomcol == int(len(self.choices)-1):
				mu[:,0] = np.log(data_D[numercol::int(self.nChannels)])
			else:
				mu[:,0] = np.log(data_D[numercol::int(self.nChannels)]/data_D[denomcol::int(self.nChannels)])
		else:
			if denomcol == int(len(self.choices)-1):
				mu[:,0] = data_D[numercol::int(self.nChannels)]
			else:
				mu[:,0] = data_D[numercol::int(self.nChannels)]/data_D[denomcol::int(self.nChannels)]
		
		RawData = pd.DataFrame()
		RawData['ang'] = np.around(ang, decimals=5)
		RawData['mu'] = mu
		
		RawData.dropna(how='any', inplace = True)
		
		RawData = RawData.groupby('ang', as_index=False).agg({'mu':['mean','std']})
		RawData.columns = ['ang','mu','sem']
		RawData = RawData[(RawData['ang'] >= np.float(self.min_ang)) & (RawData['ang'] <= np.float(self.max_ang))]
		
		RawData.sort_values(by='ang', ascending=False)
		RawData.reset_index(drop=True)
		
		return np.around(RawData['ang'].values, decimals=5), RawData['mu'].values, RawData['sem'].values