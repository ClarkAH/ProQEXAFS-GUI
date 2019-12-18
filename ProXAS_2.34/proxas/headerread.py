#Main Algorithm
class headerread:
	def __init__(self,filename):
		"""
		Defines polynomials
		
		"""
		
	def header_read_bin(self,filename):
		#print('Reading Binary Header')
		encoder_bin = filename+'_Encoder'
		f = open(encoder_bin+'.bin', 'rb')
		f.seek(0)
		data = np.fromfile(f, dtype=np.int32, count = 2)
		self.headerSize = data[0]
		fileVersion = data[1]
		data = np.fromfile(f, dtype=np.int64, count = 1)
		DT = str(data[0])
		data = np.fromfile(f, dtype='f4', count = 2)
		self.AdcClock_Hz = data[1]
		self.samplingvar.set(float(self.AdcClock_Hz/1E6))
		self.DacClock_Hz = 1
		self.nData = int((os.path.getsize(encoder_bin+'.bin') - self.headerSize)/4)

	def header_read_bin_ch(self,filename):
		g = open(filename+'.bin', 'rb')
		data = np.fromfile(g, dtype=np.int32, count = 2)
		self.headerSize = data[0]
		data = np.fromfile(g, dtype=np.int64, count = 1)
		data = np.fromfile(g, dtype='f4', count = 2)
		self.nChannels = int(data[0])
		 
	def header_read_qex(self, filename):
		header_lines = []
		with codecs.open(filename+'.qex', 'rb', encoding='cp1252', errors='ignore') as qex_file:
			line = qex_file.readline().strip('\r\n').replace('# ', '')
			while not '_Header_End_'in line:
				header_lines.append(line)
				line = qex_file.readline().strip('\r\n').replace('# ', '')
			
		#search headerlines for keywords
		self.headerSize = np.int([x for x in header_lines if 'FileHeaderSize_byte' in x][0].split(': ')[1])
		self.nColumns = np.int([x for x in header_lines if 'AdcNumberColumnsInDataFile' in x][0].split(': ')[1])
		self.nChannels = np.int([x for x in header_lines if 'AdcNumberChannelsStored' in x][0].split(': ')[1])
		DataLineFormat = ([x for x in header_lines if 'DataLineFormat' in x][0].split(': ')[1]).split(', ')
		DataLineLabels = ([x for x in header_lines if 'DataLineLabels' in x][0].split(': ')[1]).split(', ')
	
		#Interpret DataTypes to extract number of bytes per line
		self.line_bytes = int(sum([int(format_item[1]) for format_item in [re.split(r'(\d+)', s) for s in DataLineFormat]])/8)
		
		#calculate the number of datapoints in file
		self.nLines = int((os.path.getsize(filename+'.qex') - self.headerSize)/(self.line_bytes))
		self.nData = int(self.nColumns*self.nLines)
		
		d_types = [('encoder', 'f4'), ('time', 'u2')]
		for i in range(self.nChannels):
			d_types.append(tuple(('CH '+str(i), 'f4')))
	
		self.dt = np.dtype(d_types)	
	
		return self.headerSize, self.line_bytes, self.dt, self.nData