import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)
import numpy as np
import pandas as pd
import time, psutil, ast, sys, os, codecs, re
import scipy.signal as signal
from scipy.optimize import curve_fit
from scipy.interpolate import Rbf
import numpy_indexed as npi
from scipy import special
from scipy.optimize import curve_fit
import multiprocessing as mp

################################################################################  
def spawn(options):
	options = ast.literal_eval(options)
	num_spectra, data_file, folder, file_type, n_cpus, ProcessSamvar, ProcessRefvar, NormUsevar, InterpUsevar, AutoAlignUsevar, edgestepvar, Filtervar, bf_scroll_scale, eminvar, emaxxanesvar, emaxexafsvar, constkUsevar, kstepvar, estepvar, sam_numer, sam_denom, sam_logvar, ref_numer, ref_denom, ref_logvar, FlatUsevar, updownvar, edgestepwidth = options

	print('Processing '+options[1])
	procs = list()
	
	if n_cpus > num_spectra:
		n_cpus = num_spectra
		
	rootsfile = rootsfile_check(folder)
	calibrated = calibfile_check(folder)
	normalized = normfile_check(folder, NormUsevar)
	if 0.5*psutil.cpu_count() < n_cpus:
		logical = np.asarray(range(1, 2*int(n_cpus - psutil.cpu_count()/2), 2))
		print('Core Affinity Logical:', logical)
		physical = np.asarray(range(0, 2*int(psutil.cpu_count()/2), 2))
		print('Core Affinity Physical:', physical)
		core_affinity_list = np.concatenate((physical, logical))
	else:
		core_affinity_list = np.asarray(range(0, 2*int(n_cpus), 2))
		print('Core Affinity Physical:', core_affinity_list)
	sys.stdout.flush()    
	
	if (rootsfile == True) & (calibrated == True):
		if ((float(NormUsevar) == 1) & (normalized == True)) or (float(NormUsevar) == 0):
			q = mp.Queue()
			for i in range(int(options[0])):
				q.put(i)
			for cpu_i in range(n_cpus):
				d = dict(affinity = int(core_affinity_list[cpu_i].item()), q = q, ID = cpu_i, options = options)
				p = mp.Process(target=wrapper_targetFunc, kwargs=d)
				p.start()
				procs.append(p)
			for p in procs:
				p.join()
			print('Returned to Main')
	sys.stdout.flush()
	sys.exit(0)

################################################################################  	
def rootsfile_check(folder):
	try:
		global roots_file_data
		roots_file_data = (pd.read_csv(folder+'/roots.dat', sep='\t', header=None)).values
		return True
	except:
		print('Encoder Analysis Required')
		return False

################################################################################  
def calibfile_check(folder):
	try:
		global calibration_values
		calibration_values = (pd.read_csv(folder+'/calibration.dat', sep='\t', header=None)).values
		return True
	except:
		print('Cannot Find Calibration Parameter File')
		return False

################################################################################          
def normfile_check(folder, NormUsevar):
	try:
		global normalisation_values
		normalisation_values = (pd.read_csv(folder+'/Normalisation.dat', sep='\t', header=None)).values
		return True
	except:
		if NormUsevar == 1:
			print('Cannot Find Normalisation Parameter File')
		return False

################################################################################  
def wrapper_targetFunc(affinity, q, ID, options):
	proc = psutil.Process()  # get self pid
	proc.cpu_affinity([affinity])
	print('Process-'+str(ID))
	sys.stdout.flush()
	etd = time.time()
	targetFunc(q, etd, ID, options)
	print('Process-'+str(ID)+' Task Complete')
	if ID == 0:
		print('Percent Complete =',100)
		print('Returning to Main')
	sys.stdout.flush()
	os._exit(0)
 
################################################################################   
def targetFunc(q, etd, ID, options):
	global num_spectra, data_file, folder, file_type, n_cpus, ProcessSamvar, ProcessRefvar, NormUsevar, InterpUsevar, AutoAlignUsevar, edgestepvar, Filtervar, bf_scroll_scale, eminvar, emaxxanesvar, emaxexafsvar, constkUsevar, kstepvar, estepvar, sam_numer, sam_denom, sam_logvar, ref_numer, ref_denom, ref_logvar, FlatUsevar, updownvar, edgestepwidth
	num_spectra, data_file, folder, file_type, n_cpus, ProcessSamvar, ProcessRefvar, NormUsevar, InterpUsevar, AutoAlignUsevar, edgestepvar, Filtervar, bf_scroll_scale, eminvar, emaxxanesvar, emaxexafsvar, constkUsevar, kstepvar, estepvar, sam_numer, sam_denom, sam_logvar, ref_numer, ref_denom, ref_logvar, FlatUsevar, updownvar, edgestepwidth = options
		   
	spectra = np.asarray(range(int(num_spectra)))
		
	rootsfile = rootsfile_check(folder)
	calibrated = calibfile_check(folder)
	normalized = normfile_check(folder, NormUsevar)
	
	if rootsfile == True:
		global minr, min_ang, max_ang
		minr = roots_file_data[2]
		min_ang = roots_file_data[0]
		max_ang = roots_file_data[1]
	
	if ID == 0:
		j=0
		
	if float(Filtervar) == 1:
		Wn = float(bf_scroll_scale)
		N  = 3    # Filter order
		B, A = signal.butter(N, Wn, output='ba')
		Wn_interp = float(0.25*bf_scroll_scale)
		global B_interp, A_interp
		B_interp, A_interp = signal.butter(N, Wn_interp, output='ba')
		
	if float(InterpUsevar) == 1:\
		xnew = (pd.read_csv(folder+'/xnew.dat', sep='\t', header = None)).iloc[:,0].values
		
	time.sleep(0.5*ID/n_cpus)
	
	if (calibrated == True) & (rootsfile == True):
		while not q.empty():
			i = q.get()
			root_i = spectra[i]
		
			if '.bin' in file_type:
				ang,mu_sam,mu_ref,flag = data_read_bin(root_i, minr)
				
			elif '.qex' in file_type:
				ang,mu_sam,mu_ref,time_points,AdcClock,flag = data_read_qex(root_i, minr)
				time_points = (np.cumsum(time_points)-time_points[0])/AdcClock
				
			energy = 1239.852/(float(calibration_values[3])*np.sin((ang+(calibration_values[1]-calibration_values[0]))*np.pi/180))
			
			print('updownvar :', updownvar, flag)
			
			if flag == float(updownvar):
				if energy[0] > energy[1]:
					energy = energy[::-1]
					mu_sam = mu_sam[::-1]
					mu_ref = mu_ref[::-1]
			
				if float(Filtervar) == 1:
					if float(ProcessSamvar) == 1:
						mu_sam = signal.filtfilt(B,A,mu_sam)
					if float(ProcessRefvar) == 1:
						mu_ref = signal.filtfilt(B,A,mu_ref)
					
				if (normalized == True) & (float(NormUsevar) == 1):
					if float(ProcessSamvar) == 1:
						mu_sam = normalize_data(energy, mu_sam)
					if float(ProcessRefvar) == 1:
						mu_ref = normalize_data(energy, mu_ref)
						
				if (float(AutoAlignUsevar) == 1) & (float(NormUsevar) == 1) & (normalized == True) & (float(ProcessRefvar) == 1):
					popt, pcov = curve_fit(cumgauss, energy, mu_ref, (float(e0), float(edgestepwidth)))
					ec = float(edgestepvar) - popt[0]	
					energy = energy + ec
						
				if float(InterpUsevar) == 1:
					energy_input = energy
					if float(ProcessSamvar) == 1:
						energy, mu_sam = interpolate(energy_input, mu_sam, xnew)
					if float(ProcessRefvar) == 1:
						energy, mu_ref = interpolate(energy_input, mu_ref, xnew)
					
				df = pd.DataFrame()
				df['E'] = energy
				if float(ProcessSamvar) == 1:
					df['mu_sam'] = mu_sam
				if float(ProcessRefvar) == 1:
					df['mu_ref'] = mu_ref
				if 'qex' in file_type:
					df['time'] = time_points
				
				if updownvar == 0:
					df.to_csv(folder+'/Export/Individual_Up/'+str(i)+'.dat', sep='\t', header=True, index=False)
				else:
					df.to_csv(folder+'/Export/Individual_Down/'+str(i)+'.dat', sep='\t', header=True, index=False)

			qs = q.qsize()
			
			print('Percent Complete:',str(100*(int(options[0])-qs)/int(options[0])))
			if ID == 0:
				tps = (time.time()-etd)/(int(options[0])-qs)
				print('Time Per Spectrum (s):',round(tps, 4))
				print('Estimated Time Remaining:',time.strftime("%H:%M:%S", time.gmtime(qs*tps)))
			sys.stdout.flush()
	
################################################################################  
def data_read_qex(root_i, minr):
	headerSize, line_bytes, dt, nData, nChannels, AdcClock = header_read_qex(data_file)

	minr = int(minr)
	root = int(roots_file_data[root_i+3])
	
	qex_file = codecs.open(data_file+'.qex', 'rb', encoding='cp1252') 

	qex_file.seek(int(headerSize)+(int(line_bytes)*root))

	data = np.fromfile(qex_file, dtype=dt, count = minr)
	
	ang = data['encoder']
	
	if ang[0] > ang[-1]:
		flag = 0
	else:
		flag = 1
	
	mu_sam = np.zeros((minr, 1))
	if float(ProcessSamvar) == 1:
		try:
			col_num = float(sam_denom.split(' ')[1])
		except:
			try:
				col_num = nChannels
			except:
				print('failed to extract column number')
		if col_num > (nChannels-1):
			if float(sam_logvar) == 1:
				mu_sam = np.log(data[sam_numer])
			else:
				mu_sam = data[sam_numer]
		else:
			if float(sam_logvar) == 1:
				mu_sam = np.log(data[sam_numer]/data[sam_denom])
			else:
				mu_sam = data[sam_numer]/data[sam_denom]
	
	mu_ref = np.zeros((minr, 1))	
	if float(ProcessRefvar) == 1:
		try:
			col_num = float(ref_denom.split(' ')[1])
		except:
			try:
				col_num = nChannels
			except:
				print('failed to extract column number')
	
		if col_num > (nChannels-1):
			if float(ref_logvar) == 1:
				mu_ref = np.log(data[ref_numer])
			else:
				mu_ref = data[ref_numer]
		else:
			if float(ref_logvar) == 1:
				mu_ref = np.log(data[ref_numer]/data[ref_denom])
			else:
				mu_ref = data[ref_numer]/data[ref_denom]
		
	RawData = pd.DataFrame()
	RawData['ang'] = np.around(ang, decimals=5)
	RawData['mu_sam'] = mu_sam
	RawData['mu_ref'] = mu_ref
	RawData['time'] = data['time']
	
	RawData.dropna(how='any', inplace = True)
	
	RawData = RawData.groupby('ang', as_index=False).mean()
	RawData = RawData[(RawData['ang'] >= np.float(min_ang)) & (RawData['ang'] <= np.float(max_ang))]
	
	RawData.sort_values(by='ang', ascending=False)
	RawData.reset_index(drop=True)
	
	return np.around(RawData['ang'].values, decimals=5), RawData['mu_sam'].values, RawData['mu_ref'].values, RawData['time'].values, AdcClock, flag
	
################################################################################   
def header_read_qex(filename):
	header_lines = []
	with codecs.open(filename+'.qex', 'rb', encoding='cp1252', errors='ignore') as qex_file:
		line = qex_file.readline().strip('\r\n').replace('# ', '')
		while not '_Header_End_'in line:
			header_lines.append(line)
			line = qex_file.readline().strip('\r\n').replace('# ', '')
		
	#search headerlines for keywords
	headerSize = np.int([x for x in header_lines if 'FileHeaderSize_byte' in x][0].split(': ')[1])
	nColumns = np.int([x for x in header_lines if 'AdcNumberColumnsInDataFile' in x][0].split(': ')[1])
	nChannels = np.int([x for x in header_lines if 'AdcNumberChannelsStored' in x][0].split(': ')[1])
	DataLineFormat = ([x for x in header_lines if 'DataLineFormat' in x][0].split(': ')[1]).split(', ')
	DataLineLabels = ([x for x in header_lines if 'DataLineLabels' in x][0].split(': ')[1]).split(', ')
	AdcClock = np.int([x for x in header_lines if 'AdcClock_Hz' in x][0].split(': ')[1])

	#Interpret DataTypes to extract number of bytes per line
	line_bytes = int(sum([int(format_item[1]) for format_item in [re.split(r'(\d+)', s) for s in DataLineFormat]])/8)
	
	#calculate the number of datapoints in file
	nLines = int((os.path.getsize(filename+'.qex') - headerSize)/(line_bytes))
	nData = int(nColumns*nLines)
	
	d_types = [('encoder', DataLineFormat[0]), ('time', DataLineFormat[1])]
	for i in range(nChannels):
		d_types.append(tuple(('CH '+str(i), DataLineFormat[i+2])))

	dt = np.dtype(d_types)	

	return headerSize, line_bytes, dt, nData, nChannels, AdcClock
		
################################################################################  
def data_read_bin(root_i, minr):
	headerSize, nData = header_read_bin(data_file)
	ch_headerSize, nChannels, d_types = header_read_bin_ch(data_file)
	
	root = roots_file_data[root_i+3]
	minr = int(minr)
	dt = np.dtype(d_types)
	
	g = codecs.open(data_file+'.bin', 'rb', encoding='cp1252')
	g.seek(int(ch_headerSize+(4*nChannels*root)))
	data = (np.fromfile(g, dtype=dt, count = int(minr)))
	
	encoder_bin = data_file+'_Encoder'
	f = codecs.open(encoder_bin+'.bin', 'rb', encoding='cp1252')
	f.seek(int(headerSize+(4*root)))
	ang = np.around(np.fromfile(f, dtype='f4', count = int(minr)), decimals=5)

	if ang[0] > ang[-1]:
		flag = 0
	else:
		flag = 1
	
	mu_sam = np.zeros((minr, 1))
	if float(ProcessSamvar) == 1:
		try:
			col_num = float(sam_denom.split(' ')[1])
		except:
			try:
				col_num = nChannels
			except:
				print('failed to extract column number')
		if col_num > (nChannels-1):
			if float(sam_logvar) == 1:
				mu_sam = np.log(data[sam_numer])
			else:
				mu_sam = data[sam_numer]
		else:
			if float(sam_logvar) == 1:
				mu_sam = np.log(data[sam_numer]/data[sam_denom])
			else:
				mu_sam = data[sam_numer]/data[sam_denom]
	
	mu_ref = np.zeros((minr, 1))	
	if float(ProcessRefvar) == 1:
		try:
			col_num = float(ref_denom.split(' ')[1])
		except:
			try:
				col_num = nChannels
			except:
				print('failed to extract column number')
		if col_num > (nChannels-1):
			if float(ref_logvar) == 1:
				mu_ref = np.log(data[ref_numer])
			else:
				mu_ref = data[ref_numer]
		else:
			if float(ref_logvar) == 1:
				mu_ref = np.log(data[ref_numer]/data[ref_denom])
			else:
				mu_ref = data[ref_numer]/data[ref_denom]
		
	RawData = pd.DataFrame()
	RawData['ang'] = np.around(ang, decimals=5)
	RawData['mu_sam'] = mu_sam
	RawData['mu_ref'] = mu_ref
	
	RawData.dropna(how='any', inplace = True)
	
	RawData = RawData.groupby('ang', as_index=False).mean()
	RawData = RawData[(RawData['ang'] >= np.float(min_ang)) & (RawData['ang'] <= np.float(max_ang))]
	
	RawData.sort_values(by='ang', ascending=False)
	RawData.reset_index(drop=True)
	
	return np.around(RawData['ang'].values, decimals=5), RawData['mu_sam'].values, RawData['mu_ref'].values, flag
			
################################################################################     
def header_read_bin(filename):
	encoder_bin = filename+'_Encoder'
	f = open(encoder_bin+'.bin', 'rb')
	f.seek(0)
	data = np.fromfile(f, dtype=np.int32, count = 2)
	headerSize = data[0]
	data = np.fromfile(f, dtype=np.int64, count = 1)
	data = np.fromfile(f, dtype='f4', count = 2)
	AdcClock_Hz = data[1]
	DacClock_Hz = 1
	nData = int((os.path.getsize(encoder_bin+'.bin') - headerSize)/4)
	
	return headerSize, nData

################################################################################ 
def header_read_bin_ch(filename):
	g = open(filename+'.bin', 'rb')
	data = np.fromfile(g, dtype=np.int32, count = 2)
	ch_headerSize = data[0]
	data = np.fromfile(g, dtype=np.int64, count = 1)
	data = np.fromfile(g, dtype='f4', count = 2)
	nChannels = int(data[0])
	d_types = []
	
	for i in range(nChannels):
		d_types.append(tuple(('CH '+str(i), 'f4')))
	
	return ch_headerSize, nChannels, d_types
	
################################################################################ 
def normalize_data(energy, mu):   	
	linepre_y = regression('pre-edge',np.asarray(energy), np.asarray(mu))        
	linepost_y = regression('post-edge',np.asarray(energy), np.asarray(mu))
	
	index = min(range(len(energy)), key=lambda i: abs(energy[i]-float(e0)))
	normdivisor = (linepost_y-linepre_y)[index]
	
	mu = (mu-linepre_y)/normdivisor
	linepost_y = (linepost_y-linepre_y)/normdivisor
	linepre_y = np.zeros(len(linepre_y))
					
	if float(FlatUsevar) == 1:
		flat_correction = 1-linepost_y
		if energy[index] > energy[-1]:
			flat_correction[index:int(len(energy))] = 0
		else:
			flat_correction[0:index] = 0
		mu = mu + flat_correction
		return mu
	else:
		return mu
			
################################################################################                 
def regression(region,x,y): 
	global pre1,pre2,post1,post2,order_pre,order_post,e0
	pre1,pre2,post1,post2,order_pre,order_post,e0 = normalisation_values
	
	if region == 'pre-edge':
		xllim = float(e0) + float(pre1)
		xhlim = float(e0) + float(pre2)
		try:
			order = str(order_pre[0])
		except:
			order = str(order_pre)

	if region == 'post-edge':
		xllim = float(e0) + float(post1)
		xhlim = float(e0) + float(post2)
		try:
			order = str(order_post[0])
		except:
			order = str(order_post)
			
	xregion = x[(x >= xllim) & (x <= xhlim)]
	yregion = y[(x >= xllim) & (x <= xhlim)]
	
	def constant(x, a):
		return a*np.ones(len(x))

	def linear(x, a, b):
		return (a*x) + b
	
	def quadratic(x, a, b, c):
		return (a*x) + (b*x**2) + c
	
	def cubic(x, a, b, c, d):
		return (a*x) + (b*x**2) + (c*x**3) + d
	
	def victoreen(x, a, b, c):
		f = 1.23986*10**4
		return ((a*f**3)/(x**3)) - ((b*f**4)/(x**4)) + c
	
	if order == '0':
		popt, pcov = curve_fit(constant, xregion, yregion)
		yfit = constant(x, *popt)

	if order == '1':
		popt, pcov = curve_fit(linear, xregion, yregion)
		yfit = linear(x, *popt)
		
	if order == '2':
		popt, pcov = curve_fit(quadratic, xregion, yregion)
		yfit = quadratic(x, *popt)
		
	if order == '3':
		popt, pcov = curve_fit(cubic, xregion, yregion)
		yfit = cubic(x, *popt)
		
	if order == 'V':
		popt, pcov = curve_fit(victoreen, xregion, yregion)
		yfit = victoreen(x, *popt)

	return yfit
	
################################################################################  	
def interpolate(x,y,xnew):
	windowloc=200
	x = np.asarray(x)
	y = np.asarray(y)
	
	if x[0] > x[-1]:
		x=x[::-1]
		y=y[::-1]
	
	ynew = np.zeros(len(xnew))
	
	#Localised Radial Basis functions by looping through subset of data
	localisation = int(np.ceil(len(x)/windowloc))
	factor = np.ceil(len(x)/localisation)
	for idloc in range(localisation):
		#localisation allows for some overlapping to improve the consistency upon reforming total data
		#determines if the sub-data is the last sub-data group
		if idloc == localisation-1:
			start = int((idloc*factor)-5)
			end = int(len(x))
			#determines if the sub-data is the first sub-data group
		elif idloc == 0:
			start = 0
			end = (int(1*factor)+5)
		else:
			start = int((idloc*factor)-5)
			end = int(((idloc+1)*factor)+5)
	
		ydata_loop = y[start:end]
		xdata_loop = x[start:end]
	
		#selects the region of x data to extract the new y data
		xnew_cut_idx = npi.indices(xnew, xnew[(xnew >= xdata_loop[0]) & (xnew < xdata_loop[-1])])
		xnew_cut=xnew[xnew_cut_idx]

		#caries out the radial basis function interpolation on the localised sub-data
		rbf = Rbf(xdata_loop, ydata_loop, function='linear', smooth=0)
		ynew[xnew_cut_idx] = rbf(xnew_cut)
	
	return [xnew, ynew]	
	
################################################################################  	
def cumgauss(x, mu, sigma):
    return 0.5 * (1 + special.erf((x-mu)/(np.sqrt(2)*sigma)))
	
################################################################################  
if __name__ == '__main__':
	pyfile,options = sys.argv
	spawn(options)