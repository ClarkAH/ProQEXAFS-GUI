import warnings
warnings.filterwarnings("ignore",category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
import psutil, os, peakutils, time, sys, ast, codecs, re, itertools
import numpy as np     
import multiprocessing as mp
from scipy.signal import savgol_filter	
from scipy.signal import find_peaks
import matplotlib.pyplot as plt	
		
def spawn(options):
	options = ast.literal_eval(options)
	procs = list()
	n_cpus = int(options[0])
	buffer_points = np.asarray(range(int(options[1])+1)) 
	
	job_indexes = jobDiv(buffer_points, n_cpus)
	lock = mp.Lock()
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
	for cpu_i in range(n_cpus):
		d = dict(affinity = int(core_affinity_list[cpu_i].item()), work = job_indexes[cpu_i], ID = cpu_i, options = options, lock=lock)
		p = mp.Process(target=wrapper_targetFunc, kwargs=d)
		p.start()
		procs.append(p)
	for p in procs:
		p.join()
	print('Returned to Main')
	sys.stdout.flush()
	sys.exit(0)		

def wrapper_targetFunc(affinity, work, ID, options, lock):
	waiting = False
	proc = psutil.Process()  # get self pid
	proc.cpu_affinity([affinity])
	print('Process-'+str(ID)+' queueing ',work)
	sys.stdout.flush()
	etd = time.time()
	targetFunc(work, etd, ID, options, lock)
	lock.acquire()
	print('Process-'+str(ID)+' Task Complete')
	if ID == 0:
		print('Percent Complete =',100)
		print('Returning to Main')
		waiting = True
	if waiting == True:
		print('Waiting')
	sys.stdout.flush()
	lock.release()
	os._exit(0)
	
def jobDiv(inputArray, numCPUs):
	jobs = []
	arrayLength = len(inputArray)
	jobRange = int(arrayLength / numCPUs)
	extra = arrayLength - (jobRange * numCPUs)
	prevEnd = 0
	for c in range(numCPUs):
		endIdx = (c * jobRange) + jobRange - 1
		if (c == (numCPUs-1)):
			endIdx += extra
		startIdx = prevEnd
		if ( (c > 0) and (startIdx+1 < arrayLength) ):
			startIdx += 1
		jobs.append( (startIdx, endIdx) )
		prevEnd = endIdx
	return jobs
	
def targetFunc(work, etd, ID, options, lock):
	ncpus, buffer_points, headerSize, buffer, encoder_bin, nData, resample_factor, enconder_sampling, encoder_smooth, file_type = options
	roots = []
	angs = []
	
	time.sleep(ID/ncpus)   
	
	if 'bin' in file_type:
		f = open(encoder_bin+'.bin', 'rb')
		f.seek(0)
	
	if 'qex' in file_type:
		headerSize, line_bytes, dt, nData, nLines = header_read_qex(encoder_bin+'.qex')
		f = open(encoder_bin+'.qex', 'rb')
		#encoder = np.zeros(nLines)
		encoder_2 = np.zeros(np.int(np.floor(nLines/resample_factor)))
		#time_data = np.zeros(nLines)
		f.seek(0)
	
	for i in range(work[0], work[1]+1):
		start_time = etd
	
		if 'bin' in file_type:
			f.seek(int(headerSize+(i*buffer)))
			if ((i+1)*buffer) > nData*4:
				data_E = np.fromfile(f, dtype='f4', count = -1)
			else:
				data_E = np.fromfile(f, dtype='f4', count = int(buffer/4))
			y = data_E[0::resample_factor]
	
		if 'qex' in file_type:
			f.seek(int(headerSize))
			
			for k in range(int(np.ceil(nLines/(resample_factor*80000)))):
				if ((k+1)*(resample_factor*80000)) > nLines:
					data_E = np.fromfile(f, dtype=dt, count = -1)
					#encoder[k*8000000:nLines+1] = data_E['encoder']
					encoder_2[np.int(k*(resample_factor*80000)/resample_factor):np.int(np.floor(nLines/resample_factor)+1)] = data_E['encoder'][0::resample_factor]
					#time_data[k*8000000:nLines+1] = data_E['time']
				else:
					data_E = np.fromfile(f, dtype=dt, count = int(resample_factor*80000))
					#encoder[k*8000000:(k+1)*8000000] = data_E['encoder']
					encoder_2[np.int(k*(resample_factor*80000)/resample_factor):np.int((k+1)*(resample_factor*80000)/resample_factor)] = data_E['encoder'][0::resample_factor]
					#time_data[k*8000000:(k+1)*8000000] = data_E['time']
					
				lock.acquire()
				print('data read :', (100*(k+1))/int(np.ceil(nLines/(resample_factor*80000))), ' %')
				sys.stdout.flush()
				lock.release()

				
		#takes derivate of encoder data and applies Savitzky-Golay filter
		if 'bin' in file_type:
			try:
				c1_derivative = savgol_filter(np.gradient(savgol_filter(y, int(encoder_smooth), 2)), int(encoder_smooth), 2)
			except:
				c1_derivative = np.gradient(y)
				
		if 'qex' in file_type:
			##find turning points of encoder by peaks in number of averaged data points per encoder value
			##the minimum distance between peaks is defined by considering the number of encoder values in an oscillation (actually allowing for some 'short' oscillations)
			
			try:
				#smoothed_encoder = savgol_filter(encoder[0::resample_factor], int(encoder_smooth), 2)
				smoothed_encoder = savgol_filter(encoder_2, int(encoder_smooth), 2)
				c1_derivative = savgol_filter(np.gradient(smoothed_encoder), int(encoder_smooth), 2)
			except:
				#c1_derivative = np.gradient(encoder[0::resample_factor])\
				c1_derivative = np.gradient(encoder_2)
				
			print('print data :', 'derivative caluculated')
				
			indicies = np.asarray(np.where(np.diff(np.sign(c1_derivative)))[0]+1)
			
			#fig, axs = plt.subplots(2, 1)
			
			#axs[0].plot(smoothed_encoder)
			#axs[0].plot(indicies, smoothed_encoder[indicies], 'kx')
			#axs[0].set_xlabel('Encoder Datapoint')
			#axs[0].set_ylabel('Monochromator Angle')
			#
			#axs[1].plot(c1_derivative)
			#axs[1].plot(indicies, c1_derivative[indicies], 'kx')
			#axs[1].set_xlabel('Encoder Datapoint')
			#axs[1].set_ylabel('Derivative Monochromator Angle')
			#plt.show()
			#plt.close()
			
			print('print data :', 'indicies found')

			#roots_ang = list(encoder[0::resample_factor][indicies])
			roots_ang = list(encoder_2[indicies])
			
			#if (i == 0) & (ID == 0):
			#	min_diff = np.min(np.diff(indicies))
			#	#print('print data :', min_diff)
			#	#print('print data :', peaks)
			#	indicies2 = indicies[:-1] + min_diff
			#	#print('print data :', peaks2)
			#	sys.stdout.flush()
			#	
			#	roots_ang2 = list(encoder[0::resample_factor][indicies2])
			#	
			#	roots_check = list(itertools.chain(roots_ang, roots_ang2))
			#	
			#	#print('print data :', roots_check)
			#	sys.stdout.flush()
			#	
			#	B = [i for i in roots_check  if i >= np.average(roots_check)]
			#	C = [i for i in roots_check if i <= np.average(roots_check)]
			#	
			#	max_ang = np.min(B)
			#	min_ang = np.max(C)
			#    
			#	lock.acquire()
			#	print('min_ang :', min_ang)
			#	print('max_ang :', max_ang)
			#	sys.stdout.flush()
			#	lock.release()

		#finds the roots of the encoder file used for splitting the raw absorption data into spectra
		if 'bin' in file_type:
			indicies = np.asarray((np.where(np.diff(np.sign(c1_derivative)))[0]+1)*resample_factor)
			roots_ang = list(data_E[indicies])
		#if 'qex' in file_type:	
		#	indicies = np.asarray((np.where(np.diff(np.sign(c1_derivative)))[0]+1))
		#	roots_ang = list(y[indicies])
		
		if 'bin' in file_type:
			roots_loop = indicies+(i*buffer/4)
		if 'qex' in file_type:
			roots_loop = indicies*resample_factor
		
		roots = np.concatenate((roots, roots_loop), axis=0)
		angs = np.concatenate((angs, roots_ang), axis=0)
		sys.stdout.flush()
		
		if ID == 0:
			tpl = time.time()-start_time
			lock.acquire()
			print('Percent Complete:',str(100*i/(work[1]-work[0]+1)))
			print('Time Per Loop (s):',tpl)
			print('Estimated Time Remaining:',time.strftime("%H:%M:%S", time.gmtime(int(tpl*(work[1]-i-1)))))
			sys.stdout.flush()
			lock.release()

	lock.acquire()			
	print('angs_'+str(ID)+'_'+str(list(angs)))
	sys.stdout.flush()
	print('roots_'+str(ID)+'_'+str(list(roots)))
	sys.stdout.flush()
	lock.release()
	
def header_read_qex(filename):
	header_lines = []
	with codecs.open(filename, 'rb', encoding='cp1252', errors='ignore') as qex_file:
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
	
	#Interpret DataTypes to extract number of bytes per line
	line_bytes = int(sum([int(format_item[1]) for format_item in [re.split(r'(\d+)', s) for s in DataLineFormat]])/8)
	
	#calculate the number of datapoints in file
	nLines = int((os.path.getsize(filename) - headerSize)/(line_bytes))
	nData = int(nColumns*nLines)
	
	d_types = [('encoder', DataLineFormat[0]), ('time', DataLineFormat[1])]
	for i in range(nChannels):
		d_types.append(tuple(('CH '+str(i), DataLineFormat[i+2])))
	
	dt = np.dtype(d_types)	
	
	return headerSize, line_bytes, dt, nData, nLines
			
if __name__ == '__main__':
	pyfile,options = sys.argv
	spawn(options)
