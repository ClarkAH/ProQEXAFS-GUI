import warnings
warnings.filterwarnings("ignore",category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
import psutil, os, peakutils, time, sys, ast, codecs, re
import numpy as np     
import multiprocessing as mp
from scipy.signal import savgol_filter	
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
            if ((i+1)*buffer) > nLines*line_bytes:
                data_E = np.fromfile(f, dtype=dt, count = -1)
            else:
                data_E = np.fromfile(f, dtype=dt, count = int(buffer))
            y = data_E['encoder']

        disk_size = data_E.nbytes
			
		#takes derivate of encoder data and applies Savitzky-Golay filter
        c1_derivative = savgol_filter(np.gradient(savgol_filter(y, int(encoder_smooth), 2)), int(encoder_smooth), 2)
		
		#finds the roots of the encoder file used for splitting the raw absorption data into spectra
        if 'bin' in file_type:
            indicies = np.asarray((np.where(np.diff(np.sign(c1_derivative)))[0]+1)*resample_factor)
        if 'qex' in file_type:	
            indicies = np.asarray((np.where(np.diff(np.sign(c1_derivative)))[0]+1))[0::2] 
        roots_ang = list(y[indicies])
        if 'bin' in file_type:
            roots_loop = indicies+(i*buffer/4)
        if 'qex' in file_type:
            roots_loop = indicies+(i*buffer)
		
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
