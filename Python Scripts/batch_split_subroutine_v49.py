import psutil, os, peakutils, time, sys, ast
import numpy as np     
import multiprocessing as mp
from scipy.signal import savgol_filter		
		
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
    ncpus, buffer_points, headerSize, buffer, encoder_bin, nData, resample_factor, enconder_sampling, encoder_smooth = options
    roots = []
    angs = []
	
    time.sleep(ID/ncpus)    
    f = open(encoder_bin+'.bin', 'rb')
    f.seek(0)
    for i in range(work[0], work[1]+1):
        start_time = etd
        #print('loop iteration:',int(i))
        #print('file seek:',int(headerSize+(i*buffer)))
        f.seek(int(headerSize+(i*buffer)))
        if ((i+1)*buffer) > nData*4:
            data_E = np.fromfile(f, dtype='f4', count = -1)
        else:
            data_E = np.fromfile(f, dtype='f4', count = int(buffer/4))

        disk_size = data_E.nbytes
        #print('disk size:', disk_size)

        #x = np.asarray(range(len(data_E)))
        #x = x[0::resample_factor]
        y = data_E[0::resample_factor]
	
        #highs = peakutils.peak.indexes(y, min_dist=int(0.3*(enconder_sampling*1E6/float(mono_freq))/resample_factor))
        #lows = peakutils.peak.indexes(-y, min_dist=int(0.3*(enconder_sampling*1E6/float(mono_freq))/resample_factor))
		
		#takes derivate of encoder data and applies Savitzky-Golay filter
        c1_derivative = savgol_filter(np.gradient(savgol_filter(y, int(encoder_smooth), 2)), int(encoder_smooth), 2)
		
		#finds the roots of the encoder file used for splitting the raw absorption data into spectra
        indicies = np.asarray((np.where(np.diff(np.sign(c1_derivative)))[0]+1)*resample_factor)
        #indicies = resample_factor*(np.asarray(sorted(np.concatenate((highs, lows), axis=0))))
        roots_ang = list(data_E[indicies])
        roots_loop = indicies+(i*buffer/4)
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
			
if __name__ == '__main__':
    pyfile,options = sys.argv
    spawn(options)
