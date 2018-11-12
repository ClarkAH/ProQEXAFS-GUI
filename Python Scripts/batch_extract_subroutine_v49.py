import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)
import numpy as np
import pandas as pd
import time, psutil, ast, sys, os
from scipy.signal import savgol_filter
import scipy.signal as signal
from scipy.optimize import curve_fit
from scipy.interpolate import Rbf
import numpy_indexed as npi
#from lmfit.models import StepModel
from scipy import special
from scipy.optimize import curve_fit
import multiprocessing as mp

def spawn(options):
    options = ast.literal_eval(options)
    procs = list()
    n_cpus = int(options[2])
    #job_indexes = jobDiv(np.asarray(range(int(options[0]))), n_cpus)
    rootsfile = rootsfile_check(options[1])
    calibrated = calibfile_check(options[1])
    normalized = normfile_check(options[1])
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
	
    if (rootsfile == True) and (calibrated == True) and (normalized == True):
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
        
def rootsfile_check(folder):
    try:
        global calibroots
        calibroots = (pd.read_csv(folder+'/roots.dat', sep='\t', header=None)).values
        return True
    except:
        print('Encoder Analysis Required')
        return False

def calibfile_check(folder):
    try:
        global calibration_values
        calibration_values = (pd.read_csv(folder+'/calibration.dat', sep='\t', header=None)).values
        return True
    except:
        print('Cannot Find Calibration Parameter File')
        return False
        
def normfile_check(folder):
    try:
        global normalisation_values
        normalisation_values = (pd.read_csv(folder+'/Normalisation.dat', sep='\t', header=None)).values
        return True
    except:
        print('Cannot Find Normalisation Parameter File')
        return False
    
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
    
def targetFunc(q, etd, ID, options):
    global core_hole, batch_edgestepvar, batch_toggle, xnew, job_indexes, rootsfile, calibrated, normalized, spectra, n_cpus, sam_process, ref_process, export_path, folder, data_bin, batch_sam_numer, batch_sam_denom, batch_sam_logvar, batch_ref_numer, batch_ref_denom ,batch_ref_logvar, batch_updownvar, batch_filtervar, batch_firstcalibvar, batch_bf_scroll_scale, batch_emaxxanesvar, Sav_Gol, windowev
    spectra, folder, n_cpus, sam_process, ref_process, export_path, data_bin, batch_sam_numer, batch_sam_denom, batch_sam_logvar, batch_ref_numer, batch_ref_denom ,batch_ref_logvar, batch_updownvar, batch_filtervar, batch_edgestepvar, batch_bf_scroll_scale, batch_emaxxanesvar, Sav_Gol, windowev, core_hole = options

    spectra=np.asarray(range(int(spectra)))
    xnew_df= pd.read_csv(folder+'/xnew_info.dat', sep='\t', header = None)
    xnew = xnew_df.iloc[:,0].values
    
    rootsfile = rootsfile_check(folder)
    calibrated = calibfile_check(folder)
    normalized = normfile_check(folder)
    if ID == 0:
        j=0
    time.sleep(0.5*ID/n_cpus)    
    while not q.empty():
        i = q.get()
        df = pd.DataFrame()
        job = spectra[i]
        x1,y1,x2,y2,x3,y3 = batch_extract(job)
        if (len(x1) == 1) and (len(x2) == 1):
            skip = True
        elif (len(x1) == 2) and (len(x2) == 2):
            skip = True
        else:
            skip = False
            
        if skip == False:
            df['E'] = x1
            if sam_process == 1:
                df['mu_sam'] = y1 
            if ref_process == 1:
                df['mu_ref'] = y2
            df.to_csv(export_path+'/Individual/'+str(i)+'.dat', sep='\t', header=True, index=False)
            if ID == 0:
                j=j+1
                tps = (time.time()-etd)/(j)
                qs = q.qsize()
                print('Percent Complete:',str(100*(int(options[0])-qs)/int(options[0])))
                print('Time Per Spectrum (s):',round(tps/(2*n_cpus), 4))
                print('Estimated Time Remaining:',time.strftime("%H:%M:%S", time.gmtime(qs*tps/(2*n_cpus))))
                sys.stdout.flush()
                
def batch_extract(i):
    try:
        if (rootsfile == True) and (calibrated == True) and (normalized == True):
            global batch_toggle
            #i = int(ds_scroll_scale.get())
            if sam_process == 1:
                batch_toggle = 0
                ang_sam, mu_sam = spectrumread(calibroots[i+3],calibroots[2])
                energy_sam = np.transpose(1239.852/(calibration_values[3]*np.sin((ang_sam+(calibration_values[1]-calibration_values[0]))*np.pi/180))).tolist()
                if energy_sam[-1] < energy_sam[0]:
                    energy_sam = energy_sam[::-1]
                    mu_sam = mu_sam[::-1]

            if ref_process == 1:
                batch_toggle = 1
                ang_ref, mu_ref = spectrumread(calibroots[i+3],calibroots[2])
                energy_ref = np.transpose(1239.852/(calibration_values[3]*np.sin((ang_ref+(calibration_values[1]-calibration_values[0]))*np.pi/180))).tolist()
                if energy_ref[-1] < energy_ref[0]:
                    mu_ref = mu_ref[::-1]
                    energy_ref = energy_ref[::-1]
            if angl < angf:
                direction = 0
            else:
                direction = 1
            
            if direction == batch_updownvar:
                if batch_filtervar == 1:
                    Wn = float(batch_bf_scroll_scale) # Cutoff frequency
                    N  = 3    # Filter order
                    B, A = signal.butter(N, Wn, output='ba')
                    if sam_process == 1:
                        mu_sam_filtered = signal.filtfilt(B,A,mu_sam)
                    if ref_process == 1:
                        mu_ref_filtered = signal.filtfilt(B,A,mu_ref)
                else:
                    if sam_process == 1:
                        mu_sam_filtered = mu_sam
                    if ref_process == 1:
                        mu_ref_filtered = mu_ref

                pre1,pre2,post1,post2,order_pre,order_post,e0 = normalisation_values
                e0 = float(e0[0])

                if Sav_Gol == True:
                    #print('Savitzky-Golay window in eV ='+str(windowev))
                    if sam_process == 1:
                        idxWin = npi.indices(np.asarray(energy_sam), np.asarray(energy_sam)[(np.asarray(energy_sam) >= e0) & (np.asarray(energy_sam) <= (e0 + windowev))])
                    else:
                        idxWin = npi.indices(np.asarray(energy_ref), np.asarray(energy_ref)[(np.asarray(energy_ref) >= e0) & (np.asarray(energy_ref) <= (e0 + windowev))])
                    global window
                    window = len(idxWin)
                    #print('Savitzky-Golay window in datapoints ='+str(window))
                    #makes sure the window is odd in length
                    if (window % 2) == 0:
                        window = window - 1
                    if window < 3:
                        window = 3
                    global poly
                    poly = 2

                if normalized == True:
                    if sam_process == 1:
                        mu_sam = normalize(energy_sam, mu_sam)
                        mu_sam_norm = normalize(energy_sam, mu_sam_filtered)
                    if ref_process == 1:
                        mu_ref_norm = normalize(energy_ref, mu_ref_filtered)
                
                        #fits step function to the reference edge step to allow for recalibration of all spectra
                        #step_mod = StepModel(form='erf', prefix='step_')
                        #pars = step_mod.make_params(sigma=2, amplitude=1, center=float(e0))
                        #mod = step_mod
                        #out = mod.fit(mu_ref_norm, pars, x=energy_ref)
						
                        popt, pcov = curve_fit(cumgauss, energy_ref, mu_ref_norm, (e0, core_hole))
							
						
                if (i == 0) or (i==1):
                    if ref_process == 1:
                            #ec = batch_edgestepvar - out.best_values['step_center']
                            ec = batch_edgestepvar - popt[0]
                    else:
                        ec = 0
                else:
                    if ref_process == 1:
                        #ec = batch_edgestepvar - out.best_values['step_center']
                        #print('step postion = ',out.best_values['step_center'])
                        ec = batch_edgestepvar - popt[0]
                    else:
                        ec = 0

                if sam_process == 1:
                    energy_sam = np.asarray(energy_sam) + ec
                    xinterp_sam, yinterp_sam = interpolate(energy_sam, mu_sam_norm, xnew)
                if ref_process == 1:
                    energy_ref = np.asarray(energy_sam)
                    xinterp_ref, yinterp_ref = interpolate(energy_ref, mu_ref_norm, xnew)

                if np.min(yinterp_sam) < -0.4:
                    return([0],[0],[0],[0],[0],[0])						
                elif (sam_process == 1) and (ref_process == 1):
                    return(xinterp_sam, yinterp_sam, xinterp_ref, yinterp_ref, energy_sam, mu_sam)
                elif (sam_process == 0) and (ref_process == 1):
                    return([0], [0], xinterp_ref, yinterp_ref, [0], [0])
                elif (sam_process == 1) and (ref_process == 0):
                    return(xinterp_sam, yinterp_sam, [0], [0], energy_sam, mu_sam)
            else:
                return([0],[0],[0],[0],[0],[0])
    except:
        print('failed spectrum ', i)
        return([0,0],[0,0],[0,0],[0,0],[0,0],[0,0])
        
def spectrumread(rooti,minr):
    encoder_bin = data_bin+'_Encoder'
    
    minr = int(minr)
        
    g = open(data_bin+'.bin', 'rb')
    
    g.seek(0)

    data = np.fromfile(g, dtype=np.int32, count = 2)

    headerSize = data[0]

    data = np.fromfile(g, dtype=np.int64, count = 1)

    data = np.fromfile(g, dtype='f4', count = 2)

    nChannels = data[0]
    #AdcClock_Hz = data[1]
    #DacClock = 1
    
    choices = []
    for i in range(int(nChannels)):
        choices.append('CH '+str(int(i)))
    choices.append('1')

    g.seek(int(headerSize+(4*nChannels*rooti)))
    data_D = np.fromfile(g, dtype='f4', count = int(nChannels*minr))
    mu_r = np.zeros((minr, 1))
    
    if batch_toggle == 0:
        numervar = batch_sam_numer
        denomvar = batch_sam_denom    
        logvar = batch_sam_logvar
    else:
        numervar = batch_ref_numer
        denomvar = batch_ref_denom
        logvar = batch_ref_logvar
    
    for j in range(int(nChannels)+1):
        if numervar == choices[j]:
            numercol = int(j)
        if denomvar == choices[j]:
            denomcol = int(j)

    if logvar == 1:
        if denomcol == int(len(choices)-1):
            mu_r[:,0] = np.log(data_D[numercol::int(nChannels)])
        else:
            mu_r[:,0] = np.log(data_D[numercol::int(nChannels)]/data_D[denomcol::int(nChannels)])
    else:
        if denomcol == int(len(choices)-1):
            mu_r[:,0] = data_D[numercol::int(nChannels)]
        else:
            mu_r[:,0] = data_D[numercol::int(nChannels)]/data_D[denomcol::int(nChannels)]

    nans = np.argwhere(np.isnan(mu_r[:,0]))
    mu = np.zeros((minr-len(nans), 1))
    mu[:,0] = np.delete(mu_r[:,0], nans)

    f = open(encoder_bin+'.bin', 'rb')
    f.seek(int(headerSize+(4*rooti)))
    ang_r = np.fromfile(f, dtype='f4', count = int(minr))
    ang = np.delete(ang_r, nans)

    RawData = pd.DataFrame()
    RawData['ang'] = ang
    RawData['mu'] = mu
    
    global angf
    angf = ang[0]
    global angl
    angl = ang[-1]
    
    RawData=RawData[(RawData['ang'] >= calibroots[0][0]) & (RawData['ang'] <= calibroots[1][0])]  
    RawData = RawData.groupby('ang', as_index=False).mean()
    
    return RawData['ang'].values, RawData['mu'].values

def normalize(energy, mu):
    pre1,pre2,post1,post2,order_pre,order_post,e0 = normalisation_values
    e0 = float(e0[0])
    Flatvar = 1
    lineprey = regression('pre-edge',np.asarray(energy), np.asarray(mu), e0)     
    lineposty = regression('post-edge',np.asarray(energy), np.asarray(mu), e0)
    index = min(range(len(energy)), key=lambda i: abs(energy[i]-e0))
    normdivisor = (lineposty-lineprey)[index]
    mu = (mu-lineprey)/normdivisor
    lineposty = (lineposty-lineprey)/normdivisor
    lineprey = np.zeros(len(lineprey))
                    
    if Flatvar == 1:
        flat_correction = 1-lineposty
        if energy[index] > energy[-1]:
            flat_correction[index:int(len(energy))] = 0
        else:
            flat_correction[0:index] = 0
        mu = mu + flat_correction
        return mu
    else:
        return mu
    
def regression(region,x,y,calib_energy): 
    pre1,pre2,post1,post2,order_pre,order_post,e0 = normalisation_values
    if region == 'pre-edge':
        xllim = calib_energy + float(pre1[0])
        xhlim = calib_energy + float(pre2[0])
        try:
            order = str(int(order_pre[0]))
        except:
            order = str(order_pre[0])
            
    if region == 'post-edge':
        xllim = calib_energy + float(post1[0])
        xhlim = calib_energy + float(post2[0])
        try:
            order = str(int(order_post[0]))
        except:
            order = str(order_post[0])
            
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
	
def cumgauss(x, mu, sigma):
    return 0.5 * (1 + special.erf((x-mu)/(np.sqrt(2)*sigma)))
 
def interpolate(x,y,xnew):
    windowloc=200
    x = np.asarray(x)
    y = np.asarray(y)

    if Sav_Gol == True:
        idxSG = npi.indices(x, x[x <= batch_emaxxanesvar])
        y_xanes = y[0:np.max(idxSG)]
        y_exafs = y[np.max(idxSG):,]
        y_xsg = savgol_filter(y_xanes, window, poly)
        if len(y_exafs) > window*3:
            y_esg = savgol_filter(y_exafs, 3*window, poly)
        elif len(y_exafs) > window:
            y_esg = savgol_filter(y_exafs, window, poly)
        else:
            y_esg - y_exafs

    y = np.concatenate((y_xsg, y_esg), axis=0)

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

if __name__ == '__main__':
    pyfile,options = sys.argv
    spawn(options)
