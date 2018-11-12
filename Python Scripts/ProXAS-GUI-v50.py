import warnings
warnings.filterwarnings("ignore",category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
import numpy as np
import pandas as pd
import matplotlib as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import tkinter
from tkinter.filedialog import askopenfile
from tkinter import ttk
import time, os, psutil, subprocess, sys, shutil, ast
from tkinter import END
from scipy.signal import savgol_filter
from matplotlib.widgets import Button
import scipy.signal as signal
from scipy.optimize import curve_fit
from xraydb import XrayDB
from scipy.interpolate import Rbf
import numpy_indexed as npi
import gc
#from lmfit.models import StepModel
from scipy import special
from scipy.optimize import curve_fit
import itertools
import __main__ as main

print('ProXAS GUI')
class App(tkinter.Tk):
    def __init__(self,*args,**kwargs):
       tkinter.Tk.__init__(self,*args,**kwargs)
       self.resizable(True, True)
       
       self.notebook = ttk.Notebook()
       self.add_tab()
       self.notebook.grid(row=0)
       self.configure(background='white')
  
    def add_tab(self):
        tab = DataRead(self.notebook)
        tab2 = DataProcess(self.notebook)
        tab3 = PostProcess(self.notebook)
        self.notebook.add(tab,text="DataRead")
        self.notebook.add(tab2,text="DataProcess")
        self.notebook.add(tab3,text="PostProcess")
  
################################################################################	
################################################################################	
################################################################################	
class DataRead(tkinter.Frame):
    def __init__(self,name,*args,**kwargs):
        self.frame_dataread = tkinter.Frame.__init__(self,*args,**kwargs)
################################################################################

        labelfile = tkinter.Label(self, text='Load Data',anchor='w')
        labelfile.grid(column=1, row=2, columnspan=1, rowspan=1, sticky='S')
		
        self.buttonplus= tkinter.Button(self, text='+', width=4, height=1,command=self.fileplus)
        self.buttonplus.grid(column=2, row=2, columnspan=1, sticky='S')

        self.buttonminus= tkinter.Button(self, text='-', width=4, height=1,command=self.fileminus)
        self.buttonminus.grid(column=3, row=2, columnspan=1, sticky='S') 

        self.listbox = tkinter.Listbox(self)
        self.listbox.grid(column=1, row=3, columnspan = 3, rowspan=10)
        self.listbox['width'] = 30
        self.listbox['height'] = 40
        self.listbox.bind('<<ListboxSelect>>', self.onselect)
        
        self.progress_bar = ttk.Progressbar(self, orient='horizontal', mode='determinate', maximum=100, style="red.Horizontal.TProgressbar")
        self.progress_bar.grid(column=1, row=99, columnspan=99)
        self.progress_bar['length'] = 1000
		
        self.labelText = tkinter.StringVar()
        self.labelText.set('Time estimate')
        self.labeltime = tkinter.Label(self, textvariable=self.labelText)
        self.labeltime.grid(column=4, row=98, columnspan=3, rowspan=1, sticky='W')
        
        self.files_long = []
        self.edge_info = []
        self.edge_type = []
        self.edge_jump = 0
        
################################################################################
        
        self.plotcolsFrame = tkinter.Frame(self)
        self.plotcolsFrame.grid(row=3, column=7, rowspan=4, columnspan=3)
        
        plotcolslabel = tkinter.Label(self.plotcolsFrame, text='Plotting Columns',anchor='w')
        plotcolslabel.grid(column=0, row=0, columnspan=3, rowspan=1, sticky='N')
        
        numer = tkinter.Label(self.plotcolsFrame, text='Numerator',anchor='w')
        numer.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='W')
        
        denom = tkinter.Label(self.plotcolsFrame, text='Denomonator',anchor='w')
        denom.grid(column=0, row=2, columnspan=1, rowspan=1, sticky='W')
        
        self.sam_numer = tkinter.StringVar(self.plotcolsFrame)
        self.sam_numer.set('CH 0')
        self.sam_numerMenu = tkinter.OptionMenu(self.plotcolsFrame, self.sam_numer, 'CH 0','CH 1','CH 2','CH 3','1')
        self.sam_numerMenu.grid(column=1, row=1, columnspan=1) 
        
        self.sam_denom = tkinter.StringVar(self.plotcolsFrame)
        self.sam_denom.set('CH 1')
        self.sam_denomMenu = tkinter.OptionMenu(self.plotcolsFrame, self.sam_denom, 'CH 0','CH 1','CH 2','CH 3','1')
        self.sam_denomMenu.grid(column=1, row=2, columnspan=1) 
        
        self.ref_numer = tkinter.StringVar(self.plotcolsFrame)
        self.ref_numer.set('CH 1')
        self.ref_numerMenu = tkinter.OptionMenu(self.plotcolsFrame, self.ref_numer, 'CH 0','CH 1','CH 2','CH 3','1')
        self.ref_numerMenu.grid(column=2, row=1, columnspan=1) 
        
        self.ref_denom = tkinter.StringVar(self.plotcolsFrame)
        self.ref_denom.set('CH 2')
        self.ref_denomMenu = tkinter.OptionMenu(self.plotcolsFrame, self.ref_denom, 'CH 0','CH 1','CH 2','CH 3','1')
        self.ref_denomMenu.grid(column=2, row=2, columnspan=1) 
        
        log = tkinter.Label(self.plotcolsFrame, text='Log',anchor='w')
        log.grid(column=0, row=3, columnspan=1, rowspan=1, sticky='W') 
        
        self.sam_logvar = tkinter.DoubleVar()
        self.sam_logvar.set(1)
        self.sam_logCheck = tkinter.Checkbutton(self.plotcolsFrame, variable=self.sam_logvar)
        self.sam_logCheck.grid(column=1, row=3)
        
        self.ref_logvar = tkinter.DoubleVar()
        self.ref_logvar.set(1)
        self.ref_logCheck = tkinter.Checkbutton(self.plotcolsFrame, variable=self.ref_logvar)
        self.ref_logCheck.grid(column=2, row=3)
        
        filterlabel = tkinter.Label(self.plotcolsFrame, text='Filter',anchor='w')
        filterlabel.grid(column=0, row=5, columnspan=1, rowspan=1, sticky='W') 
        
        self.Filtervar = tkinter.DoubleVar()
        self.Filtervar.set(0)
        self.FilterCheck = tkinter.Checkbutton(self.plotcolsFrame, variable=self.Filtervar)
        self.FilterCheck.grid(column=2, row=5)
        
        plottoggle = tkinter.Label(self.plotcolsFrame, text='Plot',anchor='w')
        plottoggle.grid(column=0, row=4, columnspan=1, rowspan=1, sticky='W')
        
        toggle = [('',0),('',1)]
        self.radio_plot_buttons = []
        self.togglevar = tkinter.DoubleVar()
        self.togglevar.set(0)
        
        for idx, (text, mode) in enumerate(toggle):
            self.radio_plot_buttons.append(tkinter.Radiobutton(self.plotcolsFrame, text=text, variable=self.togglevar, value=mode))
            self.radio_plot_buttons[-1].grid(row=4, column=idx+1) 
        
################################################################################
        
        self.EncoderFrame = tkinter.Frame(self, width=300)
        self.EncoderFrame.grid(row=1,column=7,rowspan=3,sticky='EW')
				
        self.buttonsplit= tkinter.Button(self.EncoderFrame, text='Analyze Encoder', width=30, height=1,command=self.split)
        self.buttonsplit.grid(column=0, row=0, columnspan=3, rowspan=1)
		
        #self.buttonabort= tkinter.Button(self.EncoderFrame, text='Abort', width=30, height=1,command=self.abortsplit)
        #self.buttonabort.grid(column=0, row=4, columnspan=3, rowspan=1)
        
        refactor = tkinter.Label(self.EncoderFrame, text='Refactor',anchor='w')
        refactor.grid(column=0, row=1, columnspan=2, rowspan=1, sticky='N')
        self.splitrefacvar = tkinter.IntVar()
        self.splitrefacvar.set(1000)
        self.splitrefactorentry = tkinter.Entry(self.EncoderFrame, textvariable=self.splitrefacvar)
        self.splitrefactorentry.grid(column=2, row=1, columnspan=1) 
        
        mono = tkinter.Label(self.EncoderFrame, text='Encoder Smoothing',anchor='w')
        mono.grid(column=0, row=2, columnspan=2, rowspan=1, sticky='W')
        
        self.encoder_smoothvar = tkinter.DoubleVar()
        self.encoder_smoothvar.set(251)
        self.encoder_smoothentry = tkinter.Entry(self.EncoderFrame, textvariable=self.encoder_smoothvar)
        self.encoder_smoothentry.grid(column=2, row=2, columnspan=1)    
        
        sampling = tkinter.Label(self.EncoderFrame, text='Sampling Frequency (MHz)',anchor='w')
        sampling.grid(column=0, row=3, columnspan=2, rowspan=1, sticky='W')
        
        self.samplingvar = tkinter.IntVar()
        self.samplingvar.set(2)
        self.samplingentry = tkinter.Entry(self.EncoderFrame, textvariable=self.samplingvar)
        self.samplingentry.grid(column=2, row=3, columnspan=1)   
		
        CPUNum = tkinter.Label(self.EncoderFrame, text='Threads',anchor='w')
        CPUNum.grid(column=0, row=4, columnspan=1, rowspan=1, sticky='N')
        self.CPUNumvar = tkinter.IntVar()
        self.CPUNumvar.set(1)
        self.CPUNumentry = tkinter.Entry(self.EncoderFrame, textvariable=self.CPUNumvar)
        self.CPUNumentry.grid(column=1, row=4, columnspan=2)

        self.abort = False	
        self.active = False	
        self.active_val	= 0 	
        
################################################################################
        
        self.CalibFrame = tkinter.Frame(self, width=300)
        self.CalibFrame.grid(row=5,column=7,rowspan=3,sticky='EW')

        self.buttoncalib= tkinter.Button(self.CalibFrame, text='Calibrate', width=30, height=1,command=self.calibrate)
        self.buttoncalib.grid(column=0, row=0, columnspan=3) 
		
        mono = tkinter.Label(self.CalibFrame, text='Mono Crystal',anchor='w')
        mono.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='N')
		
        self.calibmono = tkinter.StringVar(self.CalibFrame)
        self.calibmono.set('Si 111')
        self.calibmonoMenu = tkinter.OptionMenu(self.CalibFrame, self.calibmono, 'Si 111', 'Si 311')
        self.calibmonoMenu.grid(column=1, row=1, columnspan=1) 
        
        calibener = tkinter.Label(self.CalibFrame, text='Energy (eV)',anchor='w')
        calibener.grid(column=0, row=2, columnspan=1, rowspan=1, sticky='N')
        self.calibenervar = tkinter.DoubleVar()
        self.calibenerentry = tkinter.Entry(self.CalibFrame, textvariable=self.calibenervar)
        self.calibenerentry.grid(column=1, row=2, columnspan=2)   
        
        calibedge = tkinter.Label(self.CalibFrame, text='Edge',anchor='w')
        calibedge.grid(column=0, row=3, columnspan=1, rowspan=1, sticky='N')
        self.calibedgevar = tkinter.StringVar()
        self.calibedgeentry = tkinter.Entry(self.CalibFrame, textvariable=self.calibedgevar)
        self.calibedgeentry.grid(column=1, row=3, columnspan=2)
        
        calibfilterlabel = tkinter.Label(self.CalibFrame, text='Derivative Filter',anchor='w')
        calibfilterlabel.grid(column=0, row=4, columnspan=2, rowspan=1, sticky='E') 
        self.CalibFiltervar = tkinter.DoubleVar()
        self.CalibFiltervar.set(1)
        self.CalibFilterCheck = tkinter.Checkbutton(self.CalibFrame, variable=self.CalibFiltervar)
        self.CalibFilterCheck.grid(column=2, row=4)
        
        calibuselabel = tkinter.Label(self.CalibFrame, text='Use Calibration',anchor='w')
        calibuselabel.grid(column=0, row=5, columnspan=2, rowspan=1, sticky='E') 
        self.CalibUsevar = tkinter.DoubleVar()
        self.CalibUsevar.set(0)
        self.CalibUseCheck = tkinter.Checkbutton(self.CalibFrame, variable=self.CalibUsevar)
        self.CalibUseCheck.grid(column=2, row=5)
		
        self.twod = tkinter.DoubleVar()
        self.twod.set(0.627120)
        
################################################################################
        
        self.optionFrame = tkinter.Frame(self, width=300)
        self.optionFrame.grid(row=7,column=7,rowspan=7,sticky='EW')
		
        self.buttonnorm= tkinter.Button(self.optionFrame, text='Normalise', width=30, height=1,command=self.norm)
        self.buttonnorm.grid(column=0, row=1, columnspan=3)
        
        normuselabel = tkinter.Label(self.optionFrame, text='Use Normalisation',anchor='w')
        normuselabel.grid(column=0, row=2, columnspan=2, rowspan=1, sticky='E') 
        self.NormUsevar = tkinter.DoubleVar()
        self.NormUsevar.set(0)
        self.NormUseCheck = tkinter.Checkbutton(self.optionFrame, variable=self.NormUsevar)
        self.NormUseCheck.grid(column=2, row=2)
        
        Flatlabel = tkinter.Label(self.optionFrame, text='Flatten',anchor='w')
        Flatlabel.grid(column=0, row=3, columnspan=2, rowspan=1, sticky='E') 
        self.Flatvar = tkinter.DoubleVar()
        self.Flatvar.set(0)
        self.FlatCheck = tkinter.Checkbutton(self.optionFrame, variable=self.Flatvar)
        self.FlatCheck.grid(column=2, row=3)
		
        ejlabel = tkinter.Label(self.optionFrame, text='Edge Jump',anchor='w')
        ejlabel.grid(column=0, row=4, columnspan=1, rowspan=1, sticky='E')

        ejvlabel = tkinter.Label(self.optionFrame, text=str(self.edge_jump),anchor='w')
        ejvlabel.grid(column=1, row=4, columnspan=2, rowspan=1, sticky='E') 		
		
        self.buttonloadnorm= tkinter.Button(self.optionFrame, text='Load Normalisation', width=30, height=1,command=self.normload)
        self.buttonloadnorm.grid(column=0, row=5, columnspan=3)

        self.buttonbatch= tkinter.Button(self.optionFrame, text='Batch Export', width=30, height=1,command=lambda :self.process())
        self.buttonbatch.grid(column=0, row=6, columnspan=3)
 
################################################################################       
        self.NormFrame = tkinter.Frame(self, width=300)
        
        self.norm_pre = tkinter.StringVar(self.NormFrame)
        self.norm_pre.set('0')
        self.norm_preMenu = tkinter.OptionMenu(self.NormFrame, self.norm_pre, '0','1','2','3','V')
        self.norm_preMenu.grid(column=2, row=0, columnspan=1) 
        
        self.norm_post = tkinter.StringVar(self.NormFrame)
        self.norm_post.set('V')
        self.norm_postMenu = tkinter.OptionMenu(self.NormFrame, self.norm_post, '0','1','2','3','V')
        self.norm_postMenu.grid(column=2, row=2, columnspan=1) 
        
################################################################################

        self.fig = plt.figure.Figure(figsize=(7, 7), facecolor='0.94')
        
        self.canvasFrame = tkinter.Frame(master=self, width=300)
        self.canvasFrame.grid(row=3,column=4,columnspan=3,rowspan=5,padx=(25, 25),sticky='EW')
		
        self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
        self.canvas.get_tk_widget().grid(column=0, columnspan=3, padx=(0, 0), row=2, rowspan=5, sticky='EW')
        
        self.ax=self.fig.add_subplot(111)
        
        x,y=[1,2,3],[4,5,6]
        self.ax.plot(x, y)
        self.ax.set_xlabel('Energy (eV)')
        self.ax.set_ylabel('Absorption')
        self.canvas.draw_idle()
        
        self.toolbarFrame = tkinter.Frame(master=self, width=300)
        self.toolbarFrame.grid(row=8,column=4,columnspan=3,padx=(25, 25),pady=(10,10),rowspan=1,sticky='W')
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
        self.toolbar.update()

        self.axis_color = 'lightgoldenrodyellow'	

################################################################################
    def close_figure(self):
        self.fig.clear()
        self.fig = plt.figure.Figure(figsize=(7, 7), facecolor='0.94')
        self.ax=self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
        self.canvas.get_tk_widget().grid(column=0, columnspan=3, padx=(25, 0), row=2, rowspan=5, sticky='EW')
        self.toolbar.destroy()
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
        self.toolbar.update()
        
################################################################################
    def rootsfile_check(self, folder):
        try:
            self.calibroots = (pd.read_csv(folder+'/roots.dat', sep='\t', header=None)).values
            return True
        except:
            print('Encoder Analysis Required')
            return False
        
################################################################################
    def calibfile_check(self, folder):
        try:
            self.calibration_values = (pd.read_csv(folder+'/calibration.dat', sep='\t', header=None)).values
            self.CalibUsevar.set(1)
            return True
        except:
            print('Cannot Find Calibration Parameter File')
            self.CalibUsevar.set(0)
            return False
        
################################################################################
    def normfile_check(self, folder):
        try:
            self.normalisation_values = (pd.read_csv(folder+'/Normalisation.dat', sep='\t', header=None)).values
            self.NormUsevar.set(1)
            return True
        except:
            print('Cannot Find Normalisation Parameter File')
            self.NormUsevar.set(0)
            return False
        
################################################################################        
    def normalize(self,energy, mu):
        pre1,pre2,post1,post2,order_pre,order_post,e0 = self.normalisation_values
        
        try:
            self.pre1_scroll.set(float(pre1))
            self.pre2_scroll.set(float(pre2))
            self.post1_scroll.set(float(post1))
            self.post2_scroll.set(float(post2))
            self.norm_pre.set(str(order_pre))
            self.norm_post.set(str(order_post))
        except:
            self.pre1_init = float(pre1[0])
            self.pre2_init = float(pre2[0])
            self.post1_init = float(post1[0])
            self.post2_init = float(post2[0])
            try:
                self.norm_pre_init = str(int(order_pre[0]))
                self.norm_post_init = str(int(order_post[0]))
            except:
                self.norm_pre_init = order_pre[0]
                self.norm_post_init = order_post[0]
			
        self.e0 = float(e0[0])
        
        lineprey = self.regression('pre-edge',np.asarray(energy), np.asarray(mu))        
        lineposty = self.regression('post-edge',np.asarray(energy), np.asarray(mu))
        
        index = min(range(len(energy)), key=lambda i: abs(energy[i]-self.e0))
        normdivisor = (lineposty-lineprey)[index]

        mu = (mu-lineprey)/normdivisor
        lineposty = (lineposty-lineprey)/normdivisor
        lineprey = np.zeros(len(lineprey))
                        
        if self.Flatvar.get() == 1:
            flat_correction = 1-lineposty
            if energy[index] > energy[-1]:
                flat_correction[index:int(len(energy))] = 0
            else:
                flat_correction[0:index] = 0
            mu = mu + flat_correction
            return mu
        else:
            return mu
################################################################################
    def plot_spectrum(self):
        self.initialize = True
        self.ax.clear()
        self.fig = plt.figure.Figure(figsize=(7, 7), facecolor='0.94')
        self.ax=self.fig.add_subplot(111)
        
        self.rootsfile = self.rootsfile_check(self.folder)
        self.calibrated = self.calibfile_check(self.folder)
        self.normalized = self.normfile_check(self.folder)
        
        g = open(self.data_bin+'.bin', 'rb')
        data = np.fromfile(g, dtype=np.int32, count = 2)
        self.headerSize = data[0]
        data = np.fromfile(g, dtype=np.int64, count = 1)
        data = np.fromfile(g, dtype='f4', count = 2)
        self.nChannels = int(data[0])
        
        def plot_data(evt):
            filter_toggle = int(self.Filtervar.get())
            
            try:
                i = int(self.ds_scroll_scale.get())
                self.Wn = float(self.bf_scroll_scale.get()) # Cutoff frequency
            except:
                i = 0
                self.Wn = 0.05
                self.initialize = True
           
            pdroots = pd.read_csv(self.folder+'/roots.dat', sep='\t', header=None)
            self.calibroots = pdroots.iloc[:,0].values
            ang, mu = self.spectrumread(self.calibroots[i+3],self.calibroots[2])
            if filter_toggle == 1:
                N  = 3    # Filter order
                B, A = signal.butter(N, self.Wn, output='ba')
                #C, D = signal.butter(N, [0.009, 0.2], btype='band', output='ba')
            
                mu_filtered = signal.filtfilt(B,A,mu)
                #sine_filter = signal.filtfilt(C,D,mu_filtered)
                mu_filtered = mu_filtered
            else:
                mu_filtered = mu
            
            if self.initialize == True:
                if self.CalibUsevar.get() == 1:
                    energy = np.transpose(1239.852/(float(self.calibration_values[3])*np.sin((ang+(self.calibration_values[1]-self.calibration_values[0]))*np.pi/180))).tolist()
                    return energy, mu, mu_filtered
                else:
                    return ang, mu, mu_filtered
            else:
                if self.CalibUsevar.get() == 1:
                    energy = np.transpose(1239.852/(float(self.calibration_values[3])*np.sin((ang+(self.calibration_values[1]-self.calibration_values[0]))*np.pi/180))).tolist()
                    if self.NormUsevar.get() == 1:
                        mu = self.normalize(energy, mu)
                        mu_filtered = self.normalize(energy, mu_filtered)
                        line.set_ydata(mu)
                        line2.set_ydata(mu_filtered)
                    else:
                        line.set_ydata(mu)
                        line2.set_ydata(mu_filtered)
                            
                    line.set_xdata(energy)
                    line2.set_xdata(energy)
                    self.ax.set_ylim([np.min(mu_filtered)-0.1, np.max(mu_filtered)+0.1])
                    if self.ax.get_xlim()[0] < 1000:
                        self.ax.set_xlim([np.min(energy), np.max(energy)])
                    self.canvas.draw_idle()
                else: 
                    line.set_ydata(mu)
                    line.set_xdata(ang)
                    line2.set_ydata(mu_filtered)
                    line2.set_xdata(ang)
                    self.ax.set_ylim([np.min(mu_filtered)-0.1, np.max(mu_filtered)+0.1])
                    if self.ax.get_xlim()[0] > 1000:
                        self.ax.set_xlim([np.max(ang), np.min(ang)])
                    self.canvas.draw_idle()

        if self.rootsfile == True:
            if self.initialize == True:
                x, mu, mu_filtered = plot_data(None)
                self.initialize = False
                [line] = self.ax.plot(x, mu)
                [line2] = self.ax.plot(x, mu_filtered, lw=2)
                self.ax.set_ylim([np.min(mu), np.max(mu)])
                if self.CalibUsevar.get() == 1:
                    self.ax.set_xlim([np.min(x), np.max(x)])
                else:
                    self.ax.set_xlim([np.max(x), np.min(x)])
                self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
                self.canvas.get_tk_widget().grid(column=0, columnspan=3, padx=(25, 0), row=2, rowspan=5, sticky='EW')
                self.toolbar.destroy()
                self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
                self.toolbar.update()
                self.ds_scroll_scale = tkinter.Scale(self.canvasFrame, orient='horizontal', from_=0, to=int(len(self.calibroots)-5), command=plot_data)
                self.ds_scroll_scale.grid(column=1, columnspan=2, row=8, stick='EW')
                self.bf_scroll_scale = tkinter.Scale(self.canvasFrame, orient='horizontal', from_=0.0005, to=0.995, digits = 4, resolution = 0.0005, command=plot_data)
                self.bf_scroll_scale.grid(column=1, columnspan=2, row=9, stick='EW')
                if self.CalibUsevar.get() == 1:
                    self.bf_scroll_scale.set(float(self.calibration_values[2]))
                else:
                    try:
                        self.bf_scroll_scale.set(float(self.calibration_values[2]))
                    except:
                        self.bf_scroll_scale.set(0.05)
            
                self.datasetlabel = tkinter.Label(self.canvasFrame, text='Dataset',anchor='w')
                self.datasetlabel.grid(column=0, row=8, columnspan=1, rowspan=1, sticky='S')
                self.Filterlabel = tkinter.Label(self.canvasFrame, text='Butterworth Filter',anchor='w')
                self.Filterlabel.grid(column=0, row=9, columnspan=1, rowspan=1, sticky='S')
                
                self.canvas.draw()
                self.update()
            else:
                print('Encoder Analysis required')

################################################################################     
    def fileplus(self):
        self.askloadfile()
		
################################################################################     
    def normload(self):
        file = askopenfile()
        dirname, file_short = os.path.split(os.path.abspath(file.name))
        shutil.copyfile(file.name,self.folder+'/'+file_short)
        self.Flatvar.set(1)
        self.plot_spectrum()       
		
################################################################################     
    def abortsplit(self):
        self.abort = True

################################################################################
    def fileminus(self):        
        selection = self.listbox.curselection()
        for i in reversed(selection):
            self.listbox.delete(i)
            del self.files_long[i]
            del self.edge_info[i]
            del self.edge_type[i]
        
        self.update()

################################################################################      
    def edge_box(self, file):
        self.box = tkinter.Toplevel()
        label0 =tkinter.Label(self.box, text='Enter Absorbing Atom')
        label0.pack()
	
        self.box.entry0 = tkinter.Entry(self.box)
        self.box.entry0.pack()
	
        label1 = tkinter.Label(self.box, text='Enter Edge type')
        label1.pack()
	
        self.box.entry1 = tkinter.Entry(self.box)
        self.box.entry1.pack()
	
        button1 = tkinter.Button(self.box, text='Submit', command=lambda:self.submit_edge(self.box, file))
        button1.pack()
	
        button2 = tkinter.Button(self.box, text='Cancel', command=lambda:self.box.destroy())
        button2.pack()
 
################################################################################          
    def submit_edge(self, box, file):
        self.central_atom = box.entry0.get()
        self.edge_type.append(box.entry0.get()+'_'+box.entry1.get())
        
        xdb = XrayDB()
        self.edge_info.append(xdb.xray_edge(box.entry0.get(), box.entry1.get()))
        print(self.edge_info[-1])
        self.core_hole = xdb.corehole_width(box.entry0.get(), box.entry1.get())
        print('core hole lifetime',self.core_hole)
        
        self.files_long.append(file)
        self.edge_submitted = True
        dirname, file_short = os.path.split(os.path.abspath(file.name))
        self.listbox.insert(END, file_short)
        
        try:
            self.calibedgevar.set(self.edge_type[self.listbox.curselection()[0]])
            self.calibenervar.set(self.edge_info[self.listbox.curselection()[0]].edge)
        except:
            self.calibedgevar.set('')
            self.calibenervar.set('')
            
        box.destroy()
        self.update()

################################################################################
    def askloadfile(self):
        self.edge_submitted = False
        file = askopenfile()
        try:
            dirname, file_short = os.path.split(os.path.abspath(file.name))
            #file_short = (str(file).split("'")[1]).split("/")[-1]
            print('Loading File')
            self.edge_box(file)
        except:
            print('Aborted Load')
    
################################################################################
    def spectrumread(self,rooti,minr):
        encoder_bin = self.data_bin+'_Encoder'
        minr = int(minr)
            
        g = open(self.data_bin+'.bin', 'rb')
        
        g.seek(0)

        data = np.fromfile(g, dtype=np.int32, count = 2)

        self.headerSize = data[0]
        #print("HeaderSize: ", self.headerSize)
        #print("FileVersion: ", fileVersion)

        data = np.fromfile(g, dtype=np.int64, count = 1)
        #print("DateTime: ", DT)

        data = np.fromfile(g, dtype='f4', count = 2)

        self.nChannels = data[0]
        self.AdcClock_Hz = data[1]
        self.DacClock = 1

        g.seek(int(self.headerSize+(4*self.nChannels*rooti)))
        data_D = (np.fromfile(g, dtype='f4', count = int(self.nChannels*minr)))
        #print(len(data_D))
        #print(len(data_D)/self.nChannels)
        #print(len(data_D[0::int(self.nChannels)]))
        #print(len(data_D[1::int(self.nChannels)]))
        #print(len(data_D[2::int(self.nChannels)]))
        toggle = int(self.togglevar.get())
        mu_r = np.zeros((minr, 1))
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
                mu_r[:,0] = np.log(data_D[numercol::int(self.nChannels)])
            else:
                mu_r[:,0] = np.log(data_D[numercol::int(self.nChannels)]/data_D[denomcol::int(self.nChannels)])
        else:
            if denomcol == int(len(self.choices)-1):
                mu_r[:,0] = data_D[numercol::int(self.nChannels)]
            else:
                mu_r[:,0] = data_D[numercol::int(self.nChannels)]/data_D[denomcol::int(self.nChannels)]
		
        nans = np.argwhere(np.isnan(mu_r[:,0]))
        mu = np.zeros((minr-len(nans), 1))
        mu[:,0] = np.delete(mu_r[:,0], nans)

        f = open(encoder_bin+'.bin', 'rb')
        f.seek(int(self.headerSize+(4*rooti)))
        ang_r = np.fromfile(f, dtype='f4', count = int(minr))
        ang = np.delete(ang_r, nans)

        RawData = pd.DataFrame()
        RawData['ang'] = ang
        RawData['mu'] = mu
        self.angf = ang[0]
        self.angl = ang[-1]
        
        RawData=RawData[(RawData['ang'] >= self.calibroots[0]) & (RawData['ang'] <= self.calibroots[1])]  
        RawData = RawData.groupby('ang', as_index=False).mean()
        
        return RawData['ang'].values, RawData['mu'].values
    
################################################################################
    def onselect(self,evt):
        try:
            data_bin_long = self.files_long[self.listbox.curselection()[0]].name
            self.active = True
            self.active_val = self.listbox.curselection()[0]
        except:
            if self.active == False:
                print('No File Selected')
        if self.active == True:
            data_bin_long = self.files_long[self.active_val].name
            self.calibedgevar.set(self.edge_type[self.active_val])
            self.calibenervar.set(self.edge_info[self.active_val].edge)
        try:
            self.NormFrame.grid_remove()
            self.EncoderFrame.grid(row=1,column=7,rowspan=3,sticky='EW')
            self.CalibFrame.grid(row=5,column=7,rowspan=3,sticky='EW')
            dirname, file_short = os.path.split(os.path.abspath(data_bin_long))
            print(file_short+' File Selected')
        except:
            dirname, file_short = os.path.split(os.path.abspath(data_bin_long))
            print(file_short+' File Selected')
        
        self.ax.clear()
        self.update()
        
        if data_bin_long.split('.')[-1] == 'bin':
            #self.data_bin = data_bin_long.split('.')[0:l]
            self.data_bin = os.path.splitext(data_bin_long)[0]
            self.edge = self.edge_type[self.listbox.curselection()[0]]
            self.folder = 'output_'+self.edge+'_'+self.data_bin.split('/')[-1]

            f = open(self.data_bin+'.bin', 'rb')
            f.seek(0)

            data = np.fromfile(f, dtype=np.int32, count = 2)
            self.headerSize = data[0]
            data = np.fromfile(f, dtype=np.int64, count = 1)
            data = np.fromfile(f, dtype='f4', count = 2)
            self.nChannels = data[0]
            self.AdcClock_Hz = data[1]
            self.samplingvar.set(float(self.AdcClock_Hz/1E6))
            
            self.choices = []
            for i in range(int(self.nChannels)):
                self.choices.append('CH '+str(int(i)))
            self.choices.append('1')
            
            self.sam_numerMenu.destroy()
            self.sam_denomMenu.destroy()
            self.ref_numerMenu.destroy()
            self.ref_denomMenu.destroy()
            
            self.sam_numer = tkinter.StringVar(self.plotcolsFrame)
            self.sam_numer.set('CH 0')
            self.sam_numerMenu = tkinter.OptionMenu(self.plotcolsFrame, self.sam_numer, *self.choices)
            self.sam_numerMenu.grid(column=1, row=1, columnspan=1) 
        
            self.sam_denom = tkinter.StringVar(self.plotcolsFrame)
            self.sam_denom.set('CH 1')
            self.sam_denomMenu = tkinter.OptionMenu(self.plotcolsFrame, self.sam_denom, *self.choices)
            self.sam_denomMenu.grid(column=1, row=2, columnspan=1) 
        
            self.ref_numer = tkinter.StringVar(self.plotcolsFrame)
            self.ref_numer.set('CH 1')
            self.ref_numerMenu = tkinter.OptionMenu(self.plotcolsFrame, self.ref_numer, *self.choices)
            self.ref_numerMenu.grid(column=2, row=1, columnspan=1) 
        
            self.ref_denom = tkinter.StringVar(self.plotcolsFrame)
            self.ref_denom.set('CH 2')
            self.ref_denomMenu = tkinter.OptionMenu(self.plotcolsFrame, self.ref_denom, *self.choices)
            self.ref_denomMenu.grid(column=2, row=2, columnspan=1)
            
            toggle = [('',0),('',1)]
            self.radio_plot_buttons = []
            for idx, (text, mode) in enumerate(toggle):
                self.radio_plot_buttons.append(tkinter.Radiobutton(self.plotcolsFrame, text=text, variable=self.togglevar, value=mode))
                self.radio_plot_buttons[-1].grid(row=4, column=idx+1) 
            
            self.clear_canvasFrame()
            self.initialize = True
            
            if not os.path.exists(self.folder):
                os.makedirs(self.folder)
                self.update()
            else:
                self.rootsfile = self.rootsfile_check(self.folder)
                if self.rootsfile == True:
                    self.plot_spectrum()
                else:
                    print('Encoder Analysis Required')
                    self.update()
        
        if data_bin_long.split('.')[-1] == 'txt':
            print('txt import not yet supported')
        
################################################################################             
    def plotroots(self,roots,angs):
        #loop finds the minimum and maximum difference between roots
        self.clear_canvasFrame()
        self.roots = roots
        self.angs = angs
        global test

        #for k in range(len(self.roots)-1):
        #    test[k] = self.roots[k+1] - self.roots[k]

        test = np.diff(self.roots)
		
        minr = int(np.min(test))
        maxr = int(np.max(test))
  
        self.ax.clear()
        self.fig = plt.figure.Figure(figsize=(7, 7))
        self.fig.subplots_adjust(bottom=0.2)
        self.ax=self.fig.add_subplot(111)
        [scatter] = self.ax.plot(range(len(test)), test, 'o')
        self.ax.set_xlabel('Spectrum Number')
        self.ax.set_ylabel('Datapoints per Spectrum')
        self.ax.set_ylim([0,maxr*1.1])
        self.ax.set_xlim([0,len(test)])
        
        self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
        self.canvas.get_tk_widget().grid(column=0, columnspan=3, padx=(25, 0), row=2, rowspan=5, sticky='EW')
        
        self.toolbar.destroy()
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
        self.toolbar.update()

        failed_button_ax = self.fig.add_axes([0.5, 0.025, 0.1, 0.04])
        failed_button = Button(failed_button_ax, 'Failed', color=self.axis_color, hovercolor='0.975')

        def failed_button_on_clicked(mouse_event):
            global test
            Ave = np.average(test)
            print(Ave)
            global store_k
            store_k = []
            for k in range(len(roots)-1):
                if test[k] < 0.9*Ave:
                    store_k.append(k-1)
                    store_k.append(k)
					
            store_k = list(set(store_k))
            test = [i for j, i in enumerate(test) if j not in store_k]
            self.angs = [i for j, i in enumerate(self.angs) if j not in store_k]
            self.roots = [i for j, i in enumerate(self.roots) if j not in store_k]
            global minr
            minr = int(np.min(test))
            print('Data Points Per Spectra:', minr)
            scatter.set_ydata(test)
            scatter.set_xdata(range(len(test)))
            self.ax.set_ylim([0,maxr*1.1])
            self.canvas.draw_idle()

        failed_button.on_clicked(failed_button_on_clicked)

        insert_button_ax = self.fig.add_axes([0.6, 0.025, 0.1, 0.04])
        insert_button = Button(insert_button_ax, 'Insert', color=self.axis_color, hovercolor='0.975')
        def insert_button_on_clicked(mouse_event):
            global test
            for i in range(len(test)):
                Ave = np.average(test)
	
                if test[i] > 1.8*Ave:
                    np.insert(self.roots, i+1, self.roots[i]+Ave)
                    np.insert(test, i+1, Ave)
                    test[i] = test[i] - Ave
	
            minr = int(np.min(test))
            maxr = int(np.max(test))
            print('Data Points Per Spectra:', minr)
            scatter.set_ydata(test)
            scatter.set_xdata(range(len(test)))
            self.ax.set_ylim([0,maxr*1.1])
	
            self.canvas.draw_idle()

        insert_button.on_clicked(insert_button_on_clicked)

        accept_button_ax = self.fig.add_axes([0.7, 0.025, 0.1, 0.04])
        accept_button = Button(accept_button_ax, 'Accept', color=self.axis_color, hovercolor='0.975')
        def accept_button_on_clicked(mouse_event):      
            print(np.average(self.angs))
            global test
            minr = int(np.min(test))
            B = [i for i in self.angs if i >= np.average(self.angs)]
            C = [i for i in self.angs if i <= np.average(self.angs)]
            if B[0] > C[0]:
                max_ang = min(B)-0.0025
                min_ang = max(C)+0.0025
            else:
                max_ang = min(C)-0.0025
                min_ang = max(B)+0.0025
                
            if min_ang > max_ang:
                self.roots = np.concatenate([[max_ang],[min_ang],[minr],self.roots])
            else:    
                self.roots = np.concatenate([[min_ang],[max_ang],[minr],self.roots])
            roots_info = pd.DataFrame(self.roots)
            roots_info.to_csv(self.folder+'/roots.dat', sep='\t', header=None, index=False)
            
            self.plot_spectrum()

        accept_button.on_clicked(accept_button_on_clicked)

        reject_button_ax = self.fig.add_axes([0.8, 0.025, 0.1, 0.04])
        reject_button = Button(reject_button_ax, 'Reject', color=self.axis_color, hovercolor='0.975')
        def reject_button_on_clicked(mouse_event):
            self.close_figure()

        reject_button.on_clicked(reject_button_on_clicked)
        
        self.canvas.draw()
        app.mainloop()
        
        return minr
    
################################################################################	 
    def clear_canvasFrame(self):
        list = self.canvasFrame.grid_slaves()
        for child in list:
            if str(type(child)) != "<class 'tkinter.Canvas'>":
                child.destroy()
		
################################################################################           
    def regression(self,region,x,y): 
        if region == 'pre-edge':
            try:
                xllim = self.e0_scroll_scale.get() + self.pre1_scroll_scale.get()
                xhlim = self.e0_scroll_scale.get() + self.pre2_scroll_scale.get()
                order = self.norm_pre.get()
            except:
                try:
                    xllim = self.e0 + self.pre1_init
                    xhlim = self.e0 + self.pre2_init
                except:
                    xllim = self.calibenervar.get() + self.pre1_init
                    xhlim = self.calibenervar.get() + self.pre2_init
                try:
                    order = self.norm_pre_init[0]
                except:
                    order = str(int(self.norm_pre_init))
            
        if region == 'post-edge':
            try:
                xllim = self.e0_scroll_scale.get() + self.post1_scroll_scale.get()
                xhlim = self.e0_scroll_scale.get() + self.post2_scroll_scale.get()
                order = self.norm_post.get()
            except:
                try:
                    xllim = self.e0 + self.post1_init
                    xhlim = self.e0 + self.post2_init
                except:
                    xllim = self.calibenervar.get() + self.post1_init
                    xhlim = self.calibenervar.get() + self.post2_init
                try:
                    order = self.norm_post_init[0]
                except:
                    order = str(int(self.norm_post_init))
                
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
################################################################################	
################################################################################		
    def split(self):
        self.abort = False
        self.progress_bar["value"] = 0
        self.progress_bar["style"] = "blue.Horizontal.TProgressbar"
        self.clear_canvasFrame()

        print('Anaylzing Encoder')
        
        self.edge = self.calibedgevar.get()
        buffer = 32000000
        print(self.data_bin)
        encoder_bin = self.data_bin+'_Encoder'
		
        if not os.path.exists('output_'+self.edge+'_'+self.data_bin.split('/')[-1]):
            os.makedirs('output_'+self.edge+'_'+self.data_bin.split('/')[-1])
            
        self.folder = 'output_'+self.edge+'_'+self.data_bin.split('/')[-1]

        f = open(encoder_bin+'.bin', 'rb')
        f.seek(0)
        data = np.fromfile(f, dtype=np.int32, count = 2)
        self.headerSize = data[0]
        fileVersion = data[1]
        print("HeaderSize: ", self.headerSize)
        print("FileVersion: ", fileVersion)
        data = np.fromfile(f, dtype=np.int64, count = 1)
        DT = str(data[0])
        print("DateTime: ", DT)
        data = np.fromfile(f, dtype='f4', count = 2)
        self.nChannels = data[0]
        self.AdcClock_Hz = data[1]
        self.samplingvar.set(float(self.AdcClock_Hz/1E6))
        self.DacClock = 1
        print("nChannels:", self.nChannels)
        print("AdcClock_Hz:", self.AdcClock_Hz)
        nData = int((os.path.getsize(encoder_bin+'.bin') - self.headerSize)/4)
        print("nData", nData)
        buffer_points=np.arange(np.ceil(4*nData/buffer))

        start_initial = time.time()
        #Cycles through Encoder file with set size of data to read in
        resample_factor = int(self.splitrefacvar.get())
        enconder_sampling = self.samplingvar.get()
        encoder_smooth = self.encoder_smoothvar.get()
        if (encoder_smooth % 2) == 0:
            encoder_smooth = encoder_smooth - 1
        ncpus = int(self.CPUNumvar.get())		
		
        if len(buffer_points) < ncpus:
            ncpus = len(buffer_points)
        print('num cpus =',ncpus)
        print('job length =',len(buffer_points))
        self.labelText.set('Starting Processess')
        self.update()
        time.sleep(2)
		
        angs2d = [None] * ncpus
        roots2d = [None] * ncpus		
        options = [ncpus, np.max(buffer_points), self.headerSize, buffer, encoder_bin, nData, resample_factor, enconder_sampling, encoder_smooth]       
        sub_f = 'batch_split_subroutine_v49.py'
        process = subprocess.Popen(['python', sub_f, str(options)], stdout=subprocess.PIPE, shell=True)
        
        while process.poll() == None:
            stdoutdata = process.stdout.readline()
            if self.labelText.get() == 'Starting Processess':
                self.labelText.set('Running Processess')
            if 'Affinity' in stdoutdata.decode('utf-8'):
                print(stdoutdata.decode('utf-8').split('\n')[0])
            if 'Percent' in stdoutdata.decode('utf-8'):
                perc = float(stdoutdata.decode('utf-8').split('\n')[0].split(' ')[-1])
                self.progress_bar["value"] = perc
                self.progress_bar.update()
            if 'Returning' in stdoutdata.decode('utf-8'):
                self.progress_bar["value"] = 100
                self.labelText.set('Finishing Up')
                self.progress_bar.update()
            if 'angs' in stdoutdata.decode('utf-8'):
                _ , ID, data = stdoutdata.decode('utf-8').split('\n')[0].split('_')
                angs2d[int(ID)] = ast.literal_eval(data)
            if 'roots' in stdoutdata.decode('utf-8'):
                _ , ID, data = stdoutdata.decode('utf-8').split('\n')[0].split('_')
                roots2d[int(ID)] = ast.literal_eval(data)
            if 'Waiting' in stdoutdata.decode('utf-8'):  
                self.update()	
            if 'Process' in stdoutdata.decode('utf-8'):  
                print(stdoutdata.decode('utf-8').split('\n')[0])					
            if 'Returned' in stdoutdata.decode('utf-8'):
                self.progress_bar["value"] = 100
                self.labelText.set('Joining Outputs')
                self.progress_bar.update()
                break
            self.update()
            #print(stdoutdata.decode('utf-8').split('\n')[0])
		
        print('Joining Output into list')
        roots = np.asarray(list(itertools.chain(*roots2d)))
        angs = np.asarray(list(itertools.chain(*angs2d)))	
        roots = roots.astype(np.int64)
        print('number of roots:',len(roots))
        self.minr = self.plotroots(roots,angs)

################################################################################	
################################################################################	
################################################################################	
    def calibrate(self):
        self.Wn = float(self.bf_scroll_scale.get())
        self.progress_bar["value"] = 0
        self.clear_canvasFrame()
        
        print('Initial Calibration')
        self.initialize = True
        self.ax.clear()
        self.fig = plt.figure.Figure(figsize=(7, 7), facecolor='0.94')
        self.fig.subplots_adjust(bottom=0.2)
        self.ax=self.fig.add_subplot(111)
        
        #data_bin_long = self.files_long[self.listbox.curselection()[0]].name
        
        #if data_bin_long.split('.')[-1] == 'bin':
        #    data_bin = data_bin_long.split('.')[0]

        if not os.path.exists(self.folder+'/roots.dat'):
            print('Encoder Splitting Required')
        else:
            pdroots = pd.read_csv(self.folder+'/roots.dat', sep='\t', header=None)
            self.calibroots = pdroots.iloc[:,0].values
        
        def plot_data(evt):
            filter_toggle = int(self.CalibFiltervar.get())

            if self.initialize == False:
                i = int(self.ds_scroll_scale.get())+1
                self.Wn = float(self.bf_scroll_scale.get()) # Cutoff frequency
            else:
                try:
                    self.Wn = float(self.bf_scroll_scale.get()) # Cutoff frequency
                    i = int(self.ds_scroll_scale.get())
                except:
                    self.Wn = self.Wn
                    i = 0
                
            ang, mu = self.spectrumread(self.calibroots[i+3],self.calibroots[2]) 
            
            if filter_toggle == 1:
                N  = 3    # Filter order
                B, A = signal.butter(N, self.Wn, output='ba')
                #C, D = signal.butter(N, [0.009, 0.2], btype='band', output='ba')
            
                mu_filtered = signal.filtfilt(B,A,mu)
                #sine_filter = signal.filtfilt(C,D,mu_filtered)
                mu_filtered = mu_filtered
            else:
                mu_filtered = mu
            
            if self.initialize == True:
                return ang, mu, mu_filtered
            else:
                if filter_toggle == 1:
                    ply = -np.gradient(mu_filtered)
                    line.set_ydata(ply)
                    line.set_xdata(ang)
                else:
                    ply = -np.gradient(mu)
                    line.set_ydata(ply)
                    line.set_xdata(ang)
                self.calibline.set_xdata(self.cp_scroll_scale.get())
                self.ax.autoscale(axis='y', tight=False)
                self.canvas.draw_idle()

        if self.initialize == True:
		
            if self.calibmono.get() == 'Si 111':
                self.twod.set(0.627120)
                print('Mono Crystal = Si 111')
            elif self.calibmono.get() == 'Si 311':
                self.twod.set(0.320267)
                print('Mono Cyrstal = Si 311')
				
            ang, mu, mu_filtered = plot_data(0)
            filter_toggle = int(self.CalibFiltervar.get())
            if filter_toggle == 1:
                print('Plotting Filtered')
                ply = -np.gradient(mu_filtered)
            else:
                ply = -np.gradient(mu)
				
            angmax = ang[np.where(ply == np.max(ply))]
            print(angmax[0])
            print(np.max(ply))
            [line] = self.ax.plot(ang, ply) 
            self.calibline = self.ax.axvline(ang[np.argmax(-np.gradient(mu))], ls='--', c='orange')
            self.ax.autoscale(axis='y', tight=False)
            self.ax.invert_xaxis()
            self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
            self.canvas.get_tk_widget().grid(column=0, columnspan=3, padx=(25, 0), row=2, rowspan=5, sticky='EW')
            self.toolbar.destroy()
            self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
            self.toolbar.update()
            self.ds_scroll_scale = tkinter.Scale(self.canvasFrame, orient='horizontal', from_=0, to=int(len(self.calibroots)-5), command=plot_data)
            self.ds_scroll_scale.grid(column=1, columnspan=2, row=8, stick='EW')
            self.bf_scroll_scale = tkinter.Scale(self.canvasFrame, orient='horizontal', from_=0.0005, to=0.995, digits = 5, resolution = 0.0005, command=plot_data)
            self.bf_scroll_scale.grid(column=1, columnspan=2, row=9, stick='EW')
            self.bf_scroll_scale.set(self.Wn)
            
            self.cp_scroll_scale = tkinter.Scale(self.canvasFrame, orient='horizontal', from_=np.max(ang), to=np.min(ang), digits = 6, resolution = 0.00001, command=plot_data)
            self.cp_scroll_scale.grid(column=1, columnspan=2, row=10, stick='EW')
            self.cp_scroll_scale.set(angmax[0])
            
            self.datasetlabel = tkinter.Label(self.canvasFrame, text='Dataset',anchor='w')
            self.datasetlabel.grid(column=0, row=8, columnspan=1, rowspan=1, sticky='S')
            self.Filterlabel = tkinter.Label(self.canvasFrame, text='Butterworth Filter',anchor='w')
            self.Filterlabel.grid(column=0, row=9, columnspan=1, rowspan=1, sticky='S')
            self.CPlabel = tkinter.Label(self.canvasFrame, text='Calibration Point',anchor='w')
            self.CPlabel.grid(column=0, row=10, columnspan=1, rowspan=1, sticky='S')
            
            self.initialize = False
            
            accept_button_ax = self.fig.add_axes([0.7, 0.025, 0.1, 0.04])
            accept_button = Button(accept_button_ax, 'Accept', color=self.axis_color, hovercolor='0.975')
            def accept_button_on_clicked(mouse_event):
                self.thetac = np.arcsin(1239.852/(self.twod.get()*self.calibenervar.get()))*180/np.pi
                calibvalue = [self.cp_scroll_scale.get(), self.thetac, float(self.bf_scroll_scale.get()), self.twod.get(), self.calibenervar.get()]
                pdcalibvalue = pd.DataFrame(calibvalue)
                pdcalibvalue.to_csv(self.folder+'/calibration.dat', sep='\t', header=None, index=False)
                self.close_figure()
                self.clear_canvasFrame()
                self.CalibUsevar.set(1)
                self.plot_spectrum()

            accept_button.on_clicked(accept_button_on_clicked)

            reject_button_ax = self.fig.add_axes([0.8, 0.025, 0.1, 0.04])
            reject_button = Button(reject_button_ax, 'Reject', color=self.axis_color, hovercolor='0.975')
            def reject_button_on_clicked(mouse_event):
                self.close_figure()

            reject_button.on_clicked(reject_button_on_clicked)
        
            self.canvas.draw()
            self.update()
            app.mainloop()

################################################################################	
################################################################################	
################################################################################	
    def norm(self):
        self.Wn = float(self.bf_scroll_scale.get())
        self.progress_bar["value"] = 0
        print('Normalising Data')
        self.clear_canvasFrame()

        self.initialize = True
        self.ax.clear()
        self.fig = plt.figure.Figure(figsize=(7, 7), facecolor='0.94')
        self.fig.subplots_adjust(bottom=0.2)
        self.ax=self.fig.add_subplot(111)
        
        #data_bin_long = self.files_long[self.listbox.curselection()[0]].name
        
        #if data_bin_long.split('.')[-1] == 'bin':
        #    data_bin = data_bin_long.split('.')[0]

        if not os.path.exists(self.folder+'/roots.dat'):
            print('Encoder Splitting Required')
        if not os.path.exists(self.folder+'/calibration.dat'):
            print('Calibration Required')
        else:
            self.calibroots = (pd.read_csv(self.folder+'/roots.dat', sep='\t', header=None)).iloc[:,0].values
            self.calibration_values = (pd.read_csv(self.folder+'/calibration.dat', sep='\t', header=None)).values
            
            def plot_data(evt):
                filter_toggle = int(self.Filtervar.get())
            
                if self.initialize == False:
                    i = int(self.ds_scroll_scale.get())
                    self.Wn = float(self.bf_scroll_scale.get()) # Cutoff frequency
                else:
                    try:
                        self.Wn = float(self.bf_scroll_scale.get()) # Cutoff frequency
                        i = int(self.ds_scroll_scale.get())
                    except:
                        self.Wn = self.Wn
                        i = 0
                
                ang, mu = self.spectrumread(self.calibroots[i+3],self.calibroots[2])
                
                energy = np.transpose(1239.852/(float(self.calibration_values[3])*np.sin((ang+(self.calibration_values[1]-self.calibration_values[0]))*np.pi/180))).tolist()
                    
                if filter_toggle == 1:
                    N  = 3    # Filter order
                    B, A = signal.butter(N, self.Wn, output='ba')
                    #C, D = signal.butter(N, [0.009, 0.2], btype='band', output='ba')
            
                    mu_filtered = signal.filtfilt(B,A,mu)
                    #sine_filter = signal.filtfilt(C,D,mu_filtered)
                    mu_filtered = mu_filtered
                else:
                    mu_filtered = mu
            
                if self.initialize == True:
                    return energy, mu, mu_filtered
                else:
                    if filter_toggle == 1:
                        ply = mu_filtered
                    else:
                        ply = mu
                
                    self.pre1line.set_xdata(self.pre1_scroll_scale.get()+self.e0_scroll_scale.get())
                    self.pre2line.set_xdata(self.pre2_scroll_scale.get()+self.e0_scroll_scale.get())
                    self.post1line.set_xdata(self.post1_scroll_scale.get()+self.e0_scroll_scale.get())
                    self.post2line.set_xdata(self.post2_scroll_scale.get()+self.e0_scroll_scale.get())
                    
                    lineprey = self.regression('pre-edge', np.asarray(energy), np.asarray(ply))
                    lineposty = self.regression('post-edge', np.asarray(energy), np.asarray(ply))
                    index = min(range(len(energy)), key=lambda i: abs(energy[i]-self.e0_scroll_scale.get()))
                    self.edge_jump = (lineposty-lineprey)[index]
                    ejvlabel = tkinter.Label(self.optionFrame, text=str(self.edge_jump),anchor='w')
                    ejvlabel.grid(column=1, row=4, columnspan=2, rowspan=1, sticky='E') 
                    
                    if self.ShowNormvar.get() == 1:
                        normdivisor = (lineposty-lineprey)[index]
                        ply = (ply-lineprey)/normdivisor
                        lineposty = (lineposty-lineprey)/normdivisor
                        
                        lineprey = np.zeros(len(lineprey))
                        
                        if self.Flatvar.get() == 1:
                            
                            flat_correction = 1-lineposty
                            if energy[index] > energy[-1]:
                                flat_correction[index:int(len(energy))] = 0
                            else:
                                flat_correction[0:index] = 0
                            
                            ply = ply + flat_correction
                            lineposty = np.ones(len(lineposty))
                            linepre.set_data(energy, lineprey)
                            linepost.set_data(energy, lineposty)
                            line.set_ydata(ply)
                            line.set_xdata(energy) 
                        
                            scat.set_data(energy[index], ply[index])

                            self.ax.set_ylim(-0.1, 0.1+np.max(ply))
                        else:
                            linepre.set_data(energy, lineprey)
                            linepost.set_data(energy, lineposty)
                            line.set_ydata(ply)
                            line.set_xdata(energy) 
                        
                            scat.set_data(energy[index], ply[index])

                            self.ax.set_ylim(-0.1, 0.1+np.max(ply))
                    else:
                        linepre.set_data(energy, lineprey)
                        linepost.set_data(energy, lineposty)
                        line.set_ydata(ply)
                        line.set_xdata(energy) 
                        
                        scat.set_data(energy[index], ply[index])
                        self.ax.set_ylim(np.min(np.concatenate((lineprey,lineposty,ply), axis=0))-0.1, 0.1+np.max(np.concatenate((lineprey,lineposty,ply), axis=0)))
                     
                    self.canvas.draw_idle()
                    
            if self.initialize == True:
                self.EncoderFrame.grid_remove()
                self.CalibFrame.grid_remove()
                
                energy, mu, mu_filtered = plot_data(0)
                filter_toggle = int(self.CalibFiltervar.get())
                if filter_toggle == 1:
                    ply = mu_filtered
                else:
                    ply = mu
                [line] = self.ax.plot(energy, ply) 
                self.pre1line = self.ax.axvline(x=np.min(energy)+5, ls='--', c='g')
                self.pre2line = self.ax.axvline(x=self.calibenervar.get()-30, ls='--', c='g')
                self.post1line = self.ax.axvline(x=self.calibenervar.get()+100, ls='--', c='r')
                self.post2line = self.ax.axvline(x=np.max(energy)-10, ls='--', c='r')    

                self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
                self.canvas.get_tk_widget().grid(column=0, columnspan=3, padx=(25, 0), row=2, rowspan=5, sticky='EW')
                self.toolbar.destroy()
                self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
                self.toolbar.update()
                self.ds_scroll_scale = tkinter.Scale(self.canvasFrame, orient='horizontal', from_=0, to=int(len(self.calibroots)-5), command=plot_data)
                self.ds_scroll_scale.grid(column=1, columnspan=2, row=8, stick='EW')
                self.bf_scroll_scale = tkinter.Scale(self.canvasFrame, orient='horizontal', from_=0.0005, to=0.995, digits = 5, resolution = 0.0005, command=plot_data)
                self.bf_scroll_scale.grid(column=1, columnspan=2, row=9, stick='EW')
                self.bf_scroll_scale.set(self.Wn)
                
                self.datasetlabel = tkinter.Label(self.canvasFrame, text='Dataset',anchor='w')
                self.datasetlabel.grid(column=0, row=8, columnspan=1, rowspan=1, sticky='S')
                self.Filterlabel = tkinter.Label(self.canvasFrame, text='Butterworth Filter',anchor='w')
                self.Filterlabel.grid(column=0, row=9, columnspan=1, rowspan=1, sticky='S')
                
                self.NormFrame.grid(row=5,column=7,rowspan=3,sticky='EW')
                
                self.pre1_scroll_scale = tkinter.Scale(self.NormFrame, orient='horizontal', from_=np.min(energy)-self.calibenervar.get()+5, to=0, digits=4, resolution=0.1, length=150, command=plot_data)
                self.pre1_scroll_scale.grid(column=1, columnspan=1, row=0, stick='EW')
                self.pre1_scroll_scale.set(np.min(energy)-self.calibenervar.get())
                self.pre2_scroll_scale = tkinter.Scale(self.NormFrame, orient='horizontal', from_=np.min(energy)-self.calibenervar.get()+5, to=0, digits=4, resolution=0.1, length=150, command=plot_data)
                self.pre2_scroll_scale.grid(column=1, columnspan=1, row=1, stick='EW')
                self.pre2_scroll_scale.set(-15)
                self.post1_scroll_scale = tkinter.Scale(self.NormFrame, orient='horizontal', from_=0, to=np.max(energy)-self.calibenervar.get()-10, digits=4, resolution=0.1, length=150, command=plot_data)
                self.post1_scroll_scale.grid(column=1, columnspan=1, row=2, stick='EW')
                self.post1_scroll_scale.set(100)
                self.post2_scroll_scale = tkinter.Scale(self.NormFrame, orient='horizontal', from_=0, to=np.max(energy)-self.calibenervar.get()-10, digits=5, resolution=0.1, length=150, command=plot_data)
                self.post2_scroll_scale.grid(column=1, columnspan=1, row=3, stick='EW')
                self.post2_scroll_scale.set(np.max(energy)-self.calibenervar.get()-10)
                
                pre1label = tkinter.Label(self.NormFrame, text='Pre1',anchor='w')
                pre1label.grid(column=0, row=0, columnspan=1, rowspan=1, sticky='S')
                pre2label = tkinter.Label(self.NormFrame, text='Pre2',anchor='w')
                pre2label.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='S')
                post1label = tkinter.Label(self.NormFrame, text='Post1',anchor='w')
                post1label.grid(column=0, row=2, columnspan=1, rowspan=1, sticky='S')
                post2label = tkinter.Label(self.NormFrame, text='Post2',anchor='w')
                post2label.grid(column=0, row=3, columnspan=1, rowspan=1, sticky='S')
                
                self.e0_scroll_scale = tkinter.Scale(self.NormFrame, orient='horizontal', from_=self.calibenervar.get()-20, to=self.calibenervar.get()+20, digits=6, resolution=0.1, length=150, command=plot_data)
                self.e0_scroll_scale.grid(column=1, columnspan=1, row=4, stick='EW')
                self.e0_scroll_scale.set(self.calibenervar.get())
                
                e0label = tkinter.Label(self.NormFrame, text='E0',anchor='w')
                e0label.grid(column=0, row=4, columnspan=1, rowspan=1, sticky='W')
                
                index = min(range(len(energy)), key=lambda i: abs(energy[i]-self.calibenervar.get()))
                [scat] = self.ax.plot(energy[index], ply[index], marker='x', linewidth=2, markersize=12)
        
                self.ShowNormvar = tkinter.DoubleVar()
                self.ShowNormvar.set(0)
                self.ShowNormCheck = tkinter.Checkbutton(self.NormFrame, variable=self.ShowNormvar)
                self.ShowNormCheck.grid(column=2, row=7)
            
                shownormlabel = tkinter.Label(self.NormFrame, text='Show Norm',anchor='w')
                shownormlabel.grid(column=0, row=7, columnspan=2, rowspan=1, sticky='E')

                [linepre] = self.ax.plot(energy,self.regression('pre-edge', np.asarray(energy), np.asarray(ply)), c='g')
                [linepost] = self.ax.plot(energy,self.regression('post-edge', np.asarray(energy), np.asarray(ply)), c='r')
                self.initialize = False
            
                accept_button_ax = self.fig.add_axes([0.7, 0.025, 0.1, 0.04])
                accept_button = Button(accept_button_ax, 'Accept', color=self.axis_color, hovercolor='0.975')
                def accept_button_on_clicked(mouse_event):
                    NormValues = [self.pre1_scroll_scale.get(), self.pre2_scroll_scale.get(),self.post1_scroll_scale.get(),self.post2_scroll_scale.get(),self.norm_pre.get(),self.norm_post.get(),self.e0_scroll_scale.get()]
                    pdnormvalue = pd.DataFrame(NormValues)
                    pdnormvalue.to_csv(self.folder+'/Normalisation.dat', sep='\t', header=None, index=False)
                    self.close_figure()
                    self.clear_canvasFrame()
                    self.NormFrame.grid_remove()
                    self.EncoderFrame.grid(row=1,column=7,rowspan=3,sticky='EW')
                    self.CalibFrame.grid(row=5,column=7,rowspan=3,sticky='EW')
                    self.NormUsevar.set(1)
                    self.onselect(None)

                accept_button.on_clicked(accept_button_on_clicked)

                reject_button_ax = self.fig.add_axes([0.8, 0.025, 0.1, 0.04])
                reject_button = Button(reject_button_ax, 'Return', color=self.axis_color, hovercolor='0.975')
                def reject_button_on_clicked(mouse_event):
                    self.close_figure()
                    self.clear_canvasFrame()
                    self.NormFrame.grid_remove()
                    self.EncoderFrame.grid(row=2,column=7,rowspan=3,sticky='EW')
                    self.CalibFrame.grid(row=5,column=7,rowspan=3,sticky='EW')
                    self.onselect(None)

                reject_button.on_clicked(reject_button_on_clicked)
        
                self.canvas.draw()
                self.update()
                app.mainloop()
        
################################################################################	
################################################################################	
################################################################################	
    def process(self):
        print('Processing')
        i = 0
        print("started!")
        self.progress_bar["value"] = 0
        self.progress_bar["style"] = "blue.Horizontal.TProgressbar"
        
        self.rootsfile = self.rootsfile_check(self.folder)
        self.calibrated = self.calibfile_check(self.folder)
        if self.NormUsevar.get() == 1:
            self.normalized = self.normfile_check(self.folder)
        
        if (self.rootsfile == True) and (self.calibrated == True):
            if self.NormUsevar.get() == 1:
                if not os.path.exists(self.folder+'/Extract_Norm_'+str(int(self.togglevar.get()))):
                    os.makedirs(self.folder+'/Extract_Norm_'+str(int(self.togglevar.get())))
                export_path = self.folder+'/Extract_Norm_'+str(int(self.togglevar.get()))
            else:
                if not os.path.exists(self.folder+'/Extract_Raw_'+str(int(self.togglevar.get()))):
                    os.makedirs(self.folder+'/Extract_Raw_'+str(int(self.togglevar.get())))
                export_path = self.folder+'/Extract_Raw_'+str(int(self.togglevar.get()))
            try:
                self.calibroots = (pd.read_csv(self.folder+'/roots.dat', sep='\t', header=None)).iloc[:,0].values
                self.calibration_values = (pd.read_csv(self.folder+'/calibration.dat', sep='\t', header=None)).values
                for i in range(len(self.calibroots)-5):
                    ang, mu = self.spectrumread(self.calibroots[i+3],self.calibroots[2])
                    energy = np.transpose(1239.852/(float(self.calibration_values[3])*np.sin((ang+(self.calibration_values[1]-self.calibration_values[0]))*np.pi/180))).tolist()
                    if self.NormUsevar.get() == 1:     
                        if energy[0] > energy[-1]:
                            energy=energy[::-1]
                            mu=mu[::-1]

                        export = pd.DataFrame()
                        export['Ang'] = ang	
                        export['E'] = energy					
                        if int(self.Filtervar.get()) == 1:
                            N  = 3    # Filter order
                            B, A = signal.butter(N, float(self.calibration_values[2]), output='ba')
                            #C, D = signal.butter(N, [0.009, 0.2], btype='band', output='ba')
                        
                            mu_filtered = signal.filtfilt(B,A,mu)
                            #sine_filter = signal.filtfilt(C,D,mu_filtered)
                            mu_filtered = mu_filtered
                        else:
                            mu_filtered = mu   
                        mu = self.normalize(energy, mu_filtered)							
                        export[str(i)] = mu
                        self.progress_bar["value"] = 100*i/(len(self.calibroots)-5)
                        self.progress_bar.update()
                    
                        export.to_csv(export_path+'/'+self.data_bin.split('/')[-1]+'_'+str(i)+'.dat', sep='\t', index=False)
                self.progress_bar["value"] = 0
                print('Export complete')
            except:
                print('Batch Export Aborted')
                self.progress_bar["style"] = "red.Horizontal.TProgressbar"
                self.progress_bar["value"] = 100
                self.progress_bar.update()
        else:
            print('Batch Export Aborted')
            self.progress_bar["style"] = "red.Horizontal.TProgressbar"
            self.progress_bar["value"] = 100
            self.progress_bar.update()
            
        
################################################################################	
################################################################################	
################################################################################	
class DataProcess(tkinter.Frame):
    def __init__(self,name,*args,**kwargs):
        tkinter.Frame.__init__(self,*args,**kwargs)
################################################################################
        labelfile = tkinter.Label(self, text='Load Data',anchor='w')
        labelfile.grid(column=1, row=2, columnspan=1, rowspan=1, sticky='S')
		
        self.buttonplus= tkinter.Button(self, text='+', width=4, height=1,command=self.fileplus)
        self.buttonplus.grid(column=2, row=2, columnspan=1, sticky='S')

        self.buttonminus= tkinter.Button(self, text='-', width=4, height=1,command=self.fileminus)
        self.buttonminus.grid(column=3, row=2, columnspan=1, sticky='S') 

        self.listbox = tkinter.Listbox(self)
        self.listbox.grid(column=1, row=3, columnspan = 3, rowspan=10)
        self.listbox['width'] = 30
        self.listbox['height'] = 40
        self.listbox.bind('<<ListboxSelect>>', self.onselect)
        
        self.progress_bar = ttk.Progressbar(self, orient='horizontal', mode='determinate', maximum=100)
        self.progress_bar.grid(column=1, row=99, columnspan=99)
        self.progress_bar['length'] = 1000
        
        self.labelText = tkinter.StringVar()
        self.labelText.set('Time estimate')
        self.labeltime = tkinter.Label(self, textvariable=self.labelText)
        self.labeltime.grid(column=4, row=98, columnspan=3, rowspan=1, sticky='W')
		
        self.tpstextvar = tkinter.StringVar()
        self.tpstextvar.set('Time Per Spectrum')
        self.labeltps = tkinter.Label(self, textvariable=self.tpstextvar)
        self.labeltps.grid(column=7, row=98, columnspan=4, rowspan=1, sticky='W')
        
        self.files_long = []
        self.edge_info = []
        self.edge_type = []
        self.edge_jump = 0
################################################################################

        self.fig = plt.figure.Figure(figsize=(7, 7), facecolor='0.94')
        
        self.canvasFrame = tkinter.Frame(master=self, width=300)
        self.canvasFrame.grid(row=3,column=4,columnspan=3,rowspan=5,padx=(25, 25),sticky='EW')
		
        self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
        self.canvas.get_tk_widget().grid(column=0, columnspan=3, padx=(25, 0), row=2, rowspan=5, sticky='EW')
        
        self.ax=self.fig.add_subplot(111)
        
        x,y=[1,2,3],[4,5,6]
        self.ax.plot(x, y)
        self.ax.set_xlabel('Energy (eV)')
        self.ax.set_ylabel('Absorption')
        self.canvas.draw_idle()
        
        self.toolbarFrame = tkinter.Frame(master=self, width=300)
        self.toolbarFrame.grid(row=8,column=4,columnspan=3,padx=(25, 25),pady=(10,10),sticky='W')
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
        self.toolbar.update()

        self.axis_color = 'lightgoldenrodyellow'	
        
################################################################################ 
        self.plotcolsFrame = tkinter.Frame(self)
        self.plotcolsFrame.grid(row=1, column=7, rowspan=4, columnspan=3)
        
        plotcolslabel = tkinter.Label(self.plotcolsFrame, text='Plotting Columns',anchor='w')
        plotcolslabel.grid(column=0, row=0, columnspan=3, rowspan=1, sticky='N')
        
        numer = tkinter.Label(self.plotcolsFrame, text='Numerator',anchor='w')
        numer.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='W')
        
        denom = tkinter.Label(self.plotcolsFrame, text='Denomonator',anchor='w')
        denom.grid(column=0, row=2, columnspan=1, rowspan=1, sticky='W')
        
        self.sam_numer = tkinter.StringVar(self.plotcolsFrame)
        self.sam_numer.set('CH 0')
        self.sam_numerMenu = tkinter.OptionMenu(self.plotcolsFrame, self.sam_numer, 'CH 0','CH 1','CH 2','CH 3','1')
        self.sam_numerMenu.grid(column=1, row=1, columnspan=1) 
        
        self.sam_denom = tkinter.StringVar(self.plotcolsFrame)
        self.sam_denom.set('CH 1')
        self.sam_denomMenu = tkinter.OptionMenu(self.plotcolsFrame, self.sam_denom, 'CH 0','CH 1','CH 2','CH 3','1')
        self.sam_denomMenu.grid(column=1, row=2, columnspan=1) 
        
        self.ref_numer = tkinter.StringVar(self.plotcolsFrame)
        self.ref_numer.set('CH 1')
        self.ref_numerMenu = tkinter.OptionMenu(self.plotcolsFrame, self.ref_numer, 'CH 0','CH 1','CH 2','CH 3','1')
        self.ref_numerMenu.grid(column=2, row=1, columnspan=1) 
        
        self.ref_denom = tkinter.StringVar(self.plotcolsFrame)
        self.ref_denom.set('CH 2')
        self.ref_denomMenu = tkinter.OptionMenu(self.plotcolsFrame, self.ref_denom, 'CH 0','CH 1','CH 2','CH 3','1')
        self.ref_denomMenu.grid(column=2, row=2, columnspan=1) 
        
        log = tkinter.Label(self.plotcolsFrame, text='Log',anchor='w')
        log.grid(column=0, row=3, columnspan=1, rowspan=1, sticky='W') 
        
        self.sam_logvar = tkinter.DoubleVar()
        self.sam_logvar.set(1)
        self.sam_logCheck = tkinter.Checkbutton(self.plotcolsFrame, variable=self.sam_logvar)
        self.sam_logCheck.grid(column=1, row=3)
        
        self.ref_logvar = tkinter.DoubleVar()
        self.ref_logvar.set(1)
        self.ref_logCheck = tkinter.Checkbutton(self.plotcolsFrame, variable=self.ref_logvar)
        self.ref_logCheck.grid(column=2, row=3)
        
        filterlabel = tkinter.Label(self.plotcolsFrame, text='Filter',anchor='w')
        filterlabel.grid(column=0, row=6, columnspan=1, rowspan=1, sticky='W') 
        
        self.Filtervar = tkinter.DoubleVar()
        self.Filtervar.set(0)
        self.FilterCheck = tkinter.Checkbutton(self.plotcolsFrame, variable=self.Filtervar)
        self.FilterCheck.grid(column=2, row=6)
        
        process = tkinter.Label(self.plotcolsFrame, text='Process',anchor='w')
        process.grid(column=0, row=4, columnspan=1, rowspan=1, sticky='W') 
        
        self.sam_processvar = tkinter.DoubleVar()
        self.sam_processvar.set(1)
        self.sam_processCheck = tkinter.Checkbutton(self.plotcolsFrame, variable=self.sam_processvar)
        self.sam_processCheck.grid(column=1, row=4)
        
        self.ref_processvar = tkinter.DoubleVar()
        self.ref_processvar.set(1)
        self.ref_processCheck = tkinter.Checkbutton(self.plotcolsFrame, variable=self.ref_processvar)
        self.ref_processCheck.grid(column=2, row=4)
        
        self.togglevar = tkinter.IntVar()
        self.togglevar.set(0)
        
        
################################################################################          
        self.InterpFrame = tkinter.Frame(self, width=300)
        self.InterpFrame.grid(row=4,column=7,rowspan=5,sticky='EW')
        
        Emin = tkinter.Label(self.InterpFrame, text='Emin',anchor='w')
        Emin.grid(column=0, row=0, columnspan=1, rowspan=1, sticky='N')
        self.eminvar = tkinter.DoubleVar()
        self.eminentry = tkinter.Entry(self.InterpFrame, textvariable=self.eminvar)
        self.eminentry.grid(column=1, row=0, columnspan=2)   

        EmaxXANES = tkinter.Label(self.InterpFrame, text='Emax XANES',anchor='w')
        EmaxXANES.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='N')
        self.emaxxanesvar = tkinter.DoubleVar()
        self.emaxxanesentry = tkinter.Entry(self.InterpFrame, textvariable=self.emaxxanesvar)
        self.emaxxanesentry.grid(column=1, row=1, columnspan=2)   
        
        EmaxEXAFS = tkinter.Label(self.InterpFrame, text='Emax EXAFS',anchor='w')
        EmaxEXAFS.grid(column=0, row=2, columnspan=1, rowspan=1, sticky='N')
        self.emaxexafsvar = tkinter.DoubleVar()
        self.emaxexafsentry = tkinter.Entry(self.InterpFrame, textvariable=self.emaxexafsvar)
        self.emaxexafsentry.grid(column=1, row=2, columnspan=2)   
        
        Estep = tkinter.Label(self.InterpFrame, text='E Step',anchor='w')
        Estep.grid(column=0, row=3, columnspan=1, rowspan=1, sticky='N')
        self.estepvar = tkinter.DoubleVar()
        self.estepvar.set(0.25)
        self.estepentry = tkinter.Entry(self.InterpFrame, textvariable=self.estepvar)
        self.estepentry.grid(column=1, row=3, columnspan=2)
        
        constkuselabel = tkinter.Label(self.InterpFrame, text='Use Constant K',anchor='w')
        constkuselabel.grid(column=0, row=4, columnspan=2, rowspan=1, sticky='E') 
        self.constkUsevar = tkinter.DoubleVar()
        self.constkUsevar.set(0)
        self.constkUseCheck = tkinter.Checkbutton(self.InterpFrame, variable=self.constkUsevar)
        self.constkUseCheck.grid(column=2, row=4)
        
        Kstep = tkinter.Label(self.InterpFrame, text='K Step',anchor='w')
        Kstep.grid(column=0, row=5, columnspan=1, rowspan=1, sticky='N')
        self.kstepvar = tkinter.DoubleVar()
        self.kstepvar.set(0.025)
        self.kstepentry = tkinter.Entry(self.InterpFrame, textvariable=self.kstepvar)
        self.kstepentry.grid(column=1, row=5, columnspan=2)
        
        firstcaliblabel = tkinter.Label(self.InterpFrame, text='First Calibration',anchor='w')
        firstcaliblabel.grid(column=0, row=6, columnspan=2, rowspan=1, sticky='E') 
        self.firstcalibvar = tkinter.IntVar()
        self.firstcalibvar.set(0)
        self.firstcalibCheck = tkinter.Checkbutton(self.InterpFrame, variable=self.firstcalibvar)
        self.firstcalibCheck.grid(column=2, row=6)
        
        firstcalib = tkinter.Label(self.InterpFrame, text='Edge Step',anchor='w')
        firstcalib.grid(column=0, row=7, columnspan=1, rowspan=1, sticky='N')
        self.edgestepvar = tkinter.DoubleVar()
        self.edgestepentry = tkinter.Entry(self.InterpFrame, textvariable=self.edgestepvar)
        self.edgestepentry.grid(column=1, row=7, columnspan=2)
        
        toggle = [('',0),('',1)]
        self.radio_updown_buttons = []
        self.updownvar = tkinter.IntVar()
        self.updownvar.set(0)
        
        for idx, (text, mode) in enumerate(toggle):
            self.radio_updown_buttons.append(tkinter.Radiobutton(self.InterpFrame, text=text, variable=self.updownvar, value=mode))
            self.radio_updown_buttons[-1].grid(row=9, column=idx+1) 
        
        directionlabel = tkinter.Label(self.InterpFrame, text='Direction',anchor='w')
        directionlabel.grid(column=0, row=9, columnspan=1, rowspan=1, sticky='N')
        uplabel = tkinter.Label(self.InterpFrame, text='Up',anchor='w')
        uplabel.grid(column=1, row=8, columnspan=1, rowspan=1, sticky='N')
        downlabel = tkinter.Label(self.InterpFrame, text='Down',anchor='w')
        downlabel.grid(column=2, row=8, columnspan=1, rowspan=1, sticky='N')
        
        self.buttonexport= tkinter.Button(self.InterpFrame, text='Batch Export', width=30, height=1,command=self.batch_export)
        self.buttonexport.grid(column=0, row=10, columnspan=3, rowspan=1)
		
        #self.buttonabort= tkinter.Button(self.InterpFrame, text='Abort Export', width=30, height=1,command=self.abortexport)
        #self.buttonabort.grid(column=0, row=11, columnspan=3, rowspan=1)
        
        Advancedlabel = tkinter.Label(self.InterpFrame, text='Advanced Parameters',anchor='w')
        Advancedlabel.grid(column=0, row=12, columnspan=3, rowspan=1, sticky='N')
        
        #ThreadNum = tkinter.Label(self.InterpFrame, text='Threads',anchor='w')
        #ThreadNum.grid(column=0, row=13, columnspan=1, rowspan=1, sticky='N')
        #self.ThreadNumvar = tkinter.IntVar()
        #self.ThreadNumvar.set(psutil.cpu_count()-2)
        #self.ThreadNumentry = tkinter.Entry(self.InterpFrame, textvariable=self.ThreadNumvar)
        #self.ThreadNumentry.grid(column=1, row=13, columnspan=2)
        
        CPUNum = tkinter.Label(self.InterpFrame, text='Threads',anchor='w')
        CPUNum.grid(column=0, row=14, columnspan=1, rowspan=1, sticky='N')
        self.CPUNumvar = tkinter.IntVar()
        self.CPUNumvar.set(psutil.cpu_count()-2)
        self.CPUNumentry = tkinter.Entry(self.InterpFrame, textvariable=self.CPUNumvar)
        self.CPUNumentry.grid(column=1, row=14, columnspan=2)
        
        self.windowev = 1
        self.Sav_Gol = True
        self.abort = False
        self.active = False
        self.active_val = 0
        
################################################################################    
    def rootsfile_check(self, folder):
        try:
            self.calibroots = (pd.read_csv(folder+'/roots.dat', sep='\t', header=None)).values
            return True
        except:
            print('Encoder Analysis Required')
            return False
        
################################################################################
    def calibfile_check(self, folder):
        try:
            self.calibration_values = (pd.read_csv(folder+'/calibration.dat', sep='\t', header=None)).values
            return True
        except:
            print('Cannot Find Calibration Parameter File')
            return False
        
################################################################################
    def normfile_check(self, folder):
        try:
            self.normalisation_values = (pd.read_csv(folder+'/Normalisation.dat', sep='\t', header=None)).values
            return True
        except:
            print('Cannot Find Normalisation Parameter File')
            return False

################################################################################        
    def fileplus(self):
        self.askloadfile()

################################################################################        
    def abortexport(self):
        self.abort = True
        
################################################################################
    def askloadfile(self):
        self.edge_submitted = False
        file = askopenfile()
        try:
            dirname, file_short = os.path.split(os.path.abspath(file.name))
            #file_short = (str(file).split("'")[1]).split("/")[-1]
            print('Loading File')
            self.edge_box(file)
        except:
            print('Aborted Load')
        
################################################################################      
    def edge_box(self, file):
        self.box = tkinter.Toplevel()
        label0 =tkinter.Label(self.box, text='Enter Absorbing Atom')
        label0.pack()
	
        self.box.entry0 = tkinter.Entry(self.box)
        self.box.entry0.pack()
	
        label1 = tkinter.Label(self.box, text='Enter Edge type')
        label1.pack()
	
        self.box.entry1 = tkinter.Entry(self.box)
        self.box.entry1.pack()
	
        button1 = tkinter.Button(self.box, text='Submit', command=lambda:self.submit_edge(self.box, file))
        button1.pack()
	
        button2 = tkinter.Button(self.box, text='Cancel', command=lambda:self.box.destroy())
        button2.pack()
 
################################################################################          
    def submit_edge(self, box, file):
        self.central_atom = box.entry0.get()
        self.edge_type.append(box.entry0.get()+'_'+box.entry1.get())
        
        xdb = XrayDB()
        self.edge_info.append(xdb.xray_edge(box.entry0.get(), box.entry1.get()))
        print(self.edge_info[-1])
        self.core_hole = xdb.corehole_width(box.entry0.get(), box.entry1.get())
        print('core hole lifetime',self.core_hole)
        
        self.files_long.append(file)
        self.edge_submitted = True
        dirname, file_short = os.path.split(os.path.abspath(file.name))
        self.listbox.insert(END, file_short)
            
        box.destroy()
        self.update()

################################################################################
    def fileminus(self):        
        selection = self.listbox.curselection()
        for i in reversed(selection):
            self.listbox.delete(i)
            del self.files_long[i]
            del self.edge_info[i]
            del self.edge_type[i]
        
        self.update()
            
################################################################################
    def onselect(self,evt):
        self.batch = False
        try:
            data_bin_long = self.files_long[self.listbox.curselection()[0]].name
            self.active = True
            self.active_val = self.listbox.curselection()[0]
        except:
            if self.active == False:
                print('No File Selected')
        if self.active == True:
            data_bin_long = self.files_long[self.active_val].name
        
        if data_bin_long.split('.')[-1] == 'bin':
            self.data_bin = os.path.splitext(data_bin_long)[0]
            self.edge = self.edge_type[self.listbox.curselection()[0]]
            self.folder = 'output_'+self.edge+'_'+self.data_bin.split('/')[-1]
            
            self.rootsfile = DataRead(tkinter.Frame).rootsfile_check(self.folder)
            self.calibrated = DataRead(tkinter.Frame).calibfile_check(self.folder)
            self.normalized = DataRead(tkinter.Frame).normfile_check(self.folder)
            
            f = open(self.data_bin+'.bin', 'rb')
            f.seek(0)

            data = np.fromfile(f, dtype=np.int32, count = 2)
            self.headerSize = data[0]
            data = np.fromfile(f, dtype=np.int64, count = 1)
            data = np.fromfile(f, dtype='f4', count = 2)
            self.nChannels = data[0]
            
            self.choices = []
            for i in range(int(self.nChannels)):
                self.choices.append('CH '+str(int(i)))
            self.choices.append('1')
            
            self.sam_numerMenu.destroy()
            self.sam_denomMenu.destroy()
            self.ref_numerMenu.destroy()
            self.ref_denomMenu.destroy()
            
            self.sam_numer = tkinter.StringVar(self.plotcolsFrame)
            self.sam_numer.set('CH 0')
            self.sam_numerMenu = tkinter.OptionMenu(self.plotcolsFrame, self.sam_numer, *self.choices)
            self.sam_numerMenu.grid(column=1, row=1, columnspan=1) 
        
            self.sam_denom = tkinter.StringVar(self.plotcolsFrame)
            self.sam_denom.set('CH 1')
            self.sam_denomMenu = tkinter.OptionMenu(self.plotcolsFrame, self.sam_denom, *self.choices)
            self.sam_denomMenu.grid(column=1, row=2, columnspan=1) 
        
            self.ref_numer = tkinter.StringVar(self.plotcolsFrame)
            self.ref_numer.set('CH 1')
            self.ref_numerMenu = tkinter.OptionMenu(self.plotcolsFrame, self.ref_numer, *self.choices)
            self.ref_numerMenu.grid(column=2, row=1, columnspan=1) 
        
            self.ref_denom = tkinter.StringVar(self.plotcolsFrame)
            self.ref_denom.set('CH 2')
            self.ref_denomMenu = tkinter.OptionMenu(self.plotcolsFrame, self.ref_denom, *self.choices)
            self.ref_denomMenu.grid(column=2, row=2, columnspan=1)
            
            if (self.rootsfile == True) and (self.calibrated == True) and (self.normalized == True):
                self.calibroots = (pd.read_csv(self.folder+'/roots.dat', sep='\t', header=None)).values
                self.calibration_values = (pd.read_csv(self.folder+'/calibration.dat', sep='\t', header=None)).values
                self.normalisation_values = (pd.read_csv(self.folder+'/Normalisation.dat', sep='\t', header=None)).values
                
                self.ds_scroll_scale = tkinter.Scale(self.canvasFrame, orient='horizontal', from_=0, to=int(len(self.calibroots)-6), command=self.data_extract)
                self.ds_scroll_scale.grid(column=1, columnspan=2, row=8, stick='EW')
                self.bf_scroll_scale = tkinter.Scale(self.canvasFrame, orient='horizontal', from_=0.0005, to=0.9955, digits = 5, resolution = 0.0005, command=self.data_extract)
                self.bf_scroll_scale.grid(column=1, columnspan=2, row=9, stick='EW')
                self.bf_scroll_scale.set(float(self.calibration_values[2]))
                
                self.datasetlabel = tkinter.Label(self.canvasFrame, text='Dataset',anchor='w')
                self.datasetlabel.grid(column=0, row=8, columnspan=1, rowspan=1, sticky='S')
                self.Filterlabel = tkinter.Label(self.canvasFrame, text='Butterworth Filter',anchor='w')
                self.Filterlabel.grid(column=0, row=9, columnspan=1, rowspan=1, sticky='S')
                self.update()
            
################################################################################                 
    def regression(self,region,x,y,calib_energy): 
        
        if region == 'pre-edge':
            xllim = calib_energy + self.pre1_init
            xhlim = calib_energy + self.pre2_init
            try:
                order = self.norm_pre_init[0]
            except:
                order = str(int(self.norm_pre_init))
            
        if region == 'post-edge':
            xllim = calib_energy + self.post1_init
            xhlim = calib_energy + self.post2_init
            try:
                order = self.norm_post_init[0]
            except:
                order = str(int(self.norm_post_init))
                
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
    def interpolate(self,x,y,xnew):
        windowloc=200
        x = np.asarray(x)
        y = np.asarray(y)
	
        if self.Sav_Gol == True:
            if self.batch == False:
                idxSG = npi.indices(x, x[x <= self.emaxxanesvar.get()])
            else:
                idxSG = npi.indices(x, x[x <= self.batch_emaxxanesvar])
            y_xanes = y[0:np.max(idxSG)]
            y_exafs = y[np.max(idxSG):,]
            y_xsg = savgol_filter(y_xanes, self.window, self.poly)
            if len(y_exafs) > self.window*3:
                y_esg = savgol_filter(y_exafs, 3*self.window, self.poly)
            elif len(y_exafs) > self.window:
                y_esg = savgol_filter(y_exafs, self.window, self.poly)
            else:
                y_esg - y_exafs
		
        y = np.concatenate((y_xsg, y_esg), axis=0)
        
        startinterp = time.time()
	
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
	
        #pid = os.getpid()
        #py = psutil.Process(pid)
        #memoryUse = py.memory_info()[0]/2.**30  # memory use in MB...I think
        #print('memory use interp (MB):', 1000*(memoryUse-memoryUsei))
        if self.batch == False:
            print('total interp time:',time.time()-startinterp)
        return [xnew, ynew]	
    
################################################################################	
    def xnew_gen(self):
        if self.constkUsevar.get() == 1:
            print('Processing Exafs region with constant K step')
            k_min = round(np.sqrt((self.emaxxanesvar.get() - float(self.normalisation_values[6][0]))*0.262468426033256),1)
            k_max = round(np.sqrt((self.emaxexafsvar.get() - float(self.normalisation_values[6][0]))*0.262468426033256),1)
		
            k_grid = np.arange(k_min, k_max, self.kstepvar.get())
		
            e_grid = (((k_grid)**2)/0.262468426033256) + float(self.normalisation_values[6][0])
            #e_grid = (((k_grid)**2)/0.262468426033256) + 8979
            idxSG = npi.indices(e_grid, e_grid[e_grid >= self.emaxxanesvar.get()])
            e_grid = e_grid[np.min(idxSG):np.max(idxSG)]
			
            e_xanes = np.arange(self.eminvar.get(), self.emaxxanesvar.get()+self.estepvar.get(), self.estepvar.get())
            xnew = np.concatenate((e_xanes, e_grid), axis=0)
        else:
            print('Using constant Estep throughout')
            xnew = np.arange(self.eminvar.get(), self.emaxexafsvar.get(), self.estepvar.get())
        return xnew
    
################################################################################
    def data_extract(self, evt):
        self.initialize = True
        
        if self.batch == False:
            if (self.rootsfile == True) and (self.calibrated == True) and (self.normalized == True):
                xnew = self.xnew_gen()
                i = int(self.ds_scroll_scale.get())
                if self.sam_processvar.get() == 1:
                    self.togglevar.set(0)
                    ang_sam, mu_sam = self.spectrumread(self.calibroots[i+3],self.calibroots[2])
                    energy_sam = np.transpose(1239.852/(float(self.calibration_values[3])*np.sin((ang_sam+(self.calibration_values[1]-self.calibration_values[0]))*np.pi/180))).tolist()
                    if energy_sam[-1] < energy_sam[0]:
                        energy_sam = energy_sam[::-1]
                        mu_sam = mu_sam[::-1]              
                
                if self.ref_processvar.get() == 1:
                    self.togglevar.set(1)
                    ang_ref, mu_ref = self.spectrumread(self.calibroots[i+3],self.calibroots[2])
                    energy_ref = np.transpose(1239.852/(float(self.calibration_values[3])*np.sin((ang_ref+(self.calibration_values[1]-self.calibration_values[0]))*np.pi/180))).tolist()
                    if energy_ref[-1] < energy_ref[0]:
                        energy_ref = energy_ref[::-1]
                        mu_ref = mu_ref[::-1]
            
                if self.Filtervar.get() == 1:
                    Wn = float(self.bf_scroll_scale.get()) # Cutoff frequency
                    N  = 3    # Filter order
                    B, A = signal.butter(N, Wn, output='ba')
                    if self.sam_processvar.get() == 1:
                        mu_sam_filtered = signal.filtfilt(B,A,mu_sam)
                    if self.ref_processvar.get() == 1:
                        mu_ref_filtered = signal.filtfilt(B,A,mu_ref)
                else:
                    if self.sam_processvar.get() == 1:
                        mu_sam_filtered = mu_sam
                    if self.ref_processvar.get() == 1:
                        mu_ref_filtered = mu_ref
                
                pre1,pre2,post1,post2,order_pre,order_post,e0 = self.normalisation_values
                self.e0 = float(e0[0])
            
                if self.Sav_Gol == True:
                    print('Savitzky-Golay window in eV ='+str(self.windowev))
                    if self.sam_processvar.get() == 1:
                        idxWin = npi.indices(np.asarray(energy_sam), np.asarray(energy_sam)[(np.asarray(energy_sam) >= self.e0) & (np.asarray(energy_sam) <= (float(self.e0) + self.windowev))])
                    else:
                        idxWin = npi.indices(np.asarray(energy_ref), np.asarray(energy_ref)[(np.asarray(energy_ref) >= float(self.e0)) & (np.asarray(energy_ref) <= (float(self.e0) + self.windowev))])
                    self.window = len(idxWin)
                    print('Savitzky-Golay window in datapoints ='+str(self.window))
                    #makes sure the window is odd in length
                    if (self.window % 2) == 0:
                        self.window = self.window - 1
                    if self.window < 3:
                        self.window = 3
                    self.poly = 2
            
                if self.sam_processvar.get() == 1:
                    mu_sam_norm = self.normalize(energy_sam, mu_sam)
                    xinterp_sam, yinterp_sam = self.interpolate(energy_sam, mu_sam_norm, xnew)
                if self.ref_processvar.get() == 1:
                    mu_ref_norm = self.normalize(energy_ref, mu_ref_filtered)
                    xinterp_ref, yinterp_ref = self.interpolate(energy_ref, mu_ref_norm, xnew)

                if self.initialize == True:
                    self.ax.clear()
                    self.fig = plt.figure.Figure(figsize=(7, 7), facecolor='0.94')
                    self.ax=self.fig.add_subplot(111)
                    
                    self.ax.set_xlabel('Energy (eV)')
                    self.ax.set_ylabel('Absorption')
            
                    if self.sam_processvar.get() == 1:
                        [line1]  = self.ax.plot(energy_sam, mu_sam_norm)
                        [line2] = self.ax.plot(xinterp_sam, yinterp_sam)
                    if self.ref_processvar.get() == 1:
                        [line3] = self.ax.plot(xinterp_ref, yinterp_ref)
                
                    self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
                    self.canvas.get_tk_widget().grid(column=0, columnspan=3, padx=(25, 0), row=2, rowspan=5, sticky='EW')
                    self.canvas.draw_idle()
        
                    self.toolbar.destroy()
                    self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
                    self.toolbar.update()
                    
                    self.canvas.draw()
                    
                else:
                    if self.sam_processvar.get() == 1:
                        line1.set_data(energy_sam, mu_sam)
                        line2.set_data(xinterp_sam, yinterp_sam)
                    if self.ref_processvar.get() == 1:
                        line3.set_data(xinterp_ref, yinterp_ref)
                    self.canvas.draw_idle()                    
                    
            else:
                print('Missing required preprocessing files')
                        
  ################################################################################  
    def spectrumread(self,rooti,minr):
        encoder_bin = self.data_bin+'_Encoder'
        
        minr = int(minr)
            
        g = open(self.data_bin+'.bin', 'rb')
        
        g.seek(0)

        data = np.fromfile(g, dtype=np.int32, count = 2)

        self.headerSize = data[0]
        #print("HeaderSize: ", self.headerSize)
        #print("FileVersion: ", fileVersion)

        data = np.fromfile(g, dtype=np.int64, count = 1)
        #print("DateTime: ", DT)

        data = np.fromfile(g, dtype='f4', count = 2)

        self.nChannels = data[0]
        self.AdcClock_Hz = data[1]
        self.DacClock = 1

        g.seek(int(self.headerSize+(4*self.nChannels*rooti)))
        data_D = (np.fromfile(g, dtype='f4', count = int(self.nChannels*minr)))
        #print(len(data_D))
        #print(len(data_D)/self.nChannels)
        #print(len(data_D[0::int(self.nChannels)]))
        #print(len(data_D[1::int(self.nChannels)]))
        #print(len(data_D[2::int(self.nChannels)]))
        mu_r = np.zeros((minr, 1))
        
        if self.batch == False:
            toggle = int(self.togglevar.get())
            
            if toggle == 0:
                numervar = str(self.sam_numer.get())
                denomvar = str(self.sam_denom.get())
                logvar = int(self.sam_logvar.get())
            else:
                numervar = str(self.ref_numer.get())
                denomvar = str(self.ref_denom.get())
                logvar = int(self.ref_logvar.get())
        else:
            if self.batch_toggle == 0:
                numervar = self.batch_sam_numer
                denomvar = self.batch_sam_denom    
                logvar = self.batch_sam_logvar
            else:
                numervar = self.batch_ref_numer
                denomvar = self.batch_ref_denom
                logvar = self.batch_ref_logvar
        
        for j in range(int(self.nChannels)+1):
            if numervar == self.choices[j]:
                numercol = int(j)
            if denomvar == self.choices[j]:
                denomcol = int(j)

        if logvar == 1:
            if denomcol == int(len(self.choices)-1):
                mu_r[:,0] = np.log(data_D[numercol::int(self.nChannels)])
            else:
                mu_r[:,0] = np.log(data_D[numercol::int(self.nChannels)]/data_D[denomcol::int(self.nChannels)])
        else:
            if denomcol == int(len(self.choices)-1):
                mu_r[:,0] = data_D[numercol::int(self.nChannels)]
            else:
                mu_r[:,0] = data_D[numercol::int(self.nChannels)]/data_D[denomcol::int(self.nChannels)]
		
        nans = np.argwhere(np.isnan(mu_r[:,0]))
        mu = np.zeros((minr-len(nans), 1))
        mu[:,0] = np.delete(mu_r[:,0], nans)

        f = open(encoder_bin+'.bin', 'rb')
        f.seek(int(self.headerSize+(4*rooti)))
        ang_r = np.fromfile(f, dtype='f4', count = int(minr))
        ang = np.delete(ang_r, nans)

        RawData = pd.DataFrame()
        RawData['ang'] = ang
        RawData['mu'] = mu
        
        self.angf = ang[0]
        self.angl = ang[-1]
        
        RawData=RawData[(RawData['ang'] >= self.calibroots[0][0]) & (RawData['ang'] <= self.calibroots[1][0])]  
        RawData = RawData.groupby('ang', as_index=False).mean()
        
        return RawData['ang'].values, RawData['mu'].values
    
################################################################################        
    def normalize(self,energy, mu):
        pre1,pre2,post1,post2,order_pre,order_post,e0 = self.normalisation_values
        
        self.pre1_init = float(pre1[0])
        self.pre2_init = float(pre2[0])
        self.post1_init = float(post1[0])
        self.post2_init = float(post2[0])
        self.norm_pre_init = order_pre[0]
        self.norm_post_init = order_post[0]
        self.e0 = float(e0[0])
        
        self.Flatvar = 1
        
        lineprey = self.regression('pre-edge',np.asarray(energy), np.asarray(mu), self.e0)        
        lineposty = self.regression('post-edge',np.asarray(energy), np.asarray(mu), self.e0)
        
        index = min(range(len(energy)), key=lambda i: abs(energy[i]-self.e0))
        normdivisor = (lineposty-lineprey)[index]
        self.edge_jump_norm = normdivisor

        mu = (mu-lineprey)/normdivisor
        lineposty = (lineposty-lineprey)/normdivisor
        lineprey = np.zeros(len(lineprey))
                        
        if self.Flatvar == 1:
            flat_correction = 1-lineposty
            if energy[index] > energy[-1]:
                flat_correction[index:int(len(energy))] = 0
            else:
                flat_correction[0:index] = 0
            mu = mu + flat_correction
            return mu
        else:
            return mu
			
################################################################################
    def cumgauss(self, x, mu, sigma):
        return 0.5 * (1 + special.erf((x-mu)/(np.sqrt(2)*sigma)))
        
################################################################################   
    def batch_export(self):
        self.abort = False
        self.batch = True
        start_batch = time.time()
        sam_process = self.sam_processvar.get()
        ref_process = self.ref_processvar.get()
        self.batch_sam_numer = str(self.sam_numer.get())
        self.batch_sam_denom = str(self.sam_denom.get())
        self.batch_sam_logvar = int(self.sam_logvar.get())
        self.batch_ref_numer = str(self.ref_numer.get())
        self.batch_ref_denom = str(self.ref_denom.get())
        self.batch_ref_logvar = int(self.ref_logvar.get())
        self.batch_updownvar = self.updownvar.get()
        self.batch_filtervar = self.Filtervar.get()
        self.batch_firstcalibvar = self.firstcalibvar.get()
        self.batch_bf_scroll_scale = self.bf_scroll_scale.get()
        self.batch_emaxxanesvar = self.emaxxanesvar.get()
        self.batch_edgestepvar = self.edgestepvar.get()
        
        def batch_extract(i):
            try:
                if (self.rootsfile == True) and (self.calibrated == True) and (self.normalized == True):
                    #i = int(self.ds_scroll_scale.get())
                    if sam_process == 1:
                        self.batch_toggle = 0
                        ang_sam, mu_sam = self.spectrumread(self.calibroots[i+3],self.calibroots[2])
                        energy_sam = np.transpose(1239.852/(float(self.calibration_values[3])*np.sin((ang_sam+(self.calibration_values[1]-self.calibration_values[0]))*np.pi/180))).tolist()
                        if energy_sam[-1] < energy_sam[0]:
                            energy_sam = energy_sam[::-1]
                            mu_sam = mu_sam[::-1]
                    
                    if ref_process == 1:
                        self.batch_toggle = 1
                        ang_ref, mu_ref = self.spectrumread(self.calibroots[i+3],self.calibroots[2])
                        energy_ref = np.transpose(1239.852/(float(self.calibration_values[3])*np.sin((ang_ref+(self.calibration_values[1]-self.calibration_values[0]))*np.pi/180))).tolist()
                        if energy_ref[-1] < energy_ref[0]:
                            mu_ref = mu_ref[::-1]
                            energy_ref = energy_ref[::-1]
                    if self.angl < self.angf:
                        direction = 0
                    else:
                        direction = 1
                    
                    if direction == self.batch_updownvar:
                        if self.batch_filtervar == 1:
                            Wn = float(self.batch_bf_scroll_scale) # Cutoff frequency
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
                        pre1,pre2,post1,post2,order_pre,order_post,e0 = self.normalisation_values
                        self.e0 = float(e0[0])
                
                        if self.Sav_Gol == True:
                            #print('Savitzky-Golay window in eV ='+str(self.windowev))
                            if sam_process == 1:
                                idxWin = npi.indices(np.asarray(energy_sam), np.asarray(energy_sam)[(np.asarray(energy_sam) >= self.e0) & (np.asarray(energy_sam) <= (self.e0 + self.windowev))])
                            else:
                                idxWin = npi.indices(np.asarray(energy_ref), np.asarray(energy_ref)[(np.asarray(energy_ref) >= self.e0) & (np.asarray(energy_ref) <= (self.e0 + self.windowev))])
                            self.window = len(idxWin)
                            #print('Savitzky-Golay window in datapoints ='+str(self.window))
                            #makes sure the window is odd in length
                            if (self.window % 2) == 0:
                                self.window = self.window - 1
                            if self.window < 3:
                                self.window = 3
                            self.poly = 2
                        
                        if self.normalized == True:
                            if sam_process == 1:
                                mu_sam = self.normalize(energy_sam, mu_sam)
                                mu_sam_norm = self.normalize(energy_sam, mu_sam_filtered)
                                self.edge_jump = self.edge_jump_norm
                            if ref_process == 1:
                                mu_ref_norm = self.normalize(energy_ref, mu_ref_filtered)
                        
                                #fits step function to the reference edge step to allow for recalibration of all spectra
                                #step_mod = StepModel(form='erf', prefix='step_')
                                #pars = step_mod.make_params(sigma=2, amplitude=1, center=float(self.e0))
                                #mod = step_mod
                                #out = mod.fit(mu_ref_norm, pars, x=energy_ref)
												
                                popt, pcov = curve_fit(self.cumgauss, energy_ref, mu_ref_norm, (self.e0, self.core_hole))
                                print('Optimised Parameters', popt)
                    
                        if (i == 0) or (i==1):
                            if ref_process == 1:
                                if self.batch_firstcalibvar == 1:
                                    self.batch_edgestepvar = popt[0]
                                    self.batch_firstcalibvar = 0
                                    ec = 0
                                else:
                                    ec = self.batch_edgestepvar - popt[0]
                            else:
                                self.batch_edgestepvar = 0
                                ec = 0
                        else:
                            if ref_process == 1:
                                ec = self.batch_edgestepvar - popt[0]
                                print('step postion = ',popt[0])
                            else:
                                ec = 0
                                self.batch_edgestepvar = 0
                        if ref_process == 1:
                            print('Energy Correction ',ec)
                        
                        if sam_process == 1:
                            energy_sam = np.asarray(energy_sam) + ec
                            xinterp_sam, yinterp_sam = self.interpolate(energy_sam, mu_sam_norm, self.xnew)
                        if ref_process == 1:
                            energy_ref = np.asarray(energy_sam)
                            xinterp_ref, yinterp_ref = self.interpolate(energy_ref, mu_ref_norm, self.xnew)
                        
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
                return([0,0],[0,0],[0,0],[0,0],[0,0],[0,0])
            
        if not os.path.exists(self.folder+'/Interp_Norm/Individual'):
            os.makedirs(self.folder+'/Interp_Norm/Individual')
        export_path = self.folder+'/Interp_Norm'

        b=np.asarray(range(int(len(self.calibroots)-5)))
        self.std = time.time()
        if self.sam_processvar.get() == 1:
            global df_tot_sam
            df_tot_sam = [None] * len(b)
        if self.ref_processvar.get() == 1:
            global df_tot_ref
            df_tot_ref = [None] * len(b)
        
        self.xnew = self.xnew_gen()
        df_xnew = pd.DataFrame()
        df_xnew['xnew'] = self.xnew
        df_xnew.to_csv(self.folder+'/xnew_info.dat', sep='\t', header=False, index=False)
        
        skip = True
        
        job_init=0
        while skip == True:
            x1,y1,x2,y2,x3,y3 = batch_extract(job_init)
            if (len(x1) == 1) and (len(x2) == 1):
                skip = True
                b = b[1::]
            elif (len(x1) == 2) and (len(x2) == 2):
                skip = True
                b = b[1::]
            else:
                if self.sam_processvar.get() == 1: 
                    df_tot_sam[0] = y1
                if self.ref_processvar.get() == 1:
                    df_tot_ref[0] = y2
                skip = False
					
                self.ax.clear()
                self.fig = plt.figure.Figure(figsize=(7, 7), facecolor='0.94')
                self.ax=self.fig.add_subplot(111)
            
                if self.sam_processvar.get() == 1:
                    [line1]  = self.ax.plot(x3, y3)
                    [line2] = self.ax.plot(x1, y1)
                if self.ref_processvar.get() == 1:
                    [line3] = self.ax.plot(x2, y2)
						
                self.ax.set_xlabel('Energy (eV)')
                self.ax.set_ylabel('Absorption')
            
                self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
                self.canvas.get_tk_widget().grid(column=0, columnspan=3, padx=(25, 0), row=2, rowspan=5, sticky='EW')
                self.canvas.draw_idle()
        
                self.toolbar.destroy()
                self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
                self.toolbar.update()
            
                self.canvas.draw()
                self.update()
                
            job_init = job_init + 1
                      
        num_cpus = self.CPUNumvar.get()
        print('num cpus = ',num_cpus)
        print('job length = ',int(len(self.calibroots)-5))
        self.labelText.set('Starting Processess')
        self.update()
        time.sleep(0.1)
        options = [int(len(self.calibroots)-5), self.folder, num_cpus, sam_process, ref_process, export_path, self.data_bin, self.batch_sam_numer, self.batch_sam_denom, self.batch_sam_logvar, self.batch_ref_numer, self.batch_ref_denom ,self.batch_ref_logvar, self.batch_updownvar, self.batch_filtervar, self.batch_edgestepvar, self.batch_bf_scroll_scale, self.batch_emaxxanesvar, self.Sav_Gol, self.windowev, self.core_hole]       
        sub_f = 'batch_extract_subroutine_v49.py'
        start_loop = time.time()
        process = subprocess.Popen(['python', sub_f, str(options)], stdout=subprocess.PIPE, shell=True)
        
        while process.poll() == None:
            stdoutdata = process.stdout.readline()
            if 'Affinity' in stdoutdata.decode('utf-8'):
                print(stdoutdata.decode('utf-8').split('\n')[0])
            if 'Percent' in stdoutdata.decode('utf-8'):
                perc = float(stdoutdata.decode('utf-8').split('\n')[0].split(' ')[-1])
                self.progress_bar["value"] = perc
                self.progress_bar.update()
            if 'Estimated' in stdoutdata.decode('utf-8'):
                self.labelText.set(stdoutdata.decode('utf-8').split('\n')[0])
            if 'Spectrum' in stdoutdata.decode('utf-8'):
                self.tpstextvar.set(stdoutdata.decode('utf-8').split('\n')[0])
            if 'Returning' in stdoutdata.decode('utf-8'):
                self.progress_bar["value"] = 100
                self.labelText.set('Finishing Up')
                self.progress_bar.update()
            if 'Returned' in stdoutdata.decode('utf-8'):
                self.progress_bar["value"] = 100
                self.labelText.set('Forming Matrices')
                self.progress_bar.update()
                break
            #print(stdoutdata.decode('utf-8').split('\n')[0])

        end_loop = time.time()
        print('Total Subprocess Time = ', end_loop - start_loop)			
        i=0
        none_vals = []
		
        self.tpstextvar.set(' ')
        self.update()
        
        if sam_process == 1:
            df_sam = pd.DataFrame()
        if ref_process == 1:
            df_ref = pd.DataFrame()
        
        for i in range(int(len(self.calibroots)-5)):
            try:
                data_read = pd.read_csv(self.folder+'/Interp_Norm/Individual/'+str(int(i))+'.dat', sep = '\t', header=0)
                if sam_process == 1:
                    df_sam[str(i)] = data_read['mu_sam']
                if ref_process == 1:
                    df_ref[str(i)] = data_read['mu_ref']
                i = i+1
            except:
                none_vals.append(i)
                i = i+1
            self.progress_bar["value"] = 100*(i+1)/(len(self.calibroots)-5)
            self.progress_bar.update() 
            self.update()
         
        self.labelText.set('Saving Matrices')
        self.update()		 
		 
        try:
            if sam_process == 1:
                df_sam.insert(0, 'E', data_read['E'])
                df_sam.to_csv(export_path+'/'+self.data_bin.split('/')[-1]+'_sam_matrix.dat', sep='\t', index=False)
                print('sample matrix saved')
                del df_sam
                del df_tot_sam
        except:
            print('failed to save matrix')
        
        try:            
            if ref_process == 1:
                df_ref.insert(0, 'E', data_read['E'])
                df_ref.to_csv(export_path+'/'+self.data_bin.split('/')[-1]+'_ref_matrix.dat', sep='\t', index=False)
                print('reference matrix saved')
                del df_ref
                del df_tot_ref
        except:
            print('failed to save matrix')
                
        gc.collect()
        print('Export Complete')
        end_batch = time.time()
        print('Total Process Time = ', end_batch - start_batch)
        self.labelText.set('Time estimate')
        self.progress_bar["value"] = 0
        self.progress_bar.update()

        self.edgestepvar.set(self.batch_edgestepvar)
        self.firstcalibvar.set(0)

        main_f = str(os.path.basename(main.__file__).split('.')[0])
        
        outf = open(export_path+'/parameters.dat','w')
        outf.write(str(main_f)+'\n')
        outf.write('file name '+str(self.data_bin)+'\n')
        outf.write('folder '+str(self.folder)+'\n')
        outf.write('\n')
        outf.write('Calibration Energy '+str(self.e0)+'\n')
        outf.write('Normalisation Parameters'+'\n')
        outf.write('Edge Jump '+str(self.edge_jump)+'\n')
        outf.write('Pre-edge 1 '+str(self.pre1_init)+'\n')
        outf.write('Pre-edge 2 '+str(self.pre2_init)+'\n')
        outf.write('order pre-edge '+str(self.norm_pre_init)+'\n')
        outf.write('Post-edge 1 '+str(self.post1_init)+'\n')
        outf.write('Post-edge 2 '+str(self.post2_init)+'\n')
        outf.write('order post-edge '+str(self.norm_post_init)+'\n')
        outf.write('\n')
        #outf.write('Savitzky-Golay window '+str(self.windowev)+'\n')
        #outf.write('Savitsky-Golay polynomial '+str(self.poly)+'\n')
        outf.write('Edge-step '+str(self.edgestepvar.get())+'\n')
        outf.write('Estep '+str(self.estepvar.get())+'\n')
        if self.constkUsevar.get() == 1:	
            outf.write('Kstep '+str(self.kstepvar.get())+'\n')
        outf.write('Emin '+str(self.eminvar.get())+'\n')
        outf.write('Emax XANES '+str(self.emaxxanesvar.get())+'\n')
        outf.write('Emax EXAFS '+str(self.emaxexafsvar.get())+'\n')
        outf.write('\n')
        outf.write(str(self.sam_processvar.get())+' '+str(self.ref_processvar.get())+'\n')
        outf.write('Butterworth Filter value '+str(self.calibration_values[2][0]))
        outf.close()
        
        self.batch = False

################################################################################	
################################################################################	
################################################################################	
class PostProcess(tkinter.Frame):
    def __init__(self,name,*args,**kwargs):
        self.frame_dataread = tkinter.Frame.__init__(self,*args,**kwargs)
################################################################################

        labelfile = tkinter.Label(self, text='Load Data',anchor='w')
        labelfile.grid(column=1, row=2, columnspan=1, rowspan=1, sticky='S')
		
        self.buttonplus= tkinter.Button(self, text='+', width=4, height=1,command=self.fileplus)
        self.buttonplus.grid(column=2, row=2, columnspan=1, sticky='S')

        self.buttonminus= tkinter.Button(self, text='-', width=4, height=1,command=self.fileminus)
        self.buttonminus.grid(column=3, row=2, columnspan=1, sticky='S') 

        self.listbox = tkinter.Listbox(self)
        self.listbox.grid(column=1, row=3, columnspan = 3, rowspan=10)
        self.listbox['width'] = 30
        self.listbox['height'] = 40
        self.listbox.bind('<<ListboxSelect>>', self.onselect)
        
        self.progress_bar = ttk.Progressbar(self, orient='horizontal', mode='determinate', maximum=100, style="red.Horizontal.TProgressbar")
        self.progress_bar.grid(column=1, row=99, columnspan=99)
        self.progress_bar['length'] = 1000
        
        self.files_long = []
        self.initialize = True

################################################################################

        self.fig = plt.figure.Figure(figsize=(7, 7), facecolor='0.94')
        
        self.canvasFrame = tkinter.Frame(master=self, width=300)
        self.canvasFrame.grid(row=3,column=4,columnspan=3,rowspan=5,padx=(25, 0),sticky='EW')
		
        self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
        self.canvas.get_tk_widget().grid(column=0, columnspan=3, padx=(0, 0), row=2, rowspan=5, sticky='EW')
        
        self.ax=self.fig.add_subplot(111)
        
        x,y=[1,2,3],[4,5,6]
        self.ax.plot(x, y)
        self.ax.set_xlabel('Energy (eV)')
        self.ax.set_ylabel('Absorption')
        self.canvas.draw_idle()
        
        self.toolbarFrame = tkinter.Frame(master=self, width=300)
        self.toolbarFrame.grid(row=8,column=4,columnspan=3,padx=(25, 25),pady=(10,10),sticky='W')
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
        self.toolbar.update()

        self.axis_color = 'lightgoldenrodyellow'	
        
################################################################################          
        self.AverageFrame = tkinter.Frame(self, width=300)
        self.AverageFrame.grid(row=1,column=7,rowspan=6,sticky='EW', padx=25)
        
        NSUM = tkinter.Label(self.AverageFrame, text='Nsum',anchor='w')
        NSUM.grid(column=0, row=1, columnspan=2, rowspan=1, sticky='N')
        self.nsumvar = tkinter.StringVar()
        self.nsumvar.set('5')
        self.nsumentry = tkinter.Entry(self.AverageFrame, textvariable=self.nsumvar)
        self.nsumentry.grid(column=2, row=1, columnspan=2)   
        
        self.buttonaverage= tkinter.Button(self.AverageFrame, text='Average', width=30, height=1,command=self.average)
        self.buttonaverage.grid(column=0, row=6, columnspan=3, rowspan=1)
        
################################################################################          
        self.MatrixCutFrame = tkinter.Frame(self, width=300)
        self.MatrixCutFrame.grid(row=3,column=7,rowspan=5,sticky='EW', padx=25)
        
        Elow = tkinter.Label(self.AverageFrame, text='Emin',anchor='w')
        Elow.grid(column=0, row=3, columnspan=2, rowspan=1, sticky='N')
        self.Elowvar = tkinter.DoubleVar()
        self.Elowentry = tkinter.Entry(self.AverageFrame, textvariable=self.Elowvar)
        self.Elowentry.grid(column=2, row=3, columnspan=2)   
        
        Ehigh = tkinter.Label(self.AverageFrame, text='Emax',anchor='w')
        Ehigh.grid(column=0, row=4, columnspan=2, rowspan=1, sticky='N')
        self.Ehighvar = tkinter.DoubleVar()
        self.Ehighentry = tkinter.Entry(self.AverageFrame, textvariable=self.Ehighvar)
        self.Ehighentry.grid(column=2, row=4, columnspan=2)   
        
        self.buttoncut= tkinter.Button(self.MatrixCutFrame, text='Matrix Slice', width=30, height=1,command=self.matrix_cut)
        self.buttoncut.grid(column=0, row=4, columnspan=3, rowspan=1)
        
        mcr = tkinter.Label(self.AverageFrame, text='MCR Output',anchor='w')
        mcr.grid(column=0, row=5, columnspan=1, rowspan=1, sticky='W') 
        
        self.mcrvar = tkinter.DoubleVar()
        self.mcrvar.set(0)
        self.mcrCheck = tkinter.Checkbutton(self.AverageFrame, variable=self.mcrvar)
        self.mcrCheck.grid(column=2, row=5)
          
################################################################################
    def fileminus(self):        
        i = self.listbox.curselection()[0]
        self.listbox.delete(i)
        del self.files_long[i]
        self.update()
################################################################################        
    def fileplus(self):
        self.askloadfile()
        
################################################################################
    def askloadfile(self):
        self.edge_submitted = False
        file = askopenfile()
        try:
            self.files_long.append(file)
            dirname, file_short = os.path.split(os.path.abspath(file.name))
            #file_short = (str(file).split("'")[1]).split("/")[-1]
            print('Loading File')
            self.listbox.insert(END, file_short)
        except:
            print('Aborted Load')
################################################################################
    def onselect(self,evt):
        self.initialize = True
        print(self.files_long[self.listbox.curselection()[0]].name)
        data_long = self.files_long[self.listbox.curselection()[0]].name
        dirname, filename = os.path.split(os.path.abspath(data_long))
        self.folder = dirname
        self.file = os.path.splitext(filename)[0]
        
        def plot(evt):
            i = self.ds_scroll_scale.get()
            line1.set_ydata(self.data.iloc[:,(i+1)])
            self.lowline.set_xdata(self.Ehighvar.get())
            self.highline.set_xdata(self.Elowvar.get())
            self.canvas.draw_idle()
        
        if data_long.split('.')[-1] == 'dat':
            self.data = pd.read_csv(data_long, sep='\t', header=0)
            cols = self.data.columns.values.tolist()
            
            self.Elowvar.set(min(self.data['E']))
            self.Ehighvar.set(max(self.data['E']))

            self.col_max = cols[-1]
			
            print(self.col_max)
            
            if self.initialize == True:
                self.ax.clear()
                self.fig = plt.figure.Figure(figsize=(7, 7), facecolor='0.94')
                self.ax=self.fig.add_subplot(111)
                    
                self.ax.set_xlabel('Energy (eV)')
                self.ax.set_ylabel('Absorption')
            
                [line1]  = self.ax.plot(self.data['E'], self.data.iloc[:,1])
                [line2] = self.ax.plot(self.data['E'], self.data[str(int(self.col_max))])
                
                self.lowline = self.ax.axvline(x=self.Elowvar.get(), ls='--', c='g')
                self.highline = self.ax.axvline(x=self.Ehighvar.get(), ls='--', c='g')
                
                self.datasetlabel = tkinter.Label(self.canvasFrame, text='Dataset',anchor='w')
                self.datasetlabel.grid(column=0, row=8, columnspan=1, rowspan=1, sticky='S')
                self.ds_scroll_scale = tkinter.Scale(self.canvasFrame, orient='horizontal', from_=0, to=int(self.data.shape[1]-2), command=plot)
                self.ds_scroll_scale.grid(column=1, columnspan=2, row=8, stick='EW')
                self.ds_scroll_scale.set(0)

                self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
                self.canvas.get_tk_widget().grid(column=0, columnspan=3, padx=(25, 0), row=2, rowspan=5, sticky='EW')
                self.canvas.draw_idle()
        
                self.toolbar.destroy()
                self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
                self.toolbar.update()
                    
                self.canvas.draw()
                self.initialize = False
################################################################################	 
    def clear_canvasFrame(self):
        list = self.canvasFrame.grid_slaves()
        for child in list:
            if str(type(child)) != "<class 'tkinter.Canvas'>":
                child.destroy()
                
################################################################################
    def close_figure(self):
        self.fig.clear()
        self.fig = plt.figure.Figure(figsize=(7, 7), facecolor='0.94')
        self.ax=self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
        self.canvas.get_tk_widget().grid(column=0, columnspan=3, padx=(25, 0), row=2, rowspan=5, sticky='EW')
        self.toolbar.destroy()
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
        self.toolbar.update()
                
################################################################################
    def average(self):
        self.initialize = True
        self.data_sum = pd.DataFrame()
        cols = self.data.columns.values.tolist()
        cols = np.float_(cols[1:])
        if self.nsumvar.get() == 'all':
            Nsum = int(np.max(cols))
        else:
            Nsum = int(self.nsumvar.get())
						
        Nave = int(np.floor((np.max(cols))/Nsum))
        print('Number of Averaged Spectra =',Nave)
        idxs = list(map(int, cols))
        for j in range(int(Nave)):
            idx = npi.indices(np.asarray(idxs), np.asarray(idxs)[(np.asarray(idxs) > (j*Nsum)) & (np.asarray(idxs) <= (j*Nsum)+Nsum)])
            if len(idx) > 0:
                self.data_sum[str(j)] = self.data.iloc[:,idx[0]+1]
                for i in range(len(idx)-1):
                    self.data_sum[str(j)] = self.data.iloc[:,idx[i+1]+1] + self.data_sum[str(j)]
                self.data_sum[str(j)] = self.data_sum[str(j)]/len(idx)
				
            self.progress_bar["value"] = 100*(j+1)/Nave
            self.progress_bar.update()	
        self.data_sum.insert(0, 'E', self.data['E'])				
				
        self.progress_bar["value"] = 0
        self.progress_bar.update()       
        def plot_average(evt):
            i = self.ds_scroll_scale.get()
            line2.set_ydata(self.data_sum.iloc[:,i+1])
            self.canvas.draw_idle()
        
        if self.initialize == True:
            self.ax.clear()
            self.fig = plt.figure.Figure(figsize=(7, 7), facecolor='0.94')
            self.ax=self.fig.add_subplot(111)
                    
            self.ax.set_xlabel('Energy (eV)')
            self.ax.set_ylabel('Absorption')
            
            [line1] = self.ax.plot(self.data['E'], self.data.iloc[:,self.data.shape[1]-1])
            [line2] = self.ax.plot(self.data_sum['E'], self.data_sum.iloc[:,1])
			
            self.datasetlabel = tkinter.Label(self.canvasFrame, text='Dataset',anchor='w')
            self.datasetlabel.grid(column=0, row=8, columnspan=1, rowspan=1, sticky='S')
            self.ds_scroll_scale = tkinter.Scale(self.canvasFrame, orient='horizontal', from_=0, to=int(self.data_sum.shape[1]-2), command=plot_average)
            self.ds_scroll_scale.grid(column=1, columnspan=2, row=8, stick='EW')
            self.ds_scroll_scale.set(0)

            self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
            self.canvas.get_tk_widget().grid(column=0, columnspan=3, padx=(25, 0), row=2, rowspan=5, sticky='EW')
            self.canvas.draw_idle()
            
            self.toolbar.destroy()
            self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
            self.toolbar.update()
                    
            self.initialize = False
            
            accept_button_ax = self.fig.add_axes([0.7, 0.025, 0.1, 0.04])
            accept_button = Button(accept_button_ax, 'Accept', color=self.axis_color, hovercolor='0.975')
            def accept_button_on_clicked(mouse_event):
                self.data_sum.to_csv(self.folder+'/'+self.file+'_average_'+str(Nsum)+'.dat', sep='\t', index=False)
                print('Averaged Data Matrix Saved')
                if self.Ehighvar.get() < max(self.data['E']):
                    self.data_sum = self.data_sum[(self.data_sum.iloc[:,0] >= self.Elowvar.get()) & (self.data_sum.iloc[:,0] <= self.Ehighvar.get())]
                    self.data_sum.to_csv(self.folder+'/'+self.file+'_average_'+str(Nsum)+'_cut.dat', sep='\t', index=False)
                if self.mcrvar.get() == 1:
                    self.data_mcr = self.data_sum
                    self.data_mcr.drop(columns=self.data_mcr.columns.values.tolist()[0], inplace=True)
                    self.data_mcr = self.data_mcr.T
                    self.data_mcr.to_csv(self.folder+'/'+self.file+'_average_'+str(Nsum)+'_mcr.dat', sep='\t', index=False, header=None)
                    print('MCR file created')
    
            accept_button.on_clicked(accept_button_on_clicked)
            
            close_button_ax = self.fig.add_axes([0.8, 0.025, 0.1, 0.04])
            close_button = Button(close_button_ax, 'Close', color=self.axis_color, hovercolor='0.975')
            def close_button_on_clicked(mouse_event):
                    self.close_figure()
                    self.clear_canvasFrame()
                    self.onselect(None)
    
            close_button.on_clicked(close_button_on_clicked)

            self.canvas.draw()
            self.update()
            app.mainloop()
################################################################################
    def matrix_cut(self):
        
        self.initialize = True
        
        def plot_update(evt):
            i = self.ds_scroll_scale.get()
            line1.set_ydata(self.data.iloc[:,i+1])
            self.canvas.draw_idle()
        
        if self.initialize== True:
            self.initialize = False
            self.data = self.data[(self.data.iloc[:,0] >= self.Elowvar.get()) & (self.data.iloc[:,0] <= self.Ehighvar.get())]
            self.data.to_csv(self.folder+'/'+self.file+'_cut.dat', sep='\t', index=False)
            print('Matrix Sliced')
        
            self.ax.clear()
            
            [line1] = self.ax.plot(self.data['E'], self.data.iloc[:,1])
            self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
            self.canvas.get_tk_widget().grid(column=0, columnspan=3, padx=(25, 0), row=2, rowspan=5, sticky='EW')
            self.canvas.draw_idle()
            
            self.toolbar.destroy()
            self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
            self.toolbar.update()
			
            if self.data.shape[1] > 2:
                [line2] = self.ax.plot(self.data['E'], self.data.iloc[:,self.data.shape[1]-1])
			
            self.datasetlabel = tkinter.Label(self.canvasFrame, text='Dataset',anchor='w')
            self.datasetlabel.grid(column=0, row=8, columnspan=1, rowspan=1, sticky='S')
            self.ds_scroll_scale = tkinter.Scale(self.canvasFrame, orient='horizontal', from_=0, to=int(self.data.shape[1]-2), command=plot_update)
            self.ds_scroll_scale.grid(column=1, columnspan=2, row=8, stick='EW')
            self.ds_scroll_scale.set(0)
            if self.mcrvar.get() == 1:
                self.data_mcr = self.data
                self.data_mcr.drop(columns=self.data_mcr.columns.values.tolist()[0])
                self.data_mcr = self.data_mcr.T
                self.data_mcr.to_csv(self.folder+'/'+self.file+'_mcr.dat', sep='\t', index=False, header=False)
                print('MCR file created')
            
            self.canvas.draw()
            self.update()
            app.mainloop()

################################################################################
if __name__ == "__main__":
    app = App()
    app.title('ProXAS Gui')
    app.mainloop()