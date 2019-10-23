import warnings
import warnings
warnings.filterwarnings("ignore",category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
import numpy as np
import pandas as pd
import matplotlib as plt
plt.use('WXAGG')
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import tkinter, time, os, psutil, subprocess, sys, shutil, ast, codecs, re, larch, gc, peakutils.peak, itertools
from tkinter.filedialog import askopenfile
from tkinter import ttk
from tkinter import END, MULTIPLE
from scipy.signal import savgol_filter
from matplotlib.widgets import Button
import scipy.signal as signal
from scipy.optimize import curve_fit
from scipy.interpolate import Rbf
import numpy_indexed as npi
from scipy import special
from scipy.optimize import curve_fit
from larch import Interpreter
from larch_plugins.xafs import xftf, ftwindow, pre_edge, autobk, xftr
from scipy.signal import savgol_filter
from matplotlib import gridspec
import matplotlib.ticker as ticker
from larch import ValidateLarchPlugin, parse_group_args
from larch_plugins.xafs import set_xafsGroup
from larch_plugins.xafs.cauchy_wavelet import cauchy_wavelet
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from scipy import stats
from scipy.signal import find_peaks

from pymcr.simplisma import svd, simplisma
from pymcr.mcr import McrAls
from pymcr.regressors import OLS, NNLS
from pymcr.constraints import ConstraintNonneg, ConstraintNorm
from pymcr.lcf import lcf

from proxas.interpolation import interpolate
from proxas.functions import polynomials
from proxas.headerread import headerread
from proxas.dataread import dataread
from proxas.xraydb import XrayDB
import __main__ as main

print('ProXAS GUI')
class App(tkinter.Tk):
	def __init__(self,*args,**kwargs):
		tkinter.Tk.__init__(self,*args,**kwargs)
		self.resizable(True, True)
		
		self.add_tab()

	def add_tab(self):
		App.notebook = ttk.Notebook()
		App.notebook.grid(row=0)
		tab = DataExtract(App.notebook)
		tab2 = PostProcess(App.notebook)
		tab3 = FourierTransform(App.notebook)
		tab4 = Plotting(App.notebook)
		tab5 = Analysis(App.notebook)
		App.notebook.add(tab,text="Data Extract")
		App.notebook.add(tab2,text="Post Processing")
		App.notebook.add(tab3,text="Fourier Transform")
		App.notebook.add(tab4,text="Plotting")
		App.notebook.add(tab5, text="Analysis")
		
	def get_selected(evt):
		return App.notebook.tab(App.notebook.select(), "text")
		
################################################################################		
#DATAEXTRACTTAB	
################################################################################	
class DataExtract(tkinter.Frame):
	def __init__(self,name,*args,**kwargs):
		self.frame_dataread = tkinter.Frame.__init__(self,*args,**kwargs)
				
		self.screen_width = self.winfo_screenwidth()
		self.screen_height = self.winfo_screenheight()
		
		self.progress_bar = ttk.Progressbar(self, orient='horizontal', mode='determinate', maximum=100, style="red.Horizontal.TProgressbar")
		self.progress_bar.grid(column=0, row=99, columnspan=99)
		self.progress_bar['length'] = 1000
		
		self.labelText = tkinter.StringVar()
		self.labelText.set('Time estimate')
		self.labeltime = tkinter.Label(self, textvariable=self.labelText)
		self.labeltime.grid(column=1, row=98, columnspan=2, rowspan=1, sticky='W')
		
		self.tpstextvar = tkinter.StringVar()
		self.tpstextvar.set('Time Per Spectrum')
		self.labeltps = tkinter.Label(self, textvariable=self.tpstextvar)
		self.labeltps.grid(column=3, row=98, columnspan=1, rowspan=1, sticky='W')
		
		self.frame1 = tkinter.LabelFrame(self, width=0.25*0.8333*self.screen_height, height=0.75*self.screen_height)
		self.frame2 = tkinter.LabelFrame(self, width=0.7*0.8333*self.screen_height, height=0.7*0.75*self.screen_height)
		self.frame3 = tkinter.LabelFrame(self, width=0.7*0.8333*self.screen_height, height=0.3*0.75*self.screen_height, borderwidth = 0, highlightthickness = 0)
		self.frame4 = tkinter.LabelFrame(self, width=0.25*0.8333*self.screen_height, height=0.75*self.screen_height, borderwidth = 0, highlightthickness = 0)
		self.frame5 = tkinter.LabelFrame(self, width=0.7*0.8333*self.screen_height, height=0.3*0.75*self.screen_height, borderwidth = 0, highlightthickness = 0)
		self.frame1.grid(row=0, column=0, rowspan=4, columnspan=1, padx=2, sticky='N')
		self.frame2.grid(row=0, column=1, rowspan=2, columnspan=2, padx=2, pady=15)
		self.frame3.grid(row=2, column=1, rowspan=2, columnspan=1, padx=2, pady=15)
		self.frame4.grid(row=0, column=3, rowspan=4, columnspan=1, padx=2, pady=15, sticky='N')
		self.frame5.grid(row=2, column=2, rowspan=2, columnspan=1, padx=2, pady=15)
		
################################################################################
#FRAME1
################################################################################		
		labelfile = tkinter.Label(self.frame1, text='Load Data',anchor='w')
		labelfile.grid(column=0, row=0, columnspan=3, rowspan=1, sticky='S')
		
		self.buttonplus= tkinter.Button(self.frame1, text='+', width=4, height=1,command=self.fileplus)
		self.buttonplus.grid(column=0, row=1, columnspan=1, sticky='E')
		
		self.buttonminus= tkinter.Button(self.frame1, text='-', width=4, height=1,command=self.fileminus)
		self.buttonminus.grid(column=1, row=1, columnspan=1, sticky='W') 
		
		self.listbox = tkinter.Listbox(self.frame1)
		self.listbox.grid(column=0, row=2, columnspan = 2)
		self.listbox['width'] = 50
		self.listbox['height'] = 53
		self.listbox.bind('<<ListboxSelect>>', self.onselect)
		
		self.files_long = []
		self.edge_info = []
		self.edge_type = []
		self.edge_info_calib = []
		self.edge_type_calib = []
		self.edge_jump = 0
		self.import_type = []
		
################################################################################
#FRAME2
################################################################################
		def FigureCanvas(self):
			self.fig = plt.figure.Figure(figsize=(5.75, 5.75))
			
			self.canvasFrame = tkinter.Frame(master=self.frame2, padx=16, background="white")
			self.canvasFrame.grid(row=0,column=0)
			self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
			self.canvas.get_tk_widget().grid(column=0, row=0, sticky='EW')
			
			self.ax=self.fig.add_subplot(111)
			self.canvas.draw_idle()
			
			self.toolbarFrame = tkinter.Frame(master=self.frame2)
			self.toolbarFrame.grid(row=1,column=0, sticky='W')
			self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
			self.toolbar.update()
			
			self.axis_color = 'lightgoldenrodyellow'	
			self.update()
		
		FigureCanvas(self)

################################################################################
#FRAME3
################################################################################
		self.plotcolsFrame = tkinter.Frame(master=self.frame3)
		self.plotcolsFrame.grid(row=0, column=0, rowspan=1, columnspan=1)
		
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
		self.sam_logCheck = tkinter.Checkbutton(self.plotcolsFrame, variable=self.sam_logvar, command=self.plot_update_trigger)
		self.sam_logCheck.grid(column=1, row=3)
		
		self.ref_logvar = tkinter.DoubleVar()
		self.ref_logvar.set(1)
		self.ref_logCheck = tkinter.Checkbutton(self.plotcolsFrame, variable=self.ref_logvar, command=self.plot_update_trigger)
		self.ref_logCheck.grid(column=2, row=3)
		
		filterlabel = tkinter.Label(self.plotcolsFrame, text='Filter',anchor='w')
		filterlabel.grid(column=0, row=5, columnspan=1, rowspan=1, sticky='W') 
		
		self.Filtervar = tkinter.DoubleVar()
		self.Filtervar.set(1)
		self.FilterCheck = tkinter.Checkbutton(self.plotcolsFrame, variable=self.Filtervar, command=self.plot_update_trigger)
		self.FilterCheck.grid(column=2, row=5)
		
		plottoggle = tkinter.Label(self.plotcolsFrame, text='Plot',anchor='w')
		plottoggle.grid(column=0, row=4, columnspan=1, rowspan=1, sticky='W')
		
		toggle = [('',0),('',1)]
		self.radio_plot_buttons = []
		self.togglevar = tkinter.DoubleVar()
		self.togglevar.set(0)
		
		for idx, (text, mode) in enumerate(toggle):
			self.radio_plot_buttons.append(tkinter.Radiobutton(self.plotcolsFrame, text=text, variable=self.togglevar, value=mode, command=self.plot_update_trigger))
			self.radio_plot_buttons[-1].grid(row=4, column=idx+1) 
			
		self.plot_update_button=tkinter.Button(self.plotcolsFrame, text='Update Plot', command=self.plot_update_trigger).grid(column=0, row=6)
		
################################################################################
#FRAME4
################################################################################		
		def raise_frame(frame):
			frame.tkraise()
			
		self.f1 = tkinter.Frame(master=self.frame4)
		f2 = tkinter.Frame(master=self.frame4)
		
		for frame in (self.f1, f2):
			frame.grid(row=0, column=0, sticky='news')
		
		
################################################################################
#FRAME4-1
################################################################################
		self.ToProcessFrame = tkinter.LabelFrame(master=self.f1, width=300, height=0.75*self.screen_height, pady=21, padx=123)
		tkinter.Button(self.ToProcessFrame, text='Go to Data Export', command=lambda:raise_frame(f2)).grid(row=0, column=0, columnspan=3)

################################################################################
#Encoder Analysis Frame		
		self.EncoderFrame = tkinter.LabelFrame(master=self.f1, width=300, height=0.75*self.screen_height, pady=5, padx=5)
		self.EncoderFrame.grid(row=1, column=0)
			
		self.split_accept_button=tkinter.Button(self.EncoderFrame, text='Accept', state='disabled', width=14, command=self.accept_split).grid(column=0, row=5, stick='E')
		self.split_insert_button=tkinter.Button(self.EncoderFrame, text='Insert', state='disabled', width=14, command=self.insert_split).grid(column=1, row=5, stick='EW')
		self.split_delete_button=tkinter.Button(self.EncoderFrame, text='Delete', state='disabled', width=14, command=self.delete_split).grid(column=2, row=5, stick='W')
		
		self.buttonsplit= tkinter.Button(self.EncoderFrame, text='Analyze Encoder', width=30, height=1,command=self.split)
		self.buttonsplit.grid(column=0, row=0, columnspan=3, rowspan=1)
		
		refactor = tkinter.Label(self.EncoderFrame, text='Refactor',anchor='w')
		refactor.grid(column=0, row=1, columnspan=2, rowspan=1, sticky='N')
		self.splitrefacvar = tkinter.IntVar()
		self.splitrefacvar.set(100)
		self.splitrefactorentry = tkinter.Entry(self.EncoderFrame, textvariable=self.splitrefacvar)
		self.splitrefactorentry.grid(column=2, row=1, columnspan=1) 
		
		encoder_smooth_label = tkinter.Label(self.EncoderFrame, text='Encoder Smoothing',anchor='w')
		encoder_smooth_label.grid(column=0, row=2, columnspan=2, rowspan=1, sticky='W')
		
		self.encoder_smoothvar = tkinter.DoubleVar()
		self.encoder_smoothvar.set(51)
		self.encoder_smoothentry = tkinter.Entry(self.EncoderFrame, textvariable=self.encoder_smoothvar)
		self.encoder_smoothentry.grid(column=2, row=2, columnspan=1)    
		
		sampling = tkinter.Label(self.EncoderFrame, text='Sampling Frequency (MHz)',anchor='w')
		sampling.grid(column=0, row=3, columnspan=2, rowspan=1, sticky='W')
		
		self.samplingvar = tkinter.IntVar()
		self.samplingvar.set(2)
		self.samplingentry = tkinter.Entry(self.EncoderFrame, textvariable=self.samplingvar)
		self.samplingentry.grid(column=2, row=3, columnspan=1)   
		
################################################################################
#Calibration Frame	
		self.CalibrationFrame = tkinter.LabelFrame(master=self.f1, width=300, height=0.75*self.screen_height, pady=5, padx=28)
		
		self.buttoncalib= tkinter.Button(self.CalibrationFrame, text='Calibrate', width=30, height=1,command=self.calibrate)
		self.buttoncalib.grid(column=0, row=0, columnspan=3) 
		
		mono = tkinter.Label(self.CalibrationFrame, text='Mono Crystal',anchor='w')
		mono.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='N')
		
		self.calibmono = tkinter.StringVar(self.CalibrationFrame)
		self.calibmono.set('Si 111')
		self.calibmonoMenu = tkinter.OptionMenu(self.CalibrationFrame, self.calibmono, 'Si 111', 'Si 311')
		self.calibmonoMenu.grid(column=1, row=1, columnspan=1) 
		
		calibener = tkinter.Label(self.CalibrationFrame, text='Energy (eV)',anchor='w')
		calibener.grid(column=0, row=2, columnspan=1, rowspan=1, sticky='N')
		self.calibenervar = tkinter.DoubleVar()
		self.calibenerentry = tkinter.Entry(self.CalibrationFrame, textvariable=self.calibenervar)
		self.calibenerentry.grid(column=1, row=2, columnspan=2)   
		
		calibedge = tkinter.Label(self.CalibrationFrame, text='Edge',anchor='w')
		calibedge.grid(column=0, row=3, columnspan=1, rowspan=1, sticky='N')
		self.calibedgevar = tkinter.StringVar()
		self.calibedgeentry = tkinter.Entry(self.CalibrationFrame, textvariable=self.calibedgevar)
		self.calibedgeentry.grid(column=1, row=3, columnspan=2)
		
		calibfilterlabel = tkinter.Label(self.CalibrationFrame, text='Derivative Filter',anchor='w')
		calibfilterlabel.grid(column=0, row=4, columnspan=2, rowspan=1, sticky='E') 
		self.CalibFiltervar = tkinter.DoubleVar()
		self.CalibFiltervar.set(1)
		self.CalibFilterCheck = tkinter.Checkbutton(self.CalibrationFrame, variable=self.CalibFiltervar)
		self.CalibFilterCheck.grid(column=2, row=4)
		
		calibuselabel = tkinter.Label(self.CalibrationFrame, text='Use Calibration',anchor='w')
		calibuselabel.grid(column=0, row=5, columnspan=2, rowspan=1, sticky='E') 
		self.CalibUsevar = tkinter.DoubleVar()
		self.CalibUsevar.set(0)
		self.CalibUseCheck = tkinter.Checkbutton(self.CalibrationFrame, variable=self.CalibUsevar, command=self.plot_update_trigger)
		self.CalibUseCheck.grid(column=2, row=5)
		
		blackmanfilterlabel = tkinter.Label(self.CalibrationFrame, text='BlackmanHarris Filter',anchor='w')
		blackmanfilterlabel.grid(column=0, row=6, columnspan=2, rowspan=1, sticky='E') 
		self.BlackmanHarrisFiltervar = tkinter.DoubleVar()
		self.BlackmanHarrisFiltervar.set(1)
		self.BlackmanHarrisFiltervarCheck = tkinter.Checkbutton(self.CalibrationFrame, variable=self.BlackmanHarrisFiltervar)
		self.BlackmanHarrisFiltervarCheck.grid(column=2, row=6)
		
		blackmanfilterwindowlabel = tkinter.Label(self.CalibrationFrame, text='Width',anchor='w')
		blackmanfilterwindowlabel.grid(column=0, row=7, columnspan=1, rowspan=1, sticky='N')
		self.blackmanfilterwindowvar = tkinter.DoubleVar()
		self.blackmanfilterwindowvar.set(20)
		self.blackmanfilterwindowentry = tkinter.Entry(self.CalibrationFrame, textvariable=self.blackmanfilterwindowvar)
		self.blackmanfilterwindowentry.grid(column=1, row=7, columnspan=2)   
		
		self.twod = tkinter.DoubleVar()
		self.twod.set(0.627120)
		
		self.calib_accept_button=tkinter.Button(self.CalibrationFrame, text='Accept', state='disabled', width=14, command=self.accept_calib).grid(column=1, row=8, stick='E')
		self.calib_cancel_button=tkinter.Button(self.CalibrationFrame, text='Cancel', state='disabled', width=14, command=self.cancel_calib).grid(column=2, row=8, stick='W')
		
################################################################################
#Normalization Frame	
		self.NormalisationFrame = tkinter.LabelFrame(master=self.f1, width=300, height=0.75*self.screen_height, pady=5, padx=11)
		
		self.buttonnorm= tkinter.Button(self.NormalisationFrame, text='Normalise', width=30, height=1,command=self.normalize)
		self.buttonnorm.grid(column=0, row=1, columnspan=2)
        
		normuselabel = tkinter.Label(self.NormalisationFrame, text='Use Normalization',anchor='w')
		normuselabel.grid(column=0, row=2, columnspan=2, rowspan=1, sticky='E') 
		self.NormUsevar = tkinter.DoubleVar()
		self.NormUsevar.set(0)
		self.NormUseCheck = tkinter.Checkbutton(self.NormalisationFrame, variable=self.NormUsevar, command=self.plot_update_trigger)
		self.NormUseCheck.grid(column=2, row=2)
		
		Flatlabel = tkinter.Label(self.NormalisationFrame, text='Flatten',anchor='w')
		Flatlabel.grid(column=0, row=3, columnspan=2, rowspan=1, sticky='E') 
		self.Flatvar = tkinter.DoubleVar()
		self.Flatvar.set(1)
		self.FlatCheck = tkinter.Checkbutton(self.NormalisationFrame, variable=self.Flatvar, command=self.plot_update_trigger)
		self.FlatCheck.grid(column=2, row=3)
		
		ejlabel = tkinter.Label(self.NormalisationFrame, text='Edge Jump',anchor='w')
		ejlabel.grid(column=0, row=10, columnspan=1, rowspan=1, sticky='E')
		
		self.evjlabel = tkinter.StringVar()
		self.evjlabel.set('0')
		self.ejv = tkinter.Label(self.NormalisationFrame, textvariable=self.evjlabel,anchor='w')
		self.ejv.grid(column=1, row=10, columnspan=2, rowspan=1, sticky='E') 		
		
		self.buttonloadnorm= tkinter.Button(self.NormalisationFrame, text='Load Normalization', width=14, height=1,command=self.normload)
		self.buttonloadnorm.grid(column=2, row=1, columnspan=1)
		
		self.norm_pre = tkinter.StringVar(self.NormalisationFrame)
		self.norm_pre.set('1')
		self.norm_preMenu = tkinter.OptionMenu(self.NormalisationFrame, self.norm_pre, '0','1','2','3','4','V', command=self.plot_update_trigger2)
		self.norm_preMenu.grid(column=2, row=5, columnspan=1) 
        
		self.norm_post = tkinter.StringVar(self.NormalisationFrame)
		self.norm_post.set('V')
		self.norm_postMenu = tkinter.OptionMenu(self.NormalisationFrame, self.norm_post, '0','1','2','3','4','V', command=self.plot_update_trigger2)
		self.norm_postMenu.grid(column=2, row=7, columnspan=1) 
		
		self.pre1_scroll_scale = tkinter.Scale(self.NormalisationFrame, orient='horizontal', from_=-150, to=0, state='disabled', digits=4, resolution=0.1, length=125)
		self.pre1_scroll_scale.grid(column=1, columnspan=1, row=5, stick='EW')
		self.pre1_scroll_scale.set(-100)
		self.pre2_scroll_scale = tkinter.Scale(self.NormalisationFrame, orient='horizontal', from_=-150, to=0, state='disabled', digits=4, resolution=0.1, length=125)
		self.pre2_scroll_scale.grid(column=1, columnspan=1, row=6, stick='EW')
		self.pre2_scroll_scale.set(-50)
		self.post1_scroll_scale = tkinter.Scale(self.NormalisationFrame, orient='horizontal', from_=0, to=1000, state='disabled', digits=4, resolution=0.1, length=125)
		self.post1_scroll_scale.grid(column=1, columnspan=1, row=7, stick='EW')
		self.post1_scroll_scale.set(100)
		self.post2_scroll_scale = tkinter.Scale(self.NormalisationFrame, orient='horizontal', from_=0, to=1000, state='disabled', digits=5, resolution=0.1, length=125)
		self.post2_scroll_scale.grid(column=1, columnspan=1, row=8, stick='EW')
		self.post2_scroll_scale.set(1000)
		
		self.e0_scroll_scale = tkinter.DoubleVar()
		self.e0_entry = tkinter.Entry(self.NormalisationFrame, textvariable=self.e0_scroll_scale, state='normal')
		self.e0_entry.grid(column=1, columnspan=1, row=9, stick='EW')
		self.e0_scroll_scale.set(float(0))
		
		pre1label = tkinter.Label(self.NormalisationFrame, text='Pre1',anchor='w')
		pre1label.grid(column=0, row=5, columnspan=1, rowspan=1, sticky='S')
		pre2label = tkinter.Label(self.NormalisationFrame, text='Pre2',anchor='w')
		pre2label.grid(column=0, row=6, columnspan=1, rowspan=1, sticky='S')
		post1label = tkinter.Label(self.NormalisationFrame, text='Post1',anchor='w')
		post1label.grid(column=0, row=7, columnspan=1, rowspan=1, sticky='S')
		post2label = tkinter.Label(self.NormalisationFrame, text='Post2',anchor='w')
		post2label.grid(column=0, row=8, columnspan=1, rowspan=1, sticky='S')
		e0label = tkinter.Label(self.NormalisationFrame, text='E0',anchor='w')
		e0label.grid(column=0, row=9, columnspan=1, rowspan=1, sticky='S')
		
		self.norm_accept_button=tkinter.Button(self.NormalisationFrame, text='Accept', state='disabled', width=14, command=self.accept_norm).grid(column=1, row=11, stick='E')
		self.norm_cancel_button=tkinter.Button(self.NormalisationFrame, text='Cancel', state='disabled', width=14, command=self.cancel_norm).grid(column=2, row=11, stick='W')
		
		self.update_e0_button=tkinter.Button(self.NormalisationFrame, text='Update E0').grid(column=2, row=9)
		
################################################################################
#Interpolation Frame
		self.InterpolationFrame = tkinter.LabelFrame(master=self.f1, width=300, height=0.75*self.screen_height, pady=5, padx=64)
		
		Interpuselabel = tkinter.Label(self.InterpolationFrame, text='Use Interpolation',anchor='w')
		Interpuselabel.grid(column=0, row=0, columnspan=2, rowspan=1, sticky='E') 
		self.InterpUsevar = tkinter.DoubleVar()
		self.InterpUsevar.set(0)
		self.InterpUseCheck = tkinter.Checkbutton(self.InterpolationFrame, variable=self.InterpUsevar, command=self.plot_update_trigger)
		self.InterpUseCheck.grid(column=2, row=0)

		Emin = tkinter.Label(self.InterpolationFrame, text='Emin (eV)',anchor='w')
		Emin.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='N')
		self.eminvar = tkinter.DoubleVar()
		self.eminentry = tkinter.Entry(self.InterpolationFrame, textvariable=self.eminvar)
		self.eminentry.grid(column=1, row=1, columnspan=2)   
		
		EmaxXANES = tkinter.Label(self.InterpolationFrame, text='Emax XANES (eV)',anchor='w')
		EmaxXANES.grid(column=0, row=2, columnspan=1, rowspan=1, sticky='N')
		self.emaxxanesvar = tkinter.DoubleVar()
		self.emaxxanesentry = tkinter.Entry(self.InterpolationFrame, textvariable=self.emaxxanesvar)
		self.emaxxanesentry.grid(column=1, row=2, columnspan=2)   
		
		EmaxEXAFS = tkinter.Label(self.InterpolationFrame, text='Emax EXAFS (eV)',anchor='w')
		EmaxEXAFS.grid(column=0, row=3, columnspan=1, rowspan=1, sticky='N')
		self.emaxexafsvar = tkinter.DoubleVar()
		self.emaxexafsentry = tkinter.Entry(self.InterpolationFrame, textvariable=self.emaxexafsvar)
		self.emaxexafsentry.grid(column=1, row=3, columnspan=2)   
		
		Estep = tkinter.Label(self.InterpolationFrame, text='E Step (eV)',anchor='w')
		Estep.grid(column=0, row=4, columnspan=1, rowspan=1, sticky='N')
		self.estepvar = tkinter.DoubleVar()
		self.estepvar.set(0.25)
		self.estepentry = tkinter.Entry(self.InterpolationFrame, textvariable=self.estepvar)
		self.estepentry.grid(column=1, row=4, columnspan=2)
		
		constkuselabel = tkinter.Label(self.InterpolationFrame, text='Use Constant k',anchor='w')
		constkuselabel.grid(column=0, row=5, columnspan=2, rowspan=1, sticky='E') 
		self.constkUsevar = tkinter.DoubleVar()
		self.constkUsevar.set(1)
		self.constkUseCheck = tkinter.Checkbutton(self.InterpolationFrame, variable=self.constkUsevar, command=self.plot_update_trigger)
		self.constkUseCheck.grid(column=2, row=5)
		
		Kstep = tkinter.Label(self.InterpolationFrame, text='k Step (Å\u207B\u00B9)',anchor='w')
		Kstep.grid(column=0, row=6, columnspan=1, rowspan=1, sticky='N')
		self.kstepvar = tkinter.DoubleVar()
		self.kstepvar.set(0.025)
		self.kstepentry = tkinter.Entry(self.InterpolationFrame, textvariable=self.kstepvar)
		self.kstepentry.grid(column=1, row=6, columnspan=2)
		
################################################################################		
		#self.PadFrame = tkinter.LabelFrame(master=self.f1, width=310, height=0.75*self.screen_height, pady=25, padx=166)
		#self.PadFrame.grid(row=7, column=0)
		#
		#padlabel = tkinter.Label(self.PadFrame, text='    ',anchor='w')
		#padlabel.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='N')
		
################################################################################
		
		self.windowev = 1
		self.Sav_Gol = False
		self.abort = False
		self.active = False
		self.active_val = 0
		self.trigger = 0
		self.initialize_filter_scroll = 0
		
################################################################################	
#FRAME4-2
################################################################################
		tkinter.Button(f2, text='Go to DataRead', command=lambda:raise_frame(self.f1)).grid(row=0, column=0)

################################################################################		
		self.exportFrame = tkinter.LabelFrame(master=f2, width=300, height=0.75*self.screen_height, pady=5, padx=90)
		self.exportFrame.grid(row=1, column=0)
		
		toggle = [('',0),('',1),('',2)]
		self.radio_updown_buttons = []
		self.updownvar = tkinter.IntVar()
		self.updownvar.set(0)
		
		for idx, (text, mode) in enumerate(toggle):
			self.radio_updown_buttons.append(tkinter.Radiobutton(self.exportFrame, text=text, variable=self.updownvar, value=mode))
			self.radio_updown_buttons[-1].grid(row=1, column=idx+1) 
		
		directionlabel = tkinter.Label(self.exportFrame, text='Direction',anchor='w')
		directionlabel.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='N')
		uplabel = tkinter.Label(self.exportFrame, text='Up',anchor='w')
		uplabel.grid(column=1, row=0, columnspan=1, rowspan=1, sticky='N')
		downlabel = tkinter.Label(self.exportFrame, text='Down',anchor='w')
		downlabel.grid(column=2, row=0, columnspan=1, rowspan=1, sticky='N')
		bothlabel = tkinter.Label(self.exportFrame, text='Both',anchor='w')
		bothlabel.grid(column=3, row=0, columnspan=1, rowspan=1, sticky='N')
		
		ProcessSamlabel = tkinter.Label(self.exportFrame, text='Export Sample',anchor='w')
		ProcessSamlabel.grid(column=0, row=2, columnspan=2, rowspan=1, sticky='E') 
		self.ProcessSamUsevar = tkinter.DoubleVar()
		self.ProcessSamUsevar.set(1)
		self.ProcessSamUseCheck = tkinter.Checkbutton(self.exportFrame, variable=self.ProcessSamUsevar)
		self.ProcessSamUseCheck.grid(column=2, row=2)
		
		ProcessReflabel = tkinter.Label(self.exportFrame, text='Export Reference',anchor='w')
		ProcessReflabel.grid(column=0, row=3, columnspan=2, rowspan=1, sticky='E') 
		self.ProcessRefUsevar = tkinter.DoubleVar()
		self.ProcessRefUsevar.set(1)
		self.ProcessRefUseCheck = tkinter.Checkbutton(self.exportFrame, variable=self.ProcessRefUsevar)
		self.ProcessRefUseCheck.grid(column=2, row=3)
		
		ExportMCRlabel = tkinter.Label(self.exportFrame, text='Export MCR',anchor='w')
		ExportMCRlabel.grid(column=0, row=4, columnspan=2, rowspan=1, sticky='E') 
		self.ExportMCRUsevar = tkinter.DoubleVar()
		self.ExportMCRUsevar.set(0)
		self.ExportMCRUseCheck = tkinter.Checkbutton(self.exportFrame, variable=self.ExportMCRUsevar)
		self.ExportMCRUseCheck.grid(column=2, row=4)
		
################################################################################

		self.alignFrame = tkinter.LabelFrame(master=f2, width=300, height=0.75*self.screen_height, pady=5, padx=62)
		self.alignFrame.grid(row=2, column=0)
		
		AutoAlignuselabel = tkinter.Label(self.alignFrame, text='Auto-align',anchor='w')
		AutoAlignuselabel.grid(column=0, row=0, columnspan=2, rowspan=1, sticky='E') 
		self.AutoAlignUsevar = tkinter.DoubleVar()
		self.AutoAlignUsevar.set(0)
		self.AutoAlignUseCheck = tkinter.Checkbutton(self.alignFrame, variable=self.AutoAlignUsevar)
		self.AutoAlignUseCheck.grid(column=2, row=0)
		
		firstcaliblabel = tkinter.Label(self.alignFrame, text='Initial Alignment',anchor='w')
		firstcaliblabel.grid(column=0, row=1, columnspan=2, rowspan=1, sticky='E') 
		self.firstcalibvar = tkinter.DoubleVar()
		self.firstcalibvar.set(0)
		self.firstcalibCheck = tkinter.Checkbutton(self.alignFrame, variable=self.firstcalibvar)
		self.firstcalibCheck.grid(column=2, row=1)
		
		firstcalib = tkinter.Label(self.alignFrame, text='Alignment Energy',anchor='w')
		firstcalib.grid(column=0, row=2, columnspan=1, rowspan=1, sticky='N')
		self.edgestepvar = tkinter.DoubleVar()
		self.edgestepentry = tkinter.Entry(self.alignFrame, textvariable=self.edgestepvar)
		self.edgestepentry.grid(column=1, row=2, columnspan=2)
		
################################################################################		
		self.BatchFrame = tkinter.LabelFrame(master=f2, width=300, height=0.75*self.screen_height, pady=5, padx=65)
		self.BatchFrame.grid(row=3, column=0)
		
		self.buttonexport= tkinter.Button(self.BatchFrame, text='Batch Export', width=30, height=1,command=self.batch_export)
		self.buttonexport.grid(column=0, row=0, columnspan=3, rowspan=1)
		
		Advancedlabel = tkinter.Label(self.BatchFrame, text='Advanced Parameters',anchor='w')
		Advancedlabel.grid(column=0, row=1, columnspan=3, rowspan=1, sticky='N')
		
		CPUNum = tkinter.Label(self.BatchFrame, text='Processes',anchor='w')
		CPUNum.grid(column=0, row=2, columnspan=1, rowspan=1, sticky='N')
		self.CPUNumvar = tkinter.IntVar()
		self.CPUNumvar.set(psutil.cpu_count()/2)
		self.CPUNumentry = tkinter.Entry(self.BatchFrame, textvariable=self.CPUNumvar)
		self.CPUNumentry.grid(column=1, row=2, columnspan=2)

################################################################################		
		self.PadFrame = tkinter.LabelFrame(master=f2, width=310, height=0.75*self.screen_height, pady=168, padx=166)
		self.PadFrame.grid(row=4, column=0)
		
		padlabel = tkinter.Label(self.PadFrame, text='    ',anchor='w')
		padlabel.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='N')
		
################################################################################
			
		raise_frame(self.f1)
		
################################################################################	
#FRAME5
################################################################################
		def raise_cp_frame(frame):
			frame.tkraise()
			
		self.d1 = tkinter.Frame(master=self.frame5)
		self.c1 = tkinter.Frame(master=self.frame5)
		#self.n1 = tkinter.Frame(master=self.frame5)
		
		for frame in (self.d1, self.c1):
			frame.grid(row=0, column=0, sticky='news')
			
################################################################################			
		self.scrollFrame = tkinter.Frame(master=self.d1)
		self.scrollFrame.grid(row=0, column=0, rowspan=1, columnspan=1)
		tkinter.Label(self.scrollFrame, width=57).grid(column=1, row=0)

		self.ds_scroll_scale = tkinter.Scale(self.scrollFrame, orient='horizontal', from_=0, to=1, state='disabled', width=18)
		self.ds_scroll_scale.grid(column=1, columnspan=1, row=0, stick='EW')
		self.datasetlabel = tkinter.Label(self.scrollFrame, text='Dataset',anchor='w')
		self.datasetlabel.grid(column=0, row=0, columnspan=1, rowspan=1, sticky='S')
		
		self.bf_scroll_scale = tkinter.Scale(self.scrollFrame, orient='horizontal', from_=0.0005, to=0.5, digits = 4, resolution = 0.0005, state='normal', width=18)
		self.bf_scroll_scale.grid(column=1, columnspan=1, row=1, stick='EW')
		self.bflabel = tkinter.Label(self.scrollFrame, text='Filter',anchor='w')
		self.bflabel.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='S')
		
		self.cp_scroll_scale = tkinter.Scale(self.scrollFrame, orient='horizontal', from_=0, to=1, digits = 6, resolution = 0.00005, state='disabled', width=18)
		self.cp_scroll_scale.grid(column=1, columnspan=1, row=2, stick='EW')
		self.cplabel = tkinter.Label(self.scrollFrame, text='Calibration',anchor='w')
		self.cplabel.grid(column=0, row=2, columnspan=1, rowspan=1, sticky='S')
		
################################################################################		
		raise_cp_frame(self.d1)
		
################################################################################	
#FUNCTIONS
################################################################################
		
################################################################################  
	def rootsfile_check(self, folder):
		try:
			self.calibroots = (pd.read_csv(folder+'/roots.dat', sep='\t', header=None)).values
			return True
		except:
			#print('Encoder Analysis Required')
			return False
        
################################################################################
	def calibfile_check(self, folder):
		try:
			self.calibration_values = (pd.read_csv(folder+'/calibration.dat', sep='\t', header=None)).values
			return True
		except:
			#print('Cannot Find Calibration Parameter File')
			self.CalibUsevar.set(0)
			return False
        
################################################################################
	def normfile_check(self, folder):
		try:
			self.normalisation_values = (pd.read_csv(folder+'/normalisation.dat', sep='\t', header=None)).values
			return True
		except:
			#print('Cannot Find Normalisation Parameter File')
			self.NormUsevar.set(0)
			return False
    		   
################################################################################	   
	def fileplus(self):
		self.edge_submitted = False
		file = askopenfile(filetypes = (("QEXAFS Data","*.bin;*.qex"),("all files", "*.*")))
		try:
			dirname, file_short = os.path.split(os.path.abspath(file.name))
			print('Loading File')
			self.edge_box(file)
		except:
			print('Aborted Load') 

################################################################################      
	def edge_box(self,file):

		self.box = tkinter.Toplevel()
		self.box.label0 =tkinter.Label(self.box, text='Enter Absorbing Atom')
		self.box.label0.grid(column=0, row=0, columnspan=2)
		
		self.box.entry0 = tkinter.Entry(self.box)
		self.box.entry0.grid(column=0, row=1, columnspan=2)
		
		self.box.label1 = tkinter.Label(self.box, text='Enter Edge type')
		self.box.label1.grid(column=0, row=2, columnspan=2)
		
		self.box.entry1 = tkinter.StringVar(self.box)
		self.box.entry1.set('K')
		self.box.entry1Menu = tkinter.OptionMenu(self.box, self.box.entry1, 'K','L1','L2','L3')
		self.box.entry1Menu.grid(column=0, row=3, columnspan=2)
		
		self.box.button1 = tkinter.Button(self.box, text='Submit', command=lambda:self.submit_edge(self.box, file))
		self.box.button1.grid(column=0, row=4, columnspan=1)
		
		self.box.button2 = tkinter.Button(self.box, text='Cancel', command=lambda:self.box.destroy())
		self.box.button2.grid(column=1, row=4, columnspan=1)
		
		self.box.diffcaliblabel = tkinter.Label(self.box, text='Different Calibration?',anchor='w')
		self.box.diffcaliblabel.grid(column=2, row=4, columnspan=2, rowspan=1, sticky='E') 
		
		self.box.diffcalibUsevar = tkinter.IntVar()
		self.box.diffcalibUsevar.set(0)	
		self.box.diffcalibUseCheck = tkinter.Checkbutton(self.box, variable=self.box.diffcalibUsevar, command=lambda:self.naccheck(self.box, self.box.diffcalibUsevar))
		self.box.diffcalibUseCheck.grid(column=4, row=4)
		
		self.box.label2 =tkinter.Label(self.box, text='Calibration Element')
		
		self.box.entry2 = tkinter.Entry(self.box)
		
		self.box.label3 = tkinter.Label(self.box, text='Enter Edge type')

		self.box.entry3 = tkinter.StringVar(self.box)
		self.box.entry3.set('K')
		self.box.entry3Menu = tkinter.OptionMenu(self.box, self.box.entry3, 'K','L1','L2','L3')

################################################################################# 
#TO IMPLEMENT
################################################################################# 
	def naccheck(self, box, var):
		if var.get() == 0:
			box.entry3Menu.grid_remove()
			box.label3.grid_remove()
			box.label2.grid_remove()
			box.entry2.grid_remove()
		else:
			box.entry3Menu.grid(column=2, row=3, columnspan=3)
			box.label3.grid(column=2, row=2, columnspan=3)
			box.label2.grid(column=2, row=0, columnspan=3)
			box.entry2.grid(column=2, row=1, columnspan=3)
 
################################################################################          
	def submit_edge(self, box, file):
		self.central_atom = box.entry0.get()
		self.edge_type.append(box.entry0.get()+'_'+box.entry1.get())
		
		xdb = XrayDB()
		if box.diffcalibUsevar.get() == 0:	
			self.edge_info_calib.append(0)
			self.core_hole_calib = 0
			self.edge_type_calib.append(0)
			self.diff_calib = False
		else:
			self.edge_info_calib.append(xdb.xray_edge(box.entry2.get(), box.entry3.get()))
			self.core_hole_calib = xdb.corehole_width(box.entry0.get(), box.entry1.get())
			self.edge_type_calib.append(box.entry2.get()+'_'+box.entry3.get())
			self.diff_calib = True
			print(self.edge_info_calib[-1])
			
		self.edge_info.append(xdb.xray_edge(box.entry0.get(), box.entry1.get()))
		print(self.edge_info[-1])

		self.core_hole = xdb.corehole_width(box.entry0.get(), box.entry1.get())
		print('core hole lifetime',self.core_hole)
		
		self.files_long.append(file)
		self.import_type.append(os.path.splitext(file.name)[-1])
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
			del self.edge_info_calib[i]
			del self.edge_type_calib[i]
			del self.import_type[i]
		
		self.update()
		
################################################################################
	def onselect(self, evt):
	
		self.f1.tkraise()
	
		if 'Data Extract' in App.get_selected(self):
			try:
				data_file_long = self.files_long[self.listbox.curselection()[0]].name
				self.active = True
				self.active_val = self.listbox.curselection()[0]
				
				self.data_file = os.path.splitext(data_file_long)[0]
				self.edge = self.edge_type[self.listbox.curselection()[0]]
				print('Edge Type :', self.edge)
				self.folder = 'output_'+self.edge+'_'+self.data_file.split('/')[-1]
				print('Folder :', self.folder)
				self.file_type = self.import_type[self.listbox.curselection()[0]]
				print('File Type :', self.file_type)
				try:
					if self.diff_calib == False:
						self.calibedgevar.set(self.edge_type[self.active_val])
						self.calibenervar.set(self.edge_info[self.active_val].edge)
						self.e0_scroll_scale.set(self.calibenervar.get())
					else:
						self.calibedgevar.set(self.edge_type_calib[self.active_val])
						self.calibenervar.set(self.edge_info_calib[self.active_val].edge)
						self.e0_scroll_scale.set(self.edge_info[self.active_val].edge)
						
				except:
					self.calibedgevar.set('')
					self.calibenervar.set('')
			except:
				if self.active == False:
					print('No File Selected')
					
			self.clear_figure()
			
			if '.bin' in self.file_type:
				self.header_read_bin(self.data_file)
				self.header_read_bin_ch(self.data_file)
				
			if '.qex' in self.file_type:	
				self.header_read_qex(self.data_file)
				
			print('Data Collection Start Date and Time', self.DataStartTime[0])
			
			self.Fbutter = 1.5*(self.e0_scroll_scale.get() - (1239.852/(np.sin(((np.arcsin(1239.852/(self.twod.get()*self.e0_scroll_scale.get()))*180/np.pi) + 0.00005)*np.pi/180) * self.twod.get())))
			self.Fbutter = round(0.00005*round(self.Fbutter/0.00005), 2)
			
			print('Estimated Butterworth Filter Value :', self.Fbutter)
					
			self.choices = []
			for i in range(int(self.nChannels)):
				self.choices.append('CH '+str(int(i)))
			self.choices.append('1')
			
			self.sam_numerMenu.destroy()
			self.sam_denomMenu.destroy()
			self.ref_numerMenu.destroy()
			self.ref_denomMenu.destroy()
			
			self.sam_numer = tkinter.StringVar(self.plotcolsFrame)
			self.sam_numerMenu = tkinter.OptionMenu(self.plotcolsFrame, self.sam_numer, *self.choices)
			self.sam_numerMenu.grid(column=1, row=1, columnspan=1) 
			
			self.sam_denom = tkinter.StringVar(self.plotcolsFrame)
			self.sam_denomMenu = tkinter.OptionMenu(self.plotcolsFrame, self.sam_denom, *self.choices)
			self.sam_denomMenu.grid(column=1, row=2, columnspan=1) 
			
			self.ref_numer = tkinter.StringVar(self.plotcolsFrame)
			self.ref_numerMenu = tkinter.OptionMenu(self.plotcolsFrame, self.ref_numer, *self.choices)
			self.ref_numerMenu.grid(column=2, row=1, columnspan=1) 
			
			self.ref_denom = tkinter.StringVar(self.plotcolsFrame)
			self.ref_denomMenu = tkinter.OptionMenu(self.plotcolsFrame, self.ref_denom, *self.choices)
			self.ref_denomMenu.grid(column=2, row=2, columnspan=1)
			
			self.scroll_scale_init = 0
			
			if '.bin' in self.file_type:
				self.sam_numer.set('CH 0')
				self.sam_denom.set('CH 1')
				try:
					self.ref_denom.set('CH 2')
					self.ref_numer.set('CH 1')
				except:
					print('Not enough channels in file for reference')
			if '.qex' in self.file_type:
				self.sam_numer.set('CH 0')
				self.sam_denom.set('CH 1')
				try:
					self.ref_denom.set('CH 2')
					self.ref_numer.set('CH 1')
				except:
					print('Not enough channels in file for reference')
					
			if self.active == True:
				self.InterpUsevar.set(0)
				self.rootsfile = self.rootsfile_check(self.folder)
				self.calibrated = self.calibfile_check(self.folder)
				self.normalized = self.normfile_check(self.folder)
				self.d1.tkraise()
				self.normframe_init = 0
				self.initialize_ds_scroll = 0
				
				if (self.rootsfile == True):
					self.CalibrationFrame_visible()
					data_file_long = self.files_long[self.active_val].name
					
					if self.calibrated == False:
						self.plot_spectrum()
						self.NormalisationFrame_invisible()
						self.InterpolationFrame_invisible()
						self.ToProcessFrame_invisible()
						
					else:
						self.CalibUsevar.set(1)
						self.NormalisationFrame_visible()
						self.InterpolationFrame_visible()
						self.ToProcessFrame_visible()
						if self.normalized == False:
							self.plot_spectrum()
						else:
							self.NormUsevar.set(1)
							self.plot_spectrum()
				else:
					self.CalibrationFrame_invisible()	
					self.NormalisationFrame_invisible()
					self.InterpolationFrame_invisible()
					self.ToProcessFrame_invisible()
	
################################################################################					
	def clear_figure(self):
		self.fig.clear()
		self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
		self.canvas.get_tk_widget().grid(column=0, row=0, sticky='EW')
		
		self.ax=self.fig.add_subplot(111)
		self.toolbar.destroy()
		self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
		self.toolbar.update()
		
################################################################################				
	def update_figure(self,x,y, labelx, labely):
		self.clear_figure()
		[self.line] = self.ax.plot(x, y, linewidth=0.75)
		self.ax.set_xlabel(labelx)
		self.ax.set_ylabel(labely)
		self.canvas.draw_idle()

################################################################################	
	def ToProcessFrame_invisible(self):
		self.ToProcessFrame.grid_remove()

################################################################################	
	def ToProcessFrame_visible(self):
		self.ToProcessFrame.grid(row=5, column=0)
		
################################################################################	
	def CalibrationFrame_invisible(self):
		self.CalibrationFrame.grid_remove()
		
################################################################################
	def CalibrationFrame_visible(self):
		self.CalibrationFrame.grid(row=2, column=0)
		
################################################################################				
	def NormalisationFrame_invisible(self):
		self.NormalisationFrame.grid_remove()
		
################################################################################
	def NormalisationFrame_visible(self):
		self.NormalisationFrame.grid(row=3, column=0)
		
################################################################################
	def normframe_load_values(self):
		self.pre1,self.pre2,self.post1,self.post2,self.order_pre,self.order_post,self.e0 = self.normalisation_values
		
		self.pre1_scroll_scale = tkinter.Scale(self.NormalisationFrame, orient='horizontal', from_=np.min(self.energy)-float(self.e0)+5, to=0, state='normal', digits=4, resolution=0.1, length=125)
		self.pre1_scroll_scale.grid(column=1, columnspan=1, row=5, stick='EW')
		self.pre1_scroll_scale.set(float(self.pre1))
		self.pre2_scroll_scale = tkinter.Scale(self.NormalisationFrame, orient='horizontal', from_=np.min(self.energy)-float(self.e0)+5, to=0, state='normal', digits=4, resolution=0.1, length=125)
		self.pre2_scroll_scale.grid(column=1, columnspan=1, row=6, stick='EW')
		self.pre2_scroll_scale.set(float(self.pre2))
		self.post1_scroll_scale = tkinter.Scale(self.NormalisationFrame, orient='horizontal', from_=0, to=np.max(self.energy)-float(self.e0)-10, state='normal', digits=4, resolution=0.1, length=125)
		self.post1_scroll_scale.grid(column=1, columnspan=1, row=7, stick='EW')
		self.post1_scroll_scale.set(float(self.post1))
		self.post2_scroll_scale = tkinter.Scale(self.NormalisationFrame, orient='horizontal', from_=0, to=np.max(self.energy)-float(self.e0)-10, state='normal', resolution=0.1, length=125)
		self.post2_scroll_scale.grid(column=1, columnspan=1, row=8, stick='EW')
		self.post2_scroll_scale.set(float(self.post2))
		self.e0_scroll_scale = tkinter.DoubleVar()
		self.e0_entry = tkinter.Entry(self.NormalisationFrame, textvariable=self.e0_scroll_scale, state='normal')
		self.e0_entry.grid(column=1, columnspan=1, row=9, stick='EW')
		self.e0_scroll_scale.set(float(self.e0))
		self.update_e0_button=tkinter.Button(self.NormalisationFrame, text='Update E0').grid(column=2, row=9)
		
################################################################################				
	def InterpolationFrame_invisible(self):
		self.InterpolationFrame.grid_remove()
		
################################################################################
	def InterpolationFrame_visible(self):
		self.InterpolationFrame.grid(row=4, column=0)
		
################################################################################
	def InterpolationFrame_initialise(self):
		vbf = 1
		
################################################################################
	def Split_buttons_visible(self):
		self.split_accept_button=tkinter.Button(self.EncoderFrame, text='Accept', state='normal', width=14, command=self.accept_split).grid(column=0, row=5, stick='E')
		self.split_insert_button=tkinter.Button(self.EncoderFrame, text='Insert', state='normal', width=14, command=self.insert_split).grid(column=1, row=5, stick='EW')
		self.split_delete_button=tkinter.Button(self.EncoderFrame, text='Delete', state='normal', width=14, command=self.delete_split).grid(column=2, row=5, stick='W')
		
################################################################################
	def Split_buttons_invisible(self):
		self.split_accept_button=tkinter.Button(self.EncoderFrame, text='Accept', state='disabled', width=14, command=self.accept_split).grid(column=0, row=5, stick='E')
		self.split_insert_button=tkinter.Button(self.EncoderFrame, text='Insert', state='disabled', width=14, command=self.insert_split).grid(column=1, row=5, stick='EW')
		self.split_delete_button=tkinter.Button(self.EncoderFrame, text='Delete', state='disabled', width=14, command=self.delete_split).grid(column=2, row=5, stick='W')

################################################################################
#FOR FUTURE IMPLEMENTATION
################################################################################
	def header_read(self, filename):
		if 'bin' in self.file_type:
			self.header_read_bin(self.data_file)
			self.header_read_bin_ch(self.data_file)
		elif 'qex' in self.file_type:
			self.header_read_qex(self.data_file)
			
################################################################################
#FOR FUTURE IMPLEMENTATION
################################################################################
	def data_read(self, root_i, minr):
		if 'bin' in self.file_type:
			values = self.data_read_bin(self, root_i, minr)
			return values
		elif 'qex' in self.file_type:
			values = self.data_read_qex(self, root_i, minr)
			return values		
		
################################################################################     
	def header_read_bin(self,filename):
		#print('Reading Binary Header')
		encoder_bin = filename+'_Encoder'
		f = open(encoder_bin+'.bin', 'rb')
		f.seek(0)
		data = np.fromfile(f, dtype=np.int32, count = 2)
		self.headerSize = data[0]
		fileVersion = data[1]
		#print("HeaderSize: ", self.headerSize)
		data = np.fromfile(f, dtype=np.int64, count = 1)
		self.DataStartTime = str(data)
		data = np.fromfile(f, dtype='f4', count = 2)
		self.AdcClock_Hz = data[1]
		self.samplingvar.set(float(self.AdcClock_Hz/1E6))
		self.DacClock_Hz = 1
		#print("AdcClock_Hz:", self.AdcClock_Hz)
		self.nData = int((os.path.getsize(encoder_bin+'.bin') - self.headerSize)/4)
		#print("nData", self.nData)

################################################################################ 
#FOLD THIS INTO HEADER_READ_BIN
################################################################################ 
	def header_read_bin_ch(self,filename):
		g = open(filename+'.bin', 'rb')
		data = np.fromfile(g, dtype=np.int32, count = 2)
		self.headerSize = data[0]
		data = np.fromfile(g, dtype=np.int64, count = 1)
		data = np.fromfile(g, dtype='f4', count = 2)
		self.nChannels = int(data[0])
		
################################################################################   
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
		self.DataStartTime = ([x for x in header_lines if 'StartOfMeasurement' in x][0].split(': ')[1]).split(', ')
	
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
		
################################################################################ 
#IMPLEMENT ERRORS
################################################################################ 
	def data_read_qex(self, root_i, minr):
		headerSize, line_bytes, dt, nData = self.header_read_qex(self.data_file)
	
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
		if self.BlackmanHarrisFiltervar.get() == 1:
			RawData = RawData.rolling(int(self.blackmanfilterwindowvar.get()), win_type='blackmanharris', min_periods=1, center=True).mean()
			RawData.dropna(how='any', inplace = True)
		
		#RawData.columns = ['ang','mu','sem']
		RawData = RawData[(RawData['ang'] >= np.float(self.min_ang)) & (RawData['ang'] <= np.float(self.max_ang))]
		
		RawData.sort_values(by='ang', ascending=False)
		RawData.reset_index(drop=True)
		
		#print(RawData.shape, np.min(RawData['ang']), np.max(RawData['ang']))
		return np.around(RawData['ang'].values, decimals=6), RawData['mu'].values
		
################################################################################  
	def data_read_bin(self, root_i, minr):
		self.header_read_bin(self.data_file)
		self.header_read_bin_ch(self.data_file)
		
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
		if self.BlackmanHarrisFiltervar.get() == 1:
			RawData = RawData.rolling(int(self.blackmanfilterwindowvar.get()), win_type='blackmanharris', min_periods=1, center=True).mean()
			RawData.dropna(how='any', inplace = True)

		RawData = RawData[(RawData['ang'] >= np.float(self.min_ang)) & (RawData['ang'] <= np.float(self.max_ang))]
		
		RawData.sort_values(by='ang', ascending=False)
		RawData.reset_index(drop=True)
		
		return np.around(RawData['ang'].values, decimals=5), RawData['mu'].values, RawData['sem'].values
		
################################################################################ 
	def plot_add_line(self, x, y):
		[line] = self.ax.plot(x, y) 
		self.canvas.draw_idle()
		return line

################################################################################ 		
	def plot_add_markerlines(self):
		ylims = self.ax.get_ylim()
		[self.pre1line] = self.ax.plot([self.e0_scroll_scale.get()+self.pre1_scroll_scale.get(), self.e0_scroll_scale.get()+self.pre1_scroll_scale.get()], ylims, ls='--', c='r')
		[self.pre2line] = self.ax.plot([self.e0_scroll_scale.get()+self.pre2_scroll_scale.get(), self.e0_scroll_scale.get()+self.pre2_scroll_scale.get()], ylims, ls='--', c='r')
		[self.post1line] = self.ax.plot([self.e0_scroll_scale.get()+self.post1_scroll_scale.get(), self.e0_scroll_scale.get()+self.post1_scroll_scale.get()], ylims, ls='--', c='r')
		[self.post2line] = self.ax.plot([self.e0_scroll_scale.get()+self.post2_scroll_scale.get(), self.e0_scroll_scale.get()+self.post2_scroll_scale.get()], ylims, ls='--', c='r')

		self.canvas.draw_idle()

################################################################################ 			
	def plot_update_trigger(self):
		self.plot_spectrum()
		
################################################################################ 			
	def plot_update_trigger2(self, evt):
		self.plot_spectrum()

################################################################################  
	def plot_spectrum(self):
		self.rootsfile = self.rootsfile_check(self.folder)
		self.calibrated = self.calibfile_check(self.folder)
		self.normalized = self.normfile_check(self.folder)
		
		self.initialize_plot = 0
		
		def plot_data(evt):
			if self.trigger != 3:
				root_i = int(self.ds_scroll_scale.get())
			else:
				root_i = self.align_root_i
				self.ds_scroll_scale.set(root_i)
			
			if '.bin' in self.file_type:
				try:
					ang,mu,std = self.data_read_bin(root_i, self.minr)
				except:
					print('Please check column selection, unable to read selected columns')
				
			elif '.qex' in self.file_type:
				try:
					ang,mu = self.data_read_qex(root_i, self.minr)
				except:
					print('Please check column selection, unable to read selected columns')
			
			if self.Filtervar.get() == 1:
				#nyq = 0.5 * len(ang) * 2
				self.Wn = float(self.bf_scroll_scale.get())
				#high = 2 * nyq * float(self.bf_scroll_scale.get())/nyq
				N  = 3    # Filter order
				B, A = signal.butter(N, self.Wn, output='ba')
				mu_raw = mu
				mu = signal.filtfilt(B,A,mu)
			
			if self.trigger == 0:
				if (self.calibrated == False) and (self.Filtervar.get() == 0):
					mu_raw = mu	
			
				if (self.calibrated == True)  & (self.CalibUsevar.get()==1):
					self.calibration_values = (pd.read_csv(self.folder+'/calibration.dat', sep='\t', header=None)).values
					self.energy = (1239.852/(float(self.calibration_values[3])*np.sin((ang+(self.calibration_values[1]-self.calibration_values[0]))*np.pi/180)))
					if self.Filtervar.get() == 1:
						energy_raw = self.energy
						
					if (self.NormUsevar.get() == 0):
						self.InterpolationFrame_initialise()
					
					if (self.normalized == True) & (self.NormUsevar.get() == 1):
						self.normalisation_values = pd.read_csv(self.folder+'/normalisation.dat', sep='\t', header=None).values
						self.pre1,self.pre2,self.post1,self.post2,self.order_pre,self.order_post,self.e0 = self.normalisation_values
						if self.normframe_init == 0:
							self.normframe_init = 1
							self.normframe_load_values()
							self.InterpolationFrame_initialise()
						
						mu = self.normalize_data(self.energy, mu)
						if '.bin' in self.file_type:
							if self.Filtervar.get() == 1:
								error_estimation = np.empty(len(mu))
								error_estimation.fill(np.std(mu[int(-np.floor(len(mu)/5)):]))
								std = error_estimation
							else:
								std = (std / self.edge_jump) * 1000 / len(std)
						if self.Filtervar.get()==1:
							mu_raw = self.normalize_data(self.energy, mu_raw)
							
					if self.Filtervar.get() == 0:
						mu_raw = mu	
						energy_raw = self.energy
					
					if self.InterpUsevar.get() == 1:
						self.xnew = self.xnew_gen()
						self.energy, mu = interpolate.localised(self, self.energy, mu, self.xnew)
						if '.bin' in self.file_type:
							self.energy, std = interpolate.localised(self, self.energy, std, self.xnew)
					
					if self.initialize_plot == 0:
						self.initialize_plot = 1
						if (self.normalized == True) & (self.NormUsevar.get() == 1):
							self.update_figure(energy_raw, mu_raw, 'Energy (eV)', 'Normalised Absorption')
						else: 
							self.update_figure(energy_raw, mu_raw, 'Energy (eV)', 'Absorption')
							
						if (self.Filtervar.get() == 1) or (self.InterpUsevar.get() == 1):
							self.filter_line = self.plot_add_line(self.energy, mu)
					else:
						self.line.set_data(energy_raw, mu_raw)
						if (self.Filtervar.get() == 1) or (self.InterpUsevar.get() == 1):
							self.filter_line.set_data(self.energy, mu)
						self.canvas.draw_idle()

				else:
					if self.Filtervar.get() == 0:
							mu_raw = mu
					if self.initialize_plot == 0:
						self.initialize_plot = 1
						self.update_figure(ang, mu_raw, 'Monochromator Encoder Angle (°)', 'Absorption')
						if self.Filtervar.get() == 1:
							self.filter_line = self.plot_add_line(ang, mu)
					else:
						self.line.set_data(ang, mu_raw)
						if self.Filtervar.get() == 1:
							self.filter_line.set_data(ang, mu)
					self.canvas.draw_idle()
				
			if self.trigger == 1:
				if self.CalibFiltervar.get() == 1:
					self.Wn = float(self.bf_scroll_scale.get())
					N  = 3    # Filter order
					self.B, self.A = signal.butter(N, self.Wn, output='ba')
					mu = signal.filtfilt(self.B,self.A,mu)
			
				deriv_data = -np.gradient(mu)
				self.update_figure(ang, deriv_data, 'Monochromator Encoder Angle (°)', 'Absorption Derivative')
				self.cp_scroll_scale.set(ang[np.argmax(deriv_data)])
				self.calibline = self.ax.axvline(ang[np.argmax(deriv_data)], ls='--', c='orange')
				self.ax.autoscale(axis='y', tight=False)
				self.canvas.draw_idle()
			
			if self.trigger == 2:
				
				self.energy = (1239.852/(float(self.calibration_values[3])*np.sin((ang+(self.calibration_values[1]-self.calibration_values[0]))*np.pi/180)))
				
				if self.Filtervar.get() == 1:
					self.Wn = float(self.bf_scroll_scale.get())
					N  = 3    # Filter order
					B, A = signal.butter(N, self.Wn, output='ba')
					mu = signal.filtfilt(B,A,mu)
				
				self.normalisation_values = [self.pre1_scroll_scale.get(), self.pre2_scroll_scale.get(),self.post1_scroll_scale.get(),self.post2_scroll_scale.get(),self.norm_pre.get(),self.norm_post.get(),self.e0_scroll_scale.get()]
				self.pre1,self.pre2,self.post1,self.post2,self.order_pre,self.order_post,self.e0 = self.normalisation_values
				
				linepre_y = self.regression('pre-edge', np.asarray(self.energy), np.asarray(mu))
				linepost_y = self.regression('post-edge', np.asarray(self.energy), np.asarray(mu))
				
				if self.initialize_plot == 0:
					self.initialize_plot = 1
					self.update_figure(self.energy, mu, 'Energy (eV)', 'Absorption')
					self.line_pre = self.plot_add_line(self.energy, linepre_y)
					self.line_post = self.plot_add_line(self.energy, linepost_y)
					self.plot_add_markerlines()
					self.canvas.draw_idle()	
				else:
					self.line.set_data(self.energy, mu)
					self.line_pre.set_data(self.energy, linepre_y)
					self.line_post.set_data(self.energy, linepost_y)
					ylims = self.ax.get_ylim()
					self.pre1line.set_xdata([self.e0_scroll_scale.get()+self.pre1_scroll_scale.get(), self.e0_scroll_scale.get()+self.pre1_scroll_scale.get()])
					self.pre2line.set_xdata([self.e0_scroll_scale.get()+self.pre2_scroll_scale.get(), self.e0_scroll_scale.get()+self.pre2_scroll_scale.get()])
					self.post1line.set_xdata([self.e0_scroll_scale.get()+self.post1_scroll_scale.get(), self.e0_scroll_scale.get()+self.post1_scroll_scale.get()])
					self.post2line.set_xdata([self.e0_scroll_scale.get()+self.post2_scroll_scale.get(), self.e0_scroll_scale.get()+self.post2_scroll_scale.get()])
					self.canvas.draw_idle()

				index = min(range(len(self.energy)), key=lambda i: abs(self.energy[i]-self.e0_scroll_scale.get()))
				self.edge_jump = np.around((linepost_y-linepre_y)[index], decimals=2)
				self.evjlabel.set(str(self.edge_jump))
								
			if self.trigger == 3:
				if self.Filtervar.get() == 1:
					self.Wn = float(self.bf_scroll_scale.get())
					N  = 3    # Filter order
					B, A = signal.butter(N, self.Wn, output='ba')
					mu_raw = mu
					mu = signal.filtfilt(B,A,mu)

				self.calibration_values = (pd.read_csv(self.folder+'/calibration.dat', sep='\t', header=None)).values
				self.energy = (1239.852/(float(self.calibration_values[3])*np.sin((ang+(self.calibration_values[1]-self.calibration_values[0]))*np.pi/180)))
				if self.Filtervar.get() == 1:
					energy_raw = self.energy
				
				self.normalisation_values = pd.read_csv(self.folder+'/normalisation.dat', sep='\t', header=None).values
				self.pre1,self.pre2,self.post1,self.post2,self.order_pre,self.order_post,self.e0 = self.normalisation_values
				if self.normframe_init == 0:
					self.normframe_init = 1
					self.normframe_load_values()
					
				mu = self.normalize_data(self.energy, mu)
				if self.Filtervar.get()==1:
					mu_raw = self.normalize_data(energy_raw, mu_raw)
					
				if self.InterpUsevar.get() == 1:
					self.xnew = self.xnew_gen()
					if self.Filtervar.get() == 0:
						energy_raw = self.energy
						mu_raw = mu
					self.energy, mu = interpolate.localised(self, self.energy, mu, self.xnew)
					
				if np.min(mu) < -0.2:
					pass_var = False
					popt = [0,0]
					return popt, pass_var
				else:
					pass_var = True
					popt, pcov = curve_fit(self.cumgauss, self.energy, mu, (float(self.e0), float(self.core_hole)))
					print('Alignment Parameters', popt)
					self.edgestepvar.set(popt[0])
					self.edgestepwidth = popt[1]
					
					if self.Filtervar.get() == 1:
						self.update_figure(energy_raw, mu_raw, 'Energy (eV)', 'Normalised Absorption')
						self.plot_add_line(self.energy, mu)
					else:
						self.update_figure(self.energy, mu, 'Energy (eV)', 'Normalised Absorption')
					yfit = self.cumgauss(self.energy, popt[0], popt[1])
					self.plot_add_line(self.energy, yfit)
					
					return pass_var
					
		if self.trigger != 3:
			if (self.rootsfile == True):		
				roots_file_data = pd.read_csv(self.folder+'/roots.dat', sep='\t', header=None).values
				self.minr = roots_file_data[2]
				self.max_ang = roots_file_data[1]
				self.min_ang = roots_file_data[0]
				
				if self.initialize_ds_scroll==0:
					self.ds_scroll_scale = tkinter.Scale(self.scrollFrame, orient='horizontal', from_=0, to=len(roots_file_data)-5, state='normal', width=18, command=plot_data)
					self.ds_scroll_scale.grid(column=1, columnspan=1, row=0, stick='EW')
					self.datasetlabel = tkinter.Label(self.scrollFrame, text='Dataset',anchor='w')
					self.datasetlabel.grid(column=0, row=0, columnspan=1, rowspan=1, sticky='S')
					self.initialize_ds_scroll = 1

				if self.trigger == 2:
					
					def scroll_bar_update(*args):
						self.pre1_scroll_scale = tkinter.Scale(self.NormalisationFrame, orient='horizontal', from_=np.min(self.energy)-self.e0_scroll_scale.get()+5, to=0, state='normal', digits=4, resolution=0.1, length=125, command=plot_data)
						self.pre2_scroll_scale = tkinter.Scale(self.NormalisationFrame, orient='horizontal', from_=np.min(self.energy)-self.e0_scroll_scale.get()+5, to=0, state='normal', digits=4, resolution=0.1, length=125, command=plot_data)
						self.post1_scroll_scale = tkinter.Scale(self.NormalisationFrame, orient='horizontal', from_=0, to=np.max(self.energy)-self.e0_scroll_scale.get()-10, state='normal', digits=4, resolution=0.1, length=125, command=plot_data)
						self.post2_scroll_scale = tkinter.Scale(self.NormalisationFrame, orient='horizontal', from_=0, to=np.max(self.energy)-self.e0_scroll_scale.get()-10, state='normal', resolution=0.1, length=125, command=plot_data)
						self.e0_entry = tkinter.Entry(self.NormalisationFrame, textvariable=self.e0_scroll_scale, state='normal')
						self.update()
						
						self.update_e0_button=tkinter.Button(self.NormalisationFrame, text='Update E0', command=callback).grid(column=2, row=9)
						
						self.pre1_scroll_scale.grid(column=1, columnspan=1, row=5, stick='EW')
						self.pre2_scroll_scale.grid(column=1, columnspan=1, row=6, stick='EW')
						self.post1_scroll_scale.grid(column=1, columnspan=1, row=7, stick='EW')
						self.post2_scroll_scale.grid(column=1, columnspan=1, row=8, stick='EW')
				
						self.pre1_scroll_scale.set(np.min(self.energy)-self.e0_scroll_scale.get()+5)
						self.pre2_scroll_scale.set(-50)
						self.post1_scroll_scale.set(100)
						self.post2_scroll_scale.set(np.max(self.energy)-self.e0_scroll_scale.get()-10)
						
						self.scroll_scale_init = 1
														
					def callback(*args):
						if self.e0_scroll_scale.get() > np.min(self.energy):
							scroll_bar_update(None)
							plot_data(None)
				
					if self.scroll_scale_init == 0:
						scroll_bar_update(None)
					
				if self.initialize_filter_scroll==0:
					self.bf_scroll_scale = tkinter.Scale(self.scrollFrame, orient='horizontal', from_=0.0005, to=0.5, digits = 4, resolution = 0.0005, state='normal', width=18, command=plot_data)
					self.bf_scroll_scale.grid(column=1, columnspan=1, row=1, stick='EW')
					self.bf_scroll_scale.set(self.Fbutter)
					self.initialize_filter_scroll = 1
	
				plot_data(None)
				
		else:
			pass_var = plot_data(None)
			return pass_var
		
################################################################################	
#SPLIT			CHANGE ALGORITHM?
################################################################################		
	def split(self):
		self.abort = False
		self.progress_bar["value"] = 0
		self.progress_bar["style"] = "blue.Horizontal.TProgressbar"
		
		self.Split_buttons_visible()
		
		def clear_canvasFrame():
			x,y=[],[]
			self.update_figure(x,y,'Spectrum Number','Number of Points Per Spectrum')
		
		clear_canvasFrame()
		
		print('Anaylzing Encoder')
		
		self.edge_calib = self.calibedgevar.get()
		
		if 'bin' in self.file_type:
			buffer = 32000000
		if 'qex' in self.file_type:
			buffer = self.line_bytes * self.nLines
		
		print(self.data_file)
		
		if not os.path.exists('output_'+self.edge+'_'+self.data_file.split('/')[-1]):
			os.makedirs('output_'+self.edge+'_'+self.data_file.split('/')[-1])
			
		self.folder = 'output_'+self.edge+'_'+self.data_file.split('/')[-1]
		
		if '.bin' in self.file_type :
			self.header_read_bin(self.data_file)
		elif '.qex' in self.file_type:
			self.header_read_qex(self.data_file)

		if 'bin' in self.file_type:
			buffer_points=np.arange(np.ceil(4*self.nData/buffer))
		if 'qex' in self.file_type:
			buffer_points=np.arange(np.ceil(self.line_bytes*self.nLines/buffer))
		
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
		print('number of processor cores available =',int(self.CPUNumvar.get()))
		print('job length =',len(buffer_points))
		self.labelText.set('Starting Processess')
		self.update()
		time.sleep(0.1)
		
		angs2d = [None] * ncpus
		roots2d = [None] * ncpus
		
		if 'bin' in self.file_type:
			encoder_bin = self.data_file+'_Encoder'
		if 'qex' in self.file_type:
			encoder_bin = self.data_file
				
			#f = open(encoder_bin+'.qex', 'rb')
			#f.seek(self.headerSize)
			#data_E = np.fromfile(f, dtype=self.dt, count = int(buffer))

		options = [ncpus, np.max(buffer_points), self.headerSize, buffer, encoder_bin, self.nData, resample_factor, enconder_sampling, encoder_smooth, self.file_type]       
		sub_f = 'batch_split_subroutine_v2.2.py'
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
			if 'max_ang' in stdoutdata.decode('utf-8'):
				self.max_ang_ID0 = float(stdoutdata.decode('utf-8').split('\n')[0].split(' ')[-1])
			if 'min_ang' in stdoutdata.decode('utf-8'):
				self.min_ang_ID0 = float(stdoutdata.decode('utf-8').split('\n')[0].split(' ')[-1])
			if 'print data' in stdoutdata.decode('utf-8'):
				print(stdoutdata.decode('utf-8').split('\n')[0])
			if 'data read' in stdoutdata.decode('utf-8'):
				self.labelText.set('Reading Data')
				perc = float(stdoutdata.decode('utf-8').split('\n')[0].split(' ')[3])
				self.progress_bar["value"] = perc
				if perc == 100:
					self.labelText.set('Searching Peaks')
				self.progress_bar.update()
			if 'Min Encoder Values' in stdoutdata.decode('utf-8'):
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
		self.angs = np.asarray(list(itertools.chain(*angs2d)))	
		self.roots = roots.astype(np.int64)
		print('number of roots:',len(roots))
		#print(angs)
		self.progress_bar["value"] = 0
		self.labelText.set('Plotting')
		
		if 'qex' in self.file_type:
			self.root_test = np.diff(self.roots)
			Ave = np.average(self.root_test)
			
			if self.root_test[0] < 0.95*Ave:
				self.roots = np.delete(self.roots, 0)
				self.angs = np.delete(self.angs, 0)
				
			if self.root_test[-1] < 0.95*Ave:
				self.roots = np.delete(self.roots, -1)
				self.angs = np.delete(self.angs, -1)

		self.minr = self.plotroots(self.angs)
		
################################################################################             
	def plotroots(self,angs):
		#finds the minimum and maximum difference between roots

		self.root_test = np.diff(self.roots)
		x=range(len(self.root_test))
		
		def roots_figure_plot(x, y, xlabel, ylabel):
			self.clear_figure()
			self.ax.scatter(x, y)
			self.ax.set_ylim([0,np.max(self.root_test)*1.1])
			self.ax.set_xlim([0,len(self.root_test)])
			self.canvas.draw_idle()
		
		roots_figure_plot(x, self.root_test, 'Spectrum Number', 'Number of Points Per Spectrum')
		minr = int(np.min(self.root_test))
		return minr
		
################################################################################  
	def insert_split(self):
		Ave = np.average(self.root_test)
		roots_new = []
		
		for i in range(len(self.root_test)):
			roots_new.append(self.roots[i])
			if self.root_test[i] > 1.5*Ave:
				roots_new.append(self.roots[i]+Ave)
				
		self.roots = roots_new
		self.plotroots(self.angs)
		
################################################################################  
#IMPLEMENTATION REQUIRED
################################################################################  
	def delete_split(self):
		Ave = np.average(self.root_test)
		
		print(Ave)
		
		if self.root_test[0] < 0.98*Ave:
			self.roots = np.delete(self.roots, 0)
			self.angs = np.delete(self.angs, 0)
			
		if self.root_test[-1] < 0.98*Ave:
			self.roots = np.delete(self.roots, -1)
			self.angs = np.delete(self.angs, -1)
			
		if self.root_test[0] > 2.5*Ave:
			self.roots = np.delete(self.roots, 0)
			self.angs = np.delete(self.angs, 0)
			
		if self.root_test[-1] > 2.5*Ave:
			self.roots = np.delete(self.roots, -1)
			self.angs = np.delete(self.angs, -1)
			
		self.plotroots(self.angs)

################################################################################  
	def accept_split(self):
		minr = int(np.min(self.root_test))
		B = [i for i in self.angs if i >= np.average(self.angs)]
		C = [i for i in self.angs if i <= np.average(self.angs)]
		max_ang = np.min(B)
		min_ang = np.max(C)
		
		#if 'qex' in self.file_type:
		#	if max_ang > self.max_ang_ID0:
		#		max_ang = self.max_ang_ID0
		#	if min_ang < self.min_ang_ID0:
		#		min_ang = self.min_ang_ID0
		
		print('max angle :',max_ang)
		print('min angle :',min_ang)
		 
		self.roots = np.concatenate([[min_ang],[max_ang],[minr],self.roots])
		roots_info = pd.DataFrame(self.roots)
		roots_info.to_csv(self.folder+'/roots.dat', sep='\t', header=None, index=False)
		
		print('Internal Roots File Created')
		
		self.Split_buttons_invisible()
		
		self.rootsfile = self.rootsfile_check(self.folder)
		if (self.rootsfile == True):
			self.CalibrationFrame_visible()
			data_file_long = self.files_long[self.active_val].name
		else:
			self.CalibrationFrame_invisible()	
		
		self.plot_spectrum()
		
################################################################################	
#CALIBRATE			CHECK 2D VALUES
################################################################################		
	def calibrate(self):
		print('Calibrating')
		self.trigger = 1
		
		roots_file_data = pd.read_csv(self.folder+'/roots.dat', sep='\t', header=None).values
		
		if self.calibmono.get() == 'Si 111':
			self.twod.set(0.627120)
			print('Mono Crystal = Si 111')
		elif self.calibmono.get() == 'Si 311':
			self.twod.set(0.320267)
			print('Mono Cyrstal = Si 311')
			
		self.calib_accept_button=tkinter.Button(self.CalibrationFrame, text='Accept', state='normal', width=14, command=self.accept_calib).grid(column=1, row=8, stick='E')
		self.calib_cancel_button=tkinter.Button(self.CalibrationFrame, text='Cancel', state='normal', width=14, command=self.cancel_calib).grid(column=2, row=8, stick='W')
			
		self.cp_scroll_scale = tkinter.Scale(self.scrollFrame, orient='horizontal', from_=float(self.min_ang), to=float(self.max_ang), digits = 7, resolution = 0.000025, state='normal', width=18, command=self.plot_cp_line)
		self.cp_scroll_scale.grid(column=1, columnspan=1, row=2, stick='EW')
			
		self.plot_spectrum()

################################################################################
	def plot_cp_line(self,evt):
		self.calibline.set_xdata(self.cp_scroll_scale.get())
		self.ax.autoscale(axis='y', tight=False)
		self.canvas.draw_idle()

################################################################################	
	def accept_calib(self):
		self.d1.tkraise()
		self.calib_accept_button=tkinter.Button(self.CalibrationFrame, text='Accept', state='disabled', width=14, command=self.accept_calib).grid(column=1, row=8, stick='E')
		self.calib_cancel_button=tkinter.Button(self.CalibrationFrame, text='Cancel', state='disabled', width=14, command=self.cancel_calib).grid(column=2, row=8, stick='W')
		
		self.thetac = np.arcsin(1239.852/(self.twod.get()*self.calibenervar.get()))*180/np.pi
		calibvalue = [self.cp_scroll_scale.get(), self.thetac, float(self.bf_scroll_scale.get()), self.twod.get(), self.calibenervar.get()]
		pdcalibvalue = pd.DataFrame(calibvalue)
		pdcalibvalue.to_csv(self.folder+'/calibration.dat', sep='\t', header=None, index=False)
		self.CalibUsevar.set(1)
		
		self.trigger = 0
		self.NormalisationFrame_visible()
		self.InterpolationFrame_visible()
		self.ToProcessFrame_visible()
		self.plot_spectrum()
		
################################################################################	
	def cancel_calib(self):
		self.d1.tkraise()
		self.calib_accept_button=tkinter.Button(self.CalibrationFrame, text='Accept', state='disabled', width=14, command=self.accept_calib).grid(column=1, row=8, stick='E')
		self.calib_cancel_button=tkinter.Button(self.CalibrationFrame, text='Cancel', state='disabled', width=14, command=self.cancel_calib).grid(column=2, row=8, stick='W')
		
		self.trigger = 0
		self.plot_spectrum()
		
################################################################################	
#NORMALIZE
################################################################################
	def normalize(self):
		roots_file_data = pd.read_csv(self.folder+'/roots.dat', sep='\t', header=None).values
		self.d1.tkraise()
		
		self.norm_accept_button=tkinter.Button(self.NormalisationFrame, text='Accept', state='normal', width=14, command=self.accept_norm).grid(column=1, row=11, stick='E')
		self.norm_cancel_button=tkinter.Button(self.NormalisationFrame, text='Cancel', state='normal', width=14, command=self.cancel_norm).grid(column=2, row=11, stick='W')
		
		self.trigger = 2
		self.plot_spectrum()
		
################################################################################     
	def normload(self):
		file = askopenfile()
		dirname, file_short = os.path.split(os.path.abspath(file.name))
		shutil.copyfile(file.name,self.folder+'/'+file_short)
		self.NormUsevar.set(1)
		self.Flatvar.set(1)
		self.plot_spectrum() 
			
################################################################################ 
	def normalize_data(self,energy, mu):   	
		linepre_y = self.regression('pre-edge',np.asarray(energy), np.asarray(mu))        
		linepost_y = self.regression('post-edge',np.asarray(energy), np.asarray(mu))
		
		index = min(range(len(energy)), key=lambda i: abs(energy[i]-float(self.e0)))
		normdivisor = (linepost_y-linepre_y)[index]
		self.edge_jump = normdivisor
		
		mu = (mu-linepre_y)/normdivisor
		linepost_y = (linepost_y-linepre_y)/normdivisor
		linepre_y = np.zeros(len(linepre_y))
						
		if self.Flatvar.get() == 1:
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
	def regression(self,region,x,y): 
		
		if region == 'pre-edge':
			xllim = float(self.e0) + float(self.pre1)
			xhlim = float(self.e0) + float(self.pre2)
			try:
				order = str(self.order_pre[0])
			except:
				order = str(self.order_pre)
	
		if region == 'post-edge':
			xllim = float(self.e0) + float(self.post1)
			xhlim = float(self.e0) + float(self.post2)
			try:
				order = str(self.order_post[0])
			except:
				order = str(self.order_post)
		
		xregion = x[(x >= xllim) & (x <= xhlim)]
		yregion = y[(x >= xllim) & (x <= xhlim)]
		
		if (order == '0') or (order == '0.0'):
			popt, pcov = curve_fit(polynomials.constant, xregion, yregion, p0=np.ones(1))
			yfit = polynomials.constant(x, *popt)
	
		if (order == '1') or (order == '1.0'):
			popt, pcov = curve_fit(polynomials.linear, xregion, yregion, p0=np.ones(2))
			yfit = polynomials.linear(x, *popt)
			
		if (order == '2') or (order == '2.0'):
			popt, pcov = curve_fit(polynomials.quadratic, xregion, yregion, p0=np.ones(3))
			yfit = polynomials.quadratic(x, *popt)
			
		if (order == '3') or (order == '3.0'):
			popt, pcov = curve_fit(polynomials.cubic, xregion, yregion, p0=np.ones(4))
			yfit = polynomials.cubic(x, *popt)
			
		if (order == '4') or (order == '4.0'):
			popt, pcov = curve_fit(polynomials.quartic, xregion, yregion, p0=np.ones(5))
			yfit = polynomials.quartic(x, *popt)
			
		if order == 'V':
			popt, pcov = curve_fit(polynomials.victoreen, xregion, yregion, p0=np.ones(3))
			yfit = polynomials.victoreen(x, *popt)
	
		return yfit
		
################################################################################	
	def accept_norm(self):
		self.d1.tkraise()
		self.norm_accept_button=tkinter.Button(self.NormalisationFrame, text='Accept', state='disabled', width=14, command=self.accept_norm).grid(column=1, row=11, stick='E')
		self.norm_cancel_button=tkinter.Button(self.NormalisationFrame, text='Cancel', state='disabled', width=14, command=self.cancel_norm).grid(column=2, row=11, stick='W')
		
		NormValues = [self.pre1_scroll_scale.get(), self.pre2_scroll_scale.get(),self.post1_scroll_scale.get(),self.post2_scroll_scale.get(),self.norm_pre.get(),self.norm_post.get(),self.e0_scroll_scale.get()]
		pdnormvalue = pd.DataFrame(NormValues)
		pdnormvalue.to_csv(self.folder+'/normalisation.dat', sep='\t', header=None, index=False)
		
		print(NormValues)
		
		self.trigger = 0
		self.NormUsevar.set(1)
		self.plot_spectrum()

################################################################################	
	def cancel_norm(self):
		self.d1.tkraise()
		self.norm_accept_button=tkinter.Button(self.NormalisationFrame, text='Accept', state='disabled', width=14, command=self.accept_norm).grid(column=1, row=11, stick='E')
		self.norm_cancel_button=tkinter.Button(self.NormalisationFrame, text='Cancel', state='disabled', width=14, command=self.cancel_norm).grid(column=2, row=11, stick='W')
		
		self.trigger = 0
		self.plot_spectrum()
		
################################################################################	
	def xnew_gen(self):
		if self.constkUsevar.get() == 1:
			k_min = round(np.sqrt((self.emaxxanesvar.get() - float(self.e0_scroll_scale.get()))*0.262468426033256),1)
			k_max = round(np.sqrt((self.emaxexafsvar.get() - float(self.e0_scroll_scale.get()))*0.262468426033256),1)
	
			k_grid = np.arange(k_min, k_max, self.kstepvar.get())
	
			e_grid = (((k_grid)**2)/0.262468426033256) + float(self.e0_scroll_scale.get())
			idxSG = npi.indices(e_grid, e_grid[e_grid >= self.emaxxanesvar.get()])
			e_grid = e_grid[np.min(idxSG):np.max(idxSG)]
	
			e_xanes = np.arange(self.eminvar.get(), self.emaxxanesvar.get()+self.estepvar.get(), self.estepvar.get())
			xnew = np.concatenate((e_xanes, e_grid), axis=0)
		else:
			print('Using constant Estep throughout')
			xnew = np.arange(self.eminvar.get(), self.emaxexafsvar.get(), self.estepvar.get())
		return xnew
		
################################################################################			
#BATCH EXPORT			ADD MORE INFORMATION TO PARAMTERS FILE	
################################################################################	
	def batch_export(self):
		start_batch = time.time()
		print('Beginning Export')
		self.edgestepwidth = self.core_hole
		
		self.rootsfile = self.rootsfile_check(self.folder)
		self.calibrated = self.calibfile_check(self.folder)
		self.normalized = self.normfile_check(self.folder)	
		
		if self.updownvar.get() == 0:
			if not os.path.exists(self.folder+'/Export/Individual_Up'):
				os.makedirs(self.folder+'/Export/Individual_Up')
		elif self.updownvar.get() == 1:
			if not os.path.exists(self.folder+'/Export/Individual_Down'):
				os.makedirs(self.folder+'/Export/Individual_Down')	
		else:
			if not os.path.exists(self.folder+'/Export/Individual_Both'):
				os.makedirs(self.folder+'/Export/Individual_Both')
		
		num_cpus = self.CPUNumvar.get()
		print('number of processor cores available = ',num_cpus)
		print('job length = ',int(len(self.calibroots)-5))
		
		self.labelText.set('Starting Processess')
		self.update()
		time.sleep(0.1)
		
		proceed_var = False
		
		if self.InterpUsevar.get() == 1:
			df_xnew = pd.DataFrame()
			df_xnew['xnew'] = self.xnew
			df_xnew.to_csv(self.folder+'/xnew.dat', sep='\t', header=False, index=False)
		
		if (self.AutoAlignUsevar.get() == 1) & (self.firstcalibvar.get() == 1) & (self.calibrated == True) & (self.normalized == True) & (self.NormUsevar.get()==1) & (self.CalibUsevar.get()==1):
			find_autoalign_energy = True
			self.align_root_i = 0
			toggle_store = self.togglevar.get()
			self.trigger = 3
		elif (self.AutoAlignUsevar.get() == 1) & (self.firstcalibvar.get() == 0) & (self.calibrated == True) & (self.normalized == True) & (self.NormUsevar.get()==1) & (self.CalibUsevar.get()==1):
			auto_align_energy = self.edgestepvar.get()
			if auto_align_energy == 0:
				print('Auto-align energy is zero, performing initial aligment')
				self.firstcalibvar.set(1)
				find_autoalign_energy = True
				self.align_root_i = 0
				toggle_store = self.togglevar.get()
				self.trigger = 3
			else:
				print('Proceeding with previous auto-align energy')
				proceed_var = True
		elif (self.AutoAlignUsevar.get() == 1):
			print('Conditions not met for auto-alignment, proceeding without')
			self.AutoAlignUsevar.set(0)
			self.firstcalibvar.set(0)
			if (self.calibrated == True) & (self.CalibUsevar.get()==1) & (self.NormUsevar.get()==0):
				print('Exporting without auto-alignment and normalisation')
				proceed_var = True
			elif (self.calibrated == True) & (self.CalibUsevar.get()==1) & (self.NormUsevar.get()==1):
				print('Exporting without auto-alignment')
				proceed_var = True
			else:
				print('Conditions for export not met, check for calibration')
				proceed_var = False
		else:
			self.firstcalibvar.set(0)
			if (self.calibrated == True) & (self.CalibUsevar.get()==1) & (self.NormUsevar.get()==0):
				print('Exporting without auto-alignment and normalisation')
				proceed_var = True
			elif (self.calibrated == True) & (self.CalibUsevar.get()==1) & (self.NormUsevar.get()==1):
				print('Exporting without auto-alignment')
				proceed_var = True
			else:
				print('Conditions for export not met, check for calibration')
				proceed_var = False
		
		if (self.AutoAlignUsevar.get() == 1) & (self.firstcalibvar.get() == 1) & (self.calibrated == True) & (self.normalized == True) & (self.NormUsevar.get()==1) & (self.CalibUsevar.get()==1):
			self.togglevar.set(1)
			while find_autoalign_energy == True:
				if self.align_root_i > int(len(self.calibroots)-5):
					print('Failed to find alingment energy')
					self.trigger = 0
					self.togglevar.set(int(toggle_store))
					self.AutoAlignUsevar.set(0)
					self.firstcalibvar.set(0)
					proceed_var = False
					break
				else:
					pass_var_check = self.plot_spectrum()
					if pass_var_check == True:
						find_autoalign_energy = False
						self.ds_scroll_scale.set(self.align_root_i)
						self.trigger = 0
						self.togglevar.set(int(toggle_store))
						self.firstcalibvar.set(0)
						auto_align_energy = self.edgestepvar.get()
						proceed_var = True
					else:
						self.align_root_i += 1 
		
		if proceed_var == True:
			self.update()
			options=[int(len(self.calibroots)-5), self.data_file, self.folder, self.file_type, num_cpus, self.ProcessSamUsevar.get(), self.ProcessRefUsevar.get(), self.NormUsevar.get(), self.InterpUsevar.get(), self.AutoAlignUsevar.get(), self.edgestepvar.get(), self.Filtervar.get(), self.bf_scroll_scale.get(), self.eminvar.get(), self.emaxxanesvar.get(), self.emaxexafsvar.get(), self.constkUsevar.get(), self.kstepvar.get(), self.estepvar.get(), str(self.sam_numer.get()), str(self.sam_denom.get()), int(self.sam_logvar.get()), str(self.ref_numer.get()), str(self.ref_denom.get()), int(self.ref_logvar.get()), int(self.Flatvar.get()), self.updownvar.get(), self.edgestepwidth, self.BlackmanHarrisFiltervar.get(), self.blackmanfilterwindowvar.get()]
		   
			sub_f = 'batch_extract_subroutine_v2.9.1.py'
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
				if 'ID' in stdoutdata.decode('utf-8'):
					print(stdoutdata.decode('utf-8').split('\n')[0])
				#print(stdoutdata.decode('utf-8').split('\n')[0])

		end_loop = time.time()
		print('Total Subprocess Time = ', end_loop - start_loop)	
		
		i=0
		none_vals = []
		
		self.tpstextvar.set(' ')
		self.update()
		
		if self.ProcessSamUsevar.get() == 1:
			df_sam = pd.DataFrame()
			if '.bin' in self.file_type:
				df_sam_error = pd.DataFrame()
		if self.ProcessRefUsevar.get() == 1:
			df_ref = pd.DataFrame()
			if '.bin' in self.file_type:
				df_ref_error = pd.DataFrame()
		
		for i in range(int(len(self.calibroots)-5)):
			try:
				if self.updownvar.get() == 0:
					data_read = pd.read_csv(self.folder+'/Export/Individual_Up/'+str(int(i))+'.dat', sep = '\t', header=0)
				if self.updownvar.get() == 1:
					data_read = pd.read_csv(self.folder+'/Export/Individual_Down/'+str(int(i))+'.dat', sep = '\t', header=0)
				if self.updownvar.get() == 2:
					data_read = pd.read_csv(self.folder+'/Export/Individual_Both/'+str(int(i))+'.dat', sep = '\t', header=0)
				if self.ProcessSamUsevar.get() == 1:
					df_sam[str(i)] = data_read['mu_sam']
					#if '.bin' in self.file_type:
					#	df_sam_error[str(i)] = data_read['std_sam']
				if self.ProcessRefUsevar.get() == 1:
					df_ref[str(i)] = data_read['mu_ref']
					#if '.bin' in self.file_type:
					#	df_ref_error[str(i)] = data_read['std_ref']
				i = i+1
			except:
				none_vals.append(i)
				i = i+1
			self.progress_bar["value"] = 100*(i+1)/(len(self.calibroots)-5)
			self.progress_bar.update() 
			self.update()
		
		self.labelText.set('Saving Matrices')
		self.update()		
		
		if int(self.updownvar.get()) == 0:
			updownstring = 'Up'
		elif int(self.updownvar.get()) == 1:
			updownstring = 'Down'
		else:
			updownstring = 'Both'
		
		try:
			if self.ProcessSamUsevar.get() == 1:
				df_sam.insert(0, 'E', data_read['E'])
				df_sam.to_csv(self.folder+'/Export/'+self.data_file.split('/')[-1]+'_sam_matrix_'+updownstring+'.dat', sep='\t', index=False)
				#if '.bin' in self.file_type:
				#	df_sam_error.insert(0, 'E', data_read['E'])
				#	df_sam_error.to_csv(self.folder+'/Export/'+self.data_file.split('/')[-1]+'_sam_matrix_'+updownstring+'.error', sep='\t', index=False)
				print('Sample matrix saved')
		except:
			print('failed to save sample matrix')
		
		try:            
			if self.ProcessRefUsevar.get() == 1:
				df_ref.insert(0, 'E', data_read['E'])
				df_ref.to_csv(self.folder+'/Export/'+self.data_file.split('/')[-1]+'_ref_matrix_'+updownstring+'.dat', sep='\t', index=False)
				#if '.bin' in self.file_type:
				#	df_ref_error.insert(0, 'E', data_read['E'])
				#	df_ref_error.to_csv(self.folder+'/Export/'+self.data_file.split('/')[-1]+'_ref_matrix_'+updownstring+'.error', sep='\t', index=False)
				print('Reference matrix saved')
				del df_ref
		except:
			print('failed to save reference matrix')
			   
		try:
			if self.ProcessSamUsevar.get() == 1:
				if self.ExportMCRUsevar.get() == 1:
					df_sam.drop(df_sam[df_sam['E'] > float(self.emaxxanesvar.get())].index, inplace=True)
					df_sam.drop(columns=['E'], inplace=True)
					df_sam = df_sam.T
					df_sam.to_csv(self.folder+'/Export/'+self.data_file.split('/')[-1]+'_MCR_matrix_'+updownstring+'.dat', sep='\t', index=False, header=None)
					print('MCR matrix saved')
				del df_sam
		except:
			print('failed to save MCR matrix')
			del df_sam
				
		main_f = str(os.path.basename(main.__file__).split('.')[0])
		
		outf = open(self.folder+'/Export/parameters.txt','w')
		outf.write(str(main_f)+'\n')
		outf.write('file name '+str(self.data_file)+'\n')
		outf.write('folder '+str(self.folder)+'\n')
		outf.write('\n')
		outf.write('Calibration Energy '+str(self.calibenervar.get())+'\n')
		if (self.diff_calib == True) and (self.NormUsevar.get() == 0):
			outf.write('Sample Edge Energy '+str(self.e0_scroll_scale.get())+'\n')
		outf.write('\n')
		if self.NormUsevar.get() == 1:
			outf.write('Processed with Normalisation'+'\n')
			outf.write('Sample Edge Energy '+str(self.e0_scroll_scale.get())+'\n')
			outf.write('Edge Jump '+str(self.edge_jump)+'\n')
			outf.write('Pre-edge 1 '+str(self.pre1)+'\n')
			outf.write('Pre-edge 2 '+str(self.pre2)+'\n')
			outf.write('order pre-edge '+str(self.order_pre)+'\n')
			outf.write('Post-edge 1 '+str(self.post1)+'\n')
			outf.write('Post-edge 2 '+str(self.post2)+'\n')
			outf.write('order post-edge '+str(self.order_post)+'\n')
			outf.write('\n')
		if self.AutoAlignUsevar.get() == 1:
			outf.write('Processed with Automatic Reference Alignment'+'\n')
			outf.write('Alignment Energy '+str(self.edgestepvar.get())+'\n')
			outf.write('\n')
		if self.InterpUsevar.get() == 1:
			outf.write('Processed Using Radial Basis Function Interpolation'+'\n')
			outf.write('Estep '+str(self.estepvar.get())+'\n')
			if self.constkUsevar.get() == 1:	
				outf.write('Kstep '+str(self.kstepvar.get())+'\n')
			outf.write('Emin '+str(self.eminvar.get())+'\n')
			outf.write('Emax XANES (eV) '+str(self.emaxxanesvar.get())+'\n')
			outf.write('Emax EXAFS (eV) '+str(self.emaxexafsvar.get())+'\n')
			outf.write('\n')
		if self.Filtervar.get() == 1:
			outf.write('Processed Using Butterworth Noise Filtering'+'\n')
			outf.write('Butterworth Filter value '+str(self.calibration_values[2][0]))
		outf.close()
		
		gc.collect()
		print('Export Complete')
		end_batch = time.time()
		print('Total Process Time = ', end_batch - start_batch)
		self.labelText.set('Time estimate')
		self.progress_bar["value"] = 0
		self.progress_bar.update()
		
################################################################################
	def cumgauss(self, x, mu, sigma):
		return 0.5 * (1 + special.erf((x-mu)/(np.sqrt(2)*sigma)))	
		
#################################################################################			
##POSTPROCESSTAB			
#################################################################################	
class PostProcess(tkinter.Frame):
	def __init__(self,name,*args,**kwargs):
		self.frame_dataread = tkinter.Frame.__init__(self,*args,**kwargs)
				
		self.screen_width = self.winfo_screenwidth()
		self.screen_height = self.winfo_screenheight()
		
		self.progress_bar = ttk.Progressbar(self, orient='horizontal', mode='determinate', maximum=100, style="red.Horizontal.TProgressbar")
		self.progress_bar.grid(column=0, row=99, columnspan=99)
		self.progress_bar['length'] = 1000
		
		self.labelText = tkinter.StringVar()
		self.labelText.set('Time estimate')
		self.labeltime = tkinter.Label(self, textvariable=self.labelText)
		self.labeltime.grid(column=1, row=98, columnspan=2, rowspan=1, sticky='W')
		
		self.tpstextvar = tkinter.StringVar()
		self.tpstextvar.set('Time Per Spectrum')
		self.labeltps = tkinter.Label(self, textvariable=self.tpstextvar)
		self.labeltps.grid(column=3, row=98, columnspan=1, rowspan=1, sticky='W')
		
		self.frame1 = tkinter.LabelFrame(self, width=0.25*0.8333*self.screen_height, height=0.75*self.screen_height)
		self.frame2 = tkinter.LabelFrame(self, width=0.5*0.8333*self.screen_height, height=0.7*0.75*self.screen_height)
		self.frame3 = tkinter.LabelFrame(self, width=0.25*0.8333*self.screen_height, height=0.3*0.75*self.screen_height, borderwidth = 0, highlightthickness = 0)
		self.frame4 = tkinter.LabelFrame(self, width=0.25*0.8333*self.screen_height, height=0.75*self.screen_height, borderwidth = 0, highlightthickness = 0)
		#self.frame5 = tkinter.LabelFrame(self, width=0.25*0.8333*self.screen_height, height=0.3*0.75*self.screen_height)
		self.frame1.grid(row=0, column=0, rowspan=4, columnspan=1, padx=2, sticky='N')
		self.frame2.grid(row=0, column=1, rowspan=2, columnspan=2, padx=2, pady=15)
		self.frame3.grid(row=2, column=1, rowspan=2, columnspan=1, padx=2, pady=15, sticky='N')
		self.frame4.grid(row=0, column=3, rowspan=4, columnspan=1, padx=2, pady=15, sticky='N')
		#self.frame5.grid(row=2, column=2, rowspan=2, columnspan=1, padx=2, pady=15)
		
################################################################################
#FRAME1
################################################################################		
		labelfile = tkinter.Label(self.frame1, text='Load Data',anchor='w')
		labelfile.grid(column=0, row=0, columnspan=3, rowspan=1, sticky='S')
		
		self.buttonplus= tkinter.Button(self.frame1, text='+', width=4, height=1,command=self.fileplus)
		self.buttonplus.grid(column=0, row=1, columnspan=1, sticky='E')
		
		self.buttonminus= tkinter.Button(self.frame1, text='-', width=4, height=1,command=self.fileminus)
		self.buttonminus.grid(column=1, row=1, columnspan=1, sticky='W') 
		
		self.listbox2 = tkinter.Listbox(self.frame1)
		self.listbox2.grid(column=0, row=2, columnspan = 2)
		self.listbox2['width'] = 50
		self.listbox2['height'] = 53
		self.listbox2.bind('<<ListboxSelect>>', self.onselect)
		
		self.files_long_list = []
		self.edge_info = []
		self.edge_type = []
		self.edge_jump = 0
		self.import_type = []
		
################################################################################
#FRAME2
################################################################################
		def FigureCanvas(self):
			self.fig = plt.figure.Figure(figsize=(5.75, 5.75))
			
			self.canvasFrame = tkinter.Frame(master=self.frame2, padx=16, background="white")
			self.canvasFrame.grid(row=0,column=0)
			self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
			self.canvas.get_tk_widget().grid(column=0, row=0, sticky='EW')
			
			self.ax=self.fig.add_subplot(111)
			self.canvas.draw_idle()
			
			self.toolbarFrame = tkinter.Frame(master=self.frame2)
			self.toolbarFrame.grid(row=1,column=0, sticky='W')
			self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
			self.toolbar.update()
			
			self.axis_color = 'lightgoldenrodyellow'	
			self.update()
		
		FigureCanvas(self)
		
################################################################################	
#FRAME3
################################################################################
		def raise_frame(frame):
			frame.tkraise()
			
		self.d1 = tkinter.Frame(master=self.frame3)
		self.c1 = tkinter.Frame(master=self.frame3)
		
		for frame in (self.d1, self.c1):
			frame.grid(row=0, column=0, sticky='news')
			
################################################################################			
		self.scrollFrame = tkinter.Frame(master=self.d1)
		self.scrollFrame.grid(row=0, column=0, rowspan=1, columnspan=2, pady=5)
		tkinter.Label(self.scrollFrame, width=15).grid(column=0, row=0, columnspan=1)
		tkinter.Label(self.scrollFrame, width=83).grid(column=1, row=0, columnspan=2)

		self.ds_scroll_scale = tkinter.Scale(self.scrollFrame, orient='horizontal', from_=0, to=1, state='disabled', width=18)
		self.ds_scroll_scale.grid(column=1, columnspan=2, row=1, stick='EW')
		self.datasetlabel = tkinter.Label(self.scrollFrame, text='Dataset',anchor='w')
		self.datasetlabel.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='S')
		
		self.bf_scroll_scale = tkinter.Scale(self.scrollFrame, orient='horizontal', from_=0.0005, to=0.5, digits = 4, resolution = 0.0005, state='normal', width=18)
		#self.bf_scroll_scale.grid(column=1, columnspan=1, row=1, stick='EW')
		self.bflabel = tkinter.Label(self.scrollFrame, text='Filter',anchor='w')
		#self.bflabel.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='S')
		
		raise_frame(self.d1)
		
################################################################################
#FRAME4
################################################################################		
			
		self.f1 = tkinter.Frame(master=self.frame4)
		f2 = tkinter.Frame(master=self.frame4)
		
		for frame in (self.f1, f2):
			frame.grid(row=0, column=0, sticky='news')
		
################################################################################
#FRAME4-1
################################################################################

################################################################################
#Average Frame		
		self.AverageFrame = tkinter.LabelFrame(master=self.f1, width=300, height=0.75*self.screen_height, pady=5, padx=73)
		self.AverageFrame.grid(row=1, column=0) 
		
		averageuselabel = tkinter.Label(self.AverageFrame, text='Use Averaged Data',anchor='w')
		averageuselabel.grid(column=0, row=1, columnspan=2, rowspan=1, sticky='E') 
		self.AverageUsevar = tkinter.DoubleVar()
		self.AverageUsevar.set(0)
		self.AverageUsevarCheck = tkinter.Checkbutton(self.AverageFrame, variable=self.AverageUsevar, command=self.plot_update_trigger)
		self.AverageUsevarCheck.grid(column=2, row=1)
		
		NSUM = tkinter.Label(self.AverageFrame, text='Nsum',anchor='w')
		NSUM.grid(column=0, row=2, columnspan=2, rowspan=1, sticky='N')
		self.nsumvar = tkinter.StringVar()
		self.nsumvar.set('5')
		self.nsumentry = tkinter.Entry(self.AverageFrame, textvariable=self.nsumvar)
		self.nsumentry.grid(column=2, row=2, columnspan=1)  
		
		Nstart = tkinter.Label(self.AverageFrame, text='N start',anchor='w')
		Nstart.grid(column=0, row=3, columnspan=2, rowspan=1, sticky='N')
		self.nstartvar = tkinter.StringVar()
		self.nstartentry = tkinter.Entry(self.AverageFrame, textvariable=self.nstartvar)
		self.nstartentry.grid(column=2, row=3, columnspan=1) 

		Nend = tkinter.Label(self.AverageFrame, text='N end',anchor='w')
		Nend.grid(column=0, row=4, columnspan=2, rowspan=1, sticky='N')
		self.nendvar = tkinter.StringVar()
		self.nendentry = tkinter.Entry(self.AverageFrame, textvariable=self.nendvar)
		self.nendentry.grid(column=2, row=4, columnspan=1)  		
		
		Estart = tkinter.Label(self.AverageFrame, text='Emin (eV)',anchor='w')
		Estart.grid(column=0, row=5, columnspan=2, rowspan=1, sticky='N')
		self.estartvar = tkinter.StringVar()
		self.estartentry = tkinter.Entry(self.AverageFrame, textvariable=self.estartvar)
		self.estartentry.grid(column=2, row=5, columnspan=1)  	
		
		Eend = tkinter.Label(self.AverageFrame, text='Emax (eV)',anchor='w')
		Eend.grid(column=0, row=6, columnspan=2, rowspan=1, sticky='N')
		self.eendvar = tkinter.StringVar()
		self.eendentry = tkinter.Entry(self.AverageFrame, textvariable=self.eendvar)
		self.eendentry.grid(column=2, row=6, columnspan=1)  	
		
		self.buttonaverage= tkinter.Button(self.AverageFrame, text='Average', width=30, height=1, command=self.average_trigger)
		self.buttonaverage.grid(column=0, row=7, columnspan=3, rowspan=1)	
		
		self.saveaverage= tkinter.Button(self.AverageFrame, text='Save Average', width=30, height=1, command=self.saveaverage_trigger)
		self.saveaverage.grid(column=0, row=8, columnspan=3, rowspan=1)	
	
################################################################################
	
		self.active = False
		self.active_val = 0
		self.trigger = 0
		self.initialize_filter_scroll = 0
		
################################################################################
#Normalization Frame	
		self.NormalisationFrame = tkinter.LabelFrame(master=self.f1, width=300, height=0.75*self.screen_height, pady=5, padx=15)
		
		self.buttonnorm= tkinter.Button(self.NormalisationFrame, text='Normalise', width=30, height=1,command=self.normalize)
		self.buttonnorm.grid(column=0, row=1, columnspan=2)
        
		normuselabel = tkinter.Label(self.NormalisationFrame, text='Use Normalization',anchor='w')
		normuselabel.grid(column=0, row=2, columnspan=2, rowspan=1, sticky='E') 
		self.NormUsevar = tkinter.DoubleVar()
		self.NormUsevar.set(0)
		self.NormUseCheck = tkinter.Checkbutton(self.NormalisationFrame, variable=self.NormUsevar, command=self.plot_update_trigger)
		self.NormUseCheck.grid(column=2, row=2)
		
		Flatlabel = tkinter.Label(self.NormalisationFrame, text='Flatten',anchor='w')
		Flatlabel.grid(column=0, row=3, columnspan=2, rowspan=1, sticky='E') 
		self.Flatvar = tkinter.DoubleVar()
		self.Flatvar.set(1)
		self.FlatCheck = tkinter.Checkbutton(self.NormalisationFrame, variable=self.Flatvar, command=self.plot_update_trigger)
		self.FlatCheck.grid(column=2, row=3)
		
		ejlabel = tkinter.Label(self.NormalisationFrame, text='Edge Jump',anchor='w')
		ejlabel.grid(column=0, row=10, columnspan=1, rowspan=1, sticky='E')
		
		self.evjlabel = tkinter.StringVar()
		self.evjlabel.set('0')
		self.ejv = tkinter.Label(self.NormalisationFrame, textvariable=self.evjlabel,anchor='w')
		self.ejv.grid(column=1, row=10, columnspan=2, rowspan=1, sticky='E') 		
		
		self.buttonloadnorm= tkinter.Button(self.NormalisationFrame, text='Load Normalization', width=14, height=1,command=self.normload)
		self.buttonloadnorm.grid(column=2, row=1, columnspan=1)
		
		self.norm_pre = tkinter.StringVar(self.NormalisationFrame)
		self.norm_pre.set('1')
		self.norm_preMenu = tkinter.OptionMenu(self.NormalisationFrame, self.norm_pre, '0','1','2','3','4','V', command=self.plot_update_trigger2)
		self.norm_preMenu.grid(column=2, row=5, columnspan=1) 
        
		self.norm_post = tkinter.StringVar(self.NormalisationFrame)
		self.norm_post.set('V')
		self.norm_postMenu = tkinter.OptionMenu(self.NormalisationFrame, self.norm_post, '0','1','2','3','4','V', command=self.plot_update_trigger2)
		self.norm_postMenu.grid(column=2, row=7, columnspan=1) 
		
		self.pre1_scroll_scale = tkinter.Scale(self.NormalisationFrame, orient='horizontal', from_=-150, to=0, state='disabled', digits=4, resolution=0.1, length=125)
		self.pre1_scroll_scale.grid(column=1, columnspan=1, row=5, stick='EW')
		self.pre1_scroll_scale.set(-100)
		self.pre2_scroll_scale = tkinter.Scale(self.NormalisationFrame, orient='horizontal', from_=-150, to=0, state='disabled', digits=4, resolution=0.1, length=125)
		self.pre2_scroll_scale.grid(column=1, columnspan=1, row=6, stick='EW')
		self.pre2_scroll_scale.set(-50)
		self.post1_scroll_scale = tkinter.Scale(self.NormalisationFrame, orient='horizontal', from_=0, to=1000, state='disabled', digits=4, resolution=0.1, length=125)
		self.post1_scroll_scale.grid(column=1, columnspan=1, row=7, stick='EW')
		self.post1_scroll_scale.set(100)
		self.post2_scroll_scale = tkinter.Scale(self.NormalisationFrame, orient='horizontal', from_=0, to=1000, state='disabled', digits=5, resolution=0.1, length=125)
		self.post2_scroll_scale.grid(column=1, columnspan=1, row=8, stick='EW')
		self.post2_scroll_scale.set(1000)
		
		self.e0_scroll_scale = tkinter.DoubleVar()
		self.e0_entry = tkinter.Entry(self.NormalisationFrame, textvariable=self.e0_scroll_scale, state='normal')
		self.e0_entry.grid(column=1, columnspan=1, row=9, stick='EW')
		self.e0_scroll_scale.set(float(0))
		
		pre1label = tkinter.Label(self.NormalisationFrame, text='Pre1',anchor='w')
		pre1label.grid(column=0, row=5, columnspan=1, rowspan=1, sticky='S')
		pre2label = tkinter.Label(self.NormalisationFrame, text='Pre2',anchor='w')
		pre2label.grid(column=0, row=6, columnspan=1, rowspan=1, sticky='S')
		post1label = tkinter.Label(self.NormalisationFrame, text='Post1',anchor='w')
		post1label.grid(column=0, row=7, columnspan=1, rowspan=1, sticky='S')
		post2label = tkinter.Label(self.NormalisationFrame, text='Post2',anchor='w')
		post2label.grid(column=0, row=8, columnspan=1, rowspan=1, sticky='S')
		e0label = tkinter.Label(self.NormalisationFrame, text='E0',anchor='w')
		e0label.grid(column=0, row=9, columnspan=1, rowspan=1, sticky='S')
		
		self.norm_accept_button=tkinter.Button(self.NormalisationFrame, text='Accept', state='disabled', width=14, command=self.accept_norm).grid(column=1, row=11, stick='E')
		self.norm_cancel_button=tkinter.Button(self.NormalisationFrame, text='Cancel', state='disabled', width=14, command=self.cancel_norm).grid(column=2, row=11, stick='W')
		
		self.update_e0_button=tkinter.Button(self.NormalisationFrame, text='Update E0').grid(column=2, row=9)
		
		self.save_norm_data_button=tkinter.Button(self.NormalisationFrame, text='Save Normalised Data', state='normal', width=28, command=self.save_norm_data).grid(column=1, row=12, columnspan = 2, stick='EW')
		
		mcrnormlabel = tkinter.Label(self.NormalisationFrame, text='Save MCR',anchor='w')
		mcrnormlabel.grid(column=0, row=13, columnspan=1, rowspan=1, sticky='E') 
		self.mcrnormvar = tkinter.DoubleVar()
		self.mcrnormvar.set(0)
		self.mcrnormCheck = tkinter.Checkbutton(self.NormalisationFrame, variable=self.mcrnormvar)
		self.mcrnormCheck.grid(column=1, row=13)
		
		self.emin_mcr = tkinter.DoubleVar()
		self.emin_mcr_entry = tkinter.Entry(self.NormalisationFrame, textvariable=self.emin_mcr, state='normal')
		self.emin_mcr_entry.grid(column=1, columnspan=1, row=14, stick='EW')
		
		self.emax_mcr = tkinter.DoubleVar()
		self.emax_mcr_entry = tkinter.Entry(self.NormalisationFrame, textvariable=self.emax_mcr, state='normal')
		self.emax_mcr_entry.grid(column=2, columnspan=1, row=14, stick='EW')
	
################################################################################
#Time Filter Frame		
		self.TimeFrame = tkinter.LabelFrame(master=self.f1, width=300, height=0.75*self.screen_height, pady=5, padx=75)
		self.TimeFrame.grid(row=2, column=0) 
		
		timefilteruselabel = tkinter.Label(self.TimeFrame, text='Use Time Filtering',anchor='w')
		timefilteruselabel.grid(column=0, row=1, columnspan=2, rowspan=1, sticky='E') 
		self.TimeFilterUsevar = tkinter.DoubleVar()
		self.TimeFilterUsevar.set(0)
		self.TimeFilterUsevarCheck = tkinter.Checkbutton(self.TimeFrame, variable=self.TimeFilterUsevar, command=self.plot_update_trigger)
		self.TimeFilterUsevarCheck.grid(column=2, row=1)
		
		self.windowlabel = tkinter.Label(self.TimeFrame, text='Window (odd)',anchor='w')
		self.windowlabel.grid(column=0, row=2, columnspan=2, rowspan=1, sticky='N')
		self.windowvar = tkinter.DoubleVar()
		self.windowvar.set('5')
		self.windowentry = tkinter.Entry(self.TimeFrame, textvariable=self.windowvar)
		self.windowentry.grid(column=2, row=2, columnspan=1)  
		
		self.orderlabel = tkinter.Label(self.TimeFrame, text='Order',anchor='w')
		self.orderlabel.grid(column=0, row=3, columnspan=2, rowspan=1, sticky='N')
		self.ordervar = tkinter.DoubleVar()
		self.ordervar.set('2')
		self.orderentry = tkinter.Entry(self.TimeFrame, textvariable=self.ordervar)
		self.orderentry.grid(column=2, row=3, columnspan=1)  	
		
################################################################################
#Region Exclude		
		self.CutFrame = tkinter.LabelFrame(master=self.f1, width=300, height=0.75*self.screen_height, pady=5, padx=77)
		self.CutFrame.grid(row=3, column=0) 
		
		cutuselabel = tkinter.Label(self.CutFrame, text='Exclude Region',anchor='w')
		cutuselabel.grid(column=0, row=1, columnspan=2, rowspan=1, sticky='E') 
		self.cutusevar = tkinter.DoubleVar()
		self.cutusevar.set(0)
		self.cutusevarCheck = tkinter.Checkbutton(self.CutFrame, variable=self.cutusevar, command=self.plot_update_trigger)
		self.cutusevarCheck.grid(column=2, row=1)
		
		eminlabel = tkinter.Label(self.CutFrame, text='Emin (eV)',anchor='w')
		eminlabel.grid(column=0, row=2, columnspan=1, rowspan=1, sticky='E') 
		self.emin_entry = tkinter.DoubleVar()
		self.emin_entry_box = tkinter.Entry(self.CutFrame, textvariable=self.emin_entry, state='normal')
		self.emin_entry_box.grid(column=1, columnspan=1, row=2, stick='EW')
		
		emaxlabel = tkinter.Label(self.CutFrame, text='Emax (eV)',anchor='w')
		emaxlabel.grid(column=0, row=3, columnspan=1, rowspan=1, sticky='E') 
		self.emax_entry = tkinter.DoubleVar()
		self.emax_entry_box = tkinter.Entry(self.CutFrame, textvariable=self.emax_entry, state='normal')
		self.emax_entry_box.grid(column=1, columnspan=1, row=3, stick='EW')
		
		self.cut= tkinter.Button(self.CutFrame, text='Save Cut Data', width=30, height=1, command=self.savecut_trigger)
		self.cut.grid(column=0, row=4, columnspan=3, rowspan=1)	

################################################################################	
#FRAME4-2
################################################################################
		tkinter.Button(f2, text='Go to DataRead', command=lambda:raise_frame(self.f1)).grid(row=0, column=0)	

		self.PadFrame = tkinter.LabelFrame(master=f2, width=300, height=0.75*self.screen_height, pady=100, padx=166)
		self.PadFrame.grid(row=5, column=0)
		
		padlabel = tkinter.Label(self.PadFrame, text='    ',anchor='w')
		padlabel.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='N')	
		
################################################################################

		raise_frame(self.f1)
			
################################################################################
	def normfile_check(self, folder):
		try:
			self.normalisation_values = (pd.read_csv(folder+'/normalisation.dat', sep='\t', header=None)).values
			return True
		except:
			#print('Cannot Find Normalisation Parameter File')
			self.NormUsevar.set(0)
			return False
			
################################################################################
	def fileminus(self):        
		i = self.listbox2.curselection()[0]
		self.listbox2.delete(i)
		del self.files_long_list[i]
		self.update()
		
################################################################################        
	def fileplus(self):
		self.askloadfile()
		
################################################################################
	def askloadfile(self):
		file = askopenfile()
		try:
			dirname, file_short = os.path.split(os.path.abspath(file.name))
			self.files_long_list.append(file)
			print('Loading File')
			self.listbox2.insert(END, file_short)
		except:
			print('Aborted Load')
			
################################################################################					
	def clear_figure(self):
		self.fig.clear()
		self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
		self.canvas.get_tk_widget().grid(column=0, row=0, sticky='EW')
		
		self.ax=self.fig.add_subplot(111)
		self.toolbar.destroy()
		self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
		self.toolbar.update()
		
################################################################################				
	def update_figure(self, labelx, labely):
		self.clear_figure()
		self.ax.set_xlabel(labelx)
		self.ax.set_ylabel(labely)
		
#################################################################################			
	def onselect(self,evt):			
		if 'Post' in App.get_selected(self):
			data_file_long = self.files_long_list[self.listbox2.curselection()[0]].name
			print(data_file_long)
			self.folder, self.data_file = os.path.split(os.path.abspath(data_file_long))
			self.file_type=(os.path.splitext(self.data_file)[-1])
			print('File Type :', self.file_type)
			
			try:
				del self.data_sum
				del self.column_names_sum
			except:
				pass
				
			self.normalized_check = self.normfile_check(self.folder)
				
			if '.dat' in self.file_type:
				print('Loading file ... ')
				self.data = pd.read_csv(data_file_long, sep='\t', header=0)
				try:
					self.data_error = pd.read_csv(data_file_long.split('.dat')[0]+'.error', sep='\t', header=0)
					self.errors_flag = True
				except:
					self.errors_flag = False
				print('File loaded')
				self.column_names = self.data.columns.values.tolist()
				self.trigger = 0
				
				try:
					parameters_file = open(self.folder+'/parameters.txt', 'r')
				except:
					try:
						parameters_file = open(self.folder+'/parameters.dat', 'r')
					except:
						print('No parameters file detected')
						
				self.filtered = 0
				self.normalised = 0
				self.interpolated = 0
				self.filter_value = 1
				self.e0_value_enter = 0
				self.emax_xanes = 0
				
				for line in parameters_file:
					if 'Butterworth Filter' in line:
						self.filter_value = line.split(' ')[-1]
						self.filtered = 1
					if 'Processed with Normalisation' in line:
						self.normalised = 1
						self.NormalisationFrame_invisible()
					if 'Processed Using Radial Basis Function Interpolation' in line:
						self.interpolated = 1
					if 'Sample Edge Energy' in line:
						try:
							e0 = float(ast.literal_eval(line.split(' ')[-1])[0])
						except:
							try:
								e0 = float(ast.literal_eval(line.split(' ')[-1]))
							except:
								print('Cannot determine e0 value from file, setting to 0')
								e0 = 0
						print('E0 energy :',e0)
						self.e0_scroll_scale.set(float(e0))
						self.e0_value_enter = 1
					if 'Calibration Energy' in line:
						try:
							calib_ener = float(ast.literal_eval(line.split(' ')[-1])[0])
						except:
							try:
								calib_ener = float(ast.literal_eval(line.split(' ')[-1]))
							except:
								print('Cannot determine calibration value from file, setting to 0')
								calib_ener = 0
						print('Calibation Energy :',calib_ener)
					if 'Emax XANES' in line:
						self.emax_xanes = float(ast.literal_eval(line.split(' ')[-1]))
				
				if self.e0_value_enter == 0:
					self.e0_scroll_scale.set(float(calib_ener))
				
				self.update()
				
				if any('normalised' in s for s in self.data_file.split('_')):
					self.normalised = 1
					
				if self.normalised == 1:
					self.NormalisationFrame_invisible()
	
				if (self.normalised == 0):
					self.NormalisationFrame_visible()
					
				self.plot_spectrum()
			
#################################################################################		
	def average_trigger(self):
		self.data_sum = pd.DataFrame()
		startval = False
		
		if self.errors_flag == True:
			self.data_sum_error = pd.DataFrame()
		idxs = list(map(int, self.column_names[1:]))

		while float(idxs[0]) != 0:
			idxs = list(map(lambda x: x - 1, idxs))
		
		if len(self.nstartvar.get()) > 0:
			if int(idxs[1]) - int(idxs[0]) == 2:
				idxs = list(np.asarray(idxs)[(np.asarray(idxs) >= 2*float(self.nstartvar.get()))])
				print('Using requested start dataset for averaging')
				startval = True
			else:
				idxs = list(np.asarray(idxs)[(np.asarray(idxs) >= float(self.nstartvar.get()))])
				print('Using requested start dataset for averaging')
				startval = True
		if len(self.nendvar.get()) > 0:
			if int(idxs[1]) - int(idxs[0]) == 2:
				idxs = list(np.asarray(idxs)[(np.asarray(idxs) <= 2*float(self.nendvar.get()))])
				print('Using requested end dataset for averaging')
			else:
				idxs = list(np.asarray(idxs)[(np.asarray(idxs) <= float(self.nendvar.get()))])
				print('Using requested end dataset for averaging')
			
		if self.nsumvar.get() == 'all':
			if int(idxs[1]) - int(idxs[0]) == 2:
				Nsum = int(np.max(idxs)-np.min(idxs)+2)
				self.sumlen = int(Nsum/2)
			else:
				Nsum = np.max(idxs)-np.min(idxs)+1
				self.sumlen = Nsum
		else:
			if int(idxs[1]) - int(idxs[0]) == 2:
				Nsum = 2*int(self.nsumvar.get())
				self.sumlen = int(Nsum/2)
			else:
				Nsum = int(self.nsumvar.get())
				self.sumlen = int(Nsum)
		
		if int(idxs[1]) - int(idxs[0]) == 2:
			Nave = int(np.floor((np.max(idxs)-np.min(idxs)+2)/Nsum))
		else:
			Nave = int((np.max(idxs)-np.min(idxs)+1)/Nsum)

		if Nave > 0:
			print('Number of Averaged Spectra =',Nave)
			for j in range(int(Nave)):
				if startval == True:
					if int(idxs[1]) - int(idxs[0]) == 2:
						idx = int(self.nstartvar.get()) + npi.indices(np.asarray(idxs), np.asarray(idxs)[(np.asarray(idxs) >= ((j*Nsum)+2*int(self.nstartvar.get()))) & (np.asarray(idxs) < ((j*Nsum)+Nsum+2*int(self.nstartvar.get())))])
						print(j, idx)
					else:
						idx = int(self.nstartvar.get()) + npi.indices(np.asarray(idxs), np.asarray(idxs)[(np.asarray(idxs) >= ((j*Nsum)+int(self.nstartvar.get()))) & (np.asarray(idxs) < ((j*Nsum)+Nsum+int(self.nstartvar.get())))])
						print(j, idx)
				else:
					idx = npi.indices(np.asarray(idxs), np.asarray(idxs)[(np.asarray(idxs) >= (j*Nsum)) & (np.asarray(idxs) < (j*Nsum)+Nsum)])
					print(j, idx)
				if len(idx) > 0:
					self.data_sum[str(j)] = (self.data.iloc[:,idx+1]).mean(axis=1)
					if self.errors_flag == True:
						self.data_sum_error[str(j)] = (self.data_error.iloc[:,idx+1]).mean(axis=1)/np.sqrt(len(idx))
			
				self.progress_bar["value"] = 100*(j+1)/Nave
				self.progress_bar.update()	
					
			if self.TimeFilterUsevar.get() == 1:
				if self.windowvar.get() % 2 == 0:
					self.windowvar.set(self.windowvar.get() + 1)
					self.update()
					
				if self.windowvar.get() < self.ordervar.get():
					print('Window must be greaters than order')
					print('Setting window to :', self.odervar.get()+1)
					self.windowvar.set(self.odervar.get()+1)
					self.update()
					
				if Nave > self.windowvar.get():
					print('Applying weak time filter')
					for i in range(self.data_sum.shape[0]):
						self.data_sum.iloc[[i]] = signal.savgol_filter(self.data_sum.iloc[[i]], int(self.windowvar.get()), int(self.ordervar.get()))
				else:
					print('Not eneough spectra to smooth')
					self.TimeFilterUsevar.set(0)
					
			self.data_sum.insert(0, 'E', self.data['E'])

			if self.cutusevar.get() == 1:
				self.data_sum.drop(self.data_sum[(self.data_sum['E'] >= self.emin_entry.get()) & (self.data_sum['E'] <= self.emax_entry.get())].index, inplace=True)

			if self.errors_flag == True:
				self.data_sum_error.insert(0, 'E', self.data['E'])	
				
			print(self.estartvar.get(), self.eendvar.get())
			if (len(self.estartvar.get()) > 0) and (len(self.eendvar.get()) > 0) and (float(self.eendvar.get()) > float(self.estartvar.get())) and (float(self.estartvar.get()) < np.max(self.data['E'].values)):
				self.data_sum = self.data_sum[(self.data_sum['E'] >= float(self.estartvar.get())) & (self.data_sum['E'] <= float(self.eendvar.get()))]
			
			self.column_names_sum = self.data_sum.columns.values.tolist()
			self.progress_bar["value"] = 0
			self.progress_bar.update()   
			self.AverageUsevar.set(1)
			self.trigger = 1
			self.plot_update_trigger()
		else:
			print('Not enough spectra to average')
		
#################################################################################	
	def plot_update_trigger(self):
		self.plot_spectrum()
		
#################################################################################	
	def plot_update_trigger2(self, *args):
		self.plot_spectrum()
		
#################################################################################	
	def saveaverage_trigger(self):
		self.data_sum
		self.data_sum.to_csv(self.folder+'/'+os.path.splitext(self.data_file)[0]+'_average_'+str(self.sumlen)+'.dat', sep='\t', index=False)
		if self.errors_flag == True:
			self.data_sum_error.to_csv(self.folder+'/'+os.path.splitext(self.data_file)[0]+'_average_'+str(self.sumlen)+'.error', sep='\t', index=False)
		
		main_f = str(os.path.basename(main.__file__).split('.')[0])
		
		outf = open(self.folder+'/parameters.txt','a')
		
		outf.write('\n')
		
		if len(self.nstartvar.get()) < 1:
			start = 0
		else:
			start = self.nstartvar.get()
		if len(self.nendvar.get()) < 1:
			end = len(self.data.columns)-2
		else:
			end = self.nendvar.get()
			
		outf.write('Data averaged'+'\n')
		outf.write('N sum '+str(self.nsumvar.get())+'\n')
		outf.write('N start '+str(start)+'\n')
		outf.write('N end '+str(end)+'\n')
		
		if self.cutusevar.get() == 1:
			outf.write('\n')
			outf.write('Data region excluded'+'\n')
			outf.write('Energy min '+str(self.emin_entry.get())+'\n')
			outf.write('Energy max '+str(self.emax_entry.get())+'\n')

		outf.close()
		
		print('Average Data Saved')
			
################################################################################ 
	def plot_add_line(self, x, y):
		self.ax.plot(x, y) 
		self.canvas.draw_idle()

################################################################################ 		
	def plot_add_markerlines(self):
		ylims = self.ax.get_ylim()
		self.pre1line = self.ax.plot([self.e0_scroll_scale.get()+self.pre1_scroll_scale.get(), self.e0_scroll_scale.get()+self.pre1_scroll_scale.get()], ylims, ls='--', c='r')
		self.pre2line = self.ax.plot([self.e0_scroll_scale.get()+self.pre2_scroll_scale.get(), self.e0_scroll_scale.get()+self.pre2_scroll_scale.get()], ylims, ls='--', c='r')
		self.post1line = self.ax.plot([self.e0_scroll_scale.get()+self.post1_scroll_scale.get(), self.e0_scroll_scale.get()+self.post1_scroll_scale.get()], ylims, ls='--', c='r')
		self.post2line = self.ax.plot([self.e0_scroll_scale.get()+self.post2_scroll_scale.get(), self.e0_scroll_scale.get()+self.post2_scroll_scale.get()], ylims, ls='--', c='r')

		self.canvas.draw_idle()
			
#################################################################################			
	def plot_spectrum(self):				
		try:
			self.data_sum
		except:
			self.AverageUsevar.set(0)
			
		self.initialize_plot = 0
	
		def plot_data(evt):
			data_i = self.ds_scroll_scale.get()
			
			if self.cutusevar.get() == 1:
				if self.AverageUsevar.get() == 0:
					self.data.drop(self.data[(self.data['E'] >= self.emin_entry.get()) & (self.data['E'] <= self.emax_entry.get())].index, inplace=True)
				
			if self.AverageUsevar.get() == 1:
				self.data_x = self.data_sum[str(self.column_names_sum[0])]
				self.data_y =  self.data_sum[str(self.column_names_sum[data_i+1])]
			else: 
				self.data_x = self.data[str(self.column_names[0])]
				self.data_y =  self.data[str(self.column_names[data_i+1])]
				
			if self.NormUsevar.get() == 1:
				self.pre1,self.pre2,self.post1,self.post2,self.order_pre,self.order_post,self.e0 = self.normalisation_values
				self.data_y = self.normalize_data(np.asarray(self.data_x), np.asarray(self.data_y))
				
			if (self.trigger == 0) or (self.trigger == 1):
				if self.initialize_plot == 0:
					self.update_figure('Energy (eV)', 'Normalised Absorption')
					[self.line] = self.ax.plot(self.data_x, self.data_y)
					self.canvas.draw_idle()
					self.initialize_plot = 1
				else:
					self.line.set_ydata(self.data_y)
					self.canvas.draw_idle()
					
			if (self.trigger == 2):
			
				self.normalisation_values = [self.pre1_scroll_scale.get(), self.pre2_scroll_scale.get(),self.post1_scroll_scale.get(),self.post2_scroll_scale.get(),self.norm_pre.get(),self.norm_post.get(),self.e0_scroll_scale.get()]
				self.pre1,self.pre2,self.post1,self.post2,self.order_pre,self.order_post,self.e0 = self.normalisation_values
				
				linepre_y = self.regression('pre-edge', np.asarray(self.data_x), np.asarray(self.data_y))
				linepost_y = self.regression('post-edge', np.asarray(self.data_x), np.asarray(self.data_y))
			
				self.update_figure('Energy (eV)', 'Absorption')
				[self.line] = self.ax.plot(self.data_x, self.data_y)
				try:
					self.line_pre.set_ydata(linepre_y)
					self.line_post.set_ydata(linepost_y)
				except:
					self.line_pre = self.plot_add_line(self.data_x, linepre_y)
					self.line_post = self.plot_add_line(self.data_x, linepost_y)
			
				if self.cutusevar.get() != 1:
					index = min(range(len(self.data_x)), key=lambda i: abs(self.data_x[i]-self.e0_scroll_scale.get()))
					self.edge_jump = np.around((linepost_y-linepre_y)[index], decimals=2)
					self.evjlabel.set(str(self.edge_jump))
								
				self.plot_add_markerlines()
					
		if self.AverageUsevar.get() == 1:
			self.ds_scroll_scale = tkinter.Scale(self.scrollFrame, orient='horizontal', from_=0, to=len(self.column_names_sum)-2, state='normal', width=18, command=plot_data)
			self.ds_scroll_scale.grid(column=1, columnspan=2, row=1, stick='EW')
			self.datasetlabel = tkinter.Label(self.scrollFrame, text='Dataset',anchor='w')
			self.datasetlabel.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='S')
		else:
			self.ds_scroll_scale = tkinter.Scale(self.scrollFrame, orient='horizontal', from_=0, to=len(self.column_names)-2, state='normal', width=18, command=plot_data)
			self.ds_scroll_scale.grid(column=1, columnspan=2, row=1, stick='EW')
			self.datasetlabel = tkinter.Label(self.scrollFrame, text='Dataset',anchor='w')
			self.datasetlabel.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='S')
			
		if self.trigger == 2:		
			def scroll_bar_update(*args):

				self.pre1_scroll_scale = tkinter.Scale(self.NormalisationFrame, orient='horizontal', from_=np.min(self.data_x)-self.e0_scroll_scale.get()+5, to=0, state='normal', digits=4, resolution=0.1, length=125, command=plot_data)
				self.pre2_scroll_scale = tkinter.Scale(self.NormalisationFrame, orient='horizontal', from_=np.min(self.data_x)-self.e0_scroll_scale.get()+5, to=0, state='normal', digits=4, resolution=0.1, length=125, command=plot_data)
				self.post1_scroll_scale = tkinter.Scale(self.NormalisationFrame, orient='horizontal', from_=0, to=np.max(self.data_x)-self.e0_scroll_scale.get()-10, state='normal', digits=4, resolution=0.1, length=125, command=plot_data)
				self.post2_scroll_scale = tkinter.Scale(self.NormalisationFrame, orient='horizontal', from_=0, to=np.max(self.data_x)-self.e0_scroll_scale.get()-10, state='normal', resolution=0.1, length=125, command=plot_data)
				self.e0_entry = tkinter.Entry(self.NormalisationFrame, textvariable=self.e0_scroll_scale, state='normal')
				self.update()
				
				self.update_e0_button=tkinter.Button(self.NormalisationFrame, text='Update E0', command=callback).grid(column=2, row=9)
				
				self.pre1_scroll_scale.grid(column=1, columnspan=1, row=5, stick='EW')
				self.pre2_scroll_scale.grid(column=1, columnspan=1, row=6, stick='EW')
				self.post1_scroll_scale.grid(column=1, columnspan=1, row=7, stick='EW')
				self.post2_scroll_scale.grid(column=1, columnspan=1, row=8, stick='EW')
			
				self.pre1_scroll_scale.set(np.min(self.data_x)-self.e0_scroll_scale.get()+5)
				self.pre2_scroll_scale.set(-50)
				self.post1_scroll_scale.set(100)
				self.post2_scroll_scale.set(np.max(self.data_x)-self.e0_scroll_scale.get()-10)
				
				self.scroll_scale_init = 1
				
													
			def callback(*args):
				if self.e0_scroll_scale.get() > np.min(self.data_x):
					scroll_bar_update(None)
					plot_data(None)
			
			if self.scroll_scale_init == 0:
				scroll_bar_update(None)

		plot_data(None)
		
################################################################################				
	def NormalisationFrame_invisible(self):
		self.NormalisationFrame.grid_remove()
		
################################################################################
	def NormalisationFrame_visible(self):
		self.NormalisationFrame.grid(row=4, column=0)
		
################################################################################
	def normload(self):
		file = askopenfile()
		dirname, file_short = os.path.split(os.path.abspath(file.name))
		shutil.copyfile(file.name,self.folder+'/'+file_short)
		self.NormUsevar.set(1)
		self.Flatvar.set(1)
		self.normalisation_values = (pd.read_csv(self.folder+'/normalisation.dat', sep='\t', header=None)).values
		self.plot_spectrum() 
		
################################################################################
	def normalize(self):
		self.norm_accept_button=tkinter.Button(self.NormalisationFrame, text='Accept', state='normal', width=14, command=self.accept_norm).grid(column=1, row=11, stick='E')
		self.norm_cancel_button=tkinter.Button(self.NormalisationFrame, text='Cancel', state='normal', width=14, command=self.cancel_norm).grid(column=2, row=11, stick='W')
		
		self.NormUsevar.set(0)
		self.scroll_scale_init = 0
		self.trigger = 2
		self.plot_spectrum()
		
################################################################################                 
	def regression(self,region,x,y): 
		
		if region == 'pre-edge':
			xllim = float(self.e0) + float(self.pre1)
			xhlim = float(self.e0) + float(self.pre2)
			try:
				order = str(self.order_pre[0])
			except:
				order = str(self.order_pre)
	
		if region == 'post-edge':
			xllim = float(self.e0) + float(self.post1)
			xhlim = float(self.e0) + float(self.post2)
			try:
				order = str(self.order_post[0])
			except:
				order = str(self.order_post)
		
		xregion = x[(x >= xllim) & (x <= xhlim)]
		yregion = y[(x >= xllim) & (x <= xhlim)]
		
		if (order == '0') or (order == '0.0'):
			popt, pcov = curve_fit(polynomials.constant, xregion, yregion, p0=np.ones(1))
			yfit = polynomials.constant(x, *popt)
	
		if (order == '1') or (order == '1.0'):
			popt, pcov = curve_fit(polynomials.linear, xregion, yregion, p0=np.ones(2))
			yfit = polynomials.linear(x, *popt)
			
		if (order == '2') or (order == '2.0'):
			popt, pcov = curve_fit(polynomials.quadratic, xregion, yregion, p0=np.ones(3))
			yfit = polynomials.quadratic(x, *popt)
			
		if (order == '3') or (order == '3.0'):
			popt, pcov = curve_fit(polynomials.cubic, xregion, yregion, p0=np.ones(4))
			yfit = polynomials.cubic(x, *popt)
			
		if (order == '4') or (order == '4.0'):
			popt, pcov = curve_fit(polynomials.quartic, xregion, yregion, p0=np.ones(5))
			yfit = polynomials.quartic(x, *popt)
			
		if order == 'V':
			popt, pcov = curve_fit(polynomials.victoreen, xregion, yregion, p0=np.ones(3))
			yfit = polynomials.victoreen(x, *popt)
	
		return yfit
		
################################################################################	
	def accept_norm(self):
		self.norm_accept_button=tkinter.Button(self.NormalisationFrame, text='Accept', state='disabled', width=14, command=self.accept_norm).grid(column=1, row=11, stick='E')
		self.norm_cancel_button=tkinter.Button(self.NormalisationFrame, text='Cancel', state='disabled', width=14, command=self.cancel_norm).grid(column=2, row=11, stick='W')
		
		NormValues = [self.pre1_scroll_scale.get(), self.pre2_scroll_scale.get(),self.post1_scroll_scale.get(),self.post2_scroll_scale.get(),self.norm_pre.get(),self.norm_post.get(),self.e0_scroll_scale.get()]
		pdnormvalue = pd.DataFrame(NormValues)
		pdnormvalue.to_csv(self.folder+'/normalisation.dat', sep='\t', header=None, index=False)
		
		print(NormValues)
		
		self.trigger = 0
		self.NormUsevar.set(1)
		self.plot_spectrum()

################################################################################	
	def cancel_norm(self):
		self.norm_accept_button=tkinter.Button(self.NormalisationFrame, text='Accept', state='disabled', width=14, command=self.accept_norm).grid(column=1, row=11, stick='E')
		self.norm_cancel_button=tkinter.Button(self.NormalisationFrame, text='Cancel', state='disabled', width=14, command=self.cancel_norm).grid(column=2, row=11, stick='W')
		
		self.trigger = 0
		self.plot_spectrum()
		
################################################################################ 
	def normalize_data(self,energy, mu):  
		linepre_y = self.regression('pre-edge',np.asarray(energy), np.asarray(mu))        
		linepost_y = self.regression('post-edge',np.asarray(energy), np.asarray(mu))
		
		index = min(range(len(energy)), key=lambda i: abs(energy[i]-float(self.e0)))
		normdivisor = (linepost_y-linepre_y)[index]
		self.edge_jump = normdivisor
		
		mu = (mu-linepre_y)/normdivisor
		linepost_y = (linepost_y-linepre_y)/normdivisor
		linepre_y = np.zeros(len(linepre_y))
						
		if self.Flatvar.get() == 1:
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
#OBSOLETE
################################################################################ 
	def savecut_trigger(self):  
	
		if self.cutusevar.get() == 1:
			if self.NormUsevar.get() == 1:
				self.normalisation_values = (pd.read_csv(self.folder+'/normalisation.dat', sep='\t', header=None)).values
				self.pre1,self.pre2,self.post1,self.post2,self.order_pre,self.order_post,self.e0 = self.normalisation_values
			
			self.export = pd.DataFrame()
			
			if self.AverageUsevar.get() == 1:
				data_length = len(self.column_names_sum)-1
			else:
				data_length = len(self.column_names)-1
				
			for i in range(data_length):
				if self.AverageUsevar.get() == 1:
					self.data_x = self.data_sum[str(self.column_names_sum[0])]
					self.data_y =  self.data_sum[str(self.column_names_sum[i+1])]
				else: 
					self.data_x = self.data[str(self.column_names[0])]
					self.data_y =  self.data[str(self.column_names[i+1])]
				
				if self.NormUsevar.get() == 1:			
					self.data_y = self.normalize_data(np.asarray(self.data_x), np.asarray(self.data_y))
				
				if i == 0:
					self.export['Energy'] = self.data_x
				
				self.export[str(i)] = self.data_y
				
			self.export = self.export[(self.export['Energy'] >= np.float(self.emin_entry.get())) & (self.export['Energy'] <= np.float(self.emax_entry.get()))]
	
			if self.NormUsevar.get() == 0:
				self.export.to_csv(self.folder+'/'+os.path.splitext(self.data_file)[0]+'_cut.dat', sep='\t', index=False)
			else:
				self.export.to_csv(self.folder+'/'+os.path.splitext(self.data_file)[0]+'_normalised_cut.dat', sep='\t', index=False)
			print('Cut Data Saved')
				
################################################################################ 
	def save_norm_data(self):  
		if self.cutusevar.get() == 1:
			if self.AverageUsevar.get() == 0:
				self.data.drop(self.data[(self.data['E'] >= self.emin_entry.get()) & (self.data['E'] <= self.emax_entry.get())].index, inplace=True)
	
		if self.NormUsevar.get() == 1:
			self.normalisation_values = (pd.read_csv(self.folder+'/normalisation.dat', sep='\t', header=None)).values
			self.pre1,self.pre2,self.post1,self.post2,self.order_pre,self.order_post,self.e0 = self.normalisation_values
			
			self.norm_export = pd.DataFrame()
			if self.errors_flag == True:
				self.norm_export_error = pd.DataFrame()
			
			if self.AverageUsevar.get() == 1:
				data_length = len(self.column_names_sum)-1
			else:
				data_length = len(self.column_names)-1
				
			for i in range(data_length):
				if self.AverageUsevar.get() == 1:
					self.data_x = self.data_sum[str(self.column_names_sum[0])]
					self.data_y =  self.data_sum[str(self.column_names_sum[i+1])]
					if self.errors_flag == True:
						self.data_y_error = self.data_sum_error[str(self.column_names_sum[i+1])]
				else: 
					self.data_x = self.data[str(self.column_names[0])]
					self.data_y =  self.data[str(self.column_names[i+1])]
					if self.errors_flag == True:
						self.data_y_error = self.data_error[str(self.column_names[i+1])]
					
				self.data_y = self.normalize_data(np.asarray(self.data_x), np.asarray(self.data_y))
				if self.errors_flag == True:
					self.data_y_error = self.data_y_error * self.edge_jump
				
				if i == 0:
					self.norm_export['Energy'] = self.data_x
					if self.errors_flag == True:
						self.norm_export_errors['Energy'] = self.data_x
				
				self.norm_export[str(i)] = self.data_y
				if self.errors_flag == True:
					self.norm_export_errors[str(i)] = self.data_y_error

			if self.AverageUsevar.get() == 1:
				self.norm_export.to_csv(self.folder+'/'+os.path.splitext(self.data_file)[0]+'_normalised_'+str(self.nsumvar.get())+'.dat', sep='\t', index=False)
				if self.errors_flag == True:
					self.norm_export_error.to_csv(self.folder+'/'+os.path.splitext(self.data_file)[0]+'_normalised_'+str(self.nsumvar.get())+'.error', sep='\t', index=False)
				print('Normalised and Averaged Data Saved, Nsum :', self.nsumvar.get())
			else:
				self.norm_export.to_csv(self.folder+'/'+os.path.splitext(self.data_file)[0]+'_normalised.dat', sep='\t', index=False)
				if self.errors_flag == True:
					self.norm_export_error.to_csv(self.folder+'/'+os.path.splitext(self.data_file)[0]+'_normalised.error', sep='\t', index=False)
				print('Normalised Data Saved')
				
			try:
				parameters_file = open(self.folder+'/parameters.txt', 'r')
			except:
				try:
					parameters_file = open(self.folder+'/parameters.dat', 'r')
				except:
					print('No parameters file detected')
					
			self.saved_average = 0
			self.saved_exclude = 0
			
			for line in parameters_file:
				if 'Data averaged' in line:
					self.saved_average = 1
				if 'Data region excluded' in line:
					self.saved_exclude = 1
			parameters_file.close()
	
			main_f = str(os.path.basename(main.__file__).split('.')[0])
		
			outf = open(self.folder+'/parameters.txt','a')
			outf.write('\n')
			
			if self.saved_average == 0:
				if self.AverageUsevar.get() == 1:
					outf.write('Data averaged'+'\n')
					outf.write('N sum '+str(self.nsumvar.get())+'\n')
					outf.write('N start '+str(self.nstartvar.get())+'\n')
					outf.write('N end '+str(self.nendvar.get())+'\n')
					outf.write('\n')
			elif self.saved_average == 1:
				if self.AverageUsevar.get() == 1:
					outf.write('Data average used with normalisation save'+'\n')
					outf.write('N sum '+str(self.nsumvar.get())+'\n')
					try:
						outf.write('N start '+str(self.nstartvar.get())+'\n')
						outf.write('N end '+str(self.nendvar.get())+'\n')
						outf.write('\n')
					except:
						outf.write('\n')
					
			if self.saved_exclude == 0:
				if self.cutusevar.get() == 1:
					outf.write('Data region excluded'+'\n')
					outf.write('Energy min '+str(self.emin_entry.get())+'\n')
					outf.write('Energy max '+str(self.emax_entry.get())+'\n')
					outf.write('\n')
			elif self.savesaved_exclude == 1:
				if self.cutusevar.get() == 1:
					outf.write('Data region excluded used with normalisation save'+'\n')
					outf.write('Energy min '+str(self.emin_entry.get())+'\n')
					outf.write('Energy max '+str(self.emax_entry.get())+'\n')
					outf.write('\n')
					
			outf.write('Data Normalised in Post Processing'+'\n')
			outf.write('Sample Edge Energy '+str(self.e0)+'\n')
			outf.write('Pre-edge 1 '+str(self.pre1)+'\n')
			outf.write('Pre-edge 2 '+str(self.pre2)+'\n')
			outf.write('order pre-edge '+str(self.order_pre)+'\n')
			outf.write('Post-edge 1 '+str(self.post1)+'\n')
			outf.write('Post-edge 2 '+str(self.post2)+'\n')
			outf.write('order post-edge '+str(self.order_post)+'\n')
			outf.write('\n')

			outf.close()
				
		if self.interpolated == 1:
			if self.mcrnormvar.get() == 1:
				#self.norm_export.drop(self.norm_export[self.norm_expot['Energy'] > float(self.emax_xanes)].index, inplace=True)
				self.norm_export.drop(self.norm_export[self.norm_export['Energy'] >= self.emax_mcr.get()].index, inplace=True)
				self.norm_export.drop(self.norm_export[self.norm_export['Energy'] <= self.emin_mcr.get()].index, inplace=True)
				mcr_energy = self.norm_export['Energy']
				mcr_energy = np.round(mcr_energy,decimals=5)
				mcr_energy.to_csv(self.folder+'/mcr_energy.dat', sep='\t', index=False, header=None)
				self.norm_export.drop(columns=['Energy'], inplace=True)
				self.norm_export = self.norm_export.T
				if self.AverageUsevar.get() == 1:
					self.norm_export.to_csv(self.folder+'/'+os.path.splitext(self.data_file)[0]+'_normalised_'+str(self.nsumvar.get())+'_mcr.dat', sep='\t',  index=False, header=None)
				else:
					self.norm_export.to_csv(self.folder+'/'+os.path.splitext(self.data_file)[0]+'_normalised_mcr.dat', sep='\t',  index=False, header=None)
				print('MCR file saved')
				
			outf = open(self.folder+'/parameters.txt','a')
			
			outf.write('MCR file Saved'+'\n')
			outf.write('Energy min '+str(self.emin_mcr.get())+'\n')
			outf.write('Energy max '+str(self.emax_mcr.get())+'\n')
			
			outf.close()
				
#################################################################################			
##POSTPROCESSTAB			
#################################################################################	
class FourierTransform(tkinter.Frame):
	def __init__(self,name,*args,**kwargs):
		self.frame_dataread = tkinter.Frame.__init__(self,*args,**kwargs)
				
		self.screen_width = self.winfo_screenwidth()
		self.screen_height = self.winfo_screenheight()
		
		self.progress_bar = ttk.Progressbar(self, orient='horizontal', mode='determinate', maximum=100, style="red.Horizontal.TProgressbar")
		self.progress_bar.grid(column=0, row=99, columnspan=99)
		self.progress_bar['length'] = 1000
		
		self.labelText = tkinter.StringVar()
		self.labelText.set('Time estimate')
		self.labeltime = tkinter.Label(self, textvariable=self.labelText)
		self.labeltime.grid(column=1, row=98, columnspan=2, rowspan=1, sticky='W')
		
		self.tpstextvar = tkinter.StringVar()
		self.tpstextvar.set('Time Per Spectrum')
		self.labeltps = tkinter.Label(self, textvariable=self.tpstextvar)
		self.labeltps.grid(column=3, row=98, columnspan=1, rowspan=1, sticky='W')
		
		self.frame1 = tkinter.LabelFrame(self, width=0.25*0.8333*self.screen_height, height=0.75*self.screen_height)
		self.frame2 = tkinter.LabelFrame(self, width=0.5*0.8333*self.screen_height, height=0.7*0.75*self.screen_height)
		self.frame3 = tkinter.LabelFrame(self, width=0.25*0.8333*self.screen_height, height=0.3*0.75*self.screen_height, borderwidth = 0, highlightthickness = 0)
		self.frame4 = tkinter.LabelFrame(self, width=0.25*0.8333*self.screen_height, height=0.3*0.75*self.screen_height, borderwidth = 0, highlightthickness = 0)
		self.frame1.grid(row=0, column=0, rowspan=4, columnspan=1, padx=2, sticky='N')
		self.frame2.grid(row=0, column=1, rowspan=2, columnspan=3, padx=2, pady=15)
		self.frame3.grid(row=2, column=1, rowspan=1, columnspan=3, padx=2, pady=15)
		self.frame4.grid(row=3, column=1, rowspan=1, columnspan=3, padx=2, pady=15)
		
################################################################################
#FRAME1
################################################################################		
		labelfile = tkinter.Label(self.frame1, text='Load Data',anchor='w')
		labelfile.grid(column=0, row=0, columnspan=3, rowspan=1, sticky='S')
		
		self.buttonplus= tkinter.Button(self.frame1, text='+', width=4, height=1,command=self.fileplus)
		self.buttonplus.grid(column=0, row=1, columnspan=1, sticky='E')
		
		self.buttonminus= tkinter.Button(self.frame1, text='-', width=4, height=1,command=self.fileminus)
		self.buttonminus.grid(column=1, row=1, columnspan=1, sticky='W') 
		
		self.listbox3 = tkinter.Listbox(self.frame1)
		self.listbox3.grid(column=0, row=2, columnspan = 2)
		self.listbox3['width'] = 50
		self.listbox3['height'] = 53
		self.listbox3.bind('<<ListboxSelect>>', self.onselect)
		
		self.files_long_ft = []
		self.import_type = []
		
################################################################################
#FRAME2
################################################################################
		def FigureCanvas(self):
			self.fig = plt.figure.Figure(figsize=(10.75, 4.75))
			self.fig.subplots_adjust(left=0.05, bottom=0.08, right=0.98, top=0.93, wspace=0.12, hspace=0.02)
			
			self.canvasFrame = tkinter.Frame(master=self.frame2, padx=0, background="white")
			self.canvasFrame.grid(row=0,column=0)
			self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
			self.canvas.get_tk_widget().grid(column=0, row=0, sticky='EW')
			
			self.ax1=self.fig.add_subplot(131)
			self.ax2=self.fig.add_subplot(132)
			self.ax3=self.fig.add_subplot(133)
			self.canvas.draw_idle()
			
			self.toolbarFrame = tkinter.Frame(master=self.frame2)
			self.toolbarFrame.grid(row=1,column=0, sticky='W')
			self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
			self.toolbar.update()
			
			self.axis_color = 'lightgoldenrodyellow'	
			self.update()
		
		FigureCanvas(self)
		
################################################################################	
#FRAME3
################################################################################
		def raise_frame(frame):
			frame.tkraise()
			
		self.d1 = tkinter.Frame(master=self.frame3)
		self.c1 = tkinter.Frame(master=self.frame3)
		
		for frame in (self.d1, self.c1):
			frame.grid(row=0, column=0, sticky='news')
			
################################################################################			
		self.scrollFrame = tkinter.Frame(master=self.d1)
		self.scrollFrame.grid(row=0, column=0, rowspan=1, columnspan=2, pady=5)
		tkinter.Label(self.scrollFrame, width=15).grid(column=0, row=0, columnspan=1)
		tkinter.Label(self.scrollFrame, width=135).grid(column=1, row=0, columnspan=1)

		self.ds_scroll_scale = tkinter.Scale(self.scrollFrame, orient='horizontal', from_=0, to=1, state='disabled', width=18)
		self.ds_scroll_scale.grid(column=1, columnspan=1, row=1, stick='EW')
		self.datasetlabel = tkinter.Label(self.scrollFrame, text='Dataset',anchor='w')
		self.datasetlabel.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='S')
		
		raise_frame(self.d1)
		
################################################################################	
#FRAME4-1
################################################################################
		
		self.xasbackFrame = tkinter.Frame(master=self.frame4)
		self.xasbackFrame.grid(row=1, column=0, rowspan=1, columnspan=1, pady=5)
		tkinter.Label(self.xasbackFrame, width=10).grid(column=0, row=0, columnspan=1)
		tkinter.Label(self.xasbackFrame, width=39).grid(column=1, row=0, columnspan=1)
		
		self.kwbackchoices = []
		for i in range(int(4)):
			self.kwbackchoices.append(str(int(i)))
		
		self.kwback = tkinter.IntVar(self.xasbackFrame)
		self.kwbackMenu = tkinter.OptionMenu(self.xasbackFrame, self.kwback, *self.kwbackchoices)
		self.kwbackMenu.grid(column=1, row=0, columnspan=1) 
		self.kwbacklabel = tkinter.Label(self.xasbackFrame, text='k weight',anchor='w')		
		self.kwbacklabel.grid(column=0, row=0, columnspan=1, rowspan=1, sticky='S')

		self.lclamp_scroll_scale = tkinter.Scale(self.xasbackFrame, orient='horizontal', from_=0, to=50, state='normal', width=18)
		self.lclamp_scroll_scale.grid(column=1, columnspan=1, row=1, stick='EW')
		self.lclamplabel = tkinter.Label(self.xasbackFrame, text='Low Clamp',anchor='w')
		self.lclamplabel.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='S')
		
		self.hclamp_scroll_scale = tkinter.Scale(self.xasbackFrame, orient='horizontal', from_=0, to=50, state='normal', width=18)
		self.hclamp_scroll_scale.grid(column=1, columnspan=1, row=2, stick='EW')
		self.hclamplabel = tkinter.Label(self.xasbackFrame, text='High Clamp',anchor='w')
		self.hclamplabel.grid(column=0, row=2, columnspan=1, rowspan=1, sticky='S')
		
		self.rbkg_scroll_scale = tkinter.Scale(self.xasbackFrame, orient='horizontal', from_=0.5, to=1.5, digits = 3, resolution = 0.01, state='normal', width=18)
		self.rbkg_scroll_scale.grid(column=1, columnspan=1, row=3, stick='EW')
		self.rbkglabel = tkinter.Label(self.xasbackFrame, text='RBKG',anchor='w')
		self.rbkglabel.grid(column=0, row=3, columnspan=1, rowspan=1, sticky='S')
		
		self.e0_scroll_scale = tkinter.Scale(self.xasbackFrame, orient='horizontal', from_=0, to=1, state='normal', width=18)
		self.e0_scroll_scale.grid(column=1, columnspan=1, row=4, stick='EW')
		self.e0label = tkinter.Label(self.xasbackFrame, text='E0 (eV)',anchor='w')
		self.e0label.grid(column=0, row=4, columnspan=1, rowspan=1, sticky='S')
			
		self.lclamp_scroll_scale.set(10)
		self.hclamp_scroll_scale.set(30)
		self.rbkg_scroll_scale.set(1)
		self.kwback.set(1)
		
################################################################################	
#FRAME4-2
################################################################################
		
		self.FTFrame = tkinter.Frame(master=self.frame4)
		self.FTFrame.grid(row=1, column=1, rowspan=1, columnspan=1, pady=5)
		tkinter.Label(self.FTFrame, width=10).grid(column=0, row=0, columnspan=1)
		tkinter.Label(self.FTFrame, width=39).grid(column=1, row=0, columnspan=1)
		
		self.kwftchoices = []
		for i in range(int(4)):
			self.kwftchoices.append(str(int(i)))
		
		self.kwft = tkinter.IntVar(self.FTFrame)
		self.kwftMenu = tkinter.OptionMenu(self.FTFrame, self.kwft, *self.kwftchoices)
		self.kwftMenu.grid(column=1, row=0, columnspan=1)
		self.kwftlabel = tkinter.Label(self.FTFrame, text='k weight',anchor='w')		
		self.kwftlabel.grid(column=0, row=0, columnspan=1, rowspan=1, sticky='S')

		self.kmin_scroll_scale = tkinter.Scale(self.FTFrame, orient='horizontal', from_=0, to=15, state='normal', width=18)
		self.kmin_scroll_scale.grid(column=1, columnspan=1, row=1, stick='EW')
		self.kminlabel = tkinter.Label(self.FTFrame, text='k min (Å\u207B\u00B9)',anchor='w')
		self.kminlabel.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='S')
		
		self.kmax_scroll_scale = tkinter.Scale(self.FTFrame, orient='horizontal', from_=0, to=15, state='normal', width=18)
		self.kmax_scroll_scale.grid(column=1, columnspan=1, row=2, stick='EW')
		self.kmaxlabel = tkinter.Label(self.FTFrame, text='k max (Å\u207B\u00B9)',anchor='w')
		self.kmaxlabel.grid(column=0, row=2, columnspan=1, rowspan=1, sticky='S')
		
		self.dk_scroll_scale = tkinter.Scale(self.FTFrame, orient='horizontal', from_=0, to=4, state='normal', width=18)
		self.dk_scroll_scale.grid(column=1, columnspan=1, row=3, stick='EW')
		self.dklabel = tkinter.Label(self.FTFrame, text='dk (Å\u207B\u00B9)',anchor='w')
		self.dklabel.grid(column=0, row=3, columnspan=1, rowspan=1, sticky='S')
		
		self.label2 = tkinter.Label(self.FTFrame, text='',anchor='w').grid(column=0, row=4, columnspan=1, rowspan=1, sticky='S')
			
		self.kmin_scroll_scale.set(2.5)
		self.kmax_scroll_scale.set(10)
		self.dk_scroll_scale.set(1.5)
		self.kwft.set(2)
		
################################################################################	
#FRAME4-3
################################################################################
		
		self.BFTFrame = tkinter.Frame(master=self.frame4)
		self.BFTFrame.grid(row=1, column=2, rowspan=1, columnspan=1, pady=5)
		tkinter.Label(self.BFTFrame, width=10).grid(column=0, row=0, columnspan=1)
		tkinter.Label(self.BFTFrame, width=39).grid(column=1, row=0, columnspan=1)

		self.rmin_scroll_scale = tkinter.Scale(self.BFTFrame, orient='horizontal', from_=0, to=10, state='normal', width=18)
		self.rmin_scroll_scale.grid(column=1, columnspan=1, row=0, stick='EW')
		self.rminlabel = tkinter.Label(self.BFTFrame, text='R min (Å)',anchor='w')
		self.rminlabel.grid(column=0, row=0, columnspan=1, rowspan=1, sticky='S')
		
		self.rmax_scroll_scale = tkinter.Scale(self.BFTFrame, orient='horizontal', from_=0, to=10, state='normal', width=18)
		self.rmax_scroll_scale.grid(column=1, columnspan=1, row=1, stick='EW')
		self.rmaxlabel = tkinter.Label(self.BFTFrame, text='R max (Å)',anchor='w')
		self.rmaxlabel.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='S')
		
		self.dR_scroll_scale = tkinter.Scale(self.BFTFrame, orient='horizontal', from_=0, to=4, state='normal', width=18)
		self.dR_scroll_scale.grid(column=1, columnspan=1, row=2, stick='EW')
		self.dRlabel = tkinter.Label(self.BFTFrame, text='dR (Å)',anchor='w')
		self.dRlabel.grid(column=0, row=2, columnspan=1, rowspan=1, sticky='S')
		
		self.label2 = tkinter.Label(self.BFTFrame, text='',anchor='w').grid(column=0, row=3, columnspan=1, rowspan=1, sticky='S')
			
		self.rmin_scroll_scale.set(1)
		self.rmax_scroll_scale.set(3)
		self.dR_scroll_scale.set(0)
		
		self.buttonbatch= tkinter.Button(self.BFTFrame, text='Export', width=30, height=1, command=self.batch_trigger)
		self.buttonbatch.grid(column=1, row=4, columnspan=1, rowspan=1)	
		
		self.initialize_plot = 0
		
		self.update()
			
################################################################################		
	def normfile_check(self, folder):
		try:
			#self.normalisation_values = (pd.read_csv(folder+'/normalisation.dat', sep='\t', header=None)).values
			return True
		except:
			#print('Cannot Find Normalisation Parameter File')
			#self.NormUsevar.set(0)
			return False
		
################################################################################
	def fileminus(self):        
		i = self.listbox3.curselection()[0]
		self.listbox3.delete(i)
		del self.files_long_ft[i]
		self.update()
		
################################################################################        
	def fileplus(self):
		self.askloadfile()
		
################################################################################
	def askloadfile(self):
		file = askopenfile(filetypes = (("dat files","*.dat"),("all files", "*.*")))
		try:
			dirname, file_short = os.path.split(os.path.abspath(file.name))
			self.files_long_ft.append(file)
			print('Loading File')
			self.listbox3.insert(END, file_short)
		except:
			print('Aborted Load')
			
################################################################################					
	def clear_figure(self):
		self.fig.clear()
		
		self.fig = plt.figure.Figure(figsize=(10.75, 4.75))
		self.fig.subplots_adjust(left=0.05, bottom=0.08, right=0.98, top=0.93, wspace=0.15, hspace=0.02)
		
		self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
		self.canvas.get_tk_widget().grid(column=0, row=0, sticky='EW')
		
		self.ax1=self.fig.add_subplot(131)
		self.ax2=self.fig.add_subplot(132)
		self.ax3=self.fig.add_subplot(133)
		
		self.toolbar.destroy()
		self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
		self.toolbar.update()
			
#################################################################################			
	def onselect(self,evt):			
		if 'Fourier' in App.get_selected(self):
			data_file_long = self.files_long_ft[self.listbox3.curselection()[0]].name
			print(data_file_long)
			self.folder, self.data_file = os.path.split(os.path.abspath(data_file_long))
			self.file_type=(os.path.splitext(self.data_file)[-1])
			print('File Type :', self.file_type)
			
			self.normalized_check = self.normfile_check(self.folder)
				
			if ('.dat' in self.file_type) or ('.txt' in self.file_types):
				print('Loading file ... ')
				self.data = pd.read_csv(data_file_long, sep='\t', header=0)
				try:
					self.data_error = pd.read_csv(data_file_long.split('.dat')[0]+'.error', sep='\t', header=0)
					self.errors_flag = True
				except:
					self.errors_flag = False
				self.data.dropna(how='any', inplace = True)
				print('File loaded')
				self.column_names = self.data.columns.values.tolist()
				if self.errors_flag == True:
					self.error_column_names = self.data_error.columns.values.tolist()
				self.trigger = 0
				self.Kstep = 0.05
				
				try:
					parameters_file = open(self.folder+'/parameters.txt', 'r')
					
					self.normalised = 0
					
					for line in parameters_file:
						if 'Processed with Normalisation' in line:
							self.normalised = 1
						#if 'Kstep' in line:
						#	self.Kstep = float(ast.literal_eval(line.split(' ')[-1]))
				except:
					print('No parameters file detected')
					self.normalised = 1
				
				self.update()
				
				if any('normalised' in s for s in self.data_file.split('_')):
					self.normalised = 1
					
				self.trigger = 0
				self.plot_spectrum()
			
#################################################################################			
	def plot_spectrum(self):				
		
		if self.trigger !=1:
			self.initialize_ft = 0
			self.initialize_plot = 0
	
		mylarch = Interpreter()
		global larch_data
		larch_data =  larch.Group(name='larch_data')
		
		def plot_data(evt):
			if self.trigger == 0:
				self.data_i = self.ds_scroll_scale.get()
				
			mass_e = 0.5109989461 * 10**6
			hbarc = 0.19732697
				
			self.data_x = self.data[str(self.column_names[0])]
			self.data_y =  self.data[str(self.column_names[self.data_i+1])]
			
			def adjustErrbarxy(self, errobj, x, y, y_error):
				ln, (erry_top, erry_bot), (bars,) = errobj
				x_base = x
				y_base = y
			
				yerr_top = y_base + y_error
				yerr_bot = y_base - y_error
			
				erry_top.set_xdata(x_base)
				erry_bot.set_xdata(x_base)
				erry_top.set_ydata(yerr_top)
				erry_bot.set_ydata(yerr_bot)
				
				new_segments_y = [np.array([[x, yt], [x,yb]]) for x, yt, yb in zip(x_base, yerr_top, yerr_bot)]
				bars.set_segments(new_segments_y)
				ln.set_data(x,y)
				return errobj
			
			if self.initialize_plot == 0:
				
				self.clear_figure()
				[self.lineE] = self.ax1.plot(self.data_x, self.data_y)
				
				larch_data.energy = self.data_x
				larch_data.mu = self.data_y
				
				autobk(larch_data.energy, larch_data.mu, e0=self.e0_scroll_scale.get(), rbkg=self.rbkg_scroll_scale.get(), group=larch_data, clamp_lo=self.lclamp_scroll_scale.get(), clamp_hi=self.hclamp_scroll_scale.get(), kweight=self.kwback.get(), dk=self.dk_scroll_scale.get(), kstep=self.Kstep, _larch=mylarch)
				[self.lineBack] = self.ax1.plot(larch_data.energy, larch_data.bkg, linewidth=2, color='red')
				
				xftf(larch_data.k, larch_data.chi, kmin=self.kmin_scroll_scale.get(), kmax=self.kmax_scroll_scale.get(), dk=self.dk_scroll_scale.get(), window='hanning', kweight=self.kwft.get(), group=larch_data, _larch=mylarch)
				xftr(larch_data.r, larch_data.chir, group=larch_data, qmax_out=self.kmax_scroll_scale.get(), rmin=self.rmin_scroll_scale.get(), rmax=self.rmax_scroll_scale.get(), window='hanning', dr=self.dR_scroll_scale.get(), _larch=mylarch)
		
				self.kwin = ftwindow(larch_data.k, xmin=self.kmin_scroll_scale.get(), xmax=self.kmax_scroll_scale.get(), dx=self.dk_scroll_scale.get(), dx2=self.dk_scroll_scale.get(), window='hanning')
				
				self.EKerror_data['E'] = larch_data.energy - self.e0_scroll_scale.get()
				self.EKerror_data['K'] = np.sqrt((2*mass_e*self.EKerror_data['E'].values)/((hbarc**2)*(10**8)))
				
				self.EKerror_data['kw'] = larch_data.chie * (self.EKerror_data['K']**self.kwft.get())
				if self.errors_flag == True:
					self.EKerror_data['kw_error']  = self.data_error[str(self.error_column_names[self.data_i+1])] * (self.EKerror_data['K']**self.kwft.get())
				
				self.EKerror_data = self.EKerror_data[(self.EKerror_data['E'] >= 0)].reset_index(drop=True)
				if self.errors_flag == True:
					self.linek_errorbar = self.ax2.errorbar(self.EKerror_data['K'], self.EKerror_data['kw'], yerr=self.EKerror_data['kw_error'], ecolor='r', capsize=3, elinewidth=2, markeredgewidth=1)
				else:
					[self.linek] = self.ax2.plot(self.EKerror_data['K'], self.EKerror_data['kw'], linewidth=2, color='red')
				[self.line_kmin] = self.ax2.plot([self.kmin_scroll_scale.get(), self.kmin_scroll_scale.get()], [-1,1], linewidth = 2, color='k', ls='--')
				[self.line_kmax] = self.ax2.plot([self.kmax_scroll_scale.get(), self.kmax_scroll_scale.get()], [-1,1], linewidth = 2, color='k', ls='--')
				[self.line_z] = self.ax2.plot([0, np.max(larch_data.k)], [0,0], linewidth = 1, color='black')
				[self.line_kwin] = self.ax2.plot(larch_data.k, self.kwin, linewidth = 2, color = 'green')
				
				self.rwin = ftwindow(larch_data.r, xmin=self.rmin_scroll_scale.get(), xmax=self.rmax_scroll_scale.get(), dx=self.dR_scroll_scale.get(), dx2=self.dR_scroll_scale.get(), window='hanning')
				
				[self.liner] = self.ax3.plot(larch_data.r, larch_data.chir_mag, linewidth=2, color='red')
				[self.line_rwin] = self.ax3.plot(larch_data.r, self.rwin, linewidth = 2, color = 'green')
				
				self.kwback_store = self.kwback.get()
				self.kft_store = self.kwft.get()
				self.kmax_store = self.kmax_scroll_scale.get()
				self.kmin_store = self.kmin_scroll_scale.get()
				
				self.canvas.draw_idle()
				self.initialize_plot = 1
			else:
				self.EKerror_data = pd.DataFrame()
				larch_data.energy = self.data_x
				larch_data.mu = self.data_y
				
				autobk(larch_data.energy, larch_data.mu, e0=self.e0_scroll_scale.get(), rbkg=self.rbkg_scroll_scale.get(), group=larch_data, clamp_lo=self.lclamp_scroll_scale.get(), clamp_hi=self.hclamp_scroll_scale.get(), kweight=self.kwback.get(),kstep=self.Kstep, _larch=mylarch)
				
				xftf(larch_data.k, larch_data.chi, kmin=self.kmin_scroll_scale.get(), kmax=self.kmax_scroll_scale.get(), dk=self.dk_scroll_scale.get(), window='hanning', kweight=self.kwft.get(), group=larch_data, _larch=mylarch)
				xftr(larch_data.r, larch_data.chir, group=larch_data, qmax_out=self.kmax_scroll_scale.get(), rmin=self.rmin_scroll_scale.get(), rmax=self.rmax_scroll_scale.get(), window='hanning', dr=self.dR_scroll_scale.get(), _larch=mylarch)
		
				self.kwin = ftwindow(larch_data.k, xmin=self.kmin_scroll_scale.get(), xmax=self.kmax_scroll_scale.get(), dx=self.dk_scroll_scale.get(), dx2=self.dk_scroll_scale.get(), window='hanning')
				self.rwin = ftwindow(larch_data.r, xmin=self.rmin_scroll_scale.get(), xmax=self.rmax_scroll_scale.get(), dx=self.dR_scroll_scale.get(), dx2=self.dR_scroll_scale.get(), window='hanning')
				
				self.lineE.set_ydata(self.data_y)
				self.lineBack.set_ydata(larch_data.bkg)
				
				self.EKerror_data['E'] = larch_data.energy - self.e0_scroll_scale.get()
				self.EKerror_data['K'] = np.sqrt((2*mass_e*self.EKerror_data['E'].values)/((hbarc**2)*(10**8)))
				
				self.EKerror_data['kw'] = larch_data.chie * (self.EKerror_data['K']**self.kwft.get())
				if self.errors_flag == True:
					self.EKerror_data['kw_error']  = self.data_error[str(self.error_column_names[self.data_i+1])] * (self.EKerror_data['K']**self.kwft.get())
					
				self.EKerror_data = self.EKerror_data[(self.EKerror_data['E'] >= 0)].reset_index(drop=True)
				
				if self.errors_flag == True:
					self.linek_errorbar = adjustErrbarxy(self, self.linek_errorbar, self.EKerror_data['K'], self.EKerror_data['kw'], self.EKerror_data['kw_error'])
				else:
					self.linek.set_data(self.EKerror_data['K'], self.EKerror_data['kw'])
				self.line_kmin.set_xdata([self.kmin_scroll_scale.get(), self.kmin_scroll_scale.get()])
				self.line_kmax.set_xdata([self.kmax_scroll_scale.get(), self.kmax_scroll_scale.get()])
				self.line_z.set_xdata([0, np.max(larch_data.k)])
				self.line_kwin.set_data(larch_data.k, self.kwin)
				
				self.liner.set_data(larch_data.r, larch_data.chir_mag)
				self.line_rwin.set_data(larch_data.r, self.rwin)
				
				if self.kwback_store != self.kwback.get():
					self.kwback_store = self.kwback.get()
					self.ax2.relim(visible_only=False)
					self.ax2.autoscale_view(tight=None, scalex=True, scaley=True)
					self.ax3.relim(visible_only=False)
					self.ax3.autoscale_view(tight=None, scalex=True, scaley=True)
				if self.kft_store != self.kwft.get():
					self.kft_store = self.kwft.get()
					self.ax2.relim(visible_only=False)
					self.ax2.autoscale_view(tight=None, scalex=True, scaley=True)
					self.ax3.relim(visible_only=False)
					self.ax3.autoscale_view(tight=None, scalex=True, scaley=True)
				if self.kmax_store != self.kmax_scroll_scale.get():
					self.kmax_store = self.kmax_scroll_scale.get()
					self.ax3.relim(visible_only=False)
					self.ax3.autoscale_view(tight=None, scalex=True, scaley=True)
				if self.kmin_store != self.kmin_scroll_scale.get():
					self.kmin_store = self.kmin_scroll_scale.get()
					self.ax3.relim(visible_only=False)
					self.ax3.autoscale_view(tight=None, scalex=True, scaley=True)
				
				self.canvas.draw_idle()
				
				if self.trigger != 0:
					if self.data_i == 0:
						self.dfkw['k'] = larch_data.k
						self.dfk['k'] = larch_data.k
						self.dfrr['r'] = larch_data.r
						self.dfri['r'] = larch_data.r
						self.dfrm['r'] = larch_data.r
						self.dfq['q'] = larch_data.q
						self.chiE['E'] = larch_data.energy - self.e0_scroll_scale.get()
						self.chiKw['k'] = np.sqrt((2*mass_e*self.chiE['E'].values)/((hbarc**2)*(10**8)))
						if self.errors_flag == True:
							self.chiKwerror['k'] = np.sqrt((2*mass_e*self.chiE['E'].values)/((hbarc**2)*(10**8)))
					
					self.dfkw[str(self.data_i)] = larch_data.chi * (self.dfk['k']**self.kwft.get())
					self.dfk[str(self.data_i)] = larch_data.chi
					self.dfrr[str(self.data_i)] = larch_data.chir_re
					self.dfri[str(self.data_i)] = larch_data.chir_im
					self.dfrm[str(self.data_i)] = larch_data.chir_mag
					self.dfq[str(self.data_i)] = larch_data.chiq_re
					self.chiE[str(self.data_i)] = larch_data.chie
					self.chiKw[str(self.data_i)] = larch_data.chie * (self.chiKw['k']**self.kwft.get())
					if self.errors_flag == True:
						self.chiKwerror[str(self.data_i)] = self.data_error[str(self.error_column_names[self.data_i+1])] * (self.chiKw['k']**self.kwft.get())
					
				return 'updated'
		
		if self.trigger == 0:	
		
			self.EKerror_data = pd.DataFrame()
			
			self.ds_scroll_scale.set(0)
			if len(self.column_names) == 2:
				self.ds_scroll_scale = tkinter.Scale(self.scrollFrame, orient='horizontal', from_=0, to=1, state='disabled', width=18, command=plot_data)
			else:
				self.ds_scroll_scale = tkinter.Scale(self.scrollFrame, orient='horizontal', from_=0, to=len(self.column_names)-2, state='normal', width=18, command=plot_data)
			self.ds_scroll_scale.grid(column=1, columnspan=2, row=1, stick='EW')
			self.datasetlabel = tkinter.Label(self.scrollFrame, text='Dataset',anchor='w')
			self.datasetlabel.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='S')
		
			if self.initialize_ft == 0:
				self.data_i = self.ds_scroll_scale.get()
				self.data_x = self.data[str(self.column_names[0])]
				self.data_y =  self.data[str(self.column_names[self.data_i+1])]
				larch_data.energy = self.data_x
				larch_data.mu = self.data_y
				autobk(larch_data.energy, larch_data.mu, rbkg=self.rbkg_scroll_scale.get(), group=larch_data, clamp_lo=self.lclamp_scroll_scale.get(), clamp_hi=self.hclamp_scroll_scale.get(), kweight=self.kwback.get(), kstep=self.Kstep, _larch=mylarch)
				xftf(larch_data.k, larch_data.chi, kmin=self.kmin_scroll_scale.get(), kmax=self.kmax_scroll_scale.get(), dk=self.dk_scroll_scale.get(), window='hanning', kweight=self.kwft.get(), group=larch_data, _larch=mylarch)
				xftr(larch_data.r, larch_data.chir, group=larch_data, qmax_out=self.kmax_scroll_scale.get(), rmin=self.rmin_scroll_scale.get(), rmax=self.rmax_scroll_scale.get(), window='hanning', dr=self.dR_scroll_scale.get(), _larch=mylarch)
				self.initialize_ft = 1
			
			if self.initialize_plot == 0:	
				self.kwback = tkinter.IntVar(self.xasbackFrame)
				self.kwbackMenu = tkinter.OptionMenu(self.xasbackFrame, self.kwback, *self.kwbackchoices, command=plot_data)
				self.kwbackMenu.grid(column=1, row=0, columnspan=1) 
				self.kwbacklabel = tkinter.Label(self.xasbackFrame, text='k weight',anchor='w')		
				self.kwbacklabel.grid(column=0, row=0, columnspan=1, rowspan=1, sticky='S')
		
				self.lclamp_scroll_scale = tkinter.Scale(self.xasbackFrame, orient='horizontal', from_=0, to=50, state='normal', width=18, command=plot_data)
				self.lclamp_scroll_scale.grid(column=1, columnspan=1, row=1, stick='EW')
				self.lclamplabel = tkinter.Label(self.xasbackFrame, text='Low Clamp',anchor='w')
				self.lclamplabel.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='S')
				
				self.hclamp_scroll_scale = tkinter.Scale(self.xasbackFrame, orient='horizontal', from_=0, to=50, state='normal', width=18, command=plot_data)
				self.hclamp_scroll_scale.grid(column=1, columnspan=1, row=2, stick='EW')
				self.hclamplabel = tkinter.Label(self.xasbackFrame, text='High Clamp',anchor='w')
				self.hclamplabel.grid(column=0, row=2, columnspan=1, rowspan=1, sticky='S')
				
				self.rbkg_scroll_scale = tkinter.Scale(self.xasbackFrame, orient='horizontal', from_=0.5, to=1.5, digits = 3, resolution = 0.01, state='normal', width=18, command=plot_data)
				self.rbkg_scroll_scale.grid(column=1, columnspan=1, row=3, stick='EW')
				self.rbkglabel = tkinter.Label(self.xasbackFrame, text='RBKG',anchor='w')
				self.rbkglabel.grid(column=0, row=3, columnspan=1, rowspan=1, sticky='S')
				
				self.e0_scroll_scale = tkinter.Scale(self.xasbackFrame, orient='horizontal', from_=larch_data.e0-20, to=larch_data.e0+20, digits = 6, resolution = 0.01, state='normal', width=18, command=plot_data)
				self.e0_scroll_scale.grid(column=1, columnspan=1, row=4, stick='EW')
				self.e0label = tkinter.Label(self.xasbackFrame, text='E0 (eV)',anchor='w')
				self.e0label.grid(column=0, row=4, columnspan=1, rowspan=1, sticky='S')
				
				self.lclamp_scroll_scale.set(10)
				self.hclamp_scroll_scale.set(30)
				self.rbkg_scroll_scale.set(1)
				self.kwback.set(1)
				self.e0_scroll_scale.set(larch_data.e0)
				
				self.kwft = tkinter.IntVar(self.FTFrame)
				self.kwftMenu = tkinter.OptionMenu(self.FTFrame, self.kwft, *self.kwftchoices, command=plot_data)
				self.kwftMenu.grid(column=1, row=0, columnspan=1)
				self.kwftlabel = tkinter.Label(self.FTFrame, text='k weight',anchor='w')		
				self.kwftlabel.grid(column=0, row=0, columnspan=1, rowspan=1, sticky='S')
		
				self.kmin_scroll_scale = tkinter.Scale(self.FTFrame, orient='horizontal', from_=0, to=np.max(larch_data.k), digits = 4, resolution = 0.01, state='normal', width=18, command=plot_data)
				self.kmin_scroll_scale.grid(column=1, columnspan=1, row=1, stick='EW')
				self.kminlabel = tkinter.Label(self.FTFrame, text='k min (Å\u207B\u00B9)',anchor='w')
				self.kminlabel.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='S')
				
				self.kmax_scroll_scale = tkinter.Scale(self.FTFrame, orient='horizontal', from_=0, to=np.max(larch_data.k), digits = 4, resolution = 0.01, state='normal', width=18, command=plot_data)
				self.kmax_scroll_scale.grid(column=1, columnspan=1, row=2, stick='EW')
				self.kmaxlabel = tkinter.Label(self.FTFrame, text='k max (Å\u207B\u00B9)',anchor='w')
				self.kmaxlabel.grid(column=0, row=2, columnspan=1, rowspan=1, sticky='S')
				
				self.dk_scroll_scale = tkinter.Scale(self.FTFrame, orient='horizontal', from_=0, to=4, state='normal', digits = 3, resolution = 0.01, width=18, command=plot_data)
				self.dk_scroll_scale.grid(column=1, columnspan=1, row=3, stick='EW')
				self.dklabel = tkinter.Label(self.FTFrame, text='dk (Å\u207B\u00B9)',anchor='w')
				self.dklabel.grid(column=0, row=3, columnspan=1, rowspan=1, sticky='S')
				
				self.kmin_scroll_scale.set(2.5)
				self.kmax_scroll_scale.set(10)
				self.dk_scroll_scale.set(1.5)
				self.kwft.set(2)
				
				self.rmin_scroll_scale = tkinter.Scale(self.BFTFrame, orient='horizontal', from_=0, to=np.max(larch_data.r), digits = 3, resolution = 0.01, state='normal', width=18, command=plot_data)
				self.rmin_scroll_scale.grid(column=1, columnspan=1, row=0, stick='EW')
				self.rminlabel = tkinter.Label(self.BFTFrame, text='R min (Å)',anchor='w')
				self.rminlabel.grid(column=0, row=0, columnspan=1, rowspan=1, sticky='S')
				
				self.rmax_scroll_scale = tkinter.Scale(self.BFTFrame, orient='horizontal', from_=0, to=np.max(larch_data.r), digits = 3, resolution = 0.01, state='normal', width=18, command=plot_data)
				self.rmax_scroll_scale.grid(column=1, columnspan=1, row=1, stick='EW')
				self.rmaxlabel = tkinter.Label(self.BFTFrame, text='R max (Å)',anchor='w')
				self.rmaxlabel.grid(column=0, row=1, columnspan=1, rowspan=1, sticky='S')
				
				self.dR_scroll_scale = tkinter.Scale(self.BFTFrame, orient='horizontal', from_=0, to=4, state='normal', digits = 3, resolution = 0.01, width=18, command=plot_data)
				self.dR_scroll_scale.grid(column=1, columnspan=1, row=2, stick='EW')
				self.dRlabel = tkinter.Label(self.BFTFrame, text='dR (Å)',anchor='w')
				self.dRlabel.grid(column=0, row=2, columnspan=1, rowspan=1, sticky='S')
	
				self.rmin_scroll_scale.set(1)
				self.rmax_scroll_scale.set(3)
				self.dR_scroll_scale.set(0)
			
			plot_data(None)
			
		if self.trigger == 1:
			self.dfrr = pd.DataFrame()
			self.dfri = pd.DataFrame()
			self.dfrm = pd.DataFrame()
			self.dfk = pd.DataFrame()
			self.dfkw = pd.DataFrame()
			self.dfq = pd.DataFrame()
			self.chiE = pd.DataFrame()
			self.chiKw = pd.DataFrame()
			if self.errors_flag == True:
				self.chiKwerror = pd.DataFrame()
			
			if len(self.column_names) == 2:
				ds_range = 0
				self.data_i = 0
				self.trigger = 1
				plot_data(None)
			else:
				ds_range = len(self.column_names)-1
			
				for i in range(ds_range):
					self.data_i = i
					self.ds_scroll_scale.set(i)
					plot_data(None)
					self.update()

			if not os.path.exists(self.folder+'/FT/'+str(os.path.splitext(self.data_file)[0])):
				os.makedirs(self.folder+'/FT/'+str(os.path.splitext(self.data_file)[0]))
				
			#self.rmcr = pd.concat([self.dfrr, self.dfri], axis=1)
			#self.rmcr.to_csv(self.folder+'/FT/'+str(os.path.splitext(self.data_file)[0])+'/'+str(os.path.splitext(self.data_file)[0])+'_R-space_ri_mcr.dat', sep='\t', index=False)
			self.dfrm.to_csv(self.folder+'/FT/'+str(os.path.splitext(self.data_file)[0])+'/'+str(os.path.splitext(self.data_file)[0])+'_R-space_magnitude.dat', sep='\t', index=False)
			self.dfkw.to_csv(self.folder+'/FT/'+str(os.path.splitext(self.data_file)[0])+'/'+str(os.path.splitext(self.data_file)[0])+'_K'+str(self.kwft.get())+'-space.dat', sep='\t', index=False)
			self.dfk.to_csv(self.folder+'/FT/'+str(os.path.splitext(self.data_file)[0])+'/'+str(os.path.splitext(self.data_file)[0])+'_K-space.dat', sep='\t', index=False)
			#self.dfq.to_csv(self.folder+'/FT/'+str(os.path.splitext(self.data_file)[0])+'/'+str(os.path.splitext(self.data_file)[0])+'_Q-space.dat', sep='\t', index=False)
			
			if self.errors_flag == True:
				self.chiKwerror = self.chiKwerror[(self.chiE['E'] >= 0)].reset_index(drop=True)
			self.chiKw = self.chiKw[(self.chiE['E'] >= 0)].reset_index(drop=True)
			self.chiE = self.chiE[(self.chiE['E'] >= 0)].reset_index(drop=True)
			self.chiE.to_csv(self.folder+'/FT/'+str(os.path.splitext(self.data_file)[0])+'/'+str(os.path.splitext(self.data_file)[0])+'_chiE.dat', sep='\t', index=False)	
			self.chiKw.to_csv(self.folder+'/FT/'+str(os.path.splitext(self.data_file)[0])+'/'+str(os.path.splitext(self.data_file)[0])+'_chikw.dat', sep='\t', index=False)	
			if self.errors_flag == True:
				self.chiKwerror.to_csv(self.folder+'/FT/'+str(os.path.splitext(self.data_file)[0])+'/'+str(os.path.splitext(self.data_file)[0])+'_chikw.error', sep='\t', index=False)	
				
			self.trigger = 0
		
#################################################################################			
	def batch_trigger(self):	
		self.trigger = 1
		self.plot_spectrum()
		
#################################################################################			
##PlottingTAB			
#################################################################################	
class Plotting(tkinter.Frame):
	def __init__(self,name,*args,**kwargs):
		self.frame_dataread = tkinter.Frame.__init__(self,*args,**kwargs)
				
		self.screen_width = self.winfo_screenwidth()
		self.screen_height = self.winfo_screenheight()
		
		self.progress_bar = ttk.Progressbar(self, orient='horizontal', mode='determinate', maximum=100, style="red.Horizontal.TProgressbar")
		self.progress_bar.grid(column=0, row=99, columnspan=99)
		self.progress_bar['length'] = 1000
		
		self.labelText = tkinter.StringVar()
		self.labelText.set('Time estimate')
		self.labeltime = tkinter.Label(self, textvariable=self.labelText)
		self.labeltime.grid(column=1, row=98, columnspan=2, rowspan=1, sticky='W')
		
		self.tpstextvar = tkinter.StringVar()
		self.tpstextvar.set('Time Per Spectrum')
		self.labeltps = tkinter.Label(self, textvariable=self.tpstextvar)
		self.labeltps.grid(column=3, row=98, columnspan=1, rowspan=1, sticky='W')
		
		self.frame1 = tkinter.LabelFrame(self, width=0.25*0.8333*self.screen_height, height=0.75*self.screen_height)
		self.frame2 = tkinter.LabelFrame(self, width=0.5*0.8333*self.screen_height, height=0.7*0.75*self.screen_height)
		self.frame3 = tkinter.LabelFrame(self, width=0.25*0.8333*self.screen_height, height=0.75*self.screen_height, borderwidth = 0, highlightthickness = 0)
		self.frame1.grid(row=0, column=0, rowspan=4, columnspan=1, padx=2, sticky='N')
		self.frame2.grid(row=0, column=1, rowspan=2, columnspan=2, padx=2, pady=15)
		self.frame3.grid(row=0, column=3, rowspan=4, columnspan=1, padx=2, pady=15, sticky='N')
		
################################################################################
#FRAME1
################################################################################		
		labelfile = tkinter.Label(self.frame1, text='Load Data',anchor='w')
		labelfile.grid(column=0, row=0, columnspan=3, rowspan=1, sticky='S')
		
		self.buttonplus= tkinter.Button(self.frame1, text='+', width=4, height=1,command=self.fileplus)
		self.buttonplus.grid(column=0, row=1, columnspan=1, sticky='E')
		
		self.buttonminus= tkinter.Button(self.frame1, text='-', width=4, height=1,command=self.fileminus)
		self.buttonminus.grid(column=1, row=1, columnspan=1, sticky='W') 
		
		self.listbox4 = tkinter.Listbox(self.frame1, selectmode=MULTIPLE)
		self.listbox4.grid(column=0, row=3, columnspan = 3)
		self.listbox4['width'] = 50
		self.listbox4['height'] = 53
		
		self.listboxbtn = ttk.Button(self.frame1, text="Plot Files", command=self.onplot)
		self.listboxbtn.grid(column=2, row=1, columnspan = 1)
		
		self.files_long_plot = []
		
################################################################################
#FRAME2
################################################################################
		def FigureCanvas(self):
			self.fig = plt.figure.Figure(figsize=(5.75, 5.75))
			
			self.canvasFrame = tkinter.Frame(master=self.frame2, padx=16, background="white")
			self.canvasFrame.grid(row=0,column=0)
			self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
			self.canvas.get_tk_widget().grid(column=0, row=0, sticky='EW')
			
			self.ax=self.fig.add_subplot(111)
			self.canvas.draw_idle()
			
			self.toolbarFrame = tkinter.Frame(master=self.frame2)
			self.toolbarFrame.grid(row=1,column=0, sticky='W')
			self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
			self.toolbar.update()
			
			self.axis_color = 'lightgoldenrodyellow'	
			self.update()
		
		FigureCanvas(self)

################################################################################
#FRAME3
################################################################################		
		tkinter.Label(self.frame3, width=50).grid(column=0, row=0, columnspan=3)
		
################################################################################
	def fileminus(self):        
		i = self.listbox4.curselection()[0]
		self.listbox4.delete(i)
		del self.files_long_plot[i]
		self.update()
		
################################################################################        
	def fileplus(self):
		self.askloadfile()
		
################################################################################
	def askloadfile(self):
		file = askopenfile(filetypes = (("dat files","*.dat"),("all files", "*.*")))
		try:
			dirname, file_short = os.path.split(os.path.abspath(file.name))
			self.files_long_plot.append(file)
			print('Loading File')
			self.listbox4.insert(END, file_short)
		except:
			print('Aborted Load')
			
################################################################################					
	def clear_figure(self):
		self.fig.clear()
		self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
		self.canvas.get_tk_widget().grid(column=0, row=0, sticky='EW')
		
		self.ax=self.fig.add_subplot(111)
		self.toolbar.destroy()
		self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
		self.toolbar.update()
		
################################################################################	
	def onplot(self):
		reslist = list()
		selection = self.listbox4.curselection()
		self.dataframe_collection = {}
		self.line_collection = {}
		self.sliders = {}

		self.clear_figure()
		
		for widget in self.frame3.winfo_children():
			if str(type(widget)) ==  "<class 'tkinter.Scale'>":
				widget.destroy()
			
		self.update()
		
		for i in selection:
			dirname, file_short = os.path.split(os.path.abspath(self.files_long_plot[i].name))
			reslist.append(file_short)
			self.dataframe_collection[str(i)] = pd.read_csv(self.files_long_plot[i].name, sep='\t')
			print(file_short)
			try:
				data_x = self.dataframe_collection[str(i)][str(self.dataframe_collection[str(i)].columns['E'])]
			except:
				try:
					data_x = self.dataframe_collection[str(i)][str(self.dataframe_collection[str(i)].columns['Energy'])]
				except:
					data_x = self.dataframe_collection[str(i)].iloc[:,0]
			
			data_y = self.dataframe_collection[str(i)].iloc[:,1]
			
			self.line_collection[str(i)] = self.ax.plot(data_x, data_y, label=file_short)
			
			self.sliders[str(i)] = tkinter.Scale(self.frame3, orient='horizontal', from_=0, to=len(self.dataframe_collection[str(i)].columns)-2, state='normal', width=18, sliderlength=10)
			self.sliders[str(i)].grid(column=0, columnspan=3, row=i, stick='EW')
			
		self.ax.legend(loc='lower right')
		self.canvas.draw_idle()
			
		self.store_state = []
		self.sliders_dict_keys = list(self.sliders.keys())

		for j in range(len(self.sliders)):
			self.store_state.append(0)
			self.sliders[self.sliders_dict_keys[j]].config(command = lambda iwid = j: self.plot_update_trigger(j))
			
################################################################################		
	def plot_update_trigger(self, *args):
			
		for j in range (len(self.sliders)):
			if self.sliders[self.sliders_dict_keys[j]].get() != self.store_state[j]:
				data_x = self.dataframe_collection[self.sliders_dict_keys[j]].iloc[:,0]
				data_y = self.dataframe_collection[self.sliders_dict_keys[j]].iloc[:,int(self.sliders[self.sliders_dict_keys[j]].get())+1]
				self.line_collection[self.sliders_dict_keys[j]][0].set_data(data_x, data_y)
				self.store_state[j] = self.sliders[self.sliders_dict_keys[j]].get()
				
				self.canvas.draw_idle()

#################################################################################			
##PlottingTAB			
#################################################################################	
class Analysis(tkinter.Frame):
	def __init__(self,name,*args,**kwargs):
		self.frame_dataread = tkinter.Frame.__init__(self,*args,**kwargs)
				
		self.screen_width = self.winfo_screenwidth()
		self.screen_height = self.winfo_screenheight()
		
		self.progress_bar = ttk.Progressbar(self, orient='horizontal', mode='determinate', maximum=100, style="red.Horizontal.TProgressbar")
		self.progress_bar.grid(column=0, row=99, columnspan=99)
		self.progress_bar['length'] = 1000
		
		self.labelText = tkinter.StringVar()
		self.labelText.set('Time estimate')
		self.labeltime = tkinter.Label(self, textvariable=self.labelText)
		self.labeltime.grid(column=1, row=98, columnspan=2, rowspan=1, sticky='W')
		
		self.tpstextvar = tkinter.StringVar()
		self.tpstextvar.set('Time Per Spectrum')
		self.labeltps = tkinter.Label(self, textvariable=self.tpstextvar)
		self.labeltps.grid(column=3, row=98, columnspan=1, rowspan=1, sticky='W')
		
		self.frame1 = tkinter.LabelFrame(self, width=0.25*0.8333*self.screen_height, height=0.75*self.screen_height)
		self.frame2 = tkinter.LabelFrame(self, width=0.5*0.8333*self.screen_height, height=0.7*0.75*self.screen_height)
		self.frame3 = tkinter.LabelFrame(self, width=0.25*0.8333*self.screen_height, height=0.75*self.screen_height, borderwidth = 0, highlightthickness = 0)
		self.frame1.grid(row=0, column=0, rowspan=4, columnspan=1, padx=2, sticky='N')
		self.frame2.grid(row=0, column=1, rowspan=2, columnspan=2, padx=2, pady=15)
		self.frame3.grid(row=0, column=3, rowspan=4, columnspan=1, padx=2, pady=15, sticky='N')
		
################################################################################
#FRAME1
################################################################################		
		labelfile = tkinter.Label(self.frame1, text='Load Data',anchor='w')
		labelfile.grid(column=0, row=0, columnspan=3, rowspan=1, sticky='S')
		
		self.buttonplus= tkinter.Button(self.frame1, text='+', width=4, height=1,command=self.fileplus)
		self.buttonplus.grid(column=0, row=1, columnspan=1, sticky='E')
		
		self.buttonminus= tkinter.Button(self.frame1, text='-', width=4, height=1,command=self.fileminus)
		self.buttonminus.grid(column=1, row=1, columnspan=1, sticky='W') 
		
		self.listbox5 = tkinter.Listbox(self.frame1)
		self.listbox5.grid(column=0, row=3, columnspan = 3)
		self.listbox5['width'] = 50
		self.listbox5['height'] = 53
		self.listbox5.bind('<<ListboxSelect>>', self.onselect)
		
		self.files_long_analysis = []
		
################################################################################
#FRAME2
################################################################################
		def FigureCanvas(self):
			self.fig = plt.figure.Figure(figsize=(5.75, 5.75))
			
			self.canvasFrame = tkinter.Frame(master=self.frame2, padx=16, background="white")
			self.canvasFrame.grid(row=0,column=0)
			self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
			self.canvas.get_tk_widget().grid(column=0, row=0, sticky='EW')
			
			self.ax=self.fig.add_subplot(111)
			self.canvas.draw_idle()
			
			self.toolbarFrame = tkinter.Frame(master=self.frame2)
			self.toolbarFrame.grid(row=1,column=0, sticky='W')
			self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
			self.toolbar.update()
			
			self.axis_color = 'lightgoldenrodyellow'	
			self.update()
		
		FigureCanvas(self)

################################################################################
#FRAME3
################################################################################		
		self.AnalzeFrame = tkinter.Frame(master=self.frame3)
		self.AnalzeFrame.grid(row=0, column=0, rowspan=1, columnspan=2, pady=5)
		
		tkinter.Label(self.AnalzeFrame, width=50).grid(column=0, row=0, columnspan=3)
		
		self.buttonsvd= tkinter.Button(self.AnalzeFrame, text='SVD', width=60, height=1,command=self.svd)
		self.buttonsvd.grid(column=0, row=1, columnspan=3, rowspan=1)
		
		svd_numer_label = tkinter.Label(self.AnalzeFrame, text='Nmax SVD',anchor='w')
		svd_numer_label.grid(column=0, row=2, columnspan=1, rowspan=1, sticky='W')
		self.svd_number = tkinter.IntVar()
		self.svd_number_entry = tkinter.Entry(self.AnalzeFrame, textvariable=self.svd_number, state='normal')
		self.svd_number_entry.grid(column=1, columnspan=2, row=2, stick='EW')
		
		self.buttonsimplisma= tkinter.Button(self.AnalzeFrame, text='SIMPLISMA', width=40, height=1,command=self.simplisma)
		self.buttonsimplisma.grid(column=0, row=3, columnspan=2, rowspan=1)
		
		self.buttonsimplismasave= tkinter.Button(self.AnalzeFrame, text='Save', width=20, height=1,command=self.simplisma_save)
		self.buttonsimplismasave.grid(column=2, row=3, columnspan=1, rowspan=1)
		
		simplisma_numer_label = tkinter.Label(self.AnalzeFrame, text='N Components',anchor='w')
		simplisma_numer_label.grid(column=0, row=4, columnspan=1, rowspan=1, sticky='W')
		self.simplisma_number = tkinter.IntVar()
		self.simplisma_number_entry = tkinter.Entry(self.AnalzeFrame, textvariable=self.simplisma_number, state='normal')
		self.simplisma_number_entry.grid(column=1, columnspan=2, row=4, stick='EW')
		
		self.buttonloadcomponents= tkinter.Button(self.AnalzeFrame, text='Load Components', width=20, height=1,command=self.loadcomponents)
		self.buttonloadcomponents.grid(column=0, row=5, columnspan=1, rowspan=1)
		
		self.buttonlcf= tkinter.Button(self.AnalzeFrame, text='LCF', width=20, height=1,command=self.lcf)
		self.buttonlcf.grid(column=1, row=5, columnspan=1, rowspan=1)
		
		self.buttonlcfsave= tkinter.Button(self.AnalzeFrame, text='Save LCF', width=20, height=1,command=self.lcfsave)
		self.buttonlcfsave.grid(column=2, row=5, columnspan=1, rowspan=1)
		
		constrain_label = tkinter.Label(self.AnalzeFrame, text='LCF sum to unity',anchor='w')
		constrain_label.grid(column=0, row=6, columnspan=2, rowspan=1, sticky='W')
		self.constrainusevar = tkinter.IntVar()
		self.constrainusevar.set(0)
		self.constrainusevarCheck = tkinter.Checkbutton(self.AnalzeFrame, variable=self.constrainusevar)
		self.constrainusevarCheck.grid(column=2, row=6)
		
		self.buttonmcr= tkinter.Button(self.AnalzeFrame, text='MCR', width=40, height=1,command=self.mcr)
		self.buttonmcr.grid(column=0, row=10, columnspan=2, rowspan=1)
		
		self.buttonmcrsave= tkinter.Button(self.AnalzeFrame, text='Save MCR', width=20, height=1,command=self.savemcrresults)
		self.buttonmcrsave.grid(column=2, row=10, columnspan=1, rowspan=1)
		
		constrain_nn_label = tkinter.Label(self.AnalzeFrame, text='Non negativity',anchor='w')
		constrain_nn_label.grid(column=0, row=11, columnspan=2, rowspan=1, sticky='W')
		self.constrain_nn_usevar = tkinter.IntVar()
		self.constrain_nn_usevar.set(1)
		self.constrain_nn_usevarCheck = tkinter.Checkbutton(self.AnalzeFrame, variable=self.constrain_nn_usevar)
		self.constrain_nn_usevarCheck.grid(column=2, row=11)
		
		constrain_us_label = tkinter.Label(self.AnalzeFrame, text='MCR sum to unity',anchor='w')
		constrain_us_label.grid(column=0, row=12, columnspan=2, rowspan=1, sticky='W')
		self.constrain_us_usevar = tkinter.IntVar()
		self.constrain_us_usevar.set(1)
		self.constrain_us_usevarCheck = tkinter.Checkbutton(self.AnalzeFrame, variable=self.constrain_us_usevar)
		self.constrain_us_usevarCheck.grid(column=2, row=12)
		
################################################################################
	def fileminus(self):        
		i = self.listbox5.curselection()[0]
		self.listbox5.delete(i)
		del self.files_long_analysis[i]
		self.update()
		
################################################################################        
	def fileplus(self):
		self.askloadfile()
		
################################################################################
	def askloadfile(self):
		file = askopenfile(filetypes = (("dat files","*.dat"),("all files", "*.*")))
		try:
			dirname, file_short = os.path.split(os.path.abspath(file.name))
			self.files_long_analysis.append(file)
			print('Loading File')
			self.listbox5.insert(END, file_short)
		except:
			print('Aborted Load')
			
################################################################################	
	def onselect(self, evt):	
		if 'Analysis' in App.get_selected(self):
			data_file_long = self.files_long_analysis[self.listbox5.curselection()[0]].name
			print(data_file_long)
			self.folder, self.data_file = os.path.split(os.path.abspath(data_file_long))
			self.file_type=(os.path.splitext(self.data_file)[-1])
			print('File Type :', self.file_type)
			
			if ('.dat' in self.file_type) or ('.txt' in self.file_types):
				print('Loading file ... ')
				self.data = pd.read_csv(data_file_long, sep='\t', header=0)
				self.data.dropna(how='any', inplace = True)
				print('File loaded')
				self.column_names = self.data.columns.values.tolist()
				
			if 'E' in self.column_names:
				self.energy = self.data['E']
				self.data.drop(['E'], axis=1, inplace=True)
				self.D = np.asarray((self.data).values).T
				
			elif 'Energy' in self.column_names:
				self.energy = self.data['Energy']
				self.data.drop(['Energy'], axis=1, inplace=True)
				self.D = np.asarray((self.data).values).T
				
			elif ('mcr' in self.data_file) or ('MCR' in self.data_file):
				self.D = np.asarray((self.data).values)
				
			else: 
				self.energy = self.data[self.column_names[0]]
				self.data.drop([self.column_names[0]], axis=1, inplace=True)
				self.D = np.asarray((self.data).values).T
			
################################################################################
	def loadcomponents(self):
		file = askopenfile(filetypes = (("dat files","*.dat"),("all files", "*.*")))
		try:
			dirname, file_short = os.path.split(os.path.abspath(file.name))
			print('Loading File')
			self.components = pd.read_csv(file.name, sep='\t', header=0)
			self.components_column_names = self.components.columns.values.tolist()
			
			self.clear_figure()
			self.ax=self.fig.add_subplot(1, 1, 1)
			for i in range(len(self.components_column_names)-1):
				self.ax.plot(self.components[self.components_column_names[0]], self.components[self.components_column_names[i+1]])
		
			if 'Energy' in self.components_column_names:
				self.energy = self.components['Energy']
				self.S = np.asarray(self.components.drop(['Energy'], axis=1).values)
			elif 'energy' in self.components_column_names:
				self.energy = self.components['energy']
				self.S = np.asarray(self.components.drop(['energy'], axis=1).values)
			elif 'e' in self.components_column_names:
				self.energy = self.components['e']
				self.S = np.asarray(self.components.drop(['e'], axis=1).values)
			elif 'E' in self.components_column_names:
				self.energy = self.components['E']
				self.S = np.asarray(self.components.drop(['E'], axis=1).values)
			else:
				self.S = np.asarray(self.components.drop(['E'], axis=1).values)
			self.canvas.draw_idle()
			self.update()
		except:
			print('Aborted Load')
			
################################################################################					
	def clear_figure(self):
		self.fig.clear()
		self.canvas = FigureCanvasTkAgg(self.fig, self.canvasFrame)
		self.canvas.get_tk_widget().grid(column=0, row=0, sticky='EW')
		
		self.toolbar.destroy()
		self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
		self.toolbar.update()
		
################################################################################				
	def update_figure(self):
		self.clear_figure()
		
################################################################################	
	def svd(self):
		print('Performing SVD')
		
		eigens, explained_variance_ratio = svd.svd(self.D, np.int(self.svd_number.get()))
		i = np.asarray(range(len(eigens)))+1
		
		self.clear_figure()
		self.ax1=self.fig.add_subplot(2, 1, 1)
		self.ax1.plot(i, eigens,  'o-')
		self.ax1.set_ylabel('eigenvalues')
		self.ax1.set_xlabel('components')
		self.ax2=self.fig.add_subplot(2, 1, 2)
		self.ax2.plot(i, np.cumsum(explained_variance_ratio),  'o-')
		self.ax2.set_ylabel('explained variance')
		self.ax2.set_xlabel('components')
		self.canvas.draw_idle()
		
		self.update()
		
################################################################################	
	def simplisma(self):
		print('Performing SIMPLISMA')
		
		self.S, self.C_u, self.C_c, self.LOF = simplisma.pure(self.D.T, self.simplisma_number.get(), 5, True, True)
		
		self.clear_figure()
		self.ax1=self.fig.add_subplot(4, 1, 1)
		self.ax1.plot(self.S)
		
		self.ax2=self.fig.add_subplot(4, 1, 2)
		for i in range(self.simplisma_number.get()):
			self.ax2.plot(self.C_u[:,i])
		self.ax2.set_ylabel('Unconstrained')
		
		self.ax3=self.fig.add_subplot(4, 1, 3)
		for i in range(self.simplisma_number.get()):
			self.ax3.plot(self.C_c[:,i])
		self.ax3.set_ylabel('Constrained')
		
		self.ax4=self.fig.add_subplot(4, 1, 4)
		self.ax4.plot(self.LOF)
		self.ax4.set_ylabel('LOF (%)')

		self.canvas.draw_idle()
		
		self.update()
		
		print()
		
################################################################################			
	def simplisma_save(self):
		print('saving components')
		
		save_s = pd.DataFrame(self.S)
		try:
			save_s.insert(0, 'Energy', self.energy)
			save_s.to_csv(self.folder+'\simplisma_components_'+str(self.simplisma_number.get())+'.dat', sep='\t', index=False)
		except:
			save_s.to_csv(self.folder+'\simplisma_components_'+str(self.simplisma_number.get())+'.dat', sep='\t', index=False)
		
################################################################################	
	def lcf(self):
		print('Performing LCF')
		
		self.C_u = np.dot(self.D, np.linalg.pinv(self.S.T))
		
		if self.constrainusevar.get() == 1:
			constrain = True
			self.C_c, self.LOF = lcf.lcf(self.D.T, self.S, constrain)
		else:
			constrain = False
			self.C_u, self.LOF = lcf.lcf(self.D.T, self.S, constrain)
		
		self.clear_figure()
		if constrain == False:
			self.ax1=self.fig.add_subplot(3, 1, 1)
		if constrain == True:
			self.ax1=self.fig.add_subplot(4, 1, 1)
		self.ax1.plot(self.S)
		
		if constrain == False:
			self.ax2=self.fig.add_subplot(3, 1, 2)
		if constrain == True: 
			self.ax2=self.fig.add_subplot(4, 1, 2)
		for i in range(self.S.shape[1]):
			self.ax2.plot(self.C_u[:,i])
		self.ax2.set_ylabel('Unconstrained')
		
		if constrain == True:
			self.ax3=self.fig.add_subplot(4, 1, 3)
			for i in range(self.S.shape[1]):
				self.ax3.plot(self.C_c[:,i])
			self.ax3.set_ylabel('Constrained')
			
		if constrain == True:
			self.ax4=self.fig.add_subplot(4, 1, 4)
			self.ax4.plot(self.LOF)
			self.ax4.set_ylabel('LOF (%)')
			
		if constrain == False:
			self.ax3=self.fig.add_subplot(3, 1, 3)
			self.ax3.plot(self.LOF)
			self.ax3.set_ylabel('LOF (%)')
	
		self.canvas.draw_idle()
		
		self.update()
		
		print()
		
################################################################################	
	def lcfsave(self):
		print('Save LCF')
		
		if self.constrainusevar.get() == 1:
			Concentrations_df = pd.DataFrame(self.C_c)
		else:
			Concentrations_df = pd.DataFrame(self.C_u)

		#Concentrations_df['LOF'] = self.LOF
		
		Concentrations_df.to_csv(self.folder+'\lcf_concetrations.dat', sep='\t', index=True)
		
		print('File Saved')
	
################################################################################	
	def mcr(self):
		print('Performing MCR')
		
		if (self.constrain_nn_usevar.get() == 1) & (self.constrain_us_usevar.get() == 1):
			self.mcrals = McrAls(max_iter=25, st_regr='NNLS', c_regr='NNLS', 
                c_constraints=[ConstraintNonneg(), ConstraintNorm()])
			print('contraining norm and non negativity')
		if (self.constrain_nn_usevar.get() == 1) & (self.constrain_us_usevar.get() == 0):
			self.mcrals = McrAls(max_iter=25, st_regr='NNLS', c_regr='NNLS', 
                c_constraints=[ConstraintNonneg()])
			print('contraining non negativity')
		if (self.constrain_nn_usevar.get() == 0) & (self.constrain_us_usevar.get() == 1):
			self.mcrals = McrAls(max_iter=25, st_regr='NNLS', c_regr='NNLS', 
                c_constraints=[ConstraintNorm()])
			print('contraining norm')
		if (self.constrain_nn_usevar.get() == 0) & (self.constrain_us_usevar.get() == 0):
			self.mcrals = McrAls(max_iter=25, st_regr='NNLS', c_regr='NNLS', 
                c_constraints=[])
			print('no contraints')

		self.mcrals.fit(self.D, ST=self.S.T, verbose=True)
		print('\nFinal MSE: {:.7e}'.format(self.mcrals.err[-1]))
		
		self.CS = np.dot(self.mcrals.C_opt_, self.mcrals.ST_opt_).T
		
		self.LOF = np.zeros(np.shape(self.D.T)[1])
		for i in range(np.shape(self.D.T)[1]):
			self.LOF[i] = 100*np.sqrt(np.sum(self.D.T[:,i]-self.CS[:,i])**2/np.sum(self.D.T[:,i])**2)

		self.clear_figure()
		self.ax1=self.fig.add_subplot(3, 1, 1)
		try:
			self.ax1.plot(self.energy, self.mcrals.ST_opt_.T)
			self.ax1.set_xlabel('Energy (eV)')
		except:
			self.ax1.plot(self.mcrals.ST_opt_.T)
		self.ax1.set_ylabel('Normalised Absorption')

		self.ax2=self.fig.add_subplot(3, 1, 2)
		self.ax2.plot(self.mcrals.C_opt_)
		self.ax2.set_ylabel('Component Fraction')
		self.ax2.set_xlabel('Dataset')
		
		self.ax3=self.fig.add_subplot(3, 1, 3)
		self.ax3.plot(self.LOF)
		self.ax3.set_ylabel('Lack of Fit (%)')
		self.ax3.set_xlabel('Dataset')
		
		self.canvas.draw_idle()
		
		self.update()
		
################################################################################	
	def savemcrresults(self):
		print('Saving MCR Results')
		
		Spectra_df = pd.DataFrame(self.mcrals.ST_opt_.T)
		Concentrations_df = pd.DataFrame(self.mcrals.C_opt_)
		Concentrations_df['LOF'] = self.LOF
		
		Concentrations_df.to_csv(self.folder+'\mcr_concetrations.dat', sep='\t', index=True)
		
		try:
			Spectra_df.insert(0, 'Energy', self.energy)
			Spectra_df.to_csv(self.folder+'\mcr_spectra.dat', sep='\t', index=False)
		except:
			Spectra_df.to_csv(self.folder+'\mcr_spectra.dat', sep='\t', index=False)
			
		print('Files Saved')
				
################################################################################
if __name__ == "__main__":
	app = App()
	app.title('ProXAS Gui')
	app.mainloop()
