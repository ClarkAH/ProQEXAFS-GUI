#Adam Clark UCL 2016
#adam.clark.14@ucl.ac.uk

import tkinter, subprocess, time, sys, shutil, os, wx, matplotlib,  math, random, numpy as np, pandas as pd
from distutils.spawn import find_executable
matplotlib.use('wXAgg')
import pylab as plt
import matplotlib.backends.tkagg as tkagg
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2TkAgg)
from tkinter import filedialog
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox
import warnings, larch
warnings.simplefilter(action='ignore', category=FutureWarning)
from larch import Interpreter
from larch_plugins.xafs import xftf, ftwindow, pre_edge, autobk, xftr
from scipy.signal import savgol_filter
from matplotlib import gridspec
import matplotlib.ticker as ticker
from larch import ValidateLarchPlugin, parse_group_args
from larch.utils import complex_phase
from larch_plugins.xafs import set_xafsGroup
from larch_plugins.xafs.cauchy_wavelet import cauchy_wavelet
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

class ftgui_tk(tkinter.Tk):
	def __init__(self,parent):
		tkinter.Tk.__init__(self,parent)
		self.parent = parent
		self.initialize()
	
	def initialize(self):
		self.grid()

################################################################################ 
		
		self.labelvar = tkinter.StringVar()
		labelstatus = tkinter.Label(self, textvariable=self.labelvar,
			anchor='w')
		labelstatus.grid(column=4, row=8, columnspan=1, rowspan=2, sticky='EW')
		
		self.errortype = tkinter.StringVar()
		labelerror = tkinter.Label(self, textvariable=self.errortype,
			anchor='w')
		labelerror.grid(column=0, row=25, columnspan=3, sticky='EW')
		
################################################################################
		
		buttonload = tkinter.Button(self, text=u'Load File', width=15, height=1,
			command=self.loadfile)
		buttonload.grid(column=0, row=1, columnspan=1)
		
		labele0entry = tkinter.Label(self, anchor='w', 
			text=u'E0')
		labele0entry.grid(column=0, row=3, columnspan=2, sticky='EW')
		
		self.e0var = tkinter.StringVar()
		e0entry = tkinter.Entry(self, textvariable=self.e0var)
		e0entry.grid(column=1, row=3, sticky='EW')
		
		buttoneplot = tkinter.Button(self, text=u'Plot E', width=15, height=1,
			command=self.plote)
		buttoneplot.grid(column=0, row=2, columnspan=1)
		
		buttonbplot = tkinter.Button(self, text=u'Process', width=15, height=1,
			command=self.plotback)
		buttonbplot.grid(column=1, row=2, columnspan=1)
		
		buttonsmap = tkinter.Button(self, text=u'Surface Map', width=15, height=1,
			command=self.smap)
		buttonsmap.grid(column=2, row=2, columnspan=1)
		
		buttonfitgui = tkinter.Button(self, text=u'FITGUI', width=15, height=1,
			command=self.fitgui)
		buttonfitgui.grid(column=2, row=1, columnspan=1)
		
################################################################################

		self.grid_columnconfigure(0, weight=1)
		self.resizable(True, False)
	
################################################################################

	def goplot(self):
		print('deleting plot')
		self.plot()	
		
################################################################################	

	def clicked(self, event):
		
		mx=389/(self.xlim[1]-self.xlim[0])
		cx=69-(mx*self.xlim[0])
		
		my=389/(self.ylim[0]-self.ylim[1])
		cy=43-(my*self.ylim[1])
		
		x = (self.canvas.canvasx(event.x)-cx)/mx
		y = (self.canvas.canvasy(event.y)-cy)/my
		
		self.labelvar.set('%.3f' % x)
		
################################################################################	
		
	def plote(self):
	
		root = tkinter.Tk()
		root.wm_title("Embedding in TK")
		fig = plt.Figure()
		
		axis_color = 'lightgoldenrodyellow'
		
		canvas = FigureCanvasTkAgg(fig, root)
		canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
		
		toolbar = NavigationToolbar2TkAgg(canvas, root)
		toolbar.update()
		
		ax=fig.add_subplot(111)
		fig.subplots_adjust(bottom=0.25)
		
		[line] = ax.plot(data_matrix['E'], data_matrix['0'], linewidth=2, color='red')
		
		ax.set_ylim([0,1.1*dm_max])
		
		ax_time = fig.add_axes([0.12, 0.1, 0.78, 0.03])
		s_time = Slider(ax_time, 'Dataset', 0, data_matrix.shape[1]-2, valinit=0, valstep=1, valfmt='%1.0f')
		
		def update(val):
			line.set_ydata(data_matrix[str(int(s_time.val))])
			fig.canvas.draw_idle()
			
		s_time.on_changed(update)
		
		play_button_ax = fig.add_axes([0.85, 0.025, 0.1, 0.04])
		play_button = Button(play_button_ax, 'Play', color=axis_color, hovercolor='0.975')
		def play_button_on_clicked(mouse_event):
			for i in range(data_matrix.shape[1]-1):
				s_time.set_val(i)
				line.set_ydata(data_matrix[str(int(s_time.val))])
				fig.canvas.draw()
				slt = 1/data_matrix.shape[1]
				plt.pause(slt)
				
		play_button.on_clicked(play_button_on_clicked)
		plt.tight_layout()
		
		tkinter.mainloop()

################################################################################	
	def fitgui(self):
		print('Launching FITGUI')	
		proc = subprocess.Popen(['python', 'EXAFS_FIT-13.py'])
		
################################################################################	
	def plotback(self):	
		root = tkinter.Tk()
		root.wm_title("Embedding in TK")
		fig = plt.Figure(figsize=(15,8))
		
		axis_color = 'lightgoldenrodyellow'
		
		canvas = FigureCanvasTkAgg(fig, root)
		canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
		
		toolbar = NavigationToolbar2TkAgg(canvas, root)
		toolbar.update()
		
		ax=fig.add_subplot(131)
		fig.subplots_adjust(bottom=0.45)
		
		mylarch = Interpreter()
		global sam
		sam =  larch.Group(name='sam')
		sam.energy = data_matrix['E']
		sam.mu = data_matrix['0']
		global sam_e0
		sam_e0 = float(self.e0var.get())
		vkw = 1
		
		kmin = 2.9
		kmax = 9
		fkw = 2
		
		pre_edge(sam, group=sam, e0=sam_e0, pre1=np.min(sam.energy)-sam_e0, pre2=-40, norm1=sam_e0+150, norm2=np.max(sam.energy) - sam_e0, nnorm=0, _larch=mylarch)
		autobk(sam.energy, sam.flat, rbkg=1, group=sam, clamp_lo=10, clamp_hi=30, kweight=vkw, dk=0.1, _larch=mylarch)
		
		xftf(sam.k, sam.chi, kmin=kmin, kmax=kmax, dk=2, window='hanning', kweight=fkw, group=sam, _larch=mylarch)
		xftr(sam.r, sam.chir, group=sam, qmax_out=kmax, rmin=1.2, rmax=6, window='hanning', dr=0, _larch=mylarch)
		
		kwin = ftwindow(sam.k, xmin=kmin, xmax=kmax, dx=1.5, dx2=1.5, window='hanning')
		rwin = ftwindow(sam.r, xmin=0.5, xmax=6, dx=0, dx2=0, window='hanning')
		
		[line] = ax.plot(sam.energy, sam.flat, linewidth=2, color='red')
		[line2] = ax.plot(sam.energy, sam.bkg, linewidth=2, color='black')
		
		ax.set_ylim([0,1.1*dm_max])
		
		ax_lclamp = fig.add_axes([0.1, 0.35, 0.2, 0.03])
		s_lclamp = Slider(ax_lclamp, 'Low Clamp', 0, 50, valinit=10)
		
		ax_hclamp = fig.add_axes([0.1, 0.3, 0.2, 0.03])
		s_hclamp = Slider(ax_hclamp, 'High Clamp', 0, 50, valinit=30)
		
		ax_rbkg = fig.add_axes([0.1, 0.25, 0.2, 0.03])
		s_rbkg = Slider(ax_rbkg, 'RBKG', 0.5, 1.5, valinit=1)
		
		ax_e0 = fig.add_axes([0.1, 0.2, 0.2, 0.03])
		s_e0 = Slider(ax_e0, 'E0', sam_e0-10, sam_e0+10, valinit=sam_e0)
		
		ax_kmin = fig.add_axes([0.4, 0.35, 0.2, 0.03])
		s_kmin = Slider(ax_kmin, 'K min', 0, np.max(sam.k), valinit=3)
		
		ax_kmax = fig.add_axes([0.4, 0.3, 0.2, 0.03])
		s_kmax = Slider(ax_kmax, 'K max', 0, np.max(sam.k), valinit=10)
		
		ax_dk = fig.add_axes([0.4, 0.25, 0.2, 0.03])
		s_dk = Slider(ax_dk, 'dK', 0, 4, valinit=1.5)
		
		ax_rmin = fig.add_axes([0.7, 0.35, 0.2, 0.03])
		s_rmin = Slider(ax_rmin, 'R min', 0, np.max(sam.r), valinit=0.5)
		
		ax_rmax = fig.add_axes([0.7, 0.3, 0.2, 0.03])
		s_rmax = Slider(ax_rmax, 'R max', 0, np.max(sam.r), valinit=6)
		
		ax_dr = fig.add_axes([0.7, 0.25, 0.2, 0.03])
		s_dr = Slider(ax_dr, 'dR', 0, 4, valinit=0)
		
		rax = fig.add_axes([0.1, 0.05, 0.05, 0.1], facecolor=axis_color)
		radio = RadioButtons(rax, ('1', '2', '3'))
		
		ray = fig.add_axes([0.4, 0.05, 0.05, 0.1], facecolor=axis_color)
		radio2 = RadioButtons(ray, ('1', '2', '3'),  active=1)
		
		ax2=fig.add_subplot(132)
		fig.subplots_adjust(bottom=0.45)
		
		[linek] = ax2.plot(sam.k, sam.chi*(sam.k**fkw), linewidth=2, color='red')
		
		[line_kmin] = ax2.plot([3, 3], [-1,1], linewidth = 2, color='k', ls='--')
		[line_kmax] = ax2.plot([10, 10], [-1,1], linewidth = 2, color='k', ls='--')
		[line_z] = ax2.plot([0, np.max(sam.k)], [0,0], linewidth = 1, color='black')
		[line_kwin] = ax2.plot(sam.k, kwin, linewidth = 2, color = 'green')
		
		ax3=fig.add_subplot(133)
		fig.subplots_adjust(bottom=0.45)
		
		[liner] = ax3.plot(sam.r, sam.chir_mag, linewidth=2, color='red')
		[line_rwin] = ax3.plot(sam.r, rwin, linewidth = 2, color = 'green')

		ax_time = fig.add_axes([0.1, 0.95, 0.8, 0.03])
		s_time = Slider(ax_time, 'Dataset', 0, data_matrix.shape[1]-2, valinit=0, valstep=1, valfmt='%1.0f')
		
		def update(val):
			global sam_e0
			sam_e0 = s_e0.val
			lcmp = s_lclamp.val
			hcmp = s_hclamp.val
			vrbkg = s_rbkg.val
			fkw = int(radio2.value_selected)
			vkw = int(radio.value_selected)
			kminv = s_kmin.val
			kmaxv = s_kmax.val
			dkv = s_dk.val
			rminv = s_rmin.val
			rmaxv = s_rmax.val
			drv = s_dr.val
			
			autobk(sam.energy, sam.flat, rbkg=vrbkg, group=sam, clamp_lo=lcmp, clamp_hi=hcmp, kweight=vkw, kstep = 0.025, dk=0.1, _larch=mylarch)
			xftf(sam.k, sam.chi, kmin=kminv, kmax=kmaxv, dk=dkv, window='hanning', kweight=fkw, group=sam, _larch=mylarch)
			xftr(sam.r, sam.chir, group=sam, qmax_out=kmax, rmin=rminv, rmax=rmaxv, window='hanning', dr=drv, _larch=mylarch)
			kwin = ftwindow(sam.k, xmin=kminv, xmax=kmaxv, dx=dkv, dx2=dkv, window='hanning')
			rwin = ftwindow(sam.r, xmin=rminv, xmax=rmaxv, dx=drv, dx2=drv, window='hanning')
			
			linek.set_xdata(sam.k)
			linek.set_ydata(sam.chi*(sam.k**fkw))
			liner.set_ydata(sam.chir_mag)
			line2.set_ydata(sam.bkg)
			line.set_ydata(sam.flat)
			line_kmax.set_xdata([kmaxv,kmaxv])
			line_kmin.set_xdata([kminv,kminv])
			line_kwin.set_ydata(kwin)
			line_kwin.set_xdata(sam.k)
			line_rwin.set_ydata(rwin)
			
			ax.relim()
			ax.autoscale_view()
			ax2.relim()
			ax2.autoscale_view()
			ax3.relim()
			ax3.autoscale_view()
			fig.canvas.draw()
			
		def updateds(val):
			sam_e0 = s_e0.val
			sam.mu = data_matrix[str(int(s_time.val))]
			pre_edge(sam, group=sam, e0=sam_e0, pre1=np.min(sam.energy)-sam_e0, pre2=-40, norm1=sam_e0+150, norm2=np.max(sam.energy) - sam_e0, nnorm=0, _larch=mylarch)
			update(1)
			
		def submit(val):
			global sam_e0
			sam_e0 = s_e0.val
			pre_edge(sam, group=sam, e0=sam_e0, pre1=np.min(sam.energy)-sam_e0, pre2=-40, norm1=sam_e0+150, norm2=np.max(sam.energy) - sam_e0, nnorm=0, _larch=mylarch)
			update(1)
			
		s_e0.on_changed(submit)
		s_time.on_changed(updateds)
		radio.on_clicked(update)
		radio2.on_clicked(update)
		s_lclamp.on_changed(update)
		s_hclamp.on_changed(update)
		s_rbkg.on_changed(update)
		s_lclamp.on_changed(update)
		s_kmax.on_changed(update)
		s_kmin.on_changed(update)
		s_dk.on_changed(update)
		s_rmax.on_changed(update)
		s_rmin.on_changed(update)
		s_dr.on_changed(update)
		
		WT_button_ax = fig.add_axes([0.75, 0.025, 0.05, 0.04])
		WT_button = Button(WT_button_ax, 'Wavelet', color=axis_color, hovercolor='0.975')
		def wt_button_on_clicked(mouse_event):
			root2 = tkinter.Tk()
			root2.wm_title("Embedding in TK")
			fig2 = plt.Figure(figsize=(15,6))

			canvas2 = FigureCanvasTkAgg(fig2, root2)
			canvas2.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
			toolbar = NavigationToolbar2TkAgg(canvas2, root2)
			toolbar.update()
			
			kminv = s_kmin.val
			kmaxv = s_kmax.val
			rminv = s_rmin.val
			rmaxv = s_rmax.val
			dkv = s_dk.val
			kwin = ftwindow(sam.k, xmin=kminv, xmax=kmaxv, dx=dkv, dx2=dkv, window='hanning')
			global cwsam
			cwsam =  larch.Group(name='cwsam')
			cwsam.k = sam.k
			cwsam.chi = sam.chi * kwin
			cwsam.r = sam.r

			kwwt = int(radio2.value_selected)
			cauchy_wavelet(cwsam, kweight=kwwt, rmax_out = rmaxv, _larch=mylarch)

			levels_c1 = np.linspace(np.min(cwsam.wcauchy_mag.T),np.max(cwsam.wcauchy_mag.T),25)
			levels_c2 = np.linspace(np.min(cwsam.wcauchy_re.T),np.max(cwsam.wcauchy_re.T),25)
			
			gs = gridspec.GridSpec(1, 2, width_ratios=[1,1]) 
			
			ax0 = fig2.add_subplot(gs[0])
			cax1 = make_axes_locatable(ax0).append_axes("right", size="5%", pad="2%")
			
			c1 = ax0.contourf(cwsam.r,cwsam.k,cwsam.wcauchy_mag.T,levels_c1,cmap='coolwarm')
			cb1 = plt.colorbar(c1, cax=cax1)
			cb1_ticks = np.linspace(np.min(levels_c1), np.max(levels_c1), num=11, endpoint=True)
			cb1.set_ticks(cb1_ticks) 
			cb1.draw_all() 
			
			
			ax1 = fig2.add_subplot(gs[1])
			cax2 = make_axes_locatable(ax1).append_axes("right", size="5%", pad="2%")
			
			c2 = ax1.contourf(cwsam.r,cwsam.k,cwsam.wcauchy_re.T,levels_c2,cmap='coolwarm')
			cb2 = plt.colorbar(c2, cax=cax2)
			cb2_ticks = np.linspace(np.min(levels_c2), np.max(levels_c2), num=11, endpoint=True)
			cb2.set_ticks(cb2_ticks) 
			cb2.draw_all() 
			ax0.set_ylim([kminv,kmaxv+(0.5*dkv)])
			ax1.set_ylim([kminv,kmaxv+(0.5*dkv)])
			ax0.set_xlim([rminv,rmaxv])
			ax1.set_xlim([rminv,rmaxv])
			ax0.set_title("Wavelet Transform Magnitude")
			ax1.set_title("Wavelet Transform Real Component")
			ax0.set_ylabel('k')
			ax0.set_xlabel('R')
			ax1.set_ylabel('k')
			ax1.set_xlabel('R')
			
			ax_time2 = fig2.add_axes([0.1, 0.95, 0.8, 0.03])
			s_time2 = Slider(ax_time2, 'Dataset', 0, data_matrix.shape[1]-2, valinit=int(s_time.val), valstep=1, valfmt='%1.0f')
			
			fig2.subplots_adjust(left=0.05, bottom=0.12, right=0.92, top=0.90, wspace=0.23, hspace=0.2)
			
			def updateds2(val):
				sam_e0 = s_e0.val
				lcmp = s_lclamp.val
				hcmp = s_hclamp.val
				vrbkg = s_rbkg.val
				fkw = int(radio2.value_selected)
				vkw = int(radio.value_selected)
				kminv = s_kmin.val
				kmaxv = s_kmax.val
				dkv = s_dk.val
				rminv = s_rmin.val
				rmaxv = s_rmax.val
				drv = s_dr.val
				rminv = s_rmin.val
				rmaxv = s_rmax.val
				
				cwsam.mu = data_matrix[str(int(s_time2.val))]
				cwsam.energy = sam.energy
				pre_edge(cwsam, group=cwsam, e0=sam_e0, pre1=np.min(cwsam.energy)-sam_e0, pre2=-40, norm1=sam_e0+150, norm2=np.max(cwsam.energy) - sam_e0, nnorm=0, _larch=mylarch)
			
				autobk(cwsam.energy, cwsam.flat, rbkg=vrbkg, group=cwsam, clamp_lo=lcmp, clamp_hi=hcmp, kweight=vkw, kstep = 0.025, dk=0.1, _larch=mylarch)
				xftf(cwsam.k, cwsam.chi, kmin=kminv, kmax=kmaxv, dk=dkv, window='hanning', kweight=fkw, group=cwsam, _larch=mylarch)
				xftr(cwsam.r, cwsam.chir, group=cwsam, qmax_out=kmax, rmin=rminv, rmax=rmaxv, window='hanning', dr=drv, _larch=mylarch)
						
				kwin = ftwindow(cwsam.k, xmin=kminv, xmax=kmaxv, dx=dkv, dx2=dkv, window='hanning')
				cwsam.chi = cwsam.chi * kwin

				kwwt = int(radio2.value_selected)
				cauchy_wavelet(cwsam, kweight=kwwt, rmax_out = rmaxv, _larch=mylarch)
				
				ax0.cla()
				ax1.cla()
				cax1.clear()
				cax2.clear()
				
				levels_c1 = np.linspace(np.min(cwsam.wcauchy_mag.T),np.max(cwsam.wcauchy_mag.T),25)
				levels_c2 = np.linspace(np.min(cwsam.wcauchy_re.T),np.max(cwsam.wcauchy_re.T),25)
				c1 = ax0.contourf(cwsam.r,cwsam.k,cwsam.wcauchy_mag.T,levels_c1,cmap='coolwarm')
				c2 = ax1.contourf(cwsam.r,cwsam.k,cwsam.wcauchy_re.T,levels_c2,cmap='coolwarm')
				
				cb1 = plt.colorbar(c1, cax=cax1)
				cb2 = plt.colorbar(c2, cax=cax2)

				ax0.set_ylim([kminv-(0.5*dkv),kmaxv+(0.5*dkv)])
				ax1.set_ylim([kminv-(0.5*dkv),kmaxv+(0.5*dkv)])
				ax0.set_xlim([rminv,rmaxv])
				ax1.set_xlim([rminv,rmaxv])
				ax0.set_title("Wavelet Transform Magnitude")
				ax1.set_title("Wavelet Transform Real Component")
				ax0.set_ylabel('k')
				ax0.set_xlabel('R')
				ax1.set_ylabel('k')
				ax1.set_xlabel('R')
				fig2.canvas.draw_idle()
				
			def updateds3(val):
				s_time2.set_val(s_time.val)
				updateds2(1)
			
			s_time.on_changed(updateds3)
			s_time2.on_changed(updateds2)
			s_e0.on_changed(updateds2)
			s_time.on_changed(updateds2)
			radio.on_clicked(updateds2)
			radio2.on_clicked(updateds2)
			s_lclamp.on_changed(updateds2)
			s_hclamp.on_changed(updateds2)
			s_rbkg.on_changed(updateds2)
			s_lclamp.on_changed(updateds2)
			s_kmax.on_changed(updateds2)
			s_kmin.on_changed(updateds2)
			s_dk.on_changed(updateds2)
			s_rmax.on_changed(updateds2)
			s_rmin.on_changed(updateds2)
			
			play_wtbutton_ax = fig2.add_axes([0.85, 0.025, 0.05, 0.04])
			play_wtbutton = Button(play_wtbutton_ax, 'Export', color=axis_color, hovercolor='0.975')
			def play_wtbutton_on_clicked(mouse_event):
				dfrr = pd.DataFrame()
				dfri = pd.DataFrame()
				dfrm = pd.DataFrame()
				dfk = pd.DataFrame()
				dfq = pd.DataFrame()
				if not os.path.exists(dname+'/WT'):
					os.makedirs(dname+'/WT')
				for i in range(data_matrix.shape[1]-1):
					print(i)
					s_time.set_val(i)
					root.update()
					slt = 1/(data_matrix.shape[1]-1)
					wtnumber = "%05d" % float(i)
					fig2.savefig(dname+'/WT/WT_transform_'+wtnumber+'.png', bbox_inches='tight', dpi=400)
					
					if i == 0:
						dfk['k'] = sam.k
						dfrr['r'] = sam.r
						dfri['r'] = sam.r
						dfrm['r'] = sam.r
						dfq['q'] = sam.q
					
					dfk[str(i)] = sam.chi
					dfrr[str(i)] = sam.chir_re
					dfri[str(i)] = sam.chir_im
					dfrm[str(i)] = sam.chir_mag
					dfq[str(i)] = sam.chiq_re
					
					time.sleep(slt)
					
				rmcr = pd.concat([dfrr, dfri], axis=1)
				rmcr.to_csv(dname+'/R-space_ri_mcr.dat', sep='\t', index=False)
				dfrm.to_csv(dname+'/R-space_magnitude.dat', sep='\t', index=False)
				dfk.to_csv(dname+'/K-space.dat', sep='\t', index=False)
				dfq.to_csv(dname+'/Q-space.dat', sep='\t', index=False)
					
				command = 'ffmpeg'
				if find_executable(command) is not None:
					savedPath = os.getcwd()
					os.chdir(dname+'/WT')
					subprocess.call(['ffmpeg.exe', '-r', '10', '-i', 'WT_transform_'+'%05d.png', '-c:v', 'libx264', '-strict', '-2', '-pix_fmt', 'yuv420p', '-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2', '-f', 'mp4', 'WT_transform.mp4',])
					os.chdir(savedPath)
				else:
					print('ffmpeg is required!')
					print('https://ffmpeg.org/')
					print('Unzip folder and add the ffmpeg/bin folder to the system enviroment path')

			play_wtbutton.on_clicked(play_wtbutton_on_clicked)
			
			tkinter.mainloop()

		WT_button.on_clicked(wt_button_on_clicked)
		
		play_button_ax = fig.add_axes([0.85, 0.025, 0.05, 0.04])
		play_button = Button(play_button_ax, 'Export', color=axis_color, hovercolor='0.975')
		def play_button_on_clicked(mouse_event):
			dfrr = pd.DataFrame()
			dfri = pd.DataFrame()
			dfrm = pd.DataFrame()
			dfk = pd.DataFrame()
			dfq = pd.DataFrame()
			for i in range(data_matrix.shape[1]-1):
				print(i)
				s_time.set_val(i)
				sam_e0 = s_e0.val
				lcmp = s_lclamp.val
				hcmp = s_hclamp.val
				vrbkg = s_rbkg.val
				fkw = int(radio2.value_selected)
				vkw = int(radio.value_selected)
				kminv = s_kmin.val
				kmaxv = s_kmax.val
				dkv = s_dk.val
				rminv = s_rmin.val
				rmaxv = s_rmax.val
				drv = s_dr.val
				
				sam.mu = data_matrix[str(int(s_time.val))]
				pre_edge(sam, group=sam, e0=sam_e0, pre1=np.min(sam.energy)-sam_e0, pre2=-40, norm1=sam_e0+150, norm2=np.max(sam.energy) - sam_e0, nnorm=0, _larch=mylarch)
				
				autobk(sam.energy, sam.flat, rbkg=vrbkg, group=sam, clamp_lo=lcmp, clamp_hi=hcmp, kweight=vkw, kstep = 0.025, dk=0.1, _larch=mylarch)
				xftf(sam.k, sam.chi, kmin=kminv, kmax=kmaxv, dk=dkv, window='hanning', kweight=fkw, group=sam, _larch=mylarch)
				xftr(sam.r, sam.chir, group=sam, qmax_out=kmax, rmin=rminv, rmax=rmaxv, window='hanning', dr=drv, _larch=mylarch)
				kwin = ftwindow(sam.k, xmin=kminv, xmax=kmaxv, dx=dkv, dx2=dkv, window='hanning')
				rwin = ftwindow(sam.r, xmin=rminv, xmax=rmaxv, dx=drv, dx2=drv, window='hanning')
				
				linek.set_xdata(sam.k)
				linek.set_ydata(sam.chi*(sam.k**fkw))
				liner.set_ydata(sam.chir_mag)
				line2.set_ydata(sam.bkg)
				line.set_ydata(sam.flat)
				line_kmax.set_xdata([kmaxv,kmaxv])
				line_kmin.set_xdata([kminv,kminv])
				line_kwin.set_ydata(kwin)
				line_kwin.set_xdata(sam.k)
				line_rwin.set_ydata(rwin)
				
				ax.relim()
				ax.autoscale_view()
				ax2.relim()
				ax2.autoscale_view()
				ax3.relim()
				ax3.autoscale_view()
				fig.canvas.draw()
				
				if i == 0:
					dfk['k'] = sam.k
					dfrr['r'] = sam.r
					dfri['r'] = sam.r
					dfrm['r'] = sam.r
					dfq['q'] = sam.q
					
				dfk[str(i)] = sam.chi
				dfrr[str(i)] = sam.chir_re
				dfri[str(i)] = sam.chir_im
				dfrm[str(i)] = sam.chir_mag
				dfq[str(i)] = sam.chiq_re
				
				root.update()
				slt = 1/(data_matrix.shape[1]-1)
				time.sleep(slt)
			
			#dfrr.to_csv(dname+'/R-space_real.dat', sep='\t', index=False)
			#dfri.to_csv(dname+'/R-space_imaginary.dat', sep='\t', index=False)
			
			rmcr = pd.concat([dfrr, dfri], axis=1)
			rmcr.to_csv(dname+'/R-space_ri_mcr.dat', sep='\t', index=False)
			dfrm.to_csv(dname+'/R-space_magnitude.dat', sep='\t', index=False)
			dfk.to_csv(dname+'/K-space.dat', sep='\t', index=False)
			dfq.to_csv(dname+'/Q-space.dat', sep='\t', index=False)
			
			print('FT data saved')
				
		play_button.on_clicked(play_button_on_clicked)
		
		tkinter.mainloop()
		
################################################################################	
	def loadfile(self):
		self.filename =  filedialog.askopenfilename(initialdir = "/",title = "Select file",filetypes = (("dat files","*.dat"),("all files","*.*")))
		print (self.filename)
		print('Reading Data ...')
		global data_matrix
		data_matrix = pd.read_csv(self.filename, sep='\s+', header=0)
		global dm_max
		dm_max = np.max(((data_matrix.max()).values)[1:])
		global fname
		fname = (self.filename.split('/')[-1]).split('.')[0]
		global dname
		dname = '/'.join(str(s) for s in (self.filename.split('/')[:-1]))
		print(dname)
		print(fname)
		
		if not os.path.exists(dname):
			os.makedirs(dname)
		
		print('Data Read')
		
################################################################################
	def smap(self):
		
		print('Plotting Surface Map')
		root = tkinter.Tk()
		root.wm_title("Embedding in TK")
		fig = plt.Figure()
		
		canvas = FigureCanvasTkAgg(fig, root)
		canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
		
		toolbar = NavigationToolbar2TkAgg(canvas, root)
		toolbar.update()
		
		gs = gridspec.GridSpec(2, 2, height_ratios=[1,3],width_ratios=[3, 1.5]) 
		ax2 = fig.add_subplot(gs[2])
		
		x = data_matrix[data_matrix.columns[0]]
		Z = data_matrix.drop(data_matrix.columns[0], axis=1)
		y = range(Z.shape[1])
		
		if np.max(x) < 11:
			ax2.set_xlabel(r'R ($\AA$)')
		else:
			ax2.set_xlabel('Energy (eV)')
		ax2.set_ylabel('Dataset')

		levels = np.linspace(np.min(Z.min()),np.max(Z.max()),25)

		c1 = ax2.contourf(x,y,Z.T,levels,cmap='coolwarm')
		plt.colorbar(c1, ax=ax2)
		ax2.set_ylim([0,Z.shape[1]-1])

		ax3 = fig.add_subplot(gs[3])
		[c21] = ax3.plot(Z.iloc[[0]].T,y, label=str(float('%1.2f'% x[int(0)])))
		[c22] = ax3.plot(Z.iloc[[Z.shape[0]-1]].T,y, label=str(float('%1.2f' % x[Z.shape[0]-1])))
		ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
		plt.setp(ax3.xaxis.get_majorticklabels(), rotation=90)
		ax3.set_ylim([0,Z.shape[1]-1])
		
		inner = gridspec.GridSpecFromSubplotSpec(1, 2,subplot_spec=gs[0], wspace=0, hspace=0, width_ratios=[3, 0.8])
		ax0 = fig.add_subplot(inner[0])
		[c31] = ax0.plot(x,Z[Z.columns[0]], label='0')
		[c32] = ax0.plot(x,Z[Z.columns[Z.shape[1]-1]], label=str(Z.shape[1]-1))
		yll, yhl = ax0.get_ylim()
		ax0ys = np.round(yhl/5, decimals=1)
		
		ax0.yaxis.set_ticks(np.arange(0, np.max(Z.max()), ax0ys))
		ax0.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
		
		inner2 = gridspec.GridSpecFromSubplotSpec(4, 2,subplot_spec=gs[1], wspace=0, hspace=0, height_ratios=[1, 1, 1, 1], width_ratios=[3, 1])
		ax11 = fig.add_subplot(inner2[0])
		s_ds1 = Slider(ax11, '', 0, Z.shape[1]-1, valinit=0, valstep=1, valfmt='%1.0f')
		
		ax12 = fig.add_subplot(inner2[2])
		s_ds2 = Slider(ax12, '', 0, Z.shape[1]-1, valinit=Z.shape[1]-1, valstep=1, valfmt='%1.0f')
		
		ax13 = fig.add_subplot(inner2[4])
		s_e1 = Slider(ax13, '', 0, Z.shape[0]-1, valinit=0, valstep=1)
		s_e1.valtext.set_visible(False)
		s_e1_tx = ax13.text(0.5, 0, str(float('%1.2f' % x[int(s_e1.val)])), verticalalignment='bottom', horizontalalignment='right')
		
		ax14 = fig.add_subplot(inner2[6])
		s_e2 = Slider(ax14, '',0, Z.shape[0]-1, valinit=Z.shape[0]-1, valstep=1)
		s_e2.valtext.set_visible(False)
		s_e2_tx = ax14.text(0.5, 0, str(float('%1.2f' % x[int(s_e2.val)])), verticalalignment='bottom', horizontalalignment='right')
		
		[ds1_line] = ax2.plot([np.min(x), np.max(x)], [0,0], linewidth = 1, color='k')
		[ds2_line] = ax2.plot([np.min(x), np.max(x)], [Z.shape[1],Z.shape[1]], linewidth = 1, color='k')
		[e1_line] = ax2.plot([np.min(x), np.min(x)], [0,Z.shape[1]], linewidth = 1, color='k')
		[e2_line] = ax2.plot([np.max(x), np.max(x)], [0,Z.shape[1]], linewidth = 1, color='k')
		
		c31max = np.max(Z[Z.columns[0]].values)
		c32max = np.max(Z[Z.columns[Z.shape[1]-1]].values)
		ylimh = 1.1*np.max([c31max, c32max])
		
		ax3.xaxis.set_ticks(np.arange(0, np.max(Z.max()), ax0ys))
		ax0.set_xlim([np.min(x), np.max(x)])
		ax0.set_ylim(bottom=0, top=ylimh)
		ax3ll, ax3hl = ax0.get_ylim()
		ax3.set_xlim([ax3ll, ax3hl])
		
		legend = ax0.legend(loc='best')
		legend2 = ax3.legend(loc='upper center')

		def update(val):
			c21.set_xdata(Z.iloc[[int(s_e1.val)]].values)
			c22.set_xdata(Z.iloc[[int(s_e2.val)]].values)
			e1_line.set_xdata([x[int(s_e1.val)],x[int(s_e1.val)]])
			e2_line.set_xdata([x[int(s_e2.val)],x[int(s_e2.val)]])
			
			c31.set_ydata(Z[Z.columns[int(s_ds1.val)]].values)
			c32.set_ydata(Z[Z.columns[int(s_ds2.val)]].values)
			ds1_line.set_ydata([s_ds1.val,s_ds1.val])
			ds2_line.set_ydata([s_ds2.val,s_ds2.val])
			
			s_e2_tx.set_text(str(float('%1.2f' % x[int(s_e2.val)])))
			s_e1_tx.set_text(str(float('%1.2f' % x[int(s_e1.val)])))
			
			c31max = np.max(Z[Z.columns[int(s_ds1.val)]].values)
			c32max = np.max(Z[Z.columns[int(s_ds2.val)]].values)
			
			ylimh = 1.1*np.max([c31max, c32max])
			ax0ys = np.round(ylimh/5, decimals=1)
			ax0.yaxis.set_ticks(np.arange(0, np.max(Z.max()), ax0ys))
			ax0.set_ylim(bottom=0, top=ylimh)
			
			legend.get_texts()[0].set_text(str(int(s_ds1.val)))
			legend.get_texts()[1].set_text(str(int(s_ds2.val)))
			
			legend2.get_texts()[0].set_text(str(float('%1.2f' % x[int(s_e1.val)])))
			legend2.get_texts()[1].set_text(str(float('%1.2f' % x[int(s_e2.val)])))
			
			ax0.autoscale_view()
			fig.canvas.draw_idle()
		
		s_ds1.on_changed(update)
		s_ds2.on_changed(update)
		s_e1.on_changed(update)
		s_e2.on_changed(update)
		
		fig.subplots_adjust(left=0.12, bottom=0.12, right=0.9, top=0.95, wspace=0.3, hspace=0.2)
		fig.canvas.draw()
		tkinter.mainloop()

################################################################################


if __name__ == "__main__":
		app = ftgui_tk(None)
		app.title('FT Gui')
		app.mainloop()
		