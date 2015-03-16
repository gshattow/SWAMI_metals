#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

#import h5py as h5
import numpy as np
import pylab as plt
from random import sample, seed, gauss
from os.path import getsize as getFileSize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator, MaxNLocator
from matplotlib.patches import Circle, PathPatch, Rectangle
import cPickle
import random
import matplotlib.colors as colors
import matplotlib.cm as cm

# ================================================================================
# Basic variables
# ================================================================================

# Set up some basic attributes of the run

SIM = 'MM'
inner_edge = 1.
nbins =100
fHI_d = 0.01
fb = 0.17
frame_depth = 3
numsnaps = 64


if ((SIM == 'SML') | (SIM == 'RSM') | (SIM == 'SSM')) :
	box_side_h1 = 20.
	pcl_mass_h1 = 0.001616 * 1.0e10
	
if ((SIM == 'G3') | (SIM == 'AG3') | (SIM == 'SG3')) :
	box_side_h1 = 50.
	pcl_mass_h1 = 0.010266 * 1.0e10
	npcls = 464.**3 
	 
if ((SIM == 'RG4') | (SIM == 'AG4') | (SIM == 'SG4')) :
	box_side_h1 = 100.
	pcl_mass_h1 = 0.1093 * 1.0e10
	 
if (SIM == 'LRG') :
	box_side_h1 = 248.
	pcl_mass_h1 = 1.0 * 1.0e10
	
if (SIM == 'MM') :
	box_side_h1 = 62.5
	pcl_mass_h1 = 0.086 * 1.0e10
	npcls = 270.**3 

min_bin = 0
vdepth = 1 #km/s
Zsun = 0.02

font = {'family':'serif','size':20, 'serif':'Times New Roman'}

stars_color = ['#006600', '#009900', '#00ff00']
coldgas_color = ['#000099', '#0033ff', '#66ccff']
hotgas_color = ['#330066', '#6600cc', '#cc99ff'] # '#660099',
igm_color = ['#660033', '#cc0099', '#ff66ff'] #990066
WHIM_color = ['#990000', '#cc3300', '#ff6633']
DM_color = ['#000000', '#666666', '#cccccc']

line_cycle = ['-', '--', ':']

colors = [stars_color, coldgas_color, hotgas_color, igm_color, WHIM_color, DM_color]
name_cycle = ['Stars', 'Cold Gas', 'Hot Gas', 'IGM', 'WHIM']


matplotlib.rcdefaults()
plt.rc('axes', color_cycle=[
    'k',
    'b',
    'r',
    '#00FF00',
    'm',
    '0.5',
    ], labelsize='x-large')
plt.rc('xtick', labelsize='x-large')
plt.rc('ytick', labelsize='x-large')
plt.rc('lines', linewidth='2.0')
plt.rc('font', **font)
plt.rc('legend', numpoints=1, fontsize='x-large')
#plt.rc('text', usetex=True)


redshifts = [127.000, 79.998, 50.000, 30.000, 19.916, 18.244, 16.725, 15.343, 14.086, 12.941, 11.897, 
	10.944, 10.073, 9.278, 8.550, 7.883, 7.272, 6.712, 6.197, 5.724, 5.289, 4.888, 4.520, 4.179, 3.866, 3.576, 
		3.308, 3.060, 2.831, 2.619, 2.422, 2.239, 2.070, 1.913, 1.766, 1.630, 1.504, 1.386, 1.276, 1.173, 1.078, 
			0.989, 0.905, 0.828, 0.755, 0.687, 0.624, 0.564, 0.509, 0.457, 0.408, 0.362, 0.320, 0.280, 0.242, 
				0.208, 0.175, 0.144, 0.116, 0.089, 0.064, 0.041, 0.020, 0.000]

Dir = 'plots/'
OutputFormat = '_cent.png'
TRANSPARENT = False

class Results:

	def Z_Zsun_by_z(self, snapnum) :
	
		filename = SIM + '/data_files/' + 'metals_all_z_000'
		print 'reading in file:', filename
		
		zz = np.zeros(numsnaps)
		mstars = np.zeros(numsnaps)
		mmstars = np.zeros(numsnaps)
		mdiffuse = np.zeros(numsnaps)
		mmdiffuse = np.zeros(numsnaps)
		mhotgas = np.zeros(numsnaps)
		mmhotgas = np.zeros(numsnaps)
		mcoldgas = np.zeros(numsnaps)
		mmcoldgas = np.zeros(numsnaps)
		ii = 0
		for item in file(filename) :
			item = item.split()
			zz[ii] = float(item[0])
			mstars[ii] = float(item[1])
			mmstars[ii] = float(item[2])
			mdiffuse[ii] = float(item[3])
			mmdiffuse[ii] = float(item[4])
			mhotgas[ii] = float(item[5])
			mmhotgas[ii] = float(item[6])
			mcoldgas[ii] = float(item[7])
			mmcoldgas[ii] = float(item[8])
			ii = ii + 1
		

		mstars = np.where(mstars > 0., mstars, 0.1)		
		mdiffuse = np.where(mdiffuse > 0., mdiffuse, 0.1)		
		mcoldgas = np.where(mcoldgas > 0., mcoldgas, 0.1)		
		mhotgas = np.where(mhotgas > 0., mhotgas, 0.1)		

		Zstars = mmstars/mstars/Zsun
		Zdiffuse = mmdiffuse/mdiffuse/Zsun
		Zhotgas = mmhotgas/mhotgas/Zsun
		Zcoldgas = mmcoldgas/mcoldgas/Zsun


		fig = plt.figure(figsize=(12.,10))							# square box
#		fig.subplots_adjust(wspace = .0,top = .96, bottom = .1, left = .05, right = .95)
		ax = fig.add_subplot(1,1,1)
	
		plt.axis([0,5,-4,0])
#		ax.set_yscale('log')
		
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(30)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(30)

		ax.plot(zz, np.log10(Zstars), c = '#00ff00', ls = '-', lw = 4, label = 'Stars')
		ax.plot(zz, np.log10(Zcoldgas), c = 'blue', ls = '-', lw = 4, label = 'Cold Gas (ISM)')
		ax.plot(zz, np.log10(Zhotgas), c = 'purple', ls = '-', lw = 4, label = 'Hot Gas (Halo)')
		ax.plot(zz, np.log10(Zdiffuse), c = '#cc0066', ls = '-', lw = 4, label = 'Diffuse Gas (IGM&ICM)')

		ax.set_xlabel(r'z', fontsize=34) #54 for square
		ax.set_ylabel(r'$log_{10}(Z/Z_{\odot})$', fontsize=34)  #36 for stripe
		ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0)

		outputFile = Dir + SIM + '_' + 'Z_Zsun_vs_z' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()

	def Z_Zsun_by_z_by_model(self, nreheat, neject) :
	
		zz = np.zeros((nreheat, neject, numsnaps))
		mstars = np.zeros((nreheat, neject, numsnaps))
		mmstars = np.zeros((nreheat, neject, numsnaps))
		mdiffuse = np.zeros((nreheat, neject, numsnaps))
		mmdiffuse = np.zeros((nreheat, neject, numsnaps))
		mWHIM = np.zeros((nreheat, neject, numsnaps))
		mmWHIM = np.zeros((nreheat, neject, numsnaps))
		mhotgas = np.zeros((nreheat, neject, numsnaps))
		mmhotgas = np.zeros((nreheat, neject, numsnaps))
		mcoldgas = np.zeros((nreheat, neject, numsnaps))
		mmcoldgas = np.zeros((nreheat, neject, numsnaps))

		for ir in range(nreheat) :
			for ie in range(neject) :
				filename = 'W' + str(ir + 1) + str(ie + 1) + '/metals_all_z_000'
				print 'reading in file:', filename
		
				jj = 0
				for item in file(filename) :
					item = item.split()
					zz[ir][ie][jj] = float(item[0])
					mstars[ir][ie][jj] = float(item[1])
					mmstars[ir][ie][jj] = float(item[2])
					mdiffuse[ir][ie][jj] = float(item[3]) + float(item[9]) - float(item[10])
					mmdiffuse[ir][ie][jj] = float(item[4]) - float(item[11])
					mWHIM[ir][ie][jj] = float(item[10])
					mmWHIM[ir][ie][jj] = float(item[11])
					mhotgas[ir][ie][jj] = float(item[5])
					mmhotgas[ir][ie][jj] = float(item[6])
					mcoldgas[ir][ie][jj] = float(item[7])
					mmcoldgas[ir][ie][jj] = float(item[8])
					jj = jj + 1
		

		mstars = np.where(mstars > 0., mstars, 0.1)		
		mdiffuse = np.where(mdiffuse > 0., mdiffuse, 0.1)		
		mWHIM = np.where(mWHIM > 0., mWHIM, 0.1)		
		mcoldgas = np.where(mcoldgas > 0., mcoldgas, 0.1)		
		mhotgas = np.where(mhotgas > 0., mhotgas, 0.1)		

		Zstars = mmstars/mstars/Zsun
		Zdiffuse = mmdiffuse/mdiffuse/Zsun
		ZWHIM = mmWHIM/mWHIM/Zsun
		Zhotgas = mmhotgas/mhotgas/Zsun
		Zcoldgas = mmcoldgas/mcoldgas/Zsun
		
		Zall = []
		Zall.append(np.log10(Zstars))
		Zall.append(np.log10(Zcoldgas))
		Zall.append(np.log10(Zhotgas))
		Zall.append(np.log10(Zdiffuse))
		Zall.append(np.log10(ZWHIM))


		fig = plt.figure(figsize=(12.,10))							# square box
#		fig.subplots_adjust(wspace = .0,top = .96, bottom = .1, left = .05, right = .95)
		ax = fig.add_subplot(1,1,1)
	
#		plt.axis([0,5,-3.5,0.2])
		ax.set_ylim(-3.5, 0.2)
		ax.set_xlim(0.02,5)
		ax.set_xscale('log')
		
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(30)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(30)


		for ir in range(nreheat) :
			for ie in range(neject) :
				label = 'R' + str(ir + 1) + 'E' + str(ie + 1)
				if ((ir == 0) & (ie == 0)) : label = label + ', Stars'
				ax.plot(zz[ir][ie], np.log10(Zstars[ir][ie]), c = stars_color[ir], ls = line_cycle[ie], lw = 4, label = label)

		for ir in range(nreheat) :
			for ie in range(neject) :
				label = 'R' + str(ir + 1) + 'E' + str(ie + 1)
				if ((ir == 0) & (ie == 0)) : label = label + ', ColdGas'
				ax.plot(zz[ir][ie], np.log10(Zcoldgas[ir][ie]), c = coldgas_color[ir], ls = line_cycle[ie], lw = 4, label = label)

		for ir in range(nreheat) :
			for ie in range(neject) :
				label = 'R' + str(ir + 1) + 'E' + str(ie + 1)
				if ((ir == 0) & (ie == 0)) : label = label + ', HaloGas'
				ax.plot(zz[ir][ie], np.log10(Zhotgas[ir][ie]), c = hotgas_color[ir], ls = line_cycle[ie], lw = 4, label = label)

		for ir in range(nreheat) :
			for ie in range(neject) :
				label = 'R' + str(ir + 1) + 'E' + str(ie + 1)
				if ((ir == 0) & (ie == 0)) : label = label + ', IGM'
				ax.plot(zz[ir][ie], np.log10(Zdiffuse[ir][ie]), c = igm_color[ir], ls = line_cycle[ie], lw = 4, label = label)

		
		for ir in range(nreheat) :
			for ie in range(neject) :
				label = 'R' + str(ir + 1) + 'E' + str(ie + 1)
				if ((ir == 0) & (ie == 0)) : label = label + ', WHIM'
				ax.plot(zz[ir][ie], np.log10(ZWHIM[ir][ie]), c = WHIM_color[ir], ls = line_cycle[ie], lw = 4, label = label)

		
		ax.set_xlabel(r'z', fontsize=34) #54 for square
		ax.set_ylabel(r'$log_{10}(Z/Z_{\odot})$', fontsize=34)  #36 for stripe
#		ax.legend(prop = matplotlib.font_manager.FontProperties(size=20),fancybox=True,loc=0, ncol = 3)
		plt.yticks(np.linspace(-3, 0, 4))
		ax.set_xticks([0.1, 0.5,1,5])
		ax.set_xticklabels([0.1, 0.5,1,5])

		outputFile = Dir + 'W' + '_' + 'Z_Zsun_vs_z' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()

#####
		ymin = [-0.59, -0.59, -1.09, -3.5,-3.5]
		ymax = [0, 0, 0, 0, 0]
		fig = plt.figure(figsize=(10.,20))							# square box
		fig.subplots_adjust(hspace = .05,top = .96, bottom = .05, left = .15, right = .95)

		for it in range(len(name_cycle)) :

			ax = fig.add_subplot(5,1,it + 1)	
			ax.set_ylim(ymin[it], ymax[it])
			ax.set_xlim(0.02,5)
			ax.set_xscale('log')
		
			for tick in ax.xaxis.get_major_ticks():
				tick.label1.set_fontsize(30)
				if (it != len(name_cycle) - 1) :		
					tick.label1On = False
			for tick in ax.yaxis.get_major_ticks():
				tick.label1.set_fontsize(30)

			for ir in range(nreheat) :
				for ie in range(neject) :
					label = 'R' + str(ir + 1) + 'E' + str(ie + 1)
	#				if ((ir == 0) & (ie == 0)) : label = label #+ ', ' + name_cycle[it]
					ax.plot(zz[ir][ie], Zall[it][ir][ie], c = colors[it][ir], ls = line_cycle[ie], lw = 4, label = label)

			ax.set_ylabel(r'$log_{10}(Z/Z_{\odot})$', fontsize=34)  #36 for stripe
			ax.legend(prop = matplotlib.font_manager.FontProperties(size=20),fancybox=True,loc=0, ncol = 3)
			plt.text(0.022, ymin[it]*0.03, name_cycle[it], size = 34, color = colors[it][0], verticalalignment='top', horizontalalignment='left')

			if (ymin[it] == -3.5) : 	plt.yticks(np.linspace(-3, 0, 4))


			if (it == len(name_cycle) - 1) :		
				ax.set_xlabel(r'z', fontsize=34) #54 for square
				ax.set_ylabel(r'$log_{10}(Z/Z_{\odot})$', fontsize=34)  #36 for stripe
				ax.legend(prop = matplotlib.font_manager.FontProperties(size=20),fancybox=True,loc=0, ncol = 3)
				ax.set_xticks([0.1, 0.5,1,5])
				ax.set_xticklabels([0.1, 0.5,1,5])

		outputFile = Dir + 'W' + '_' + 'Z_Zsun_vs_z_sep' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()


	def Z_frac_by_z(self, snapnum) :
	
		filename = SIM + '/data_files/' + 'metals_all_z_000'
		print 'reading in file:', filename
		
		zz = np.zeros(numsnaps)
		mstars = np.zeros(numsnaps)
		mmstars = np.zeros(numsnaps)
		mdiffuse = np.zeros(numsnaps)
		mmdiffuse = np.zeros(numsnaps)
		mhotgas = np.zeros(numsnaps)
		mmhotgas = np.zeros(numsnaps)
		mcoldgas = np.zeros(numsnaps)
		mmcoldgas = np.zeros(numsnaps)
		ii = 0
		for item in file(filename) :
			item = item.split()
			zz[ii] = float(item[0])
			mstars[ii] = float(item[1])
			mmstars[ii] = float(item[2])
			mdiffuse[ii] = float(item[9])
			mmdiffuse[ii] = float(item[4])
			mhotgas[ii] = float(item[5])
			mmhotgas[ii] = float(item[6])
			mcoldgas[ii] = float(item[7])
			mmcoldgas[ii] = float(item[8])
			ii = ii + 1
		
		
		Z_tot = mmstars + mmdiffuse + mmhotgas + mmcoldgas
		Z_tot = np.where(Z_tot > 0., Z_tot, 0.1)		

		Zstars = mmstars/Z_tot
		Zdiffuse = mmdiffuse/Z_tot
		Zhotgas = mmhotgas/Z_tot
		Zcoldgas = mmcoldgas/Z_tot


		print 'Z frac: stars, diffuse, hot gas, cold gas:'
		print Zstars[snapnum]
		print Zdiffuse[snapnum]
		print Zhotgas[snapnum]
		print Zcoldgas[snapnum]


		fig = plt.figure(figsize=(12.,10))							# square box
#		fig.subplots_adjust(wspace = .0,top = .96, bottom = .1, left = .05, right = .95)
		ax = fig.add_subplot(1,1,1)
	
		plt.axis([0,5,0.0,1])
#		ax.set_yscale('log')
		
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(30)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(30)

		ax.plot(zz, Zstars, c = '#00ff00', ls = '-', lw = 4, label = 'Stars')
		ax.plot(zz, Zcoldgas, c = 'blue', ls = '-', lw = 4, label = 'Cold Gas (ISM)')
		ax.plot(zz, Zhotgas, c = 'purple', ls = '-', lw = 4, label = 'Hot Gas (Halo)')
		ax.plot(zz, Zdiffuse, c = '#cc0066', ls = '-', lw = 4, label = 'Diffuse Gas (IGM&ICM)')

		ax.set_xlabel(r'z', fontsize=34) #54 for square
		ax.set_ylabel(r'$M_Z/M_{Z}^{total}$', fontsize=34)  #36 for stripe
		ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0)

# 		axR = fig.add_subplot(1,1,1, sharex=ax, frameon=False)
# 		plt.axis([0,5,0.0,1])
# 		axR.yaxis.tick_right()
# 		axR.yaxis.set_label_position("right")
# 		axR.set_ylabel(r'M/M$_{total}$')


		outputFile = Dir + SIM + '_' + 'Z_Ztot_vs_z' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()

	def Z_frac_by_z_by_model(self, nreheat, neject) :
	
		zz = np.zeros((nreheat, neject, numsnaps))
		mstars = np.zeros((nreheat, neject, numsnaps))
		mmstars = np.zeros((nreheat, neject, numsnaps))
		mdiffuse = np.zeros((nreheat, neject, numsnaps))
		mmdiffuse = np.zeros((nreheat, neject, numsnaps))
		mWHIM = np.zeros((nreheat, neject, numsnaps))
		mmWHIM = np.zeros((nreheat, neject, numsnaps))
		mhotgas = np.zeros((nreheat, neject, numsnaps))
		mmhotgas = np.zeros((nreheat, neject, numsnaps))
		mcoldgas = np.zeros((nreheat, neject, numsnaps))
		mmcoldgas = np.zeros((nreheat, neject, numsnaps))

		for ir in range(nreheat) :
			for ie in range(neject) :
				filename = 'W' + str(ir + 1) + str(ie + 1) + '/metals_all_z_000'
				print 'reading in file:', filename
		
				jj = 0
				for item in file(filename) :
					item = item.split()
					zz[ir][ie][jj] = float(item[0])
					mstars[ir][ie][jj] = float(item[1])
					mmstars[ir][ie][jj] = float(item[2])
					mdiffuse[ir][ie][jj] = float(item[3]) + float(item[9]) - float(item[10])
					mmdiffuse[ir][ie][jj] = float(item[4]) - float(item[11])
					mWHIM[ir][ie][jj] = float(item[10])
					mmWHIM[ir][ie][jj] = float(item[11])
					mhotgas[ir][ie][jj] = float(item[5])
					mmhotgas[ir][ie][jj] = float(item[6])
					mcoldgas[ir][ie][jj] = float(item[7])
					mmcoldgas[ir][ie][jj] = float(item[8])
					jj = jj + 1
		
		
		
		Z_tot = mmstars + mmdiffuse + mmhotgas + mmcoldgas + mmWHIM
		Z_tot = np.where(Z_tot > 0., Z_tot, 0.1)		

		Zstars = mmstars/Z_tot
		Zdiffuse = mmdiffuse/Z_tot
		ZWHIM = mmWHIM/Z_tot
		Zhotgas = mmhotgas/Z_tot
		Zcoldgas = mmcoldgas/Z_tot


		print 'Z frac: stars, diffuse, hot gas, cold gas:'
		print Zstars[1][1][snapnum]
		print Zdiffuse[1][1][snapnum]
		print Zhotgas[1][1][snapnum]
		print Zcoldgas[1][1][snapnum]

		Zall = []
		Zall.append(np.log10(Zstars))
		Zall.append(np.log10(Zcoldgas))
		Zall.append(np.log10(Zhotgas))
		Zall.append(np.log10(Zdiffuse))
		Zall.append(np.log10(ZWHIM))

		logz = np.log10(1. + zz)

		fig = plt.figure(figsize=(12.,10))							# square box
#		fig.subplots_adjust(wspace = .0,top = .96, bottom = .1, left = .05, right = .95)
		ax = fig.add_subplot(1,1,1)
	
		plt.axis([0.02,5,0.0,1])
		ax.set_xscale('log')
		
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(30)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(30)

		for ir in range(nreheat) :
			for ie in range(neject) :
				label = 'R' + str(ir + 1) + 'E' + str(ie + 1)
				if ((ir == 0) & (ie == 0)) : label = label + ', Stars'
				ax.plot(zz[ir][ie], np.log10(Zstars[ir][ie]), c = stars_color[ir], ls = line_cycle[ie], lw = 4, label = label)

		for ir in range(nreheat) :
			for ie in range(neject) :
				label = 'R' + str(ir + 1) + 'E' + str(ie + 1)
				if ((ir == 0) & (ie == 0)) : label = label + ', ColdGas'
				ax.plot(zz[ir][ie], np.log10(Zcoldgas[ir][ie]), c = coldgas_color[ir], ls = line_cycle[ie], lw = 4, label = label)

		for ir in range(nreheat) :
			for ie in range(neject) :
				label = 'R' + str(ir + 1) + 'E' + str(ie + 1)
				if ((ir == 0) & (ie == 0)) : label = label + ', HaloGas'
				ax.plot(zz[ir][ie], np.log10(Zhotgas[ir][ie]), c = hotgas_color[ir], ls = line_cycle[ie], lw = 4, label = label)

		for ir in range(nreheat) :
			for ie in range(neject) :
				label = 'R' + str(ir + 1) + 'E' + str(ie + 1)
				if ((ir == 0) & (ie == 0)) : label = label + ', IGM'
				ax.plot(zz[ir][ie], np.log10(Zdiffuse[ir][ie]), c = igm_color[ir], ls = line_cycle[ie], lw = 4, label = label)

		
		for ir in range(nreheat) :
			for ie in range(neject) :
				label = 'R' + str(ir + 1) + 'E' + str(ie + 1)
				if ((ir == 0) & (ie == 0)) : label = label + ', WHIM'
				ax.plot(zz[ir][ie], np.log10(ZWHIM[ir][ie]), c = WHIM_color[ir], ls = line_cycle[ie], lw = 4, label = label)

		
		ax.set_xlabel(r'z', fontsize=34) #54 for square
		ax.set_ylabel(r'$log_{10}(M_Z/M_{Z}^{total})$', fontsize=34)  #36 for stripe
#		ax.legend(prop = matplotlib.font_manager.FontProperties(size=20),fancybox=True,loc=0, ncol = 3)
		plt.yticks(np.linspace(-3, 0, 4))
		ax.set_xticks([0.1, 0.5,1,5])
		ax.set_xticklabels([0.1, 0.5,1,5])

		outputFile = Dir + 'W' + '_' + 'Z_Ztot_vs_z' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()

#####
#		ymin = [-0.59, -0.59, -1.09, -3.5,-3.5]
		ymin = [-2, -2, -2, -2, -2]
		ymax = [-0., -0.1, -0.1, -0.1, -0.1]
		fig = plt.figure(figsize=(10.,20))							# square box
		fig.subplots_adjust(hspace = .05,top = .96, bottom = .05, left = .15, right = .95)

		for it in range(len(name_cycle)) :

			ax = fig.add_subplot(5,1,it + 1)	
			ax.set_ylim(ymin[it], ymax[it])
			ax.set_xlim(0.02,5)
			ax.set_xscale('log')
		
			for tick in ax.xaxis.get_major_ticks():
				tick.label1.set_fontsize(30)
				if (it != len(name_cycle) - 1) :		
					tick.label1On = False
			for tick in ax.yaxis.get_major_ticks():
				tick.label1.set_fontsize(30)

			for ir in range(nreheat) :
				for ie in range(neject) :
					label = 'R' + str(ir + 1) + 'E' + str(ie + 1)
	#				if ((ir == 0) & (ie == 0)) : label = label #+ ', ' + name_cycle[it]
					ax.plot(zz[ir][ie], Zall[it][ir][ie], c = colors[it][ir], ls = line_cycle[ie], lw = 4, label = label)

			if (it == 2) : ax.set_ylabel(r'$log_{10}(M_Z/M_{Z}^{total})$', fontsize=34)  #36 for stripe
			ax.legend(prop = matplotlib.font_manager.FontProperties(size=20),fancybox=True,loc=0, ncol = 3)
			plt.text(4.5, ymin[it]*0.85, name_cycle[it], size = 34, color = colors[it][0], verticalalignment='top', horizontalalignment='right')

			if (ymin[it] == -3.5) : 	plt.yticks(np.linspace(-3, 0, 4))


			if (it == len(name_cycle) - 1) :		
				ax.set_xlabel(r'z', fontsize=34) #54 for square
#				ax.set_ylabel(r'$log_{10}(Z/Z_{\odot})$', fontsize=34)  #36 for stripe
				ax.legend(prop = matplotlib.font_manager.FontProperties(size=20),fancybox=True,loc=0, ncol = 3)
				ax.set_xticks([0.1, 0.5,1,5])
				ax.set_xticklabels([0.1, 0.5,1,5])

		outputFile = Dir + 'W' + '_' + 'Z_Ztot_vs_z_sep' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()


	def M_frac_by_z(self, snapnum) :
	
		filename = SIM + '/data_files/' + 'metals_all_z_000'
		print 'reading in file:', filename
		
		zz = np.zeros(numsnaps)
		mstars = np.zeros(numsnaps)
		mmstars = np.zeros(numsnaps)
		mdiffuse = np.zeros(numsnaps)
		mmdiffuse = np.zeros(numsnaps)
		mhotgas = np.zeros(numsnaps)
		mmhotgas = np.zeros(numsnaps)
		mcoldgas = np.zeros(numsnaps)
		mmcoldgas = np.zeros(numsnaps)
		ii = 0
		for item in file(filename) :
			item = item.split()
			zz[ii] = float(item[0])
			mstars[ii] = float(item[1])
			mmstars[ii] = float(item[2])
			mdiffuse[ii] = float(item[3]) 
			mmdiffuse[ii] = float(item[4])
			mhotgas[ii] = float(item[5])
			mmhotgas[ii] = float(item[6])
			mcoldgas[ii] = float(item[7])
			mmcoldgas[ii] = float(item[8])
			ii = ii + 1
		

		M_tot = mstars + mdiffuse + mhotgas + mcoldgas
		M_tot = np.where(M_tot > 0., M_tot, 0.1)
		Mstars = mstars/M_tot
		Mdiffuse = mdiffuse/M_tot
		Mhotgas = mhotgas/M_tot
		Mcoldgas = mcoldgas/M_tot


		print 'M frac: stars, diffuse, hot gas, cold gas:'
		print Mstars[snapnum]
		print Mdiffuse[snapnum]
		print Mhotgas[snapnum]
		print Mcoldgas[snapnum]

		fig = plt.figure(figsize=(12.,10))							# square box
#		fig.subplots_adjust(wspace = .0,top = .96, bottom = .1, left = .05, right = .95)
		ax = fig.add_subplot(1,1,1)
	
		plt.axis([0,5,0.0,1])
#		ax.set_yscale('log')
		
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(30)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(30)

		ax.plot(zz, Mstars, c = '#00ff00', ls = '-', lw = 4, label = 'Stars')
		ax.plot(zz, Mcoldgas, c = 'blue', ls = '-', lw = 4, label = 'Cold Gas (ISM)')
		ax.plot(zz, Mhotgas, c = 'purple', ls = '-', lw = 4, label = 'Hot Gas (Halo)')
		ax.plot(zz, Mdiffuse, c = '#cc0066', ls = '-', lw = 4, label = 'Diffuse Gas (IGM&ICM)')

		ax.set_xlabel(r'z', fontsize=34) #54 for square
		ax.set_ylabel(r'$M_b/M_{b}^{total}$', fontsize=34)  #36 for stripe
		ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0)

# 		axR = fig.add_subplot(1,1,1, sharex=ax, frameon=False)
# 		plt.axis([0,5,0.0,1])
# 		axR.yaxis.tick_right()
# 		axR.yaxis.set_label_position("right")
# 		axR.set_ylabel(r'M/M$_{total}$')


		outputFile = Dir + SIM + '_' + 'M_Mtot_vs_z' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()

	def M_frac_by_z_by_model(self, nmodels) :
	
		zz = np.zeros((nmodels, numsnaps))
		mstars = np.zeros((nmodels, numsnaps))
		mmstars = np.zeros((nmodels, numsnaps))
		mdiffuse = np.zeros((nmodels, numsnaps))
		mmdiffuse = np.zeros((nmodels, numsnaps))
		mhotgas = np.zeros((nmodels, numsnaps))
		mmhotgas = np.zeros((nmodels, numsnaps))
		mcoldgas = np.zeros((nmodels, numsnaps))
		mmcoldgas = np.zeros((nmodels, numsnaps))

		for ii in range(nmodels) :
			filename = 'W' + str(ii) + '/data_files/' + 'metals_all_z_000'
			print 'reading in file:', filename
		
			jj = 0
			for item in file(filename) :
				item = item.split()
				zz[ii][jj] = float(item[0])
				mstars[ii][jj] = float(item[1])
				mmstars[ii][jj] = float(item[2])
				mdiffuse[ii][jj] = float(item[9])
				mmdiffuse[ii][jj] = float(item[4])
				mhotgas[ii][jj] = float(item[5])
				mmhotgas[ii][jj] = float(item[6])
				mcoldgas[ii][jj] = float(item[7])
				mmcoldgas[ii][jj] = float(item[8])
				jj = jj + 1
		
		
		M_tot = mstars + mdiffuse + mhotgas + mcoldgas
		M_tot = np.where(M_tot > 0., M_tot, 0.1)
		Mstars = mstars/M_tot
		Mdiffuse = mdiffuse/M_tot
		Mhotgas = mhotgas/M_tot
		Mcoldgas = mcoldgas/M_tot


		print 'M frac: stars, diffuse, hot gas, cold gas:'
		print Mstars[2][snapnum]
		print Mdiffuse[2][snapnum]
		print Mhotgas[2][snapnum]
		print Mcoldgas[2][snapnum]

		fig = plt.figure(figsize=(12.,10))							# square box
#		fig.subplots_adjust(wspace = .0,top = .96, bottom = .1, left = .05, right = .95)
		ax = fig.add_subplot(1,1,1)
	
		plt.axis([0,5,0.0001,1.1])
		ax.set_yscale('log')
		
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(30)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(30)


		ii = 1
		ax.plot(zz[ii], Mstars[ii], c = '#00ff00', ls = '-', lw = 4, label = 'Stars, Constant Wind')
		ax.plot(zz[ii], Mcoldgas[ii], c = 'blue', ls = '-', lw = 4, label = 'Cold Gas (ISM)')
		ax.plot(zz[ii], Mhotgas[ii], c = 'purple', ls = '-', lw = 4, label = 'Hot Gas (Halo)')
		ax.plot(zz[ii], Mdiffuse[ii], c = '#cc0066', ls = '-', lw = 4, label = 'Diffuse Gas (IGM&ICM)')

		ii = 2
		ax.plot(zz[ii], Mstars[ii], c = '#00ff00', ls = '--', lw = 4)
		ax.plot(zz[ii], Mcoldgas[ii], c = 'blue', ls = '--', lw = 4)
		ax.plot(zz[ii], Mhotgas[ii], c = 'purple', ls = '--', lw = 4)
		ax.plot(zz[ii], Mdiffuse[ii], c = '#cc0066', ls = '--', lw = 4, label = 'Momentum Driven Wind')

		ii = 0
		ax.plot(zz[ii], Mstars[ii], c = '#00ff00', ls = ':', lw = 4)
		ax.plot(zz[ii], Mcoldgas[ii], c = 'blue', ls = ':', lw = 4)
		ax.plot(zz[ii], Mhotgas[ii], c = 'purple', ls = ':', lw = 4)
		ax.plot(zz[ii], Mdiffuse[ii], c = '#cc0066', ls = ':', lw = 4, label = 'No Wind')

		ax.set_xlabel(r'z', fontsize=34) #54 for square
		ax.set_ylabel(r'$M_b/M_{b}^{total}$', fontsize=34)  #36 for stripe
		ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0)

# 		axR = fig.add_subplot(1,1,1, sharex=ax, frameon=False)
# 		plt.axis([0,5,0.0,1])
# 		axR.yaxis.tick_right()
# 		axR.yaxis.set_label_position("right")
# 		axR.set_ylabel(r'M/M$_{total}$')


		outputFile = Dir + 'W' + '_' + 'M_Mtot_vs_z' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()


	def Z_Zsun_by_rho(self, snapnum) :
	
		snap = "%03d" % snapnum
		Zfilename = SIM + '/data_files/' + SIM + '_' + 'metals_' + str(nbins) + '_' + snap + '.cpickle'
		Hfilename = SIM + '/data_files/' + SIM + '_' + 'H_' + str(nbins) + '_' + snap + '.cpickle'
		dfilename = SIM + '/data_files/' + SIM + '_' + 'delta_' + str(nbins) + '_' + snap + '.cpickle'
		
		f = open(Zfilename)
		metals = cPickle.load(f)
		f.close()
		metals = np.array(metals)
		w = np.where(metals > 0.0)[0]
		metals = metals[w]

		f = open(Hfilename)
		H = cPickle.load(f)
		f.close()
		H = np.array(H)
		H = H[w]

		f = open(dfilename)
		delta1 = cPickle.load(f)
		f.close()
		delta1 = np.array(delta1)
		delta1 = delta1[w]

		delta1 = np.where(delta1 > 0., delta1, 0.01)
		logrho = np.log10(delta1)
		
		
		Zdiffuse = metals/H/Zsun
		Zdiffuse = np.where(Zdiffuse > 0., Zdiffuse, 1.0e-8)
		
		logZ = np.log10(Zdiffuse)
#		logZMW = np.log10(Zdiffuse*H)


		print 'Z/Zsun min-max', np.min(Zdiffuse), np.max(Zdiffuse)
		print 'rho min-max', np.min(logrho), np.max(logrho)
		rhobins = np.linspace(np.min(logrho)+2., np.max(logrho), 21)
		rhox = np.linspace(np.min(logrho)+2., np.max(logrho), 20)
		zrho = []
		zrho_err = []
		for ii in range(len(rhobins) - 1) :
			w = np.where((logrho > rhobins[ii]) & (logrho < rhobins[ii + 1]))[0]
			print len(w)
			if len(w) > 0 :
				ave = np.average(Zdiffuse[w])
				stdev = np.std(Zdiffuse[w])
			else :
				ave = -0.
				stdev = 0.
			zrho.append(ave)
			zrho_err.append(stdev)

		zrho = np.log10(zrho)
		zrho_err = np.log10(zrho_err)
		rho_err = (rhox[1] - rhox[0])/2.
		
		print rhobins
		print zrho
		print rho_err

		fig = plt.figure(figsize=(12.,10))							# square box
#		fig.subplots_adjust(wspace = .0,top = .96, bottom = .1, left = .05, right = .95)
		ax = fig.add_subplot(1,1,1)
	
		plt.axis([-2.5,8,-8,2])
#		ax.set_yscale('log')
		
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(30)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(30)

		ax.scatter(logrho, logZ, s = 2, alpha = .1, c = '#cc0066')
		ax.scatter(rhox[10], zrho[10], label = 'Diffuse Gas (IGM&ICM)', s = 2, alpha = 1, c = '#cc0066',)
		ax.errorbar(rhox, zrho, yerr = zrho_err, xerr = rho_err, ls = 'None')
		ax.set_xlabel(r'$log_{10}(\rho/\bar{\rho})$', fontsize=34) #54 for square
		ax.set_ylabel(r'$log_{10}(Z/Z_{\odot})$', fontsize=34)  #36 for stripe
		ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0)

		outputFile = Dir + SIM + '_' + 'Z_Zsun_vs_rho' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()

	def Z_Zsun_by_rho_z(self, snapnum0, snapnum1, snapnum2, snapnum3) :

		zs = [snapnum0, snapnum1, snapnum2, snapnum3]

		metals = []
		H = []
		logrho = []
		logZ = []
		Zdiffuse = []
		
		zsnap = []
		snapfile = SIM + '/data_files/' + SIM + 'z.txt'
		for item in file(snapfile) :
			item = item.split()
			zsnap.append(float(item[0]))
	
		for zz in zs :
			snap = "%03d" % zz
			Zfilename = SIM + '/data_files/' + SIM + '_' + 'metals_' + str(nbins) + '_' + snap + '.cpickle'
			Hfilename = SIM + '/data_files/' + SIM + '_' + 'H_' + str(nbins) + '_' + snap + '.cpickle'
			dfilename = SIM + '/data_files/' + SIM + '_' + 'delta_' + str(nbins) + '_' + snap + '.cpickle'
		
			f = open(Zfilename)
			ms = cPickle.load(f)
			f.close()
			ms = np.array(ms)
			w = np.where(ms > 0.0)[0]
			metals.append(ms[w])

			f = open(Hfilename)
			H1 = cPickle.load(f)
			f.close()
			H1 = np.array(H1)
			H.append(H1[w])

			f = open(dfilename)
			delta1 = cPickle.load(f)
			f.close()
			delta1 = np.array(delta1)
			delta1 = np.where(delta1 > 0., delta1, 0.01)
			delta1 = np.log10(delta1)
			logrho.append(delta1[w])
		
			Zdiff = ms/H1/Zsun
			Zdiffuse.append(np.where(Zdiff > 0., Zdiff, 1.0e-8))
		
			logZ.append(np.log10(Zdiffuse))
#		logZMW = np.log10(Zdiffuse*H)


		print 'Z/Zsun min-max', np.min(Zdiffuse), np.max(Zdiffuse)
		print 'rho min-max', np.min(logrho[0]), np.max(logrho[0])
		rhobins = np.linspace(-2.5, 5., 21)
		off = (rhobins[1] - rhobins[0])/2.
		rhox = np.linspace(-2.5 + off, 5 - off, 20)
		zrho = []
		zrho_err = []
		for zz in range(len(zs)) :
			aves = []
			errs = []
			for ii in range(len(rhobins) - 1) :
				w = np.where((logrho[zz] > rhobins[ii]) & (logrho[zz] < rhobins[ii + 1]))[0]
				print len(w)
				if len(w) > 0 :
					ave = np.average(Zdiffuse[zz][w])
					stdev = np.std(Zdiffuse[zz][w])
				else :
					ave = -0.
					stdev = 0.
				aves.append(ave)
				errs.append(stdev)
			zrho.append(aves)
			zrho_err.append(errs)

		zrho = np.log10(zrho)
		zrho_err = np.log10(zrho_err)
		rho_err = (rhox[1] - rhox[0])/2.
		
		print rhobins
		print zrho
		print rho_err

		fig = plt.figure(figsize=(12.,10))							# square box
#		fig.subplots_adjust(wspace = .0,top = .96, bottom = .1, left = .05, right = .95)
		ax = fig.add_subplot(1,1,1)
	
		plt.axis([-2.5,8,-8,2])
#		ax.set_yscale('log')
		
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(30)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(30)

		for zz in range(len(zs)) :
			rdsh = "%0.3f" % zsnap[zs[zz]]
			ax.plot(rhox, zrho[zz], ls = '-', label = 'z = '+ rdsh)
#			ax.errorbar(rhox, zrho[zz], yerr = zrho_err[zz], xerr = rho_err, ls = 'None', label = 'z = '+ rdsh)
		ax.set_xlabel(r'$log_{10}(\rho/\bar{\rho})$', fontsize=34) #54 for square
		ax.set_ylabel(r'$log_{10}(Z/Z_{\odot})$', fontsize=34)  #36 for stripe
		ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0)

		outputFile = Dir + SIM + '_' + 'Z_Zsun_vs_rho_vs_z' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()

	def f_b_by_z(self, snapnum) :
	
		snap = "%03d" % snapnum
		filename = SIM + '/data_files/' + 'bound_all_z_000'
		
		zz = np.zeros(numsnaps)
		mvir = np.zeros(numsnaps)
		bound_DM = np.zeros(numsnaps)
		bound_b = np.zeros(numsnaps)
		ii = 0
		for item in file(filename) :
			item = item.split()
			zz[ii] = float(item[0])
			bound_DM[ii] = float(item[1])
			mvir[ii] = float(item[2])
			bound_b[ii] = float(item[3])
			ii = ii + 1

		filename = SIM + '/' + 'metals_all_z_000'
		mstars = np.zeros(numsnaps)
		mdiffuse = np.zeros(numsnaps)
		mhotgas = np.zeros(numsnaps)
		mcoldgas = np.zeros(numsnaps)
		ii = 0
		for item in file(filename) :
			item = item.split()
			mstars[ii] = float(item[1])
			mdiffuse[ii] = float(item[3])
			mhotgas[ii] = float(item[5])
			mcoldgas[ii] = float(item[7])
			ii = ii + 1
		
		
		fb_bound = bound_b/bound_DM
		fb_mvir = bound_b/mvir
		fb_all = (mstars + mcoldgas + mhotgas)/mvir
		f_mvir = mvir/bound_DM
		f_stars = mstars/mvir
		f_hot = mhotgas/mvir
		f_cold = mcoldgas/mvir


		fig = plt.figure(figsize=(12.,10))							# square box
#		fig.subplots_adjust(wspace = .0,top = .96, bottom = .1, left = .05, right = .95)
		ax = fig.add_subplot(1,1,1)
	
		plt.axis([0,5,0.001,1])
		ax.set_yscale('log')
		
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(30)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(30)

#		ax.plot(zz, fb_all, ls = '--', lw = 4, label = '$M_{Stars + Cold + Hot}/M_{vir}$')
		ax.plot(zz, f_stars, ls = '--', lw = 4, c = '#00FF00', label = '$M_{Stars}/M_{vir}$')
		ax.plot(zz, f_cold, ls = '--', lw = 4, c = 'blue', label = '$M_{Cold}/M_{vir}$')
		ax.plot(zz, f_hot, ls = '--', lw = 4, c = 'purple', label = '$M_{Hot}/M_{vir}$')
		ax.plot(zz, fb_bound, ls = '-', lw = 4, c = '0.5', label = '$M_b/M_{DM}$')
		ax.plot(zz, fb_mvir, ls = '-', lw = 4, c = 'k', label = '$M_b/M_{vir}$')
		ax.axhline(y=0.16, ls = ':', lw = 3)


		ax.set_xlabel(r'z', fontsize=34) #54 for square
		ax.set_ylabel(r'$f_b$', fontsize=34)  #36 for stripe
		ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0, ncol = 2)

		outputFile = Dir + SIM + '_' + 'f_b_by_z' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()

# =================================================================


#  'Main' section of code.  This if statement executes if the code is run from the 
#   shell command line, i.e. with 'python allresults.py'

if __name__ == '__main__':
	res = Results()

 	snapnum = 63




#	res.Z_frac_by_z_by_model(nreheat, neject)
## 	res.M_frac_by_z_by_model(nmodels)
#	res.Z_Zsun_by_z_by_model(nreheat, neject)
#	res.fM_Z_halos(snapnum, 1)
#	res.M_Z_halos(snapnum, 1)
# 	res.f_b_halos(snapnum, 'True', 'False', 'False')
#	res.f_b_by_z(snapnum)
# 	res.Z_Zsun_by_z(snapnum)
# 	res.Z_frac_by_z(snapnum)
# 	res.M_frac_by_z(snapnum)
#	res.Z_Zsun_by_rho(snapnum)
#	res.Z_hist(snapnum)
#	res.Z_Zsun_by_rho_z(27, 32, 40, 63)
#	res.rho_hist(snapnum, nreheat, neject)
# 	res.T_hist(snapnum)
#	res.T_vs_rho(snapnum, nreheat, neject)
# 	res.M_vs_rho_vs_reservoir(snapnum, nreheat, neject)
# 	res.M_vs_rho_vs_reservoir_resid(snapnum, nres, nreheat, neject)
# 	res.M_Z_vs_rho_vs_reservoir_resid(snapnum, nres, nreheat, neject)
# 	res.Z_vs_rho_vs_reservoir_resid(snapnum, nres, nreheat, neject)	