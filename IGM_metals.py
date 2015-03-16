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
from os.path import getsize as getFileSize


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
Hubble_h = 0.73


np.seterr(divide='ignore')

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
IGMgas_color = ['#000099', '#0033ff', '#66ccff']
ZZ_color = ['#330066', '#9933cc', '#cc99ff'] # '#660099',
Z_color = ['#660033', '#cc0099', '#ff99ff'] #990066
WHIM_color = ['#990000', '#cc3300', '#ff6633']
DM_color = ['#000000', '#666666', '#cccccc']

line_cycle = ['-', '--', ':']

colors = [stars_color, IGMgas_color, ZZ_color, Z_color, WHIM_color, DM_color]
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
OutputFormat = '.eps'
TRANSPARENT = False


class Results:

	def read_gals(self, model_name, first_file, last_file):
	
		# The input galaxy structure:
		Galdesc_full = [
			('Type'                         , np.int32),                    
			('GalaxyIndex'                  , np.int64),                    
			('HaloIndex'                    , np.int32),                    
			('FOFHaloIdx'                   , np.int32),                    
			('TreeIdx'                      , np.int32),                    
			('SnapNum'                      , np.int32),                    
			('CentralGal'                   , np.int32),                    
			('CentralMvir'                  , np.float32),                  
			('mergeType'                    , np.int32),                    
			('mergeIntoID'                  , np.int32),                    
			('mergeIntoSnapNum'             , np.int32),                    
			('dT'                           , np.int32),                    
			('Pos'                          , (np.float32, 3)),             
			('Vel'                          , (np.float32, 3)),             
			('Spin'                         , (np.float32, 3)),             
			('Len'                          , np.int32),                    
			('Mvir'                         , np.float32),                  
			('Rvir'                         , np.float32),                  
			('Vvir'                         , np.float32),                  
			('Vmax'                         , np.float32),                  
			('VelDisp'                      , np.float32),                  
			('ColdGas'                      , np.float32),                  
			('StellarMass'                  , np.float32),                  
			('BulgeMass'                    , np.float32),                  
			('HotGas'                       , np.float32),                  
			('EjectedMass'                  , np.float32),                  
			('BlackHoleMass'                , np.float32),                  
			('IntraClusterStars'            , np.float32),                  
			('MetalsColdGas'                , np.float32),                  
			('MetalsStellarMass'            , np.float32),                  
			('MetalsBulgeMass'              , np.float32),                  
			('MetalsHotGas'                 , np.float32),                  
			('MetalsEjectedMass'            , np.float32),                  
			('MetalsIntraClusterStars'      , np.float32),                  
			('SfrDisk'                      , np.float32),                  
			('SfrBulge'                     , np.float32),                  
			('SfrDiskZ'                     , np.float32),                  
			('SfrBulgeZ'                    , np.float32),                  
			('DiskRadius'                   , np.float32),                  
			('Cooling'                      , np.float32),                  
			('Heating'                      , np.float32),
			('LastMajorMerger'              , np.float32),
			('OutflowRate'                  , np.float32),
			('infallMvir'                   , np.float32),
			('infallVvir'                   , np.float32),
			('infallVmax'                   , np.float32)
			]
		names = [Galdesc_full[i][0] for i in xrange(len(Galdesc_full))]
		formats = [Galdesc_full[i][1] for i in xrange(len(Galdesc_full))]
		Galdesc = np.dtype({'names':names, 'formats':formats}, align=True)


		# Initialize variables.
		TotNTrees = 0
		TotNGals = 0
		FileIndexRanges = []

		print "Determining array storage requirements."
	
		# Read each file and determine the total number of galaxies to be read in
		goodfiles = 0
		for fnr in xrange(first_file,last_file+1):
			fname = model_name+'_'+str(fnr)  # Complete filename
	
			if getFileSize(fname) == 0:
				print "File\t%s  \tis empty!  Skipping..." % (fname)
				continue
	
			fin = open(fname, 'rb')  # Open the file
			Ntrees = np.fromfile(fin,np.dtype(np.int32),1)  # Read number of trees in file
			NtotGals = np.fromfile(fin,np.dtype(np.int32),1)[0]  # Read number of gals in file.
			TotNTrees = TotNTrees + Ntrees  # Update total sim trees number
			TotNGals = TotNGals + NtotGals  # Update total sim gals number
			goodfiles = goodfiles + 1  # Update number of files read for volume calculation
			fin.close()

		print
		print "Input files contain:\t%d trees ;\t%d galaxies ." % (TotNTrees, TotNGals)
		print

		# Initialize the storage array
		G = np.empty(TotNGals, dtype=Galdesc)

		offset = 0  # Offset index for storage array

		# Open each file in turn and read in the preamble variables and structure.
		print "Reading in files."
		for fnr in xrange(first_file,last_file+1):
			fname = model_name+'_'+str(fnr)  # Complete filename

			if getFileSize(fname) == 0:
				continue
	
			fin = open(fname, 'rb')  # Open the file
			Ntrees = np.fromfile(fin, np.dtype(np.int32), 1)  # Read number of trees in file
			NtotGals = np.fromfile(fin, np.dtype(np.int32), 1)[0]  # Read number of gals in file.
			GalsPerTree = np.fromfile(fin, np.dtype((np.int32, Ntrees)),1) # Read the number of gals in each tree
			print ":   Reading N=", NtotGals, "   \tgalaxies from file: ", fname
			GG = np.fromfile(fin, Galdesc, NtotGals)  # Read in the galaxy structures
	
			FileIndexRanges.append((offset,offset+NtotGals))
	
			# Slice the file array into the global array
			# N.B. the copy() part is required otherwise we simply point to
			# the GG data which changes from file to file
			# NOTE THE WAY PYTHON WORKS WITH THESE INDICES!
			G[offset:offset+NtotGals]=GG[0:NtotGals].copy()
		
			del(GG)
			offset = offset + NtotGals  # Update the offset position for the global array
	
			fin.close()  # Close the file


		print
		print "Total galaxies considered:", TotNGals

		# Convert the Galaxy array into a recarray
		G = G.view(np.recarray)

		w = np.where(G.StellarMass > 0.001)[0]
		print "Galaxies more massive than 10^7Msun/h:", len(w)

		print 'minimum mass is', np.min(G.StellarMass[w])
		
		print

		return G

	def f_b_halos(self, snapnum, p1, p2, p3) :
	
		plot1 = p1
		plot2 = p2
		plot3 = p3
		snap = "%03d" % snapnum
		filename = SIM + '/data_files/' + 'halos_0.000.dat'
		
		nhalos = 0
		for item in file(filename) :
			nhalos = nhalos + 1
		print nhalos, 'halos in file'
		
		type = np.zeros(nhalos)
		mvir = np.zeros(nhalos)
		mstars = np.zeros(nhalos)
		mcoldgas = np.zeros(nhalos)
		mhotgas = np.zeros(nhalos)
		other = np.zeros(nhalos)
		ii = 0
		for item in file(filename) :
			item = item.split()
			type[ii] = int(item[0])
			mvir[ii] = float(item[1])
			mstars[ii] = float(item[3])
			mcoldgas[ii] = float(item[9])
			mhotgas[ii] = float(item[7])
			other[ii] = float(item[11])
			ii = ii + 1

		baryons = mstars + mcoldgas + mhotgas + other
		fb_all = baryons/mvir
		fb_stars = mstars/mvir
		fb_cold = mcoldgas/mvir
		fb_hot = mhotgas/mvir

		f_stars = mstars/baryons
		f_cold = mcoldgas/baryons
		f_hot = mhotgas/baryons

		print np.sum(other), np.sum(baryons)
		print min(fb_all), max(fb_all)
		print min(fb_stars), max(fb_stars)
		print min(fb_cold), max(fb_cold)
		print min(fb_hot), max(fb_hot)
		
		w = np.where(fb_all > 1.)[0]
		for ii in range(len(w)) :
			print type[w][ii], "%.4f \t %.4f \t %e \t %e \t %e \t %e" % (fb_all[w][ii], fb_cold[w][ii], mvir[w][ii], mstars[w][ii], mcoldgas[w][ii], mhotgas[w][ii])
		
		
		bins = 200
		fhist, edges = np.histogram(fb_all, bins = bins)
		w = np.where(type == 0)[0]
		print len(w)
#		fb_cent = fb_all[w]
		fhist_0, edges = np.histogram(fb_all[w], bins = edges)
		print edges[15], np.sum(fhist_0[15::])
		print min(fb_all[w]), max(fb_all[w]), np.average(fb_all[w])
		w = np.where(type == 1)[0]
		print len(w)
#		fb_sat = fb_all[w]
		fhist_1, edges = np.histogram(fb_all[w], bins = edges)
		print min(fb_all[w]), max(fb_all[w]), np.average(fb_all[w])
		print edges[15], np.sum(fhist_1[15::])
		w = np.where((type == 1) & (fb_all <= 0.5))
		print np.average(fb_all[w])
		w = np.where((type == 0) & (mvir > 10**11.) & (mvir < 10**12.))[0]
		print len(w), 'mvir 11-12: ', np.average(fb_all[w])
		fhist_11, edges = np.histogram(fb_all[w], bins = edges)
		w = np.where((type == 0) & (mvir > 10**12.) & (mvir < 10**13.))[0]
		print len(w), 'mvir 12-13: ', np.average(fb_all[w])
		fhist_12, edges = np.histogram(fb_all[w], bins = edges)
		w = np.where((type == 0) & (mvir > 10**13.))[0]
		print len(w), 'mvir > 13: ', np.average(fb_all[w])
		fhist_13, edges = np.histogram(fb_all[w], bins = edges)


#		print fhist
#		print edges
		off = abs(edges[0] - edges[1])/2.
		print off
		f = np.linspace(edges[0] + off, edges[-1] - off, bins)
#		print f

		if plot1 == 'True' :
			fig = plt.figure(figsize=(12.,10))							# square box
	#		fig.subplots_adjust(wspace = .0,top = .96, bottom = .1, left = .05, right = .95)
			ax = fig.add_subplot(1,1,1)
	
			plt.axis([0,1,1,1.0e3])
			ax.set_yscale('log')
		
			for tick in ax.xaxis.get_major_ticks():
				tick.label1.set_fontsize(30)
			for tick in ax.yaxis.get_major_ticks():
				tick.label1.set_fontsize(30)

			ax.plot(f, fhist, label = 'All Galaxies')#, label = 'Diffuse Gas (IGM&ICM)', s = 2, alpha = .1)
			ax.plot(f, fhist_0, label = 'Central Galaxies', c = '#ff0066')#, label = 'Diffuse Gas (IGM&ICM)', s = 2, alpha = .1)
			ax.plot(f, fhist_1, label = 'Satellite Galaxies', c='#3399cc')#, label = 'Diffuse Gas (IGM&ICM)', s = 2, alpha = .1)
			ax.plot(f, fhist_11, label = r'$10^{11} < M_{vir} < 10^{12}M_{\odot}$', c='#6600cc', ls = ':')
			ax.plot(f, fhist_12, label = r'$10^{12} < M_{vir} < 10^{13}M_{\odot}$', c='#6600cc', ls = '--')
			ax.plot(f, fhist_13, label = r'$10^{13}M_{\odot} < M_{vir}$', c='#6600cc', ls = '-')
			ax.axvline(x = 0.17, ls = ':', lw = 4, c = 'k')

			ax.set_xlabel(r'$f_b$', fontsize=34) #54 for square
			ax.set_ylabel(r'$Number$', fontsize=34)  #36 for stripe
			ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0)

			outputFile = Dir + SIM + '_fb_hist' + OutputFormat
			plt.savefig(outputFile)  # Save the figure
			print 'Saved file to', outputFile
			plt.close()

		if plot2 == 'True' :
			fig = plt.figure(figsize=(12.,10))							# square box
	#		fig.subplots_adjust(wspace = .0,top = .96, bottom = .1, left = .05, right = .95)
			ax = fig.add_subplot(1,1,1)
	
			plt.axis([10,15,0,0.5])
	#		ax.set_xscale('log')
		
			for tick in ax.xaxis.get_major_ticks():
				tick.label1.set_fontsize(30)
			for tick in ax.yaxis.get_major_ticks():
				tick.label1.set_fontsize(30)

	#		ax.scatter(mvir, fb_all, label = 'All Galaxies')#, label = 'Diffuse Gas (IGM&ICM)', s = 2, alpha = .1)
			w = np.where(type == 0)[0]
			ax.scatter(np.log10(mvir[w]), fb_all[w], label = 'Central Galaxies', c = '#ff0066', edgecolor =  '#ff0066', alpha = 0.3)
			w = np.where(type == 1)[0]
			ax.scatter(np.log10(mvir[w]), fb_all[w], label = 'Satellite Galaxies', c='#3399cc', edgecolor = '#3399cc', alpha = 0.3)
			ax.axhline(y = 0.17, ls = ':', lw = 4, c = 'k')

			ax.set_xlabel(r'$log_{10}(M_{vir}/M_{\odot})$', fontsize=34) #54 for square
			ax.set_ylabel(r'$f_b$', fontsize=34)  #36 for stripe
			ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0)

			outputFile = Dir + SIM + '_fb_vs_Mvir' + OutputFormat
			plt.savefig(outputFile)  # Save the figure
			print 'Saved file to', outputFile
			plt.close()

		if plot3 == 'True' :
		
			fig = plt.figure(figsize=(12.,10))						
			fig.subplots_adjust(wspace = .0, hspace = 0.)
			ax = fig.add_subplot(2,2,1)
			plt.axis([10.1,15,0,1])
		
			for tick in ax.xaxis.get_major_ticks():
				tick.label1.set_fontsize(30)
				tick.label1On = False
			for tick in ax.yaxis.get_major_ticks():
				tick.label1.set_fontsize(30)

			w = np.where(type == 0)[0]
			ax.scatter(np.log10(mvir[w]), f_stars[w], label = 'Central Galaxies', c = '#ff0066', edgecolor =  '#ff0066', alpha = 0.3)
			w = np.where(type == 1)[0]
			ax.scatter(np.log10(mvir[w]), f_stars[w], label = 'Satellite Galaxies', c='#3399cc', edgecolor = '#3399cc', alpha = 0.5)

#			ax.set_xlabel(r'$log_{10}(M_{vir}/M_{\odot})$', fontsize=34) #54 for square
			ax.set_ylabel(r'$log_{10}(M_{\star}/M_{b})$', fontsize=34) #54 for square

			ax = fig.add_subplot(2,2,2)
			plt.axis([10.1,15,0,1])
		
			for tick in ax.xaxis.get_major_ticks():
				tick.label1.set_fontsize(30)
			for tick in ax.yaxis.get_major_ticks():
				tick.label1.set_fontsize(30)
				tick.label1On = False

			w = np.where(type == 0)[0]
			ax.scatter(np.log10(mvir[w]), f_cold[w], label = 'Central Galaxies', c = '#ff0066', edgecolor =  '#ff0066', alpha = 0.3)
			w = np.where(type == 1)[0]
			ax.scatter(np.log10(mvir[w]), f_cold[w], label = 'Satellite Galaxies', c='#3399cc', edgecolor = '#3399cc', alpha = 0.5)

			ax.yaxis.tick_right()
			ax.yaxis.set_label_position("right")

			ax.set_xlabel(r'$log_{10}(M_{vir}/M_{\odot})$', fontsize=34) #54 for square
			ax.set_ylabel(r'$log_{10}(M_{ColdGas}/M_{b})$', fontsize=34) #54 for square
			ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0)

			ax = fig.add_subplot(2,2,3)
			plt.axis([10.1,15,0,.999])
		
			for tick in ax.xaxis.get_major_ticks():
				tick.label1.set_fontsize(30)
			for tick in ax.yaxis.get_major_ticks():
				tick.label1.set_fontsize(30)

			w = np.where(type == 0)[0]
			ax.scatter(np.log10(mvir[w]), f_hot[w], label = 'Central Galaxies', c = '#ff0066', edgecolor =  '#ff0066', alpha = 0.3)
			w = np.where(type == 1)[0]
			ax.scatter(np.log10(mvir[w]), f_hot[w], label = 'Satellite Galaxies', c='#3399cc', edgecolor = '#3399cc', alpha = 0.5)

			ax.set_xlabel(r'$log_{10}(M_{vir}/M_{\odot})$', fontsize=34) #54 for square
			ax.set_ylabel(r'$log_{10}(M_{HotGas}/M_{b})$', fontsize=34) #54 for square

			outputFile = Dir + SIM + '_fb_frac_vs_Mvir' + OutputFormat
			plt.savefig(outputFile)  # Save the figure
			print 'Saved file to', outputFile
			plt.close()

	def rho_vs_Mvir_hist(self, snapnum, mass, ir, ie) :
	
		if mass == 'mvir' :
			mindex = 1
			Mname = 'Mvir'
			Mlabel = r'$\mathrm{M}_{vir}$'
			cmap = cm.BuPu
			cmapE = cm.BuPu
			cmapZ = cm.BuPu
			cmapZZ = cm.BuPu
			vlim = [0, 3.5]
			vlimI = [8.7, 11.2]
			vlimE = [7, 13]
			vlimEH = [6.5, 10.2]
			vlimZ = [4.5, 10.5]
			vlimZZ = [-4.5, -0.75]
			vlimZH = [3.8, 8.2]
			cmapI = cm.RdPu #cm.autumn_r

		if mass == 'mstars' :
			mindex = 3
			Mname = 'Mstars'
			Mlabel = r'M$_{\star}$'
			cmap = cm.YlGn #YlOrRd #autumn_r
			cmapE = cm.YlOrRd #cm.autumn_r
			cmapZ = cm.Greens# cm.autumn_r
			cmapZZ = cm.Purples# cm.autumn_r
			cmapZhZ = cm.Purples# cm.autumn_r
			cmapI = cm.Blues #cm.autumn_r
			vlim = [0, 3.]
			vlimE = [7.5, 12.5]
			vlimI = [8.7, 11.2]
			vlimEH = [7.0, 11.0]
			vlimZ = [4.5, 10.5]
			vlimZZ = [-4.25, -0.75]
			vlimZhZ = [-3.25, 0.25]
			vlimZH = [4., 8.5]

	
		snap = "%03d" % snapnum
		dir = Dir + mass + '_rho_' + snap + '/'
		sim = 'W' +  str(ir + 1) + str(ie + 1)
		mod = 'R' +  str(ir + 1) + 'E' + str(ie + 1)

		
		bound = np.zeros((nbins, nbins, nbins))
		unbound = np.zeros((nbins, nbins, nbins))
		diffuse = np.zeros((nbins, nbins, nbins))
		metals_grid = np.zeros((nbins, nbins, nbins))
		
		gridfile = 'MM/diffuse_files/diffuse_' + snap
		for item in file(gridfile) :
			item = item.split()
			xbin = int(item[0])
			ybin = int(item[1])
			zbin = int(item[2])
			unbound[xbin][ybin][zbin] = float(int(item[3]))
			bound[xbin][ybin][zbin] = float(int(item[4]))

		rho_grid = (unbound + bound)/(npcls)*float(nbins)**3

		diffusefile = 'millennium_' + sim + '/metals_files/metals_' + snap
		print 'reading in', diffusefile
		for item in file(diffusefile) :
			item = item.split()
			xbin = int(item[0])
			ybin = int(item[1])
			zbin = int(item[2])
			diffuse[xbin][ybin][zbin] = float(item[9])
			metals_grid[xbin][ybin][zbin] = float(item[10])
			
		diffuse = np.where(diffuse >= 0., diffuse, 0.)
		print diffuse[10][20][42]
		print metals_grid[10][20][42]
		print np.float(len(np.where(metals_grid > 0.)[0]))/100.**3 * 100., '% of cells have metals'

		w = np.where((G.Type == 0) & ( G.StellarMass > 0.001))[0]
		print len(w), 'central galaxies have > 10^7 stars'
		xbin = (G.Pos[w,0]/box_side_h1*nbins).astype(int)
		ybin = (G.Pos[w,1]/box_side_h1*nbins).astype(int)
		zbin = (G.Pos[w,2]/box_side_h1*nbins).astype(int)
		Mass = np.log10(G.StellarMass[w] * 1.0e10)
		ejected = G.EjectedMass[w] * 1.0e10
		zejected = G.MetalsEjectedMass[w] * 1.0e10
		halo = G.HotGas[w] * 1.0e10
		zhalo = G.MetalsHotGas[w] * 1.0e10
		mvir = G.Mvir[w] * 1.0e10
		stars = G.StellarMass[w] * 1.0e10
		zstars = G.MetalsStellarMass[w] * 1.0e10
		all_metals = G.MetalsStellarMass[w] + G.MetalsColdGas[w] + G.MetalsHotGas[w] + G.MetalsEjectedMass[w] + G.MetalsIntraClusterStars[w]
		all_metals = all_metals * 1.0e10
		coldgas = G.ColdGas[w] * 1.0e10
		vmax = G.Vmax[w]
		vvir = G.Vvir[w]
		rdisk = G.DiskRadius[w]
		sfr = G.SfrDisk[w]
		
		print np.sum(zejected)/np.sum(all_metals)*100., '% of metals are in the IGM'

		xyz = np.zeros((100, 100, 100))
		rho = np.zeros(len(w))
		IGM = np.zeros(len(w))
		zgrid = np.zeros(len(w))
		for gg in range(len(w)) :
			rho[gg] = np.log10(rho_grid[xbin[gg]][ybin[gg]][zbin[gg]])
			IGM[gg] = diffuse[xbin[gg]][ybin[gg]][zbin[gg]]
			zgrid[gg] = metals_grid[xbin[gg]][ybin[gg]][zbin[gg]]	
			if ((Mass[gg] >= 9.) & (Mass[gg] <= 11.) & (rho[gg] > 0.5) & (rho[gg] < 2.)) :
				xyz[xbin[gg]][ybin[gg]][zbin[gg]] = xyz[xbin[gg]][ybin[gg]][zbin[gg]] + 1	
		
		print 'cells in galaxy subset'
		wxyz = np.where(xyz > 0)[0]
		print len(wxyz), 'cells have at least one galaxy'
		wxyz = np.where(xyz > 1)[0]
		print len(wxyz), 'cells have more than one galaxy'
		wxyz = np.where(xyz > 2)[0]
		print len(wxyz), 'cells have more than two galaxies'
		print
		print


		halo = np.where(halo > 0., halo, 1.0e7)
		zhalo = np.where(halo > 0., zhalo, 1.0)
		zgrid = np.where(IGM > 0., zgrid, 0.0)
		IGM = np.where(IGM > 0., IGM, 1.0)
		IGM = np.where(IGM > ejected, IGM, ejected)
		
		
		print len(w), 'galaxies loaded'	
		w = np.where((Mass >= 9.) & (Mass <= 11.) & (IGM > 0.) & (rho > 0.5) & (rho < 2.))[0]
		w1 = np.where((Mass >= 9.) & (Mass <= 11.) & (rho > 0.5) & (rho < 2.))[0]
		w2 = np.where((Mass >= 9.) & (Mass <= 11.) & (rho > 0.5) & (rho < 2.) & (zgrid > 0.))[0]
		print 'There are', len(w), 'galaxies with'
		print np.min(Mass[w]), '< Mstars <', np.max(Mass[w])
		print np.min(rho[w]), '< rho <', np.max(rho[w])
		print float(len(w2))/float(len(w1))*100., '% of galaxies have ejected metals'
		
		
		w = np.where((Mass >= 9.) & (Mass <= 11.) & (rho > 0.5) & (rho < 2.) & (zgrid/IGM/0.02 > 1.0e-4))[0]
		print 'Their densities are', np.min(rho[w]), np.max(rho[w])
		lg_rho = rho[w]
		print 'Their IGM metallicities are', np.min(zgrid[w]/IGM[w]), np.max(zgrid[w]/IGM[w])
		print 'Their IGM metallicities are', np.log10(np.min(zgrid[w]/IGM[w])/0.02), np.log10(np.max(zgrid[w]/IGM[w])/0.02), 'log solar'
		print 'Median metallicity is', np.average(np.log10(zgrid[w]/IGM[w]/0.02)), 'log solar'
		print 'Numpy Median metallicity is', np.log10(np.median(zgrid[w]/IGM[w]/0.02)), 'log solar'
		print 'Average metallicity is', np.log10(np.average(zgrid[w]/IGM[w]/0.02)), 'log solar'
		lg_Z_IGM = zgrid[w]/IGM[w]
		lg_Z_halo = zhalo[w]/halo[w]
		lg_Z_stars = zstars[w]/stars[w]
		print 'location of the highest metallicity:'
		xb = xbin[np.where(zgrid/IGM == np.max(zgrid[w]/IGM[w]))]
		yb = ybin[np.where(zgrid/IGM == np.max(zgrid[w]/IGM[w]))]
		zb = zbin[np.where(zgrid/IGM == np.max(zgrid[w]/IGM[w]))]
		rhob = rho[np.where(zgrid/IGM == np.max(zgrid[w]/IGM[w]))]
		mb = Mass[np.where(zgrid/IGM == np.max(zgrid[w]/IGM[w]))]
		mvb = halo[np.where(zgrid/IGM == np.max(zgrid[w]/IGM[w]))]
		Zb = zgrid[np.where(zgrid/IGM == np.max(zgrid[w]/IGM[w]))]
		Ib = IGM[np.where(zgrid/IGM == np.max(zgrid[w]/IGM[w]))]
		eb = ejected[np.where(zgrid/IGM == np.max(zgrid[w]/IGM[w]))]
		
		if ir == 0 :
			edisk = 3.0
		if ir == 1 :
			edisk = 150./vmax[np.where(zgrid/IGM == np.max(zgrid[w]/IGM[w]))]
		if ir == 2 :
			cg = coldgas[np.where(zgrid/IGM == np.max(zgrid[w]/IGM[w]))]
			rd = rdisk[np.where(zgrid/IGM == np.max(zgrid[w]/IGM[w]))]
			st = stars[np.where(zgrid/IGM == np.max(zgrid[w]/IGM[w]))]
			edisk = (cg/(2.*np.pi*rd**2)/1.6e15)**(0.6) * (cg/(cg + st)/0.12)**(0.8)
			
		vv = vvir[np.where(zgrid/IGM == np.max(zgrid[w]/IGM[w]))]	
		vm = vmax[np.where(zgrid/IGM == np.max(zgrid[w]/IGM[w]))]
		if ie == 0 :
			ehalo = 0.3 * (680./vv)**2
		if ie == 1 :
			ehalo = 0.3 * (3. * vm/vv)**2
		if ie == 2 :
			ehalo = edisk * 0.1 * (680./vv)**2
			
		SFR = sfr[np.where(zgrid/IGM == np.max(zgrid[w]/IGM[w]))]
		
		print xb, yb, zb
		print 'Mass and density:'
		print mb, rhob, mvb
		print 'SFR, ehalo-edisk'
		print SFR, ehalo-edisk
		print 'M_Z, M_ejected, M_IGM'
		print Zb, eb, Ib
		
		print np.sum(zejected[w])/np.sum(zejected)*100., '% of metals in the IGM are in this range'

		print

		w = np.where(zgrid > 0.)[0]
		Zmstar = zgrid[np.where(zgrid > 0.)]/IGM[np.where(zgrid > 0.)]
		print 'For all galaxies:'
		print 'Their IGM metallicities are', np.min(Zmstar[np.where(Zmstar > 0.)]), '-', np.max(zgrid[w]/IGM[w])
		print 'Their IGM metallicities are', np.log10(np.min(Zmstar[np.where(Zmstar > 0.)])/0.02), np.log10(np.max(zgrid[w]/IGM[w])/0.02), 'log solar'
		print 'Numpy Median metallicity is', np.log10(np.median(Zmstar[np.where(Zmstar > 0.)]/0.02)), 'log solar'
		
		print
		print
# 
# 		die

		res.gal_Z_vs_IGM_Z(mod, sim, dir, lg_Z_halo, lg_Z_stars, lg_Z_IGM)


# 		w = np.where((10.0 < Mass) & ( Mass < 10.5) & (0.75 < rho) & ( rho < 1.25))[0]
# 		print 'There are', len(w), 'M* -0.5 galaxies with 0.75 < log(rho/rhobar) < 1.25'
# 		w = np.where((10.0 < Mass) & ( Mass < 10.5) & (0.75 < rho) & ( rho < 1.25) &(ejected > 0.))[0]
# 		print 'There are', len(w), 'M* -0.5 galaxies with 0.75 < log(rho/rhobar) < 1.25 and ejected IGM'
# 		Zmstar = zgrid[w]/IGM[w]
# 		print 'Their IGM metallicities are', np.min(Zmstar[np.where(Zmstar > 0.)]), np.max(zgrid[w]/IGM[w])
# 		print 'Their IGM metallicities are', np.log10(np.min(Zmstar[np.where(Zmstar > 0.)])/0.02), np.log10(np.max(zgrid[w]/IGM[w])/0.02), 'log solar'
# 		print 'Median metallicity is', np.average(np.log10(Zmstar[np.where(Zmstar > 0.)]/0.02)), 'log solar'
# 		print
# 		print 'Mvir of these halos:', np.log10(min(mvir[w])), np.log10(max(mvir[w]))
# 		print 'Med Mvir of these halos:', np.average(np.log10(mvir[w]))
# 
# 		print np.min(Mass), np.max(Mass)	
# 		print np.min(rho), np.max(rho)
		massrange = 4.3 #np.max(np.log10(mvir)) - np.min(np.log10(mvir))
# 		print 'massrange = ', massrange
		rhorange = 3.62 #np.max(rho) - np.min(rho)
		maxmass = 11.4 #np.max(Mass)
		maxrho = 3.32 # np.max(rho)
		if mass == 'mstars' : maxmass = 11.3
# 
# 
		H, yedges, xedges = np.histogram2d(Mass, rho, bins = 20, range = ((maxmass - massrange, maxmass), (maxrho - rhorange, maxrho)))
		logH = np.log10(H)
		print 'min, max log(H):', np.min(np.where(logH > -8.)), np.max(logH)
# 
# 		eH, yedges, xedges = np.histogram2d(Mass, rho, bins = 20, weights = ejected, range = ((maxmass - massrange, maxmass), (maxrho - rhorange, maxrho)))
# 		logeH = np.log10(eH)
# 		print 'min, max log(eH):', np.min(np.where(logeH > -8.)), np.max(logeH)
# 
		eZH, yedges, xedges = np.histogram2d(Mass, rho, bins = 20, weights = zejected, range = ((maxmass - massrange, maxmass), (maxrho - rhorange, maxrho)))
		logeZH = np.log10(eZH)
		print 'min, max log(eZH):', np.min(np.where(logeZH > -8.)), np.max(logeZH)
# 
		IH, yedges, xedges = np.histogram2d(Mass, rho, bins = 20, weights = IGM, range = ((maxmass - massrange, maxmass), (maxrho - rhorange, maxrho)))
		logIH = np.log10(IH)
		print 'min, max log(IH):', np.min(np.where(logIH > -8.)), np.max(logIH)

		ZH, yedges, xedges = np.histogram2d(Mass, rho, bins = 20, weights = zgrid/IGM/0.02, range = ((maxmass - massrange, maxmass), (maxrho - rhorange, maxrho)))
		logZH = np.log10(ZH)
		print 'min, max log(ZH):', np.min(np.where(logZH > -8.)), np.max(logZH)
# 
# 		hH, yedges, xedges = np.histogram2d(Mass, rho, bins = 20, weights = halo, range = ((maxmass - massrange, maxmass), (maxrho - rhorange, maxrho)))
# 		loghH = np.log10(hH)
# 		print 'min, max log(hH):', np.min(np.where(loghH > -8.)), np.max(loghH)
# 
# 		ZhH, yedges, xedges = np.histogram2d(Mass, rho, bins = 20, weights = zhalo/halo/0.02, range = ((maxmass - massrange, maxmass), (maxrho - rhorange, maxrho)))
# 		logZhH = np.log10(ZhH)
# 		print 'min, max log(ZhH):', np.min(np.where(logZhH > -8.)), np.max(logZhH)
# 
# 
# 		zIGM = zgrid/IGM/0.02
# 		zIGM = np.where(zIGM > 0., zIGM, 1.0e-8)
# 		nIGM, edges = np.histogram(rho, bins = 20, range = (maxrho - rhorange, maxrho))
# 		mh, edges = np.histogram(rho, bins = 20, weights = IGM, range = (maxrho - rhorange, maxrho))
# 		zh, edges = np.histogram(rho, bins = 20, weights = zgrid, range = (maxrho - rhorange, maxrho))
# 		Zh, edges = np.histogram(rho, bins = 20, weights = zIGM, range = (maxrho - rhorange, maxrho)) #zh/mh/0.02
# 
# 		print len(np.where(IGM < 0.)[0]), 'cells have less than 0 M'
# 		print len(np.where((IGM <= 0.) & ( zgrid > 0.))[0]), 'cells have metals but no igm'
# 
# 

		w = np.where((Mass >= 9.) & (Mass <= 11.) & (IGM > 0.) & (rho > 0.5) & (rho < 2.))[0]
		Zbins = np.linspace(-3.8, 0.2, 21)
		ZZhist, edges = np.histogram(np.log10(zgrid[w]/IGM[w]/0.02), bins = Zbins)
		off = (edges[1] - edges[0])/2.
		ZZ_x = np.linspace(edges[0] + off, edges[-1] - off, 20)

		Zhistfile = dir + sim + '_hist_ZZ.dat'
		print 'printing hist vs ZZ data to', Zhistfile
		f = open(Zhistfile, 'w')
		for rr in range(len(ZZ_x)) :
			f.write(str(ZZ_x[rr]) + '\t' + str(ZZhist[rr]) + '\n')
		f.close()

		
		Zrho, edges = np.histogram(rho, bins = xedges, weights = zgrid)
		IGMrho, edges = np.histogram(rho, bins = xedges, weights = IGM)
		ZZrho = Zrho/IGMrho
		rho_hist, rho_bins = np.histogram(rho, bins = xedges)
		off = abs(rho_bins[0] - rho_bins[1])/2.
		rho_x = np.linspace(rho_bins[0] + off, rho_bins[-1] - off, 20)
		
		rhofile = dir + sim + '_rho_Z.dat'
		print 'printing rho vs Z data to', rhofile
		f = open(rhofile, 'w')
		for rr in range(len(rhox)) :
			f.write(str(rho_x[rr]) + '\t' + str(Zrho[rr]) + '\n')
		f.close()
		
		Zrhofile = dir + sim + '_rho_ZZ.dat'
		print 'printing rho vs ZZ data to', Zrhofile
		f = open(Zrhofile, 'w')
		for rr in range(len(rhox)) :
			f.write(str(rho_x[rr]) + '\t' + str(ZZrho[rr]) + '\n')
		f.close()

		Mass_hist, Mass_bins = np.histogram(Mass, bins = yedges)
		off = abs(Mass_bins[1] - Mass_bins[0])/2.
		Mass_x = np.linspace(Mass_bins[0] + off, Mass_bins[-1] - off, 20)

		ZMstar, m_bins = np.histogram(Mass, bins = yedges, weights = zgrid)
		IGMMstar, m_bins = np.histogram(Mass, bins = yedges, weights = IGM)
		ZZMstar = ZMstar/IGMMstar

		off = abs(m_bins[0] - m_bins[1])/2.
		m_x = np.linspace(m_bins[0] + off, m_bins[-1] - off, 20)

		print np.log10(ZMstar)
		
		mfile = dir + sim + '_Mstar_Z.dat'
		print 'printing Mstar vs Z data to', mfile
		f = open(mfile, 'w')
		for rr in range(len(m_x)) :
			f.write(str(m_x[rr]) + '\t' + str(ZMstar[rr]) + '\t' + str(Mass_hist[rr]) + '\n')
		f.close()

		Zmfile = dir + sim + '_Mstar_ZZ.dat'
		print 'printing Mstar vs ZZ data to', Zmfile
		f = open(Zmfile, 'w')
		for rr in range(len(m_x)) :
			f.write(str(m_x[rr]) + '\t' + str(ZZMstar[rr]) + '\t' + str(Mass_hist[rr]) + '\n')
		f.close()

# 
# 		fig = plt.figure(figsize=(12.,10))							# square box
# 		fig.subplots_adjust(top = 0.85, bottom = .15, left = .15, right = .95)
# 		ax = fig.add_subplot(1,1,1)
# 
# 		ax.set_yscale('log')
# 		plt.axis([min(m_bins),max(m_bins),1.0e12, 1.0e15])
# 
# 		for tick in ax.xaxis.get_major_ticks():
# 			tick.label1.set_fontsize(30)
# 		for tick in ax.yaxis.get_major_ticks():
# 			tick.label1.set_fontsize(30)
# 
# 		ax.plot(m_x, IGMMstar, lw = 5, c = IGMgas_color[2], label = r'IGM($\rho$)')
# 		ax.plot(m_x, IGMMstar, lw = 5, c = IGMgas_color[0], label = r'IGM(M$_{\star}$)')
# 
# 
# 		ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0, ncol = 1)
# 
# 
# 		ax.set_xlabel(r'log$_{10}$(M$_{\star}$)', fontsize=34)  #36 for stripe
# 		ax.set_ylabel(r'M(IGM)', fontsize=34) #54 for square
# 
# 		ax2 = ax.twiny() # ax2 is responsible for "top" axis and "right" axis
# 		for tick in ax2.xaxis.get_major_ticks():
# 			tick.label2.set_fontsize(30)
# 
# 		plt.axis([min(rho_bins),max(rho_bins),1.0e12, 1.0e15])
# 		ax2.set_xlabel(r'log$_{10}$($\rho/\bar{\rho}$)' + "\n", fontsize=34)  #36 for stripe
# 		ax2.plot(rho_x, IGMrho, lw = 5, c = IGMgas_color[2])
# #		ax2.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0, ncol = 1)
# 
# 
# 		outputFile = 'IGM_hists.png'
# 		plt.savefig(outputFile)  # Save the figure
# 		print 'Saved file to', outputFile
# 		plt.close()
# 
# 		
# 		
# 		logHnorm = np.log10(H/rho_hist)
# 		logeHnorm = np.log10(eH/rho_hist)
# 		logeZHnorm = np.log10(eH/rho_hist)
# 
		res.mstars_rho_histogram(mod, sim, dir, logH, xedges, yedges, cmap, vlim, maxmass, massrange, maxrho, rhorange, rho_hist, rho_x, Mass_hist, Mass_x)
# 		res.MZ_ejected_mstars_rho_histogram(mod, sim, dir, logeZH, logH, xedges, yedges, cmapZ, vlimZH, zz)
		res.Z_Zsun_mstars_rho_histogram(mod, sim, dir, H, logH, eZH, logZH, xedges, yedges, cmapZZ, vlimZZ, zz)
		res.IGM_mstars_rho_histogram(mod, sim, dir, logIH, logH, xedges, yedges, cmapI, vlimI, zz)

	def mstars_rho_histogram(self, mod, sim, dir, logH, xedges, yedges, cmap, vlim, maxmass, massrange, maxrho, rhorange, rho_hist, rho_x, Mass_hist, Mass_x) :

		fig = plt.figure(figsize=(13.5,10))
		fig.subplots_adjust(left = .12, right = .89, top = 0.95, bottom = 0.15) 
		ax = fig.add_subplot(1,1,1)
	
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(30)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(30)

		extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]
		im = plt.imshow(logH,extent=extent,interpolation='nearest',origin='lower', cmap = cmap, vmin = vlim[0], vmax = vlim[1], aspect = 'auto')
#		ax.set_aspect(0.5)
#		plt.text(xedges[0] + 0.2, yedges[-1] - 0.2, 'z = ' + zz, size = 40, color = 'k', verticalalignment='top', horizontalalignment='left')
		
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)

		cbar = plt.colorbar(im, cax=cax)
		cbar.set_label(r'log$_{10}$(N$_{\mathrm{Galaxies}}$)')
		cbar.set_ticks([0,0.5,1,1.5,2, 2.5, 3, 3.5])

		ax.set_ylim(maxmass - massrange, maxmass)

		plt.axes(extent)
		plt.text(-1, 0.05, mod, size = 40, color = 'k', verticalalignment='bottom', horizontalalignment='right')
		ax.set_xlabel(r'log$_{10}$($\rho/\bar{\rho}$)', fontsize=34)
		ax.set_ylabel(r'log$_{10}$($\mathrm{M}_{\star}/h^{-1}$ M$_{\odot}$)', fontsize=34) 
#		ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0, ncol = 2)

		ax2 = ax.twinx()
		rho_histN = rho_hist/float(np.max(rho_hist))*massrange*0.1 + maxmass - massrange
		ax2.plot(rho_x*0.95, rho_histN, lw = 5, c = 'k')
		ax2.set_xlim(xedges[0], xedges[-1])
		ax2.set_ylim(yedges[0], yedges[-1])
		ax2.set_yticks([])

		ax3 = ax.twiny()
		Mass_hist = maxrho*0.94 - Mass_hist/float(np.max(Mass_hist))*rhorange*0.1
		ax3.plot(Mass_hist, Mass_x, lw = 5, c = 'k', ls = '--')
		ax3.set_xlim(xedges[0], xedges[-1])
#		ax3.set_ylim(yedges[0],yedges[-1])
		ax3.set_ylim(maxmass - massrange, maxmass)
		ax3.set_xticks([])

#		ax4 = ax.twiny()
#		plt.axes(extent)


		outputFile = dir + sim + '_mstars_vs_rho' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()
		
	def MZ_ejected_mstars_rho_histogram(self, mod, sim, dir, logeZH, logH, xedges, yedges, cmapZ, vlimZH, zz) :
		lg10 = r'log$_{10}$'
		h = r'$h^{-1}$'
		MZej = r'$\mathrm{M}_{\mathrm{ejected}}^{Z}$'
		label = lg10 + '(' + h + MZej + ')'

		fig = plt.figure(figsize=(13.5,10))
		fig.subplots_adjust(left = .12, right = .89, top = 0.95, bottom = 0.15) 
		ax = fig.add_subplot(1,1,1)
	
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(30)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(30)

		extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]
		im = plt.imshow(logeZH - logH,extent=extent,interpolation='nearest',origin='lower', cmap = cmapZ, vmin = vlimZH[0], vmax = vlimZH[1], aspect = 'auto')

#		plt.text(xedges[0] + 0.2, yedges[-1] - 0.2, 'z = ' + zz, size = 40, color = 'k', verticalalignment='top', horizontalalignment='left')

		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)

		cbar = plt.colorbar(im, cax=cax)
		label = lg10 + '(' + h + MZej + '/Galaxy)'
		cbar.set_label(label)
		cbar.set_ticks([4, 5, 6, 7, 8])


		ax.set_xlabel(r'log$_{10}$($\rho/\bar{\rho}$)', fontsize=34)
		ax.set_ylabel(r'log$_{10}$($\mathrm{M}_{\star}/h^{-1}$ M$_{\odot}$)', fontsize=34) 
		plt.text(-1, 0.05, mod, size = 40, color = 'k', verticalalignment='bottom', horizontalalignment='right')
#		ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0, ncol = 2)

		outputFile = dir + sim + '_' + 'MZ_ejected_' + 'Mstars' + '_vs_rho_per_halo' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()

	def Z_Zsun_mstars_rho_histogram(self, mod, sim, dir, H, logH, eZH, logZH, xedges, yedges, cmapZZ, vlimZZ, zz) :

		fig = plt.figure(figsize=(13.5,10))
		fig.subplots_adjust(left = .12, right = .89, top = 0.95, bottom = 0.15) 
		ax = fig.add_subplot(1,1,1)
	
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(30)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(30)

		extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]
		lgH = np.where(eZH > 0., logH, 100)
		Zh = np.where(H > 5, logZH, -100.)
		im = plt.imshow(Zh - lgH,extent=extent,interpolation='nearest',origin='lower', cmap = cmapZZ, vmin = vlimZZ[0], vmax = vlimZZ[1], aspect = 'auto')
#		plt.text(xedges[0] + 0.2, yedges[-1] - 0.2, 'z = ' + zz, size = 40, color = 'k', verticalalignment='top', horizontalalignment='left')
		print 'min, max Z/Zsun = ', np.min(Zh - lgH), np.max(Zh - lgH)
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)

		cbar = plt.colorbar(im, cax=cax)
		lg10 = r'log$_{10}$'
		h = r'$h^{-1}$'
		Z = r'<$\mathrm{Z}/\mathrm{Z}_{\odot}$>'
		label =  lg10 + '('  + Z + ')' 
		cbar.set_label(label)
		cbar.set_ticks([-5, -4, -3, -2, -1, 0, 1])


		ax.set_xlabel(r'log$_{10}$($\rho/\bar{\rho}$)', fontsize=34)
		ax.set_ylabel(r'log$_{10}$($\mathrm{M}_{\star}/h^{-1}$ M$_{\odot}$)', fontsize=34) 
#		ax.set_ylabel(ylab, fontsize=34) 
		print sim
		print plt.gca().get_xlim()
		plt.text(-1, 0.05, mod, size = 40, color = 'k', verticalalignment='bottom', horizontalalignment='right')


		outputFile = dir + sim + '_' + 'Z_Zsun_' + 'Mstars' + '_vs_rho' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()

	def IGM_mstars_rho_histogram(self, mod, sim, dir, logIH, logH, xedges, yedges, cmapI, vlimI, zz) :

		output_list = []
		fig = plt.figure(figsize=(13.5,10))
		fig.subplots_adjust(left = .12, right = .89, top = 0.95, bottom = 0.15) 
		ax = fig.add_subplot(1,1,1)
	
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(30)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(30)

		extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]
		im = plt.imshow(logIH - logH,extent=extent,interpolation='nearest',origin='lower', cmap = cmapI, vmin = vlimI[0], vmax = vlimI[1], aspect = 'auto')
#		plt.text(xedges[0] + 0.2, yedges[-1] - 0.2, 'z = ' + zz, size = 40, color = 'k', verticalalignment='top', horizontalalignment='left')
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)



		cbar = plt.colorbar(im, cax=cax)
		lg10 = r'log$_{10}$'
		h = r'$h^{-1}$'
		Z = r'<$\mathrm{Z}/\mathrm{Z}_{\odot}$>'
		label = lg10 + '(M' + r'$_{\mathrm{IGM}}$' + '/' + h + 'M' + r'$_{\odot}$' + ' per Galaxy)'
		cbar.set_label(label)
		cbar.set_ticks([9, 10, 11, 12])


		ax.set_xlabel(r'log$_{10}$($\rho/\bar{\rho}$)', fontsize=34)
		ax.set_ylabel(r'log$_{10}$($\mathrm{M}_{\star}/h^{-1}$ M$_{\odot}$)', fontsize=34) 
		plt.text(-1, 0.05, mod, size = 40, color = 'k', verticalalignment='bottom', horizontalalignment='right')

		outputFile = dir + sim + '_' + 'IGM_' + 'Mstars' + '_vs_rho_per_halo' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()
		output_list.append(outputFile)

# 		fig = plt.figure(figsize=(13.5,10))
# 		fig.subplots_adjust(left = .12, right = .89, top = 0.95, bottom = 0.15) 
# 		ax = fig.add_subplot(1,1,1)
# 	
# 		for tick in ax.xaxis.get_major_ticks():
# 			tick.label1.set_fontsize(30)
# 		for tick in ax.yaxis.get_major_ticks():
# 			tick.label1.set_fontsize(30)
# 
# 		extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]
# 		im = plt.imshow(logZhH - logH,extent=extent,interpolation='nearest',origin='lower', cmap = cmapZhZ, vmin = vlimZhZ[0], vmax = vlimZhZ[1], aspect = 'auto')
# 		plt.text(xedges[0] + 0.2, yedges[-1] - 0.2, 'z = ' + zz, size = 40, color = 'k', verticalalignment='top', horizontalalignment='left')
# 		print 'min, max Z/Zsun = ', np.min(Zh - lgH), np.max(Zh - lgH)
# 		divider = make_axes_locatable(ax)
# 		cax = divider.append_axes("right", size="5%", pad=0.05)
# 
# 		cbar = plt.colorbar(im, cax=cax)
# 		lg10 = r'log$_{10}$'
# 		h = r'$h^{-1}$'
# 		Z = r'<$\mathrm{Z}/\mathrm{Z}_{\odot}$>'
# 		label =  lg10 + '('  + Z + ')' 
# 		cbar.set_label(label)
# 		cbar.set_ticks([-5, -4, -3, -2, -1, 0, 1])
# 
# 
# 		ax.set_xlabel(r'log$_{10}$($\rho/\bar{\rho}$)', fontsize=34)
# 		ylab = r'log$_{10}$($h^{-1}$' + Mlabel + ')'
# 		ax.set_ylabel(ylab, fontsize=34) 
# 		print sim
# 		print plt.gca().get_xlim()
# 		plt.text(-1, 0.05, mod, size = 40, color = 'k', verticalalignment='bottom', horizontalalignment='right')
# 
# 		outputFile = dir + sim + '_' + 'Z_Zsun_hothalo_' + Mname + '_vs_rho' + OutputFormat
# 		plt.savefig(outputFile)  # Save the figure
# 		print 'Saved file to', outputFile
# 		plt.close()
# 		output_list.append(outputFile)
# 
# 		fig = plt.figure(figsize=(13.5,10))
# 		fig.subplots_adjust(left = .12, right = .89, top = 0.95, bottom = 0.15) 
# 		ax = fig.add_subplot(1,1,1)
# 	
# 		for tick in ax.xaxis.get_major_ticks():
# 			tick.label1.set_fontsize(30)
# 		for tick in ax.yaxis.get_major_ticks():
# 			tick.label1.set_fontsize(30)
# 
# 		extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]
# 		im = plt.imshow(loghH - logH,extent=extent,interpolation='nearest',origin='lower', cmap = cmapZhZ, vmin = vlimI[0], vmax = vlimI[1] + 1., aspect = 'auto')
# 		plt.text(xedges[0] + 0.2, yedges[-1] - 0.2, 'z = ' + zz, size = 40, color = 'k', verticalalignment='top', horizontalalignment='left')
# 		divider = make_axes_locatable(ax)
# 		cax = divider.append_axes("right", size="5%", pad=0.05)
# 
# 
# 
# 		cbar = plt.colorbar(im, cax=cax)
# 		label = lg10 + '(' + h + 'M' + r'$_{\mathrm{Hot Gas}}$' + '/Galaxy)'
# 		cbar.set_label(label)
# 		cbar.set_ticks([9, 10, 11, 12])
# 
# 
# 		ax.set_xlabel(r'log$_{10}$($\rho/\bar{\rho}$)', fontsize=34)
# 		ax.set_ylabel(r'log$_{10}$($h^{-1}$' + Mlabel + ')', fontsize=34) 
# 		plt.text(-1, 0.05, mod, size = 40, color = 'k', verticalalignment='bottom', horizontalalignment='right')
# #		ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0, ncol = 2)
# 
# 		outputFile = dir + sim + '_' + 'hothalo_' + Mname + '_vs_rho_per_halo' + OutputFormat
# 		plt.savefig(outputFile)  # Save the figure
# 		print 'Saved file to', outputFile
# 		plt.close()
# 		output_list.append(outputFile)
# 
# 

		print output_list, 'printed to file'

	def gal_Z_vs_IGM_Z(self, mod, sim, dir, lg_Z_halo, lg_Z_stars, lg_Z_IGM) :
		fig = plt.figure(figsize=(10.,10))
		fig.subplots_adjust(left = .15, right = .87) 
		ax = fig.add_subplot(1,1,1)
	
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(30)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(30)
			
		log_lg_Z_h = np.log10(lg_Z_halo/0.02)
		log_lg_Z_s = np.log10(lg_Z_stars/0.02)
		log_lg_Z_I = np.log10(lg_Z_IGM/0.02)

#		plt.axes([xedges[0], xedges[-1], -4.25, -0.75])
#		ax.set_xlim(xedges[0], xedges[-1])
		ax.set_xlim(-1., 0.5)
		ax.set_ylim(-2.5, -0.5)
		ax.set_xticks([ -1, -0.5, 0, 0.5])
		ax.set_yticks([ -2, -1.5, -1, -0.5])
		ax.axhline(-1.10, c = 'k', lw = 1)
		ax.axhline(-1.14, c = 'k', lw = 1)
#		ax.axhline(-1.06, c = 'k', lw = 1, ls = '--')
#		ax.axhline(-1.18, c = 'k', lw = 1, ls = '--')
		ax.axvline(-0.29, c = 'k', lw = 1)
		ax.axvline(-0.13, c = 'k', lw = 1)
		
#		print log_lg_Z_h
		w = np.where((-0.29 <= log_lg_Z_s) & (log_lg_Z_s <= -0.13) & (-1.14 <= log_lg_Z_I) & (log_lg_Z_I <= -1.10))[0]
		print len(w), 'galaxies out of', len(lg_Z_halo), '(', float(len(w))/float(len(lg_Z_halo))*100., '%) are within the errorbars'
		w = np.where((-0.29 <= log_lg_Z_s) & (log_lg_Z_s <= -0.13) & (-1.18 <= log_lg_Z_I) & (log_lg_Z_I <= -1.06))[0]
		print len(w), 'galaxies out of', len(lg_Z_halo), '(', float(len(w))/float(len(lg_Z_halo))*100., '%) are within 3 sigma'
		
		colors = np.log10(lg_Z_halo/0.02)
		print 'max halo Z:', np.max(colors)
		print 'max IGM Z:', np.max(np.log10(lg_Z_IGM/0.02))
		im = ax.scatter(log_lg_Z_s, log_lg_Z_I, c =  stars_color[1],  edgecolor = 'white', s = 100) #, cmap = cmapZhZ, vmin = vlimZhZ[0] + 2., vmax = vlimZhZ[1], s = 100)
		ax.set_xlabel(r'log$_{10}$(Z$_{\star}$/$\mathrm{Z}_{\odot}$)', fontsize=34) 
		ax.set_ylabel(r'log$_{10}$(Z$_{\mathrm{IGM}}$/$\mathrm{Z}_{\odot}$)', fontsize=34) 
		plt.text(-1, 0.05, mod, size = 40, color = 'k', verticalalignment='bottom', horizontalalignment='right')

		outputFile = dir + sim + '_' + 'IGM_metals_vs_Stellar_metals' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()

	def StellarMassFunction(self, nreheat, neject, snapnum) :

		SMF = np.zeros((nreheat, neject, 130))
		Ms = np.zeros((nreheat, neject, 130))
		snap = "%03d" % snapnum

		for ir in range(nreheat) :
			for ie in range(neject) :
				dir = 'millennium_W' + str(ir + 1) + str(ie + 1) + '/' + 'plots' + '/' 
				sim = 'W' +  str(ir + 1) + str(ie + 1)
				mod = 'R' +  str(ir + 1) + 'E' + str(ie + 1)
		
				SMfile = dir + 'SMF.txt'
				print 'reading SMF data from', SMfile
				ii = 0
				for item in file(SMfile) :
					item = item.split()
					Ms[ir][ie][ii] = float(item[0])
					SMF[ir][ie][ii] = float(item[1])
					ii = ii + 1

		Baldry = np.array([
			[7.05, 1.3531e-01, 6.0741e-02],
			[7.15, 1.3474e-01, 6.0109e-02],
			[7.25, 2.0971e-01, 7.7965e-02],
			[7.35, 1.7161e-01, 3.1841e-02],
			[7.45, 2.1648e-01, 5.7832e-02],
			[7.55, 2.1645e-01, 3.9988e-02],
			[7.65, 2.0837e-01, 4.8713e-02],
			[7.75, 2.0402e-01, 7.0061e-02],
			[7.85, 1.5536e-01, 3.9182e-02],
			[7.95, 1.5232e-01, 2.6824e-02],
			[8.05, 1.5067e-01, 4.8824e-02],
			[8.15, 1.3032e-01, 2.1892e-02],
			[8.25, 1.2545e-01, 3.5526e-02],
			[8.35, 9.8472e-02, 2.7181e-02],
			[8.45, 8.7194e-02, 2.8345e-02],
			[8.55, 7.0758e-02, 2.0808e-02],
			[8.65, 5.8190e-02, 1.3359e-02],
			[8.75, 5.6057e-02, 1.3512e-02],
			[8.85, 5.1380e-02, 1.2815e-02],
			[8.95, 4.4206e-02, 9.6866e-03],
			[9.05, 4.1149e-02, 1.0169e-02],
			[9.15, 3.4959e-02, 6.7898e-03],
			[9.25, 3.3111e-02, 8.3704e-03],
			[9.35, 3.0138e-02, 4.7741e-03],
			[9.45, 2.6692e-02, 5.5029e-03],
			[9.55, 2.4656e-02, 4.4359e-03],
			[9.65, 2.2885e-02, 3.7915e-03],
			[9.75, 2.1849e-02, 3.9812e-03],
			[9.85, 2.0383e-02, 3.2930e-03],
			[9.95, 1.9929e-02, 2.9370e-03],
			[10.05, 1.8865e-02, 2.4624e-03],
			[10.15, 1.8136e-02, 2.5208e-03],
			[10.25, 1.7657e-02, 2.4217e-03],
			[10.35, 1.6616e-02, 2.2784e-03],
			[10.45, 1.6114e-02, 2.1783e-03],
			[10.55, 1.4366e-02, 1.8819e-03],
			[10.65, 1.2588e-02, 1.8249e-03],
			[10.75, 1.1372e-02, 1.4436e-03],
			[10.85, 9.1213e-03, 1.5816e-03],
			[10.95, 6.1125e-03, 9.6735e-04],
			[11.05, 4.3923e-03, 9.6254e-04],
			[11.15, 2.5463e-03, 5.0038e-04],
			[11.25, 1.4298e-03, 4.2816e-04],
			[11.35, 6.4867e-04, 1.6439e-04],
			[11.45, 2.8294e-04, 9.9799e-05],
			[11.55, 1.0617e-04, 4.9085e-05],
			[11.65, 3.2702e-05, 2.4546e-05],
			[11.75, 1.2571e-05, 1.2571e-05],
			[11.85, 8.4589e-06, 8.4589e-06],
			[11.95, 7.4764e-06, 7.4764e-06],
			], dtype=np.float32)
		
		fig = plt.figure(figsize=(24.,8))							# square box
		fig.subplots_adjust(top = 0.98, bottom = .15, left = .06, right = .98, wspace = 0.)

		for ie in range(neject) :
			ax = fig.add_subplot(1,3,ie + 1)

			for tick in ax.xaxis.get_major_ticks():
				tick.label1.set_fontsize(30)
			for tick in ax.yaxis.get_major_ticks():
				tick.label1.set_fontsize(30)
				if (ie != 0) : tick.label1On = False

			if (ie == 0) : plt.ylabel(r'$\phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$')  # Set the y...
			plt.xlabel(r'$\log_{10} \mathrm{M}_{\star}\ (\mathrm{M}_{\odot})$')  # and the x-axis labels

			plt.yscale('log', nonposy='clip')
			plt.axis([8.1, 11.9, 2.0e-6, 9.0e-2])

			Baldry_xval = np.log10(10 ** Baldry[:, 0])
			Baldry_yvalU = (Baldry[:, 1]+Baldry[:, 2]) * Hubble_h*Hubble_h*Hubble_h
			Baldry_yvalL = (Baldry[:, 1]-Baldry[:, 2]) * Hubble_h*Hubble_h*Hubble_h

			plt.fill_between(Baldry_xval, Baldry_yvalU, Baldry_yvalL, 
				facecolor='gray', alpha=0.25, label='Baldry et al. 2008 (z=0.1)')


			for ir in range(nreheat) :
				mod = 'R' +  str(ir + 1) + 'E' + str(ie + 1)
				ax.plot(Ms[ir][ie], SMF[ir][ie], lw = 3 + ie, c = stars_color[ir], ls = line_cycle[ie], label = mod)

			ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0, ncol = 1)



		outputFile =  'SMF' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()

	def rho_hist_Z(self, nreheat, neject, snapnum) :
		
		rhox = np.zeros((nreheat, neject, 20))
		Zrho = np.zeros((nreheat, neject, 20))


		snap = "%03d" % snapnum

		for ir in range(nreheat) :
			for ie in range(neject) :
				dir = 'millennium_W' + str(ir + 1) + str(ie + 1) + '/' + mass + '_rho_' + snap + '/'
				sim = 'W' +  str(ir + 1) + str(ie + 1)
				mod = 'R' +  str(ir + 1) + 'E' + str(ie + 1)
		
				rhofile = dir + sim + '_rho_Z.dat'
				print 'reading rho vs Z data from', rhofile
				ii = 0
				for item in file(rhofile) :
					item = item.split()
					rhox[ir][ie][ii] = float(item[0])
					Zrho[ir][ie][ii] = np.log10(float(item[1]))
					ii = ii + 1
					

		fig = plt.figure(figsize=(12.,10))
		fig.subplots_adjust(left = .15, right = .97, top = 0.95, bottom = 0.15) 
		ax = fig.add_subplot(1,1,1)
	
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(30)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(30)

#		plt.text(xedges[0] + 0.2, yedges[-1] - 0.2, 'z = ' + zz, size = 40, color = 'k', verticalalignment='top', horizontalalignment='left')
		
		ax.set_ylim(8.5, 11.25)

		ax.set_xlabel(r'log$_{10}$($\rho/\bar{\rho}$)', fontsize=34)
		ax.set_ylabel(r'log$_{10}$(M$^{\mathrm{Z}}$/$h^{-1}$ M$_{\odot}$)', fontsize=34) 


		for ir in range(nreheat) :
			for ie in range(neject) :
#				print Zrho[ir][ie] 
				mod = 'R' +  str(ir + 1) + 'E' + str(ie + 1)
				ax.plot(rhox[ir][ie], Zrho[ir][ie], lw = 3 + ie, c = Z_color[ir], ls = line_cycle[ie], label = mod)

		ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0, ncol = 1)
		plt.xticks([0, 0.5, 1, 1.5, 2, 2.5, 3])
#		plt.yticks([-5, -4, -3, -2, -1])


		outputFile =  'Z_vs_rho_model' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()
		
	def rho_hist_ZZ(self, nreheat, neject, snapnum) :
		
		rhox = np.zeros((nreheat, neject, 20))
		Zrho = np.zeros((nreheat, neject, 20))


		snap = "%03d" % snapnum

		for ir in range(nreheat) :
			for ie in range(neject) :
				dir = 'millennium_W' + str(ir + 1) + str(ie + 1) + '/' + mass + '_rho_' + snap + '/'
				sim = 'W' +  str(ir + 1) + str(ie + 1)
				mod = 'R' +  str(ir + 1) + 'E' + str(ie + 1)
		
				rhofile = dir + sim + '_rho_ZZ.dat'
				print 'reading rho vs ZZ data from', rhofile
				ii = 0
				for item in file(rhofile) :
					item = item.split()
					rhox[ir][ie][ii] = float(item[0])
					Zrho[ir][ie][ii] = np.log10(float(item[1])/0.02)
					ii = ii + 1
					

		fig = plt.figure(figsize=(12.,10))
		fig.subplots_adjust(left = .15, right = .97, top = 0.95, bottom = 0.15) 
		ax = fig.add_subplot(1,1,1)
	
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(30)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(30)

#		plt.text(xedges[0] + 0.2, yedges[-1] - 0.2, 'z = ' + zz, size = 40, color = 'k', verticalalignment='top', horizontalalignment='left')
		
		ax.set_ylim(-3.5, -1)

		ax.set_xlabel(r'log$_{10}$($\rho/\bar{\rho}$)', fontsize=34)
		ax.set_ylabel(r'log$_{10}$(Z$_{\mathrm{IGM}}$/Z$_{\odot}$)', fontsize=34) 


		for ir in range(nreheat) :
			for ie in range(neject) :
#				print Zrho[ir][ie] 
				mod = 'R' +  str(ir + 1) + 'E' + str(ie + 1)
				ax.plot(rhox[ir][ie], Zrho[ir][ie], lw = 3 + ie, c = ZZ_color[ir], ls = line_cycle[ie], label = mod)

		ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0, ncol = 3)
		plt.xticks([0, 0.5, 1, 1.5, 2, 2.5, 3])
#		plt.yticks([-5, -4, -3, -2, -1])


		outputFile =  'ZZ_vs_rho_model' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()

	def hist_ZZ(self, nreheat, neject, snapnum) :
		
		ZZx = np.zeros((nreheat, neject, 80))
		ZZhist = np.zeros((nreheat, neject, 80))


		snap = "%03d" % snapnum

		for ir in range(nreheat) :
			for ie in range(neject) :
				dir = 'millennium_W' + str(ir + 1) + str(ie + 1) + '/' + mass + '_rho_' + snap + '/'
				sim = 'W' +  str(ir + 1) + str(ie + 1)
				mod = 'R' +  str(ir + 1) + 'E' + str(ie + 1)
		
				rhofile = dir + sim + '_hist_ZZ.dat'
				print 'reading hist vs ZZ data from', rhofile
				ii = 0
				for item in file(rhofile) :
					item = item.split()
					ZZx[ir][ie][ii] = float(item[0])
					ZZhist[ir][ie][ii] = float(item[1])
					ii = ii + 1
					

		fig = plt.figure(figsize=(12.,10))
		fig.subplots_adjust(left = .15, right = .97, top = 0.95, bottom = 0.15) 
		ax = fig.add_subplot(1,1,1)
	
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(30)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(30)

#		plt.text(xedges[0] + 0.2, yedges[-1] - 0.2, 'z = ' + zz, size = 40, color = 'k', verticalalignment='top', horizontalalignment='left')
		
#		ax.set_ylim(-3.5, -1)
		ax.set_yscale('log')
		ax.set_xlabel(r'log$_{10}$(Z$_{\mathrm{IGM}}$/Z$_{\odot}$)', fontsize=34) 
		ax.set_ylabel(r'N', fontsize=34)


		for ir in range(nreheat) :
			for ie in range(neject) :
#				print Zrho[ir][ie] 
				mod = 'R' +  str(ir + 1) + 'E' + str(ie + 1)
				ax.plot(ZZx[ir][ie], ZZhist[ir][ie], lw = 3 + ie, c = ZZ_color[ir], ls = line_cycle[ie], label = mod)
#				ax.plot(ZZx[ir][ie], ZZhist[ir][ie]/np.max(ZZhist[ir][ie]), lw = 3 + ie, c = ZZ_color[ir], ls = line_cycle[ie], label = mod)
				print mod, ZZx[ir][ie][np.where(ZZhist[ir][ie] == np.max(ZZhist[ir][ie]))[0]]

		ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0, ncol = 3)
#		plt.xticks([0, 0.5, 1, 1.5, 2, 2.5, 3])
		plt.xticks([-4, -3, -2, -1, 0])
		print ZZx[0][0]

		outputFile =  'ZZ_hist_model' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()

	def rho_hist_Z_cumulative(self, nreheat, neject, snapnum) :
		
		rhox = np.zeros((nreheat, neject, 20))
		Zrho = np.zeros((nreheat, neject, 20))


		snap = "%03d" % snapnum

		for ir in range(nreheat) :
			for ie in range(neject) :
				dir = 'millennium_W' + str(ir + 1) + str(ie + 1) + '/' + mass + '_rho_' + snap + '/'
				sim = 'W' +  str(ir + 1) + str(ie + 1)
				mod = 'R' +  str(ir + 1) + 'E' + str(ie + 1)
		
				rhofile = dir + sim + '_rho_Z.dat'
				print 'reading rho vs Z data from', rhofile
				ii = 0
				for item in file(rhofile) :
					item = item.split()
					rhox[ir][ie][ii] = float(item[0])
					Zrho[ir][ie][ii] = float(item[1])
					ii = ii + 1
					

		fig = plt.figure(figsize=(12.,10))
		fig.subplots_adjust(left = .15, right = .97, top = 0.95, bottom = 0.15) 
		ax = fig.add_subplot(1,1,1)
	
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(30)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(30)

#		plt.text(xedges[0] + 0.2, yedges[-1] - 0.2, 'z = ' + zz, size = 40, color = 'k', verticalalignment='top', horizontalalignment='left')
		
		ax.set_ylim(-7, -5)

		ax.set_xlabel(r'log$_{10}$($\rho/\bar{\rho}$)', fontsize=34)
		ax.set_ylabel(r'log$_{10}$($\Omega_{Z}$)', fontsize=34) 

		Zrho_c = np.zeros((nreheat, neject, len(Zrho[0][0])))
		M90 = np.zeros((nreheat, neject))
		for ir in range(nreheat) :
			for ie in range(neject) :
				for ii in range(len(Zrho[0][0])) :
					Zrho_c[ir, ie, ii] = np.sum(Zrho[ir, ie, 0:ii + 1])
				M90[ir][ie] = Zrho_c[ir, ie, -1]*0.9

		for ir in range(nreheat) :
			for ie in range(neject) :
				mn = np.where(Zrho_c[ir][ie] < M90[ir][ie])[0][-1]
#				print mn
#				print Zrho_c[ir][ie][0:mn[-1]]/Zrho_c[ir][ie][-1]
				print 'R', ir + 1, 'E', ie + 1, ':', rhox[ir][ie][mn], '-', rhox[ir][ie][mn + 1]
				print Zrho_c[ir][ie][mn]/Zrho_c[ir][ie][-1], Zrho_c[ir][ie][mn + 1]/Zrho_c[ir][ie][-1]

		for ir in range(nreheat) :
			for ie in range(neject) :
#				print Zrho_c[ir][ie] 
				mod = 'R' +  str(ir + 1) + 'E' + str(ie + 1)
				ax.plot(rhox[ir][ie], np.log10(Zrho_c[ir][ie]) - np.log10(6.774e16), lw = 3 + ie, c = Z_color[ir], ls = line_cycle[ie], label = mod)

		ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0, ncol = 1)
		plt.xticks([0, 0.5, 1, 1.5, 2, 2.5, 3])
#		plt.yticks([-5, -4, -3, -2, -1])


		print 'Omega', np.log10(Zrho_c[:, :, -1]/6.774e16)
		print '90% Omega_Z', np.log10(M90) - np.log10(6.774e16)

		outputFile =  'Z_vs_rho_cumulative' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()
		
	def mass_hist_Z(self, nreheat, neject, snapnum) :
		
		Mx = np.zeros((nreheat, neject, 20))
		ZM = np.zeros((nreheat, neject, 20))


		snap = "%03d" % snapnum

		for ir in range(nreheat) :
			for ie in range(neject) :
				dir = 'millennium_W' + str(ir + 1) + str(ie + 1) + '/' + mass + '_rho_' + snap + '/'
				sim = 'W' +  str(ir + 1) + str(ie + 1)
				mod = 'R' +  str(ir + 1) + 'E' + str(ie + 1)
		
				mfile = dir + sim + '_Mstar_Z.dat'
				print 'reading Mstar vs Z data from', mfile
				ii = 0
				for item in file(mfile) :
					item = item.split()
					Mx[ir][ie][ii] = float(item[0])
					ZM[ir][ie][ii] = (float(item[1]))
					ii = ii + 1
					

		fig = plt.figure(figsize=(12.,10))
		fig.subplots_adjust(left = .15, right = .97, top = 0.95, bottom = 0.15) 
		ax = fig.add_subplot(1,1,1)
	
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(30)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(30)

#		plt.text(xedges[0] + 0.2, yedges[-1] - 0.2, 'z = ' + zz, size = 40, color = 'k', verticalalignment='top', horizontalalignment='left')
		
		ax.set_ylim(8.5, 11.25)

		ax.set_xlabel(r'log$_{10}$(M$_{\star}$/$h^{-1}$ M$_{\odot}$)', fontsize=34)
		ax.set_ylabel(r'log$_{10}$(M$^{\mathrm{Z}}$/$h^{-1}$ M$_{\odot}$)', fontsize=34) 


		for ir in range(nreheat) :
			for ie in range(neject) :
#				print ZM[ir][ie] 
				mod = 'R' +  str(ir + 1) + 'E' + str(ie + 1)
				ax.plot(Mx[ir][ie], np.log10(ZM[ir][ie]), lw = 3 + ie, c = Z_color[ir], ls = line_cycle[ie], label = mod)

		ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0, ncol = 3)
		plt.xticks([7, 8, 9, 10, 11])
#		plt.yticks([-5, -4, -3, -2, -1])


		outputFile =  'Z_vs_Mstar_model' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()
		
	def mass_hist_Z_cumulative(self, nreheat, neject, snapnum) :
		
		Mx = np.zeros((nreheat, neject, 20))
		Nx = np.zeros((nreheat, neject, 20))
		ZM = np.zeros((nreheat, neject, 20))


		snap = "%03d" % snapnum

		for ir in range(nreheat) :
			for ie in range(neject) :
				dir = 'millennium_W' + str(ir + 1) + str(ie + 1) + '/' + mass + '_rho_' + snap + '/'
				sim = 'W' +  str(ir + 1) + str(ie + 1)
				mod = 'R' +  str(ir + 1) + 'E' + str(ie + 1)
		
				mfile = dir + sim + '_Mstar_Z.dat'
				print 'reading Mstar vs Z data from', mfile
				ii = 0
				for item in file(mfile) :
					item = item.split()
					Mx[ir][ie][ii] = float(item[0])
					ZM[ir][ie][ii] = float(item[1])
					Nx[ir][ie][ii] = float(item[2])
					ii = ii + 1
					

		fig = plt.figure(figsize=(12.,10))
		fig.subplots_adjust(left = .15, right = .97, top = 0.95, bottom = 0.15) 
		ax = fig.add_subplot(1,1,1)
	
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(30)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(30)

#		plt.text(xedges[0] + 0.2, yedges[-1] - 0.2, 'z = ' + zz, size = 40, color = 'k', verticalalignment='top', horizontalalignment='left')
		
		ax.set_ylim(-7, -5)

		ax.set_xlabel(r'log$_{10}$(M$_{\star}$/$h^{-1}$ M$_{\odot}$)', fontsize=34)
		ax.set_ylabel(r'log$_{10}$($\Omega_{Z}$)', fontsize=34) 

		ZM_c = np.zeros((nreheat, neject, len(ZM[0][0])))
		N_c = np.zeros((nreheat, neject, len(ZM[0][0])))
		M90 = np.zeros((nreheat, neject))
		for ir in range(nreheat) :
			for ie in range(neject) :
				for ii in range(len(ZM[0][0])) :
					ZM_c[ir, ie, ii] = np.sum(ZM[ir, ie, 0:ii + 1])
					N_c[ir, ie, ii] = np.sum(Nx[ir, ie, 0:ii + 1])
				M90[ir][ie] = ZM_c[ir, ie, -1]*0.9

		for ir in range(nreheat) :
			for ie in range(neject) :
				mn = np.where(ZM_c[ir][ie] < M90[ir][ie])[0][-1]
				print 'R', ir + 1, 'E', ie + 1, ':', Mx[ir][ie][mn], '-', Mx[ir][ie][mn + 1]
				print ZM_c[ir][ie][mn]/ZM_c[ir][ie][-1], ZM_c[ir][ie][mn + 1]/ZM_c[ir][ie][-1]
				print N_c[ir][ie][mn]/N_c[ir][ie][-1], N_c[ir][ie][mn + 1]/N_c[ir][ie][-1]
				if ie == 1 : 
					print ZM[ir][ie]
					print np.sum(ZM[ir][ie]), ZM_c[ir][ie][-1]

		for ir in range(nreheat) :
			for ie in range(neject) :
#				print ZM[ir][ie] 
				mod = 'R' +  str(ir + 1) + 'E' + str(ie + 1)
				print mod + '/R1E1', (ZM_c[ir][ie][-1])/(ZM_c[0][0][-1])
				ax.plot(Mx[ir][ie], np.log10(ZM_c[ir][ie]) - np.log10(6.774e16), lw = 3 + ie, c = Z_color[ir], ls = line_cycle[ie], label = mod)

		ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0, ncol = 3)
		plt.xticks([7, 8, 9, 10, 11])
#		plt.yticks([-5, -4, -3, -2, -1])

		print 'OmegaZ = '
		for ir in range (nreheat) :
			for ie in range(neject) :
				print("R%d E%d = %1.3f e-5" % (ir + 1, ie + 1, ZM_c[ir, ie, -1]/6.774e16*1.e5))
		print '90% Omega_Z', M90 - np.log10(6.774e16)

		outputFile =  'Z_vs_Mstar_cumulative' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()

	def mass_hist_ZZ(self, nreheat, neject, snapnum) :
		
		Mx = np.zeros((nreheat, neject, 20))
		ZM = np.zeros((nreheat, neject, 20))


		snap = "%03d" % snapnum

		for ir in range(nreheat) :
			for ie in range(neject) :
				dir = 'millennium_W' + str(ir + 1) + str(ie + 1) + '/' + mass + '_rho_' + snap + '/'
				sim = 'W' +  str(ir + 1) + str(ie + 1)
				mod = 'R' +  str(ir + 1) + 'E' + str(ie + 1)
		
				mfile = dir + sim + '_Mstar_ZZ.dat'
				print 'reading Mstar vs Z data from', mfile
				ii = 0
				for item in file(mfile) :
					item = item.split()
					Mx[ir][ie][ii] = float(item[0])
					ZM[ir][ie][ii] = float(item[1])/0.02
					ii = ii + 1
					

		fig = plt.figure(figsize=(12.,10))
		fig.subplots_adjust(left = .15, right = .97, top = 0.95, bottom = 0.15) 
		ax = fig.add_subplot(1,1,1)
	
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(30)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(30)
		
		ax.set_ylim(-3.5, -1)

		ax.set_xlabel(r'log$_{10}$(M$_{\star}$/$h^{-1}$ M$_{\odot}$)', fontsize=34)
		ax.set_ylabel(r'log$_{10}$(Z$_{\mathrm{IGM}}$/Z$_{\odot}$)', fontsize=34) 


		for ir in range(nreheat) :
			for ie in range(neject) :
#				print ZM[ir][ie] 
				mod = 'R' +  str(ir + 1) + 'E' + str(ie + 1)
				ax.plot(Mx[ir][ie], np.log10(ZM[ir][ie]), lw = 3 + ie, c = ZZ_color[ir], ls = line_cycle[ie], label = mod)

		ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0, ncol = 3)
		plt.xticks([7, 8, 9, 10, 11])
#		plt.yticks([-4, -3, -2, -1])


		outputFile =  'ZZ_vs_Mstar_model' + OutputFormat
		plt.savefig(outputFile)  # Save the figure
		print 'Saved file to', outputFile
		plt.close()
		
		
# =================================================================


#  'Main' section of code.  This if statement executes if the code is run from the 
#   shell command line, i.e. with 'python allresults.py'

if __name__ == '__main__':
	res = Results()
	
 	snapnum = 63

 	nmodels = 3
 	nreheat = 3
 	neject = 3
 	nres = 5

	massrange = 0.

 	snapnum = 63

	mass = 'mstars'
	Zrho = np.zeros((3,3,20))
	rhox = np.zeros(20)

	for ir in range(nreheat) :
		for ie in range(neject) :
			Dir = 'millennium_W' + str(ir + 1) + str(ie + 1) + '/'
			mname = 'model_z'
			zz = "%0.3f" % redshifts[snapnum]
			fin_base = Dir + mname + zz
			G = res.read_gals(fin_base, 0, 7)
			res.rho_vs_Mvir_hist(snapnum, mass, ir, ie)
			print
			print
			print
			
	res.hist_ZZ(nreheat, neject, snapnum)
	res.rho_hist_Z(nreheat, neject, snapnum)
	res.rho_hist_ZZ(nreheat, neject, snapnum)
	res.rho_hist_Z_cumulative(nreheat, neject, snapnum)
	res.mass_hist_Z(nreheat, neject, snapnum)
	res.mass_hist_ZZ(nreheat, neject, snapnum)
	res.mass_hist_Z_cumulative(nreheat, neject, snapnum)
	res.StellarMassFunction(nreheat, neject, snapnum)

