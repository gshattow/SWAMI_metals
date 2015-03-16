#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import numpy as np
import pylab as plt
from random import sample, seed, gauss
from os.path import getsize as getFileSize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator, MaxNLocator
from matplotlib.patches import Circle, PathPatch, Ellipse
import cPickle
from scipy import ndimage
#import mpl_toolkits.axisartist as AA
font = {'family':'serif','size':20, 'serif':'Times New Roman'}

stars_color = ['#006600', '#009900', '#00ff00']
IGMgas_color = ['#000099', '#0066cc', '#66ccff']
ZZ_color = ['#330066', '#9933cc', '#cc99ff'] # '#660099',
Z_color = ['#660033', '#cc0099', '#ff99ff'] #990066

WHIM_color = ['#990000', '#cc3300', '#ff6633']
DM_color = ['#000000', '#666666', '#cccccc']


# stars_color = ['#006600', '#009900', '#00ff00']
# ZZ_color = ['#330066', '#660099', '#9999ff']
# IGMgas_color = ['#000099', '#0033ff', '#66ccff']
# Z_color = ['#660066', '#990066', '#ff00ff'] 
# WHIM_color = ['#000000', '#666666', '#cccccc']

line_cycle = ['-', '--']

colors = [stars_color, IGMgas_color, ZZ_color, Z_color, WHIM_color]
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


Vmax = np.linspace(30, 1500, 10000)
Mvir = 10**(2.81594306*np.log10(Vmax) - 4.5037275 ) * 1.0e10
blanks = np.ones(len(Vmax))
Vvir = (Mvir/3.32e5)**(1./3.)
#Mvir = Vvir**3*3.32e5
e_disk_1 = 3.*blanks
e_disk_2 = 150./Vmax  #4./3. * .3 * (501.4 * 501.4) / (Vmax * Vmax)
e_disk_3 = 10**(np.log10(Vmax)*(-2.24387532) + 5.52911963) #6.5 * (0.5 + pow(Vmax/70., -3.5))
#e_disk_3 = np.where(Vmax < 266., e_disk_3, 10**(np.log10(266.)*(-2.24387532) + 5.52911963))

print Vmax/Vvir

e_halo_1 = 0.3 * (680./Vvir)**2 #blanks
e_halo_2 = 0.3 * (3. * Vmax/Vvir)**2   #0.1 * 4./3. * .3 * (501.4 * 501.4) / (Vmax * Vmax)
e_halo_3a = e_disk_1*0.1*(680./Vvir)**2 #0.32 * (0.5 + pow(Vmax/70., -3.5)) * (680./Vmax)**2
e_halo_3b = e_disk_2*0.1*(680./Vvir)**2 
e_halo_3c = e_disk_3*0.1*(680./Vvir)**2 

e_halo_2 = np.where(Vmax < 266., e_halo_2, 0.3 * (3. * 266./Vvir)**2)

fig = plt.figure(figsize=(12.,10))							# square box
fig.subplots_adjust(top = 0.85, bottom = .15, left = .15, right = .95)
#ax = host_subplot(111, axes_class=AA.Axes)
ax = fig.add_subplot(1,1,1)

ax.set_xscale('log')
ax.set_yscale('log')
plt.axis([30,1500,0.02,200])

for tick in ax.xaxis.get_major_ticks():
	tick.label1.set_fontsize(30)
for tick in ax.yaxis.get_major_ticks():
	tick.label1.set_fontsize(30)

ax.plot(Vmax, e_disk_1, lw = 3, c = IGMgas_color[0], label = r'$\epsilon_{disk}$(R1)')
ax.plot(Vmax, e_disk_2, lw = 3, c = IGMgas_color[1], label = r'$\epsilon_{disk}$(R2)')
ax.plot(Vmax, e_disk_3, lw = 3, c = IGMgas_color[2], label = r'$\epsilon_{disk}$(R3)')


ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0, ncol = 1)


ax.set_xlabel(r'V$_{max}$ (km/s)', fontsize=34)  #36 for stripe
ax.set_ylabel(r'$\epsilon_{disk} = \Delta $M$_{reheat}/\Delta$ M$_{\star}$', fontsize=34) #54 for square

ax2 = ax.twiny() # ax2 is responsible for "top" axis and "right" axis
for tick in ax2.xaxis.get_major_ticks():
	tick.label2.set_fontsize(30)

ax2.set_xscale('log')
plt.axis([min(Mvir),max(Mvir),0.02,200])
ax2.set_xlabel(r'M$_{vir}$ [$h^{-1}$ M$_{\odot}$]' + "\n", fontsize=34)  #36 for stripe


outputFile = 'plots/e_disk.eps'
plt.savefig(outputFile)  # Save the figure
print 'Saved file to', outputFile
plt.close()

fig = plt.figure(figsize=(12.,10))							# square box
fig.subplots_adjust(top = 0.85, bottom = .15, left = .15, right = .95)
#ax = host_subplot(111, axes_class=AA.Axes)
ax = fig.add_subplot(1,1,1)

ax.set_xscale('log')
ax.set_yscale('log')
plt.axis([30,1500,0.02,200])

for tick in ax.xaxis.get_major_ticks():
	tick.label1.set_fontsize(30)
for tick in ax.yaxis.get_major_ticks():
	tick.label1.set_fontsize(30)

width = np.log10(max(Vmax)) - np.log10(min(Vmax))
len150 = (np.log10(150) - np.log10(min(Vmax)))/width
height5 = (np.log10(5) - np.log10(0.02))/4.

ax.plot(Vmax, e_halo_1, lw = 3, c = ZZ_color[1], label = r'$\epsilon_{halo}$(E1)', ls = '-')
ax.plot(Vmax, e_halo_2, lw = 3, c = ZZ_color[1], label = r'$\epsilon_{halo}$(E2)', ls = '--')
#ax.axhline(y = 5, xmin = 0, xmax = len150, lw = 3, c = ZZ_color[1], label = r'$\epsilon_{halo}$(R2)')
#ax.axvline(x = 150., ymin = 0, ymax = height5, lw = 3, c = ZZ_color[1])
ax.plot(Vmax, e_halo_3c, lw = 3, c = ZZ_color[1], label = r'$\epsilon_{halo}$(E3)', ls = ':')

ax.legend(prop = matplotlib.font_manager.FontProperties(size=25),fancybox=True,loc=0, ncol = 1)


ax.set_xlabel(r'V$_{max}$ (km/s)', fontsize=34)  #36 for stripe
ax.set_ylabel(r'$\epsilon_{halo}$', fontsize=34) #54 for square

ax2 = ax.twiny() # ax2 is responsible for "top" axis and "right" axis
for tick in ax2.xaxis.get_major_ticks():
	tick.label2.set_fontsize(30)

ax2.plot(Mvir, e_halo_1, lw = 3, c = ZZ_color[1])

ax2.set_xscale('log')
plt.axis([min(Mvir),max(Mvir),0.02,200])
ax2.set_xlabel(r'M$_{vir}$ [$h^{-1}$ M$_{\odot}$]' + "\n", fontsize=34)  #36 for stripe

#ax2.set_xticks([0., .5*np.pi, np.pi, 1.5*np.pi, 2*np.pi])
#ax2.set_xticklabels(["$0$", r"$\frac{1}{2}\pi$",
 #                    r"$\pi$", r"$\frac{3}{2}\pi$", r"$2\pi$"])

#ax2.axis["right"].major_ticklabels.set_visible(False)

outputFile = 'plots/e_halo.eps'
plt.savefig(outputFile)  # Save the figure
print 'Saved file to', outputFile
plt.close()


fig = plt.figure(figsize=(12.,10))							# square box
fig.subplots_adjust(top = 0.85, bottom = .15, left = .15, right = .95)
#ax = host_subplot(111, axes_class=AA.Axes)
ax = fig.add_subplot(1,1,1)

ax.set_xscale('log')
ax.set_yscale('log')
plt.axis([30,1500,0.02,200])

for tick in ax.xaxis.get_major_ticks():
	tick.label1.set_fontsize(30)
for tick in ax.yaxis.get_major_ticks():
	tick.label1.set_fontsize(30)

width = np.log10(max(Vmax)) - np.log10(min(Vmax))
len150 = (np.log10(150) - np.log10(min(Vmax)))/width
height5 = (np.log10(5) - np.log10(0.02))/4.

ax.plot(Vmax, e_halo_1 - e_disk_1, lw = 3, c = Z_color[0], label = r'$\epsilon_{halo} - \epsilon_{disk}$(R1E1)', ls = '-')
ax.plot(Vmax, e_halo_2 - e_disk_1, lw = 3, c = Z_color[0], label = r'$\epsilon_{halo} - \epsilon_{disk}$(R1E2)', ls = '--')
#ax.plot(Vmax, e_halo_3a - e_disk_1, lw = 3, c = Z_color[0], label = r'$\epsilon_{halo} - \epsilon_{disk}$(R1E3)', ls = ':')

ax.plot(Vmax, e_halo_1 - e_disk_2, lw = 3, c = Z_color[1], label = r'$\epsilon_{halo} - \epsilon_{disk}$(R2E1)', ls = '-')
#ax.plot(Vmax, e_halo_2 - e_disk_2, lw = 3, c = Z_color[1], label = r'$\epsilon_{halo} - \epsilon_{disk}$(R2E2)', ls = '--')
ax.plot(Vmax, e_halo_3b - e_disk_2, lw = 3, c = Z_color[1], label = r'$\epsilon_{halo} - \epsilon_{disk}$(R2E3)', ls = ':')

ax.plot(Vmax, e_halo_1 - e_disk_3, lw = 3, c = Z_color[2], label = r'$\epsilon_{halo} - \epsilon_{disk}$(R3E1)', ls = '-')
ax.plot(Vmax, e_halo_2 - e_disk_3, lw = 3, c = Z_color[2], label = r'$\epsilon_{halo} - \epsilon_{disk}$(R3E2)', ls = '--')
ax.plot(Vmax, e_halo_3c - e_disk_3, lw = 3, c = Z_color[2], label = r'$\epsilon_{halo} - \epsilon_{disk}$(R3E3)', ls = ':')

ax.legend(prop = matplotlib.font_manager.FontProperties(size=20),fancybox=True,loc=0, ncol = 1)


ax.set_xlabel(r'V$_{max}$ (km/s)', fontsize=34)  #36 for stripe
ax.set_ylabel(r'$\epsilon_{halo} - \epsilon_{disk}$', fontsize=34) #54 for square

ax2 = ax.twiny() # ax2 is responsible for "top" axis and "right" axis
for tick in ax2.xaxis.get_major_ticks():
	tick.label2.set_fontsize(30)

ax2.plot(Mvir, e_halo_1 - e_disk_1, lw = 3, c = Z_color[0])

ax2.set_xscale('log')
plt.axis([min(Mvir),max(Mvir),0.02,200])
ax2.set_xlabel(r'M$_{vir}$ [$h^{-1}$ M$_{\odot}$]' + "\n", fontsize=34)  #36 for stripe

#ax2.set_xticks([0., .5*np.pi, np.pi, 1.5*np.pi, 2*np.pi])
#ax2.set_xticklabels(["$0$", r"$\frac{1}{2}\pi$",
 #                    r"$\pi$", r"$\frac{3}{2}\pi$", r"$2\pi$"])

#ax2.axis["right"].major_ticklabels.set_visible(False)

outputFile = 'plots/e_halo-e_disk.eps'
plt.savefig(outputFile)  # Save the figure
print 'Saved file to', outputFile
plt.close()



fig = plt.figure(figsize=(12.,13))							# square box
fig.subplots_adjust(top = 0.88, bottom = .09, left = .15, right = .95, hspace = 0.)

ax = fig.add_subplot(3,1,1)

ax.set_xscale('log')
ax.set_yscale('log')
plt.axis([30,1500,0.09,90])

for tick in ax.xaxis.get_major_ticks():
	tick.label1.set_fontsize(30)
for tick in ax.yaxis.get_major_ticks():
	tick.label1.set_fontsize(30)

ax.plot(Vmax, e_disk_1, lw = 3, c = IGMgas_color[0], label = r'$\epsilon_{disk}$(R1)')
ax.plot(Vmax, e_disk_2, lw = 3, c = IGMgas_color[1], label = r'$\epsilon_{disk}$(R2)')
ax.plot(Vmax, e_disk_3, lw = 3, c = IGMgas_color[2], label = r'$\epsilon_{disk}$(R3)')


ax.legend(prop = matplotlib.font_manager.FontProperties(size=18),fancybox=True,loc=1, ncol = 1)
ax.set_ylabel(r'$\epsilon_{disk}$', fontsize=34) #54 for square


ax2 = ax.twiny() # ax2 is responsible for "top" axis and "right" axis
for tick in ax2.xaxis.get_major_ticks():
	tick.label2.set_fontsize(30)

ax2.plot(Mvir, e_disk_1, lw = 3, c = IGMgas_color[0])

ax2.set_xscale('log')
plt.axis([min(Mvir),max(Mvir),0.09,90])
ax2.set_xlabel(r'M$_{vir}$ [$h^{-1}$ M$_{\odot}$]' + "\n", fontsize=34)  #36 for stripe

ax = fig.add_subplot(3,1,2)

ax.set_xscale('log')
ax.set_yscale('log')
plt.axis([30,1500,0.09,90])

for tick in ax.xaxis.get_major_ticks():
	tick.label1.set_fontsize(30)
for tick in ax.yaxis.get_major_ticks():
	tick.label1.set_fontsize(30)

ax.plot(Vmax, e_halo_1, lw = 3, c = ZZ_color[1], label = r'$\epsilon_{halo}$(E1)', ls = '-')
ax.plot(Vmax, e_halo_2, lw = 3, c = ZZ_color[1], label = r'$\epsilon_{halo}$(E2)', ls = '--')
ax.plot(Vmax, e_halo_3c, lw = 3, c = ZZ_color[1], label = r'$\epsilon_{halo}$(E3)', ls = ':')


ax.legend(prop = matplotlib.font_manager.FontProperties(size=18),fancybox=True,loc=1, ncol = 1)
ax.set_ylabel(r'$\epsilon_{halo}$', fontsize=34) #54 for square
ax.set_xlabel(r'V$_{max}$ (km/s)', fontsize=34)  #36 for stripe

ax = fig.add_subplot(3,1,3)

ax.set_xscale('log')
ax.set_yscale('log')
plt.axis([30,1500,0.09,90])

for tick in ax.xaxis.get_major_ticks():
	tick.label1.set_fontsize(30)
for tick in ax.yaxis.get_major_ticks():
	tick.label1.set_fontsize(30)

ax.plot(Vmax, e_halo_1 - e_disk_1, lw = 3, c = Z_color[0], label = r'R1E1', ls = '-')
ax.plot(Vmax, e_halo_2 - e_disk_1, lw = 3, c = Z_color[0], label = r'R1E2', ls = '--')
ax.plot(Vmax, e_halo_3a - e_disk_1, lw = 3, c = Z_color[0], label = r'R1E3', ls = ':')

ax.plot(Vmax, e_halo_1 - e_disk_2, lw = 3, c = Z_color[1], label = r'R2E1', ls = '-')
ax.plot(Vmax, e_halo_2 - e_disk_2, lw = 3, c = Z_color[1], label = r'R2E2', ls = '--')
ax.plot(Vmax, e_halo_3b - e_disk_2, lw = 3, c = Z_color[1], label = r'R2E3', ls = ':')

ax.plot(Vmax, e_halo_1 - e_disk_3, lw = 3, c = Z_color[2], label = r'R3E1', ls = '-')
ax.plot(Vmax, e_halo_2 - e_disk_3, lw = 3, c = Z_color[2], label = r'R3E2', ls = '--')
ax.plot(Vmax, e_halo_3c - e_disk_3, lw = 3, c = Z_color[2], label = r'R3E3', ls = ':')
# ax.plot(Vmax, e_halo_1 - e_disk_1, lw = 3, c = Z_color[0], label = r'$\epsilon_{halo} - \epsilon_{disk}$(R1E1)', ls = '-')
# ax.plot(Vmax, e_halo_2 - e_disk_2, lw = 3, c = Z_color[1], label = r'$\epsilon_{halo} - \epsilon_{disk}$(R2E2)', ls = '--')
# ax.plot(Vmax, e_halo_3c - e_disk_3, lw = 3, c = Z_color[2], label = r'$\epsilon_{halo} - \epsilon_{disk}$(R3E3)', ls = ':')


ax.legend(prop = matplotlib.font_manager.FontProperties(size=18),fancybox=True,loc=1, ncol = 3)
ax.set_ylabel(r'$\epsilon_{halo} - \epsilon_{disk}$', fontsize=34) #54 for square
ax.set_xlabel(r'V$_{max}$ (km/s)', fontsize=34)  #36 for stripe

outputFile = 'plots/e_halo+e_disk.eps'
plt.savefig(outputFile)  # Save the figure
print 'Saved file to', outputFile
plt.close()





fig = plt.figure(figsize=(12.,13))							# square box
fig.subplots_adjust(top = 0.88, bottom = .09, left = .15, right = .95, hspace = 0.)

ax = fig.add_subplot(3,1,1)

ax.set_xscale('log')
ax.set_yscale('log')
plt.axis([30,1500,0.09,90])

for tick in ax.xaxis.get_major_ticks():
	tick.label1.set_fontsize(30)
for tick in ax.yaxis.get_major_ticks():
	tick.label1.set_fontsize(30)

ax.plot(Vmax, e_halo_1 - e_disk_1, lw = 3, c = IGMgas_color[0], label = r'R1E1', ls = '-')
ax.plot(Vmax, e_halo_1 - e_disk_2, lw = 3, c = IGMgas_color[1], label = r'R2E1', ls = '-')
ax.plot(Vmax, e_halo_1 - e_disk_3, lw = 3, c = IGMgas_color[2], label = r'R3E1', ls = '-')


ax.legend(prop = matplotlib.font_manager.FontProperties(size=18),fancybox=True,loc=1, ncol = 1)
#ax.set_ylabel(r'$\Delta m_{eject}/\Delta m_{*} = \epsilon_{halo} - \epsilon_{disk}$', fontsize=34) #54 for square


ax2 = ax.twiny() # ax2 is responsible for "top" axis and "right" axis
for tick in ax2.xaxis.get_major_ticks():
	tick.label2.set_fontsize(30)

ax.plot(Mvir, e_halo_1 - e_disk_1, lw = 3, c = Z_color[0], ls = '-')

ax2.set_xscale('log')
plt.axis([min(Mvir),max(Mvir),0.09,90])
ax2.set_xlabel(r'M$_{\mathrm{vir}}$ ($h^{-1} \mathrm{M}_{\odot}$)' + "\n", fontsize=34)  #36 for stripe

ax = fig.add_subplot(3,1,2)

ax.set_xscale('log')
ax.set_yscale('log')
plt.axis([30,1500,0.09,90])

for tick in ax.xaxis.get_major_ticks():
	tick.label1.set_fontsize(30)
for tick in ax.yaxis.get_major_ticks():
	tick.label1.set_fontsize(30)

ax.plot(Vmax, e_halo_2 - e_disk_1, lw = 3, c = IGMgas_color[0], label = r'R1E2', ls = '--')
ax.plot(Vmax, e_halo_2 - e_disk_2, lw = 3, c = IGMgas_color[1], label = r'R2E2', ls = '--')
ax.plot(Vmax, e_halo_2 - e_disk_3, lw = 3, c = IGMgas_color[2], label = r'R3E2', ls = '--')



ax.legend(prop = matplotlib.font_manager.FontProperties(size=18),fancybox=True,loc=1, ncol = 1)
ax.set_ylabel(r'$\Delta \mathrm{m}_{\mathrm{eject}}/\Delta \mathrm{m}_{*} = \epsilon_{\mathrm{halo}} - \epsilon_{\mathrm{disk}}$', fontsize=34) #54 for square

ax = fig.add_subplot(3,1,3)

ax.set_xscale('log')
ax.set_yscale('log')
plt.axis([30,1500,0.09,90])

for tick in ax.xaxis.get_major_ticks():
	tick.label1.set_fontsize(30)
for tick in ax.yaxis.get_major_ticks():
	tick.label1.set_fontsize(30)


ax.plot(Vmax, e_halo_3a - e_disk_1, lw = 5, c = IGMgas_color[0], label = r'R1E3', ls = ':')
ax.plot(Vmax, e_halo_3b - e_disk_2, lw = 5, c = IGMgas_color[1], label = r'R2E3', ls = ':')
ax.plot(Vmax, e_halo_3c - e_disk_3, lw = 5, c = IGMgas_color[2], label = r'R3E3', ls = ':')


ax.legend(prop = matplotlib.font_manager.FontProperties(size=18),fancybox=True,loc=1, ncol = 1)
#ax.set_ylabel(r'$\Delta m_{eject}/\Delta m_{*} = \epsilon_{halo} - \epsilon_{disk}$', fontsize=34) #54 for square
ax.set_xlabel(r'V$_{\mathrm{max}}$ ($\mathrm{km}\ \mathrm{s}^{-1}$)', fontsize=34)  #36 for stripe

outputFile = 'plots/e_halo-e_disk3.eps'
plt.savefig(outputFile)  # Save the figure
print 'Saved file to', outputFile
plt.close()


fig = plt.figure(figsize=(12.,10))							# square box
fig.subplots_adjust(top = 0.85, bottom = .15, left = .15, right = .95)
#ax = host_subplot(111, axes_class=AA.Axes)
ax = fig.add_subplot(1,1,1)

#plt.axis([-1.85,-1.4,0.25,1.6])
plt.axis([-1.46,-1.29,0.25,1.6])

for tick in ax.xaxis.get_major_ticks():
	tick.label1.set_fontsize(30)
for tick in ax.yaxis.get_major_ticks():
	tick.label1.set_fontsize(30)

medZ = [[-1.498, -1.799, -1.498], [-1.486, -1.436, -1.637], [-1.835, -1.535, -1.666]]
NmedZ = [[-1.319, -1.372, -1.319], [-1.310, -1.300, -1.448], [-1.441, -1.341, -1.393]]
aveZ = [[-1.197, -1.243, -1.197], [-1.200, -1.180, -1.278], [-1.282, -1.239, -1.266]]
maxZ = [[0.850, 0.710, 0.850], [0.411, 0.970, 0.754], [1.565, 0.310, 0.804]]
shapes = ['o', '^', 's']
for ir in range(3) :
	for ie in range(3) :
		label = 'R'+ str(ir + 1) + 'E' + str(ie + 1)
		ax.scatter(NmedZ[ir][ie], maxZ[ir][ie], c = ZZ_color[ir], marker = (5,ie), s = 200, edgecolor = ZZ_color[ir], label = label)


ax.legend(prop = matplotlib.font_manager.FontProperties(size=20),fancybox=True,loc=0, ncol = 3)


ax.set_xlabel(r'median(log(Z/Z$_{\odot}$))', fontsize=34)  #36 for stripe
ax.set_ylabel(r'maximum(log(Z/Z$_{\odot}$))', fontsize=34) #54 for square


outputFile = 'plots/NmedZ_vs_maxZ.eps'
plt.savefig(outputFile)  # Save the figure
print 'Saved file to', outputFile
plt.close()

