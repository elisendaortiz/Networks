#!/usr/bin/python
# -*- coding: utf-8 -*-
import collections
from networkx import *
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy import stats
from math import *
import powerlaw
import pandas as pd
import sys, traceback
import os
import itertools
from natsort import natsorted, ns
import scipy.ndimage as ndi
import matplotlib.gridspec as gridspec
from pylab import *
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
rc('text.latex', preamble=r'\usepackage{cmbright}')
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]#COMMENT THIS LINE AND THE FONTS WON'T BE IN ITALICS/LATEX STYLE
from matplotlib import cm
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import LogLocator
import matplotlib.style
import matplotlib as mpl
from matplotlib.lines import Line2D
from collections import Counter
#mpl.style.use('classic')


# --- XMGARCE COLORS ----
g_white='#FFFFFF'
g_black='#000000'
g_grey='#DCDCDC'
g_grey1='#999999'
g_lightred='#FB9A99'
g_red='#E41A1C'
g_red1='#B21010'
g_lightgreen='#B2DF8A'
g_green='#00FF00'
g_green1='#9ACD32'
g_green2='#008000'
g_green3='#006400'
g_lightblue='#ADD8E6'
g_cyan='#00FFFF'
g_aquamarine='#7FFFD4'
g_turquoise='#40E0D0'
g_blue='#1E90FF'
g_blue1='#0000FF'
g_blue2='#004C99'
g_yellow='#FFFF00'
g_yellow1='#FFD700'
g_yellow2='#DAA520'
g_lightorange='#FD776F'
g_orange='#FFA500'
g_orange1='#FF7F00'
g_brown='#D2B48C'
g_brown1='#994C00'
g_brown2='#663300'
g_pink='#FF66B2'
g_pink1='#FF1493'
g_magenta='#FF00FF'
g_lightpurple='#CAB2D6'
g_purple='#984EA3'
g_violet='#9400D3'
g_indigo='#7221BC'
g_maroon='#670748'
#----Mathematica colors ----
m_blue='#5E81B5'
m_green='#8FB131'
m_mustard='#E19C24'
m_tile='#EC6235'
l_grey='#CCCCC6'
# --------------------------


# 1 HYPERBOLIC PLOT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def read_coords(file_to_read,net):
#READ nodelabels from coordfile
	nodelabels=[]
	with open('WTW-2013.inf_coord') as f:
		for line in f.readlines():
		    li = line.lstrip()
		    if not li.startswith("#"):
				data = (line.rstrip()).split()
				nodelabels.append(int(data[0]))
	iline=-1
	with open(file_to_read) as f:
		for line in f.readlines():
		    li = line.lstrip()
		    if not li.startswith("#"):
				iline=iline+1
				data = (line.rstrip()).split()
#				node=int(data[0])
				node=nodelabels[iline]
				if node in net.nodes():
					net.node[node]['r'] = float(data[2])
					net.node[node]['theta'] = float(data[1])
					net.node[node]['state']= float(data[3]) #We'll generate a file akin to .inf_coord w states for that specific time we wanna plot


				
#Create list of networks we're gonna work with
path0='./'
networks=['WTW-2013']
name=networks[0]



network_names=[]
arr = os.listdir('./')
for item in arr:       
        if '.edge' in item:
                name=item.replace('.edge','')
                network_names.append(name)
nets=natsorted(network_names, key=lambda y: y.lower())

states_files=[]
arr = os.listdir('./')
for item in arr:       
        if 'states_' in item:
                states_files.append(item)
states_files=natsorted(states_files, key=lambda y: y.lower())

for i in range(len(states_files)):
	name=nets[0]
	coord_file=states_files[i]
	time=coord_file.split('_')[1][5:]

	#READ EGDES----------------------------------------------------------------------------------
	G_full = nx.read_edgelist('./'+name+'.edge', nodetype=int)#data=(('weight',float),)
	Ncomps=number_connected_components(G_full)
	comps=sorted(nx.connected_components(G_full), key=len, reverse=True)
	G=max(nx.connected_component_subgraphs(G_full), key=len)      #To extract the GCC
	G=nx.Graph(G)                                          #To remove parallel edges from GCC (just in case, in fact there shouldn't be any)
	G.remove_edges_from(G.selfloop_edges())                #To remove self loops from GCC

	#READ COORDS AND STATES------------------------------------------------------------------------
	read_coords(coord_file,G)
	r = nx.get_node_attributes(G,'r').values()
	theta = nx.get_node_attributes(G,'theta').values()
	thetas=[]
	rs=[]
	states=[]
	for node in G:
		thetas.append(float(G.node[node]['theta']))
		rs.append(float(G.node[node]['r']))

		state=float(G.node[node]['state'])
		if state<0.0: 
			color=m_blue
		else:
			color=g_red1#m_tile
		states.append(color)
	
	rmax=max(rs)

	#MAKE HYPER PLOT----------------------------------------------------------------------------------
	#	fig, ax = plt.subplots()
	ax = subplot(111, polar=True)
	for (u,v) in G.edges:
	    t1 = G.node[u]['theta']
	    t2 = G.node[v]['theta']
	    r1 = G.node[u]['r']
	    r2 = G.node[v]['r']
	    ax.plot([t1,t2], [r1,r2], c='black', linewidth=0.2,zorder=2,alpha=0.1)
	#area=100000/(r**3)
	ax.scatter(thetas, rs, c=states, clip_on = False, edgecolors='black', s=120,zorder=3,alpha=0.6)#s=area # cmap=cm.hsv,m_blue
	ax.set_alpha(1.0)
	ax.plot(np.linspace(0, 2*np.pi, 100), np.ones(100)*rmax, color='black', linestyle='-', zorder=1)
	ax.set_rmax(rmax)
	ax.set_title('t='+time)
	plt.axis('off')
	plt.savefig('./'+name+'_t='+time+'.pdf',bbox_inches='tight')#,bbox_inches='tight')
	#plt.show()
	plt.close()






































