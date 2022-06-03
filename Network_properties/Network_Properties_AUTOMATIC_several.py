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
#mpl.style.use('classic')

network_names=[]
arr = os.listdir('./Data')
for item in arr:       
        if '.edges' in item:
                name=item.replace('.edges','')
                network_names.append(name)

nets=natsorted(network_names, key=lambda y: y.lower())

for net in nets:    
	name=net 
	if not os.path.exists('./Figures_properties/'):
		os.makedirs('./Figures_properties/')
	if not os.path.exists('./Edgelists_gcc/'):
		os.makedirs('./Edgelists_gcc/')
	if not os.path.exists('./Files_properties/'):
		os.makedirs('./Files_properties/')

	print '.....',name,'.....'

	# EXTRACT NETWORK AND GCC:
	#WARNING: Edgelists to be read  need to be CLEAN--> No multiple connections (BUT they may have self loops in order to include isolated nodes).
	G_full = nx.read_edgelist('./Data/'+name+'.edges', nodetype=int, data=(('weight',float),))#data=(('weight',float),)
	Ncomps=number_connected_components(G_full)
	comps=sorted(nx.connected_components(G_full), key=len, reverse=True)
	G=max(nx.connected_component_subgraphs(G_full), key=len)      #To extract the GCC
	G=nx.Graph(G)                                          #To remove parallel edges from GCC (just in case, in fact there shouldn't be any)
	G.remove_edges_from(G.selfloop_edges())                #To remove self loops from GCC

	
	#PRODUCE GCC EDGELIST FILE when the full network has more than one component:
	if G_full.number_of_nodes() != G.number_of_nodes():
		print 'Producing edgelist'
		outedgelist = open('./Edgelists_gcc/'+name+'.edge', "w")
		for line in nx.generate_edgelist(G, data=False):
		    print>>outedgelist, line


	#EXTRACT BASIC INFORMATION
	print 'Extracting basic properties'
	K = [val for (node, val) in G.degree()]
	avdegree=sum(K)/float(len(K))
	degree_one=K.count(1)
	clustering = nx.clustering(G).values()
	avclust=sum(clustering)/float(G.number_of_nodes()-degree_one)#This is the clustering considering all nodes except the ones with k=1.
	maxdegree=np.amax(K)  
	f = open('./Files_properties/'+name+'_properties.txt', 'w')
	print >> f, name+"  THIS IS GCC INFORMATION:"
	print >> f, "Nodes full net:              ", G_full.number_of_nodes()
	print >> f, "Nodes gcc:                   ", G.number_of_nodes()
	print >> f, "Egdes:                       ", G.number_of_edges()
	print >> f, "Av. Degree:                  ", avdegree
	print >> f, "Max Degree:                  ", maxdegree
	print >> f, "Mean clustering coefficient: ", avclust
	print >> f, "Av. shortest path length     ", nx.average_shortest_path_length(G)
	print >> f, "Diameter                     ", nx.diameter(G, e=None)
	f.close()

	# REAL NETWORK INFORMATION ----------------------------------------------------------------
	# For each element in the degree view, obtain node label and its degree
	Node_label_degrees = [node for (node, val) in G.degree()]
	k_nodelabels, degree_values = zip(*G.degree())#degrees_values=degrees sorted by nodelabel in ascending order
	######## Calculate CCP(k):
	#Obtain how many nodes there are of each degree and find P(k)
	 #Real
	counter=collections.Counter(K)
	frequency_of_degree_class=np.array(counter.values())
	degree_class=np.array(counter.keys())
	PK=zip(degree_class,frequency_of_degree_class.astype(np.float32)/G.number_of_nodes())
	PK.sort(key=lambda x: x[0])
	Pk=zip(*PK)
	CPK=np.zeros(len(Pk[0])+1)
	CCPK=np.zeros(len(Pk[0]))
	for i in range(1,len(Pk[0])+1):
	    CPK[i]=CPK[i-1]+Pk[1][i-1]
	complementary=1.-CPK
	CCPK=list(complementary)
	CCPK.pop()
	degree_class_sorted=sorted(degree_class)

	######### CLUSTERING ALL NODES AND AUTOMATIC BINNING
	clustering_dict= nx.clustering(G)
	clust_to_sort= clustering_dict.items()
	clust_sorted= sorted(clust_to_sort, key=lambda x: x[0])
	clust_nodelabels, clust_values = zip(*clust_sorted)
	#les dos coses q has de plotar una against d laltre son : (degree_values, clust_values) son aquestes perq lu de adalt ensures que el degree i el clustering son del mateix node.
	binning=np.logspace(log(0.1)/log(10),log(max(degree_values))/log(10),30,endpoint=True)#, base=10 by default, so start= 10**sth, end=10**sthelse
	biins=binning.tolist()
	bin_means, bin_edges, binnumber = stats.binned_statistic(degree_values, clust_values, statistic='mean', bins=biins)
	bin_std, bin_std_edges, bin_std_number = stats.binned_statistic(degree_values, clust_values, statistic='std', bins=biins)
	biins.pop(0)
	biins=np.array(biins) #Doing this is necessary to apply the mask that ignores the "Nan" values of the bins and hence connects all points in the plot with a line
	bin_means=np.array(bin_means)
	bin_means[bin_means == 0] = 'nan' #To mask too the number 0.0 if it appears
	bin_std=np.array(bin_std)
	bin_means_mask = np.isfinite(bin_means)

	######### KNN ALL NODES AND AUTOMATIC BINNING
	knn_dict=nx.average_neighbor_degree(G)
	knn_to_sort= knn_dict.items()
	knn_sorted= sorted(knn_to_sort, key=lambda x: x[0])
	knn_nodelabels, knn_values = zip(*knn_sorted)
	#les dos coses q has de plotar una against d laltre son : (degree_values, knn_values) son aquestes perq lu de adalt ensures que el degree i el clustering son del mateix node.
	binning_knn=np.logspace(log(0.1)/log(10),log(max(degree_values))/log(10),30,endpoint=True)#, base=10 by default, so start= 10**sth, end=10**sthelse
	biins_knn=binning_knn.tolist()
	bin_means_knn, bin_edges_knn, binnumber_knn = stats.binned_statistic(degree_values, knn_values, statistic='mean', bins=biins_knn)
	bin_std_knn, bin_std_edges_knn, bin_std_number_knn = stats.binned_statistic(degree_values, knn_values, statistic='std', bins=biins_knn)
	biins_knn.pop(0)
	biins_knn=np.array(biins_knn) #Doing this is necessary to apply the mask that ignores the "Nan" values of the bins and hence connects all points in the plot with a line
	bin_means_knn=np.array(bin_means_knn)
	bin_means_knn[bin_means_knn == 0] = 'nan' #To mask too the number 0.0 if it appears
	bin_std_knn=np.array(bin_std_knn)
	bin_means_mask_knn = np.isfinite(bin_means_knn)
	#Knn QUICK VERSION: pero no tens tant de marge amb el binning per fer q quedi mes bonic ja q l'output del dictionary ja es per classes de degree.
	#knn_dict=nx.k_nearest_neighbors(G)
	#knn_to_sort= knn_dict.items()
	#knn_sorted= sorted(knn_to_sort, key=lambda x: x[0])#we sort knn dictionary by degree to be able to plot it correctly
	#knn_degree_classes, knn_values = zip(*knn_sorted)


	#PLOT PROPERTIES TOGETHER:#############################################################################################
	#----Mathematica colors ----
	m_blue='#5E81B5'
	m_green='#8FB131'
	m_mustard='#E19C24'
	m_tile='#EC6235'
	l_grey='#CCCCC6'
	# --------------------------

	print 'Plotting data'
	fig, axes = plt.subplots(nrows=1, ncols=4, sharex=False, sharey=False, figsize = (12.96,2.7))#15.64,5.4 12.74,4.4
	gs1 = gridspec.GridSpec(1, 3)
	gs1.update(wspace=0.28, hspace=0.43)

	ax1 = plt.subplot(gs1[0,0])
	ax1.plot(degree_class_sorted,CCPK, marker='o',c=m_blue,linewidth=1.0,clip_on=True)
	ax1.set_xscale("log")#nonposx='clip'
	ax1.set_yscale("log")#nonposy='clip'
	ax1.set_ylabel(r"$\mathrm{CCP}$",fontsize=13)
	ax1.set_xlabel(r"$k$",fontsize=13)
	ax1.xaxis.set_major_locator(LogLocator(base=10))
	ax1.yaxis.set_major_locator(LogLocator(base=10))
	#ax1.legend(loc="lower left", framealpha=0.0, frameon=False, fontsize=13, numpoints=1)
	#plt.title(name, fontsize=13)#name

	ax2 = plt.subplot(gs1[0,1])
	#ax2.scatter(degree_values, clust_values,marker='o', facecolors='none', edgecolors=m_mustard, linewidth=1.0,label='raw data')
	ax2.errorbar(biins[bin_means_mask], bin_means[bin_means_mask], yerr=bin_std[bin_means_mask], c=m_blue, linewidth=1.0, fmt='-o', ecolor='black')
	ax2.set_xscale("log")#, nonposx='clip'
	ax2.set_yscale("log")#, nonposy='clip'
	ax2.set_ylabel(r"$\langle c \rangle$",fontsize=14)
	ax2.set_xlabel(r"$k$",fontsize=13)
	plt.ylim(10**(-3),10**1)
	ax2.xaxis.set_major_locator(LogLocator(base=10))
	ax2.yaxis.set_major_locator(LogLocator(base=10))
	#ax2.legend(loc="lower left", framealpha=0.0, frameon=False, fontsize=13, numpoints=1)
	#plt.title(name, fontsize=13)

	ax3 = plt.subplot(gs1[0,2])
	#ax3.plot(knn_degree_classes,knn_values, c='g', marker='o',linewidth=0,label='Real',clip_on=True)#KNN QUICK VERSION
	#ax3.scatter(degree_values, knn_values,marker='o', facecolors='none', edgecolors=m_mustard, linewidth=1.0,label='raw data')
	ax3.errorbar(biins_knn[bin_means_mask_knn], bin_means_knn[bin_means_mask_knn], yerr=bin_std_knn[bin_means_mask_knn], c=m_blue, linewidth=1.0, fmt='-o', ecolor='black')
	ax3.set_xscale("log")#, nonposy='clip'
	ax3.set_yscale("log")#, nonposy='clip'
	ax3.set_ylabel(r"$\langle k_{\mathsf{nn}}\rangle$",fontsize=14)
	ax3.set_xlabel(r"$k$",fontsize=13)
	plt.ylim(10**0,10**3)
	ax3.xaxis.set_major_locator(LogLocator(base=10))
	ax3.yaxis.set_major_locator(LogLocator(base=10))
	#ax3.legend(loc="best", framealpha=0.0, frameon=False, fontsize=13, numpoints=1)
	#plt.title(name, fontsize=13)

	fig.savefig('./Figures_properties/'+name+'_properties.eps',bbox_inches='tight')#,bbox_inches='tight'
	#plt.show()


########################################################################################################################





