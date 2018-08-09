#slope of 2 represents Larsons's law - constant column density. slope of 3 represents constant volume density.

import numpy as np
from numpy import polyfit
import matplotlib.pyplot as plt
plt.rc('font',family='serif')
import os
#####

fd_dir_cloud_d_smooth = './../Output_data_products/fractal_dimension/cloud_d_smooth/'
fd_dir_cloud_ef_smooth = './../Output_data_products/fractal_dimension/cloud_ef_smooth/'

fd_dir_cloud_d_smooth_list = os.listdir(fd_dir_cloud_d_smooth)
fd_dir_cloud_ef_smooth_list = os.listdir(fd_dir_cloud_ef_smooth)

fd_dir_cloud_d_smooth_list_txt = [s for s in fd_dir_cloud_d_smooth_list if ".dat" in s]
fd_dir_cloud_ef_smooth_list_txt = [s for s in fd_dir_cloud_ef_smooth_list if ".dat" in s]

for file in fd_dir_cloud_d_smooth_list_txt:
    
    raw_data = np.loadtxt(fd_dir_cloud_d_smooth+file)
    
    log_data = np.log10(raw_data)
    
    fit = polyfit(log_data[:,1], log_data[:,2], 1)
    area_xaxis = np.logspace(-5,5)
    fit_plot = 10**fit[1] * area_xaxis**fit[0]
    D_p = fit[0] * 2.
    
    
    fig = plt.figure(figsize=(4, 4))
    ax1 = plt.subplot2grid(shape=(1,1), loc=(0,0), colspan=1)
    ax1.scatter(raw_data[:,1], raw_data[:,2])
    ax1.plot(area_xaxis, fit_plot)
    ax1.text(0.05,0.05, r'D$_{\rm P}$ = '+str(round(D_p,3)), transform = ax1.transAxes)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"Area (pc$^{2}$)")
    ax1.set_ylabel(r"Perimeter (pc)")
    #ax1.set_xlim(1e-3, 1e3)
    #ax1.set_ylim(1e-3, 1e3)

    plt.tight_layout()
    print 'making file: '+fd_dir_cloud_d_smooth+file[:-4]+".pdf"
    plt.savefig(fd_dir_cloud_d_smooth+file[:-4]+".pdf")
    plt.close

for file in fd_dir_cloud_ef_smooth_list_txt:
    
    raw_data = np.loadtxt(fd_dir_cloud_ef_smooth+file)
    log_data = np.log10(raw_data)

    fit = polyfit(log_data[:,0], log_data[:,1], 1)
    area_xaxis = np.logspace(-5,5)
    fit_plot = 10**fit[1] * area_xaxis**fit[0]
    D_p = fit[0] * 2.
    
    
    fig = plt.figure(figsize=(4, 4))
    ax1 = plt.subplot2grid(shape=(1,1), loc=(0,0), colspan=1)
    ax1.scatter(raw_data[:,1], raw_data[:,2])
    ax1.plot(area_xaxis, fit_plot)
    #ax1.text(0.05,0.05, r'Cloud D, D$_{\rm P}$ = '+str(round(D_p_cloud_d,3)), transform = ax1.transAxes)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"Area (pc$^{2}$)")
    ax1.set_ylabel(r"Perimeter (pc)")
    ax1.set_xlim(1e-2, 1e2)
    ax1.set_ylim(1e-2, 1e2)

    plt.tight_layout()
    print 'making file: '+fd_dir_cloud_ef_smooth+file[:-4]+".pdf"
    plt.savefig(fd_dir_cloud_ef_smooth+file[:-4]+".pdf")
    plt.close()