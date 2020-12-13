#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 16:33:52 2020

@author: karenchen-wiegart
"""

import numpy as np
import matplotlib.pyplot as plt
from skimage import io
import h5py
#import timeit

# %matplotlib qt
# %matplotlib notebook
# %matplotlib inline
plt.close('all')

'''
Read fly scan files from h5 format
'''
in_dir = '/NSLS2/xf18id1/users/2020Q3/KAREN_Proposal_305052/'
scan_start = 79253
scan_end = 79444
fn = in_dir + f'fly_scan_id_{scan_start}.h5'
h5 = h5py.File(fn, 'r')
# print(list(h5.keys()))
first_prj = np.array(h5['img_tomo'][0])
dark = np.array(h5['img_dark_avg'])[0]
bkg = np.array(h5['img_bkg_avg'])[0]
first_prj = (first_prj-dark)/(bkg-dark)
h5.close()
plt.figure()
plt.imshow(first_prj)
# out_file = f'first_prj_fly_{scan_start}.tiff'
# io.imsave(in_dir+out_file, np.float32(first_prj))


# Check the x, y, z coordinates of mosaic scans
for scan_ID in range(scan_start, scan_end+1):
    fn = in_dir + f'fly_scan_id_{scan_ID}.h5'
    h5 = h5py.File(fn, 'r')    
    x = np.float32(h5['x_ini'])
    y = np.float32(h5['y_ini'])
    z = np.float32(h5['z_ini'])
    h5.close()
    print([scan_ID, x, y, z])

# Make an empty array for stitch
num_col = 8
num_row = 2
depth = 0
stitch = np.zeros([first_prj.shape[0]*num_row, first_prj.shape[1]*num_col])

# Stitch together
for i in range(num_col*num_row):
    scan_ID = scan_start +num_col*num_row*depth + i
    fn = in_dir + f'fly_scan_id_{scan_ID}.h5'
    h5 = h5py.File(fn, 'r')    
    first_prj = np.array(h5['img_tomo'][0])
    dark = np.array(h5['img_dark_avg'])[0]
    bkg = np.array(h5['img_bkg_avg'])[0]
    first_prj = (first_prj-dark)/(bkg-dark)
    n = int(i/num_col)
    m = i%num_col
    stitch[first_prj.shape[0]*n:first_prj.shape[0]*(n+1), first_prj.shape[1]*m:first_prj.shape[1]*(m+1)] = first_prj
    h5.close()

plt.figure()
plt.imshow(stitch)
# out_file = f'stitch_{scan_start}_{scan_start+num_col*num_row-1}.tiff'
# io.imsave(in_dir+out_file, np.float32(stitch))
#
#
## Stitch XZ plan after recon
#for i in range(num_col*num_row):
#    scan_ID = scan_start + i
#    fn = in_dir + f'recon_scan_{scan_ID}_bin1.h5'
#    h5 = h5py.File(fn, 'r')    
#    xz_prj = np.array(h5['img'])[:,int(first_prj.shape[1]/2),:]
#    n = int(i/num_col)
#    m = i%num_col
#    stitch[xz_prj.shape[0]*n:xz_prj.shape[0]*(n+1), xz_prj.shape[1]*m:xz_prj.shape[1]*(m+1)] = xz_prj
#    h5.close()
#
#plt.figure()
#plt.imshow(stitch)
## out_file = f'stitch_reco_{scan_start}_{scan_start+num_col*num_row-1}.tiff'
## io.imsave(in_dir+out_file, np.float32(stitch))

plt.show()
