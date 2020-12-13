#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 17:04:55 2020

@author: karenchen-wiegart
"""

import numpy as np
import matplotlib.pyplot as plt
from skimage import io
import h5py
from scipy.signal import find_peaks
# import os

# %matplotlib qt
# %matplotlib notebook
# %matplotlib inline
# plt.close('all')

'''
Read fly scan files from h5 format
'''
fly_dir = '/media/karenchen-wiegart/hard_disk/65344_65391_CNT/'
recon_dir = '/media/karenchen-wiegart/hard_disk/65344_65391_CNT/binning/'
# recon_dir = '/home/karenchen-wiegart/ChenWiegartgroup/CHLin/20200818_FXI/65286_65339_SuperP/recon_FBP_CUDA_parzen_bin2/'
scan_start = 65344
scan_end = 65391
fn = fly_dir + f'fly_scan_id_{scan_start}.h5'
h5 = h5py.File(fn, 'r')
# print(list(h5.keys()))
first_prj = np.array(h5['img_tomo'][0])
dark = np.array(h5['img_dark_avg'])[0]
bkg = np.array(h5['img_bkg_avg'])[0]
first_prj = (first_prj-dark)/(bkg-dark)  # Normalization
h5.close()
# plt.figure()
# plt.imshow(first_prj)
# out_file = f'first_prj_fly_{scan_start}.tiff'
# io.imsave(in_dir+out_file, np.float32(first_prj))

# Check the x, y, z coordinates of mosaic scans
x_stage = []
y_stage = []
z_stage = []
i=0
for scan_ID in range(scan_start, scan_end+1):
    fn = fly_dir + f'fly_scan_id_{scan_ID}.h5'
    h5 = h5py.File(fn, 'r')    
    # x = np.float32(h5['x_ini'])
    # y = np.float32(h5['y_ini'])
    # z = np.float32(h5['z_ini'])
    x = np.round(np.float32(h5['x_ini']))
    y = np.round(np.float32(h5['y_ini']))
    z = np.round(np.float32(h5['z_ini']))
    x_stage.append(x)
    y_stage.append(y)
    z_stage.append(z)
    h5.close()
    print([scan_ID, i, x, y, z])
    # print([i, x, y, z])
    i=i+1

# x_stage[5] = x_stage[2]
# x_stage[16] = x_stage[13]
# x_stage[25] = x_stage[22]

#use reconstruction coordination: x-y horizontal, z vertical
#note that in projection coordinates: x-horizontal in project, y-verticl in projection, z-along the depth of focus


#additional binning from the down-sampling
add_bin = 4 #total bin = 2*add_bin
h_FOV = int(1280/add_bin)
v_FOV = int(1080/add_bin)
h_mosaicstep = 750 # bin2 (mosaic step size is the center-to-center distance)
v_mosaicstep = 750 # bin2 (750 for CB, 125 CNT )
d_FOV = int(1280/add_bin)
depth = 125

h_size = int(h_mosaicstep/add_bin)
v_size = int(v_mosaicstep/add_bin)
d_size = int(depth/add_bin)

#define the horizontal and vertical size of the cropped array from image
v_cen =  int(v_FOV/2)
v_w = int(v_size/2)
h_cen =  int(h_FOV/2)
h_w = int(h_size/2)
d_cen = int(d_FOV/2)
d_w = int(d_size/2)

# flip_z = True

#x_stage > x in recon
#y_stage > z in recon
#z_stage > y in recon

x_shift1 = 0
x_shift2 = 0
y_shift1 = 0
y_shift2 = 0
x_pc_shift = 4

pixel_size = 40 #nm
x_dim = int((max(x_stage)-min(x_stage))*1000/pixel_size/add_bin+h_size) + (x_shift1+x_shift2)*3 + x_pc_shift*2
y_dim = int((max(z_stage)-min(z_stage))*1000/pixel_size/add_bin+h_size) + (y_shift1+y_shift2)*3
z_dim = int((max(y_stage)-min(y_stage))*1000/pixel_size/add_bin+v_size)

stitch=np.zeros([z_dim, y_dim, x_dim])
#stitch_all = np.zeros([z_size*v_FOV, y_size*h_FOV, x_size*h_FOV])

h_motorstep = 30 # in um
v_motorstep = 30 # in um


single_dir = recon_dir + 'z[100]/'
i=0
for scan_ID in range(scan_start, scan_end+1):
    fn = recon_dir + f'recon_scan_{scan_ID}_bin{add_bin}.h5'
    h5 = h5py.File(fn, 'r')
    # print(scan_ID)

    # img = io.imread(fn)
    img = np.array(h5['img'])
    img = np.flip(img, axis=1)
    # img = np.flip(img, axis=0)
    
    ### Define pc = pixel center ###
    # x_pc = int((x_stage[i]-min(x_stage))*1000/pixel_size/add_bin)
    # y_pc = int((z_stage[i]-min(z_stage))*1000/pixel_size/add_bin)
    # z_pc = int((y_stage[i]-min(y_stage))*1000/pixel_size/add_bin)  
    # stitch[z_pc : z_pc+2*v_w+1, y_pc : y_pc+2*h_w+1, x_pc : (x_pc+2*h_w+1)] = img[v_cen-v_w:v_cen+v_w+1, h_cen-h_w:h_cen+h_w+1, h_cen-h_w:h_cen+h_w+1]
    # stitch[z_pc : z_pc+2*v_w, y_pc : y_pc+2*h_w, x_pc : (x_pc+2*h_w)] = img[v_cen-v_w:v_cen+v_w, h_cen-h_w:h_cen+h_w, h_cen-h_w:h_cen+h_w]
    
    img_slice = img[v_cen-v_w:v_cen+v_w, h_cen-h_w-y_shift1:h_cen+h_w+y_shift2, h_cen-h_w-x_shift1:h_cen+h_w+x_shift2]
    x_slice = img_slice.shape[2]
    y_slice = img_slice.shape[1] 
    z_slice = img_slice.shape[0]
    
    x_pc = int((x_stage[i]-min(x_stage))*x_slice/h_motorstep)
    y_pc = int((z_stage[i]-min(z_stage))*y_slice/v_motorstep)
    z_pc = int((y_stage[i]-min(y_stage))*z_slice/h_motorstep)
    
    if y_pc == 0: x_pc = x_pc + x_pc_shift*2
    elif (0 < y_pc and y_pc < 300): x_pc = x_pc + x_pc_shift
    
    y_end = y_pc + 2*h_w + y_shift1 + y_shift2
    z_end = z_pc + 2*v_w
    x_end = x_pc + 2*h_w + x_shift1 + x_shift2
    
    # stitch[z_pc : z_pc+2*v_w, y_pc : y_pc+2*h_w, x_pc : (x_pc+2*h_w)] = img_slice
    stitch[z_pc : z_end, y_pc : y_end, x_pc : x_end] = img_slice
    
    # stitch_slice = stitch[z_pc : z_pc+2*v_w, y_pc : y_pc+2*h_w, x_pc : (x_pc+2*h_w)]
    print (f'{scan_ID}, {i}, {x_pc}+{x_end}, {y_pc}+{y_end}, {z_pc}+{z_end}')
    # print (f'{scan_ID}, {img_slice.shape}, {stitch_slice.shape}')
    i=i+1
    h5.close()
    
    # fout = single_dir + f'{scan_ID}_z[100].tiff'+++++
    # io.imsave(fout, np.float32(img[100, h_cen-h_w:h_cen+h_w, h_cen-h_w:h_cen+h_w]))

fout  = recon_dir + 'stitch_motor_04_pc.tiff' 
io.imsave(fout, np.float32(stitch))


'''
Non-overlapped
'''
# recon_xy = np.zeros([row*col, h_FOV, h_FOV])
# for i in range (row*col):
#     x1 = 0 + (i%row)*h_FOV
#     x2 = x1 + h_FOV
#     y1 = 0 + (int(i/col))*h_FOV
#     y2 = y1 + h_FOV
#     recon_xy[i] = stitch[0][y1:y2, x1:x2]
#     print(f'{x1}:{x2} and {y1}:{y2}')

# fout  = recon_dir + 'recon_xy_z[0].tiff' 
# io.imsave(fout, np.float32(recon_xy))
