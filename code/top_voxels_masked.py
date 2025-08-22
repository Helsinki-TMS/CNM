#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 12:51:32 2022

@author: pamilos1
"""

from nilearn import image
import nibabel as nib
import numpy as np
import os
import sys
import pandas as pd
from nilearn.input_data import NiftiMasker

def max3D(x):
    m=np.max(x)
    ind=np.unravel_index(x.argmax(), x.shape)
    
    return m,ind

def peak_coords(img):
    # find max value with coordinates
    a = np.array(img.dataobj)
    m,ind=max3D(a)
    
    # Find the MNI coordinates of the voxel 
    img_coords=image.coord_transform(ind[0], ind[1], ind[2], img.affine)
    return  img_coords

def divide_img_to_hemispheres(img):
    data = np.array(img.dataobj)
    midline=image.coord_transform(0,0,0, np.linalg.inv(img.affine))
    x0=int(np.round(midline[0]))
    aleft=data.copy()
    aleft[:x0,:,:]=0
    aright=data.copy()
    aright[x0:,:,:]=0
    left_img = nib.Nifti1Image(aleft, img.affine)
    right_img = nib.Nifti1Image(aright, img.affine)
    return left_img,right_img

def get_top100_voxels(data,masker):
    data_sort=np.sort(np.abs(data))
    th = data_sort[0,-100]
    #data[np.abs(data)<th]=0
    data[data<th]=0
    top100img=masker.inverse_transform(data)
    return top100img
    

def peak_coords_hems(img):
    # find max value with coordinates
    a = np.array(img.dataobj)
    affine=img.affine
    midline=image.coord_transform(0,0,0, np.linalg.inv(affine))
    x0=int(np.round(midline[0]))
    aleft=a.copy()
    aleft[:x0,:,:]=0
    aright=a.copy()
    aright[x0:,:,:]=0
    m,indleft=max3D(aleft)
    m,indright=max3D(aright)
    
    # Find the MNI coordinates of the voxel 
    img_coords_left=image.coord_transform(indleft[0], indleft[1], indleft[2], img.affine)
    img_coords_right=image.coord_transform(indright[0], indright[1], indright[2], img.affine)
    return  img_coords_left,img_coords_right
    
def top_voxels(folder):

    subs = os.listdir(folder)
    subs = [s for s in subs if os.path.isdir(folder+s)]  
    #subs=['subject1']
    D=[]
    
    for sub in subs:
        files = os.listdir(folder+sub + '/FCmaps/')
        maskfile = os.listdir(folder+sub + '/Mask/')
        masker = NiftiMasker(folder+sub+ '/Mask/' + maskfile[0])
        result_path = folder + sub  + '/top100/'
        if not os.path.exists(result_path):
            print('Creating directory ' + result_path)
            os.makedirs(result_path)
        D=[]
        for f in files:
            img=nib.load(folder+sub+'/FCmaps/' +f)
            img_left,img_right = divide_img_to_hemispheres(img)
            datal=masker.fit_transform(img_left)
            datar=masker.fit_transform(img_right)
            top100img_left = get_top100_voxels(datal, masker)
            top100img_right = get_top100_voxels(datar, masker)
            fname=os.path.basename(f)
            fname=fname.split('.')[0]
            nib.save(top100img_left,result_path + fname + '_top100_left.nii')
            nib.save(top100img_right,result_path + fname + '_top100_right.nii')
            p_coords_left = peak_coords(top100img_left)
            p_coords_right = peak_coords(top100img_right)
            d=[f,p_coords_left,p_coords_right]
            D.append(d)
    
        df=pd.DataFrame(data=D,columns=['file','peak_coords_left','peak_coords_right'])
        df.to_csv(result_path + 'peak_coords.csv',sep=',',index=False)
        
    return
    
    
if __name__ == '__main__':
    folder = sys.argv[1]
    top_voxels(folder)
    
    