#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 10:24:46 2020

@author: pamilos1
"""

import glob
from nilearn import image
import nibabel as nib
from nilearn.masking import apply_mask, unmask
import numpy as np
from scipy import stats
import pandas as pd
import os

def compute_correlations_2(images,j,mask_file):
    values=[]
    mean_vals=[]
    cnames = ['S%d' % x for x in range(1,len(images)+1)]
    for i in images:
        if j==0:
            v = i.get_fdata().flatten()
        elif (j>0 and j<3):
            v = apply_mask(i,mask_file) 
        elif j>2:
            v = apply_mask(i,mask_file) 
        values.append(v)
        if j > 2:
            v = v[abs(v)>0] 
            mean_vals.append(v.mean())
        else:
            mean_vals.append(v.mean())
        
    df=pd.DataFrame(data=np.asarray(values).T,columns=cnames)
    C = df.corr()
    return C, mean_vals

def compute_correlations(images):
    values=[]
    cnames = ['S%d' % x for x in range(1,len(images)+1)]
    for i in images:
        v = i.get_fdata().flatten()
        values.append(v)
    df=pd.DataFrame(data=np.asarray(values).T,columns=cnames)
    C = df.corr()
    return C

def make_tail(include,multiply,abs_val):
    
    if include:
        inc = '_'.join([str(elem) for elem in include])
        inc = '_' + inc
    else: inc = ''
    
    if abs_val:
        mult = '_abs'
    elif multiply:
        mult = '_m'.join([str(elem) for elem in multiply])
        mult = '_m' + mult
    else: mult=''
    
    tail ='%s%s' % (inc,mult)
    return tail

def get_sessions(fn):
    p0 = [x.split('Results')[0] for x in fn]
    l1 = [1 for s in p0 if s[-1]== '/']
    l2 = [int(s[-2]) for s in p0 if s[-1] != '/']
    ses = list(set(l1+l2))
    return ses

def corr_seeds(fn,s,include,multiply,abs_val,mask_file,savepath,tail):
    #correlations between seeds
    R = []
    for r in include:
        fn_roi = [x for x in fn if '/'+r in x]
        img = image.load_img(fn_roi)
        v = apply_mask(img,mask_file) 
        if abs_val:
            v=np.abs(v)
        elif r in multiply:
            v = -v                       
        R.append(v)
        
    mean_img = image.load_img(savepath +'S%d_mean_connectivity_map%s.nii' % (s,tail))
    meanv = apply_mask(mean_img,mask_file)
    meanv = np.reshape(meanv,v.shape)
    R.append(meanv)
    cnames=include.copy()
    cnames.append('avg')
    df=pd.DataFrame(data=np.asarray(np.squeeze(R)).T,columns=cnames)
    C = df.corr()
    return C
    
def corr_ses(fn,ses,include,multiply,abs_val,mask_file):
    #correlations between seeds
    C_roi=[]
    for r in include:
        R = []
        for s in ses:
            if s==1:
                fn_ses =[x for x in fn if '/Results' in x]
            else:
                fn_ses =[x for x in fn if 'S%d_Results' % s in x]           
            
            fn_roi = [x for x in fn_ses if '/'+r in x]
            img = image.load_img(fn_roi)
            v = apply_mask(img,mask_file) 
            if abs_val:
                v=np.abs(v)
            elif r in multiply:
                v = -v                       
            R.append(v)
            
        mean_roi = np.mean(np.asarray(R),axis=0)
        R.append(mean_roi)
        cnames=['%s_ses%d' % (r,x) for x in ses]
        cnames.append('%s_avg' % r)
        df=pd.DataFrame(data=np.asarray(np.squeeze(R)).T,columns=cnames)
        C = df.corr()
        C_roi.append(C)
    return C_roi


def compute_mean_connectivity_maps(dname=None,subj_name='subject',multiply='',include='',abs_val = False,mask_file = '/m/nbe/project/heps/eregmasks/greymask04.nii',nvox=500,th=0):
    #################################################################################################################
    # Input:
    # dname - path to the directory containing the subject folders
    # subj_name - prefix of the subject folders, default 'subject'. 
    # multiply - ROIs to be multiplied by -1, inside brackets and separated by comma, e.g.['ROI1,'ROI2']. By default no images are multiplied.
    # include - ROIs to be included in the average, inside brackets and separated by comma, e.g. ['ROI1','ROI2']. By default all ROIs are included.
    # abs_val - True/False, default False. True = use absolute values when computing the mean.
    # mask_file - path to mask to be used, default '/m/nbe/project/heps/eregmasks/greymask04.nii'
    # nvox - number of voxels to be included when thresholding, e.g. top nvox voxels will be selected. Default 500 voxels.
    # th - threshold to be used in thresholding. Default 0.
    #
    # Output:
    # Average connectivity maps and correlation matrices will be saved in 'jointconnectivity'-folder
    #################################################################################################################
    if dname is None:
        raise ValueError('Path to the folder containing the subject folders must me specified.')
    
    th_str=str(th)
    th_str=th_str.replace('.','_')
    #names of all nii files
    #fn=glob.glob(dname + subj_name + '*/*Results/FC_FunImgRWCFB/*.nii')
    fn=glob.glob(dname + subj_name + '*/*Results/FC_FunImgARCFB/*.nii')
    #include only certain ROIs
    if include:
        fn[:] = [x for x in fn if any('/'+y in x for y in include)]
    
    # list subjects    
    p1 = [sub.split(subj_name)[1] for sub in fn]
    subs = list(set([int(sub.split('/')[0][0]) for sub in p1]))
    
    
    for sub in subs:
        fn_sub_all = [x for x in fn if subj_name+str(sub) in x]
        
        p1 = [sub.split(subj_name)[1] for sub in fn_sub_all]
        subs_list = list(set([sub.split('/')[0] for sub in p1]))
        if len(subs_list[0])>1:
            subjs=list(subs_list[0].split('and'))
            subjs = [int(x) for x in subjs]
        else: subjs = [sub]
        

        F=[]
        if len(subjs)>1:
            for f in subjs:
                ff = [x for x in fn_sub_all if subj_name+'%d.nii' % f in x]
                F.append(ff)
        else:
              F.append(fn_sub_all)
        i=0    
        for fn_sub in F:
            isub=subjs[i]
            i+=1
        
            #list sessions
            ses = get_sessions(fn_sub)            
            if isub<10:
               savepath = dname + 'jointconnectivity/' + subj_name + '0%d/' % isub
            else:
               savepath = dname + 'jointconnectivity/' + subj_name + '%d/' % isub

            #savepath = '/m/home/home5/58/pamilos1/unix/psykoosi/test/' + 'jointconnectivity/' + subj_name + '%d/' % isub
            if not os.path.exists(savepath):
                os.makedirs(savepath)
            # initialize mean connectivity maps over all sessions
            M = []
            M_masked = []
            M_zscored = []
            M_topn = []
            M_th = []
            #loop over sessions
            for s in ses:
                if s==1:
                    fn_sub_ses =[x for x in fn_sub if '/Results' in x]
                else:
                    fn_sub_ses =[x for x in fn_sub if 'S%d_Results' % s in x]
                if abs_val:
                    mean_img = image.mean_img(image.math_img('np.absolute(img)',img=fn_sub_ses))
                elif multiply:
                    fn_neg = [x for x in fn_sub_ses if any('/'+y in x for y in multiply)]
                    neg_mean = image.mean_img(image.math_img("-img", img=fn_neg))
                    fn_pos = [ele for ele in fn_sub_ses if ele not in fn_neg] 
                    pos_mean=image.mean_img(fn_pos)
                    mean_img = image.mean_img([neg_mean,pos_mean])
                else:
                    mean_img = image.mean_img(fn_sub_ses)
                
                tail = make_tail(include,multiply,abs_val)
                
                nib.save(mean_img, savepath +'S%d_mean_connectivity_map%s.nii' % (s,tail))
                M.append(mean_img)
                # mask and zscore
                values_masked = apply_mask(mean_img,mask_file) 
                values_zscored = values_masked #stats.zscore(values_masked)
                mean_img_masked=unmask(values_masked,mask_file)
                nib.save(mean_img_masked, savepath +'S%d_mean_connectivity_map_masked%s.nii' % (s,tail))
                mi_zscore=unmask(values_zscored,mask_file)   
                nib.save(mi_zscore, savepath +'S%d_mean_connectivity_map_masked_zscored%s.nii' % (s,tail))
                M_masked.append(mean_img_masked)
                M_zscored.append(mi_zscore)
                
                #threshold, top nvox
                if abs_val:
                    values_sorted=np.sort(np.abs(values_zscored))
                else:
                    values_sorted=np.sort(values_zscored)
                
                top = values_sorted[-nvox::]
                thn=top[0]
                values_thn= np.copy(values_zscored)
                if abs_val:
                    values_thn[np.abs(values_thn)<thn]=0
                else:
                    values_thn[values_thn<thn]=0
                mi_thn = unmask(values_thn,mask_file)
                nib.save(mi_thn, savepath +'S%d_mean_connectivity_map_masked_zscored_top%d%s.nii' % (s,nvox,tail))
                M_topn.append(mi_thn)
                
                #threshold using th
                values_th= np.copy(values_zscored)
                if abs_val:
                    values_th[np.abs(values_th)<th]=0
                elif th >= 0:
                    values_th[values_th<th]=0
                else:
                    values_th[values_th>th]=0
                mi_th = unmask(values_th,mask_file)
                nib.save(mi_th, savepath +'S%d_mean_connectivity_map_masked_zscored_th%s%s.nii' % (s,th_str,tail))
                M_th.append(mi_th)
                
                del mean_img,mean_img_masked,mi_zscore,mi_thn,mi_th
                
                #correlations between seeds
                C = corr_seeds(fn_sub_ses,s,include,multiply,abs_val,mask_file,savepath,tail)
                C.to_csv(savepath +'S%d_seed_correlation_matrix%s_masked.txt' % (s,tail),float_format='%.6f')
                
            #correlations between sessions
            C = corr_ses(fn_sub,ses,include,multiply,abs_val,mask_file)
            i=0;
            for c in C:
                c.to_csv(savepath +'%s_correlation_matrix_masked.txt' % (include[i]),float_format='%.6f')    
                i+=1    
                    
            # mean over sessions
            M = image.mean_img(M)
            M_masked = image.mean_img(M_masked)
            M_zscored = image.mean_img(M_zscored)
            M_topn = image.mean_img(M_topn)
            M_th = image.mean_img(M_th)
            nib.save(M, savepath +'mean_connectivity_map%s.nii' % (tail))
            nib.save(M_masked, savepath +'mean_connectivity_map_masked%s.nii' % (tail))
            nib.save(M_zscored, savepath +'mean_connectivity_map_masked_zscored%s.nii' % (tail))
            nib.save(M_topn, savepath +'mean_connectivity_map_masked_zscored_top%s%s.nii' % (nvox,tail))
            nib.save(M_th, savepath +'mean_connectivity_map_masked_zscored_th%s%s.nii' % (th_str,tail))
            # correlations
            names= ['','_masked','_masked_zscored','_masked_zscored_top%d' % nvox,'_masked_zscored_th%s' % (th_str)]
            j=0
            for n in names:
                images=[]
                for s in ses:
                    img = nib.load(savepath +'S%d_mean_connectivity_map%s%s.nii' % (s,n,tail))
                    images.append(img) 
                C, mean_vals = compute_correlations_2(images,j,mask_file)
                C.to_csv(savepath +'correlation_matrix%s%s.txt' % (n,tail),float_format='%.6f')
                with open(savepath +'correlation_matrix%s%s.txt' % (n,tail), 'a') as f:
                    f.write("\nMean connectivities:\n")
                    for m in range(0,len(mean_vals)):
                        f.write('S%d: %.6f\n' % (m+1,mean_vals[m]))
                f.close()
                j+=1
        
     
     
         
         
         
         
         
         
         
         
         
