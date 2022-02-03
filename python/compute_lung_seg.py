# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: python/3.8
#     language: python
#     name: py3.8
# ---

"""
Lung MRI segmentation.
"""

import numpy as np
import collections
import os
import sys
import math
import time
import random
import logging
import pickle
from pathlib import Path
import matplotlib.pyplot as plt
#%matplotlib inline

# +
import torch
import torchvision
import scipy.io as sio # RR

Git_DIR = Path(__file__).parents[1].resolve()
sys.path.insert(1, str(Git_DIR))

import scipy


def compute_lung_seg(data, display=False):
    
    def prob_thresh(probs, device, p_thresh=0.5, params=None):
    
    	cpu_device = torch.device('cpu')

    	probs = probs.to(device=cpu_device)

    	RO = probs.shape[0]
    	E1 = probs.shape[1]
    
    	number_of_blobs = float("inf")
    	blobs = np.zeros((RO,E1))
    
    	mask = (probs > p_thresh).float()
    
    	return mask
    
    def cpad_2d(data, RO, E1):
        '''
        data: [dRO, dE1, N], padded it round center to [RO, E1, N]
        return value: (s_ro, s_e1), starting point of data in padded array
        '''

        dRO, dE1, N = data.shape

        s_ro = int((RO-dRO)/2)
        s_e1 = int((E1-dE1)/2)

        #print(data.shape, RO, E1, s_ro, s_e1)
        if(s_ro>=0):
            data_padded = np.zeros((RO, dE1, N))
            if(dRO>=RO):
                data_padded = data[s_ro:s_ro+RO,:,:]
            else:
                data_padded[s_ro:s_ro+dRO,:,:] = data
            data = data_padded
        else:
            data_padded = np.zeros((RO, dE1, N))
            if(dRO+s_ro+s_ro>RO):
                data_padded = data[-s_ro:(dRO+s_ro-1),:,:]
            else:
                data_padded = data[-s_ro:(dRO+s_ro),:,:]
            data = data_padded

        #print(data.shape)

        if(s_e1>=0):
            data_padded = np.zeros((RO, E1, N))
            if(dE1>=E1):
                data_padded = data[:,s_e1:s_e1+E1,:]
            else:
                data_padded[:,s_e1:s_e1+dE1,:] = data
            data = data_padded
        else:
            data_padded = np.zeros((RO, E1, N))
            if(dE1+s_e1+s_e1>E1):
                data_padded = data[:,-s_e1:(dE1+s_e1-1),:]
            else:
                data_padded = data[:,-s_e1:(dE1+s_e1),:]
            data = data_padded

        return data_padded, s_ro, s_e1


    # takes in data in the form of a numpy array [RO E1 N], and returns masks as a numpy array of same dimension
    device = 'cpu'
    model_file = str(Git_DIR) + "/best_lung_seg_model.pkl"
 
    model = torch.load(model_file)
    
    data = np.transpose(data, (2,0,1))
    N, orig_RO, orig_E1 = data.shape
    #print(orig_RO, orig_E1)
    RO = 384
    E1 = 384

    if torch.cuda.is_available():
        device = torch.device('cuda')
    else:
        device = torch.device('cpu')

    #print("Lung segmentation, device is ", device, file=sys.stderr)
    
    data_normalized = np.zeros((N, RO, E1), dtype='float32')
    
    for n in range(N):
        data2D, s_ro, s_e1 = cpad_2d(np.expand_dims(data[n,:,:], axis=2), RO, E1)
        data2D = data2D.squeeze()
        if np.max(data2D) != 0:
            data2D = data2D / np.max(data2D)
        data_normalized[n,:,:] = data2D


    im = np.expand_dims(data_normalized, axis=1)
    img_rs = torch.from_numpy(im.astype(np.float32)).to(device=device)
    
    
    model.to(device=device)  
    model.eval() 
    

    with torch.no_grad():
        t0 = time.time()

        scores = model(img_rs)

        probs = torch.sigmoid(scores)
        probs = probs.detach().cpu().float().numpy().astype(np.float32)

        ## Resize output mask to required iutput size 
        output = np.zeros((384, 384, N), np.float32)
        
        if display:
            fig, axes = plt.subplots(nrows=8, ncols=9, figsize=(50,50), sharex=True, sharey=True)
        
        for i in range(N):
            mask = prob_thresh(torch.from_numpy(probs[i,0,:,:]), 'cpu', p_thresh=0.5)
            output[:,:, i] = mask

            masked = np.ma.masked_where(mask == 0, mask)
            
            if display:
                y = i//9
                x = i%9
                axes[y,x].imshow(data[i,:,:], 'gray', clim=(0.0, 0.5))
                axes[y,x].grid(False)
                axes[y,x].imshow(masked, 'flag', interpolation='none', alpha=0.2)
        
        t1 = time.time()
        
        if display:
            plt.show()

        print("Mask computed in %f seconds" % (t1-t0), file=sys.stderr)
        output, _, _ = cpad_2d(output, orig_RO, orig_E1)

        return output

if __name__ == "__main__":

    mat_contents = sio.loadmat(str(sys.argv[1]));
    image = mat_contents['im']

    mask = compute_lung_seg(image, display=False)
  
    sio.savemat(str(str(sys.argv[2])), {"mask": mask})
    print("Mask saved.")





