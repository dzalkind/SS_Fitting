import os, time
import numpy as np
import pandas as pd
import xarray as xr
from itertools import islice
import capytaine as cpt

def ordered_unique(seq):
    seen      = set()
    seen_add  = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

def load_nc(fname):
    # load Capytaine results data file
    return cpt.io.xarray.merge_complex_values(xr.open_dataset(fname))

def load_wamit(fname_base):

    outputs   = {}

    fname_hst = fname_base + '.hst'
    fname_1   = fname_base + '.1'
    fname_3   = fname_base + '.3'

    if os.path.exists(fname_hst):
        outputs.update(load_wamit_hst(fname_hst))
    if os.path.exists(fname_1):
        outputs.update(load_wamit_1(fname_1))
    if os.path.exists(fname_3):
        outputs.update(load_wamit_3(fname_3))

    return outputs

def load_wamit_hst(fname_hst):
    # Load .hst: linear hydrostatic-restoring matrix from the effects of water-plane area and center of buoyancy, C

    data_hst = np.loadtxt(fname_hst)

    C = np.zeros((6,6))

    for i_data in range(len(data_hst)):
        i      = int(data_hst[i_data, 0]) - 1
        j      = int(data_hst[i_data, 1]) - 1
        C[i,j] = data_hst[i_data, 2]

    outputs      = {}
    outputs['C'] = C

    return outputs

    
def load_wamit_1(fname_1):
    # Load .1: oscillation-frequency dependent hydrodnamic-added mass and damping matrices A and B
    
    # For negative or 0 periods, B matrix is all 0s and not given.  This results in non-constant number of columns in the file.
    # taking a first pass at loading this file to determine where the number of columns changes.
    fid    = open(fname_1, 'r')
    line   = fid.readline().strip().split()
    n_col0 = len(line)
    N      = 0
    while line:
        line = fid.readline().strip().split()
        N   += 1
        if len(line) > n_col0:
            break
    fid.close()

    # Load complete data set, join
    with open(fname_1, 'r') as fid:
        data1a = np.loadtxt(islice(fid, N))
        data1b = np.loadtxt(fid)
    fid.close()

    data1a = np.concatenate((data1a, np.zeros((np.shape(data1a)[0],1))), axis=1)
    data1  = np.concatenate((data1a, data1b), axis=0)

    # Format A and B Matrices
    periods   = ordered_unique(data1[:,0])
    n_periods = len(periods)

    A = np.zeros((n_periods, 6, 6))
    B = np.zeros((n_periods, 6, 6))
    for i_data in range(len(data1[:,0])):
        i = int(data1[i_data,1])-1
        j = int(data1[i_data,2])-1
        k = periods.index(data1[i_data,0])
        A[k, i, j] = data1[i_data,3]
        B[k, i, j] = data1[i_data,4]

    outputs = {}
    outputs['periods_rad'] = periods
    outputs['A']           = A
    outputs['B']           = B

    return outputs

def load_wamit_3(fname_3):
    # Load .3: wave-frequency and direction-dependent hydrodynamic-wave-excitation vector X

    data3 = np.loadtxt(fname_3)

    # Get Matrix Sizes
    periods     = ordered_unique(data3[:,0])
    n_periods   = len(periods)
    
    direction   = ordered_unique(data3[:,1])
    n_direction = len(direction)

    # Format X
    X_mag   = np.zeros((n_periods, n_direction, 6))
    X_phase = np.zeros((n_periods, n_direction, 6))
    X_re    = np.zeros((n_periods, n_direction, 6))
    X_imag  = np.zeros((n_periods, n_direction, 6))

    for i_data in range(len(data3[:,0])):

        i = periods.index(data3[i_data,0])
        j = direction.index(data3[i_data,1])
        k = int(data3[i_data,2])-1

        X_mag[i,j,k]   = data3[i_data, 3]
        X_phase[i,j,k] = data3[i_data, 4]
        X_re[i,j,k]    = data3[i_data, 5]
        X_imag[i,j,k]  = data3[i_data, 6]

    outputs = {}
    outputs['periods_diff']   = periods
    outputs['direction_diff'] = direction
    outputs['X_mag']          = X_mag
    outputs['X_phase']        = X_phase
    outputs['X_re']           = X_re
    outputs['X_imag']         = X_imag

    return outputs



if __name__ == "__main__":

    fname_base = 'C:/Users/egaertne/WT_Codes/models/openfast-dev/r-test/glue-codes/openfast/5MW_Baseline/HydroData/marin_semi'
    
    wamit = load_wamit(fname_base)
    print(wamit.keys())




