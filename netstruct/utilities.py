# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: Omer Tzuk <omertz@post.bgu.ac.il>
"""
__version__= 1.0
__author__ = """Omer Tzuk (omertz@post.bgu.ac.il)"""
import numpy as np
import h5py

def readcsv(fname):
    with open(fname) as f:
       ncols = len(f.readline().split(','))
    area = np.loadtxt(fname, delimiter=',', usecols=[0],dtype=str)
    ind  = np.loadtxt(fname, delimiter=',', usecols=[1],dtype=str)
    data = np.loadtxt(fname, delimiter=',', usecols=range(2,ncols),dtype=str)
    return data,area,ind


def csv2h5(fname):
    if fname.endswith('.csv'):
        data_array,area,ind = readcsv(fname)
        fname = fname[:-4]
    else:
        data_array,area,ind = readcsv(fname+".csv")
    with h5py.File(fname+'.hdf5', 'w') as f:
        f.create_dataset("area", data=area, compression="gzip")
        f.create_dataset("ind", data=ind, compression="gzip")
        f.create_dataset("data", data=data_array, compression="gzip")

def saveGDmatrix(fname,area,ind,data_a,GDmatrix):
    if fname.endswith('.csv'):
        fname=fname[:-4]
    if fname.endswith(('.h5','.hdf5')):
        pass
    else:
        fname = fname+".hdf5"
    with h5py.File(fname, 'w') as f:
        f.create_dataset("area", data=area, compression="gzip")
        f.create_dataset("ind", data=ind, compression="gzip")
        f.create_dataset("data", data=data_a, compression="gzip")
        f.create_dataset("GDmatrix", data=GDmatrix, compression="gzip")

def addGDmatrix(fname,GDmatrix):
    if fname.endswith(('.h5','.hdf5')):
        pass
    else:
        fname = fname+".hdf5"
    with h5py.File(fname, 'a') as f:
        f.create_dataset("GDmatrix", data=GDmatrix, compression="gzip")

def loadGDmatrix(fname):
    if fname.endswith('.hdf5'):
        pass
    else:
        fname = fname+".hdf5"
    with h5py.File(fname, 'r') as f:
        area_a = f.get('area')
        area = np.array(area_a)
        ind_a = f.get('ind')
        ind = np.array(ind_a)
        GDmatrix_a = f.get('GDmatrix')
        GDmatrix = np.array(GDmatrix_a)
    return area,ind,GDmatrix

def loadh5(fname):
    if fname.endswith(('.h5','.hdf5')):
        pass
    else:
        fname = fname+".hdf5"
    with h5py.File(fname, 'r') as f:
        area = f.get('area')
        area = np.array(area)
        ind = f.get('ind')
        ind = np.array(ind)
        data = f.get('data')
        data = np.array(data)
    return data,area,ind

