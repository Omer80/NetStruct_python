# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: Omer Tzuk <omertz@post.bgu.ac.il>
"""
__version__= 1.0
__author__ = """Omer Tzuk (omertz@post.bgu.ac.il)"""
import numpy as np
#import h5py
import deepdish as dd

def read_data(fname):
    """ read initial data from either csv of hdf5/h5 files """
    if fname.endswith('.csv'):
        data,area,ind = readcsv(fname)
    elif fname.endswith(('.h5','.hdf5')):
        data,area,ind=readhdf5(fname)
    return data,area,ind

def readcsv(fname):
    """ read initial data from csv file """
    with open(fname) as f:
       ncols = len(f.readline().split(','))
    area = np.loadtxt(fname, delimiter=',', usecols=[0],dtype=str)
    ind  = np.loadtxt(fname, delimiter=',', usecols=[1],dtype=str)
    data = np.loadtxt(fname, delimiter=',', usecols=range(2,ncols),dtype=str)
    return data,area,ind

def readhdf5(fname):
    """ read initial data from hdf5 file """
    dict_data = dd.io.load(fname)
    data=dict_data['data']
    area=dict_data['area']
    ind =dict_data['ind']
    return data,area,ind

def csv2h5(fname):
    """ transform a csv initial file data to hdf5 format """
    data,area,ind = readcsv(fname)
    if fname.endswith('.csv'):
        fname=fname[:-4]
    data_dict = {'data':data,'area':area,'ind':ind}
    dd.io.save(fname,data_dict,compression='blosc')
#    if fname.endswith('.csv'):
#        data_array,area,ind = readcsv(fname)
#        fname = fname[:-4]
#    else:
#        data_array,area,ind = readcsv(fname+".csv")
#    with h5py.File(fname+'.hdf5', 'w') as f:
#        f.create_dataset("area", data=area, compression="gzip")
#        f.create_dataset("ind", data=ind, compression="gzip")
#        f.create_dataset("data", data=data_array, compression="gzip")

def saveGDmatrix(fname,area,ind,data_a,GDmatrix):
    if fname.endswith('.csv'):
        fname=fname[:-4]
        fname+=".hdf5"
    elif not fname.endswith(('.h5','.hdf5')):
        fname+=".hdf5"
    data_dict = {'area':area,'ind':ind,'data':data_a,'GDmatrix':GDmatrix}
    print "Writing to file ",fname
    dd.io.save(fname,data_dict,compression='blosc')
#    else:
#        fname = fname+".hdf5"
#    with h5py.File(fname, 'w') as f:
#        f.create_dataset("area", data=area, compression="gzip")
#        f.create_dataset("ind", data=ind, compression="gzip")
#        f.create_dataset("data", data=data_a, compression="gzip")
#        f.create_dataset("GDmatrix", data=GDmatrix, compression="gzip")
#
#def addGDmatrix(fname,GDmatrix):
#    if fname.endswith(('.h5','.hdf5')):
#        pass
#    else:
#        fname = fname+".hdf5"
#    with h5py.File(fname, 'a') as f:
#        f.create_dataset("GDmatrix", data=GDmatrix, compression="gzip")

def loadareas(fname):
    if fname.endswith('.hdf5'):
        pass
    else:
        fname = fname+".hdf5"
    return np.unique(dd.io.load(fname,'/area'))

def loadGDmatrix(fname):
    if fname.endswith('.hdf5'):
        pass
    else:
        fname = fname+".hdf5"
    data_dict = dd.io.load(fname)
    return data_dict['area'],data_dict['ind'],data_dict['GDmatrix']
#    with h5py.File(fname, 'r') as f:
#        area_a = f.get('area')
#        area = np.array(area_a)
#        ind_a = f.get('ind')
#        ind = np.array(ind_a)
#        GDmatrix_a = f.get('GDmatrix')
#        GDmatrix = np.array(GDmatrix_a)
#    return area,ind,GDmatrix
#
#def loadh5(fname):
#    if fname.endswith(('.h5','.hdf5')):
#        pass
#    else:
#        fname = fname+".hdf5"
#    with h5py.File(fname, 'r') as f:
#        area = f.get('area')
#        area = np.array(area)
#        ind = f.get('ind')
#        ind = np.array(ind)
#        data = f.get('data')
#        data = np.array(data)
#    return data,area,ind

