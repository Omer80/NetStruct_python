# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 10:49:02 2016

@author: Omer Tzuk <omertz@post.bgu.ac.il>
"""
import numpy as np
#import h5py as h5

class NetStruct(object):
    def __init__(self,data,area=None,ind=None):
        self.data = data
        self.npop , num_columns = data.shape
        self.nloci = int(num_columns/2)
        if area is not None:
            self.area=area
        if ind is not None:
            self.ind=ind
        self.calc_allele_freq()
        self.createMatrix()

    def calc_allele_freq(self):
        max_l = np.max(self.data)
        freq  = np.zeros((self.nloci,max_l+1))
        for l in range(self.nloci):
            locus=np.concatenate((self.data[:,l*2],self.data[:,l*2+1]))
            count=np.bincount(np.append(locus,max_l))
            count[max_l]-=1
            zeros=count[0]
            count[0]=0
            freq[l]= count/float((2*self.npop)-zeros)
        self.freq_allele = freq

    def createMatrix(self):
        A = np.zeros((self.npop,self.npop))
        for i in np.arange(self.npop):
            for j in np.arange(i+1,self.npop):
                A[i,j]=self.calc_edge(i,j)
        self.A = A + A.T

    def calc_edge(self,i,j):
        Sij = np.empty(self.nloci)
        for l in range(self.nloci):
            Iac = int(self.data[i,l*2]==self.data[j,l*2])
            Iad = int(self.data[i,l*2]==self.data[j,l*2+1])
            Ibc = int(self.data[i,l*2+1]==self.data[j,l*2])
            Ibd = int(self.data[i,l*2+1]==self.data[j,l*2+1])
            Sij[l]=(1.0/4)*((1.0-self.freq_allele[l][self.data[i,l*2]])*(Iac+Iad)+(1.0-self.freq_allele[l][self.data[i,l*2+1]])*(Ibc+Ibd))
        return (1.0/self.nloci)*np.sum(Sij)


def readfile(fname):
    with open(fname) as f:
       ncols = len(f.readline().split(','))
    area = np.loadtxt(fname, delimiter=',', usecols=[0],dtype=str)
    ind  = np.loadtxt(fname, delimiter=',', usecols=[1],dtype=int)
    data = np.loadtxt(fname, delimiter=',', usecols=range(2,ncols),dtype=int)
    return data,area,ind

def main(args):
    print args.fname
    global test
    data,area,ind=readfile(args.fname)
    test = NetStruct(data,area,ind)

def add_parser_arguments(parser):
    parser.add_argument('--csvfile', type=str, nargs='?',
                        default='wildass4pops.csv',
                        dest='fname',
                        help='output file')
    return parser


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(prog='PROG', usage='%(prog)s [options]')
    parser = add_parser_arguments(parser)
    args = parser.parse_args()
    main(args)