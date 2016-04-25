# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 10:49:02 2016

@author: Omer Tzuk <omertz@post.bgu.ac.il>
"""
__version__= 1.0
__author__ = """Omer Tzuk (omertz@post.bgu.ac.il)"""
import numpy as np
from scipy import stats
import igraph
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
        self.buildGDmatrix()

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

    def buildGDmatrix(self):
        A = np.zeros((self.npop,self.npop))
        for i in np.arange(self.npop):
            for j in np.arange(i+1,self.npop):
                A[i,j]=self.calc_edge(i,j)
        self.A = A  + A.T

    def Athreshold(self,threshold):
        return stats.threshold(self.A, threshmin=threshold,newval=0)

    def clustering_algorithm(self,graph,algorithm_num=1):
#        print graph.es["weight"]
#        print "min,max weight", min(graph.es["weight"]),max(graph.es["weight"])
        if algorithm_num==1:
            dendrogram = graph.community_label_propagation(weights="weight")
        elif algorithm_num==2:
            dendrogram = graph.community_edge_betweenness(weights="weight")
        elif algorithm_num==3:
            dendrogram = graph.community_fastgreedy(weights="weight")
        elif algorithm_num==4:
            dendrogram = graph.community_walktrap(weights="weight")
        elif algorithm_num==5:
            dendrogram = graph.community_spinglass(weights="weight")
        elif algorithm_num==6:
            dendrogram = graph.community_leading_eigenvector(weights="weight")
        else:
            print "No such option!"
            dendrogram = None
        return dendrogram

    def FindCommunities(self,threshold=0.0,algorithm_num=1):
        """ (threshold,algorithm_num)->(list_of_clusters_memberships)
        matx - a threshold for the adjacenncy matrix
        algorithm_num - an integer for choosing the clustering algorithm type.

        Choose one of the following options for community plus modularity detection
        algorithm_num chooses one of the following methods:
        1 : label.propagation.community
        2 : edge.betweenness.community
        3 : fastgreedy.community
        4 : walktrap.community
        5 : spinglass.community
        6 : leading.eigenvector.community
        """
        matx = self.Athreshold(threshold)
        graph = igraph.Graph.Weighted_Adjacency(matx.tolist(),attr="weight",mode="UPPER")
        graph["name"] = "NetStruck Weighted Adjacency with threshold={0}".format(threshold)
        graph.vs['area']=self.area
        graph.vs['ind']=self.ind
        components = graph.components()
        initial_membership_num = []
        initial_membership_num.append(0)
        graphs = {}
        for i,subgraph in enumerate(components.subgraphs()):
#            print i,subgraph.__class__,subgraph.vcount()
            if subgraph.vcount()>1:
                dendrogram=self.clustering_algorithm(subgraph,algorithm_num)
                clusters = dendrogram.as_clustering()
                membership = clusters.membership
                membership = [ x + initial_membership_num[i] for x in membership]
                initial_membership_num.append(max(membership)+1)
            else:
                initial_membership_num.append(initial_membership_num[i]+1)
                membership =  initial_membership_num[i]
            subgraph.vs["cluster"]=membership
            graphs[i]=subgraph
        return graphs

    def calc_edge(self,i,j):
        Sij = np.zeros(self.nloci)
        nzeros=0
        for l in range(self.nloci):
            a=self.data[i,l*2]
            b=self.data[i,l*2+1]
            c=self.data[j,l*2]
            d=self.data[j,l*2+1]
            if np.count_nonzero([a,b,c,d])!=4:
                nzeros+=1
            else:
                Iac = int(a==c)
                Iad = int(a==d)
                Ibc = int(b==c)
                Ibd = int(b==d)
                fa=self.freq_allele[l][self.data[i,l*2]]
                fb=self.freq_allele[l][self.data[i,l*2+1]]
                Sij[l]=(1.0/4)*((1.0-fa)*(Iac+Iad)+(1.0-fb)*(Ibc+Ibd))
        if nzeros==self.nloci:
            Sijtot=0
        else:
            Sijtot=(np.sum(Sij)/float((self.nloci-nzeros)))
        return Sijtot

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
                        help='input file')
    parser.add_argument('--outputfile', type=str, nargs='?',
                        default='test',
                        dest='outfname',
                        help='output file')
    return parser


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(prog='PROG', usage='%(prog)s [options]')
    parser = add_parser_arguments(parser)
    args = parser.parse_args()
    main(args)