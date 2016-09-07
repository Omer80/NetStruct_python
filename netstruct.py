# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 10:49:02 2016

@author: Omer Tzuk <omertz@post.bgu.ac.il>

- Create argparse and destroy class
1) Construct Genetic Distance Matrix - buildGDmatrix
2) Community Detection - comDetect
3) Strength of Assosiation - SAD
4) Significance 
"""
__version__= 1.0
__author__ = """Omer Tzuk (omertz@post.bgu.ac.il)"""
import numpy as np
from scipy import stats
import igraph
import zipfile
from collections import Counter
import time
import utilities
""" Matrix construction algorithms """
class buildGDmatrix(object):
    """ class for building the matrix from dataset """
    def __init__(self,fname,outputfname=None):
        if fname.endswith('.csv'):
            data,area,ind = utilities.readcsv(fname)
        elif fname.endswith(('.h5','.hdf5')):
            data,area,ind = utilities.loadh5(fname)
        else:
            print "Please specify full file name"
            raise IOError
        start_time = time.time()
        print "Constructing network"
        self.data = data
        self.npop , num_columns = data.shape
        self.nloci = int(num_columns/2)
        if area is not None:
            self.area=area
        if ind is not None:
            self.ind=ind
        self.calc_allele_freq()
        print "Finished calculating alleles frequencies in : %f seconds" % (time.time() - start_time)
        start_time = time.time()
        self.build()
        print "Finished building GD matrix in : %f seconds" % (time.time() - start_time)
        if fname.endswith('.csv'):
            utilities.saveGDmatrix(fname,self.area,self.ind,self.data,self.A)
        elif fname.endswith(('.h5','.hdf5')):
            if outputfname is None:
                utilities.addGDmatrix(fname,self.A)
            else:
                utilities.saveGDmatrix(outputfname,self.area,self.ind,self,data,self.A)


    def calc_edge(self,i,j):
        """ """
        Sij = np.zeros(self.nloci)
        nzeros=0
        for l in range(self.nloci):
            a=self.data[i,l*2]
            b=self.data[i,l*2+1]
            c=self.data[j,l*2]
            d=self.data[j,l*2+1]
            if np.where(np.array([a,b,c,d])=='0')[0].size!=0:
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

    def calc_allele_freq(self):
        """ """
        freq  = []
        for l in range(self.nloci):
            locus=np.concatenate((self.data[:,l*2],self.data[:,l*2+1]))
            local_freq = {}
            total = locus.size - np.where(locus=='0')[0].size
            locus = Counter(locus)
            for allele,bincounter in locus.most_common():
                local_freq[allele]=float(bincounter)/float(total)
            freq.append(local_freq)
        self.freq_allele = freq

    def build(self):
        """ """
        A = np.zeros((self.npop,self.npop))
        for i in np.arange(self.npop):
            for j in np.arange(i+1,self.npop):
                A[i,j]=self.calc_edge(i,j)
        self.A = A  + A.T
        # output to file CSV
""" Communities detection algorithms """
def clustering_algorithm(graph,algorithm_num=1):
    """ """
    if algorithm_num==1:
        clusters = graph.community_label_propagation(weights="weight")
    elif algorithm_num==2:
        dendrogram = graph.community_edge_betweenness(weights="weight")
        clusters = dendrogram.as_clustering()
    elif algorithm_num==3:
        dendrogram = graph.community_fastgreedy(weights="weight")
        clusters = dendrogram.as_clustering()
    elif algorithm_num==4:
        dendrogram = graph.community_walktrap(weights="weight")
        clusters = dendrogram.as_clustering()
    elif algorithm_num==5:
        clusters = graph.community_spinglass(weights="weight")
    elif algorithm_num==6:
        clusters = graph.community_leading_eigenvector(weights="weight")
    else:
        print ("No such option!")
        clusters = None
    return clusters
# the function should accept a data file of the matrix output of buildGDmatrix
def FindCommunities(fname,threshold=0.0,algorithm_num=1):
    """ (threshold,algorithm_num)->(list_of_clusters_memberships)
    matx - a threshold for the adjacenncy matrix
    algorithm_num - an integer for choosing the clustering algorithm type.

    Choose one of the following options for community plus modularity detection
    algorithm_num chooses one of the following methods:
    1 : label.propagation.community    -> returns VertexClustering object.
    2 : edge.betweenness.community     -> returns VertexDendrogram object.
    3 : fastgreedy.community           -> returns VertexDendrogram object.
    4 : walktrap.community             -> returns VertexDendrogram object.
    5 : spinglass.community            -> returns VertexClustering object.
    6 : leading.eigenvector.community  -> returns VertexClustering object.
    """
    area,ind,GDmatrix=utilities.loadGDmatrix(fname)
    matx=stats.threshold(GDmatrix, threshmin=threshold,newval=0)
    graph = igraph.Graph.Weighted_Adjacency(matx.tolist(),attr="weight",mode="UPPER")
    graph["name"] = "NetStruck Weighted Graph with threshold={0}".format(threshold)
    graph.vs["area"]=area
    graph.vs["name"]=map(str,ind)
    graph.vs["degree"]=matx.sum(axis=1)
    components = graph.components()
    initial_membership_num = []
    initial_membership_num.append(0)
    graphs = []
    membership_list = []
    for i,subgraph in enumerate(components.subgraphs()):
        if subgraph.vcount()>1:
            clusters=clustering_algorithm(subgraph,algorithm_num)
            membership = clusters.membership
            membership = [ x + initial_membership_num[i] for x in membership]
            initial_membership_num.append(max(membership)+1)
        else:
            initial_membership_num.append(initial_membership_num[i]+1)
            membership =  [initial_membership_num[i]]
        membership_list.append(membership)
    membership_list = organize_in_decreasing_order(membership_list)
    for i,subgraph in enumerate(components.subgraphs()):
        subgraph.vs["cluster"]=membership_list[i]
        graphs.append(subgraph)
        for vertex in subgraph.vs:
            graph.vs.find(str(vertex["name"]))["cluster"]=vertex["cluster"]
    return graph

def loopGDmat(data_fname,threshold_min, threshold_max,threshold_delta,piecharts_fname="piecharts",zipped_pickle_fname="zipped_graphs", algorithm_num=1):
    """ Produce CSV file with columns area,name,cluster for each threshold
        and fix pie charts creation
    """
    area,ind,GDmatrix=utilities.loadGDmatrix(data_fname)
    areas = np.unique(area)
    del ind
    del GDmatrix
    start_time = time.time()
    graphs=[]
    import os
    thresholds = np.arange(threshold_min, threshold_max+threshold_delta,threshold_delta)
    print zipped_pickle_fname
    with zipfile.ZipFile(zipped_pickle_fname+".zip", 'w') as myzip:
        for threshold in thresholds:
            graph = FindCommunities(data_fname,threshold,algorithm_num)
            graphs.append(graph)
            pickled_graph_name=str(int(threshold*1000))
            graph.write_picklez(fname=pickled_graph_name)
            myzip.write(pickled_graph_name)
            try:
                os.remove(pickled_graph_name)
            except OSError:
                pass
    print "Finished calculating GDmat in : %f seconds" % (time.time() - start_time)
    start_time = time.time()
    createPieCharts(areas,thresholds,graphs,piecharts_fname=piecharts_fname)
    print "Finished producing PDFs in : %f seconds" % (time.time() - start_time)

def createPieCharts(areas,thresholds,graphs,piecharts_fname="piecharts"):
    """ """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    num_columns = len(areas)+1
    num_rows    = len(thresholds)
    fig_x = num_columns * 5
    fig_y = num_rows * 5
    fig, ax = plt.subplots(num_rows,num_columns,figsize=(fig_x,fig_y))
    for i,threshold in enumerate(thresholds):
        fracs = calcPiecharts(graphs[i],areas)
        ax[i,0].set_axis_off() # Turn off axes & ticks
        ax[i,0].text(0.5,0.5,threshold,ha='center',va='center',fontsize=30, fontweight='bold')
        for j,area in enumerate(areas):
            ax[i,j+1].set_title(area,fontsize=30, fontweight='bold')
            ax[i,j+1].pie(fracs[area])
            ax[i,j+1].set_aspect('equal')
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        plt.savefig(piecharts_fname+".pdf",format="pdf")

def calcPiecharts(graph,areas):
    """ """
    fracs = {}
    mat=np.asarray(zip(graph.vs["area"],graph.vs["cluster"]))
    for area in areas:
        fracs[area]=np.bincount(map(int,mat[np.where(mat[:,0]==area)][:,1]))/float(len(map(int,mat[np.where(mat[:,0]==area)][:,1])))
    return fracs

def SADanalysis(zipped_pickle_fname,threshold,csv_fname="SADresults.csv"):
    """ (fname_matrix,fname_communities)->
    Loads one of the pickled graphs from the zip file
    Loops all of the individuals of the graph
    and calculate SA on each individual of the graph
    returns: 1) CSV file that first column area, 2nd ind, 3rd cluster,
    4th strength of association
    2) pickled graph with a new attribute of strength of association
    """
    print str(int(threshold*1000))
    print "Outputfile ",csv_fname
    with zipfile.ZipFile(zipped_pickle_fname+".zip", 'r') as myzip:
        myzip.extract(str(int(threshold*1000)))
    graph=igraph.Graph.Read_Picklez(str(int(threshold*1000)))
    import os
    try:
        os.remove(str(int(threshold*1000)))
    except OSError:
        pass
    names = graph.vs["name"]
    rows=[]
    for name in names:
        indx,min_modularity=SA(graph,name)
        graph.vs.find(name)["SOA"]=min_modularity
        rows.append((graph.vs.find(name)["area"],
                  graph.vs.find(name)["name"],
                  graph.vs.find(name)["cluster"],
                  min_modularity,indx))
    import csv
    with open(csv_fname, 'w') as csvfile:
        fieldnames = ['area', 'name','cluster','SOA','MAC']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
#                print row
            writer.writerow({'area':row[0],'name':row[1],'cluster':row[2],'SOA':row[3],'MAC':row[4]})
    return graph
#        return f_SAD #save in CSV format

def SA(graph,name):
    """ Loops on all of the clusters except the one from which name is from
    create graph object that is similar to the input, just with the ind
    belonging to a different graph, and then calculating the modularity
    of the new graph. Produces a minimum of a vector that gives
    the difference between the original modularity and each one of the
    calculated modularities (for different cluster).
    
    !!! A BUG - what happens when there is only one cluster??
    """
    clusters=np.unique(graph.vs["cluster"])
    original_cluster=graph.vs.find(name)["cluster"]
    original_modularity = graph.modularity(graph.vs["cluster"],"weight")
#        print clusters,original_cluster,original_modularity
    modularity_list = []
    for i,cluster in enumerate(clusters):
        if original_cluster!=cluster:
            graph.vs.find(name)["cluster"]=cluster
            modularity_list.append([i,original_modularity-graph.modularity(graph.vs["cluster"],"weight")])
    graph.vs.find(name)["cluster"]=original_cluster
    modularity_list = np.array(modularity_list)
#        min_modularity_list = np.amin(np.array(modularity_list))
    return modularity_list[np.argmin(modularity_list[:,1])]
#    def SignificanceTest(self,)

def organize_in_decreasing_order(membership_list):
    conlist = np.concatenate(membership_list)
    count=np.bincount(conlist)
    ordered_list = []
    for i in range(len(count)):
        ordered_list.append(np.argmax(count))
        count[np.argmax(count)]=0
#        print ordered_list
    maxvalue = max(conlist)
    locations = []
    for number in range(maxvalue+1):
        locs = []
        for i in range(len(membership_list)):
            for j in range(len(membership_list[i])):
                if membership_list[i][j]==number:
                    locs.append((i,j))
        locations.append(locs)
    for i,value in enumerate(ordered_list):
        for location in locations[value]:
#                print location
#                print location[0],location[1]
            membership_list[location[0]][location[1]]=i
    return membership_list


def writefile(graphs):
    import csv
    with open('temp.csv', 'w') as csvfile:
        fieldnames = ['area', 'name','cluster']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        tonumpy = []
        for subgraph in graphs:
            rows = zip(subgraph.vs["area"],subgraph.vs["name"],subgraph.vs["cluster"])
            for row in rows:
#                print str(row[2])
                writer.writerow({'area': row[0], 'name':row[1] ,'cluster':row[2]})
                tonumpy.append((row[0],row[1],row[2]))
    return np.asarray(tonumpy)




def main(args):
    if args.createGDmatrix is not None:
        buildGDmatrix(args.createGDmatrix,args.outfname)
    elif args.loopGDmat is not None:
        fname,threshold_min,threshold_max,threshold_delta = args.loopGDmat
        print fname,threshold_min,threshold_max,threshold_delta
        loopGDmat(fname,float(threshold_min),float(threshold_max),float(threshold_delta)
                 ,args.piecharts_fname,args.zipped_pickle_fname, args.algorithm_num)
    elif args.SAD is not None:
        input_file, threshold = args.SAD
        if args.outfname is not None:
            SADanalysis(input_file,float(threshold),args.outfname)
        else:
            SADanalysis(input_file,float(threshold))
    return 0

def add_parser_arguments(parser):
    parser.add_argument('--csvfile', type=str, nargs='?',
                        default=None,
                        dest='fname',
                        help='input file')
    parser.add_argument('--outputfile', type=str, nargs='?',
                        default=None,
                        dest='outfname',
                        help='output file')
    parser.add_argument('-c','--createGD', type=str, nargs='?',
                        default=None,
                        dest='createGDmatrix',
                        help='Create GD matrix from csv or hdf5 file')
    parser.add_argument('-l','--loopGDmat', type=str, nargs=4,
                        default=None,
                        metavar=('loopGDmat', 'threshold_min', 'threshold_max', 'threshold_delta'),
                        help='Loop GD matrix, with min, max and delta of threshold')
    parser.add_argument('-a','--algorithm_num', type=int,
                        default=1,
                        dest='algorithm_num',
                        help='Choose one of six algorithms of community detection.')
    parser.add_argument('-z','--zipped_pickle_fname', type=str, nargs='?',
                        default="zipped_pickled_graphs",
                        dest='zipped_pickle_fname',
                        help='zipped_pickle_fname file')                        
    parser.add_argument('-p','--piecharts_fname', type=str, nargs='?',
                        default="piecharts",
                        dest='piecharts_fname',
                        help='piecharts_fname file')                        
    parser.add_argument('-s','--SAD', type=str, nargs=2,
                        default=None,
                        metavar=('SADanalysis','threshold'),
                        help='SAD analysis file')                        
    return parser


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(prog='PROG', usage='%(prog)s [options]')
    parser = add_parser_arguments(parser)
    args = parser.parse_args()
    main(args)