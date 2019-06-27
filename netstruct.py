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
import numpy as np
import deepdish as dd
from scipy import stats
import igraph
import zipfile
from collections import Counter
import time

""" Matrix construction algorithms """
class buildGDmatrix(object):
    """ class for building Genetic Similarity matrix from dataset """
    def __init__(self,fname,outputfname=None):
        """ fname can end either with csv of hdf5/h5 
        the outputfile will be in hdf5 format, using deepdish.io tool
        in case outputfile has not been specified, the outputfile name will
        be the same as the input 
        (in case it was CSV it will be changed to hdf5)         
        """
        try:
            data,area,ind = read_data(fname)
        except IOError:
            print "Please specify full file name"
#        if fname.endswith('.csv'):
#            data,area,ind = readcsv(fname)
#        elif fname.endswith(('.h5','.hdf5')):
#            data,area,ind = loadh5(fname)
#        else:
#            print "Please specify full file name"
#            raise IOError
        start_time = time.time()
        print "Constructing genetic similarity matrix"
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
        print "Finished building genetic similarity matrix in : %f seconds" % (time.time() - start_time)
        if outputfname is None:
            if fname.endswith('.hdf5'):
                outputfname=fname[:-5]
            elif fname.endswith('.csv'):
                outputfname=fname[:-4]
            else:
                pass
        saveGDmatrix(outputfname,self.area,self.ind,self.data,self.A)
#        if fname.endswith('.csv'):# Returns a CSV file with the genetic similarity matrix
#            saveGDmatrix(fname,self.area,self.ind,self.data,self.A)
#        elif fname.endswith(('.h5','.hdf5')):# Returns an hdf5 file with the genetic similarity matrix
#            if outputfname is None:
#                addGDmatrix(fname,self.A)
#            else:
#                saveGDmatrix(outputfname,self.area,self.ind,self,data,self.A)


    def calc_edge(self,i,j):
        """ (genetic data, i and j)->(similarity value between individual i and individual j)
		Calculates the genetic similarity between individuals i and j, as defined by equation 3 in Greenbaum et al. 2016 """
        Sij = np.zeros(self.nloci)
        nzeros=0 # This is a counter for the missing loci in indiduals i or j (i.e. number of loci with missing information in either i or j or both)
        for l in range(self.nloci):
            a=self.data[i,l*2]
            b=self.data[i,l*2+1]
            c=self.data[j,l*2]
            d=self.data[j,l*2+1]
            if np.where(np.array([a,b,c,d])=='0')[0].size!=0:
                nzeros+=1 # If there is a "0" in any of the genotypes of i or j, the counter is increased
            else:
                Iac = int(a==c)
                Iad = int(a==d)
                Ibc = int(b==c)
                Ibd = int(b==d)
                fa=self.freq_allele[l][self.data[i,l*2]]
                fb=self.freq_allele[l][self.data[i,l*2+1]]
                Sij[l]=(1.0/4)*((1.0-fa)*(Iac+Iad)+(1.0-fb)*(Ibc+Ibd))
        if nzeros==self.nloci:
            Sijtot=0 # If i and j are missing information in all loce, we set their similarity to 0. Note - this should not happen in reasonable datasets
        else:
            Sijtot=(np.sum(Sij)/float((self.nloci-nzeros))) # The average similarity over all non-missing loci 
        return Sijtot

    def calc_allele_freq(self):
        """ (genetic data)->()
		Calculates the allele frequencies in all loci"""
        freq  = []
        for l in range(self.nloci):
            locus=np.concatenate((self.data[:,l*2],self.data[:,l*2+1])) # This creates a list of all gametes in the locus
            local_freq = {}
            total = locus.size - np.where(locus=='0')[0].size  # Number of non-missing gametes in the locus
            locus = Counter(locus) # Returns a counter struct, in which each value (gametes) from the locus appear with the number of times it repeats in the original array
            for allele,bincounter in locus.most_common(): # For each pair of gametes and number of times that it appears do
                local_freq[allele]=float(bincounter)/float(total)  # Calculation of the frequency of the specific gametes in the locus
            freq.append(local_freq)
        self.freq_allele = freq

    def build(self):
        """(genetic data)->(genetic similarity matrix)
		Builds the genetic similarity matrix """
        A = np.zeros((self.npop,self.npop))
        for i in np.arange(self.npop):
            for j in np.arange(i+1,self.npop):
                A[i,j]=self.calc_edge(i,j)
        self.A = A  + A.T  # The loop above builds the matrix for the upper half above the diagonal of the matrix. This line simply summs the matrix with its transose to create a symmetric simmilarity matrix. This works since calc(data,i,j)=calc(data,j,i)
        # output to file CSV
		
""" Communities detection algorithms """

def clustering_algorithm(graph,algorithm_num=6): #The default is algorithim 6
    """ (igraph graph object, algorithim number 1-6)->(an igraph cluster object with the detected community partition)
	Implements one of the igraph community detection algorithims, all with 'weighted networks' option"""
    algorithm_num=int(algorithm_num)
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

def FindCommunities(fname,threshold=0.0,algorithm_num=6):
    """ (threshold,algorithm_num)->(list_of_clusters_memberships)
    Implments community detection for a single threshold
	
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
    area,ind,GDmatrix=loadGDmatrix(fname)
    # matx=stats.threshold(GDmatrix, threshmin=threshold,newval=0) # Creates an adjacency matrix with a threshold applied to GDmatrix (i.e. all values below the threshold are taken to be 0)
    matx = GDmatrix
    matx[matx<threshold] = 0

    graph = igraph.Graph.Weighted_Adjacency(matx.tolist(),attr="weight",mode="UPPER")#Creates an undirected weighted igraph graph object with the adjacency matrix
    graph["name"] = "NetStruck Weighted Graph with threshold={0}".format(threshold)
    graph.vs["area"]=area
    graph.vs["name"]=map(str,ind)
    graph.vs["degree"]=matx.sum(axis=1)
    components = graph.components()
    initial_membership_num = []
    initial_membership_num.append(0)
    graphs = []
    membership_list = []
    for i,subgraph in enumerate(components.subgraphs()):# Network detection is performed for each component independintly (most igraph community detection algorithms do not work on disconnected graphs)
        if subgraph.vcount()>1:
            clusters=clustering_algorithm(subgraph,algorithm_num)
            membership = clusters.membership
            membership = [ x + initial_membership_num[i] for x in membership]
            initial_membership_num.append(max(membership)+1)
        else: # The case of singeltons
            initial_membership_num.append(initial_membership_num[i]+1)
            membership =  [initial_membership_num[i]]
        membership_list.append(membership)
    membership_list = organize_in_decreasing_order(membership_list)#sort communities from largest to smallest
    for i,subgraph in enumerate(components.subgraphs()):
        subgraph.vs["cluster"]=membership_list[i]
        graphs.append(subgraph)
        for vertex in subgraph.vs:
            graph.vs.find(str(vertex["name"]))["cluster"]=vertex["cluster"]
    return graph #returns the graph (igraph object) with an additional attribute to each node indicating its community (1 is the largest community, 2 after, and so on..)

def loopGDmat(data_fname,threshold_min, threshold_max,threshold_delta,algorithm_num,outputfname=None):
    """ This function performs community detection for several thresholds (from threshold_min to threshold_max in jumps of threshold_delta)
	It produces a ZIP file which contains:
		For each threshold, a ZIP file which contains
			1. CSV file with columns area,name,cluster for each threshold
			2. A PDF with the piecharts for each threshold
			3. A pickled (igraph format) file which contains the partitioned graph. This file is intended for the further analyses
    """
    print "Running command: ", data_fname,threshold_min, threshold_max,threshold_delta,algorithm_num,outputfname
    if outputfname is None:
        if data_fname.endswith('.hdf5'):
            outputfname=data_fname[:-5]
        else:
            outputfname=data_fname
#    area,ind,GDmatrix=loadGDmatrix(data_fname)
    areas=loadareas(data_fname)
#    areas = np.unique(area)
#    del ind
#    del area
#    del GDmatrix
    start_time = time.time()
    graphs=[]
    import os
    thresholds = np.arange(threshold_min, threshold_max+threshold_delta,threshold_delta)
    np.savetxt('thresholds.dat',map(str,map(int,(thresholds*1000))), fmt="%s")
    print "Output filename: ",outputfname
    with zipfile.ZipFile(outputfname+".zip", 'w') as myzip:
        myzip.write('thresholds.dat')
        try:
            os.remove('thresholds.dat')
        except OSError:
            pass
        for threshold in thresholds:
            graph = FindCommunities(data_fname,threshold,algorithm_num)
            graphs.append(graph)
            pickled_graph_name=str(int(threshold*1000))
            graph.write_picklez(fname="0_"+pickled_graph_name)
            myzip.write("0_"+pickled_graph_name)
            try:
                os.remove("0_"+pickled_graph_name)
            except OSError:
                pass
            graph_data = np.array([graph.vs['area'],graph.vs['name'],graph.vs['cluster']])
            np.savetxt("0_"+pickled_graph_name+".csv",graph_data.T,delimiter=',', fmt="%s")
            myzip.write("0_"+pickled_graph_name+".csv")
            try:
                os.remove("0_"+pickled_graph_name+".csv")
            except OSError:
                pass
    print "Finished community detection analysis in : %f seconds" % (time.time() - start_time)
    start_time = time.time()
    createPieCharts(areas,thresholds,graphs,piecharts_fname="piecharts_"+outputfname)#Creates and saves the PDF with piecharts
    print "Finished producing PDF with piecharts in : %f seconds" % (time.time() - start_time)

def createPieCharts(areas,thresholds,graphs,piecharts_fname="piecharts"):
    """ Creates the PDF with the piecharts for the different thresholds """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    num_columns = len(areas)+1
    num_rows    = len(thresholds)
    fig_x = num_columns * 5
    fig_y = num_rows * 5
    fig, ax = plt.subplots(num_rows,num_columns,figsize=(fig_x,fig_y))
    for i,threshold in enumerate(thresholds): #For each threshold create a piechart line
        fracs = calcPiecharts(graphs[i],areas) #This is the fraction of each community in each area
        ax[i,0].set_axis_off() # Turn off axes & ticks
        ax[i,0].text(0.5,0.5,threshold,ha='center',va='center',fontsize=30, fontweight='bold')
        for j,area in enumerate(areas):
            ax[i,j+1].set_title(area,fontsize=30, fontweight='bold')
            ax[i,j+1].pie(fracs[area])
            ax[i,j+1].set_aspect('equal')
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        plt.savefig(piecharts_fname+".pdf",format="pdf")

def calcPiecharts(graph,areas):
    """ Calculates the fraction of each community in each area """
    fracs = {}
    mat=np.asarray(zip(graph.vs["area"],graph.vs["cluster"]))
    for area in areas:
        fracs[area]=np.bincount(map(int,mat[np.where(mat[:,0]==area)][:,1]))/float(len(map(int,mat[np.where(mat[:,0]==area)][:,1])))
    return fracs

def SADanalysis(zipped_pickle_fname,threshold,csv_fname):
    """ (fname_matrix,fname_communities)->
    Loads one of the pickled graphs from the zip file
    Loops all of the individuals of the graph
    and calculate SA on each individual of the graph
    returns: 1) CSV file that first column area, 2nd ind, 3rd cluster,
    4th strength of association
    2) pickled graph with a new attribute of strength of association
    """
#    print str(int(threshold*1000))
    if csv_fname is None:
        csv_fname="SADresults.csv"
    print "Output file: ",csv_fname
    with zipfile.ZipFile(zipped_pickle_fname+".zip", 'r') as myzip:
#        myzip.extract('thresholds.dat')
#        file_thresholds = np.array(map(int,np.loadtxt('thresholds.dat')))
##        print file_thresholds
#        str_threshold=int(threshold*10000)
##        print str_threshold
#        threshold=str(file_thresholds[(np.abs(file_thresholds-str_threshold)).argmin()])
#        print "Extracting 0.%s threshold from zip file" % threshold
        myzip.extract("0_"+str(int(threshold*1000)))
    graph=igraph.Graph.Read_Picklez("0_"+str(int(threshold*1000))) #Reads the igraph graph object from the pickle-zipped file for the selected threshold
    import os
    try:
        os.remove("0_"+str(int(threshold*1000)))
    except OSError:
        pass
    names = graph.vs["name"]
    clusters=np.unique(graph.vs["cluster"])
    try:
        assert len(clusters)>1 #For a community of size 1, SAD analysis cannot be performed
    except AssertionError:
        print "Error: Cannot run SAD analysis since only one cluster (no structure) was detected"
        raise
    """ Save results into CSV file """
    rows=[]
    fieldnames = "Sampled population,Individual's ID,Assigned community,SA,Most associated alternative community"
    for name in names:
        indx,min_modularity=SA(graph,name,clusters) #For each individual, calculate the Strength of Association (SA)
        graph.vs.find(name)["SOA"]=min_modularity
        rows.append((graph.vs.find(name)["area"],
                  graph.vs.find(name)["name"],
                  graph.vs.find(name)["cluster"],
                  min_modularity,indx))
    np.savetxt(csv_fname,rows,delimiter=',',header=fieldnames, fmt="%s")
#    import csv
#    with open(csv_fname, 'w') as csvfile:
#        fieldnames = ['area', 'name','cluster','SAA','Most associated alternative community']
#        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
#        writer.writeheader()
#        for row in rows:
##                print row
#            writer.writerow({'area':row[0],'name':row[1],'cluster':row[2],'SAA':row[3],'Most associated alternative community':row[4]})
    return graph
#        return f_SAD #save in CSV format

def SA(graph,name,clusters):
    """ Returns the Strength of Association(SA) of "name".  	
	
	Loops on all of the clusters except the one from which name is from
    create graph object that is similar to the input, just with the ind
    belonging to a different graph, and then calculating the modularity
    of the new graph. Produces a minimum of a vector that gives
    the difference between the original modularity and each one of the
    calculated modularities (for different cluster).
    """
#    clusters=np.unique(graph.vs["cluster"])
    original_cluster=graph.vs.find(name)["cluster"]
    original_modularity = graph.modularity(graph.vs["cluster"],"weight")
#        print clusters,original_cluster,original_modularity
    modularity_list = []
    for i,cluster in enumerate(clusters):
        if original_cluster!=cluster:
            graph.vs.find(name)["cluster"]=cluster
            modularity_list.append([i,original_modularity-graph.modularity(graph.vs["cluster"],"weight")])
    graph.vs.find(name)["cluster"]=original_cluster
    modularity_list = np.array(modularity_list) #This is a list of length (number of clusters - 1) of the Strength of association of individual "name" to all clusters except the detected cluster of "name"
#        min_modularity_list = np.amin(np.array(modularity_list))
    return modularity_list[np.argmin(modularity_list[:,1])] #Returns the minimum of the list of association values; this is the Strength of Association (SA) of "name".
#    def SignificanceTest(self,)

def organize_in_decreasing_order(membership_list):
    """Reorganizes the communities list in decresing order from the largest community to the smallest"""
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


# Gone!



def main(args):
    if args.createGDmatrix is not None:
        buildGDmatrix(args.createGDmatrix,args.outfname)
    elif args.loopGDmat is not None:
        if args.outfname=="output_file":
            args.outfname = "Community_Detection_Results"
        input_fname,threshold_min,threshold_max,threshold_delta,algorithm_num = args.loopGDmat
        print "Processing command: ",input_fname,threshold_min,threshold_max,threshold_delta,algorithm_num,args.outfname
        loopGDmat(input_fname,float(threshold_min),float(threshold_max),float(threshold_delta)
                 ,algorithm_num,args.loopGDmat_results)
    elif args.SAD is not None:
        input_file, threshold = args.SAD
        SADanalysis(input_file,float(threshold),args.outfname)
    return 0

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

def saveGDmatrix(fname,area,ind,data_a,GDmatrix):
    if fname.endswith('.csv'):
        fname=fname[:-4]
        fname+=".hdf5"
    elif not fname.endswith(('.h5','.hdf5')):
        fname+=".hdf5"
    data_dict = {'area':area,'ind':ind,'data':data_a,'GDmatrix':GDmatrix}
    print "Writing to file ",fname
    dd.io.save(fname,data_dict,compression='blosc')

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

def add_parser_arguments(parser):
#    parser.add_argument('filename', type=str, nargs='?', default=None,
#                        dest='infname',
#                        help='input file')
    parser.add_argument('-o','--outputfile', type=str, nargs='?',
                        default=None,
                        dest='outfname',
                        help='output file')
    parser.add_argument('-b','--createGD', type=str, nargs='?',
                        default=None,
#                        action="store_true",
                        dest='createGDmatrix',
                        help='Create GD matrix from csv or hdf5 file')
    parser.add_argument('-c','--loopGDmat', type=str, nargs=5,
                        default=None,
                        metavar=('input_fname','threshold_min', 'threshold_max', 'threshold_delta','algorithm_num'),
                        help='Loop GD matrix, args are: min threshold, max threshold, delta of threshold, and algorithm type')
    parser.add_argument('-z','--loopGDmat_results', type=str, nargs='?',
                        default="Community_Detection_Results",
                        dest='loopGDmat_results',
                        help='loopGDmat results filename')                        
    parser.add_argument('-s','--SAD', type=str, nargs=2,
                        default=None,
                        metavar=('SADthreshold','threshold'),
                        help='SAD analysis file')                        
    return parser


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(prog='PROG', usage='%(prog)s [options]')
    parser = add_parser_arguments(parser)
    args = parser.parse_args()
    main(args)
