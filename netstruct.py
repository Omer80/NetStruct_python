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
import zipfile
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
        # output to file CSV

    def Athreshold(self,threshold):
        return stats.threshold(self.A, threshmin=threshold,newval=0)

    def clustering_algorithm(self,graph,algorithm_num=1):
#        print graph.es["weight"]
#        print "min,max weight", min(graph.es["weight"]),max(graph.es["weight"])
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
    def FindCommunities(self,threshold=0.0,algorithm_num=1):
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
        matx = self.Athreshold(threshold)
        self.Astar = matx.sum()
        graph = igraph.Graph.Weighted_Adjacency(matx.tolist(),attr="weight",mode="UPPER")
        graph["name"] = "NetStruck Weighted Graph with threshold={0}".format(threshold)
        graph.vs["area"]=self.area
        graph.vs["name"]=map(str,self.ind)
        graph.vs["degree"]=matx.sum(axis=1)
        components = graph.components()
        initial_membership_num = []
        initial_membership_num.append(0)
        graphs = []
        membership_list = []
        for i,subgraph in enumerate(components.subgraphs()):
#            print i,subgraph.__class__,subgraph.vcount()
            if subgraph.vcount()>1:
                clusters=self.clustering_algorithm(subgraph,algorithm_num)
#                for c in clusters:
#                    print clusters.vs[c]["name"]
#                clusters = dendrogram.as_clustering()
                membership = clusters.membership
                membership = [ x + initial_membership_num[i] for x in membership]
                initial_membership_num.append(max(membership)+1)
            else:
                initial_membership_num.append(initial_membership_num[i]+1)
                membership =  [initial_membership_num[i]]
            membership_list.append(membership)
#        print membership_list
        membership_list = self.organize_in_decreasing_order(membership_list)
#        print membership_list
        for i,subgraph in enumerate(components.subgraphs()):
#            print membership_list[i]
            subgraph.vs["cluster"]=membership_list[i]
            graphs.append(subgraph)
            for vertex in subgraph.vs:
                graph.vs.find(str(vertex["name"]))["cluster"]=vertex["cluster"]
        self.Agraph = graph
        self.Amatx  = matx
        return graph

    def loopGDmat(self,threshold_min, threshold_max,threshold_delta,piecharts_fname="piecharts",zipped_pickle_fname="zipped_graphs", algorithm_num=1):
        """ """
        graphs=[]
        import os
        thresholds = np.arange(threshold_min, threshold_max+threshold_delta,threshold_delta)
        with zipfile.ZipFile(zipped_pickle_fname+".zip", 'w') as myzip:
            for threshold in thresholds:
                graph = self.FindCommunities(threshold,algorithm_num)
                graphs.append(graph)
                graph.write_picklez(fname=str(int(threshold*1000)))
                myzip.write(str(int(threshold*1000)))
                try:
                    os.remove(str(int(threshold*1000)))
                except OSError:
                    pass
        self.createPieCharts(thresholds,graphs,piecharts_fname=piecharts_fname)
#            graphs.append(self.FindCommunities(threshold,algorithm_num))
        return thresholds,graphs

    def createPieCharts(self,thresholds,graphs,piecharts_fname="piecharts"):
        """ """
#        import datetime
#        from matplotlib.backends.backend_pdf import PdfPages
#        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        areas = np.unique(self.area)
        self.dif_areas=areas
        num_columns = len(areas)+1
        num_rows    = len(thresholds)
        fig_x = num_columns * 5
        fig_y = num_rows * 5
#        with PdfPages(piecharts_fname+'.pdf') as pdf:
        fig, ax = plt.subplots(num_rows,num_columns,figsize=(fig_x,fig_y))
        for i,threshold in enumerate(thresholds):
            fracs = self.calcPiecharts(graphs[i])
            ax[i,0].set_axis_off()                                  # Turn off axes & ticks
            ax[i,0].text(0.5,0.5,threshold,ha='center',va='center',fontsize=30, fontweight='bold')
            for j,area in enumerate(areas):
                ax[i,j+1].set_title(area,fontsize=30, fontweight='bold')
                ax[i,j+1].pie(fracs[area])#, explode=explode, labels=labels,
                #                                autopct='%1.1f%%', shadow=True, startangle=90)
                ax[i,j+1].set_aspect('equal')
            plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
            plt.savefig(piecharts_fname+".pdf",format="pdf")
#            plt.close()
            # Adding metadata for PDF file
#            d = pdf.infodict()
#            d['Title'] = 'Pie charts'
#            d['Author'] = u'Jouni K. Sepp\xe4nen'
#            d['Subject'] = 'How to create a multipage pdf file and set its metadata'
#            d['Keywords'] = 'PdfPages multipage keywords author title subject'
#            d['CreationDate'] = datetime.datetime.today()

    def calcPiecharts(self,graph):
        """ """
        fracs = {}
        mat=np.asarray(zip(graph.vs["area"],graph.vs["cluster"]))
        for area in self.dif_areas:
#            print map(int,mat[where(mat[:,0]==area)][:,1])
            fracs[area]=np.bincount(map(int,mat[np.where(mat[:,0]==area)][:,1]))/float(len(map(int,mat[np.where(mat[:,0]==area)][:,1])))
#            print fracs[area]
#        print zip(self.area,graph.vs["cluster"])
        return fracs

    def SADanalysis(self,zipped_pickle_fname,threshold,csv_fname="SADresults.csv"):
        """ (fname_matrix,fname_communities)->
        Loads one of the pickled graphs from the zip file
        Loops all of the individuals of the graph
        and calculate SA on each individual of the graph
        returns: 1) CSV file that first column area, 2nd ind, 3rd cluster,
        4th strength of association
        2) pickled graph with a new attribute of strength of association
        """
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
            min_modularity=self.SA(graph,name)
            graph.vs.find(name)["SOA"]=min_modularity
            rows.append((graph.vs.find(name)["area"],
                      graph.vs.find(name)["name"],
                      graph.vs.find(name)["cluster"],
                      min_modularity))
        import csv
        with open(csv_fname, 'w') as csvfile:
            fieldnames = ['area', 'name','cluster','SOA']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
#                print row
                writer.writerow({'area':row[0],'name':row[1],'cluster':row[2],'SOA':row[3]})
        return graph
#        return f_SAD #save in CSV format

    def SA(self,graph,name):
        """ Loops on all of the clusters except the one from which name is from
        create graph object that is similar to the input, just with the ind
        belonging to a different graph, and then calculating the modularity
        of the new graph. Produces a minimum of a vector that gives
        the difference between the original modularity and each one of the
        calculated modularities (for different cluster).
        """
        clusters=np.unique(graph.vs["cluster"])
        original_cluster=graph.vs.find(name)["cluster"]
        original_modularity = graph.modularity(graph.vs["cluster"],"weight")
#        print clusters,original_cluster,original_modularity
        modularity_list = []
        for i,cluster in enumerate(clusters):
            if graph.vs.find(name)["cluster"]!=cluster:
                graph.vs.find(name)["cluster"]=cluster
                modularity_list.append(original_modularity-graph.modularity(graph.vs["cluster"],"weight"))
        graph.vs.find(name)["cluster"]=original_cluster
#        print modularity_list
        return np.min(np.array(modularity_list))
#    def SignificanceTest(self,)

    def organize_in_decreasing_order(self,membership_list):
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


def readfile(fname):
    with open(fname) as f:
       ncols = len(f.readline().split(','))
    area = np.loadtxt(fname, delimiter=',', usecols=[0],dtype=str)
    ind  = np.loadtxt(fname, delimiter=',', usecols=[1],dtype=str)
    data = np.loadtxt(fname, delimiter=',', usecols=range(2,ncols),dtype=str)
    return data,area,ind

def main(args):
    print (args.fname)
    global test,data,area,ind
    data,area,ind=readfile(args.fname)
#    test = NetStruct(data,area,ind)

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