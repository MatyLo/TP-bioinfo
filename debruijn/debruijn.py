#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import textwrap
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file (default contigs.fasta)")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as an image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
	#print("aaa")
	with open (fastq_file,"r") as monfich:
		for i in monfich:
			yield next(monfich).strip("\n")
			next(monfich)
			next(monfich)
		pass


def cut_kmer(read, kmer_size):
    for i,_ in enumerate(read[:len(read)-kmer_size+1]):
            yield read[i:i+kmer_size]
    pass


def build_kmer_dict(fastq_file, kmer_size):
    #On crée le dictionnaire vide
    dict_kmer = {}

    #On récupère la séquence grace a read_fastq
    for read in read_fastq(fastq_file):
        for k_mer in cut_kmer(read, kmer_size):
            #Si le k_mer se trouve déjà dans la table de hachage, on ajoute 1 au compteur
            if k_mer in dict_kmer:
                dict_kmer[k_mer] += 1
            #Sinon, on ajoute une nouvelle clé à la table avec la valeur 1
            else:
                dict_kmer[k_mer] = 1

    return dict_kmer
    pass


def build_graph(kmer_dict):
    digraph = nx.DiGraph()
    for kmer, occurrence in kmer_dict.items():
        digraph.add_edge(kmer[:-1],kmer[1:],weight=occurrence)
    return digraph
    pass


def est_dans_la_liste(a,liste):
    for i in liste:
        if(i==a):
            return True
    return False 

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    #On récupère les noeuds dans path_list
    list = []
    list.append(path_list[0][0])
    for i in path_list:
        print(i)
        for j in range(1,len(i)-1):
            if(est_dans_la_liste(i[j],list)==False):
                print("false")
                list.append(i[j])
    list.append(path_list[0][len(path_list[0])-1])
    print(list)
    print("Les noeuds avant remove: ",graph.nodes())
    #Pour supprimer tous les noeuds du chemin
    if(delete_entry_node==True and delete_sink_node==True):
        graph.remove_nodes_from(list)
    #Pour supprimer tous les noeuds du chemin sauf le dernier
    elif(delete_entry_node==True and delete_sink_node==False):
        graph.remove_nodes_from(list[:-1])
    #Pour supprimer tous les noeuds du chemin sauf le premier
    elif(delete_entry_node==False and delete_sink_node==True):
        graph.remove_nodes_from(list[1:])
    #Pour supprimer tous les noeuds du chemin sauf le premier et le dernier
    elif(delete_entry_node==False and delete_sink_node==False):
        print("F-F")
        graph.remove_nodes_from(list[1:-1])
        print("Les noeuds apres remove:",graph.nodes())
    return graph


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    #Comparaison selon la frequence(poids)
    if(statistics.stdev(weight_avg_list)>0):
        #On choisit le chemin dont le poids est le plus élevé
        max1=weight_avg_list[0]
        indice1=0
        for j in range(len(weight_avg_list)):
            if(weight_avg_list[j]>max1):
                max1=weight_avg_list[j]
                indice1=j
        path_list.pop(indice1)
        remove_paths(graph,path_list,delete_entry_node,delete_sink_node)
        return graph
    else:
        #Comparaison de la longueur
        if(statistics.stdev(path_length)>0):
            max2=path_length[0]
            indice2=0
            for j in range(1,len(path_length)):
                if(path_length[j]>max2):
                    max2=path_length[j]
                    indice2=j
            path_list.pop(indice2)
            remove_paths(graph,path_list,delete_entry_node,delete_sink_node)
            return  graph
        else:
            indice3 = random.randint(0,len(path_length))
            path_list.pop(indice3)
            remove_paths(graph,path_list,delete_entry_node,delete_sink_node)
            return graph
    #return graph
    pass

def path_average_weight(graph, path):
    """Compute the weight of a path"""
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])


def solve_bubble(graph, ancestor_node, descendant_node):
    #Détermination des chemins simples entre 2 noeuds
    path_list = []
    path_length = []
    weight_list = []
    weight_avg_list = []
    for path in nx.all_simple_paths(graph,ancestor_node,descendant_node):
        print(path)
        print("taille: ",len(path))
        path_list.append(path)
        path_length.append(len(path))
        a = graph.subgraph(path).edges(data='weight')
        c=0
        for i in a:
            c+=i[2]
        #print(list(b)[0])
        print("c: ",c)
        weight = c/(len(path)-1)
        print("w: ",weight)
        weight_avg_list.append(weight)
    print(path_list)
    print(path_length)
    print("Poids moyen",weight_avg_list)
    graph = select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False)
    return graph

def simplify_bubbles(graph):
    bubble = False
    for n in graph:
        print(n)
        liste_predecesseurs = list(graph.predecessors(n))
        #print(liste_predecesseurs)
        #print("len", len(list(liste_predecesseurs)))
        if(len(list(liste_predecesseurs)) > 1):
            print("non null")
            print("liste prédecesseurs: ",liste_predecesseurs)
            for i in range(len(liste_predecesseurs)-1):
                print("i:" ,liste_predecesseurs[i])
                j = i+1
                print("j: ",liste_predecesseurs[j])
                noeud_ancetre = nx.lowest_common_ancestor(graph,liste_predecesseurs[i],liste_predecesseurs[j])
                print("noeud ancetre: ",noeud_ancetre)
                if(noeud_ancetre != None):
                    bubble = True
                    break
        if(bubble == True):
            break
    if(bubble):
        graph = simplify_bubbles(solve_bubble(graph,noeud_ancetre,n))
    return graph
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    list_start = []
    for node in graph.nodes():
        if(len(list(graph.predecessors(node)))==0):
            list_start.append(node)
    return list_start
    pass

def get_sink_nodes(graph):
    list_sink = []
    for node in graph:
        if(len(list(graph.successors(node)))==0):
            list_sink.append(node)
    return list_sink
    pass

def get_contigs(graph, starting_nodes, ending_nodes):
    tuple=[]
    for noeud_start in starting_nodes:
        for noeud_end in ending_nodes:
            if(nx.has_path(graph,noeud_start,noeud_end)==True):
                contig = noeud_start
                for path in nx.all_simple_paths(graph,noeud_start,noeud_end):
                    for node in path[1:]:
                        contig+=node[-1]
                tuple.append([contig, len(contig)])
    return tuple
    pass

def save_contigs(contigs_list, output_file):
    with open (output_file,"r") as out:
        i=0
        for tuple in contigs_list:
            print("contig_",i," len=",tuple[1])
            textwrap.fill(tuple[0],width=80)
            i+=1


def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    
    #read_fastq(args.dest)

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    #if args.graphimg_file:
        #draw_graph(graph, args.graphimg_file)


if __name__ == '__main__':
    main()
