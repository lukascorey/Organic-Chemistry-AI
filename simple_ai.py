#!/usr/bin/env python

import networkx as nx
import matplotlib.pyplot as plt
from copy import deepcopy
import time
from helpers import draw_graph, name_products, test_molecule_isomorphism, display_molecules, molecule_to_drawgraph, carbon_lst, print_path
from helpers import write_graph_to_file, load_graph_from_file
from reactions import hydrobromination, hydration, dehydration
from classes import atom, molecule, molecule_path_node
from methods import brute_force

# LOAD GRAPH INFO IF PREVIOUSLY RUN
(ad_list, graph_store) = load_graph_from_file("list_store.txt", "graph_store.txt")

# ADMINITRATIVE, IGNORE THESE BASICALLY 
reaction_dict = {"hydrobromination": hydrobromination, "hydration": hydration, "dehydration": dehydration}
carbons = carbon_lst(10)
oxygen1 = atom("oxygen", "oxygen1")
hydrogen1 = atom("hydrogen", "hydrogen1")

# HARD CODE IN RETROSYNTHESIS PROBLEM
start = molecule()
start.add_atoms(carbons[0:5] + [oxygen1, hydrogen1])
start.add_bonds([(carbons[0], carbons[1], 1), (carbons[1], carbons[2], 2),(carbons[2], carbons[3], 1), (carbons[3], carbons[4], 2)])
start.add_bonds([(carbons[3], oxygen1, 1), (oxygen1, hydrogen1, 1)])


goal = molecule()

bromine1 = atom("bromine", "bromine1")
bromine2 = atom("bromine", "bromine2")
bromine3 = atom("bromine", "bromine3")
goal.add_atoms(carbons[5:10] + [bromine1, bromine2, bromine3])
goal.add_bonds([(carbons[5], carbons[6], 1), (carbons[6], carbons[7], 1), (carbons[7], carbons[8], 1), (carbons[8], carbons[9], 1)])
goal.add_bonds([(carbons[6], bromine1, 1), (carbons[7], bromine2, 1), (carbons[8], bromine3, 1)])

# START TIMER
beg = time.time()
done = False

# can we find the path in our stored graph of molecules?
if graph_store != None:
    end_present = -1
    for i in range(len(graph_store)): 
        if (test_molecule_isomorphism(graph_store[i][0], goal)):
            end_present = i
            print("found finish molecule")
    if (end_present != -1):
        print("Found path from start to finish in stored graph")
        done = True
        end = time.time()
        molecule_to_drawgraph(start, "start")
        molecule_to_drawgraph(goal, "goal")
        for i in graph_store: 
            molecule_to_drawgraph(i[0], i[1])
# if not, have to brute force
if not done:
    print_list = brute_force(start, goal, 3, reaction_dict)
    end = time.time()
    # PRINT RESULTS 
    for i in print_list: 
        molecule_to_drawgraph(i[0], i[1])
    ad_list = []
    # CREATE ADJACENCY LIST AND GRAPH TO STORE
    for i in range(len(print_list[2:])):
        temp = []
        for j in range(i+1):
            temp.append(0)
        temp.append(1)
        for j in range(len(print_list)-i-3):
            temp.append(0)
        ad_list.append(temp)
    print("generated adjacency list: ")
    print(ad_list)
    write_graph_to_file("list_store.txt", "graph_store.txt", ad_list, print_list[2:])
# PRINT TIME TAKEN
print("took " + str(end-beg) + " seconds")