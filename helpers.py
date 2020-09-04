#! usr/bin/env/python3

import networkx as nx
import matplotlib.pyplot as plt
from classes import atom, molecule, molecule_path_node
import time
import json
import os
import pickle

prefixes = ["first", "second", "third", "fourth", "fifth", 'sixth', "seventh", "eigtht", "ninth", "tenth", "eleventh"]
prefixes = prefixes + ["twelfth", "thirteen", "fourteenth", "fifteenth", "sixteenth", "seventeenth", "eighteenth", "nineteenth", "twentieth"]

def write_graph_to_file(list_filename,graph_filename, ad_list, graphs_list):
    with open(list_filename, 'w+') as filehandle:
        json.dump(ad_list, filehandle)
    with open(graph_filename, "wb+") as filehandle:
        pickle.dump(graphs_list, filehandle)
    pass

def load_graph_from_file(list_filename, graph_filename):
    if os.path.exists(list_filename):
        with open(list_filename, 'r') as filehandle:
            x = len(filehandle.read())
        if x>10:
            with open(list_filename, "r") as filehandle:
                ad_list = json.load(filehandle)
        else: 
            return (None, None)
    else:
        return (None, None)
    print("loaded ad list as")
    print(ad_list)
    if os.path.exists(graph_filename):
        with open(graph_filename, 'rb') as filehandle: 
            graphs_list = pickle.load(filehandle)
            print("loaded in graphs list, first is ")
            print(graphs_list[0])
    else:
        return (None, None)
    return (ad_list, graphs_list)
    


def molecule_to_drawgraph(m, title):
    print(title)
    print(m.bonds)
    graph = []
    labels = []
    for i in m.bonds: 
        graph.append((i[0].label, i[1].label))
        labels.append(i[2])
    draw_graph(graph, title, labels)
    pass

def display_molecules(lst, titles):
    for i in range(len(lst)):
        molecule_to_drawgraph(lst[i], titles[i])

def print_path(lst, counter):
    numlst = []
    while lst[counter].recent_reac != "start":
        numlst.append(counter)
        counter = lst[counter].pathindex
    string = "start"
    display_molecules([lst[0].molec], ["start"])
    for i in reversed(numlst): 
        string = string + ", " + lst[i].recent_reac
        display_molecules([lst[i].molec], [string])
    return 0

def get_path(lst, counter):
    returnlist = []
    numlst = []
    while lst[counter].recent_reac != "start":
        numlst.append(counter)
        counter = lst[counter].pathindex
    string = "start"
    returnlist.append((lst[0].molec, "start"))
    for i in reversed(numlst): 
        string = string + ", " + lst[i].recent_reac
        returnlist.append((lst[i].molec, string))
    return returnlist

def make_5carbon_alkene():
    global carbon1, carbon2, carbon3, carbon4, carbon5
    carbon1 = atom("carbon", "carbon1")
    carbon2 = atom("carbon", "carbon2")
    carbon3 = atom("carbon", "carbon3")
    carbon4 = atom("carbon", "carbon4")
    carbon5 = atom("carbon", "carbon5")
    alkene = molecule()
    alkene.add_atom(carbon1)
    alkene.add_atom(carbon2)
    alkene.add_atom(carbon3)
    alkene.add_atom(carbon4)
    alkene.add_atom(carbon5)

    alkene.add_bond((carbon5, carbon4, 1))
    alkene.add_bond((carbon1, carbon2, 2))
    alkene.add_bond((carbon2, carbon3, 2))
    alkene.add_bond((carbon3, carbon4, 2))
    return alkene

def make_5carbon_alcohol():
    global oxygen1, hydrogen1 
    oxygen1 = atom("oxygen", "oxygen1")
    hydrogen1 = atom("hydrogen", "hydrogen1")
    oxygen2 = atom("oxygen", "oxygen2")
    hydrogen2 = atom("hydrogen", "hydrogen2")

    alcohol = make_5carbon_alkene()

    alcohol.replace_bond((carbon1, carbon2, 2), (carbon1, carbon2, 1))
    alcohol.replace_bond((carbon2, carbon3, 2), (carbon2, carbon3, 1))
    alcohol.replace_bond((carbon3, carbon4, 2), (carbon3, carbon4, 1))

    alcohol.add_atom(oxygen1)
    alcohol.add_atom(hydrogen1)
    alcohol.add_bond((oxygen1, hydrogen1, 1))
    alcohol.add_bond((oxygen1, carbon5, 1))

    alcohol.add_atom(oxygen2)
    alcohol.add_atom(hydrogen2)
    alcohol.add_bond((oxygen2, hydrogen2, 1))
    alcohol.add_bond((oxygen2, carbon3, 1))
    return alcohol


def carbon_lst(num):
    lst = []
    for i in range(num):
        lst.append(atom("carbon", "carbon" + str(i)))
    return lst


# https://www.udacity.com/wiki/creating-network-graphs-with-python
def draw_graph(graph, title, labels=None, 
               node_size=1200, node_color='blue', node_alpha=.8,
               node_text_size=12,
               edge_color='blue', edge_alpha=0.3, edge_tickness=2,
               edge_text_pos=0.3,
               text_font='sans-serif'):

    # create networkx graph
    G=nx.Graph()

    # add edges
    for edge in graph:
        G.add_edge(edge[0], edge[1])

    
    # 'bipartite_layout',
    # 'circular_layout',
    # 'kamada_kawai_layout',
    # 'random_layout',
    # 'rescale_layout',
    # 'shell_layout',
    # 'spring_layout',
    # 'spectral_layout',
    # 'planar_layout',
    # 'fruchterman_reingold_layout',
    # 'spiral_layout'
    graph_pos=nx.kamada_kawai_layout(G)
    

    # draw graph
    nx.draw_networkx_nodes(G,graph_pos,node_size=node_size, 
                           alpha=node_alpha, node_color=node_color)
    nx.draw_networkx_edges(G,graph_pos,width=edge_tickness,
                           alpha=edge_alpha,edge_color=edge_color)
    nx.draw_networkx_labels(G, graph_pos,font_size=node_text_size,
                            font_family=text_font)
    #nx.draw_networkx(G, label="title")

    if labels is None:
        labels = range(len(graph))


    edge_labels = dict(zip(graph, labels))
    nx.draw_networkx_edge_labels(G, graph_pos, edge_labels=edge_labels, 
                                 label_pos=edge_text_pos)

    # show graph
    plt.title(title)
    plt.show()

               
def name_products(length, word):
    name_lst = []
    for i in range(length): 
        name_lst.append(prefixes[i] + " " + word + " product")
    return name_lst


def remove_molecule_duplicates(molecule_lst):
    return_list = []
    for i in molecule_lst: 
        inalready = False
        for j in return_list: 
            if test_molecule_isomorphism(i, j):
                inalready = True
        if not inalready: 
            return_list.append(i)
    return return_list

def test_molecule_isomorphism(m1, m2):

    #simple checks
    if len(m1.atoms) == len(m2.atoms):
        atom_lst_1 = []
        atom_lst_2 = []
        for i in m1.atoms: 
            atom_lst_1.append(i.element)
        for i in m2.atoms: 
            atom_lst_2.append(i.element)
        if sorted(atom_lst_1) != sorted(atom_lst_2):
            return False    
    else: 
        return False
    
    counter1 = 0
    counter2 = 0
    for i in m1.atoms: 
        counter1 += (i.bonds)
    for i in m2.atoms: 
        counter2 += (i.bonds)
    if counter1 != counter2:
        return False    


    # try to find isomorph
    for a in range(len(m2.atoms)):
        if find_isomorphism(0, m1, [m1.atoms[0].label], a, m2, [m2.atoms[a].label]):
            return True
    return False

def find_isomorphism(m1start, m1, m1_visited_lst, m2start, m2, m2_visited_lst):
    if m1.atoms[m1start].element == m2.atoms[m2start].element:
        # go through connected atoms 
        for i in m1.atoms[m1start].connected_atoms:
            # if they arent in the visited list, then need to visit 
            if i.label not in m1_visited_lst:
                is_there_a_match = False
                # go through m2 connected atoms to see if any match
                for j in m2.atoms[m2start].connected_atoms: 
                    if j.label not in m2_visited_lst:
                        # if they are the same element 
                        if j.element == i.element:
                            # if they have the right number of bonds
                            correct_num_bonds = False
                            correct_num = 0
                            for k in m1.bonds: 
                                if (k[0] == i and k[1] == m1.atoms[m1start]) or (k[0] == m1.atoms[m1start] and k[1] == i):
                                    correct_num = k[2]
                            for k in m2.bonds: 
                                if (k[0] == j and k[1] == m2.atoms[m2start]) or (k[0] == m2.atoms[m2start] and k[1] == j):
                                    if k[2] == correct_num:
                                        correct_num_bonds = True
                            # if they have the right number of bonds and is the right element, continue checking
                            if correct_num_bonds == True:
                                # update visited list
                                m1_visited_lst.append(i.label)
                                m2_visited_lst.append(j.label)
                                # find indices of i and j to pass in correctly
                                newstart1 = 0
                                newstart2 = 0
                                for k in range(len(m1.atoms)): 
                                    if i.label == m1.atoms[k].label:
                                        newstart1 = k
                                for k in range(len(m2.atoms)):
                                    if j.label == m2.atoms[k].label:
                                        newstart2 = k
                                if (find_isomorphism(newstart1, m1, m1_visited_lst, newstart2, m2, m2_visited_lst)):
                                    is_there_a_match = True
                                m1_visited_lst.pop()
                                #m2_visited_lst.pop()
                if is_there_a_match == False: 
                    return False

        # if no connected atoms or all conected atoms visited, return true 
        return True
    # if starts are different elements 
    return False

