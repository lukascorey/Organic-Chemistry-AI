#! usr/bin/env/python3

from classes import atom, molecule, molecule_path_node

from helpers import molecule_to_drawgraph

from helpers import test_molecule_isomorphism, print_path, get_path

from copy import deepcopy
# reaction dict format is {name_of_reaction : function}
def brute_force(start, goal, max_reactions, reaction_dict):
    print_list = []
    print_list.append((start, "starting molecule"))
    print_list.append((goal, "goal molecule"))

    counter = 0
    start_node = molecule_path_node(start, 0, 0, "start")
    lst = [start_node]
    returnv = False
    while (lst[counter].reactions < max_reactions and (not (test_molecule_isomorphism((lst[counter].molec), goal)))): 
        #x = deepcopy(lst[counter].molec)
        for j in reaction_dict:
            # apply every possible reaction
            x = deepcopy(lst[counter].molec)
            products = reaction_dict[j](x)
            product_nodes = []
            for i in products:
                product_nodes.append(molecule_path_node(i, counter, lst[counter].reactions + 1, j)) 
            lst = lst + product_nodes
        counter = counter + 1
    if test_molecule_isomorphism(lst[counter].molec, goal):
        print_list = print_list + (get_path(lst, counter))
        print(lst[counter].molec.bonds) # here
        print(goal)
        returnv = True
    else: 
        while counter < len(lst):
            if test_molecule_isomorphism(lst[counter].molec, goal):
                print_list = print_list + (get_path(lst, counter))
                print(lst[counter].molec.bonds) # here 
                print(goal)
                returnv = True
                break
            else: 
                counter = counter + 1
    #print(test_molecule_isomorphism(lst[counter].molec, goal))
    print(returnv)
    return print_list
        
