# %%
# make list of reactions in the central metabolic network to draw a network in the cytoscape.
# %%
# buffering structures in the original network
from src import ReactionNetwork
import matplotlib.pyplot as plt
import numpy as np
# %%
network = ReactionNetwork.from_csv(
    '../../../data/central_metabolism/kgml/reaction_list_renamed.csv')
# %%
edge_list = []

# one reaction is converted to more than 2 edges: [sub,reac] & [reac,pro]


def conv_edge(reac):
    edge_list = []
    rname = reac[0]
    for sub in reac[1]:
        if sub == 'out':
            edge = ['outnode_'+rname, rname]
        else:
            edge = [sub, rname]
        edge_list.append(edge)

    for pro in reac[2]:
        if pro == 'out':
            edge = [rname, 'outnode_'+rname]
        else:
            edge = [rname, pro]
        edge_list.append(edge)
    return edge_list


# %%
edge_list = [['source', 'target', 'reversibility']]
for n, reac in enumerate(network.reaction_list):
    # if the reaction is reversible, flag the edge "reversible"
    if [reac[2], reac[1]] in network.reaction_list_noid:
        if network.reaction_list_noid.index([reac[2], reac[1]]) > n:
            for edge in conv_edge(reac):
                edge_list.append(edge+['reversible'])
        else:
            continue
    else:
        edge_list += conv_edge(reac)
# %%
with open('./central_metabolism_cytoscape.csv', 'w') as f:
    for edge in edge_list:
        f.write(','.join(edge))
        f.write('\n')
# %%
node_list = []
for edge in edge_list:
    node_list.append(edge[0])
    node_list.append(edge[1])
node_list = list(set(node_list))
# %%
node_table = []
for node in node_list:
    node_info = [node]
    if node.isdigit():
        node_table.append([node, 'reaction'])
    elif 'outnode' in node:
        node_table.append([node, 'reaction'])
    else:
        node_table.append([node, 'substance'])
