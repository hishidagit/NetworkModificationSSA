'''make a random subgraph, compare its index on the original network and the network with outflow of 2-Oxo'''
# %%
import numpy as np
import matplotlib.pyplot as plt
from src import ReactionNetwork
import random
# %%
# network from renamed csv
network1 = ReactionNetwork.from_csv(
    '../../../data/central_metabolism/kgml/reaction_list_renamed.csv', info=True)
reaction_list = network1.reaction_list
# %%
# list of index before and after adding outflow
index_pare = []

# make a network with new outflow
cpd_out = '2-Oxo'  # add an outflow to a cpd
if [[cpd_out], ['out']] in network1.reaction_list_noid:
    print('out flux already exists; ', cpd_out)
reaction_list_out = reaction_list.copy()
outflow = ['outflow', [cpd_out], ['out']]
reaction_list_out.append(outflow)
network2 = ReactionNetwork.ReactionNetwork(reaction_list_out)

# number of metabolites
M = network1.M

# select a random subgraph and compare its index on the original network and the network with outflow of 2-Oxo
trials = 1000
for t in range(trials):
    # choose metabolites
    number_cpds = 10
    choose_cpds_index = random.sample(
        list(range(M)), number_cpds)  # random sampling
    choose_cpds = [network1.cpd_list_noout[i] for i in choose_cpds_index]

    # reactions in subgraph of the origina  
    choose_reacs_1 = []
    for r in range(network1.R):
        if not set(network1.reaction_list_noid[r][0]).isdisjoint(choose_cpds):
            choose_reacs_1.append(r)

    # reactions in subgraph of the network with outflow of 2-Oxo
    choose_reacs_2 = []
    for r in range(network2.R):
        if not set(network2.reaction_list_noid[r][0]).isdisjoint(choose_cpds):
            choose_reacs_2.append(r)

    subg_1 = network1.make_ocompSubg(choose_cpds)
    subg_2 = network2.make_ocompSubg(choose_cpds)
    index_pare.append([network1.index_subg(subg_1),
                      network2.index_subg(subg_2)])
index_pare = np.array(index_pare)
# %%
fig = plt.figure(figsize=(8, 8), facecolor='w', tight_layout=True)

# plt.scatter(np.arange(-30,0,0.1),np.arange(-30,0,0.1),s=)
plt.plot([10, -30], [10, -30], color='black')
plt.scatter(index_pare[:, 0], index_pare[:, 1])
plt.xlabel('Index in the original network', size=20)
plt.ylabel('Index in the network with outflow of 2-Oxo', size=20)
# plt.title(f'change of index by the outflow to {cpd_out}', size=25)
plt.xlim(-22, 2)
plt.ylim(-22, 2)
plt.xticks(np.array([-20, -15, -10, -5, 0]),
           np.array([-20, -15, -10, -5, 0]), size=16)
plt.yticks(np.array([-20, -15, -10, -5, 0]),
           np.array([-20, -15, -10, -5, 0]), size=16)
plt.savefig(
    f'../../../results/central_metabolism/indexChange_out_{cpd_out}.svg')
plt.savefig(
    f'../../../results/central_metabolism/indexChange_out_{cpd_out}.png')
plt.show()
# %%
