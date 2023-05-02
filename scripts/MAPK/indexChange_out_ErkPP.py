# %%
'''compare the index of random subgraphs on the original network and network with outflow of ErkPP'''
# %%
import random
import numpy as np
import matplotlib.pyplot as plt
from src import ReactionNetwork
# %%
PATH = '../../data/MAPK/MAPKnetwork2.csv'
network1 = ReactionNetwork.from_csv(PATH)

# %%
# make a random subgraph, compare its index on network1/2
random.seed(0)
index_pare = []

cpd_out = 'ErkPP'
outflow = [f'out_{cpd_out}', [cpd_out], ['out']]
reaction_list_out = network1.reaction_list_reg.copy()
reaction_list_out.append(outflow)
network2 = ReactionNetwork.ReactionNetwork(reaction_list_out)

M = network1.M
trials = 1000
for t in range(trials):
    number_cpds = 4
    choose_cpds = random.sample(network1.cpd_list_noout, number_cpds)

    # reactions for network1
    subg_1 = network1.make_ocompSubg(choose_cpds)

    # reactions for network2
    subg_2 = network2.make_ocompSubg(choose_cpds)

    index_pare.append([network1.index_subg(subg_1),
                      network2.index_subg(subg_2)])
index_pare = np.array(index_pare)
# %%
fig = plt.figure(figsize=(8, 8), facecolor='w')
plt.plot([10, -30], [10, -30], color='black')
plt.scatter(index_pare[:, 0], index_pare[:, 1], s=200)
plt.xlabel('Index in the original network', size=20)
plt.ylabel('Index in the network with outflow of ErkPP', size=20)
# plt.title(f'change of index by the outflow to {cpd_out}', size=25)
plt.xlim(-5, 0)
plt.ylim(-5, 0)
plt.xticks(np.arange(-4, 1), np.arange(-4, 1), size=16)
plt.yticks(np.arange(-4, 1), np.arange(-4, 1), size=16)
plt.savefig('../../results/MAPK/indexChange_out_ErkPP.SVG')
plt.show()
# %%
network1.cpd_list_noout[1]
# %%
network2.compute_smat()
