"""add an outflow of one compound in the central metabolic network, and  check the index change of TCA cycle"""
# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from src import ReactionNetwork

# %%
# network from renamed csv
network = ReactionNetwork.from_csv(
    '../../../data/central_metabolism/kgml/reaction_list_renamed.csv', info=True)


# %%
cpds = ['2-Oxo', 'AcCoA', 'CisAco', 'Citrate', 'Fum', 'Hco3-', 'Isocitrate',
        'Lactate', 'Malate', 'OxaloSuc', 'Pyruvate', 'Suc', 'SucCoA']
subn = network.make_ocompSubg(cpds)
# %%
network.index_subg(subn)

# %%
# add outflow for a compound in the network
for cpd in network.cpd_list_noout:
    if [[cpd], ['out']] in network.reaction_list_noid:
        print(f'outflow of {cpd} already exists.')
        continue
    
    reaction_list_out = network.reaction_list.copy()
    outflow = [f'newoutflow_{cpd}', [cpd], ['out']]
    reaction_list_out.append(outflow)
    
    # index of TCA in the network with new outflow
    network_out = ReactionNetwork.ReactionNetwork(
        reaction_list_out, info=False)
    subg = network_out.make_ocompSubg(cpds)
    index = network_out.index_subg(subg)

    # if index is not zero, TCA becomes not a buffeirn
    print(cpd, index)


# %%
