'''sensitivity of the original network'''
# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from src.plot_numcalc import plot_allcpd_tca
from src.plot_numcalc import plot_4cpds_tca
from src import ReactionNetwork
from src import massaction
# %%
# network from renamed csv
network = ReactionNetwork.from_csv(
    '../../../data/central_metabolism/kgml/reaction_list_renamed.csv', info=True)

# %%
# perturb Oxa -> PEP (rn:R00431)
M = network.M
R = network.R
params = np.full(R, 3.0)
perturbed = '13'
# perturbed='rn:R00431_rn:R00726'
for i, reac in enumerate(network.reaction_list):
    if reac[0] == perturbed:
        perturbed_index = i
params2 = params.copy()
params2[perturbed_index] += 3.
np.random.seed(seed=1)
ini = np.random.rand(M)
N_1 = 40000
N_2 = 80000
dt = 0.01

ans1 = massaction.perturb(network, ini=ini, steps=[
                          N_1, N_2], params=params, perturbed=[perturbed, 3.])
# %%
savepath = '../../../results/central_metabolism/numcalc/numcalc_perturb_Oxa_PEP_allcpds.SVG'
plot_allcpd_tca(ans1, dt, network, [perturbed, 3.0], savepath=savepath)
savepath = '../../../results/central_metabolism/numcalc/numcalc_perturb_Oxa_PEP_allcpds.png'
plot_allcpd_tca(ans1, dt, network, [perturbed, 3.0], savepath=savepath)

# plot_allcpd_tca(ans1,network,perturbed)
# %%
cpds = ['Deoxyribose', 'PEP', 'Fum', '2-Oxo']
ymin_list = [6.34, 2.3, 0.0, 0.0]
ymax_list = [6.36, 2.7, 2.0, 0.4]
perturb_name = 'Oxa -> PEP'
savepath = '../../../results/central_metabolism/numcalc/numcalc_perturb_Oxa_PEP_4cpds.SVG'
plot_4cpds_tca(ans1, dt, cpds, network, (ymin_list,
               ymax_list), perturb_name, savepath)
savepath = '../../../results/central_metabolism/numcalc/numcalc_perturb_Oxa_PEP_4cpds.png'
plot_4cpds_tca(ans1, dt, cpds, network, (ymin_list,
               ymax_list), perturb_name, savepath)
# %%
M = network.M
R = network.R
params = np.full(R, 3.0)
perturb_name = 'Oxa -> Cit'
perturbed = '9'
# perturbed='rn:R00352'
for i, reac in enumerate(network.reaction_list):
    if reac[0] == perturbed:
        perturbed_index = i
np.random.seed(seed=1)
ini = np.random.rand(M)
N_1 = 40000
N_2 = 80000
dt = 0.01

ans1 = massaction.perturb(network, ini=ini, steps=[N_1, N_2],
                          params=params, perturbed=[perturbed, 0.6])

# plot_allcpd_tca(ans1,network,perturbed)
# %%
savepath = '../../../results/central_metabolism/numcalc/numcalc_perturb_Oxa_Cit_allcpds.SVG'
plot_allcpd_tca(ans1, dt, network, [perturbed, 3.0], savepath=savepath)
savepath = '../../../results/central_metabolism/numcalc/numcalc_perturb_Oxa_Cit_allcpds.png'
plot_allcpd_tca(ans1, dt, network, [perturbed, 3.0], savepath=savepath)
# %%
cpds = ['Deoxyribose', 'PEP', 'Fum', '2-Oxo']
ymin_list = [6.33, 2.30, 1.5, 0.29]
ymax_list = [6.35, 2.45, 1.8, 0.36]
savepath = '../../../results/central_metabolism/numcalc/numcalc_perturb_Oxa_Cit_4cpds.SVG'
plot_4cpds_tca(ans1, dt, cpds, network, (ymin_list,
               ymax_list), perturb_name, savepath)
savepath = '../../../results/central_metabolism/numcalc/numcalc_perturb_Oxa_Cit_4cpds.png'
plot_4cpds_tca(ans1, dt, cpds, network, (ymin_list,
               ymax_list), perturb_name, savepath)
# %%
# perturb Fum -> Mal (rn:R01082_2)
M = network.M
R = network.R
params = np.full(R, 3.0)
perturbed = '43'
# perturbed='rn:R01082_2'
perturb_name = 'Fum -> Mal'
for i, reac in enumerate(network.reaction_list):
    if reac[0] == perturbed:
        perturbed_index = i
np.random.seed(seed=1)
ini = np.random.rand(M)
N_1 = 40000
N_2 = 80000
dt = 0.01

ans1 = massaction.perturb(network, ini=ini, steps=[N_1, N_2],
                          params=params, perturbed=[perturbed, 3.0])

# plot_allcpd_tca(ans1,network,perturbed)
# %%
savepath = '../../../results/central_metabolism/numcalc/numcalc_perturb_Fum_Mal_allcpds.SVG'
plot_allcpd_tca(ans1, dt, network, [perturbed, 3.0], savepath=savepath)
savepath = '../../../results/central_metabolism/numcalc/numcalc_perturb_Fum_Mal_allcpds.png'
plot_allcpd_tca(ans1, dt, network, [perturbed, 3.0], savepath=savepath)
# %%
cpds = ['Deoxyribose', 'PEP', 'Fum', '2-Oxo']
ymin_list = [6.32, 2.3, 0.6, 0.29]
ymax_list = [6.36, 2.6, 1.8, 0.35]
savepath = '../../../results/central_metabolism/numcalc/numcalc_perturb_Fum_Mal_4cpds.SVG'
plot_4cpds_tca(ans1, dt, cpds, network, (ymin_list,
               ymax_list), perturb_name, savepath)
savepath = '../../../results/central_metabolism/numcalc/numcalc_perturb_Fum_Mal_4cpds.png'
plot_4cpds_tca(ans1, dt, cpds, network, (ymin_list,
               ymax_list), perturb_name, savepath)

# %%
# perturb Fum -> Oxa (rn:R01082_2)
M = network.M
R = network.R
params = np.full(R, 3.0)
perturbed = '38'
# perturbed='rn:R01066'
perturb_name = 'DR5P -> G3P'
for i, reac in enumerate(network.reaction_list):
    if reac[0] == perturbed:
        perturbed_index = i
np.random.seed(seed=1)
ini = np.random.rand(M)
N_1 = 40000
N_2 = 80000
dt = 0.01

ans1 = massaction.perturb(network, ini=ini, steps=[N_1, N_2],
                          params=params, perturbed=[perturbed, 3.0])

# plot_allcpd_tca(ans1,network,perturbed)
# %%
savepath = '../../../results/central_metabolism/numcalc/numcalc_perturb_DR5P_G3P_allcpds.SVG'
plot_allcpd_tca(ans1, dt, network, [perturbed, 3.0], savepath=savepath)
savepath = '../../../results/central_metabolism/numcalc/numcalc_perturb_DR5P_G3P_allcpds.png'
plot_allcpd_tca(ans1, dt, network, [perturbed, 3.0], savepath=savepath)
# %%
cpds = ['Deoxyribose', 'PEP', 'Fum', '2-Oxo']
ymin_list = [2.0, 2.3, 1.5, 0.25]
ymax_list = [7.0, 2.5, 1.7, 0.40]
savepath = '../../../results/central_metabolism/numcalc/numcalc_perturb_DR5P_G3P_4cpds.SVG'
plot_4cpds_tca(ans1, dt, cpds, network, (ymin_list,
               ymax_list), perturb_name, savepath)
savepath = '../../../results/central_metabolism/numcalc/numcalc_perturb_DR5P_G3P_4cpds.png'
plot_4cpds_tca(ans1, dt, cpds, network, (ymin_list,
               ymax_list), perturb_name, savepath)
# %%
