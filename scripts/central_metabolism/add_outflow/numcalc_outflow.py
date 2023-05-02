'''sensitivity of the network with outflow of 2-oxo'''
# %%
import numpy as np
from src import ReactionNetwork
from src import massaction
from src.plot_numcalc import plot_4cpds_tca
from src.plot_numcalc import plot_allcpd_tca

# %%
# network from renamed csv
network = ReactionNetwork.from_csv(
    '../../../data/central_metabolism/kgml/reaction_list_renamed.csv', info=True)
# new network with outflow of 2-Oxo
target_cpd = '2-Oxo'
outflow = [f'out_{target_cpd}', [target_cpd], ['out']]
reaction_list_new = network.reaction_list+[outflow]
network_new = ReactionNetwork.ReactionNetwork(reaction_list_new, info=True)

# %%
# perturb Oxa -> PEP (rn:R00431)
M = network_new.M
R = network_new.R
params = np.full(R, 3.0)
perturbed = '13'
for i, reac in enumerate(network_new.reaction_list):
    if reac[0] == perturbed:
        perturbed_index = i
params2 = params.copy()
params2[perturbed_index] += 3.
np.random.seed(seed=1)
ini = np.random.rand(M)
N_1 = 40000
N_2 = 80000
dt = 0.01

ans1 = massaction.perturb(network_new, ini=ini, steps=[
                          N_1, N_2], params=params, perturbed=[perturbed, 3.])

# %%
# plot and save the all results
savepath = f'../../../results/central_metabolism/numcalc/numcalc_perturb_Oxa_PEP_allcpds_out_{target_cpd}.SVG'
plot_allcpd_tca(ans1, dt, network_new, [perturbed, 3.0], savepath=savepath)
savepath = f'../../../results/central_metabolism/numcalc/numcalc_perturb_Oxa_PEP_allcpds_out_{target_cpd}.png'
plot_allcpd_tca(ans1, dt, network_new, [perturbed, 3.0], savepath=savepath)

# %%
# plot only 4 compounds
cpds = ['Deoxyribose', 'PEP', 'Fum', '2-Oxo']
ymin_list = [6.30, 2.2, 0.5, 0.07]
ymax_list = [6.34, 2.7, 1.5, 0.15]
perturb_name = f'Oxa -> PEP, out_{target_cpd}'
savepath = f'../../../results/central_metabolism/numcalc/numcalc_perturb_Oxa_PEP_4cpds_out_{target_cpd}.SVG'
plot_4cpds_tca(ans1, dt, cpds, network_new,
               (ymin_list, ymax_list), perturb_name, savepath)
savepath = f'../../../results/central_metabolism/numcalc/numcalc_perturb_Oxa_PEP_4cpds_out_{target_cpd}.png'
plot_4cpds_tca(ans1, dt, cpds, network_new,
               (ymin_list, ymax_list), perturb_name, savepath)
# plot_4cpds_tca(ans1,dt,cpds,network_new)
# %%
# perturb inflow of a-Glu (in_Glu)
M = network_new.M
R = network_new.R
params = np.full(R, 3.0)
perturb_name = 'Oxa -> Cit'
perturbed = '9'
# perturbed = 'rn:R00352'
for i, reac in enumerate(network_new.reaction_list):
    if reac[0] == perturbed:
        perturbed_index = i
np.random.seed(seed=1)
ini = np.random.rand(M)
N_1 = 40000
N_2 = 80000
dt = 0.01

ans1 = massaction.perturb(network_new, ini=ini, steps=[N_1, N_2],
                          params=params, perturbed=[perturbed, 0.6])


# %%
savepath = f'../../../results/central_metabolism/numcalc/numcalc_perturb_Oxa_Cit_allcpds_out_{target_cpd}.SVG'
plot_allcpd_tca(ans1, dt, network_new, [perturbed, 3.0], savepath=savepath)
savepath = f'../../../results/central_metabolism/numcalc/numcalc_perturb_Oxa_Cit_allcpds_out_{target_cpd}.png'
plot_allcpd_tca(ans1, dt, network_new, [perturbed, 3.0], savepath=savepath)
# %%
cpds = ['Deoxyribose', 'PEP', 'Fum', '2-Oxo']
ymin_list = [6.31, 2.2, 1.12, 0.135]
ymax_list = [6.32, 2.3, 1.17, 0.170]
savepath = f'../../../results/central_metabolism/numcalc/numcalc_perturb_Oxa_Cit_4cpds_out_{target_cpd}.SVG'
plot_4cpds_tca(ans1, dt, cpds, network_new,
               (ymin_list, ymax_list), perturb_name, savepath)
savepath = f'../../../results/central_metabolism/numcalc/numcalc_perturb_Oxa_Cit_4cpds_out_{target_cpd}.png'
plot_4cpds_tca(ans1, dt, cpds, network_new,
               (ymin_list, ymax_list), perturb_name, savepath)
# plot_4cpds_tca(ans1,dt,cpds,network_new)

# %%
# perturb Fum -> Mal (rn:R01082_2)
M = network_new.M
R = network_new.R
params = np.full(R, 3.0)
perturbed = '43'
# perturbed = 'rn:R01082_2'
perturb_name = f'Fum -> Mal, out_{target_cpd}'
for i, reac in enumerate(network_new.reaction_list):
    if reac[0] == perturbed:
        perturbed_index = i
np.random.seed(seed=1)
ini = np.random.rand(M)
N_1 = 40000
N_2 = 80000
dt = 0.01

ans1 = massaction.perturb(network_new, ini=ini, steps=[N_1, N_2],
                          params=params, perturbed=[perturbed, 3.0])

# plot_allcpd_tca(ans1,network_new,perturbed)
# %%
savepath = f'../../../results/central_metabolism/numcalc/numcalc_perturb_Fum_Mal_allcpds_out_{target_cpd}.SVG'
plot_allcpd_tca(ans1, dt, network_new, [perturbed, 3.0], savepath=savepath)
savepath = f'../../../results/central_metabolism/numcalc/numcalc_perturb_Fum_Mal_allcpds_out_{target_cpd}.png'
plot_allcpd_tca(ans1, dt, network_new, [perturbed, 3.0], savepath=savepath)
# %%
cpds = ['Deoxyribose', 'PEP', 'Fum', '2-Oxo']
ymin_list = [6.30, 2.2, 0.4, 0.13]
ymax_list = [6.33, 2.4, 1.2, 0.16]
savepath = f'../../../results/central_metabolism/numcalc/numcalc_perturb_Fum_Mal_4cpds_out_{target_cpd}.SVG'
plot_4cpds_tca(ans1, dt, cpds, network_new,
               (ymin_list, ymax_list), perturb_name, savepath)
savepath = f'../../../results/central_metabolism/numcalc/numcalc_perturb_Fum_Mal_4cpds_out_{target_cpd}.png'
plot_4cpds_tca(ans1, dt, cpds, network_new,
               (ymin_list, ymax_list), perturb_name, savepath)
# plot_4cpds_tca(ans1,dt,cpds,network_new)

# %%
# perturb Fum -> Mal (rn:R01082_2)
M = network_new.M
R = network_new.R
params = np.full(R, 3.0)
perturbed = '38'
# perturbed = 'rn:R01066'
perturb_name = f'DR5P -> G3P, out_{target_cpd}'
for i, reac in enumerate(network_new.reaction_list):
    if reac[0] == perturbed:
        perturbed_index = i
np.random.seed(seed=1)
ini = np.random.rand(M)
N_1 = 40000
N_2 = 80000
dt = 0.01

ans1 = massaction.perturb(network_new, ini=ini, steps=[N_1, N_2],
                          params=params, perturbed=[perturbed, 3.0])

# plot_allcpd_tca(ans1,network_new,perturbed)
# %%
savepath = f'../../../results/central_metabolism/numcalc/numcalc_perturb_DR5P_G3P_allcpds_out_{target_cpd}.SVG'
plot_allcpd_tca(ans1, dt, network_new, [perturbed, 3.0], savepath=savepath)
savepath = f'../../../results/central_metabolism/numcalc/numcalc_perturb_DR5P_G3P_allcpds_out_{target_cpd}.png'
plot_allcpd_tca(ans1, dt, network_new, [perturbed, 3.0], savepath=savepath)

# %%

cpds = ['Deoxyribose', 'PEP', 'Fum', '2-Oxo']
ymin_list = [2.5, 2.26, 1.10, 0.13]
ymax_list = [7.0, 2.29, 1.18, 0.17]
savepath = f'../../../results/central_metabolism/numcalc/numcalc_perturb_DR5P_G3P_4cpds_out_{target_cpd}.SVG'
plot_4cpds_tca(ans1, dt, cpds, network_new,
               (ymin_list, ymax_list), perturb_name, savepath)
savepath = f'../../../results/central_metabolism/numcalc/numcalc_perturb_DR5P_G3P_4cpds_out_{target_cpd}.png'
plot_4cpds_tca(ans1, dt, cpds, network_new,
               (ymin_list, ymax_list), perturb_name, savepath)
# plot_4cpds_tca(ans1,dt,cpds,network_new)
