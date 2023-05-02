# %%
# demo of the effect to sensitivity by addition of new outflux
# %%
import numpy as np
import matplotlib.pyplot as plt
from src import ReactionNetwork
from src import massaction
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
# %%
def plot_demo(ans1, dt, cpds, network, yrange=(), perturb_name='', savepath=False, start=5000, end=15000):
    triangle = plt.imread('../../data/triangle.png')
    # triangle=Image.open('../../../data/triangle.png')
    imagebox = OffsetImage(triangle, zoom=0.02)
    # fumarate;18, Oxaloacetate31 , D-glucose45;, deoxyribose12
    M = network.M
    cpd_index = [network.cpd_list_noout.index(cpd) for cpd in cpds]
    fig, ax = plt.subplots(1, 4, figsize=(
        16, 4), tight_layout=True, facecolor='w')
    ax_1d = ax.reshape(-1)
    for i, subax in enumerate(ax_1d):

        ymin = np.min(ans1[start:end, i])*0.8
        ymax = np.max(ans1[start:end, i])*1.2

        if i == 0:
            subax.set_ylabel('concentration', size=25)
        if i == 3:
            subax.set_xlabel('time', loc='right', size=30)

        # plot the answer
        subax.plot(np.arange(end)[start:end]*dt,
                   ans1[start:end, cpd_index[i]], linewidth=5.0)

        # equib value
        equib_before = ans1[int((start+end)/2), cpd_index[i]]
        equib_after = ans1[end, cpd_index[i]]
        subax.plot(np.array(
            [start, end])*dt, np.array([equib_before, equib_before]), linestyle='dashed')
        subax.plot(np.array(
            [start, end])*dt, np.array([equib_after, equib_after]), linestyle='dashed')

        # if equib value changes, set spine color red
        if abs(equib_before - equib_after) > 1.0e-3:
            subax.spines['left'].set_color('red')
            subax.spines['right'].set_color('red')
            subax.spines['top'].set_color('red')
            subax.spines['bottom'].set_color('red')
            subax.spines['left'].set_linewidth(3)
            subax.spines['right'].set_linewidth(3)
            subax.spines['top'].set_linewidth(3)
            subax.spines['bottom'].set_linewidth(3)

        # label ticks
        subax.set_xticks(np.array([start*dt, (start+end)*dt/2, end*dt]))
        subax.set_ylim(ymin, ymax)
        yticks_str = ['{:.2f}'.format(x) for x in [ymin, (ymin+ymax)/2, ymax]]
        subax.set_yticks([ymin, (ymin+ymax)/2, ymax], yticks_str)
        subax.tick_params(axis='y', labelsize=20)
        subax.tick_params(axis='x', labelsize=20)

        # title
        subax.set_title(network.cpd_list_noout[cpd_index[i]], fontsize=28)

        # triangle
        ab = AnnotationBbox(imagebox, ((start+end)*dt/2, ymin),
                            xybox=(0, -40.),
                            xycoords='data',
                            boxcoords="offset points", frameon=False)
        subax.add_artist(ab)

#     fig.suptitle(f'perturbation on {perturb_name}',fontsize=30)
    if savepath:
        plt.savefig(savepath)


# %%
# load the original network
network1 = ReactionNetwork.from_csv('./demoNetwork.csv')
# %%
# numcalc of original network
ini = np.random.rand(network1.M)
N1 = 10000
N2 = 20000
perturbed = '6'
perturbation = 6.0
params = np.full(network1.R, 2.5)
dt = 0.01
ans = massaction.perturb(network1, ini, (N1, N2), params,
                         (perturbed, perturbation), dt=0.01)
# %%
# plot the result
cpds = network1.cpd_list_noout
savepath = '../../results/demo_network/numcalc_demo_6.png'
plot_demo(ans, dt, cpds, network1, savepath=savepath)
savepath = '../../results/demo_network/numcalc_demo_6.svg'
plot_demo(ans, dt, cpds, network1, savepath=savepath)


# %%
# add an outflux to D
reaction_list2 = network1.reaction_list+[['outflux', ['D'], ['out']]]
network2 = ReactionNetwork.ReactionNetwork(reaction_list2)
# %%
# numcalc of new network
ini = np.random.rand(network2.M)
N1 = 10000
N2 = 20000
perturbed = '6'
perturbation = 6.0
params = np.full(network2.R, 2.5)
dt = 0.01
ans2 = massaction.perturb(network2, ini, (N1, N2),
                          params, (perturbed, perturbation), dt=0.01)
# %%
# plot the result
cpds = network2.cpd_list_noout
savepath = '../../results/demo_network/numcalc_demo_6_out_D.png'
plot_demo(ans2, dt, cpds, network2, savepath=savepath)
savepath = '../../results/demo_network/numcalc_demo_6_out_D.svg'
plot_demo(ans2, dt, cpds, network2, savepath=savepath)

# %%
