'''sensitivity analysis of MAPK network'''
# %%
from tqdm import tqdm
import csv
import numpy as np
import matplotlib.pyplot as plt
import sys
from src import ReactionNetwork
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)

# %%
PATH = '../../data/MAPK/MAPKnetwork2.csv'
network1 = ReactionNetwork.from_csv(PATH)
# %%
def make_reguMat(network):
    # reaction and its substrates
    reguMat = []
    for cpd in network.cpd_list_noout:
        row = [network.reaction_list_reg[r]
               [1].count(cpd) for r in range(network.R)]
        reguMat.append(row)
    reguMat = np.array(reguMat)
    return reguMat


def activation(x):
    return np.exp((x-0.5)*5)/(1+np.exp((x-0.5)*5))


def inhibition(x):
    return 1/(1+np.exp((x-0.5)*2))


K = 0.5


def f1(k1, RasGDP, ErkPP):
    return k1*RasGDP/(K+RasGDP)*1/(1+ErkPP)


def f2(k2, RasGTP):
    return k2*RasGTP/(K+RasGTP)


def f3(k3, Raf, RasGTP, ErkPP):
    return k3*Raf*RasGTP/(K+Raf) * 1/(1+ErkPP)


def f4(k4, RafP):
    return k4*RafP/(K+RafP)


def f5(k5, Mek, RafP, ErkPP):
    return k5*Mek*RafP/(K+Mek) * 1/(1+ErkPP)


def f6(k6, MekPP):
    return k6*MekPP/(K+MekPP)


def f7(k7, Erk, MekPP):
    return k7*Erk*MekPP/(K+Erk)


def f8(k8, ErkPP):
    return k8*ErkPP/(K+ErkPP)


def f9(k9):  # in
    return k9


def f10(k10, ErkPP):  # out
    return k10*ErkPP/(K+ErkPP)


def compute_flux(network, x, reguMat, params, activate=[], inhibit=[]):
    # activate/inhibit=[(index_substrate, index_reaction),...]
    M = network.M
    R = network.R
    flux = np.empty(R)
    flux[0] = f1(params[0], x[6], x[1])
    flux[1] = f2(params[1], x[7])
    flux[2] = f3(params[2], x[4], x[7], x[1])
    flux[3] = f4(params[3], x[5])
    flux[4] = f5(params[4], x[2], x[5], x[1])
    flux[5] = f6(params[5], x[3])
    flux[6] = f7(params[6], x[0], x[3])
    flux[7] = f8(params[7], x[1])
    if R > 8:
        flux[8] = f9(params[8])
        flux[9] = f10(params[9], x[1])
    return flux


def perturb(network, ini, steps, params, perturbed, dt=0.01):
    # numcalc of the network, petrubation to a specific reaction
    # steps=[N1, N2]
    # perturbed=[reactionName, perburbation]
    M = network.M
    R = network.R
    reguMat = make_reguMat(network)

    for i, reac in enumerate(network.reaction_list):
        if reac[0] == perturbed[0]:
            perturbed_index = i
    params2 = params.copy()
    params2[perturbed_index] += perturbed[1]
    ans = np.zeros((steps[1], M))
    ans[0] = ini

    activate = []
    inhibit = []
    for r, reac in enumerate(network.reaction_list_reg):
        if len(reac) > 3:
            activators = reac[3]
            for cpd in activators:
                m = network.cpd_list_noout.index(cpd)
                activate.append((m, r))
        if len(reac) > 4:
            inhibitors = reac[4]
            for cpd in inhibitors:
                m = network.cpd_list_noout.index(cpd)
                inhibit.append((m, r))
    # print(activate,inhibit)

    for n in tqdm(range(1, steps[1])):
        if n < steps[0]+1:
            flux = compute_flux(network, ans[n-1], reguMat, params,
                                activate=activate, inhibit=inhibit)
        else:
            flux = compute_flux(network, ans[n-1], reguMat, params2,
                                activate=activate, inhibit=inhibit)

        ans[n] = ans[n-1]+np.dot(network.stoi, flux)*dt
        if np.max(ans[n]) > 1.0e20:
            print(n)
            print('overflow')
            break
    return ans
# %%


def plot_4cpds_mapk(ans1, dt, cpds, network, yrange=(), perturb_name='', savepath=False, start=7000, end=13000):
    triangle = plt.imread('../../data/triangle.png')
    # triangle=Image.open('../../../data/triangle.png')
    imagebox = OffsetImage(triangle, zoom=0.02)
    # fumarate;18, Oxaloacetate31 , D-glucose45;, deoxyribose12
    M = network.M
    cpd_index = [network.cpd_list_noout.index(cpd) for cpd in cpds]
    fig, ax = plt.subplots(1, 4, figsize=(
        16, 4), tight_layout=True, facecolor='w')

    for i, subax in enumerate(ax):
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
        if yrange:
            ymin = yrange[0][i]
            ymax = yrange[1][i]
            subax.set_ylim(ymin, ymax)
            subax.set_yticks([ymin, (ymin+ymax)/2, ymax])
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


def plot_mapk(ans, dt, network, perturbed, savepath=False):
    X = 4
    Y = 2
    fig, ax = plt.subplots(Y, X, figsize=(
        15, 8), facecolor='w', tight_layout=True)
    ymax = np.max(ans)
    start = 7000
    end = 13000
    ax_1d = ax.reshape(-1)
    for i, subax in enumerate(ax_1d):
        if i > network.M-1:
            subax.axi('off')
            continue
        if i % X == 0:
            subax.set_ylabel('concentration', size=25)
        if i % X == 3:
            subax.set_xlabel('time', loc='right', size=30)

        ymin = np.min(ans[start:end, i])*0.8
        ymax = np.max(ans[start:end, i])*1.2
        subax.plot(np.arange(end)[start:end]*dt,
                   ans[start:end, i], linewidth=5.0)
        subax.set_ylim(ymin, ymax)
        subax.set_yticks([ymin, (ymin+ymax)/2, ymax])
        subax.set_title(network.cpd_list_noout[i], fontsize=28)
        subax.set_xticks([start*dt, (start+end)*dt/2, end*dt])
        subax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
        subax.tick_params(axis='y', labelsize=20)
        subax.tick_params(axis='x', labelsize=20)
    # fig.suptitle(f' perturbation on {perturbed[0]}',fontsize=25)
    if savepath:
        plt.savefig(savepath)
    plt.show()


# %%
params = np.full(network1.R, 3.)
# params=np.array()
params = np.array([3., 3., 3., 3., 3., 1.5, 3., 1.5])
perturbed_reaction = '1'
np.random.seed(seed=0)
ini = np.random.rand(network1.M)
N_1 = 10000
N_2 = 20000
dt = 0.01

# %%
steps = [N_1, N_2]
perturbed = (perturbed_reaction, 5.)
ans1 = perturb(network1, ini, steps, params, perturbed)

# %%
cpds = ['Ras-GTP', 'RafP', 'MekPP', 'ErkPP']
ymin = [0.0, 0.0, 0.0, 0.0]
ymax = [1.5, 0.6, 0.7, 0.5]
yrange = (ymin, ymax)
savepath = f'../../results/MAPK/numcalc/numcalc_perturb_{perturbed_reaction}.SVG'
plot_4cpds_mapk(ans1, dt, cpds, network1, yrange, savepath=savepath)
savepath = f'../../results/MAPK/numcalc/numcalc_perturb_{perturbed_reaction}.png'
plot_4cpds_mapk(ans1, dt, cpds, network1, yrange, savepath=savepath)
# %%

savepath = f'../../results/MAPK/numcalc/numcalc_perturb_{perturbed_reaction}_all.SVG'
plot_mapk(ans1, dt, network1, perturbed, savepath)
savepath = f'../../results/MAPK/numcalc/numcalc_perturb_{perturbed_reaction}_all.png'
plot_mapk(ans1, dt, network1, perturbed, savepath)
# %%
'''
add outflow to network
influx is needed
outflow to ErkPP, does signal network work?
'''
out_cpd = 'ErkPP'
in_cpd = 'Erk'
reaction_list2 = network1.reaction_list_reg.copy()
reaction_list2.append([f'inflow_{in_cpd}', ['out'], [in_cpd]])
reaction_list2.append([f'outflow_{out_cpd}', [out_cpd], ['out']])
network2 = ReactionNetwork.ReactionNetwork(reaction_list2)
# %%
perturbed_reaction = '1'
np.random.seed(seed=0)
ini = np.random.rand(network2.M)
params = np.array([3., 3., 3., 3., 3., 1.5, 3., 1.5, 0.2, 2.])
# ini[:2]+=1.0
N_1 = 10000
N_2 = 20000
dt = 0.01

# %%
steps = [N_1, N_2]
perturbed = (perturbed_reaction, 3.)
ans2 = perturb(network2, ini, steps, params, perturbed)
# %%
cpds = ['Ras-GTP', 'RafP', 'MekPP', 'ErkPP']
start = 9000
end = 11000
savepath = f'../../results/MAPK/numcalc/numcalc_perturb_{perturbed_reaction}_out_ErkPP.SVG'
plot_4cpds_mapk(ans2, dt, cpds, network2, yrange, savepath=savepath)
savepath = f'../../results/MAPK/numcalc/numcalc_perturb_{perturbed_reaction}_out_ErkPP.png'
plot_4cpds_mapk(ans2, dt, cpds, network2, yrange, savepath=savepath)

# %%
savepath = f'../../results/MAPK/numcalc/numcalc_perturb_{perturbed_reaction}_out_ErkPP_all.SVG'
plot_mapk(ans2, dt, network2, perturbed, savepath)
savepath = f'../../results/MAPK/numcalc/numcalc_perturb_{perturbed_reaction}_out_ErkPP_all.png'
plot_mapk(ans2, dt, network2, perturbed, savepath)
