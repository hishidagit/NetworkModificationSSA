# %%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
# functions to plot numerical calculation results

# %%
def plot_mapk(ans, dt, network, perturbed, savepath=False):
    X = 4
    Y = 2
    fig, ax = plt.subplots(Y, X, figsize=(60, 40), facecolor='w')
    ymax = np.max(ans)
    start = 20000
    end = 60000
    ax_1d = ax.reshape(-1)
    for i, subax in enumerate(ax_1d):
        if i > network.M-1:
            subax.axi('off')
            continue
        ymin = np.min(ans[start:end, i])*0.8
        ymax = np.max(ans[start:end, i])*1.2
        subax.plot(np.arange(end)[start:end]*dt,
                   ans[start:end, i], linewidth=5.0)
        subax.set_ylim(ymin, ymax)
        subax.set_yticks([ymin, (ymin+ymax)/2, ymax])
        subax.set_title(network.cpd_list_noout[i], fontsize=40)

    fig.suptitle(f' perturbation on {perturbed[0]}', fontsize=50)
    if savepath:
        plt.savefig(savepath)
    plt.show()


def plot_4cpds_tca(ans1, dt, cpds, network, yrange=(), perturb_name='', savepath=False, start=20000, end=60000):
    triangle = plt.imread('../../../data/triangle.png')
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
            # set yticks as x.xxx
            yticks_str = ['{:.2f}'.format(x)
                          for x in [ymin, (ymin+ymax)/2, ymax]]
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


def plot_allcpd_tca(ans, dt, network, perturbed, savepath=False):
    X = 8
    Y = 6
    fig, ax = plt.subplots(Y, X, figsize=(
        30, 20), facecolor='w', tight_layout=True)
    ymax = np.max(ans)
    start = 20000
    end = 60000
    ax_1d = ax.reshape(-1)
    for i, subax in enumerate(ax_1d):
        if i > network.M-1:
            subax.axi('off')
            continue
        if i % 8 == 0:
            subax.set_ylabel('concentration', size=25)
        if i % 8 == 7:
            subax.set_xlabel('time', loc='right', size=30)
        ymin = np.min(ans[start:end, i])*0.8
        ymax = np.max(ans[start:end, i])*1.2

        # plot answer
        subax.plot(np.arange(end)[start:end]*dt,
                   ans[start:end, i], linewidth=5.0)

        # subtitle
        subax.set_title(network.cpd_list_noout[i], fontsize=30)

        # label ticks
        subax.set_xticks(np.array([start*dt, (start+end)*dt/2, end*dt]))
        subax.set_ylim(ymin, ymax)
        # set yticks as x.xxx
        yticks_str = ['{:.2f}'.format(x) for x in [ymin, (ymin+ymax)/2, ymax]]
        subax.set_yticks([ymin, (ymin+ymax)/2, ymax], yticks_str)
        subax.tick_params(axis='y', labelsize=20)
        subax.tick_params(axis='x', labelsize=20)

    # fig.suptitle(f' perturbation on {perturbed[0]}',fontsize=50)
    if savepath:
        plt.savefig(savepath)
    plt.show()
