import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def get_plot(result, mtype):
    
    mpl.rcParams['text.usetex'] = True
    mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

    fontsize = 15
    labelsize = 18
    markersize = 8

    fig, ax = plt.subplots(1, 3, figsize = (17, 5), constrained_layout = True)




    pars = {'m': 5}
    index = 0
    result1 = result
    result1 = result1.loc[result1['m'] == 5]
    result1 = result1.loc[result1['s'] == 5]

    types = ['sparse']
    mks = ['x']

    estimators = ['adele', 'mrlasso']
    cols = ['orange', 'blue']
    mks = ['x', 'o']

    if mtype == "l2.":
        ylab = r'$\|\hat\beta_0-\beta_0\|_2$'
    elif mtype == "error.":
        ylab = r'$E_{\beta \sim F}\Big[E_{x,y\sim P_{\beta}}\big[(y - x^\top \hat \beta_0)^2\big]\Big]$'
    else:
        ylab = ""



    for est, col, mk in zip(estimators, cols, mks):
        typ = 'sparse'
        result2 = result1
        msr = mtype + est + '.' + typ
        mean, std = result2[msr]['mean'], result2[msr]['std']
        x = result2['n']
        ax[index].errorbar(x, mean, std, linestyle = '-', color = col, marker = mk, markersize = markersize)
    ax[index].set_xscale('log')
    ax[index].set_yscale('log')
    ax[index].set_xlabel('$n_k$', fontsize = labelsize)
    ax[index].set_ylabel(ylab, labelpad = 0, fontsize = labelsize)
    ax[index].set_xticklabels(x, fontsize = fontsize)
    ax[index].set_xticks(x)
    ax[index].tick_params(axis = 'y', labelsize = fontsize)



    result1 = result
    result1 = result1.loc[result['m'] != 5]
    result1 = result1.loc[result['s'] == 5]
    index = 1

    types = ['sparse']
    mks = ['x']

    estimators = ['adele', 'mrlasso']
    cols = ['orange', 'blue']
    mks = ['x', 'o']


    for est, col, mk in zip(estimators, cols, mks):
        typ = 'sparse'
        result2 = result1
        msr = mtype + est + '.' + typ

        mean, std = result2[msr]['mean'], result2[msr]['std']
        x = result2['m']
        ax[index].errorbar(x, mean, std, linestyle = '-', color = col, marker = mk, markersize = markersize)
    ax[index].set_xscale('log')
    ax[index].set_yscale('log')
    ax[index].set_xlabel('$m$', fontsize = labelsize)
    ax[index].set_xticks(x)
    ax[index].set_xticklabels(x, fontsize = fontsize)
    ax[index].tick_params(axis = 'y', labelsize = fontsize)


    result1 = result
    result1 = result1.loc[result['m'] == 5]
    result1 = result1.loc[result['s'] != 5]
    index = 2

    types = ['sparse']
    mks = ['x']

    estimators = ['adele', 'mrlasso']
    cols = ['orange', 'blue']
    mks = ['x', 'o']


    for est, col, mk in zip(estimators, cols, mks):
        typ = 'sparse'
        result2 = result1
        msr = mtype + est + '.' + typ

        mean, std = result2[msr]['mean'], result2[msr]['std']
        x = result2['s']
        ax[index].errorbar(x, mean, std, linestyle = '-', color = col, marker = mk, markersize = markersize)
    ax[index].set_xscale('log')
    ax[index].set_yscale('log')
    ax[index].set_xlabel(r'$s\big(\beta_{0, \text{ADELE}}\big)$', fontsize = labelsize)
    ax[index].set_xticks(x)
    ax[index].set_xticklabels(x, fontsize = fontsize)
    ax[index].tick_params(axis = 'y', labelsize = fontsize)




    estimators = ['ADELE', 'MrLasso']
    cols = ['orange', 'blue']
    mks = ['x', 'o']

    labels = []
    lines = []


    for est, col, mk in zip(estimators, cols, mks):
        labels.append(est)

        lines.append(Line2D([0], [0], linestyle = '-', color = col, marker = mk, markersize = markersize))

    ax[0].legend(lines, labels, fontsize = fontsize)
    ax[1].legend(lines, labels, fontsize = fontsize)
    ax[2].legend(lines, labels, fontsize = fontsize)
    
    return fig, ax