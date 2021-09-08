import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

outlier_dist = 10

def get_plot(result, mtype):
    
    
    alpha_shade = 0.1
    col_shade = 'g'
    
    global outlier_dist
    mpl.rcParams['text.usetex'] = True
    mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

    fontsize = 15
    labelsize = 18
    markersize = 8

    fig, ax = plt.subplots(2, 2, figsize = (8, 6), constrained_layout = True)




    pars = {'m': 5}
    index = 0,0
    result1 = result
    result1 = result1.loc[result1['eta.pre'] < 0]
    result1 = result1.loc[result1['m'] == 5]
    result1 = result1.loc[result1['s'] == 5]
    result1 = result1.loc[result1['outlier.dist'] == outlier_dist]

    types = ['sparse']
    mks = ['x']

    estimators = ['adele', 'mrlasso']
    cols = ['orange', 'blue']
    mks = ['x', 'o']

    if mtype == "l2.":
        ylab = r'$\|\hat\beta_0-\beta_0\|_2$'
    elif mtype == "error.":
        ylab = r'$E_{\beta}\Big[E_{(x,y)}\big[(y - x^\top \hat \beta_0)^2\big]\Big]$'
    else:
        ylab = ""



    for est, col, mk in zip(estimators, cols, mks):
        typ = 'sparse'
        result2 = result1
        msr = mtype + est + '.' + typ
        mean, std = result2[msr]['mean'], result2[msr]['std']
        x = result2['n']
        ax[index].errorbar(x, mean, std, linestyle = '-', color = col, marker = mk, markersize = markersize)
        
        
        if mtype == "l2.":
            msr = 'eta'
            mean, std = result2[msr]['mean'], result2[msr]['std']
            ax[index].errorbar(x, mean, std, linestyle = ':', color = col, marker = mk, markersize = markersize)
            ax[index].fill_between(x, 6, 1, alpha = alpha_shade, color = col_shade)
            
            
    ax[index].set_xscale('log')
    ax[index].set_yscale('log')
    ax[index].set_xlabel('$n_k$', fontsize = labelsize)
    ax[index].set_ylabel(ylab, labelpad = 0, fontsize = labelsize)
    ax[index].set_xticklabels(x, fontsize = fontsize)
    ax[index].set_xticks(x)
    ax[index].tick_params(axis = 'y', labelsize = fontsize)



    result1 = result
    
    result1 = result1.loc[result1['eta.pre'] < 0]
    result1 = result1.loc[result['m'] != 5]
    index = 0, 1

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
        
        if mtype == "l2.":
            msr = 'eta'
            mean, std = result2[msr]['mean'], result2[msr]['std']
            ax[index].errorbar(x, mean, std, linestyle = ':', color = col, marker = mk, markersize = markersize)
            ax[index].fill_between(x, 6, 1, alpha = alpha_shade, color = col_shade)
            
    ax[index].set_xscale('log')
    ax[index].set_yscale('log')
    ax[index].set_xlabel('$m$', fontsize = labelsize)
    ax[index].set_xticks(x)
    ax[index].set_xticklabels(x, fontsize = fontsize)
    ax[index].tick_params(axis = 'y', labelsize = fontsize)


    result1 = result
    result1 = result1.loc[result1['eta.pre'] < 0]
    result1 = result1.loc[result1['s'] != 5]
    result1 = result1.loc[result1['outlier.dist'] == outlier_dist]
    index = 1, 0

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
        
        if mtype == "l2.":
            msr = 'eta'
            mean, std = result2[msr]['mean'], result2[msr]['std']
            ax[index].errorbar(x, mean, std, linestyle = ':', color = col, marker = mk, markersize = markersize)
            ax[index].fill_between(x, 6, 1, alpha = alpha_shade, color = col_shade)
            
    ax[index].set_xscale('log')
    ax[index].set_yscale('log')
    ax[index].set_xlabel(r'$s$', fontsize = labelsize)
    ax[index].set_xticks(x)
    ax[index].set_ylabel(ylab, labelpad = 0, fontsize = labelsize)
    ax[index].set_xticklabels(x, fontsize = fontsize)
    ax[index].tick_params(axis = 'y', labelsize = fontsize)
    

    result1 = result
    result1 = result1.loc[result1['eta.pre'] > 0]
    index = 1, 1

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
        x = result2['eta.pre']
        ax[index].errorbar(x, mean, std, linestyle = '-', color = col, marker = mk, markersize = markersize)
   
    
    ax[index].set_xlabel(r'$\eta$', fontsize = labelsize)
    ax[index].set_xticks(x)
    ax[index].set_xticklabels(x, fontsize = fontsize)
    ax[index].tick_params(axis = 'y', labelsize = fontsize)




    estimators = ['SHIR', 'MrLasso']
    cols = ['orange', 'blue']
    mks = ['x', 'o']

    

    for i in [0, 1]:
        for j in [0, 1]:
            labels = []
            lines = []


            for est, col, mk in zip(estimators, cols, mks):
                labels.append(est)
                lines.append(Line2D([0], [0], linestyle = '-', color = col, marker = mk, markersize = markersize))
            if i*j == 0 and mtype == "l2.":
                lines.append(Line2D([0], [0], linestyle = ':', color = 'b', marker = 'o', markersize = markersize))
                labels.append(r'$\hat\eta$')
            ax[i, j].legend(lines, labels, fontsize = fontsize)
            
    
    
    return fig, ax










def get_plot2(result, mtype):
    
    global outlier_dist
    mpl.rcParams['text.usetex'] = True
    mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

    fontsize = 15
    labelsize = 18
    markersize = 8

    fig, ax = plt.subplots(1, 3, figsize = (14, 3.5), constrained_layout = True)




    pars = {'m': 5}
    index = 0
    result1 = result
    result1 = result1.loc[result1['eta.pre'] < 0]
    result1 = result1.loc[result1['m'] == 5]
    result1 = result1.loc[result1['s'] == 5]
    result1 = result1.loc[result1['outlier.dist'] == outlier_dist]

    types = ['sparse']
    mks = ['x']

    estimators = ['adele', 'mrlasso']
    cols = ['orange', 'blue']
    mks = ['x', 'o']

    if mtype == "l2.":
        ylab = r'$\|\hat\beta_0-\beta_0\|_2$'
    elif mtype == "error.":
        ylab = r'$E_{\beta}\Big[E_{(x,y)}\big[(y - x^\top \hat \beta_0)^2\big]\Big]$'
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
    result1 = result1.loc[result1['eta.pre'] < 0]
    result1 = result1.loc[result['m'] != 5]
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
    result1 = result1.loc[result1['eta.pre'] < 0]
    result1 = result1.loc[result1['m'] == 5]
    result1 = result1.loc[result1['s'] != 5]
    result1 = result1.loc[result1['outlier.dist'] == outlier_dist]
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
    ax[index].set_xlabel(r'$s$', fontsize = labelsize)
    ax[index].set_xticks(x)
    ax[index].set_xticklabels(x, fontsize = fontsize)
    ax[index].tick_params(axis = 'y', labelsize = fontsize)
    

   



    estimators = ['SHIR', 'MrLasso']
    cols = ['orange', 'blue']
    mks = ['x', 'o']

    labels = []
    lines = []


    for est, col, mk in zip(estimators, cols, mks):
        labels.append(est)

        lines.append(Line2D([0], [0], linestyle = '-', color = col, marker = mk, markersize = markersize))
    for i in [0, 1, 2]:
        ax[i].legend(lines, labels, fontsize = fontsize)
    
    
    return fig, ax