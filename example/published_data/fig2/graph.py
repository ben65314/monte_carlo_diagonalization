#!/usr/bin/python3
import csv
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import sys;
from matplotlib.ticker import AutoMinorLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
import matplotlib as mpl

from matplotlib.legend_handler import HandlerBase
from matplotlib.colors import LinearSegmentedColormap

# --- Helper: Truncate colormap and make it discrete
def truncate_colormap_discrete(cmap, minval=0.0, maxval=1.0, n_colors=5):
    colors = cmap((np.linspace(minval, maxval, n_colors))**(0.2))
    return LinearSegmentedColormap.from_list(f'{cmap.name}_trunc', colors, N=n_colors)

# --- Custom Handler for discrete colormap patches
class HandlerDiscreteColormap(HandlerBase):
    def __init__(self, n_colors=5, **kwargs):
        self.n_colors = n_colors
        super().__init__(**kwargs)

    def create_artists(self, legend, orig_handle, xdescent, ydescent,
                       width, height, fontsize, trans):
        cmap = orig_handle.get_cmap()
        blocks = []
        block_width = width / self.n_colors
        for i in range(self.n_colors):
            color = cmap(i / (self.n_colors - 1))
            rect = mpl.patches.Rectangle(
                (xdescent + i * block_width - 0.5, ydescent),
                block_width + 1, height,
                transform=trans,
                facecolor=color,
                edgecolor='none', linewidth=0
            )
            blocks.append(rect)
        return blocks

if (len(sys.argv) == 1): exit
def reduce_data(data,start=0,keep=1):
    n = len(data[0]) #Initial 
    k = int((1-keep) * n) #Reduce
    indices_to_remove = np.linspace(start, n - 1, k, dtype=int)
    
    # Remove the selected indices
    mask = np.ones(n, dtype=bool)
    mask[indices_to_remove] = False
    
    return_data = []
    for d in data:
        print('d:',len(d))
        return_data.append(np.array(d)[mask])
    return return_data


# Set font family and text bubbles
hfont = {'fontname':'Times'}
plt.rcParams['text.usetex'] = True
plt.rc('font', family='serif',size=14)
plt.rc('text', usetex=True)
plt.rc('legend', fontsize=24)
size_text = 14

# Color maps for each
cmap_red   = plt.get_cmap("Reds")
cmap_yellow= plt.get_cmap("YlOrBr")
cmap_green = plt.get_cmap("Greens")
cmap_blue  = plt.get_cmap("Blues") 
cmaps = [cmap_red,cmap_yellow,cmap_green,cmap_blue]

# Markers
markers = ["o","s","D","^","v","*","P","X"]
zorders = [10+i for i in range(len(markers))]
start_marker = 8
diff_size_marker = 10
markers_size = [start_marker + diff_size_marker*i for i in range(len(markers))]
markers_size = [12,32,24,32,32,48,32,32]

sites = [16,12,8,4]

# Sample data
q_matrix_files = []
for j in sites:
    q_matrix_files.append('./'+str(j)+'s/new_fund.txt')

q_ghost_files = []
for j in sites:
    q_ghost_files.append('./'+str(j)+'s/full_fund.txt')

# Figure
fig, ax = plt.subplots(figsize=(6,4),constrained_layout=True)
symbol_legend = []
color_legend = []

# Ghost outline (gray)
for i in range(len(q_ghost_files)):

    cumul_weight = []
    cumul_nU = []
    n_states = []
    data_file_name = q_ghost_files[i]

    # READING
    with open(data_file_name) as file:
        reader = csv.reader(file)
        for lines in reader:
            # Ignore comments
            if(lines[0].startswith('#')):
                continue;
            line = lines[0].split()
            cumul_weight.append(float(line[1]))
            cumul_nU.append(int(line[2]))
            n_states.append(float(line[0]))

    combined = sorted(zip(n_states,cumul_weight,cumul_nU))
    n_states, cumul_weight, cumul_nU = zip(*combined)
    n_states = list(n_states)
    cumul_weight = list(cumul_weight)
    cumul_nU = list(cumul_nU)

    #Ghost outline (gray)
    ax.plot(n_states,cumul_weight,color='k',alpha=0.2,lw=9,zorder=-100)
    
# Data points
for i in range(len(q_matrix_files)):
    #Color map used
    max_color, min_color = 0.95, 0.006
    if i == 1 :
        max_color, min_color = 0.5, 0.012
    colorsa = (np.linspace(min_color,max_color,num=(int(sites[i]/2)+1)))**(0.2);
    colors = [cmaps[i](j) for j in colorsa]

    cumul_weight = []
    cumul_nU = []
    n_states = []
    data_file_name = q_matrix_files[i]
    
    # READING
    with open(data_file_name) as file:
        reader = csv.reader(file)
        for lines in reader:
            # Ignore comments
            if(lines[0].startswith('#')):
                continue;
            line = lines[0].split()
            cumul_weight.append(float(line[1]))
            cumul_nU.append(int(line[2]))
            n_states.append(float(line[0]))

    combined = sorted(zip(n_states,cumul_weight,cumul_nU))
    n_states, cumul_weight, cumul_nU = zip(*combined)
    n_states = list(n_states)
    cumul_weight = list(cumul_weight)
    cumul_nU = list(cumul_nU)
    
    # Seperate data according to its nU value
    cw = [[] for l in range(max(cumul_nU)+1)]#[[]]*(max(cumul_nU)+1)
    cs = [[] for l in range(max(cumul_nU)+1)]#[[]]*(max(cumul_nU)+1)
    for j, nu in enumerate(cumul_nU):
        cw[nu].append(cumul_weight[j])
        cs[nu].append(n_states[j])
    
    # Scatter plot
    for m, (c_weight, c_states) in enumerate(zip(cw,cs)):
        using_c = colors[m]
        if m!=-1:
          col = 'black'
        else:
          col = colors[m]
        sl = plt.scatter(c_states,c_weight,facecolor=using_c,marker=markers[m],s=markers_size[m],zorder=zorders[m]+10*i, edgecolors=col,lw=0.02)

    # Notify nU weight markers (Used for the nU legend)
    if i == 0:
        for j,k in enumerate(markers):
            col = 'black'
            sl = plt.scatter([5,5],[1,1],marker=k,label=j,s=markers_size[j],facecolor='white',edgecolors=col,lw=0.05)
            symbol_legend.append(sl)

    x_c = np.linspace(10000,len(colors)+100001,len(colors))
    y_c = [0.5]*len(colors)
    cm2 = truncate_colormap_discrete(cmaps[i],min_color,max_color,len(colors))
    ax.scatter(x_c,y_c,c=colors,cmap=cm2)

    # Dummy handle with colormap info
    gradient_patch = mpl.cm.ScalarMappable(cmap=cm2)
    color_legend.append(gradient_patch)




ax.set_xscale('log')
ax.set_xlim(3e-9, 1.7)  # Can't use 0 on log scale
ax.set_ylim(0, 1.05)
ax.tick_params(which='minor', length=4, color='black')
ax.grid(True)
#plt.xscale('log')

# Affichage
# Add to legend
sites_label = [str(i) for i in sites]
legend_sites = ax.legend(color_legend, sites_label, loc='lower right',
            handler_map={type(color_legend[0]): HandlerDiscreteColormap(int(sites[3]/2+1)),
            type(color_legend[1]): HandlerDiscreteColormap(int(sites[2]/2+1)),
            type(color_legend[2]): HandlerDiscreteColormap(int(sites[1]/2+1)),
            type(color_legend[3]): HandlerDiscreteColormap(int(sites[0]/2+1))},
            bbox_to_anchor=(1,0.),
            bbox_transform=ax.transAxes,
            fontsize=size_text,
            ncol=1,columnspacing=0.2,
            labelspacing=0.1,
            handletextpad=0.2,
            handlelength=1.,
            borderpad=0.2, frameon=True, 
            title=r'$N_c$')
legend_sites._legend_box.align = "center"

ax.add_artist(legend_sites)
legend_markers = ax.legend(handles=symbol_legend, loc='upper left',
        bbox_to_anchor=(0.,1),
        bbox_transform=ax.transAxes,
        fontsize=size_text,
        ncol=1,columnspacing=0.2,
        labelspacing=0.1,
        handletextpad=0.2,
        handlelength=1,
        borderpad=0.2, frameon=True,
        title=r'$n_U$')
ax.add_artist(legend_markers)
legend_markers._legend_box.align = "center"

plt.xticks([1e-8, 1e-6, 1e-4, 1e-2, 1])
ax.set_xlabel(r'$f$')
ax.set_ylabel(r'$w_f$')
plt.savefig("energy_conv_2D.pdf")

