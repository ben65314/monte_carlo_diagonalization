import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import FancyBboxPatch
import matplotlib
import numpy as np
import matplotlib as mpl
import csv

import sys
import os
# add the parent directory (where stack_axes.py lives) to the path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from stack_axes import make_stacked_axes

def ax_to_fig(ax, x, y):
    fig = ax.figure
    display_coords = ax.transAxes.transform((x, y))
    f = fig.transFigure.inverted().transform(display_coords)
    return f[0], f[1]

def reader(file_name):
    """
    READER

        """
    x = []
    y = []
    y_both = []
    with open(file_name) as file:
        reader = csv.reader(file,delimiter='\t')
        for lines in reader:
            if(lines[0].startswith('#') or lines[0].startswith('PARAMETERS') or lines[0].startswith('NEXT_SITE')):
                continue
            x.append(float(lines[0]))
            #y.append(float(lines[1])+float(lines[2]))
            y.append(float(lines[1])+float(lines[2]))
            y_both.append([float(lines[1]),float(lines[2])])
    return x,y,y_both

##Curves information
#limits
xlims_min = 0#-7.5
xlims_max = 7.5
#position of # of states abs position
pos_x_name = 0.005
pos_y_name = 0.975

#data
#JUST NEED TO CHANGE THOSE VALUES, works up to 4 graphs
data_perc = '05'
data_site = 16
data_methods = ['MCD', 'BFS-Boltzmann', 'DFS', 'Monte-Carlo']
data_methods_label = ['MCD', 'Variant 1', 'Variant 2', 'Variant 3']
data_beta = [[0.1, 0.2, 0.3], # MCD
             [0.1, 0.2, 0.3], # BFS-Boltzmann
             [0.1, 0.2, 0.3], # DFS
             [0.1, 0.2, 0.3]] # Monte-Carlo
data_mu_label = ['2','1.6','1.2','0.8','0.4']
ylims = [0.45,0.45,0.45,0.45,0.6,0.6]


#lines esthetics
line_style = ['solid']*5
line_width = [1.5]*5
lwidth = 1.5 #if j!=0 else 5

#Colors #states
cmap = plt.get_cmap('gnuplot')
colorsa = np.linspace(0,.90,num=4);
colors = [cmap(i) for i in colorsa]

betas = ['ED', 0.1, 0.2, 0.3]
betas_label = ['ED   ', r'$\beta:0.1$', 0.2, 0.3]
color_dict = {}
for i,b in enumerate(betas):
    color_dict[b] = colors[i]


#legend labels
labels = ['$N_{1.00}$','$N_{0.05}$','Ncst']
letters = ['a) ','b) ','c) ','d) ','e) ','f) ']
handles = []

#Set font family
hfont = {'fontname':'Times'}
plt.rc('font', family='serif',size=14)
plt.rc('text', usetex=True)
plt.rc('legend', fontsize=14)
size_text = 16


# Sample data
q_matrix_files =[]
for i,method in enumerate(data_methods):
    q_matrix_files.append([])
    for j in data_beta[i]:
        q_matrix_files[-1].append('./'+method+'/dos_green05_'+str(j)+'.txt')

# Create subplots with shared x-axis
fig,axes = make_stacked_axes(n_panels=4,panel_height_cm=3.2,left=0.01,right=0.99,top_margin_cm=0.1)

#BOX PARAMS
ss_x = 0.9975#xlims_max*0.50
ss_y = 0.98#*np.array(ylims)
ss_spacing = 0.14 #Was .15
ss_size_text = 14 #Was 16
ha_pos_title = 'center'
ha_pos = 'right'
va_pos = 'top'

box_h = 0.22
box_pad = 0.01

num_slider = [-0.0075,-0.005,-0.005,-0.005,0]

nf = [0.965,0.93,0.90,0.85]
handle1 = handle2 = None
for i,ax in enumerate(axes):
    q_files = q_matrix_files[i]
    print(q_files)

    # Plot references
    x,y,y2 = reader("./dos_green100.txt")
    handle, = ax.plot(x, np.array(y),linestyle=line_style[0],lw=2,color='k')

    for j in range(len(q_matrix_files[i])):

        x,y,y_both = reader(q_files[j])
        print(q_files[j][18:-4])
        handle, = ax.plot(x, np.array(y),linestyle=line_style[j],lw=lwidth,color=color_dict[data_beta[i][j]],label=labels[j],alpha=0.65)
        ax.set_ylim(0, float(ylims[i]))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax.set_xlim(xlims_min,xlims_max)
        ax.tick_params(axis='x',labelsize=size_text+1)
        ax.set_yticks([0, ylims[i]])  # Tick only at 0 and 1 and 2

        #Vert line at 0
        ax.axvline(c='k',ls='--',lw=0.75)

        #Remove y tick labels
        ax.set_yticklabels([])
        if i == 0:  # Add legend only to the first plot
            handles.append(handle)
            #ax.legend(loc='upper right')

    # Optionally add individual titles or labels
    ax.text(pos_x_name,pos_y_name,letters[i],fontsize=size_text,ha='left',va='top', transform=ax.transAxes)
    ax.text(ss_x,ss_y - 0*ss_spacing, data_methods_label[i], fontsize=size_text,ha='right',va='top', transform=ax.transAxes,c='k')

# legend
handles_color = []
for i,b in enumerate(betas):
    hdl, = axes[0].plot(20,20, color = color_dict[b], label=betas_label[i],alpha=0.65)
    handles_color.append(hdl)
data_methods = ['MCD', 'BFS-Boltzmann', 'DFS', 'Monte-Carlo']

# Set x-axis label on the bottom subplot
legend = fig.legend(handles=handles_color,loc='center',bbox_to_anchor=(0.5,0.96),fontsize=size_text,ncol=10,columnspacing=0.7,labelspacing=0.1,handletextpad=0.1,handlelength=1,borderpad=0.1, frameon=True,markerfirst=False)
axes[-1].set_xlabel('$\omega$',fontsize = size_text)

# Shared legend placed outside the plots (bottom center)
plt.savefig('method_variation.pdf',bbox_inches='tight', bbox_extra_artists=[legend])

plt.show()

