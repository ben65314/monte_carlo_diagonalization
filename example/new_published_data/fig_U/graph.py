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

def reader(file_name):
    x = []
    y = []
    y_both= []
    with open(file_name) as file:
        reader = csv.reader(file,delimiter='\t')
        for lines in reader:
            if(lines[0].startswith('#') or lines[0].startswith('PARAMETERS') or lines[0].startswith('NEXT_SITE')):
                continue
            x.append(float(lines[0]))
            y.append(float(lines[1])+float(lines[2]))
            y_both.append([float(lines[1]),float(lines[2])])
    return x,y,y_both

##Curves information
#limits
xlims_min = -9
xlims_max = 9
#position of # of states abs position
pos_x_name = 0.005
pos_y_name = 0.975

#data
#JUST NEED TO CHANGE THOSE VALUES, works up to 4 graphs
data_perc = ['100','50','25','15','05']
data_u = ['4','8','12']
ylims = [0.7,0.4,0.4,0.4,0.4]
#legend labels
labels = []
letters = ['a) ','b) ','c) ','d) ','e) ']
for j in data_perc:
    labels.append('{:.2f}'.format(float(j)/100))
labels[0] = '$f$: {:.2f}'.format(float(data_perc[0])/100)
handles = []

#lines esthetics
line_style = ['solid']*5
line_width = [1.5]*5

#Colors #states
cmap = plt.get_cmap('gnuplot')
colorsa = np.linspace(0,0.90,num=len(data_perc));
colors = [cmap(i) for i in colorsa]



#Set font family
hfont = {'fontname':'Times'}
plt.rc('font', family='serif',size=14)
plt.rc('text', usetex=True)
plt.rc('legend', fontsize=14)
size_text = 16


# Sample data
q_matrix_files =[]
for j in data_u:
    q_matrix_files.append([])
    for i in data_perc:
        q_matrix_files[-1].append('./u'+j+'/dos_green'+i+'.txt')

# Create subplots with shared x-axis
fig,axes = make_stacked_axes(n_panels=3,panel_height_cm=3.2,left=0.01,right=0.99,top_margin_cm=0.1)


for i,ax in enumerate(axes):
    q_files = q_matrix_files[i]
    print(q_files)
    for j in range(len(q_matrix_files[i])):

        x,y,_ = reader(q_files[j])
        handle, = ax.plot(x, np.array(y),linestyle=line_style[j],lw=line_width[j],color=colors[j],label=labels[j])
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

    # Optionally add individual titles or labels
    ax.text(pos_x_name,pos_y_name,letters[i] + r'$U=$ ' + data_u[i],fontsize=size_text,ha='left',va='top', transform=ax.transAxes)

# Set x-axis label on the bottom subplot only
legend = fig.legend(handles=handles,loc='center',bbox_to_anchor=(0.5,0.96),fontsize=size_text,ncol=5,
        columnspacing=.5,labelspacing=0.1,handletextpad=0.1,handlelength=1,borderpad=0.1, frameon=True,markerfirst=False)
axes[-1].set_xlabel('$\omega$',fontsize = size_text)

# Shared legend placed outside the plots (bottom center)
plt.savefig('var_U.pdf', bbox_inches='tight', bbox_extra_artists=[legend])

plt.show()

