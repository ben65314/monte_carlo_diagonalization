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
xlims_min = -6
xlims_max = 10
#position of # of states abs position
pos_x_name = 0.005
pos_y_name = 0.975

#data
data_perc = ['100','Ncst']
data_site = 16
data_mu = ['2','1_6','1_2','0_8','0_4','0']
data_mu_label = ['2.0','1.6','1.2','0.8','0.4', '0.0']
ylims = [0.35,0.35,0.45,0.55,0.8,0.6]
#legend labels
labels = ['ED',r'$N_{f}$']
letters = ['a) ','b) ','c) ','d) ','e) ','f) ']
handles = []


ref_sector = [r'$N_e=16$'+'\t'+'$S_z=0$', r'$N_e=16$'+'\t'+'$S_z=0$', r'$N_e=14$'+'\t'+'$S_z=0$', r'$N_e=12$'+'\t'+'$S_z=0$', r'$N_e=10$'+'\t'+'$S_z=0$', r'$N_e=10$'+'\t'+'$S_z=0$']

ref_number = [r'$N_\mathrm{tot}$=165,636,900', r'$N_\mathrm{tot}$=165,636,900', r'$N_\mathrm{tot}$=130,873,600', r'$N_\mathrm{tot}$=64,128,064', r'$N_\mathrm{tot}$=19,079,424', r'$N_\mathrm{tot}$=19,079,424']
per_Ncst = [r'$f=0.05$',r'$f=0.05$',r'$f=0.06$',r'$f=0.13$',r'$f=0.43$',r'$f=0.43$']
Ncst = "8,281,845"

#lines esthetics
line_style = ['solid']*5
line_width = [1.5]*5

#Colors #states
cmap = plt.get_cmap('gnuplot')
colorsa = np.linspace(0,0.85,num=len(data_perc));
colors = [cmap(i) for i in colorsa]

#Set font family
hfont = {'fontname':'Times'}
plt.rc('font', family='serif',size=14)
plt.rc('text', usetex=True)
plt.rc('legend', fontsize=14)
size_text = 16

# Sample data
q_matrix_files =[]
for j in data_mu:
    q_matrix_files.append([])
    for i in data_perc:
        if (i == '05') :
            continue
        q_matrix_files[-1].append('./mu'+j+'/dos_green'+i+'.txt')

# Create subplots with shared x-axis
fig,axes = make_stacked_axes(n_panels=6,panel_height_cm=3.2,left=0.01,right=0.99,top_margin_cm=0.1)

#BOX PARAMS
ss_x = 0.9975#xlims_max*0.50
ss_y = 0.975#*np.array(ylims)
ss_spacing = 0.15 #Was .15
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
    for j in range(len(q_matrix_files[i])):
        alpha = 1 if j!=0 else 0.75
        lwidth = 2 #if j!=0 else 5

        x,y,y_both = reader(q_files[j])
        handle, = ax.plot(x, np.array(y)+0.01,linestyle=line_style[j],lw=lwidth,color=colors[j],label=labels[j],alpha=alpha)
        if j != 0:
            ##if i != 2 : 
            handle1, = ax.plot(x, np.array(y_both)[:,0],linestyle='--',lw=1.,label=r"$G^-(\omega)$",color='red', alpha=.75)
            handle2, = ax.plot(x, np.array(y_both)[:,1],linestyle='--',lw=1.,label=r"$G^+(\omega)$",color='blue', alpha=.75)

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
            #pass
            handles.append(handle)
            #ax.legend(loc='upper right')
    # Optionally add individual titles or labels
    ax.text(pos_x_name,pos_y_name,letters[i] + r'$\mu=$ ' + data_mu_label[i],fontsize=size_text,ha='left',va='top', transform=ax.transAxes)

    #Sectors
    ax.text(ss_x,ss_y - 0*ss_spacing, ref_sector[i], fontsize=ss_size_text, ha='right', va='top', transform=ax.transAxes,c=colors[0])
    #Tot number
    ax.text(ss_x,ss_y - 1*ss_spacing, ref_number[i], fontsize=ss_size_text, ha='right', va='top', transform=ax.transAxes,c=colors[0])
    # f
    ax.text(ss_x,ss_y - 2*ss_spacing, per_Ncst[i], fontsize=ss_size_text, ha='right', va='top', transform=ax.transAxes,c=colors[1])

handles.append(handle1)
handles.append(handle2)

# Set x-axis label on the bottom subplot only
legend = fig.legend(handles=handles,loc='center',bbox_to_anchor=(0.5,0.93),fontsize=size_text,ncol=1,columnspacing=0.5,labelspacing=0.1,handletextpad=0.3,handlelength=1,borderpad=0.2, frameon=True,markerfirst=True)
axes[-1].set_xlabel('$\omega$',fontsize = size_text)

# Shared legend placed outside the plots (bottom center)
plt.savefig('mu_variation.pdf', bbox_inches='tight', bbox_extra_artists=[legend])

plt.show()

