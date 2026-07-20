
#!/usr/bin/python3
import csv
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec
import numpy as np

import scipy.stats as stats

import sys
import os
# add the parent directory (where stack_axes.py lives) to the path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from stack_axes import make_stacked_axes


plt.style.use('tableau-colorblind10')

def chi(exp_value, obs_value, dx):
    chi2 = 0
    for exp,obs in zip(exp_value,obs_value):
        chi2 += dx*(obs-exp)**2
    return chi2

def reader(file_name):
    print(file_name)
    x = []
    y = []
    y_both= []
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
xlims_min = -7
xlims_max = 7
#position of # of states abs position
pos_x_name = 0.005
pos_y_name = 0.975

#data
#JUST NEED TO CHANGE THOSE VALUES, works up to 4 graphs
data_perc = ['100','05']
data_r = ['1', '09999', '0999', '099', '09']
data_r_float_w = ['1.0000', '0.9999', '0.999', '0.99', '0.9']
data_r_float = ['0.00000', '0.0001', '0.001', '0.01', '0.1']
#exp_label_10 = [r'$10^{-\infty}$',r'$10^{-3}$',r'$10^{-2}$',r'$10^{-1}$']
data_site = ['16','16']

ylims = [0.44,0.44]



#lines esthetics
line_style = ['solid']*6
line_width = [1.5]*6
markers = ['.','x']
m_size = [70,40]

#Colors #states
cmap = plt.get_cmap('viridis')
colorsa = np.linspace(0,1,num=len(data_r)-1);
colors = [cmap(i) for i in colorsa]
colors.insert(0,'k')
colors[-1]='goldenrod'

#legend labels
labels = []
letters = ['a) ','b) ','c) ','d) ']
for j in data_perc:
    labels.append('{:.2f}'.format(float(j)/100))
labels.append(r'$\chi^2(w_t)$')
handles = []

#Set font family
hfont = {'fontname':'Times'}
plt.rc('font', family='serif',size=14)
plt.rc('text', usetex=True)
plt.rc('legend', fontsize=14)
size_text = 16

# Sample data
q_matrix_files =[]
for i in data_perc:
    q_matrix_files.append([])
    for j in data_r:
        q_matrix_files[-1].append('./'+'16s_u8/dos_green_'+i+'_'+j+'.txt')

#print(q_matrix_files)
# Create subplots with shared x-axis
fig,axes = make_stacked_axes(n_panels=2,panel_height_cm=3.2,left=0.01,right=0.99,top_margin_cm=0.1)
gs = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[1, 1])

# Create inset axis (for c)) positioned relative to the lower plot
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
rect = 0.42,0.40,0.16,0.30
inset_ax = fig.add_axes(rect)
axes.append(inset_ax)  # c), 2, figsize=(6, 3),constrained_layout=True)
#Shared x-axes
axes[1].sharex(axes[0])
axes[0].label_outer()   # hides x labels on top row

#BOX PARAMS
ss_x = 0.93#xlims_max*0.50
ss_y = 1.0#*np.array(ylims)
ss_spacing = 0.15 #Was .15
ss_size_text = 12 #Was 16
ha_pos_title = 'center'
ha_pos = 'right'
va_pos = 'top'

box_h = 0.22
box_pad = 0.01

num_slider = [-0.0075,-0.005,-0.005,-0.005,0]

A_omega = [r"$A(\omega)$",r"$A(\omega)_{reduced}$"]

#SPECTRAL FUNCTIONS
for i,ax in enumerate(axes):

    if i < 2 :
        q_files = q_matrix_files[i]
        #print(q_files)
        x_ref,y_ref,_ = reader(q_files[0])

        for j in range(len(q_matrix_files[i])):
            x,y,_ = reader(q_files[j])
            #SPECTRAL FUNCTIONS
            handle, = ax.plot(x, np.array(y),linestyle=line_style[j],lw=line_width[j],color=colors[j],label=data_r_float_w[j])

            ax.set_ylim(0, float(ylims[i]))
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ax.set_xlim(xlims_min,xlims_max)
            ax.tick_params(axis='x',labelsize=size_text+1)
            ax.set_yticks([0, ylims[i]])  # Tick only at 0 and 1 and 2

            #Vert line at 0
            ax.axvline(c='k',ls='--',lw=0.75)

            ax.scatter(xlims_max*0.95,ylims[i]*0.9,marker=markers[i],color='dimgray',s=50)

            #Remove y tick labels
            ax.set_yticklabels([])
            if i == 0 :#and j!=0:  # Add legend only to the first plot
                handles.append(handle)
                #ax.legend(loc='upper right')
            #CHI SQUARE CONVERGENCE
            if j > 0 :
                ks = stats.kstest(y,y_ref)
                chi_2 = chi(y_ref,y,x[1]-x[0])
                print(chi_2)


                axes[2].scatter(float(data_r_float[j]),chi_2,color=colors[j], marker=markers[i],s=m_size[i])

    # Optionally add individual titles a), b), ...
    if i == 2 : shift = 0.025
    else :shift = 0
    f = ["$f$=","$f$=",""]
    ax.text(pos_x_name+shift,pos_y_name,letters[i] + f[i]+ labels[i],fontsize=size_text,ha='left',va='top', transform=ax.transAxes)

#Reverse handles
h0, = axes[0].plot(0, 1,label='$w_t$',alpha=0)
reverse_handle = []
for i in handles:
    reverse_handle.insert(0,i)
reverse_handle.insert(0,h0)
handles = reverse_handle

# Set x-axis label on the bottom subplot only
leg = fig.legend(handles=handles,loc='upper center',
        bbox_to_anchor=(0.525,1.01),
        #bbox_to_anchor=(0.725,0.16),
        #bbox_to_anchor=(0.66,0.5828),
        fontsize=size_text,ncol=6,
        columnspacing=0.5,labelspacing=0.1,
        handletextpad=0.2,handlelength=0.65,borderpad=0.1,
        frameon=True, markerfirst=True)
axes[1].set_xlabel('$\omega$',fontsize = size_text)

#axes[2].set_xlabel('1-$w_t$',fontsize = size_text-4,labelpad=0)
#axes[2].set_ylabel('$\chi^2$',fontsize = size_text-4,rotation=0)
axes[2].set_facecolor((1, 1, 1, 0.8))  # RGBA, alpha=0.5
axes[2].set_xscale('log')
axes[2].tick_params(
    axis='both',        # both x and y
    which='both',       # major and minor
    bottom=False,       # remove bottom ticks
    top=False,          # remove top ticks
    left=False,         # remove left ticks
    right=False,        # remove right ticks
    labelbottom=False,  # remove x labels
    labelleft=False     # remove y labels
)
axes[2].set_xlim(2e-1,5e-5)
axes[2].set_ylim(-0.005,0.15)
axes[2].axhline(0,c='k',lw=0.7,linestyle='--')

plt.savefig('further_truncation.pdf')

plt.show()
