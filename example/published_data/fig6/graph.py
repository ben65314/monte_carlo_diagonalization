
#!/usr/bin/python3
import csv
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import FancyBboxPatch
from matplotlib.ticker import AutoMinorLocator
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
import numpy as np

import scipy.stats as stats

from matplotlib.legend_handler import HandlerBase
from matplotlib.colors import LinearSegmentedColormap


def chi(exp_value, obs_value, dx):
    chi2 = 0
    for exp,obs in zip(exp_value,obs_value):
        chi2 += dx*(obs-exp)**2
    return chi2

def cum_curve(exp_values,obs_values,x=None,name = ''):
    cum_sum_exp = np.cumsum(exp_values)
    cum_sum_obs = np.cumsum(obs_values)

    #if(x!=None):

    #    fig,ax = plt.subplots(1, 1,figsize=(6,3),constrained_layout=True)
    #    ax.plot(x,cum_sum_exp)
    #    ax.plot(x,cum_sum_obs)

    #    fig.savefig('csum'+name+'.png')
    chi2 = chi(cum_sum_exp,cum_sum_obs)
    return chi2


def norm(data):
    return (data)/(max(data)-min(data))

def reader(file_name):
    #print(file_name)
    x = []
    y = []
    #y = []
    with open(file_name) as file:
        reader = csv.reader(file,delimiter='\t')
        for lines in reader:
            if(lines[0].startswith('#') or lines[0].startswith('PARAMETERS') or lines[0].startswith('NEXT_SITE')):
                continue
            x.append(float(lines[0]))
            #y.append(float(lines[1])+float(lines[2]))
            y.append(float(lines[1])+float(lines[2]))
            #print
    return x,y

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

ylims = [0.35,0.35,0.45,0.45,0.45,0.45]



#lines esthetics
line_style = ['solid']*6
#line_style = ['solid','dashed','dashdot',(0,(1,1)),'solid']
#line_width = [2.6,2.2,1.8,1.4,1]
line_width = [1.5]*6
markers = ['.','x']
m_size = [70,40]
#line_width = [2,2,2,2,1]
#Colors #states
cmap = plt.get_cmap('viridis')
colorsa = np.linspace(0.2,0.9,num=len(data_r)-1);
colors = [cmap(i) for i in colorsa]
colors.insert(0,'k') 


#legend labels
labels = []
letters = ['a) ','b) ','c) ','d) ']
for j in data_perc:
    labels.append('{:.2f}'.format(float(j)/100))
labels.append('')
handles = []

#Set font family
hfont = {'fontname':'Times'}
plt.rc('font', family='serif',size=24)
plt.rc('text', usetex=True)
plt.rc('legend', fontsize=24)
#mpl.rcParams['text.latex.preamble'] = r'\boldmath'
#plt.rcParams["font.family"] = "Computer Modern"
size_text = 14

# Sample data
q_matrix_files =[]
for i in data_perc:
    q_matrix_files.append([])
    for j in data_r:
        q_matrix_files[-1].append('./'+'16s_u8/dos_green_'+i+'_'+j+'.txt')

#print(q_matrix_files)
# Create subplots with shared x-axis
fig = plt.figure(figsize=(6,3),constrained_layout=True)
L, R, T, B = 0.002, 0.002, 0.0, 0.0
#fig.get_layout_engine().set(w_pad=0.1,h_pad=0.03,wspace=0.,rect=(L, B, 1 - L - R, 1 - B - T))
#mpl.layout_engine.ConstrainedLayoutEngine.set(,w_pad=0.0)
#mpl.layout_engine.ConstrainedLayoutEngine(w_pad=0.0)
#mpl.layout_engine.ConstrainedLayoutEngine.execute(fig)

gs = gridspec.GridSpec(2, 2, figure=fig,width_ratios=[3, 1])
axes = [fig.add_subplot(gs[0,0])]
axes.append(fig.add_subplot(gs[1,0]))
axes.append(fig.add_subplot(gs[:,1]))
#fig, axes = plt.subplots(2, 2, figsize=(6, 3),constrained_layout=True)
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
        x_ref,y_ref = reader(q_files[0])

        for j in range(len(q_matrix_files[i])):
            x,y = reader(q_files[j])
            #SPECTRAL FUNCTIONS
            handle, = ax.plot(x, np.array(y),linestyle=line_style[j],lw=line_width[j],color=colors[j],label=data_r_float_w[j])

            ax.set_ylim(0, float(ylims[i]))
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ax.set_xlim(xlims_min,xlims_max)
            ax.tick_params(axis='x',labelsize=size_text)
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
                #print(ks)
                #chi_2 = ks[0]
                #print(chi_2)


                axes[2].scatter(float(data_r_float[j]),chi_2,color=colors[j], marker=markers[i],s=m_size[i])
                # Move y-axis to the right
                axes[2].yaxis.set_label_position("right")
                axes[2].yaxis.tick_right()
                #ax.set_ylim(0, float(ylims[i]))
                #ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                #ax.set_xlim(xlims_min,xlims_max)
                #axes[2].set_ylim(bottom=-0.4,top = 15.4)
                axes[2].tick_params(axis='x',labelsize=size_text)
                axes[2].tick_params(axis='y', labelsize=size_text)
                axes[2].set_xscale('log')
                #axes[2].set_yscale('log')
                axes[2].minorticks_off()
                axes[2].axhline(0,c='k',lw=0.7,linestyle='--')
                #axes[2].set_xlim(float(data_r_float[-1])-0.02,1.02)
                #axes[2].set_xticks([10e-4,10e-2,])  # Tick only at 0 and 1 and 2
                #axes[2].set_xticklabels([data_r_fl])  # Tick only at 0 and 1 and 2

    # Optionally add individual titles a), b), ...
    if i == 2 : shift = 0.025 
    else :shift = 0
    f = ["$f$=","$f$=",""]
    ax.text(pos_x_name+shift,pos_y_name,letters[i] + f[i]+ labels[i],fontsize=size_text,ha='left',va='top', transform=ax.transAxes)

#Reverse handles
reverse_handle = []
for i in handles:
    reverse_handle.insert(0,i)
handles = reverse_handle

# Set x-axis label on the bottom subplot only
fig.legend(handles=handles,loc='upper left',
        title='$w_t$',title_fontsize=size_text,
        bbox_to_anchor=(0.70,1),
        #bbox_to_anchor=(0.725,0.16),
        #bbox_to_anchor=(0.66,0.5828),
        fontsize=size_text,ncol=1,
        columnspacing=0.5,labelspacing=0.1,
        handletextpad=0.2,handlelength=0.65,borderpad=0.2,
        frameon=True, markerfirst=True)
axes[1].set_xlabel('$\omega$',fontsize = size_text)
axes[2].set_xlabel('1-$w_t$',fontsize = size_text)
axes[2].set_ylabel('$\chi^2$',fontsize = size_text,rotation=0)


# Shared legend placed outside the plots (bottom center)
#plt.tight_layout()
plt.savefig('further_truncation.pdf')

#plt.show()

