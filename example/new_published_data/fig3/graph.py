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


#Each name tag is label numberOfSites_Percentage
dns4 = {
        "4_100" : '36         ',
        "4_50"  : '18         ',
        "4_25"  : '9          ',
        "4_15"  : '6          ',
        "4_10"  : '4          ',
        "4_05"  : '2          ',
        "8_100" : '4 900      ',
        "8_50"  : '2 450      ',
        "8_25"  : '1 225      ',
        "8_15"  : '735        ',
        "8_10"  : '490        ',
        "8_05"  : '245        ',
        "12_100": '853 776    ',
        "12_50" : '426 888    ',
        "12_25" : '213 444    ',
        "12_15" : '128 067    ',
        "12_10" : '85 378     ',
        "12_05" : '42 689     ',
        "16_100": '165 636 900',
        "16_50" : '82 818 450 ',
        "16_25" : '41 409 225 ',
        "16_15" : '24 845 535 ',
        "16_10" : '16 563 690 ',
        "16_05" : '8 281 845  '
        }

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
            #y.append(float(lines[1])+float(lines[2]))
            y.append(float(lines[1])+float(lines[2]))
            y_both.append([float(lines[1]),float(lines[2])])
    return x,y,y_both

##Curves information
#limits
xlims_min = -7
xlims_max = 7
#position of letters of subfig
pos_x_name = 0.005
pos_y_name = 0.975

#data 
data_perc = ['100','50','25','15','05']
data_site = ['4','8','12','16']
ylims = [2.5,0.8,0.5,0.45]
#legend labels
labels = []
letters = ['a) ','b) ','c) ','d) ']
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


# Get data files
q_matrix_files =[]
for j in data_site:
    q_matrix_files.append([])
    for i in data_perc:
        q_matrix_files[-1].append('./'+j+'s_u8/dos_green'+i+'.txt')

# Create subplots with shared x-axis
fig,axes = make_stacked_axes(n_panels=4,panel_height_cm=3.2,left=0.01,right=0.99,top_margin_cm=0.1)

#BOX PARAMS
ss_x = 0.93#xlims_max*0.50
ss_y = 0.975#*np.array(ylims)
ss_spacing = 0.12 #Was .15
ss_size_text = 14 #Was 16
ha_pos_title = 'center'
ha_pos = 'right'
va_pos = 'top'

box_h = 0.22
box_pad = 0.01

num_slider = [-0.0075,-0.005,-0.005,-0.005,0]

nf = [0.96,0.92,0.89,0.835]

for i,ax in enumerate(axes):
    q_files = q_matrix_files[i]
    #print(q_files)
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



        if j == 0:
            ax.text(nf[i],
                    ss_y+box_h/2-box_pad-(j+1)*ss_spacing,
                    r'$N_f=$',ha=ha_pos,va=va_pos,
                    fontsize=ss_size_text,transform=ax.transAxes,
                    color=colors[j])

        ax.text(ss_x+0.073+num_slider[i],
                ss_y+box_h/2-box_pad-(j+1)*ss_spacing,
                dns4[data_site[i]+'_'+data_perc[j]],ha=ha_pos,va=va_pos,
                fontsize=ss_size_text,transform=ax.transAxes,
                color=colors[j])
        #Remove y tick labels
        ax.set_yticklabels([])
        if i == 0:  # Add legend only to the first plot
            handles.append(handle)
            #ax.legend(loc='upper right')
    # Optionally add individual titles or labels
    #ax.set_ylabel(graph_name[i],fontsize=14)
    ax.text(pos_x_name,pos_y_name,letters[i] + r'$N_c=$ ' + data_site[i],fontsize=size_text,ha='left',va='top', transform=ax.transAxes)
    #ax.text(pos_x_name,pos_y_name,letters[i] + data_site[i] + ' sites',fontsize=size_text,ha='left',va='top', transform=ax.transAxes)
    #ax.text(pos_x_name+0.05,pos_y_name-0.025,labels[i],fontsize=size_text-1,ha='left',va='top', transform=ax.transAxes, color=colors[i])



# Set x-axis label on the bottom subplot only
legend = fig.legend(handles=handles,loc='center',bbox_to_anchor=(0.5,0.96),
                    fontsize=size_text,ncol=5,columnspacing=0.7,labelspacing=0.1,
                    handletextpad=0.1,handlelength=1,borderpad=0.1, frameon=True,
                    markerfirst=False)
#fig.get_layout_engine().set(rect=[0, 0, 1, 1])  # top 5% reserved for legend
axes[-1].set_xlabel('$\omega$',fontsize = size_text)

# Shared legend placed outside the plots (bottom center)
#plt.tight_layout()
plt.savefig('1D_4s_to_16s.pdf',bbox_inches='tight',bbox_extra_artists=[legend])

plt.show()

