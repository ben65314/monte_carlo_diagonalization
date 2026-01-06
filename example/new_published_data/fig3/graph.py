import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import FancyBboxPatch
import matplotlib
import numpy as np
import matplotlib as mpl
import csv

from scipy.signal import argrelextrema as extrema

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

def integrate(y, dx):
    sum = 0.0
    for i in range(len(y)-1):
        sum += (y[i]+y[i+1])*dx/2
    return sum

def norm(data):
    return (data)/(max(data)-min(data))

def reader(file_name):
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
    return x,y

def gapWidth(x_data,y_data):
    n_points = len(y_data)

    left_x = x_data[:n_points//2]
    right_x = x_data[n_points//2:]
    
    left_y = y_data[:n_points//2]
    right_y = y_data[n_points//2:]

    left_index = extrema(np.array(left_y),np.greater)
    right_index = extrema(np.array(right_y),np.greater)

    max_left = left_index[0][-1]
    max_right = right_index[0][0]

    max_x_left = left_x[max_left]
    max_x_right = right_x[max_right]

    return max_x_left, max_x_right
##Curves information
#limits
xlims_min = -7
xlims_max = 7
#position of # of states abs position
pos_x_name = 0.005
pos_y_name = 0.975

#data
#JUST NEED TO CHANGE THOSE VALUES, works up to 4 graphs 
data_perc = ['100','50','25','15','05']
data_site = ['4','8','12','16']
ylims = [2.5,0.8,0.5,0.45]



#lines esthetics
line_style = ['solid']*5
#line_style = ['solid','dashed','dashdot',(0,(1,1)),'solid']
#line_width = [2.6,2.2,1.8,1.4,1]
line_width = [1.5]*5
#line_width = [2,2,2,2,1]
#Colors #states
cmap = plt.get_cmap('gnuplot')
colorsa = np.linspace(0,0.90,num=len(data_perc));
colors = [cmap(i) for i in colorsa]


#legend labels
labels = []
letters = ['a) ','b) ','c) ','d) ']
for j in data_perc:
    labels.append('{:.2f}'.format(float(j)/100))
labels[0] = '$f$: {:.2f}'.format(float(data_perc[0])/100)
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
for j in data_site:
    q_matrix_files.append([])
    for i in data_perc:
        q_matrix_files[-1].append('./'+j+'s_u8/dos_green'+i+'.txt')

# Create subplots with shared x-axis
fig, axes = plt.subplots(len(data_site), 1, sharex=True, figsize=(6, 2+len(data_site)),constrained_layout=True)
#plt.subplots_adjust(left=.02,right=.98,bottom=.16,top=.90)  #Change top space
#plt.subplots_adjust(hspace=.1)  # Reduce vertical space between plots
#plt.subplots_adjust(left=0, bottom=0, right=1, top=0, wspace=0, hspace=0)

#BOX PARAMS
ss_x = 0.93#xlims_max*0.50
ss_y = 1#*np.array(ylims)
ss_spacing = 0.12 #Was .15
ss_size_text = 12 #Was 16
ha_pos_title = 'center'
ha_pos = 'right'
va_pos = 'top'

box_h = 0.22
box_pad = 0.01

num_slider = [-0.0075,-0.005,-0.005,-0.005,0]

nf = [0.965,0.93,0.90,0.85]

for i,ax in enumerate(axes):
    q_files = q_matrix_files[i]
    #print(q_files)
    for j in range(len(q_matrix_files[i])):

        x,y = reader(q_files[j])
        handle, = ax.plot(x, np.array(y),linestyle=line_style[j],lw=line_width[j],color=colors[j],label=labels[j])
        #print(integrate(y,x[1]-x[0]))
        ax.set_ylim(0, float(ylims[i]))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax.set_xlim(xlims_min,xlims_max)
        ax.tick_params(axis='x',labelsize=size_text)
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
fig.legend(handles=handles,loc='center',bbox_to_anchor=(0.5,0.96),fontsize=size_text,ncol=5,
        columnspacing=0.5,labelspacing=0.1,handletextpad=0.1,handlelength=1,borderpad=0.1, frameon=True,markerfirst=False)
axes[-1].set_xlabel('$\omega$',fontsize = size_text)

# Shared legend placed outside the plots (bottom center)
#plt.tight_layout()
plt.savefig('1D_4s_to_16s.pdf')

#plt.show()

