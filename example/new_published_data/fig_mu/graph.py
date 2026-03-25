import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import FancyBboxPatch
import matplotlib
import numpy as np
import matplotlib as mpl
import csv

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
xlims_max = 12
#position of # of states abs position
pos_x_name = 0.005
pos_y_name = 0.975

#data
#JUST NEED TO CHANGE THOSE VALUES, works up to 4 graphs
data_perc = ['100','Ncst']
data_site = 16
data_mu = ['2','1_6','1_2','0_8','0_4']
data_mu_label = ['2','1.6','1.2','0.8','0.4']
ylims = [0.35,0.35,0.6,0.6,0.6,0.6]


ref_number = ['165,636,900', '165,636,900', '130,873,600', '64,128,064', '19,079,424']
per_5 = ['8,281,845', '8,281,845', '6,543,680', '3,206,403', '953,972']
per_Ncst = ['(f=0.05)','(f=0.05)','(f=0.06)','(f=0.13)','(f=0.43)']
Ncst = "8,281,845"



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
labels = ['$N_{1.00}$','$N_{0.05}$','Ncst']
letters = ['a) ','b) ','c) ','d) ','e) ','f) ']
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
for j in data_mu:
    q_matrix_files.append([])
    for i in data_perc:
        if (i == '05') :
            continue
        q_matrix_files[-1].append('./mu'+j+'/average/av_green_'+i+'.txt')

# Create subplots with shared x-axis
fig, axes = plt.subplots(len(data_mu), 1, sharex=True, figsize=(6, 2+len(data_mu)),constrained_layout=True)
#plt.subplots_adjust(left=.02,right=.98,bottom=.16,top=.90)  #Change top space
#plt.subplots_adjust(hspace=.1)  # Reduce vertical space between plots
#plt.subplots_adjust(left=0, bottom=0, right=1, top=0, wspace=0, hspace=0)

#BOX PARAMS
ss_x = 0.9975#xlims_max*0.50
ss_y = 0.99#*np.array(ylims)
ss_spacing = 0.14 #Was .15
ss_size_text = 12 #Was 16
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
            handle1, = ax.plot(x, np.array(y_both)[:,0],linestyle='--',lw=1.,label=r"$G^-(\omega)$",color='red', alpha=.75)
            handle2, = ax.plot(x, np.array(y_both)[:,1],linestyle='--',lw=1.,label=r"$G^+(\omega)$",color='blue', alpha=.75)
        ax.set_ylim(0, float(ylims[i]))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax.set_xlim(xlims_min,xlims_max)
        ax.tick_params(axis='x',labelsize=size_text)
        ax.set_yticks([0, ylims[i]])  # Tick only at 0 and 1 and 2

        #Vert line at 0
        ax.axvline(c='k',ls='--',lw=0.75)

        #Gap width
        #if i == 0 :
        #    c=0.64
        #else:
        #    c=0.74
        #if j == 0:
        #    l,r = gapWidth(x,y)
        #    ax.axvline(x=l,ymin=c-0.02, ymax=(c+0.02), c='gray', lw=1,zorder=-10)
        #    ax.axvline(x=r,ymin=c-0.02, ymax=(c+0.02), c='gray', lw=1,zorder=-10)
        #    y_h = (c)*ylims[i]
        #    print(y_h)
        #    ax.axhline(y=y_h,xmin=(l-xlims_min)/(2*xlims_max),xmax=(r-xlims_min)/(2*xlims_max),c='gray',lw=1)



        #Remove y tick labels
        ax.set_yticklabels([])
        if i == 0:  # Add legend only to the first plot
            
            handles.append(handle)
            #ax.legend(loc='upper right')
    # Optionally add individual titles or labels
    #ax.set_ylabel(graph_name[i],fontsize=14)
    ax.text(pos_x_name,pos_y_name,letters[i] + r'$\mu=$ ' + data_mu_label[i],fontsize=size_text,ha='left',va='top', transform=ax.transAxes)
    ax.text(ss_x,ss_y, ref_number[i], fontsize=size_text,ha='right',va='top', transform=ax.transAxes,c=colors[0])
    #ax.text(ss_x,ss_y - ss_spacing, per_5[i], fontsize=size_text,ha='right',va='top', transform=ax.transAxes,c=colors[1])
    ax.text(ss_x,ss_y - 1*ss_spacing, Ncst + per_Ncst[i], fontsize=size_text,ha='right',va='top', transform=ax.transAxes,c=colors[1])
    #ax.text(pos_x_name+0.05,pos_y_name-0.025,labels[i],fontsize=size_text-1,ha='left',va='top', transform=ax.transAxes, color=colors[i])


handles.append(handle1)
handles.append(handle2)

# Set x-axis label on the bottom subplot only
fig.legend(handles=handles,loc='center',bbox_to_anchor=(0.45,0.96),fontsize=size_text-2,ncol=2,columnspacing=0.5,labelspacing=0.1,handletextpad=0.1,handlelength=1,borderpad=0.1, frameon=True,markerfirst=False)
axes[-1].set_xlabel('$\omega$',fontsize = size_text)

# Shared legend placed outside the plots (bottom center)
#plt.tight_layout()
##plt.savefig('1D_4s_to_16s.svg')
plt.savefig('mu_variation.pdf')

#plt.show()

