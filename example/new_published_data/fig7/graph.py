import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import FancyBboxPatch
import matplotlib
import numpy as np
import matplotlib as mpl
import csv

from scipy.signal import argrelextrema as extrema

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
mh_beta = ['02','02','02']
data_perc = ['10e-5','10e-7','10e-10']
data_perc__ = ['0.00001']
site = ['24','28','32']
data_dim = ['1D']
#data_site = ['16']
ylims = [0.6,0.8,0.5,0.45]



#lines esthetics
line_style = ['solid']*5
#line_style = ['solid','dashed','dashdot',(0,(1,1)),'solid']
#line_width = [2.6,2.2,1.8,1.4,1]
line_width = [1.5]*5
#line_width = [2,2,2,2,1]
#Colors #states
cmap = plt.get_cmap('gnuplot')
colorsa = np.linspace(0.2,0.90,num=len(data_perc));
colors = [cmap(i) for i in colorsa]


#legend labels
labels = []
letters = ['a) ','b) ','c) ','d) ']
for i,j in zip(letters,data_perc):
    labels.append('{:.2f}'.format(float(j)/100))
handles = []
labels = []
labels.append(r'$N_c=24$ ($f\sim10^{-5}$)')
labels.append(r'$N_c=28$ ($f\sim10^{-7}$)')
labels.append(r'$N_c=32$ ($f\sim10^{-10}$)')

#Set font family
hfont = {'fontname':'Times'}
plt.rc('font', family='serif',size=24)
plt.rc('text', usetex=True)
plt.rc('legend', fontsize=24)
#mpl.rcParams['text.latex.preamble'] = r'\boldmath'
#plt.rcParams["font.family"] = "Computer Modern"
size_text = 14
size_text_gap = 9

#Compare Width
c_width = ['width/green_4.txt','width/green_8.txt','width/green_12.txt','width/green_16.txt']

# Sample data
a = ['_02','_05']
q_matrix_files =[]
for j in data_dim:
    q_matrix_files.append([])
    for i,k,l in zip(data_perc,mh_beta,site):
        file = './large_cluster_data/dos_'+j+'_'+l+'green'+i + '_' + k+'.txt'
        q_matrix_files[-1].append(file)
        c_width.append(file)

#Find gap width
left_w, right_w = [], []
for w in c_width : 
    x,y = reader(w)
    a,b = gapWidth(x,y)
    left_w.append(a)
    right_w.append(b)


# Create subplots with shared x-axis
n_subplots = len(data_dim)
fig, axes = plt.subplots(n_subplots, 1, sharex=True, figsize=(6, 1.6+n_subplots),constrained_layout=True)
#plt.subplots_adjust(left=.02,right=.98,bottom=.16,top=.90)  #Change top space
#plt.subplots_adjust(hspace=.1)  # Reduce vertical space between plots
#plt.subplots_adjust(left=0, bottom=0, right=1, top=0, wspace=0, hspace=0)
if n_subplots == 1 :
    axes = [axes]

#BOX PARAMS
ss_x = 0.93#xlims_max*0.50
ss_y = 0.98#*np.array(ylims)
ss_spacing = 0.10
ss_size_text = 12
ha_pos_title = 'center'
ha_pos = 'right'
va_pos = 'top'

box_h = 0.22
box_pad = 0.01

num_slider = [-0.0075,-0.00075,-0.0005,-0.0012]


#print(q_matrix_files)
for i,ax in enumerate(axes):
    q_files = q_matrix_files[i]
    #print(q_files)
    for j in range(len(q_matrix_files[i])):
        #print(f"i:{i}\tj:{j}")
        x,y = reader(q_files[j])
        handle, = ax.plot(x, np.array(y),linestyle=line_style[j],lw=line_width[j],color=colors[j],label=labels[j],zorder=-j)

        ax.set_ylim(0, float(ylims[i]))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax.set_xlim(xlims_min,xlims_max)
        ax.tick_params(axis='x',labelsize=size_text)
        ax.set_yticks([0, float(ylims[i])])  # Tick only at 0 and 1 and 2
        
        #Vert line at 0
        ax.axvline(c='k',ls='--',lw=0.75)
        #ax.text(ss_x+0.073+num_slider[i],
        #        ss_y+box_h/2-box_pad-(j+2)*ss_spacing,
        #        r'$\beta=$'+f'0.{int(mh_beta[j])}',ha=ha_pos,va=va_pos,
        #        fontsize=ss_size_text,transform=ax.transAxes,
        #        color=colors[j])
        
        #Remove y tick labels
        ax.set_yticklabels([])
        if i == 0:  # Add legend only to the first plot
            handles.append(handle)
            #ax.legend(loc='upper right')

    #Add vertical lines
    ns=['4','8','12','16','24','28','32']
    clr = ['dimgray']*4
    for c in colors:
        clr.append(c)

    shift = [0,0,0,0,0.0,0,-0.0]
    space_line = 0.06
    for k,(l,r) in enumerate(zip(left_w,right_w)):
        #if ns[k] == '24':
        #    m = 4
        #else:
        m=k
        y_h_min = (0.95-space_line*m+shift[k])
        y_h_max = (0.99-space_line*m+shift[k])
        tempC = clr[k]
        #print(tempC)
        var = 0.80
        #try:
        #    tempC = (clr[k][0]*var, clr[k][1]*var, clr[k][2]*var, clr[k][3]*var)
        #except:
        #    pass
        #print(tempC)

        ax.axvline(x=l,ymin=y_h_min, ymax=y_h_max, c=tempC, lw=0.76,zorder=10)
        ax.axvline(x=r,ymin=y_h_min, ymax=y_h_max, c=tempC,lw=0.76,zorder=10)
        y_h = (0.97-space_line*m+shift[k])*ylims[i]
        #print(y_h)
        ax.axhline(y=y_h,xmin=(l-xlims_min)/(2*xlims_max),xmax=(r-xlims_min)/(2*xlims_max),c=tempC,lw=0.75,zorder=5)
        #ax.axhline(y=y_h,xmin=(l-xlims_min)/(2*xlims_max),xmax=(-1.6-xlims_min)/(2*xlims_max),c=tempC,lw=0.75,zorder=5)
        #ax.axhline(y=y_h,xmin=(-0.4-xlims_min)/(2*xlims_max),xmax=(r-xlims_min)/(2*xlims_max),c=tempC,lw=0.75,zorder=5)
        text_site = r'$N_c$='+ns[k]
        if k > 0:
            text_site = '\t'+ns[k]
        if m%1==0 :
            ax.text(l-0.1,y_h-0.0035,text_site,fontsize=size_text_gap,va='center',ha='right',c='k')
        else:
            ax.text(r+0.1,y_h-0.0035,text_site,fontsize=size_text_gap,va='center',ha='left',c='k')
        #    ax.text(r,y_h,r'$N_S$='+ns[k],fontsize=ss_size_text,va='bottom',ha='left')

    # Optionally add individual titles or labels
    #ax.set_ylabel(graph_name[i],fontsize=14)    
    #ax.text(pos_x_name,pos_y_name,letters[i] + "$f$=" + data_perc__[i],fontsize=size_text,ha='left',va='top', transform=ax.transAxes)
    #ax.text(pos_x_name,pos_y_name,letters[i],fontsize=size_text,ha='left',va='top', transform=ax.transAxes)
    #ax.text(ss_x+0.073+num_slider[i],
    #            ss_y+box_h/2-box_pad-(1)*ss_spacing,'73 124 592',ha=ha_pos,va=va_pos,
    #            fontsize=ss_size_text,transform=ax.transAxes,
    #            color='k')
    #ax.text(pos_x_name+0.05,pos_y_name-0.025,labels[i],fontsize=size_text-1,ha='left',va='top', transform=ax.transAxes, color=colors[i])


# Set x-axis label on the bottom subplot only
fig.legend(handles=handles,loc='upper right',bbox_to_anchor=(1,1),fontsize=size_text-4,ncol=1,
        columnspacing=0.5,labelspacing=0.1,handletextpad=0.1,handlelength=1,borderpad=0.1, frameon=True)
axes[-1].set_xlabel('$\omega$',fontsize = size_text)

# Shared legend placed outside the plots (bottom center)
#plt.tight_layout()
plt.savefig('large_clusters.pdf')

#plt.show()

