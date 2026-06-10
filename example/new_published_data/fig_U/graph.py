import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import FancyBboxPatch
import matplotlib
import numpy as np
import matplotlib as mpl
import csv

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
xlims_min = -9
xlims_max = 9
#position of # of states abs position
pos_x_name = 0.005
pos_y_name = 0.975

#data
#JUST NEED TO CHANGE THOSE VALUES, works up to 4 graphs 
data_perc = ['100','05']#['100','50','25','15','05']
data_perc_label = ['100','05']#['100','50','25','15','05']
perc_label = ['ED',r'$f=0.05$',r'$U\beta=0.2$']#['100','50','25','15','05']
data_site = 16
data_u = ['4','8','12']
ylims = [0.6,0.4,0.4,0.4,0.4]



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
letters = ['a) ','b) ','c) ','d) ','e) ']
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
for j in data_u:
    q_matrix_files.append([])
    for i in data_perc_label:
        q_matrix_files[-1].append('./u'+j+'/dos_green'+i+'.txt')

# Create subplots with shared x-axis
fig, axes = plt.subplots(len(data_u), 1, sharex=True, figsize=(6, 2+len(data_u)),constrained_layout=True)
#plt.subplots_adjust(left=.02,right=.98,bottom=.16,top=.90)  #Change top space
#plt.subplots_adjust(hspace=.1)  # Reduce vertical space between plots
#plt.subplots_adjust(left=0, bottom=0, right=1, top=0, wspace=0, hspace=0)
#new_axes = []
#for i in axes:
#    for j in i:
#        new_axes.append(j)
#axes=new_axes
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
    print(q_files)
    for j in range(len(q_matrix_files[i])):

        x,y = reader(q_files[j])
        handle, = ax.plot(x, np.array(y),linestyle=line_style[j],lw=line_width[j],color=colors[j],label=perc_label[j])
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
    ax.text(pos_x_name,pos_y_name,letters[i] + r'$U/t=$ ' + data_u[i],fontsize=size_text,ha='left',va='top', transform=ax.transAxes)
    #ax.text(pos_x_name,pos_y_name,letters[i] + data_site[i] + ' sites',fontsize=size_text,ha='left',va='top', transform=ax.transAxes)
    #ax.text(pos_x_name+0.05,pos_y_name-0.025,labels[i],fontsize=size_text-1,ha='left',va='top', transform=ax.transAxes, color=colors[i])



# Set x-axis label on the bottom subplot only
legend = fig.legend(handles=handles,loc='center',bbox_to_anchor=(0.5,0.975),fontsize=size_text,ncol=5,
        columnspacing=2.,labelspacing=0.1,handletextpad=0.1,handlelength=1,borderpad=0.1, frameon=True,markerfirst=False)
fig.get_layout_engine().set(rect=[0, 0, 1, 0.955])  # top 5% reserved for legend
axes[-1].set_xlabel('$\omega$',fontsize = size_text)

# Shared legend placed outside the plots (bottom center)
#plt.tight_layout()
##plt.savefig('1D_4s_to_16s.svg')
plt.savefig('var_U.pdf', bbox_inches='tight', bbox_extra_artists=[legend])

plt.show()

