import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import FancyBboxPatch
import matplotlib
import numpy as np
import matplotlib as mpl
import csv

dns0 = {
        "16_100": '19 079 424',
        "16_43" : '8 281 845 ',
        "16_05" : '953 972   '    
        }

def norm(data):
    return (data)/(max(data)-min(data))

def reader(file_name):
    x = []
    y = [[],[],[]]
    #y = []
    with open(file_name) as file:
        reader = csv.reader(file,delimiter='\t')
        for lines in reader:
            if(lines[0].startswith('#') or lines[0].startswith('PARAMETERS') or lines[0].startswith('NEXT_SITE')):
                continue
            x.append(float(lines[0]))
            #y.append(float(lines[1])+float(lines[2]))
            y[0].append(float(lines[1])+float(lines[2]))
            y[1].append(float(lines[1]))
            y[2].append(float(lines[2]))
    return x,y

##Curves information
graph_name = [r'$A(\omega)$',r'$A^+(\omega)$',r'$A^-(\omega)$']
#limits
ylims = [1,1,1,1]
xlims_min = -3
xlims_max = 12
#position of # of states abs position
pos_x_name = 0.005
pos_y_name = 0.975

#data
#JUST NEED TO CHANGE THOSE VALUES, works up to 4 graphs 
data_perc = ['100','43','05']
system_size = str(16)



#lines esthetics
line_style = ['--','-','-']
#line_style = ['solid','dashed','dashdot',(0,(1,1)),'solid']
#line_width = [2.6,2.2,1.8,1.4,1]
line_width = [2.5,1.75,1.75]
#line_width = [2,2,2,2,1]
#Colors #states
cmap = plt.get_cmap('gnuplot')
colorsa = np.linspace(0,0.90,num=len(data_perc));
#colors_states = [cmap(i) for i in colorsa]
colors_states = ['k' for i in colorsa]


#Colors warm cold
cmap = plt.get_cmap('coolwarm')
colorsa = 0.1
colors = ['k',cmap(0+colorsa),cmap(1-colorsa)]
alpha = [1,0.75,0.75]
#legend labels
labels = []
letters = ['a) ','b) ','c) ','d) ']
for i,j in zip(letters,data_perc):
    labels.append('{:.2f}'.format(float(j)/100))
handles = []

#Set font family
hfont = {'fontname':'Times'}
plt.rc('font', family='serif',size=20)
plt.rc('text', usetex=True)
plt.rc('legend', fontsize=20)
#mpl.rcParams['text.latex.preamble'] = r'\boldmath'
#plt.rcParams["font.family"] = "Computer Modern"
size_text = 14


# Sample data
q_matrix_files =[]
for i in data_perc:
    q_matrix_files.append('./'+system_size+'s_u8/dos_green'+i+'.txt')

# Create subplots with shared x-axis
fig, axes = plt.subplots(len(data_perc), 1, sharex=True, figsize=(6, 2+len(data_perc)),constrained_layout=True)
#plt.subplots_adjust(left=.02,right=.98,bottom=.16,top=.90)  #Change top space
#plt.subplots_adjust(hspace=.1)  # Reduce vertical space between plots
#plt.subplots_adjust(left=0, bottom=0, right=1, top=0, wspace=0, hspace=0)

#BOX PARAMS
ss_x = 0.92#xlims_max*0.50
ss_y = 0.96#*np.array(ylims)
ss_spacing = 0.1
ss_size_text = 12
ha_pos_title = 'center'
ha_pos = 'right'
va_pos = 'top'

box_h = 0.22
box_pad = 0.01

num_slider = [0.0005,-0.0005,-0.0005,-0.0012]

nf = [0.86,0.87,0.89,0.85]

for i,ax in enumerate(axes):
    q_files = q_matrix_files[i]
    #print(q_files)
    x,y = reader(q_files)
    for j in range(3):
        handle, = ax.plot(x, np.array(y[j]),linestyle=line_style[j],lw=line_width[j],color=colors[j],label=graph_name[j],alpha=alpha[j])

        ax.set_ylim(0, float(ylims[i]))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax.set_xlim(xlims_min,xlims_max)
        ax.tick_params(axis='x',labelsize=size_text)
        ax.set_yticks([0, ylims[i]])  # Tick only at 0 and 1 and 2
        
        #Vert line at 0
        ax.axvline(c='k',ls='--',lw=0.75)
        
        if j == 0 :
            ax.text(nf[i],
                    ss_y+box_h/2-box_pad-(j+1)*ss_spacing,
                    r'$N_f=$',ha=ha_pos,va=va_pos,
                    fontsize=ss_size_text,transform=ax.transAxes,
                    color=colors[j])
            ax.text(ss_x+0.073+num_slider[i],
                    ss_y+box_h/2-box_pad-(j+1)*ss_spacing,
                    dns0[system_size+'_'+data_perc[i]],ha=ha_pos,va=va_pos,
                    fontsize=ss_size_text,transform=ax.transAxes,
                    color=colors_states[i])
        #Remove y tick labels
        ax.set_yticklabels([])
        if i == 0:  # Add legend only to the first plot
            handles.append(handle)
            #ax.legend(loc='upper right')
    # Optionally add individual titles or labels
    #ax.set_ylabel(graph_name[i],fontsize=14)    
    ax.text(pos_x_name,pos_y_name,letters[i] + "$f$=" + labels[i],fontsize=size_text,ha='left',va='top', transform=ax.transAxes)
    #ax.text(pos_x_name+0.05,pos_y_name-0.025,labels[i],fontsize=size_text,ha='left',va='top', transform=ax.transAxes, color=colors_states[i])


# Set x-axis label on the bottom subplot only
fig.legend(handles=handles,loc='center',bbox_to_anchor=(0.5,0.95),fontsize=size_text,ncol=4,
        columnspacing=0.5,labelspacing=0.1,handletextpad=0.1,handlelength=1.5,borderpad=0, frameon=False)
axes[-1].set_xlabel('$\omega$',fontsize = size_text)

# Shared legend placed outside the plots (bottom center)
#plt.tight_layout()
plt.savefig('1D_'+system_size+'s_mu0.pdf')

#plt.show()
