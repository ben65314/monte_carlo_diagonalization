import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import FancyBboxPatch
import matplotlib.ticker as ticker
import matplotlib
import numpy as np
import matplotlib as mpl
import csv

dns4 = {
        "4_100" : f'{36:>2,}'.replace(',',' '),
        "4_50"  : f'{18:>2,}'.replace(',',' '),
        "4_25"  : f'{9:>2,}'.replace(',',' '),
        "4_15"  : f'{6:>2,}'.replace(',',' '),
        "4_10"  : f'{4:>2,}'.replace(',',' '),
        "4_05"  : f'{2:>2,}'.replace(',',' '),
        "8_100" : f'{4_900:>4,}'.replace(',',' '),
        "8_50"  : f'{2_450:>4,}'.replace(',',' '),
        "8_25"  : f'{1_225:>4,}'.replace(',',' '),
        "8_15"  : f'{735:>4,}'.replace(',',' '),
        "8_10"  : f'{490:>4,}'.replace(',',' '),
        "8_05"  : f'{245:>4,}'.replace(',',' '),
        "12_100": f'{853_776:>6,}'.replace(',',' '),
        "12_50" : f'{426_888:>6,}'.replace(',',' '),
        "12_25" : f'{213_444:>6,}'.replace(',',' '),
        "12_15" : f'{128_067:>6,}'.replace(',',' '),
        "12_10" : f'{85_378:>6,}'.replace(',',' '),
        "12_05" : f'{42_689:>6,}'.replace(',',' '),
        "16_100": f'{165_636_900:>9,}'.replace(',',' '),
        "16_50" : f'{82_818_450:>9,}'.replace(',',' '),
        "16_25" : f'{41_409_225:>9,}'.replace(',',' '),
        "16_15" : f'{24_845_535:>9,}'.replace(',',' '),
        "16_10" : f'{16_563_690:>9,}'.replace(',',' '),
        "16_05" : f'{8_281_845:>9,}'.replace(',',' ')    
        }

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
data_site = ['4','8','12','16']
ylims = [2.5,1.3,0.65,0.85]



#lines esthetics
line_style = ['solid']*10
#line_style = ['solid','dashed','dashdot',(0,(1,1)),'solid']
#line_width = [2.6,2.2,1.8,1.4,1]
line_width = [1.5]*10
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
for j in data_site:
    q_matrix_files.append([])
    for i in data_perc:
        q_matrix_files[-1].append('./'+j+'s_u8/dos_green'+i+'.txt')

# Create subplots with shared x-axis
width = 6 #in inches
height = 2 + len(data_site) #in inches
fig, axes = plt.subplots(len(data_site), 1, sharex=True, figsize=(width,height),constrained_layout=True)
#plt.subplots_adjust(left=.02,right=.98,bottom=.16,top=.90)  #Change top space
#plt.subplots_adjust(hspace=.1)  # Reduce vertical space between plots
#plt.subplots_adjust(left=0, bottom=0, right=1, top=0, wspace=0, hspace=0)

#BOX PARAMS
ss_x = 0.93#xlims_max*0.50
ss_y = 1#*np.array(ylims)
ss_spacing = 0.12 #was .15
ss_size_text = 12 #Was 16
ha_pos_title = 'center'
ha_pos = 'right'
va_pos = 'top'

box_h = 0.22
box_pad = 0.01

num_slider = [-0.0075,-0.005,-0.005,-0.005]
nf = [r'$N_f=$',r'$N_f=$',r'$N_f=$',r'$N_f=$']
nf = [0.965,0.93,0.90,0.85]


for i,ax in enumerate(axes):
    q_files = q_matrix_files[i]
    #print(q_files)
    for j in range(len(q_matrix_files[i])):
        x,y = reader(q_files[j])
        handle, = ax.plot(x, np.array(y),linestyle=line_style[j],lw=line_width[j],color=colors[j],label=labels[j])

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
    #ax.text(pos_x_name+0.05,pos_y_name-0.025,labels[i],fontsize=size_text-1,ha='left',va='top', transform=ax.transAxes, color=colors[i])


# Set x-axis label on the bottom subplot only
fig.legend(handles=handles,loc='center',bbox_to_anchor=(0.5,0.96),fontsize=size_text,ncol=6,
        columnspacing=0.5,labelspacing=0.1,handletextpad=0.1,handlelength=1,borderpad=0.1, frameon=True, markerfirst=False)
axes[-1].set_xlabel('$\omega$',fontsize = size_text)
axes[-1].xaxis.set_major_locator(ticker.MultipleLocator(2))

# Shared legend placed outside the plots (bottom center)
#plt.tight_layout()
#plt.savefig('1D_'+"AA"+'s_mu0_wN.svg')
plt.savefig('2D_mu4.pdf')

#plt.show()


#ADD LATTICE 
#fig.savefig('./cover.svg', transparent=False)

#xPtLat = 15
#yPtLat = 30
#yPtDiff = 122
#width_svg = str(width*90)
#height_svg = str(height*90)
#print(width_svg,'\t',height_svg)
#sc.Figure(width_svg,height_svg,
#        sc.SVG('./images/cover.svg').scale(1.25),
#        sc.SVG("./images/2_2_lattice.svg").scale(0.14).move(xPtLat,yPtLat),
#        sc.SVG("./images/4_2_lattice.svg").scale(0.14).move(xPtLat,yPtLat+1*yPtDiff),
#        sc.SVG("./images/4_3_lattice.svg").scale(0.14).move(xPtLat,yPtLat+2*yPtDiff),
#        sc.SVG("./images/4_4_lattice.svg").scale(0.14).move(xPtLat,yPtLat+3*yPtDiff)
#        ).save("2D_4s_to_16s_wl.svg")



