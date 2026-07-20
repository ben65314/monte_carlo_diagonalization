import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import FancyBboxPatch
import matplotlib.ticker as ticker
import matplotlib
import numpy as np
import matplotlib as mpl
import csv

from PIL import Image
import svgutils.compose as sc

import sys
import os
# add the parent directory (where stack_axes.py lives) to the path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from stack_axes import make_stacked_axes

dns4 = {
        "4_100" : f'{         36:>10,}'.replace(',',' '),
        "4_50"  : f'{         18:>10,}'.replace(',',' '),
        "4_25"  : f'{          9:>10,}'.replace(',',' '),
        "4_15"  : f'{          6:>10,}'.replace(',',' '),
        "4_10"  : f'{          4:>10,}'.replace(',',' '),
        "4_05"  : f'{          2:>10,}'.replace(',',' '),
        "8_100" : f'{      4_900:>10,}'.replace(',',' '),
        "8_50"  : f'{      2_450:>10,}'.replace(',',' '),
        "8_25"  : f'{      1_225:>10,}'.replace(',',' '),
        "8_15"  : f'{        735:>10,}'.replace(',',' '),
        "8_10"  : f'{        490:>10,}'.replace(',',' '),
        "8_05"  : f'{        245:>10,}'.replace(',',' '),
        "12_100": f'{    853_776:>10,}'.replace(',',' '),
        "12_50" : f'{    426_888:>10,}'.replace(',',' '),
        "12_25" : f'{    213_444:>10,}'.replace(',',' '),
        "12_15" : f'{    128_067:>10,}'.replace(',',' '),
        "12_10" : f'{     85_378:>10,}'.replace(',',' '),
        "12_05" : f'{     42_689:>10,}'.replace(',',' '),
        "16_100": f'{165_636_900:>10,}'.replace(',',' '),
        "16_50" : f'{ 82_818_450:>10,}'.replace(',',' '),
        "16_25" : f'{ 41_409_225:>10,}'.replace(',',' '),
        "16_15" : f'{ 24_845_535:>10,}'.replace(',',' '),
        "16_10" : f'{ 16_563_690:>10,}'.replace(',',' '),
        "16_05" : f'{  8_281_845:>10,}'.replace(',',' ')    
        }

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
data_perc = ['100','50','25','15','05']
data_site = ['4','8','12','16']
ylims = [2.5,1.3,0.65,0.85]
#legend labels
labels = []
letters = ['a) ','b) ','c) ','d) ','e) ']
for j in data_perc:
    labels.append('{:.2f}'.format(float(j)/100))
labels[0] = '$f$: {:.2f}'.format(float(data_perc[0])/100)
handles = []


#lines esthetics
line_style = ['solid']*10
line_width = [1.5]*10

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


# Get data file
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
ss_spacing = 0.12 #was .15
ss_size_text = 14 #Was 16
ha_pos = 'right'
va_pos = 'top'

box_h = 0.22
box_pad = 0.01

num_slider = [-0.0075,-0.005,-0.005,-0.005,0]
nf = [r'$N_f=$',r'$N_f=$',r'$N_f=$',r'$N_f=$']
nf = [0.96,0.92,0.89,0.835]


for i,ax in enumerate(axes):
    q_files = q_matrix_files[i]
    #print(q_files)
    for j in range(len(q_matrix_files[i])):
        x,y = reader(q_files[j])
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
        ax.text(ss_x + 0.073 + num_slider[i],
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
    ax.text(pos_x_name,pos_y_name,letters[i] + r'$N_c=$ ' + data_site[i],fontsize=size_text,ha='left',va='top', transform=ax.transAxes)


# Set x-axis label on the bottom subplot only
fig.legend(handles=handles,loc='center',bbox_to_anchor=(0.5,0.96),fontsize=size_text,ncol=6,
        columnspacing=0.7,labelspacing=0.1,handletextpad=0.1,handlelength=1,borderpad=0.1, frameon=True, markerfirst=False)
axes[-1].set_xlabel('$\omega$',fontsize = size_text)
axes[-1].xaxis.set_major_locator(ticker.MultipleLocator(2))

# Shared legend placed outside the plots (bottom center)
plt.show()


#ADD LATTICE 
fig.savefig('./images/cover.pdf', transparent=False)

import fitz
def place(page, pdf_path, x, y, scale, page_num=0):
    src = fitz.open(pdf_path)
    r = src[page_num].rect
    target = fitz.Rect(x, y, x + r.width*scale, y + r.height*scale)
    page.show_pdf_page(target, src, page_num)
    src.close()

width, height = fig.get_size_inches()
xPtLat = 20
yPtLat = 20
yPtDiff = 99
width_svg = str(width*90)
height_svg = str(height*90)


width_pdf = width * 72
height_pdf = height * 72

doc = fitz.open()
page = doc.new_page(width=width_pdf, height=height_pdf)

place(page, './images/cover.pdf', 0, 0, 1)
place(page, './images/2_2_lattice.pdf', xPtLat, yPtLat, 0.14)
place(page, './images/4_2_lattice.pdf', xPtLat, yPtLat + 1*yPtDiff, 0.14)
place(page, './images/4_3_lattice.pdf', xPtLat, yPtLat + 2*yPtDiff, 0.14)
place(page, './images/4_4_lattice.pdf', xPtLat, yPtLat + 3*yPtDiff, 0.14)

doc.save("2D_mu4.pdf")
doc.close()
