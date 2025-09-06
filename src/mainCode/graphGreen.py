#!/usr/bin/python3
import csv
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import sys;

COMPUTE_AVERAGE = True
is_both_data = False

colors = ['black','red','green','orange','blue','yellow','cyan','lime','orange','blue','black','red','green','orange','blue']
def fg(x,a,b,c,d):
    return a*(x**3)+b*(x**2)+c*x+d

def f(x,a,b,c,d):
    return a*np.log(b*(x+d))+c

if (len(sys.argv) == 1): exit

#Number of files to compare
n_files = len(sys.argv) - 1

#File name
data_file = sys.argv[1]
#dataTAG = ['Reticle = 1','Reticle = 5','Reticle = 20','Reticle = 50']

#Recherche du fichier ayant les données nécessaire
#Reference
for i in range(n_files):
    linestyle = '-'
    data = [[],[],[]];
    both_data = [];
    index = 1;
    omegaRead = False
    param = []
    #data.append([])
    data_file_name = sys.argv[1+i]
    with open(data_file_name) as file:
        reader = csv.reader(file,delimiter='\t')
        for lines in reader:
            #Ignore les lignes de commentaire contenues dans le fichier 
            if(lines[0].startswith('#')):
                continue;
            if(lines[0].startswith('NEXT_SITE')):
                if is_both_data : 
                    plt.plot(data[0],data[1],lw=1,ls='dashed', c=colors[i%15])
                    plt.plot(data[0],data[2],lw=1,ls='dotted', c=colors[i%15])
                
                combined = np.array(data[1])+np.array(data[2]) - (index-1) 
                if index == 1 :
                    plt.plot(data[0],combined,lw=1, c=colors[i%15],label = data_file_name)
                else:
                    plt.plot(data[0],combined,lw=1, c=colors[i%15])

                data[1].clear()
                data[2].clear()
                index+=1
                omegaRead = True
                ##data = [[],[]]
                continue
            if(lines[0].startswith('PARAMETERS')):
                param.append(lines[1])
                param.append(lines[2]) 
                param.append(lines[3])
                param.append(lines[4])
                param.append(lines[5])

            else :
                if not omegaRead :
                    data[0].append(float(lines[0]))
                data[1].append(float(lines[1]))
                if len(lines) > 2 : data[2].append(float(lines[2]))
                else : data[2].append(index-1)

       
title = "#Sites = " +str(param[0])+"\tU = "+str(param[1])+"\t$\mu$ = "+str(param[2])+"\tN = "+str(param[3])+"\t$S_z$ = "+str(param[4])
plt.title(title)

#ax.legend()
plt.grid()
#Affichage
plt.legend()
plt.show()




