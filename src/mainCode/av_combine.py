#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import sys;
import csv

#File name
data_file = sys.argv[1]
out_file = "av_" + data_file

x = []
y_data = [[]] 
y_data2 = [[]] 
index = 0
xRead = False

param = []

#Reading data
with open(data_file) as file:
    reader = csv.reader(file,delimiter='\t')
    for lines in reader:
        #Ignores comments
        if(lines[0].startswith('#')):
            continue;
        elif(lines[0].startswith('NEXT_SITE')):
            index += 1
            y_data.append([])
            y_data2.append([])
            xRead  = True


        elif(lines[0].startswith('PARAMETERS')):
            param.append(int(lines[1]))
            param.append(float(lines[2])) 
            param.append(float(lines[3]))
            param.append(int(lines[4]))
            param.append(int(lines[5]))
            param.append(float(lines[6]))

        else :
            if not xRead :
                x.append(float(lines[0]))

            y_data[index].append(float(lines[1]))
            y_data2[index].append(float(lines[2]))

    y_data.pop(-1)
    y_data2.pop(-1)

max_shift = -len(y_data)*(len(y_data)-1)/2
av_data = np.array([max_shift for i in range(len(x))])
av_data2 = np.array([max_shift for i in range(len(x))])
for i,(y,y2) in enumerate(zip(y_data,y_data2)):
    av_data += np.array(y)
    av_data2 += np.array(y2)
    plt.plot(x,y)
av_data/=len(y_data)
av_data2/=len(y_data2)

print(param)
lineString = "PARAMETERS\t" + str(param[0]) + "\t" + str(param[1]) + "\t" + str(param[2]) + "\t" + str(param[3]) + "\t" + str(param[4]) + "\t" + str(param[5]) + "\n#OMEGA\tVALUE\n";
file = open(out_file,'w')
for x,y,y2 in zip(x,av_data,av_data2):
    lineString += str(x)+"\t"+str(y)+"\t"+str(y2)+"\n"
lineString += 'NEXT_SITE'
file.write(lineString)

