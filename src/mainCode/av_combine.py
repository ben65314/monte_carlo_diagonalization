#!/usr/bin/python3
import numpy as np
import sys;
import csv

#File name
data_file = sys.argv[1]
out_file = "av_" + data_file

x = []
y_data_e = [[]] 
y_data_h = [[]] 
index = 0
xRead = False

param = []

#Reading data
with open(data_file) as file:
    reader = csv.reader(file,delimiter='\t')
    for lines in reader:
        #Ignores comments
        if lines[0].startswith('#') :
            continue;
        elif lines[0].startswith('NEXT_SITE') :
            index += 1
            y_data_e.append([])
            y_data_h.append([])
            xRead  = True


        elif lines[0].startswith('PARAMETERS') :
            param.append(int(lines[1]))
            param.append(float(lines[2])) 
            param.append(float(lines[3]))
            param.append(int(lines[4]))
            param.append(int(lines[5]))
            param.append(float(lines[6]))

        else :
            if not xRead :
                x.append(float(lines[0]))

            y_data_e[index].append(float(lines[1]))
            y_data_h[index].append(float(lines[2]))

    y_data_e.pop(-1)
    y_data_h.pop(-1)

max_shift = -len(y_data_e)*(len(y_data_e) - 1) / 2
av_data_e = np.array([max_shift for i in range(len(x))])
av_data_h = np.array([max_shift for i in range(len(x))])

# Sums the average of the curves
for i, (y_e, y_h) in enumerate(zip(y_data_e, y_data_h)):
    av_data_e += np.array(y_e)
    av_data_h += np.array(y_h)

av_data_e /= len(y_data_e)
av_data_h /= len(y_data_h)

# Writing
lineString = ("PARAMETERS\t" + str(param[0]) + "\t" + str(param[1]) + "\t" 
            + str(param[2]) + "\t" + str(param[3]) + "\t" + str(param[4]) 
            + "\t" + str(param[5]) + "\n#OMEGA\tVALUE\n")

file = open(out_file,'w')
for x, y_e, y_h in zip(x, av_data_e, av_data_h):
    lineString += str(x) + "\t" + str(y_e) + "\t" + str(y_h) + "\n"
lineString += 'NEXT_SITE'
file.write(lineString)

