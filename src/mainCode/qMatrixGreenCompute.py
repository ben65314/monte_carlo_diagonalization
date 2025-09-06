#!/usr/bin/python3
import csv
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import sys;

set_omega = False
fitting = False;
fileNameQM = "qMatrices.txt";
if (len(sys.argv) == 4):
    min_omega = float(sys.argv[2])
    max_omega = float(sys.argv[3])
    set_omega = True 
    fileNameQM = sys.argv[1]
elif len(sys.argv) > 2 :
    print("Number of given arguments incorrect")
elif (len(sys.argv) == 2):
    fileNameQM = sys.argv[1];
#elif (len(sys.argv) == 3):
#    min_omega = float(sys.argv[1])
#    max_omega = float(sys.argv[2])
#    print(type(min_omega))
#    set_omega = True 

def gF_electrons(omega, fundE, QMatrixE, eigenE, nu_mu):
    temp = 0;
    z = omega + eta*1j;

    #Electron part
    for l in range(len(eigenE)):
        denom_e = z + fundE - eigenE[l];
        temp += (QMatrixE[nu_mu][l]**2)/denom_e;
    return (temp.imag/-np.pi);

def gF_holes(omega, fundE, QMatrixH, eigenH, nu_mu):
    temp = 0;
    z = omega + eta*1j;

    #Hole part
    for l in range(len(eigenH)):
        denom_h = z - fundE + eigenH[l];
        temp += (QMatrixH[nu_mu][l]**2)/denom_h;
    return (temp.imag/-np.pi);

def greenFunction(omega, fundE, QMatrixE, QMatrixH, eigenE, eigenH, nu_mu):

    temp = gF_electrons(omega,fundE,QMatrixE,eigenE,nu_mu)\
            + gF_holes(omega,fundE,QMatrixH,eigenH,nu_mu)
    return temp;



#Reading data
evE = []
QMatrixE = [[]]
evH = [];
QMatrixH = [[]];
index = 0;

with open(fileNameQM) as file:
    reader = csv.reader(file,delimiter='\t')
    for lines in reader:
        #Nothing
        if (len(lines) == 0) : continue
        #Ignore les lignes de commentaire contenues dans le fichier 
        if(lines[0].startswith('#')):
            index += 1
            continue;
        else :
            if index == 0 :
                #Systems parameters
                sites = int(lines[0]);
                u = float(lines[1]);
                mu = float(lines[2]);
                N = int(lines[3]);
                S = int(lines[4]);
                fundE = float(lines[5]);
                eta = float(lines[6]);
                size = float(lines[7]);

                for k in range(sites-1): 
                    QMatrixE.append([]);
                    QMatrixH.append([]);
            elif index == 1 : 
                evE.append(float(lines[0]));
                for i in range(sites):
                    QMatrixE[i].append(complex(lines[i+1]))
            elif index == 2 : 
                evH.append(float(lines[0]));
                for i in range(sites):
                    QMatrixH[i].append(complex(lines[i+1]))

#Green points x-axis
nbrPoints = 2000;
if not set_omega :
    omegaMin = -3* np.pi;
    omegaMax = 3 * np.pi
else : 
    omegaMin = min_omega
    omegaMax = max_omega

nu_mu = np.arange(0,sites,1,dtype=int);

omega = np.linspace(omegaMin,omegaMax,nbrPoints,dtype = float, endpoint=True);

#Receiving data ref ###########################33
data = [[]];
corrected_fundE = [];
cov = [];
for i in range(len(nu_mu)-1): data.append([]);
if fitting : 
    data = [[]];
    lookingIndex = 0;
    ref_fundE = 0;

    for k in range(sites-1): data.append([]);
    with open(refFileName) as file:
        reader = csv.reader(file,delimiter='\t')
        for lines in reader:
            #Ignore les lignes de commentaire contenues dans le fichier 
            if(lines[0].startswith('#')):
                continue;
            elif(lines[0].startswith('NEXT_SITE')):
                lookingIndex += 1;

            elif(lines[0].startswith('PARAMETERS')):
                ref_fundE = float(lines[6]);
            else :
                data[lookingIndex].append(float(lines[1]) - lookingIndex);
    
    #Creating Green average ref
    averageData = np.zeros(len(data[0]));
    for data_set in data:
        for data_point in range(len(data_set)):
            averageData[data_point] += data_set[data_point];
    guess_fundE = ref_fundE#* (np.random.rand() * 0.1 + 1);

    #print(guess_fundE);
    #CURVE FIT ON EACH INDIVIDUAL CURVE
    #for i in range(sites) :
    #    popt, pcov = opt.curve_fit(lambda x, a: greenFunction(x, a, QMatrixE, QMatrixH, evE, evH, nu_mu[i]),omega,data[i], p0=guess_fundE, maxfev = 1000);
    #    corrected_fundE.append(popt[0]);
    #    cov.append(pcov[0][0]);
    #print(corrected_fundE);
    #for i in range(sites):
    #    if (np.sqrt(cov[i]) > 0.01 or ((corrected_fundE[i] - mean_fundE)*100/mean_fundE > 0.1)):
    #        print("WARNING MAY BE WRONG");
    #    #print(f"The standard deviation of the fundamental energy fitting of cruve {i} is {np.sqrt(cov[i])}. The deviation with the mean is {(corrected_fundE[i] - mean_fundE)*100/mean_fundE}%")
    #CURVE FIT ON THE AVERAGE CURVE (ACTUALLY ITS ONLY FITTED ON CURVE 0)
    popt, pcov = opt.curve_fit(lambda x, a: greenFunctionAverage(x, a, QMatrixE, QMatrixH, evE, evH, nu_mu,sites), omega, data[0], p0=guess_fundE, maxfev = 10000);
    corrected_fundE.append(popt[0]);
    cov.append(pcov[0][0]);

    mean_fundE = -33.90#np.mean(corrected_fundE);


    print(f" --- Green function fitting completed --- using mean fundE = {mean_fundE}");
############################################################

#Writting File
if fitting : fileName = "mod_greenFunctionValue"+".txt";
else : fileName = "greenFunctionValueQM"+".txt";
file = open(fileName,'w')
parametersGreenFileString = "PARAMETERS\t" + str(sites) + "\t" + str(u) + "\t" + str(mu) + "\t" + str(N) + "\t" + str(S) + "\t" + f"{fundE:.12f}";
#file.write("PARAMETERS\t4\t8\t4\t4\t0\n");
file.write(parametersGreenFileString);
file.write("\n#OMEGA\tVALUE\n");

for i in range(len(nu_mu)):
    for j in omega:
        ##mod_fundE = corrected_fundE[i] if fitting else fundE
        mod_fundE = mean_fundE if fitting else fundE
        tempe = gF_electrons(j, mod_fundE, QMatrixE, evE, nu_mu[i]);
        temph = gF_holes(j, mod_fundE, QMatrixH, evH, nu_mu[i]);
        lineString = str(j)+"\t"+str(tempe + nu_mu[i])+"\t"+str(temph + nu_mu[i])+"\n"; 
        file.write(lineString);
    file.write("NEXT_SITE\n");



print(f" --- Green function A_ii calculated and written in '{fileName}'")
