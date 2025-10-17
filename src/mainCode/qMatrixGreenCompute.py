#!/usr/bin/python3
import csv
import matplotlib.pyplot as plt
import numpy as np
import sys;

set_omega = False
file_name_qm = "qMatrices.txt";
if (len(sys.argv) == 4):
    # Custom bounds.
    min_omega = float(sys.argv[2])
    max_omega = float(sys.argv[3])
    set_omega = True 
    # File name given.
    file_name_qm = sys.argv[1]

elif (len(sys.argv) == 2):
    # File name given.
    file_name_qm = sys.argv[1];
elif len(sys.argv) > 2 :
    print("Number of given arguments incorrect")

def green_electrons(omega, fundE, q_matrix_e, eigen_e, nu_mu, eta):
    ######################################################
    # Computes the electron part of the green function at eneergy omega.
    #
    # Parameters
    # ----------
    # omega      : (float) energy value.
    # fundE      : (float) fundamental energy of the system.
    # q_matrix_e : (float-array) Q-matrix pre-computed with BL.
    # eigen_e    : (float-array) energies of the BL approximation.
    # nu_mu      : (int) which site to compute.
    # eta        : (float) eta value for Lorentzian.
    #
    # Returns
    # -------
    # g_omega_float : (float) green density weight
    ######################################################
    g_omega = 0;
    z = omega + eta*1j;

    #Electron part
    for l in range(len(eigen_e)):
        denom_e = z + fundE - eigen_e[l]
        g_omega += (q_matrix_e[nu_mu][l]**2) / denom_e

    g_omega_float = g_omega.imag / -np.pi
    return g_omega_float;

def green_holes(omega, fundE, q_matrix_h, eigen_h, nu_mu, eta):
    ######################################################
    # Computes the hole part of the green function at eneergy omega.
    #
    # Parameters
    # ----------
    # omega      : (float) energy value.
    # fundE      : (float) fundamental energy of the system.
    # q_matrix_h : (float-array) Q-matrix pre-computed with BL.
    # eigen_h    : (float-array) energies of the BL approximation.
    # nu_mu      : (int) which site to compute.
    # eta        : (float) eta value for Lorentzian.
    #
    # Returns
    # -------
    # g_omega_float : (float) green density weight
    ######################################################
    g_omega = 0;
    z = omega + eta*1j;

    #Hole part
    for l in range(len(eigen_h)):
        denom_h = z - fundE + eigen_h[l]
        g_omega += (q_matrix_h[nu_mu][l]**2) / denom_h

    g_omega_float = g_omega.imag / -np.pi
    return g_omega_float

def green_function(omega, fundE, q_matrix_e, q_matrix_h, eigen_e,
                   eigen_h, nu_mu, eta):
    ######################################################
    # Computes the green function at eneergy omega.
    #
    # Parameters
    # ----------
    # omega      : (float) energy value.
    # fundE      : (float) fundamental energy of the system.
    # q_matrix_e : (float-array) electron Q-matrix pre-computed with BL.
    # q_matrix_h : (float-array) hole Q-matrix pre-computed with BL.
    # eigen_e    : (float-array) electron energies of the BL approximation.
    # eigen_h    : (float-array) hole energies of the BL approximation.
    # nu_mu      : (int) which site to compute.
    # eta        : (float) eta value for Lorentzian.
    #
    # Returns
    # -------
    # g : (float) green density weight
    ######################################################

    g = green_electrons(omega, fundE, q_matrix_e, eigen_e, nu_mu, eta)\
        + green_holes(omega, fundE, q_matrix_h, eigen_h, nu_mu, eta)
    return g;


#Reading data
ev_e = []
qm_e = [[]]
ev_h = [];
qm_h= [[]];
index = 0;

with open(file_name_qm) as file:
    reader = csv.reader(file, delimiter='\t')
    for lines in reader:
        #Nothing
        if (len(lines) == 0) : continue
        #Ignore les lignes de commentaire contenues dans le fichier 
        if (lines[0].startswith('#')):
            index += 1
            continue;
        else :
            if index == 0 :
                #Systems parameters
                sites = int(lines[0])
                u = float(lines[1])
                mu = float(lines[2])
                N = int(lines[3])
                S = int(lines[4])
                fundE = float(lines[5])
                eta = float(lines[6])
                size = float(lines[7])

                for k in range(sites-1): 
                    qm_e.append([])
                    qm_h.append([])

            elif index == 1 : 
                ev_e.append(float(lines[0]))
                for i in range(sites):
                    qm_e[i].append(complex(lines[i+1]))
            elif index == 2 : 
                ev_h.append(float(lines[0]));
                for i in range(sites):
                    qm_h[i].append(complex(lines[i+1]))

#Green points x-axis
nbr_points = 2000;
if not set_omega :
    omega_min = -3* np.pi;
    omega_max = 3 * np.pi
else : 
    omega_min = min_omega
    omega_max = max_omega

nu_mu = np.arange(0, sites, 1, dtype=int);

omega = np.linspace(omega_min, omega_max, nbr_points, 
                    dtype=float, endpoint=True);

#Writting File
file_name = "greenFunctionValueQM"+".txt";
file = open(file_name, 'w')
parametersGreenFileString = ("PARAMETERS\t" + str(sites) + "\t" + str(u) 
                          + "\t" + str(mu) + "\t" + str(N) + "\t" + str(S) 
                          + "\t" + f"{fundE:.12f}")
#file.write("PARAMETERS\t4\t8\t4\t4\t0\n");
file.write(parametersGreenFileString);
file.write("\n#OMEGA\tVALUE\n");

for i in range(len(nu_mu)):
    for j in omega:
        mod_fundE = fundE
        tempe = green_electrons(j, mod_fundE, qm_e, ev_e, nu_mu[i], eta);
        temph = green_holes(j, mod_fundE, qm_h, ev_h, nu_mu[i], eta);
        line_string = (str(j) + "\t" + str(tempe + nu_mu[i]) + "\t" 
                    + str(temph + nu_mu[i])+"\n")
        file.write(line_string);
    file.write("NEXT_SITE\n");



print(f" --- Green function A_ii calculated and written in '{file_name}'")
