#!/usr/bin/env python
# coding: utf-8
# In[30]:
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import os
import re
# import pandas as pd
# import sklearn
# from sklearn.linear_model import LinearRegression
from scipy.integrate import odeint

def HSE_energy(structure_path):
    # with open(path+'/HSE/OUTCAR') as f:
    #     content = f.readlines()
    
    os.chdir(structure_path + "/HSE/")
    os.system("(grep sigma OUTCAR | tail -n 1) > HSEenergy.txt")
    with open("HSEenergy.txt", 'r') as f:
        HSE = f.readline()
        # print(HSE)
        HSE = float(HSE[(len(HSE)-16):].strip())
    
    
    # energy = subprocess.run(["grep sigma OUTCAR | tail -n 1", "OUTCAR"], check=True, capture_output=True, text=True).stdout
    
    # for line in content:
    #     if ('energy' in line) and ('without' in line) and ('entropy' in line) and ('sigma' in line):
    #         energy_zero = float(line.split()[-1])
    return HSE
    
def PBE_energy(structure_path):
    # with open(path+'/PBE/OUTCAR') as f:
    # content = f.readlines()
    
    os.chdir(structure_path + "/PBE/")
    os.system("(grep sigma OUTCAR | tail -n 1) > PBEenergy.txt")
    with open("PBEenergy.txt", 'r') as f:
        PBE = f.readline()
        # print(PBE)
        PBE = float(PBE[(len(PBE)-16):].strip())
        
        
    # energy = subprocess.run(["grep sigma OUTCAR | tail -n 1", "OUTCAR"], check=True, capture_output=True, text=True).stdout
    
    # for line in content: 
    #       if('energy' in line) and ('without' in line) and ('entropy' in line) and ('sigma' in line):
    #           energy_zero = float(line.split()[-1])
    return PBE
    
def shomate(t,molecule):
    ''' shomate
            Computes the standard entropy of methane based on temperature using Shomate equation    
        Parameters:
            t = temperature (K)
            molecule = molecule of choice (methane, oxygen, hydrogen, carbon monoxide, carbon dioxide, water, formaldehyde, methanol, etc)
    '''
    
    shomate_params = {"A": 0, "B": 0, "C": 0, "D": 0, "E": 0, "F": 0, "G": 0, "H": 0}
    if molecule == ("methane"):
        if t in range(273,1300):
            shomate_params = {"A": -0.703029, "B": 108.4773, "C": -42.52157, "D": 5.862788, "E": 0.678565, "F": -76.84376, "G": 158.7163, "H": -74.87310}
        elif t in range(1300, 6000):
            shomate_params = {"A": 85.81217, "B": 11.26467, "C": -2.114146, "D": 0.138190, "E": -26.42221, "F": -153.5327, "G": 224.4143, "H": -74.87310}
    elif molecule == ("oxygen"):
        if t in range(100,700):
            shomate_params = {"A": 31.32234, "B": -20.23531, "C": 57.86644, "D": -36.50624, "E": -0.007374, "F": -8.903471, "G": 246.7945, "H": 0}
        elif t in range(700, 2000):
            shomate_params = {"A": 30.03235, "B": 8.772972, "C": -3.988133, "D": 0.788313, "E": -0.741599, "F": -11.32468, "G": 236.1663, "H": 0}
        elif t in range(2000, 6000):
            shomate_params = {"A": 20.91111, "B": 10.72071, "C": -2.020498, "D": 0.146449, "E": 9.245722, "F": 5.337651, "G": 237.6185, "H": 0}
    elif molecule == ("water"):
        if t in range(500,1700):
            shomate_params = {"A": 30.092, "B": 6.832514, "C": 6.793435, "D": -2.53448, "E": 0.082139, "F": -250.881, "G": 223.3967, "H": -241.8264}
    elif molecule == ("CO2"):
        if t in range(298,1200):
            shomate_params = {"A": 24.99735, "B": 55.18696, "C": -33.69137, "D": 7.948387, "E": -0.136638, "F": -403.6075, "G": 228.2431, "H": -393.5224}
    elif molecule == ("CO"):
        if t in range(298,1300):
            shomate_params = {"A": 25.56759, "B": 6.09613, "C": 4.054656, "D": -2.671301, "E": 0.131021, "F": -118.0089, "G": 227.3665, "H": -110.5271}
    elif molecule == ("formaldehyde"):
        if t in range(298,1200):
            shomate_params = {"A": 5.193767, "B": 93.23249, "C": -44.85457, "D": 7.882279, "E": 0.551175, "F": -119.3591, "G": 202.4663, "H": -115.8972}
    else:
        return 0
    
    t = t/1000
    S = shomate_params["A"]*np.log(t) + shomate_params["B"]*t + shomate_params["C"]*(t*t)/2 + shomate_params["D"]*(t*t*t)/3 - shomate_params["E"]/(2*(t*t)) + shomate_params["G"] #J/(mol*K). Note that log() is natural log in python
    
    # print("J/mol K:")
    # print(S)
    #conversion to eV/K
    S = S*6.242e18/6.022e23
    # S = S/96.487/1000
    # print("eV/K:")
    # print(S)
    
    return S

def freq(structure_path):
    os.chdir(structure_path + "/vtst/")
    os.system('/projects/academic/mdupuis2/software/vtst/vtsttools/vtstscripts/dymmatrix.pl')
    
    freq = []
    with open('freq.dat') as f:
        content = f.readlines()
        for line in content:
            line = line.strip()
            nums = re.findall(r"[-+]?\d*\.\d+|\d+", line)
            freq.append(float(nums[0]))
    freq = np.array(freq)
    l = len(freq)
    for k in range(l):
        if freq[k] < 100:
            freq[k] = 100
    freq = freq * c
    zpe = .5*h*freq
    ZPE = np.sum(zpe)
    qvib = 1/(1-np.exp((-h*freq)/(kB*T))) #check if this should be freq or omega
    Qvib = np.prod(qvib)
    # Fvib = ZPE + kB*T*np.log(1/Qvib)
    # print(Fvib)
    Svib = kB*np.log(Qvib)
   
    return [ZPE, Svib]
    

kB = 8.617E-5 # eV/K
h = 4.136E-15 # eV*s
c = 2.998E10 # cm/s
T = 298.15 # K

project_path = "/projects/academic/ericwalk/ruthbell/Ads_Ag-SSZ13/"

# Dictionary of structures. The value is a list to contain HSE, PBE, ZPE, Svib, and G (free energy), in that order
structures = {"Blank_Ag-SSZ-13":[],
              "C2H4_1":[],
              "CO2_1_Alternate":[],
              "H2O_1":[],
              "O_1":[],
              "CO2_gas":[],
              "O2_gas":[],
              "COS-C":[],
              "COS-S":[],
              "COS_gas":[],
              "SO2-O":[],
              "SO2-S":[],
              "SO2_gas":[],
              "O2_1":[],
              "DMS":[],
              "DMS_gas":[],
              "H2S-S":[],
              "H2S_gas":[],
            }
# Dictionary of gases. The value is a list to contain HSE, PBE, ZPE, S_shomate, and G (free energy), in that order
gases = {"O2_gas":[],
         "CO2_gas":[]
        }

# for s in structures.keys():
#     structures[s]

# quit()



# # Structure Energies

# # Blank_Ag-SSZ-13 - GOOD - ASK IF WE SHOULD DO FREQ JOB FOR THE BLANK STRUCTURE
# structures["Blank_Ag-SSZ-13"].append(HSE_energy(project_path + "Blank_Ag-SSZ-13")) # Blank_Ag-SSZ-13 HSE energy
# structures["Blank_Ag-SSZ-13"].append(PBE_energy(project_path + "Blank_Ag-SSZ-13")) # Blank_Ag-SSZ-13 PBE energy
# # structures["Blank_Ag-SSZ-13"].append(freq(project_path + "Blank_Ag-SSZ-13")[0]) # Blank_Ag-SSZ-13 ZPE energy
# # structures["Blank_Ag-SSZ-13"].append(freq(project_path + "Blank_Ag-SSZ-13")[1]) # Blank_Ag-SSZ-13 entropy
# # structures["Blank_Ag-SSZ-13"].append(structures["Blank_Ag-SSZ-13"][0] + structures["Blank_Ag-SSZ-13"][1] - structures["Blank_Ag-SSZ-13"][2]*T) # Blank_Ag-SSZ-13 free energy = HSE + ZPE + S*T

# print("Blank_Ag-SSZ-13: HSE, PBE") # , ZPE, S, G") 
# print(structures["Blank_Ag-SSZ-13"])
# print()

# # #C2H4_1 - GOOD
# # structures["C2H4_1"].append(HSE_energy(project_path + "C2H4_1")) # C2H4_1 HSE energy 
# # structures["C2H4_1"].append(PBE_energy(project_path + "C2H4_1")) # C2H4_1 PBE energy
# # structures["C2H4_1"].append(freq(project_path + "C2H4_1")[0]) # C2H4_1 ZPE energy
# # structures["C2H4_1"].append(freq(project_path + "C2H4_1")[1]) # C2H4_1 entropy
# # structures["C2H4_1"].append(structures["C2H4_1"][0] + structures["C2H4_1"][2] - structures["C2H4_1"][3]*T) # C2H4_1 free energy = HSE + ZPE + S*T
# # structures["C2H4_1"].append(structures["C2H4_1"][1] + structures["C2H4_1"][2] - structures["C2H4_1"][3]*T) # C2H4_1 free energy = PBE + ZPE + S*T

# # print("C2H4_1: HSE, PBE, ZPE, S, HSEG, PBEG")
# # print(structures["C2H4_1"])
# # print()

# #CO2_1_Alternate - GOOD
# structures["CO2_1_Alternate"].append(HSE_energy(project_path + "CO2_1_Alternate")) # CO2_1_Alternate HSE energy
# structures["CO2_1_Alternate"].append(PBE_energy(project_path + "CO2_1_Alternate")) # CO2_1_Alternate PBE energy
# structures["CO2_1_Alternate"].append(freq(project_path + "CO2_1_Alternate")[0]) # CO2_1_Alternate ZPE energy
# structures["CO2_1_Alternate"].append(freq(project_path + "CO2_1_Alternate")[1]) # CO2_1_Alternate entropy
# structures["CO2_1_Alternate"].append(structures["CO2_1_Alternate"][0] + structures["CO2_1_Alternate"][2] - structures["CO2_1_Alternate"][3]*T) # CO2_1_Alternate free energy = HSE + ZPE + S*T
# structures["CO2_1_Alternate"].append(structures["CO2_1_Alternate"][1] + structures["CO2_1_Alternate"][2] - structures["CO2_1_Alternate"][3]*T) # CO2_1_Alternate free energy = PBE + ZPE + S*T

# print("CO2_1_Alternate: HSE, PBE, ZPE, S, HSEG, PBEG")
# print(structures["CO2_1_Alternate"])
# print()

# # CO2_gas
# structures["CO2_gas"].append(HSE_energy(project_path + "CO2_gas")) # CO2 gas HSE energy
# structures["CO2_gas"].append(PBE_energy(project_path + "CO2_gas")) # CO2_gas PBE energy
# structures["CO2_gas"].append(freq(project_path + "CO2_gas")[0]) # CO2 gas ZPE energy
# structures["CO2_gas"].append(freq(project_path + "CO2_gas")[1]) # CO2 gas S_shomate
# structures["CO2_gas"].append(structures["CO2_gas"][0] + structures["CO2_gas"][2] - structures["CO2_gas"][3]*T) # CO2 gas free energy = HSE + ZPE + S*T
# structures["CO2_gas"].append(structures["CO2_gas"][1] + structures["CO2_gas"][2] - structures["CO2_gas"][3]*T) # CO2_gas free energy = PBE + ZPE + S*T

# print("CO2_gas: HSE, PBE, ZPE, S, HSEG, PBEG")
# print(structures["CO2_gas"])
# print()

# # #H2O_1 - GOOD
# # structures["H2O_1"].append(HSE_energy(project_path + "H2O_1")) # H2O_1 HSE energy
# # structures["H2O_1"].append(PBE_energy(project_path + "H2O_1")) # H2O_1 PBE energy
# # structures["H2O_1"].append(freq(project_path + "H2O_1")[0]) # H2O_1 ZPE energy
# # structures["H2O_1"].append(freq(project_path + "H2O_1")[1]) # H2O_1 entropy
# # structures["H2O_1"].append(structures["H2O_1"][0] + structures["H2O_1"][2] - structures["H2O_1"][3]*T) # H2O_1 free energy = HSE + ZPE + S*T
# # structures["H2O_1"].append(structures["H2O_1"][1] + structures["H2O_1"][2] - structures["H2O_1"][3]*T) # H2O_1 free energy = PBE + ZPE + S*T

# print("H2O_1: HSE, PBE, ZPE, S, HSEG, PBEG")
# print(structures["H2O_1"])
# print()

# #O_1 - GOOD
# structures["O_1"].append(HSE_energy(project_path + "O_1")) # O_1 HSE energy
# structures["O_1"].append(PBE_energy(project_path + "O_1")) # O_1 PBE energy
# structures["O_1"].append(freq(project_path + "O_1")[0]) # O_1 ZPE energy
# structures["O_1"].append(freq(project_path + "O_1")[1]) # O_1 entropy
# structures["O_1"].append(structures["O_1"][0] + structures["O_1"][2] - structures["O_1"][3]*T) # O_1 free energy = HSE + ZPE + S*T
# structures["O_1"].append(structures["O_1"][1] + structures["O_1"][2] - structures["O_1"][3]*T) # O_1 free energy = PBE + ZPE + S*T

# print("O_1: HSE, PBE, ZPE, S, HSEG, PBEG")
# print(structures["O_1"])
# print()


# #O2_gas
# structures["O2_gas"].append(HSE_energy(project_path + "O2_gas")) # O2 gas HSE energy
# structures["O2_gas"].append(PBE_energy(project_path + "O2_gas")) # O2_gas PBE energy
# structures["O2_gas"].append(freq(project_path + "O2_gas")[0]) # O2 gas ZPE energy
# structures["O2_gas"].append(freq(project_path + "O2_gas")[1]) # O2 gas S_shomate
# structures["O2_gas"].append(structures["O2_gas"][0] + structures["O2_gas"][2] - structures["O2_gas"][3]*T) # O2 gas free energy = HSE + ZPE + S*T
# structures["O2_gas"].append(structures["O2_gas"][1] + structures["O2_gas"][2] - structures["O2_gas"][3]*T) # O2_gas free energy = PBE + ZPE + S*T

# print("O2_gas: HSE, PBE, ZPE, S, HSEG, PBEG")
# print(structures["O2_gas"])
# print()

# #O2_1
# structures["O2_1"].append(HSE_energy(project_path + "O2_1")) # O2_1 HSE energy
# structures["O2_1"].append(PBE_energy(project_path + "O2_1")) # O2_1 PBE energy
# structures["O2_1"].append(freq(project_path + "O2_1")[0]) # O2_1 ZPE energy
# structures["O2_1"].append(freq(project_path + "O2_1")[1]) # O2_1 entropy
# structures["O2_1"].append(structures["O2_1"][0] + structures["O2_1"][2] - structures["O2_1"][3]*T) # O2_1 free energy = HSE + ZPE + S*T
# structures["O2_1"].append(structures["O2_1"][1] + structures["O2_1"][2] - structures["O2_1"][3]*T) # O2_1 free energy = PBE + ZPE + S*T

# print("O2_1: HSE, PBE, ZPE, S, HSEG, PBEG")
# print(structures["O2_1"])
# print()

# #SO2-O
# structures["SO2-O"].append(HSE_energy(project_path + "SO2-O")) # SO2-O HSE energy
# structures["SO2-O"].append(PBE_energy(project_path + "SO2-O")) # SO2-O PBE energy
# structures["SO2-O"].append(freq(project_path + "SO2-O")[0]) # SO2-O ZPE energy
# structures["SO2-O"].append(freq(project_path + "SO2-O")[1]) # SO2-O entropy
# structures["SO2-O"].append(structures["SO2-O"][0] + structures["SO2-O"][2] - structures["SO2-O"][3]*T) # SO2-O free energy = HSE + ZPE + S*T
# structures["SO2-O"].append(structures["SO2-O"][1] + structures["SO2-O"][2] - structures["SO2-O"][3]*T) #SO2-O free energy = PBE + ZPE + S*T

# print("SO2-O: HSE, PBE, ZPE, S, HSEG, PBEG")
# print(structures["SO2-O"])
# print()

# #SO2-S
# structures["SO2-S"].append(HSE_energy(project_path + "SO2-S")) # SO2-S HSE energy
# structures["SO2-S"].append(PBE_energy(project_path + "SO2-S")) # SO2-S PBE energy
# structures["SO2-S"].append(freq(project_path + "SO2-S")[0]) # SO2-S ZPE energy
# structures["SO2-S"].append(freq(project_path + "SO2-S")[1]) # SO2-S entropy
# structures["SO2-S"].append(structures["SO2-S"][0] + structures["SO2-S"][2] - structures["SO2-S"][3]*T) # SO2-S free energy = HSE + ZPE + S*T
# structures["SO2-S"].append(structures["SO2-S"][1] + structures["SO2-S"][2] - structures["SO2-S"][3]*T) #SO2-S free energy = PBE + ZPE + S*T

# print("SO2-S: HSE, PBE, ZPE, S, HSEG, PBEG")
# print(structures["SO2-S"])
# print()

# #SO2_gas
# structures["SO2_gas"].append(HSE_energy(project_path + "SO2_gas")) # SO2 gas HSE energy
# structures["SO2_gas"].append(PBE_energy(project_path + "SO2_gas")) # SO2_gas PBE energy
# structures["SO2_gas"].append(freq(project_path + "SO2_gas")[0]) # SO2 gas ZPE energy
# structures["SO2_gas"].append(freq(project_path + "SO2_gas")[1]) # SO2 gas S_shomate
# structures["SO2_gas"].append(structures["SO2_gas"][0] + structures["SO2_gas"][2] - structures["SO2_gas"][3]*T) # SO2 gas free energy = HSE + ZPE + S*T
# structures["SO2_gas"].append(structures["SO2_gas"][1] + structures["SO2_gas"][2] - structures["SO2_gas"][3]*T) # SO2_gas free energy = PBE + ZPE + S*T

# print("SO2_gas: HSE, PBE, ZPE, S, HSEG, PBEG")
# print(structures["SO2_gas"])
# print()

# #COS-C
# structures["COS-C"].append(HSE_energy(project_path + "COS-C")) # COS-C HSE energy
# structures["COS-C"].append(PBE_energy(project_path + "COS-C")) # COS-C PBE energy
# structures["COS-C"].append(freq(project_path + "COS-C")[0]) # COS-C ZPE energy
# structures["COS-C"].append(freq(project_path + "COS-C")[1]) # COS-C entropy
# structures["COS-C"].append(structures["COS-C"][0] + structures["COS-C"][2] - structures["COS-C"][3]*T) # COS-C free energy = HSE + ZPE + S*T
# structures["COS-C"].append(structures["COS-C"][1] + structures["COS-C"][2] - structures["COS-C"][3]*T) #COS-C free energy = PBE + ZPE + S*T

# print("COS-C: HSE, PBE, ZPE, S, HSEG, PBEG")
# print(structures["COS-C"])
# print()

# #COS-S
# structures["COS-S"].append(HSE_energy(project_path + "COS-S")) # COS-S HSE energy
# structures["COS-S"].append(PBE_energy(project_path + "COS-S")) # COS-S PBE energy
# structures["COS-S"].append(freq(project_path + "COS-S")[0]) # COS-S ZPE energy
# structures["COS-S"].append(freq(project_path + "COS-S")[1]) # COS-S entropy
# structures["COS-S"].append(structures["COS-S"][0] + structures["COS-S"][2] - structures["COS-S"][3]*T) # COS-S free energy = HSE + ZPE + S*T
# structures["COS-S"].append(structures["COS-S"][1] + structures["COS-S"][2] - structures["COS-S"][3]*T) # COS-S free energy = PBE + ZPE + S*T

# print("COS-S: HSE, PBE, ZPE, S, HSEG, PBEG")
# print(structures["COS-S"])
# print()

# #COS_gas
# structures["COS_gas"].append(HSE_energy(project_path + "COS_gas")) # COS_gas HSE energy
# structures["COS_gas"].append(PBE_energy(project_path + "COS_gas")) # COS_gas PBE energy
# structures["COS_gas"].append(freq(project_path + "COS_gas")[0]) # COS_gas ZPE energy
# structures["COS_gas"].append(freq(project_path + "COS_gas")[1]) # COS_gas S_shomate
# structures["COS_gas"].append(structures["COS_gas"][0] + structures["COS_gas"][2] - structures["COS_gas"][3]*T) # COS_gas free energy = HSE + ZPE + S*T
# structures["COS_gas"].append(structures["COS_gas"][1] + structures["COS_gas"][2] - structures["COS_gas"][3]*T) # COS_gas free energy = PBE + ZPE + S*T

# print("COS_gas: HSE, PBE, ZPE, S, HSEG, PBEG")
# print(structures["COS_gas"])
# print()

# #DMS
# structures["DMS"].append(HSE_energy(project_path + "DMS")) # DMS HSE energy
# structures["DMS"].append(PBE_energy(project_path + "DMS")) # DMS PBE energy
# structures["DMS"].append(freq(project_path + "DMS")[0]) # DMS ZPE energy
# structures["DMS"].append(freq(project_path + "DMS")[1]) # DMS entropy
# structures["DMS"].append(structures["DMS"][0] + structures["DMS"][2] - structures["DMS"][3]*T) # DMS free energy = HSE + ZPE + S*T
# structures["DMS"].append(structures["DMS"][1] + structures["DMS"][2] - structures["DMS"][3]*T) # DMS free energy = PBE + ZPE + S*T

# print("DMS: HSE, PBE, ZPE, S, HSEG, PBEG")
# print(structures["DMS"])
# print()

# #DMS_gas
# structures["DMS_gas"].append(HSE_energy(project_path + "DMS_gas")) # DMS_gas HSE energy
# structures["DMS_gas"].append(PBE_energy(project_path + "DMS_gas")) # DMS_gas PBE energy
# structures["DMS_gas"].append(freq(project_path + "DMS_gas")[0]) # DMS_gas ZPE energy
# structures["DMS_gas"].append(freq(project_path + "DMS_gas")[1]) # DMS_gas entropy
# structures["DMS_gas"].append(structures["DMS_gas"][0] + structures["DMS_gas"][2] - structures["DMS_gas"][3]*T) # DMS_gas free energy = HSE + ZPE + S*T
# structures["DMS_gas"].append(structures["DMS_gas"][1] + structures["DMS_gas"][2] - structures["DMS_gas"][3]*T) # DMS_gas free energy = PBE + ZPE + S*T

# print("DMS_gas: HSE, PBE, ZPE, S, HSEG, PBEG")
# print(structures["DMS_gas"])
# print()

# #H2S-S
# structures["H2S-S"].append(HSE_energy(project_path + "H2S-S")) # H2S-S HSE energy
# structures["H2S-S"].append(PBE_energy(project_path + "H2S-S")) # H2S-S PBE energy
# structures["H2S-S"].append(freq(project_path + "H2S-S")[0]) # H2S-S ZPE energy
# structures["H2S-S"].append(freq(project_path + "H2S-S")[1]) # H2S-S entropy
# structures["H2S-S"].append(structures["H2S-S"][0] + structures["H2S-S"][2] - structures["H2S-S"][3]*T) # H2S-S free energy = HSE + ZPE + S*T
# structures["H2S-S"].append(structures["H2S-S"][1] + structures["H2S-S"][2] - structures["H2S-S"][3]*T) # H2S-S free energy = PBE + ZPE + S*T

# print("H2S-S: HSE, PBE, ZPE, S, HSEG, PBEG")
# print(structures["H2S-S"])
# print()

# #H2S_gas  
# structures["H2S_gas"].append(HSE_energy(project_path + "H2S_gas")) # H2S_gas HSE energy
# structures["H2S_gas"].append(PBE_energy(project_path + "H2S_gas")) # H2S_gas PBE energy
# structures["H2S_gas"].append(freq(project_path + "H2S_gas")[0]) # H2S_gas ZPE energy
# structures["H2S_gas"].append(freq(project_path + "H2S_gas")[1]) # H2S_gas entropy
# structures["H2S_gas"].append(structures["H2S_gas"][0] + structures["H2S_gas"][2] - structures["H2S_gas"][3]*T) # H2S_gas free energy = HSE + ZPE + S*T
# structures["H2S_gas"].append(structures["H2S_gas"][1] + structures["H2S_gas"][2] - structures["H2S_gas"][3]*T) # H2S_gas free energy = PBE + ZPE + S*T

# print("H2S_gas: HSE, PBE, ZPE, S, HSEG, PBEG")
# print(structures["H2S_gas"])
# print()
#-----------------------------------------

# quit()

#---------------------lateral interaction microkinetic model

def microkinetic_model(y, t, CO2_1_Alternate_gas_frac, O2_gas_frac, SO2_gas_frac, DMS_gas_frac, COS_gas_frac, H2S_gas_frac, G_1, G_2, G_3, G_4, G_5, G_6, f1, f2, f3, f4, f5, f6): 
    theta_CO2_1_Alternate, theta_O2, theta_SO2, theta_DMS, theta_COS, theta_H2S, theta_vacant = y
    #Equilibrium Constants
    K_1 = np.exp(-(G_1)/(kB*T))
    K_2 = np.exp(-(G_2)/(kB*T))
    K_3 = np.exp(-(G_3)/(kB*T))
    K_4 = np.exp(-(G_4)/(kB*T))
    K_5 = np.exp(-(G_5)/(kB*T))
    K_6 = np.exp(-(G_6)/(kB*T))
    
    #print(G_1,G_2,G_3,G_4,G_5,G_6) 
    #Reverse Rate Constants
    r1 = f1/K_1
    r2 = f2/K_2
    r3 = f3/K_3
    r4 = f4/K_4
    r5 = f5/K_5
    r6 = f6/K_6
    #print(r1,r2,r3,r4,r5,r6)
    #print(CO2_1_Alternate_gas_frac, O2_gas_frac, SO2_gas_frac, DMS_gas_frac, COS_gas_frac, H2S_gas_frac)
    #Reaction Rates
    r_rxn1 = f1*theta_vacant*CO2_1_Alternate_gas_frac - r1*theta_CO2_1_Alternate
    r_rxn2 = f2*theta_vacant*O2_gas_frac - r2*theta_O2
    r_rxn3 = f3*theta_vacant*SO2_gas_frac - r3*theta_SO2
    r_rxn4 = f4*theta_vacant*DMS_gas_frac - r4*theta_DMS
    r_rxn5 = f5*theta_vacant*COS_gas_frac - r5*theta_COS
    r_rxn6 = f6*theta_vacant*H2S_gas_frac - r6*theta_H2S
    #Microkinetic Model Equations (Mass Balances)
    d_theta_CO2_1_Alternate_dt = r_rxn1 
    d_theta_O2_dt = r_rxn2
    d_theta_SO2_dt = r_rxn3
    d_theta_DMS_dt = r_rxn4
    d_theta_COS_dt = r_rxn5
    d_theta_H2S_dt = r_rxn6
    d_theta_vacant_dt = -r_rxn1 - r_rxn2 - r_rxn3 - r_rxn4 - r_rxn5 - r_rxn6 
    dydt = [d_theta_CO2_1_Alternate_dt, d_theta_O2_dt, d_theta_SO2_dt, d_theta_DMS_dt, d_theta_COS_dt, d_theta_H2S_dt, d_theta_vacant_dt]
    return dydt
    
# CO2_Alternate_fractions = np.zeros(len(M_1))
# H2O_fractions = np.zeros(len(M_1))
y0 = y=np.array([0, 0, 0, 0, 0, 0, 1])
t = np.linspace(0, 1E10, int(1E4)) #1E10,
#for i in range(len(M_1)):  #M_1 holds all slopes for C2H4.
f1 = 6.2014895E5 #adsorption rate constant ethylene on zeolite                    # from collision theory
f2 = 7.7382354E5 #adsorption rate water from collision theory
f3 = 7.7382354E5
f4 = 7.7382354E5
f5 = 7.7382354E5
f6 = 7.7382354E5
G_1 = -0.4510 #CO2
G_2 = -0.1213 #O2
G_3 = -0.57039 #SO2
G_4 = -0.8319 #DMS
G_5 = -0.7561 #COS
G_6 = -1.1473 #H2S
theta_individual = odeint(microkinetic_model, y0, t, args=(0.75, 0.05, 0.05, 0.05, 0.05, 0.05, G_1, G_2,  G_3, G_4, G_5, G_6, f1, f2, f3, f4, f5, f6)) 
#print(str(theta_individual[-1,:]))
CO2_1_Alternate_fractions = theta_individual[-1,0]
O2_fractions = theta_individual[-1,1]
SO2_fractions = theta_individual[-1,2]
DMS_fractions = theta_individual[-1,3]
COS_fractions = theta_individual[-1,4]
H2S_fractions = theta_individual[-1,5]
vacant_fractions = theta_individual[-1,6]
adsorption_fractions = np.column_stack((CO2_1_Alternate_fractions,O2_fractions,SO2_fractions,DMS_fractions,COS_fractions,H2S_fractions))
#print('adsadsorption_fractions.shape: ' + str(adsorption_fractions.shape))

CO2_1_Alternate_fractions = CO2_1_Alternate_fractions
O2_fractions = O2_fractions
vacant_fractions = vacant_fractions

print("CO2_1_Alternate_fractions",CO2_1_Alternate_fractions)
print("O2_fractions",O2_fractions)
print("SO2_fractions",SO2_fractions)
print("DMS_fractions",DMS_fractions)
print("COS_fractions",COS_fractions)
print("H2S_fractions",H2S_fractions)
print("vacant_fractions",vacant_fractions)

theta_CO2_1_Alternate = CO2_1_Alternate_fractions

#-----------------------------------------

#Plot
preset_bins = np.linspace(0,1,10)
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)

fig, ax = plt.subplots()
ax.hist([CO2_1_Alternate_fractions, O2_fractions], bins=preset_bins, label=[r'$\theta_{CO2_1_Alternate}$',r'$\theta_{O2}$'])
ax.set_title('Microkinetic Model Adsorption Fractions',  fontsize=18)
ax.set_ylabel('Count', fontsize=18)
ax.set_xlabel('Surface Fraction', fontsize=18)
ax.set_xlim(-.02,1.02)
ax.set_ylim(0,2000)
plt.legend(loc='upper right', fontsize=16)
#fig.savefig('ThetaMicMod.png', dpi=220)

# barwidth=2
# fig = plt.figure(figsize=(5,4.5))  #5,4.5
# ax = fig.add_subplot(111)
# plt.bar(theta_CO2_1_Alternate, color='b', width=barwidth, label=r'$\theta_{CO2_1_Alternate}$') 
# plt.bar(theta_O2, bottom=theta_C2H4, color='g', width=barwidth, label=r'$\theta_{O2}$')
# plt.bar(theta_vacant, color='k', width=barwidth, label=r'$\theta_{vacant}$')
# ax.set_xlabel('Count', fontsize=12)
# ax.set_ylabel('Adsorption Site Coverage (%)', fontsize=12)
# ax.legend(loc='lower left', bbox_to_anchor=(-0.01,1.0,1,0.2), ncol=6, fontsize='large')
# fig.savefig('AdFrac_barchart_MKM.png', dpi=220)

quit() # Caitlin's old code is below
# Gas phase energies inputted manually because we have the values from previous projects

# # O2
# gases["O2_gas"].append(-17.00568799) # O2 gas HSE energy
# gases["O2_gas"].append(0.1280465295225186) # O2 gas ZPE energy
# gases["O2_gas"].append(shomate(T,"oxygen")) # O2 gas S_shomate
# gases["O2_gas"].append(gases["O2_gas"][0] + gases["O2_gas"][1] - gases["O2_gas"][2]*T) # O2 gas free energy = HSE + ZPE + S*T

# print("O2_gas: HSE, ZPE, S, G")
# print(gases["O2_gas"])
# print()

# # CO2
# gases["CO2_gas"].append(-32.16876137) # CO2 gas HSE energy
# gases["CO2_gas"].append(0.3378522930999075) # CO2 gas ZPE energy
# gases["CO2_gas"].append(shomate(T,"CO2")) # CO2 gas S_shomate
# gases["CO2_gas"].append(gases["CO2_gas"][0] + gases["CO2_gas"][1] - gases["CO2_gas"][2]*T) # CO2 gas free energy = HSE + ZPE + S*T

# print("CO2_gas: HSE, ZPE, S, G")
# print(gases["CO2_gas"])
# print()
# 
#########get slopes

#C2H4
C2H4_1_egy = np.genfromtxt('BEEF_vdw_C2H4-1.csv') + 1.43 - 1.34
#C2H4_2_egy = np.genfromtxt('BEEF_vdw_C2H4-2.csv') + 2*(1.43 - 1.34)
#C2H4_3_egy = np.genfromtxt('BEEF_vdw_C2H4-3.csv') + 3*(1.43 - 1.34)
C2H4_2_egy = np.genfromtxt('BEEF_vdw_C2H4-2.csv') + 0.528 - 0.301
C2H4_3_egy = np.genfromtxt('BEEF_vdw_C2H4-3.csv') + 0.279 - (-0.00861)
C2H4_free_np = np.genfromtxt('BEEF_vdw_ensemble_energies_C2H4_free.csv') + 1.3886 - 0.67770
zeolite_np = np.genfromtxt('BEEF_vdw_Zeolite.csv')



#y_energies is delta G
y_energies = np.zeros((len(C2H4_1_egy),3))
y_energies[:,0] = C2H4_1_egy - C2H4_free_np - zeolite_np
y_energies[:,1] = C2H4_2_egy - C2H4_1_egy - C2H4_free_np 
y_energies[:,2] = C2H4_3_egy - C2H4_2_egy - C2H4_free_np
X = np.array([0,1/3,2/3])
X = X.reshape((-1,1))
M_1 = np.zeros(len(C2H4_1_egy))
for i in range(len(C2H4_1_egy)): #len(C2H4_1_egy) 
  y = y_energies[i,:]
  y = y-y[0] # has to be through origin
  	#y = y.reshape((1,-1))
  	#reg = LinearRegression().fit(X, y)
  	#m = reg.coef_  # will be one entry in delta_G_m_2
  m1, _, _, _ = np.linalg.lstsq(X, y)
  M_1[i] = m1

C2H4_eng = np.concatenate(y_energies).reshape((2000,3))

#------------------------------------------------------------------------
#PairPlot
#C2H4pd_eng = pd.DataFrame(np.array(C2H4_eng), columns=['$G_{rxn1}$ (eV)', '$G_{rxn2}$ (eV)', '$G_{rxn3}$ (eV)'])
#C2H4_PairPlot = sns.pairplot(C2H4pd_eng)
#C2H4_PairPlot.set(xlim=(-3,2))
#C2H4_PairPlot.set(ylim=(-3,2))
#C2H4_PairPlot.fig.suptitle("C2H4 Energies")
#C2H4_PairPlot.savefig('PairPlotC2H4.png')

#-----------------------------------------------------------------------
#print(C2H4_ads - C2H4_free - zeolite)  #same as y[:,0]
#--------------------------------------------------------------------------
#H2O
H2O_1_egy = np.genfromtxt('BEEF_vdw_H2O-1.csv') + 0.682 - 0.629
#H2O_2_egy = np.genfromtxt('BEEF_vdw_H2O-2.csv') + 2*(0.682 - 0.629)
#H2O_3_egy = np.genfromtxt('BEEF_vdw_H2O-3.csv') + 3*(0.682 - 0.629)
H2O_2_egy = np.genfromtxt('BEEF_vdw_H2O-2.csv') + 1.35 - 1.24
H2O_3_egy = np.genfromtxt('BEEF_vdw_H2O-3.csv') + 0.0896 - (-0.073)
H2O_free_np = np.genfromtxt('BEEF_vdw_ensemble_energies_H2O_free.csv') + 0.6034 - 0.58352
zeolite_np = np.genfromtxt('BEEF_vdw_Zeolite.csv')


#y_energies is delta G
y2_energies = np.zeros((len(H2O_1_egy),3))
y2_energies[:,0] = H2O_1_egy - H2O_free_np - zeolite_np
y2_energies[:,1] = H2O_2_egy - H2O_1_egy - H2O_free_np
y2_energies[:,2] = H2O_3_egy - H2O_2_egy -  H2O_free_np
X = np.array([0,1/3,2/3])
X = X.reshape((-1,1))
M_2 = np.zeros(len(H2O_1_egy))
for i in range(len(H2O_1_egy)): 
  y2 = y2_energies[i,:]
  y2 = y2-y2[0] # has to be through origin
  	#y = y.reshape((1,-1))
  	#reg = LinearRegression().fit(X, y)
  	#m = reg.coef_  # will be one entry in delta_G_m_2
  m2, _, _, _ = np.linalg.lstsq(X, y2)
  M_2[i] = m2
#----------------------------------------------------------------
#print(H2O_ads - H2O_free - zeolite) #same as y[:,0]
#----------------------------------------------------------------
H2O_eng = np.concatenate(y2_energies).reshape((2000,3))
#-----------------------------------------------------------------------
#PairPlot
#H2Opd_eng = pd.DataFrame(np.array(H2O_eng), columns=['$G_{rxn1}$ (eV)', '$G_{rxn2}$ (eV)', '$G_{rxn3}$ (eV)'])
#H2O_PairPlot = sns.pairplot(H2Opd_eng)
#H2O_PairPlot.set(xlim=(-2,2.5))
#H2O_PairPlot.set(ylim=(-2,2.5))
#H2O_PairPlot.fig.suptitle("H2O Energies")
#H2O_PairPlot.savefig('PairPlotH2O.png')

#---------------------lateral interaction microkinetic model

def microkinetic_model(y, t, C2H4_gas_frac,H2O_gas_frac,G_1,f1,m_1,G_2,f2,m_2):
    theta_C2H4, theta_H2O, theta_vacant = y
    K1 = np.exp(-(G_1+theta_C2H4*m_1)/(kB*298.15))
    K2 = np.exp(-(G_2+theta_H2O*m_2)/(kB*298.15))
    r1 = f1/K1
    r2 = f2/K2
    d_theta_C2H4_dt = f1*C2H4_gas_frac*theta_vacant - r1*theta_C2H4
    d_theta_H2O_dt = f2*H2O_gas_frac*theta_vacant - r2*theta_H2O
    d_theta_vacant_dt = -d_theta_C2H4_dt-d_theta_H2O_dt
    dydt = [d_theta_C2H4_dt, d_theta_H2O_dt, d_theta_vacant_dt]
    return dydt
C2H4_fractions = np.zeros(len(M_1))
H2O_fractions = np.zeros(len(M_1))
y0 = y=np.array([0, 0, 1])
t = np.linspace(0,1E5,1E4)
for i in range(len(M_1)):  #M_1 holds all slopes for C2H4.
    f1 = 6.2014895E5 #adsorption rate constant ethylene on zeolite                    # from collision theory
    f2 = 7.7382354E5 #adsorption rate water from collision theory
    G_1 = y_energies[i,0]
    m_1 = M_1[i]
    G_2 = y2_energies[i,0]
    m_2 = M_2[i]
    theta_individual = odeint(microkinetic_model, y0, t, args=(6.0E-4, 0.06, G_1,f1,m_1,G_2,f2,m_2))
    #print(str(theta_individual[-1,:]))
    C2H4_fractions[i] = theta_individual[-1,0]
    H2O_fractions[i] = theta_individual[-1,1]
adsorption_fractions = np.column_stack((C2H4_fractions,H2O_fractions))
print('adsadsorption_fractions.shape: ' + str(adsorption_fractions.shape))

#C2H4_fractions.flatten()
#H2O_fractions.flatten()

C2H4_fractions = np.reshape(C2H4_fractions, (2000,1))
H2O_fractions = np.reshape(H2O_fractions, (2000,1)) 


print(C2H4_fractions)
print(C2H4_fractions.shape)
print(H2O_fractions)
print(H2O_fractions.shape)


#------------------
#Plot
#preset_bins = np.linspace(0,1,10)
#matplotlib.rc('xtick', labelsize=14)
#matplotlib.rc('ytick', labelsize=14)

#fig, ax = plt.subplots()
#ax.hist([C2H4_fractions, H2O_fractions], bins=preset_bins, label=[r'$\theta_{C2H4}$',r'$\theta_{H2O}$'])
#ax.set_title('Microkinetic Model Adsorption Fractions',  fontsize=18)
#ax.set_ylabel('Count', fontsize=18)
#ax.set_xlabel('Surface Fraction', fontsize=18)
#ax.set_xlim(-.02,1.02)
#ax.set_ylim(0,2000)
#plt.legend(loc='upper right', fontsize=16)
#fig.savefig('ThetaMicMod.png', dpi=220)


#------------------
#Langmuir

K1_L = np.exp(-(C2H4_ads - C2H4_free - zeolite)/(kB*T))
K2_L = np.exp(-(H2O_ads - H2O_free - zeolite)/(kB*T))

P_C2H4 = 0.0006      # Can be 0.0006, 0.00045, 0.0003, 0.0002
P_H2O = 0.06
theta_C2H4 = np.array(K1_L*P_C2H4/(1 + K1_L*P_C2H4 + K2_L*P_H2O))
theta_H2O = np.array(K2_L*P_H2O/(1 + K1_L*P_C2H4 + K2_L*P_H2O))
#theta_C2H4.flatten()
#theta_H2O.flatten()

print(theta_C2H4)
print(theta_C2H4.shape)
print(theta_H2O)
print(theta_H2O.shape)


#-------------------------------------------------------------------
#Plot 

#matplotlib.rc('xtick', labelsize=14)     
#matplotlib.rc('ytick', labelsize=14)

#fig, ax = plt.subplots()
#ax.hist([theta_C2H4, theta_H2O], bins=preset_bins, label=[r'$\theta_{C2H4}$',r'$\theta_{H2O}$'])
#ax.set_title('Langmuir Isotherm  Model Adsorption Fractions', fontsize=18)
#ax.set_ylabel('Count', fontsize=18)
#ax.set_xlabel('Surface Fraction', fontsize=18)
#ax.set_xlim(-.02,1.02)
#ax.set_ylim(0,2000)
#plt.legend(loc='upper right', fontsize=16)
#fig.savefig('ThetaLang.png', dpi=220)


#-------------
#PairPlot
#ads_frac_pd = pd.DataFrame(adsorption_fractions, columns=[r'$\theta_{C_2H_4}$',r'$\theta_{H_2O}$'])
#ads_PairPlot = sns.pairplot(ads_frac_pd)
#ads_PairPlot.set(xlim=(0,1))
#ads_PairPlot.set(ylim=(0,1))
#ads_PairPlot.savefig('PairPlotAdsFrac.png')


#-------------
#Violin Plot

#C2H4_fractions_df = pd.DataFrame(C2H4_fractions, columns=['6e-4','4.5e-4','3e-4','2e-4'])   
#print(C2H4_fractions_df)

#C2H4_fractions    MicModel
#H2O_fractions     MicModel
#theta_C2H4        Lang
#theta_H2O         Lang

#stack1 = np.hstack((C2H4_fractions,H2O_fractions))
#stack2 = np.hstack((stack1, theta_C2H4))
#stack3 = np.hstack((stack2, theta_H2O))
#ALL_Frac = np.reshape(stack3, (2000,4))

# 0 = MicMod    1 = Langmuir
a = np.array([[0],[1]])
ModelType = np.repeat(a, 4000)
ModelType = np.reshape(ModelType, (8000,1)) 
print(ModelType)
print(ModelType.shape)

# 0 = C2H4     1 = H2O
b = np.array([[0],[1]])
b1 = np.repeat(b, 2000)
Adsorbate = np.vstack((b1,b1))
Adsorbate = np.reshape(Adsorbate, (8000,1))
print(Adsorbate)
print(Adsorbate.shape)

stack1 = np.vstack((C2H4_fractions,H2O_fractions))
stack2 = np.vstack((stack1, theta_C2H4))
stack3 = np.vstack((stack2, theta_H2O))
ALL_Frac = np.reshape(stack3, (8000,1))
print(ALL_Frac)
print(ALL_Frac.shape)

dstack1 = np.hstack((ALL_Frac,ModelType))
Data = np.hstack((dstack1,Adsorbate))
Data = np.reshape(Data, (8000,3))

Fractions = pd.DataFrame(Data, columns=['Adsorption Fraction','Model Type','Adsorbate'])
print('Langmuir_fractions',Fractions)
Fractions.to_excel('Langmuir_fractions.xls')

sns.set(style="whitegrid", font_scale=1)
fig,ax = plt.subplots()
sns.violinplot(x="Model Type", y="Adsorption Fraction", hue="Adsorbate", inner=None, split=True, palette={0: "r", 1: "b"}, data = Fractions, bw=.2, cut=1, linewidth=1) 
sns.despine(left=True, bottom=True)
ax.set_ylabel('Adsorption Fraction')
ax.set_xlabel('Model Type')
fig.savefig('violinSurfFrac.png',dpi=220)




#fig,ax = plt.subplots()
#sns.violinplot(data=C2H4_fractions_beta_df,palette="Set3", bw=.2, cut=1, linewidth=1)
#sns.despine(left=True, bottom=True)
#ax.set_title('$C_2H_4$')
#ax.set_ylabel('Adsorption Fractions')
#ax.set_xlabel('$C_2H_4$ Concentration')
#fig.savefig('violinC2H4.png',dpi=220)

#print('H2O_fractions_beta',H2O_fractions_beta)
#H2O_fractions_beta_df = pd.DataFrame(H2O_fractions_beta, columns=['6e-4','4.5e-4','3e-4','2e-4'])
#print(H2O_fractions_beta_df)
#fig,ax = plt.subplots()
#sns.violinplot(data=H2O_fractions_beta_df,palette="Set3", bw=.2, cut=1, linewidth=1)
#sns.despine(left=True, bottom=True)
#ax.set_title('$H_2O$')
#ax.set_ylabel('Adsorption Fractions')
#ax.set_xlabel('$C_2H_4$ Concentration')
#fig.savefig('violinH2O.png',dpi=220)


