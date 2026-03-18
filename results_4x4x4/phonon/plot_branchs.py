import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParamsDefault
import numpy as np
def Heavi(x,delta):
    if x<=delta/2:
        return 1
    else :
        return 0
import math
def Gaussian(x,delta):
    return math.exp(-0.5*x*x/(delta*delta))/np.sqrt(2*3.1459*delta*delta)


# Specify the filename
filename = 'Nb_modes.freq'

# Read all lines from the file
with open(filename, 'r') as f:
    lines = f.readlines()

# Find the &plot line to extract nbnd and nks
plot_line_index = None
for idx, line in enumerate(lines):
    if '&plot' in line:
        plot_line_index = idx
        match = re.search(r'nbnd=\s*(\d+),\s*nks=\s*(\d+)', line)
        nbnd, nks = map(int, match.groups())
        break

# Extract the data lines after the &plot line
data_lines = [l.strip() for l in lines[plot_line_index+1:] if l.strip()]

# Initialize lists to collect frequencies per branch
branches = {f'branch_{i+1}': [] for i in range(nbnd)}
qpoints = []

# Each q-point has 2 lines: coordinates and frequencies
for k in range(nks):
    coord_line = data_lines[2 * k]
    freq_line = data_lines[2 * k + 1]

    # Parse q-point coordinates
    qx, qy, qz = map(float, coord_line.split())
    qpoints.append((qx, qy, qz))

    # Parse frequencies
    freqs = list(map(float, freq_line.split()))
    for i, f_val in enumerate(freqs):
        branches[f'branch_{i+1}'].append(f_val)

branch1_freqs = branches['branch_1']  # this is already a Python list
Phonon_1 = [f * 0.02998 for f in branch1_freqs]  # in THz

branch2_freqs = branches['branch_2']  # this is already a Python list
Phonon_2 = [f * 0.02998 for f in branch2_freqs]  # in THz

branch3_freqs = branches['branch_3']  # this is already a Python list
Phonon_3 = [f * 0.02998 for f in branch3_freqs]  # in THz
# I now how to extract so the next is to copy the sma code for the DOS
NQ=len(branch1_freqs) #Number of K points used on the Calcualtions.
V=3# Number of branch for the Phonons
frequency=[]
for i in range(0,NQ):
    frequency.append(i*0.001)

plt.plot(frequency, Phonon_1,'o-')
plt.plot(frequency, Phonon_2, linewidth=1, alpha=0.5, color='b')
plt.plot(frequency, Phonon_3, linewidth=1, alpha=0.5, color='g')


plt.xlabel("Wavevector")
plt.ylabel(" Frequency [THz]")
plt.legend()
plt.show()


print("Maximum Frequencies for Acoutic Phonons TF TS and L",max(Phonon_1),max(Phonon_2),max(Phonon_3))
# print("Maximum Frequencies for Optical Phonons TF TS and L",max(Phonon_4),max(Phonon_5),max(Phonon_6))

# Defining the Density of States vs the Frequency
# Defining the same number of frequency points as q points in a ra
# The density function depends of w.

minW=0.0
maxW=5.2
Ste=(maxW-minW)/(len(frequency)) #+100 WAS WORKING GOOD
print(Ste)
#Finding the range
W=np.arange(minW, maxW, Ste)

DeltaW= 0.1 #THz
Density=[]
Density_1=[]
Density_2=[]
Density_3=[]
# Density_4=[]
# Density_5=[]
# Density_6=[]
# # Density_7=[]
# # Density_8=[]
# # Density_9=[]
# # Density_10=[]
# # Density_11=[]
# # Density_12=[]
# # Density_13=[]
# # Density_14=[]
# # Density_15=[]
# # Density_16=[]
# # Density_17=[]
# # Density_18=[]
# # Density_19=[]
# # Density_20=[]
# # Density_21=[]
# # Density_22=[]
# # Density_23=[]
# # Density_24=[]
# # Density_25=[]
# # Density_26=[]
# # Density_27=[]
# # Density_28=[]
# # Density_29=[]
# # Density_30=[]
# # A1=0
# # A2=0
# # A3=0
# # A4=0
# # A5=0
# # A6=0
# # for i in range(0,len(W)):
# #     for j in range(0,len(frequency)):
# #         A1+=Heavi(abs(Phonon_1[j]-W[i]),DeltaW)
# #         A2+=Heavi(abs(Phonon_2[j]-W[i]),DeltaW)
# #         A3+=Heavi(abs(Phonon_3[j]-W[i]),DeltaW)
# #         A4+=Heavi(abs(Phonon_4[j]-W[i]),DeltaW)
# #         A5+=Heavi(abs(Phonon_5[j]-W[i]),DeltaW)
# #         A6+=Heavi(abs(Phonon_6[j]-W[i]),DeltaW)
# #     Density.append(A1+A2+A3+A4+A5+A6)
# #     Density_1.append(A1)
# #     Density_2.append(A2)
# #     Density_3.append(A3)
# #     A1=0
# #     A2=0
# #     A3=0
# #     A4=0
# #     A5=0
# #     A6=0
A1=0
A2=0
A3=0
# A4=0
# A5=0
# A6=0
# # A7=0
# # A8=0
# # A9=0
# # A10=0
# # A11=0
# # A12=0
# # A13=0
# # A14=0
# # A15=0
# # A16=0
# # A17=0
# # A18=0
# # A19=0
# # A20=0
# # A21=0
# # A22=0
# # A23=0
# # A24=0
# # A25=0
# # A26=0
# # A27=0
# # A28=0
# # A29=0
# # A30=0
#
for i in range(0,len(W)):
    for j in range(0,len(frequency)):
        A1+=Gaussian(abs(Phonon_1[j]-W[i]),DeltaW)
        A2+=Gaussian(abs(Phonon_2[j]-W[i]),DeltaW)
        A3+=Gaussian(abs(Phonon_3[j]-W[i]),DeltaW)
#         A4+=Gaussian(abs(Phonon_4[j]-W[i]),DeltaW)
#         A5+=Gaussian(abs(Phonon_5[j]-W[i]),DeltaW)
#         A6+=Gaussian(abs(Phonon_6[j]-W[i]),DeltaW)
#         # A7+=Gaussian(abs(Phonon_7[j]-W[i]),DeltaW)
#         # A8+=Gaussian(abs(Phonon_8[j]-W[i]),DeltaW)
#         # A9+=Gaussian(abs(Phonon_9[j]-W[i]),DeltaW)
#         # A10+=Gaussian(abs(Phonon_10[j]-W[i]),DeltaW)
#         # A11+=Gaussian(abs(Phonon_11[j]-W[i]),DeltaW)
#         # A12+=Gaussian(abs(Phonon_12[j]-W[i]),DeltaW)
#         # A13+=Gaussian(abs(Phonon_13[j]-W[i]),DeltaW)
#         # A14+=Gaussian(abs(Phonon_14[j]-W[i]),DeltaW)
#         # A15+=Gaussian(abs(Phonon_15[j]-W[i]),DeltaW)
#         # A16+=Gaussian(abs(Phonon_16[j]-W[i]),DeltaW)
#         # A17+=Gaussian(abs(Phonon_17[j]-W[i]),DeltaW)
#         # A18+=Gaussian(abs(Phonon_18[j]-W[i]),DeltaW)
#         # A19+=Gaussian(abs(Phonon_19[j]-W[i]),DeltaW)
#         # A20+=Gaussian(abs(Phonon_20[j]-W[i]),DeltaW)
#         # A21+=Gaussian(abs(Phonon_21[j]-W[i]),DeltaW)
#         # A22+=Gaussian(abs(Phonon_22[j]-W[i]),DeltaW)
#         # A23+=Gaussian(abs(Phonon_23[j]-W[i]),DeltaW)
#         # A24+=Gaussian(abs(Phonon_24[j]-W[i]),DeltaW)
#         # A25+=Gaussian(abs(Phonon_25[j]-W[i]),DeltaW)
#         # A26+=Gaussian(abs(Phonon_26[j]-W[i]),DeltaW)
#         # A27+=Gaussian(abs(Phonon_27[j]-W[i]),DeltaW)
#         # A28+=Gaussian(abs(Phonon_28[j]-W[i]),DeltaW)
#         # A29+=Gaussian(abs(Phonon_29[j]-W[i]),DeltaW)
#         # A30+=Gaussian(abs(Phonon_30[j]-W[i]),DeltaW)
    Density.append(A1+A2+A3)

    Density_1.append(A1)
    Density_2.append(A2)
    Density_3.append(A3)
#     Density_4.append(A4)
#     Density_5.append(A5)
#     Density_6.append(A6)
#     # Density_7.append(A7)
#     # Density_8.append(A8)
#     # Density_9.append(A9)
#     # Density_10.append(A10)
#     # Density_11.append(A11)
#     # Density_12.append(A12)
#     # Density_13.append(A13)
#     # Density_14.append(A14)
#     # Density_15.append(A15)
#     # Density_16.append(A16)
#     # Density_17.append(A17)
#     # Density_18.append(A18)
#     # Density_19.append(A19)
#     # Density_20.append(A20)
#     # Density_21.append(A21)
#     # Density_22.append(A22)
#     # Density_23.append(A23)
#     # Density_24.append(A24)
#     # Density_25.append(A25)
#     # Density_26.append(A26)
#     # Density_27.append(A27)
#     # Density_28.append(A28)
#     # Density_29.append(A29)
#     # Density_30.append(A30)
    A1=0
    A2=0
    A3=0
#     A4=0
#     A5=0
#     A6=0
#     # A7=0
#     # A8=0
#     # A9=0
#     # A10=0
#     # A11=0
#     # A12=0
#     # A13=0
#     # A14=0
#     # A15=0
#     # A16=0
#     # A17=0
#     # A18=0
#     # A19=0
#     # A20=0
#     # A21=0
#     # A22=0
#     # A23=0
#     # A24=0
#     # A25=0
#     # A26=0
#     # A27=0
#     # A28=0
#     # A29=0
#     # A30=0
#
# #Checking Heaviside
# print(len(Density))
# print(len(W))
plt.plot(W,Density,label='Total DOS')
plt.plot(W,Density_1,label='TF')
plt.plot(W,Density_2,label='TS')
plt.plot(W,Density_3,label='LA')
# plt.plot(W,Density_4,label='TF')
# plt.plot(W,Density_5,label='TS')
# plt.plot(W,Density_6,label='LA')
# #plt.yscale("log")
plt.show()
import csv
Radial_posTS='DOS_Full_q_real.txt'
with open(Radial_posTS, 'a') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(W,Density_1,Density_2,Density_3))
