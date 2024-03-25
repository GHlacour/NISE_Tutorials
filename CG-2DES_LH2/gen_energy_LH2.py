# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 09:30:21 2020

@author: kaizhong
"""

import numpy as np
import random
import math
import re
import itertools


num_step = 1000000
lambda1 = lambda2 = 1/150

delta_t = 2
sigma_D = 256 #disorder magnitudes for B850
sigma_A  = 169 #disorder magnitudes for B80
sampleNo = 100000
C = 2.99792458e-5
mu = 0  #average value equal zero
part_a_1 = math.exp(-delta_t*lambda1)
part_b_1 = ((1-(part_a_1**2))**(0.5))

part_a_2 = math.exp(-delta_t*lambda2)
part_b_2 = ((1-(part_a_2**2))**(0.5))

num_D = 18
num_A = 9

def Noise_data(Initial_value,sigma,part_a, part_b):
    def itfunc(a, b):
        return a*part_a + b
    mu = 0

    all_G = np.random.normal(mu, sigma, num_step)
    all_G *= part_b
    all_G[0] = Initial_value
    input_noise = [*itertools.accumulate(all_G, itfunc)]
    return input_noise

noise_total1  = []
noise_total2  = []

for i in range(num_D):
    Initial_value_delta_t2 =  np.random.normal(mu,sigma_D)
    noise2 = Noise_data(Initial_value_delta_t2,sigma_D,part_a_1, part_b_1)
    noise_total1.append(noise2)

for i in range(num_A):
    Initial_value_delta_t1 =  np.random.normal(mu,sigma_A)
    noise1 = Noise_data(Initial_value_delta_t1,sigma_A,part_a_1, part_b_1)
    noise_total2.append(noise1)


inputfile = open(r"input_eng.txt" ,"r")
dataline = inputfile.readlines()
#line_num = 200
sum_int  = 0
delta_t  = 3
input_ham = []
coor = re.split(r"\s+", dataline[0].strip())
np.array(coor,dtype=float)
count = 0
with open ('flu_energy.txt','a') as j:
    j.write(str(0))
    j.write("  ")

    j.write(" ".join(str(item) for item in coor))
    j.write("\n")
    for m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18, m19, m20, m21, m22, m23, m24, m25, m26 in zip(noise_total1[0],noise_total2[0],noise_total1[1],noise_total1[2],noise_total2[1],noise_total1[3],noise_total1[4],noise_total2[2],noise_total1[5],noise_total1[6],noise_total2[3],noise_total1[7],noise_total1[8],noise_total2[4],noise_total1[9],noise_total1[10],noise_total2[5],noise_total1[11],noise_total1[12],noise_total2[6],noise_total1[13],noise_total1[14],noise_total2[7],noise_total1[15],noise_total1[16],noise_total2[8],noise_total1[17]):
        coor = re.split(r"\s+", dataline[0].strip())
        np.array(coor,dtype=float)
        #print(type(coor[1]))        
        coor[0]   = float(coor[0]  ) + float(m0 )
        coor[27]  = float(coor[27] ) + float(m1 )
        coor[53]  = float(coor[53] ) + float(m2 )
        coor[78]  = float(coor[78] ) + float(m3 )
        coor[102] = float(coor[102]) + float(m4 )
        coor[125] = float(coor[125]) + float(m5 )
        coor[147] = float(coor[147]) + float(m6 )
        coor[168] = float(coor[168]) + float(m7 )
        coor[188] = float(coor[188]) + float(m8 )
        coor[207] = float(coor[207]) + float(m9 )
        coor[225] = float(coor[225]) + float(m10)
        coor[242] = float(coor[242]) + float(m11)
        coor[258] = float(coor[258]) + float(m12)
        coor[273] = float(coor[273]) + float(m13)
        coor[287] = float(coor[287]) + float(m14)
        coor[300] = float(coor[300]) + float(m15)
        coor[312] = float(coor[312]) + float(m16)
        coor[323] = float(coor[323]) + float(m17)
        coor[333] = float(coor[333]) + float(m18)
        coor[342] = float(coor[342]) + float(m19)
        coor[350] = float(coor[350]) + float(m20)
        coor[357] = float(coor[357]) + float(m21)
        coor[363] = float(coor[363]) + float(m22)
        coor[368] = float(coor[368]) + float(m23)
        coor[372] = float(coor[372]) + float(m24)
        coor[375] = float(coor[375]) + float(m25)
        coor[377] = float(coor[377]) + float(m26)
        count +=1
        j.write(str(0))
        j.write("  ")
        j.write(" ".join(str(item) for item in coor))
        j.write("\n")
j.close()
