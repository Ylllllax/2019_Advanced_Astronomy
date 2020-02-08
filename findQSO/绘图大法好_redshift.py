#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt


# In[2]:


# WISE vega星等转流量

def flux_wise_w1(mag_vega):
    mag_ab = mag_vega + 2.699
    F = 10**(-(mag_ab+48.60)/2.5) * 10**23
    return F

def flux_wise_w2(mag_vega):
    mag_ab = mag_vega + 3.339
    F = 10**(-(mag_ab+48.60)/2.5) * 10**23
    return F

def flux_wise_w3(mag_vega):
    mag_ab = mag_vega + 5.174
    F = 10**(-(mag_ab+48.60)/2.5) * 10**23
    return F

def flux_wise_w4(mag_vega):
    mag_ab = mag_vega + 6.620
    F = 10**(-(mag_ab+48.60)/2.5) * 10**23
    return F


# In[3]:


# 2MASS vega星等转流量

def flux_2mass_J(mag_vega):
    mag_ab = mag_vega + 0.91
    F = 10**(-(mag_ab+48.60)/2.5) * 10**23
    return F

def flux_2mass_H(mag_vega):
    mag_ab = mag_vega + 1.39
    F = 10**(-(mag_ab+48.60)/2.5) * 10**23
    return F

def flux_2mass_K(mag_vega):
    mag_ab = mag_vega + 1.85
    F = 10**(-(mag_ab+48.60)/2.5) * 10**23
    return F


# In[ ]:


# 导入数据：HERSCHEL flux 单位：mJy；WISE 2MASS vega星等
data = pd.read_csv("flux_with_redshift.csv")
index = len(data)

for num in range(0, index):
    print ('The echo is', num, '/', index)
    flux = []
    nu = []
    wavelength = []

    
    # 波长
    c = 3 * 10**14  # um/s
    wavelength_um = [70, 100, 160, 250, 350, 500, 3.4, 4.2, 12, 22, 1.23, 1.66, 2.15]
    nu_all = [c/i for i in wavelength_um]
    wavelength_all = [i*10000 for i in wavelength_um]


    # 流量数据存在时提取
    if not np.isnan(data['herf_70'][num]):
        flux.append(data['herf_70'][num]/1000)
        nu.append(nu_all[0])
        wavelength.append(wavelength_all[0])
    if not np.isnan(data['herf_100'][num]):
        flux.append(data['herf_100'][num]/1000)
        nu.append(nu_all[1])
        wavelength.append(wavelength_all[1])
    if not np.isnan(data['herf_160'][num]):
        flux.append(data['herf_160'][num]/1000)
        nu.append(nu_all[2])
        wavelength.append(wavelength_all[2])
    if not np.isnan(data['herf_250'][num]):
        flux.append(data['herf_250'][num]/1000)
        nu.append(nu_all[3])
        wavelength.append(wavelength_all[3])
    if not np.isnan(data['herf_350'][num]):
        flux.append(data['herf_350'][num]/1000)
        nu.append(nu_all[4])
        wavelength.append(wavelength_all[4])
    if not np.isnan(data['herf_500'][num]):
        flux.append(data['herf_500'][num]/1000)
        nu.append(nu_all[5])
        wavelength.append(wavelength_all[5])
    if not np.isnan(data['wisevega_1'][num]):
        flux.append(flux_wise_w1(data['wisevega_1'][num]))
        nu.append(nu_all[6])
        wavelength.append(wavelength_all[6])
    if not np.isnan(data['wisevega_2'][num]):
        flux.append(flux_wise_w2(data['wisevega_2'][num]))
        nu.append(nu_all[7])
        wavelength.append(wavelength_all[7])
    if not np.isnan(data['wisevega_3'][num]):
        flux.append(flux_wise_w3(data['wisevega_3'][num]))
        nu.append(nu_all[8])
        wavelength.append(wavelength_all[8])
    if not np.isnan(data['wisevega_4'][num]):
        flux.append(flux_wise_w4(data['wisevega_4'][num]))
        nu.append(nu_all[9])
        wavelength.append(wavelength_all[9])
    if not np.isnan(data['2mvega_j'][num]):
        flux.append(flux_2mass_J(data['2mvega_j'][num]))
        nu.append(nu_all[10])
        wavelength.append(wavelength_all[10])
    if not np.isnan(data['2mvega_h'][num]):
        flux.append(flux_2mass_H(data['2mvega_h'][num]))
        nu.append(nu_all[11])
        wavelength.append(wavelength_all[11])
    if not np.isnan(data['2mvega_k'][num]):
        flux.append(flux_2mass_K(data['2mvega_k'][num]))
        nu.append(nu_all[12])
        wavelength.append(wavelength_all[12])
    nuF = np.multiply(np.array(nu), np.array(flux))
    if len(nuF) <= 1:
        continue
    midnum = math.ceil(len(nuF)/2)
    norm_nuF = [i/nuF[midnum] for i in nuF]



    # 去除红移
    z = data['redshift'][num]
    wavelength0 = [i/(1+z) for i in wavelength]


    # 模板QSO1
    temp_QSO1 = pd.read_csv("QSO1_template_norm.sed", sep='\s+', header=None, names=['wavelength1','Flambda1'])
    wavelength1 = temp_QSO1['wavelength1']
    Flambda1 = temp_QSO1['Flambda1']
    lambdaF1 = np.multiply(np.array(wavelength1), np.array(Flambda1))  # lambda F_lambda = nu F_nu
    a = 500
    norm_lambdaF1 = [i/lambdaF1[a] for i in lambdaF1]


    # 模板QSO2
    temp_QSO2 = pd.read_csv("Torus_template_norm.sed", sep='\s+', header=None, names=['wavelength2','Flambda2'])
    wavelength2 = temp_QSO2['wavelength2']
    Flambda2 = temp_QSO2['Flambda2']
    lambdaF2 = np.multiply(np.array(wavelength2), np.array(Flambda2))  # lambda F_lambda = nu F_nu
    a = 500
    norm_lambdaF2 = [i/lambdaF2[a] for i in lambdaF2]


    # 绘图
    fig= plt.figure(figsize=(15,6))
    
    plt.subplot(121)
    plt.loglog()
    plt.scatter(wavelength0, norm_nuF, color="r", linewidth=1)
    plt.plot(wavelength1, norm_lambdaF1, color="b", linewidth=1)
    plt.xlabel(r"$Wavelength(Angstrom)$",fontsize=20)
    plt.ylabel(r"$ Normalized \ \nu F_{\nu} $",fontsize=20)
    plt.title(r"$QSO1$",fontsize=25)
    plt.tick_params(labelsize=15)

    plt.subplot(122)
    plt.loglog()
    plt.scatter(wavelength0, norm_nuF, color="r", linewidth=1)
    plt.plot(wavelength2, norm_lambdaF2, color="b", linewidth=1)
    plt.xlabel(r"$Wavelength(Angstrom)$",fontsize=20)
    plt.ylabel(r"$ Normalized \ \nu F_{\nu} $",fontsize=20)
    plt.title(r"$QSO2$",fontsize=25)
    plt.tick_params(labelsize=15)

    fig.tight_layout()
    plt.savefig('SED/%s.jpg' % (num))
    plt.show()
    plt.close()


# In[ ]:




