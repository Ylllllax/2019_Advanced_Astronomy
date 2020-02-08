#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import scipy.optimize as optimization
import warnings


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


# In[4]:


def return_template_information(st):
    st1 = st
    f = open(st1, 'r')
    str = f.read()
    f.close()
    alist = str.split("\n")
    alist = [i for i in alist if len(i)>0]
    xlist = []
    ylist = []   
    for i in range(len(alist)):
        singlelist = alist[i].split(" ")
        singlelist = [j for j in singlelist if len(j)>0]
        xlist.append(float(singlelist[0]))
        ylist.append(float(singlelist[1]))
    for i in range(len(ylist)):
        ylist[i] = ylist[i]*xlist[i]
    return xlist, ylist
     
def choice_list(xlist, ylist, model_xlist, model_ylist):
    choice_xlist = []
    choice_ylist = []
    new_ylist = []
    new_xlist = []
    m = model_xlist[0]
    for i in range(len(xlist)):
        if xlist[i] <= model_xlist[-1]:
            new_xlist.append(xlist[i])
            new_ylist.append(ylist[i])
        else:
            pass          
    for i in new_xlist:
        for j in range(len(model_xlist)):
            if abs(model_xlist[j]-i)<abs(m-i):
                m = model_xlist[j]
            else:
                pass
        choice_xlist.append(m)
        choice_ylist.append(model_ylist[model_xlist.index(m)])
    return choice_xlist, choice_ylist, new_ylist

def return_photoz(object_xlist, object_ylist, template_xlist, template_ylist):
    min_z = ''
    min_number = 50.0 
    for redshift_i in range(0, 5000):
        redshift = float(redshift_i)/1000
        object_new_xlist = [0.0 for m in range(len(object_xlist))]
        for i in range(len(object_ylist)):
            object_new_xlist[i] = object_xlist[i]/(1+redshift)  
        choice_xlist, choice_ylist, new_ylist = choice_list(object_new_xlist, object_ylist, template_xlist, template_ylist)
        def func(x,  a):
            return a*x
        x0 = 0.0
        a = optimization.curve_fit(func, choice_ylist, new_ylist, x0)
        choice_model_new_ylist = [0.0 for j in range(len(choice_ylist))]
        for n in range(len(choice_ylist)):
            choice_model_new_ylist[n] = a[0][0]*choice_ylist[n]
        number = 0.0
        for i in range(len(choice_model_new_ylist)):
            number = number+(new_ylist[i]-choice_model_new_ylist[i])*(new_ylist[i]-choice_model_new_ylist[i])
        if number < min_number:
            min_number = number 
            min_z = redshift 
        else:
            pass
    return min_z


# In[5]:


# 导入数据：HERSCHEL flux 单位：mJy；WISE 2MASS vega星等
data = pd.read_csv("flux_without_redshift.csv")
index = len(data)

for num in range(0, index):
    warnings.filterwarnings('ignore')
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
    Sum = np.sum(nuF)
    norm_nuF = [i/Sum for i in nuF]
    
    
    object_xlist, object_ylist = wavelength, norm_nuF  #归一化处理好的候选体数据，横坐标为波段，纵坐标为vFv
    
    template_xlist_1,template_ylist_1 = return_template_information('QSO1_template_norm.sed') #st为模板的路径
    photoz_1 = return_photoz(object_xlist,object_ylist,template_xlist_1,template_ylist_1)
    data['redshift1'][num] = photoz_1
    
    template_xlist_2,template_ylist_2 = return_template_information('Torus_template_norm.sed') #st为模板的路径
    photoz_2 = return_photoz(object_xlist,object_ylist,template_xlist_2,template_ylist_2)
    data['redshift2'][num] = photoz_2
    
data.to_csv('flux_with_photoz.csv', float_format="%.6f", index=0)


# In[ ]:




