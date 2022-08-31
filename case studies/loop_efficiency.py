# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:14:57 2020

@author: bav
"""
#%%
import numpy as np
import time

mat = np.arange(21 * 841 * 1029).reshape(21, 841, 1029) * 145.001
coef = np.arange(21)
mat2 = mat * np.nan


def mult_channel(c, A):
    tmp = A.T * c
    return tmp.T


nloop = 10

# case 1: loop
start_time = time.process_time()
for k in range(nloop):
    for i in range(21):
        mat2[i, :, :] = np.exp(np.sqrt(coef[i] * mat[i, :, :]))
end_time = time.process_time()
print("1")
print(end_time - start_time)

# case 2: function
start_time = time.process_time()
for k in range(nloop):
    mat2 = np.exp(np.sqrt(mult_channel(coef, mat)))
end_time = time.process_time()
print("2")
print(end_time - start_time)

# case 2: function
mat3 = np.arange(21 * 841 * 1029).reshape(841, 1029, 21) * 145.001

start_time = time.process_time()
for k in range(nloop):
    mat4 = np.exp(np.sqrt(coef, mat3 * coef))
end_time = time.process_time()
print("3")
print(end_time - start_time)

#%%
# def mult_channel(c,A):
#    tmp = A.T*c
#    return tmp.T
#
# alb_sph = np.exp(-np.sqrt(1000.*4.* mult_channel(np.pi*bai/w, np.tile(al,(21,1,1)))))
#
##alpha = 4.*np.pi*bai/w
##absor = g*np.nan
##for i_channel in range(21):
##    alb_sph[i_channel, :,:] = np.exp(-np.sqrt(1000.*alpha[i_channel]*al))
##    alb_sph[i_channel, absor[i_channel,:,:] <=  1.e-6] = 1.0
#
##%%
#
# start_time= time.process_time()
# for i_channel in range(21):
#    rp[i_channel,:,:] = alb_sph[i_channel,:,:]**ak1
#
## derivation of snow reflectance function
#    refl[i_channel,:,:]=r0*alb_sph[i_channel,:,:]**(ak1*ak2/r0)
# end_time = time.process_time()
# print(end_time-start_time)
#
#
# start_time= time.process_time()
# rp_2 = np.power (alb_sph, ak1)
## derivation of snow reflectance function
# refl_2 =r0* np.power(alb_sph, (ak1*ak2/r0))
#
# end_time = time.process_time()
# print(end_time-start_time)
#
# np.nanmax(rp-rp_2)
# np.nanmax(refl-refl_2)
