# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 16:49:50 2020

@author: bav
"""

#%% Test for integration
def cube(x):
    return x**3

def res(a,b):
    return b**4/4-a**4/4
a = 0
b=0.3

start_time = time.time()
test1 = sl.qsimp(cube,a,b)
print(test1)
print(time.time() - start_time)


start_time = time.time()
test2 = res(a,b)
print(test2)
print(time.time() - start_time)

import scipy.integrate as Integrate
start_time = time.time()
test3 = Integrate.quad(cube, a,b)
print(test3)
print(time.time() - start_time)

# Conclusion: qsimp is a good integration function