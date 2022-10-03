# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
x = np.linspace(0.33, 2.4,200)


# u1 between 0.9 and 1.1
# al between o and 70

fig, ax = plt.subplots(3,1, sharex=True)
ax[0].plot(xa,ya)
ax[0].set_yscale('log')
ax[1].plot(x, f0 + f1 * np.exp(-x * bet) + f2 * np.exp(-x * gam))
for k in range(1, 70,10):
    print(k)
    rs = x*np.nan
    for i in range(len(x)):
        rs[i] = funp(x[i],k,1,1)
    ax[2].plot(x, rs)

# %% 
def plane_albedo_sw_approx(al, u1):
    # convert al to D as it was done in 2021
    D = al/(9.2*16/9)
    # convert u1 to cos_sza as it was done in 2021
    cos_sza = (u1*7/3-1)/2
    anka = 0.7389 - 0.1783 * cos_sza + 0.0484 * cos_sza**2.
    banka = 0.0853 + 0.0414 * cos_sza - 0.0127 * cos_sza**2.
    canka = 0.1384 + 0.0762 * cos_sza - 0.0268 * cos_sza**2.
    diam1 = 187.89 - 69.2636 * cos_sza + 40.4821 * cos_sza**2.
    diam2 = 2687.25 - 405.09 * cos_sza + 94.5 * cos_sza**2.
    return anka + banka * np.exp(-1000 * D / diam1) + canka * np.exp(-1000 * D / diam2)


def spher_albedo_sw_approx(al):
    # convert al to D as it was done in 2021
    D = al/(9.2*16/9)
    anka = 0.6420
    banka = 0.1044
    canka = 0.1773
    diam1 = 158.62
    diam2 = 2448.18
    return anka + banka * np.exp(-1000 * D / diam1) + canka * np.exp(-1000 * D / diam2)

# %% Evaluation of BBA approximation planar
u1_list = np.arange(0.9, 1.1,0.01)
plt.close('all')
for al in [1, 10, 20, 30, 50]: 

    col = plt.cm.seismic(np.arange(len(u1_list))/len(u1_list))    
    import matplotlib as mpl
    fig, ax = plt.subplots(2,1, figsize=(6,10))
    fig.suptitle('Evaluation of planar BBA integration: al = '+str(al))
    for k, u1 in enumerate(u1_list):
        rp = x*np.nan
        for i in range(len(x)):
            rp[i] = funp(x[i], al, 0, u1) / (f0 + f1 * np.exp(-x[i] * bet) + f2 * np.exp(-x[i] * gam))
        ax[0].plot(x, rp, color=col[k,:])
    
            
        _, _, bba = BBA_calc_clean(al, u1, mode = 'planar')
        ax[1].plot(u1, bba, color='tab:red', marker='o', linestyle='None', label='_no_legend_')
        bba_approx = 0.5271 + 0.3612 * np.exp(-u1 * np.sqrt(0.02350 * al))
        ax[1].plot(u1, bba_approx, color='tab:blue', marker='^', linestyle='None', label='_no_legend_')
        ax[1].plot(u1, plane_albedo_sw_approx(al, u1), 
                   color='tab:green', marker='v', linestyle='None', label='_no_legend_')
    
        alb = wls.values*np.nan
        alb[0] = np.interp(wls[0], x, rp)
        alb[5] = np.interp(wls[5], x, rp)
        alb[10] = np.interp(wls[10], x, rp)
        alb[11] = np.interp(wls[11], x, rp)
        alb[16] = np.interp(wls[16], x, rp)
        alb[20] = np.interp(wls[20], x, rp)
        
        _, _, bba_sw_pol = BBA_calc_pol(np.expand_dims(alb, 1), asol, sol_vis, sol_nir, sol_sw)
        ax[1].plot(u1, bba_sw_pol, color='tab:orange', marker='d', linestyle='None', label='_no_legend_')
  
    ax_cb = fig.add_axes([0.7, 0.8, 0.15, 0.05])
    cmap = mpl.cm.seismic
    norm = mpl.colors.Normalize(vmin=np.min(u1_list), vmax=np.max(u1_list))
    cb1 = mpl.colorbar.ColorbarBase(ax_cb, cmap=cmap,
                                    norm=norm,
                                    orientation='horizontal',
                                    label='u1 value')
    
    ax[1].plot(np.nan, np.nan, marker='o', color='tab:red', linestyle='None', label='integration')
    ax[1].plot(np.nan, np.nan, marker='^', color='tab:blue', linestyle='None', label='approximation version 2022')
    ax[1].plot(np.nan, np.nan, marker='v', color='tab:green', linestyle='None', label='approximation version 2021')
    ax[1].plot(np.nan, np.nan, marker='d', color='tab:orange', linestyle='None', label='approximation piecewise')
    ax[1].set_ylim(0.6,0.95)
    ax[1].legend()
    ax[1].set_xlabel('u1 \n escape function applied to the cos of SZA')
    ax[1].set_ylabel('Planar SW BBA (-)')
    ax[0].set_xlabel('Wavelength ($\mu$m)')
    ax[0].set_ylabel('Planar albedo (-)')
    
# %% Evaluation of BBA approximation spherical
al_list = [1, 10, 20, 30, 40, 50]
plt.close('all')
x = np.linspace(0.33,2.5,100)
col = plt.cm.seismic(np.arange(len(al_list))/len(al_list))    
import matplotlib as mpl

fig, ax = plt.subplots(1,2, figsize=(10,6))
fig.suptitle('Evaluation of spherical BBA integration')
for k, al in enumerate(al_list):
    rs = x*np.nan
    for i in range(len(x)):
        if x[i]>=0.33:
            solar_flux = (f0 + f1 * np.exp(-x[i] * bet) + f2 * np.exp(-x[i] * gam))
        else:
            solar_flux = 0
        rs[i] = funp(x[i], al, 1, u1) / solar_flux
    ax[0].plot(x, rs, color=col[k,:])
        
    _, _, bba = BBA_calc_clean(al, u1, mode = 'spherical')
    ax[1].plot(al, bba, color='tab:red', marker='o', linestyle='None', label='_no_legend_')
    
    bba_approx = 0.5271 + 0.3612 * np.exp(-np.sqrt(0.02350 * al))
    ax[1].plot(al, bba_approx, color='tab:blue', marker='^', linestyle='None', label='_no_legend_')
    ax[1].plot(al, spher_albedo_sw_approx(al), color='tab:green', marker='v', linestyle='None', label='_no_legend_')

    alb = wls.values*np.nan
    alb[0] = np.interp(wls[0], x, rs)
    alb[5] = np.interp(wls[5], x, rs)
    alb[10] = np.interp(wls[10], x, rs)
    alb[11] = np.interp(wls[11], x, rs)
    alb[16] = np.interp(wls[16], x, rs)
    alb[20] = np.interp(wls[20], x, rs)
    
    _, _, bba_sw_pol = BBA_calc_pol(np.expand_dims(alb, 1), asol, sol_vis, sol_nir, sol_sw)
    ax[1].plot(al, bba_sw_pol, color='tab:orange', marker='d', linestyle='None', label='_no_legend_')

ax_cb = fig.add_axes([0.3, 0.8, 0.15, 0.05])
cmap = mpl.cm.seismic
norm = mpl.colors.Normalize(vmin=np.min(al_list), vmax=np.max(al_list))
cb1 = mpl.colorbar.ColorbarBase(ax_cb, cmap=cmap,
                                norm=norm,
                                orientation='horizontal',
                                label='al value')

ax[1].plot(np.nan, np.nan, marker='o', color='tab:red', linestyle='None', label='integration')
ax[1].plot(np.nan, np.nan, marker='^', color='tab:blue', linestyle='None', label='approximation version 2022')
ax[1].plot(np.nan, np.nan, marker='v', color='tab:green', linestyle='None', label='approximation version 2021')
ax[1].plot(np.nan, np.nan, marker='d', color='tab:orange', linestyle='None', label='approximation piecewise')
ax[1].legend()
ax[1].set_xlabel('al')
ax[1].set_ylabel('Spherical SW BBA (-)')
ax[0].set_xlabel('Wavelength ($\mu$m)')
ax[0].set_ylabel('Spherical albedo (-)')

# %%
x = xa
al= 10

flux_out = x*np.nan
rs = x*np.nan
solar_flux = x*np.nan

for i in range(len(x)):
    y = np.interp(x[i], xa, ya)
    dega = 1000.0 * al * 4.0 * np.pi * y / x[i]
    pow = np.sqrt(dega)
    if pow >= 1.0e-6:
        rs[i] = np.exp(-pow)
    else:
        rs[i] = 1.0
    flux_out[i] = funp(x[i], 10, 1, 1)

# input reflectances
alb = wls.values*np.nan
alb[0] = np.interp(wls[0], x, rs)
alb[5] = np.interp(wls[5], x, rs)
alb[10] = np.interp(wls[10], x, rs)
alb[11] = np.interp(wls[11], x, rs)
alb[16] = np.interp(wls[16], x, rs)
alb[20] = np.interp(wls[20], x, rs)

# QUADRATIC POLYNOMIal for the range 400-709nm
_, a1, b1, c1 = quad_func(0.4, 0.56, 0.709, alb[0], alb[5], alb[10])
_, a2, b2, c2 = quad_func(0.709, 0.753, 0.865, alb[10], alb[11], alb[16])

rati = alb[16] / alb[20]
alasta = (1.02 - 0.865) / np.log(rati)
an = 1.0 / alasta
p = alb[16] * np.exp(0.865 / alasta)

def rs_fit(x):
    if x < 0.709:
        return a1 + b1 * x + c1 * x**2
    elif x < 0.865:
        return a2 + b2 * x + c2 * x**2
    else:
        return p * np.exp(-an * x)
rs_fit_v = np.vectorize(rs_fit)

plt.figure()
plt.plot(x, rs, label = 'exact')
plt.plot(x, rs_fit_v(x), label = 'fitted with polynomial and exponential functions')
plt.plot([0.4, 0.56, 0.709, 0.753, 0.865, 1.02], 
         [alb[0], alb[5], alb[10], alb[11], alb[16], alb[20]],
          marker = 'o', linestyle = 'None',
          label = '($\lambda$, rs) pairs used for the fitting')
plt.ylabel('Spherical albedo')
plt.ylabel('$\lambda (\mu m)$')
plt.legend()

print()
print('Approx 2022')
%timeit bba_approx = 0.5271 + 0.3612 * np.exp(-np.sqrt(0.02350 * al))
bba_approx = 0.5271 + 0.3612 * np.exp(-np.sqrt(0.02350 * al))
print(bba_approx)
print()

# integrating BBA
def func_integ(x):
    flux = funp(x, al,1,1)
    return flux
print('Integration')
from scipy.integrate import quad

%timeit bba_sw_clean = qsimp(func_integ, 0.33,  2.4)/sol_sw
bba_sw_clean = qsimp(func_integ, 0.33,  2.4)/sol_sw
print(bba_sw_clean)
print()

print('Integration of the piecewise approx (Alex)')
%timeit _, _, bba_sw_pol = BBA_calc_pol(np.expand_dims(alb, 1), asol, sol_vis, sol_nir, sol_sw)
_, _, bba_sw_pol = BBA_calc_pol(np.expand_dims(alb, 1), asol, sol_vis, sol_nir, sol_sw)
print(bba_sw_pol)
print()