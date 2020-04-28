# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 11:45:26 2019
These scripts are used to:
    1) Run pySICE on dat files for debugging purpose 
    2) Plot the output of both pySICE and sice.f as well as scattered plot
        This part only works when having the output of pySICE loaded in workspace
        and with the good path set to the sice.f output
@author: bav
"""

#%% ========= run from data file for debug ================
olci_data = genfromtxt('../SnowProcessor/2.2/olci_toa_newformat_polluted.dat')
water_vod = genfromtxt('../SnowProcessor/2.2/tg_water_vod.dat', delimiter='   ')
voda = water_vod[range(21),1]

ozone_vod = genfromtxt('../SnowProcessor/2.2/tg_vod.dat', delimiter='   ',skip_header=2)
tozon = ozone_vod[range(21),1]
aot = 0.1

for i_pixel in range(6):
    ns, lat, lon, sza,vza,saa,vaa,height = olci_data[i_pixel,range(8)]
    toa = olci_data[i_pixel,range(8,29)]
    ozone =  olci_data[i_pixel,29] #kg m-2
    water = olci_data[i_pixel,30]  # kg/m2

    toa_uncor = toa
#    BXXX[i], D[i], area[i], al[i], r0[i], isnow[i], conc[i], ntype[i],alb_sph[i,:], rp[i,:],refl[i,:], rp1[i], rp2[i], rp3[i], rs1[i], rs2[i], rs3[i] = \
    sl.pySICE(sza,vza,saa,vaa,height,toa,ozone,water,voda,tozon,aot)
    
    #%% Comparing sice output with Python output
import matplotlib.pyplot as plt
    
def compare_plot(var, file_name, var_name, var_id):
    x, y, var_f = output_sice_f(file_name,var_name,var_id)
    x1,y1,var_py = WriteOutput(var,(var_name+'_py'),InputFolder)
    print(np.shape(var_f))
    print(np.shape(var_py))
    
    if (np.shape(var_f)[0]>np.shape(var_py)[0]):
        var_f=var_f[range(np.shape(var_py)[0]), :]
    if (np.shape(var_f)[1]>np.shape(var_py)[1]):
        var_f=var_f[:,range(np.shape(var_py)[1])]

    x0 = min(np.nanmin(var_f),np.nanmin(var_py))
    x1 = max(np.nanmax(var_f),np.nanmax(var_py))
    fig, ax = plt.subplots()
    ax.scatter(var_f.flatten(), var_py.flatten())
    ax.plot([x0, x1], [x0, x1], 'k-', color = 'r')
    plt.title(var_name)
    plt.xlabel('Fortran')
    plt.ylabel('Python')
    plt.show()

# "../SnowProcessor/2.2/bba.dat")
#ns,ndate(3),alat, alon,rp3,rp1,rp2,rs3,rs1, rs2,isnow
#rp3[rp3==0]=np.nan
#rp1[rp1<1]=np.nan
#rp2[rp2<1]=np.nan
#rs3[rs3<1]=np.nan
compare_plot(rs3, "../SnowProcessor/2.2/bba.dat", 'rs3', 7)
compare_plot(rp1, "../SnowProcessor/2.2/bba.dat", 'rp1', 5)
compare_plot(rp2, "../SnowProcessor/2.2/bba.dat", 'rp2', 6)
compare_plot(rp3, "../SnowProcessor/2.2/bba.dat", 'rp3', 4)

# NOT GOOD

compare_plot(isnow, "../SnowProcessor/2.2/bba.dat", 'isnow', 10)
# GOOD

# ("../SnowProcessor/2.2/size.dat")
#ns,ndate(3),alat,alon,D,area,al,r0, andsi,andbi,indexs,indexi,indexd,isnow
compare_plot(D, "../SnowProcessor/2.2/size.dat", 'D', 4)
compare_plot(area, "../SnowProcessor/2.2/size.dat", 'area', 5)
compare_plot(al, "../SnowProcessor/2.2/size.dat", 'al', 6)
compare_plot(r0, "../SnowProcessor/2.2/size.dat", 'r0', 7)
# ALL GOOD

# ("../SnowProcessor/2.2/retrieved_O3.dat")
# ns,alat,alon,BXXX,totadu,deltak,sza,vza,amf
BXXX_save = BXXX
BXXX[np.where(BXXX>500)] = np.nan
BXXX[np.where(BXXX<0)] = np. nan
compare_plot(BXXX, "../SnowProcessor/2.2/retrieved_O3.dat", 'O3', 4)
# output_sice_f("../SnowProcessor/2.2/retrieved_O3.dat",'O3',4)

# ('spherical_albedo.dat' )
# ns,ndate(3),alat,alon,(answer(i),i=1,21),isnow
for i in range(15,21):
    alb_sph[alb_sph[:,i]>100,i]=np.nan
    compare_plot(alb_sph[:,i], "../SnowProcessor/2.2/spherical_albedo.dat", ('alb_sph'+str(i)), 4+i)
    input("Press Enter to continue...")
# GOOD
    
for i in range(15,21):
    rp[rp[:,i]>100,i]=np.nan
    compare_plot(rp[:,i], "../SnowProcessor/2.2/planar_albedo.dat", ('alb_pl'+str(i)), 4+i)
    input("Press Enter to continue...")
# GOOD    
    
for i in range(15,21):
    refl[refl[:,i]>100,i]=np.nan
    compare_plot(refl[:,i], "../SnowProcessor/2.2/boar.dat", ('refl'+str(i)), 4+i)
    input("Press Enter to continue...")
# GOOD

#  open( 'impurity.dat'              )
# ns,ndate(3),alat,alon,ntype,conc,bf,bm,thv,  toa(1),isnow
compare_plot(ntype, "../SnowProcessor/2.2/impurity.dat", 'ntype', 4)
compare_plot(conc, "../SnowProcessor/2.2/impurity.dat", 'conc', 5)
# GOOD
