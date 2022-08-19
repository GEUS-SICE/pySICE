                        PROGRAM SICE

c***********************************************      
c                  A. KOKHANOVSKY       
c                  RETRIEVAL OF SNOW PROPERTIES
c                  USING OLCI
c				   version 6.0     
c                  a.a.kokhanovsky@gmail.com
c                  April 15, 2022               
c***********************************************
     
      REAL TOA(21),bai(21),CALIB(21),wls(21)
      REAL CTOA(21),BCTOA(21),answ2(21),answer(21),boa(21)
      REAL TOK(21),F(21),answ2(21),cabsoz(21),panswer(21),botswer(21)
      REAL albs(21),albp(21),botka(21)
      
      COMMON sza,vza,raa, am1,cos_vza,u1,u2,co,
     c  alam,reflec,height,aot,anna,pi,rv,tauaer,
     c  taumol,gaer,foto,t620,MOLEC
      COMMON /CD/ WELT
      COMMON /ANS/ ANSW2,Panswer
      COMMON/ANB/ albs,albp
      
      EXTERNAL fun,zbrent,func1,func2,func3
c******************************************************************
c           OLCI channels
            DATA wls/400., 412.5, 442.5, 490., 510., 560., 620.,
     c            665., 673.75, 681.25, 708.75, 753.75, 761.25,
     c            764.375, 767.5, 778.75, 865., 885., 900., 940.,
     c            1020./
        
        
c       imaginary part of ice refractive index at OLCI channels      
        DATA bai/6.27E-10,5.78E-10,6.49E-10,
     c 1.08E-9,1.46E-9,3.35E-09,    
     c 8.58E-09,1.78E-08,1.95E-08,2.1E-08,3.3E-08,6.23E-08,7.1E-08,
     c 7.68E-08,8.13E-08,9.88E-08,2.4E-07,3.64E-07,4.2E-07,5.53e-07,
     c       2.25E-06/

c       ozone vertical optical density  at OLCI channels      
        DATA cabsoz/
     c  1.3782e-4, 3.0488e-4, 1.6457e-3, 8.9359e-3, 1.7505e-2,
     c  4.3471e-2, 4.4871e-2, 2.1016e-2, 1.7162e-2, 1.4663e-2,
     c  7.9830e-3, 3.8797e-3, 2.9338e-3, 2.7992e-3, 2.7297e-3,
     c  3.2560e-3, 8.9569e-4, 5.1888e-4, 6.7158e-4, 3.1278e-4,
     c  1.4088e-5/

c     The array OLCI gains calib(ih) must be changed in case
c     gains are used        
c     No gains:        
        do 897 JH=1,21
           calib(JH)=1.
897     continue

c**********************************************        
c     input files:
c     input-1: OLCI data:        
         open(2022,file='input.dat')    
c     input-2: auxiliary data         
         open(5002, file='thv.dat')
c********************************************** 

c     output files:
c      main output before postprocessing:         
       open(1004, file='snow_parameters.dat')
       open(1001,file= 'spectral_spherical_albedo.dat')
       open(1002,file= 'spectral_plane_albedo.dat')
       open(1003,file= 'spectral_BOAR.dat')

c     main output after postprocessing
c     (high quality retrieval results):
       open(7004, file='output_snow_parameters_post.dat')
       open(7001,file= 'output_spectral_spherical_albedo_post.dat')
       open(7002, file='output_spectral_plane_albedo_post.dat') 
       open(7003, file='output_spectral_BOAR_post.dat')
c************************************************************        
     
c      constants:                     
       pi=acos(-1.)
          
c    relative vertical optical density of ozone f at all OLCI channels
c    (normalized to that at 620nm)
      
          abs620= 4.4871e-2
          do 802 mkls=1,21              
          f(mkls)=cabsoz(mkls)/abs620       
 802      continue
                        
c     aerosol properties(aot=aot500nm,ANNA=Angström exponent)
c     MOLEC not used anymore
c     number of pixels to be printed (with spectral data),
c     THV0=THV for the channel R(400nm),
c     THV1=THV for the measured and simulated spectra differences,
c     THV2=THV for the grain diameter,
c          THV3=THV for the total ozone column retrieved     
       read(5002,*) aot, ANNA, MOLEC,jend,THV0, THV1,THV2,THV3

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c***********************************************
c     imaginary part of ice refractive index at 865 and 1020nm      
          akap3=2.40e-7
          akap4=2.25e-6
       alt3=4.*pi*akap3/0.865
       alt4=4.*pi*akap4/1.02
          eps=          1./(1.-sqrt(alt3/alt4))

c     START of reading and processing (OLCI)

c         cycle for each OLCI pixels          
              do 87 j=1, 10 000 000
                   ind=0
                   NPOLE=0
c     reading OLCI data
         
        read(2022,*,IOSTAT=istat)
c     pixel number_x, pixelnumber_y, longitude,latitude,
c     solar zenith angle,soalr azimuthal angle,
c     viewing zenith angle, viewing azimuth angle, 
c     21 OLCI TOA reflectances (pi*I/cos(sza)/F0,
c     height of the surface (m),TOTAL OZONE load (ECMWF)     
     c  nsx,nsy,alon,alat,sza,saa,vza,vaa,
     c  (toa(ikss),ikss=1,21),height,ozon
        if (IS_IOSTAT_END(istat)) go to 87
        if(toa(1).lt.0.2)         go to 87
        raa=vaa-saa
        isnow=5
        
c       calculation of OLCI spectral indices - NDSI, NDBI,ratka      
        sk1=toa(17)+toa(21)
        sk2=toa(1) +toa(21)
        if (sk1.gt.0.)andsi=(toa(17)-toa(21))/sk1
        if (sk2.gt.0.)andbi=(toa(1)-toa(21))/sk2
        ratka=        toa(21)/toa(1)
        
c       ozone, transfer ECMWF OLCI total ozone to OLCI total ozone column in Dobson units:        
        akozon=ozon/2.1415e-5
        
c          calibration of OLCI channels(currently not used)      
                do 9011 MKL=1,21
          toa_cal(MKL)=toa(mkl)*calib(mkl)
          toa_cal_cor(MKL)= toa_cal(mkl)
9011            continue
         
c      snow fraction:
       factor=1.
                  cos_sza=cos(sza*pi/180.)
                  cos_vza=cos(vza*pi/180.)  
                  u1=0.6*cos_sza +1./3. +sqrt(cos_sza)/3.
                  u2=0.6*cos_vza +1./3. +sqrt(cos_vza)/3.
                  sin_sza=sin(sza*pi/180.)
                  sin_vza=sin(vza*pi/180.)   
                  cos_raa=cos(raa*pi/180.)
                  inv_cos_za=1./cos_sza+1./cos_vza
                    
c       cosine of scattering angle:                  
                  co = - cos_sza*cos_vza + sin_sza*sin_vza*cos_raa
                    isnow=1
c     isnow=1 - clean 100% snow cover
                   
c******************NOT USED for 100% snow cover:**************                    
c    the case of 100% snow cover (R(400nm)>0.65):                    
                    if (toa_cal(1).ge.THV0) go to 5021
                    isnow=3
c     the case of not 100% snow cover:              
c     scaling factor for patchy snow at 400nm                     
                    psi=RINFF(cos_sza,cos_vza)
c     factor=snow fraction ( SMALLER THAN 1.0):       
                           factor=toa_cal(1)/psi
c          snow TOA corrected for snow fraction                            
        do 2011 MKL=1,21
        toa_cal_cor(MKL)=toa_cal_cor(MKL)/factor
2011    continue
c***************************************************************+
5021      continue         

c     *****************   STEP 1   ******************                  
c     two - channel retrieval of EAL and R0 ( using 865 and 1020nm)
c     valid for clean and polluted PARTIALLY SNOW COVERED PIXELS
c     atmospheric contribution is ignored at 865 amd 1020nm          
c     nonabsorbing snow reflectance
:        
          r0=toa_cal_cor(17)**eps*toa_cal_cor(21)**(1.-eps)
		  
          write(*,*) r0,toa_cal_cor(1),toa_cal_cor(2),eps,toa(1),toa(2),factor,
     c      toa_cal(1),psi
          rv=r0
      rr2=toa_cal_cor(21)
	  
c     derivation of effective absorption length - dlina (mm):
      bal = alog(rr2/r0)*alog(rr2/r0)/alt4/(u1*u2/r0)**2.
      al = 1.e-3*bal  
	  
c     diameter of grains(mm):      
      diam=al/16.
	  
c     specific surface area(m*m/kg):      
      ssa= 104.7/al
c                 ************END of STEP 1**********

c     ***********************STEP 2************************
c     calculation of ozone gaseous transmittance at 620nm
c     snow reflectance at 620nm in absence of ozone absorption:
      reska_620=exp(-sqrt( 4.*pi * 8.58e-9/ (1.e-3 * 0.620) * al))
      rozone=r0*exp(-u1*u2/r0*sqrt( DARM* al))

c     atmospheric spectral characteristics at 620nm:      
      tauaer_620=aot*(0.62/0.5)**(-anna)
      taumol_620=0.0053/0.62**(4.0932)

c      total aerosol optical thickness at 620nm:      
      tau_620=tauaer_620+taumol_620
      
c      atmospheric reflectance, total transmittance and spherical albedo      
                              refatm_620=   func1(tau_620)
                              tatm_620=     func2(tau_620)
                              albatm_620=   func3(tau_620)
                              
c        gaseous transmittance at 620nm:
          TOAOZ=refatm_620+tatm_620*rozone/(1.-reska_620*albatm_620)
          T620=toa_cal_cor(7)/TOAOZ
          
c                 ************END of STEP 2**********

c     ***********************STEP 3************************
c                   *****the most important step*****
c                         step   3.1
c          snow spherical albedo determination
c     (with account for atmospheric scattering and snow impurities)
c*************************************************************************
          do 6363 nkl=1,21
c           wavelength(micron):             
            alam=wls(nkl)/1000.
            
c     OLCI CALBRATED TOA REFLECTANCE:            
      reflec= toa_cal_cor(nkl)
      tauaer=aot*(alam/0.5)**(-anna)
      taumol=0.0053/alam**(4.0932)
      tau=tauaer+taumol

c          asymmetry parameter
                       g0=0.5263
                       g1=0.4627
                       wave0=0.4685                       
                       gaer=g0+g1*exp(-alam/wave0)
                       g=tauaer*gaer/tau
					   
c     atmospheric reflectance, transmittance and albedo
                    refatm=func1(tau)
                    tatm=  func2(tau)
                    albatm=func3(tau)
                    z=u1*u2/r0

c     retrieval of  snow spherical albedo at all OLCI channels     
c     only retrievals at NKL=1 (400nm) and NKL =4 (490nm) are used to get
c     snow impurity parameters         

c     retrieval of  snow spectral spherical albedo at OLCI channel N=NKL:
c     lower, upper boundaries for snow spherical albedo retrieval:
      x1=0.001
      x2=1.0
c     accuracy of search routine:      
      tol=1.e-10    
      answer(nkl)=zbrent(fun,x1,x2,tol)
c     the mininimization routine has been finished
      if (answer(nkl).lt.0.0) answer(nkl)=1.
 6363 continue

c        Array answer(nkl) gives first guess
c      snow spherical albedo at all OLCI channels
c*************************************************************************
c     postprocessing -1
            if (answer(1).lt.0.0.and.answer(2).gt.0.0) answer(1)=1.0
            if (answer(1).lt.0.0.and.answer(3).gt.0.0) answer(1)=1.0
            if (answer(1).lt.0.0.and.answer(4).gt.0.0) answer(1)=1.0
       if (answer(2).lt.0.0.and.answer(3).gt.0.0) answer(2)=1.0
       if (answer(2).lt.0.0.and.answer(4).gt.0.0) answer(2)=1.0
       if (answer(3).lt.0.0.and.answer(2).gt.0.0) answer(3)=1.0
       if (answer(3).lt.0.0.and.answer(4).gt.0.0) answer(3)=1.0
c            end of postprocessing -1

c     retrieved spherical albedo at 2 wavelengths
c     retrieval first guess spherical/plane albedo
c     and bottom of atmosphere reflectance
c     at all channels for 100% and partially snow covered pixels:
       do 9025 JJ=1,21
          answ2(jj)  = factor*answer(jj)
          panswer(jj)= factor*(answer(jj)**u1)
          botswer(jj)= factor*r0*answer(jj)**(u1*u2/r0) 
          
 9025  continue
c           pixel classification ( 1- clean, 2-polluted,3-partailly snow covered)
       
      if (answ2(1).gt.0.98)                    isnow=1
      if (answ2(1).le.0.98.and.factor.gt.0.99) isnow=2
      if (answ2(1).le.0.98.and.factor.le.0.99) isnow=3

9002      format (i5,3e12.4,f8.4,40e12.4)

c******************************************
c            STEP_3_2                     *
c     retrieval of properties of          *
c     snow impurities + load of impurities*                    
c******************************************
      r1=answer(1)
      r2=answer(4)
      zara=1.
      if (r1.gt.0.999)  go to 8965
      if (r2.gt. 0.999) go to 8965
      zara=  (alog(r1)/alog(r2))**2.
8965  continue
      
c     1-retrieved absorption Angström exponent (AAE):      
      powe=alog(zara)/alog(490./400.)
      
      if (powe.le.0.9)   powe=0.
      if (powe.le.0.9)   go to 7754
c     2-retrieved pollution load coefficient (PLC), 1/mm:      
      polut=((0.4)**powe) * (alog(r1)*alog(r1))/al
c     special case of soot impurities:     
      if (powe.lt.1.2) go to 1941
c     DUST IMPURITIES:      
         NPOLE=1
      
c     3- retrieved effective diameter of dust grains:      
      DEFF=39.7373-11.8195*powe+0.8325*powe*powe
      
c     4- retrieved volumetric absorption coefficient of dust impurities
c     at the wavelength 1 micron      (1/mm)      
      absor1=10.916-2.0831*powe+0.5441*powe*powe
      
c     mass absorption coefficient (MAC) of dust in snow(    cm**3/g/mm) 
c      density of dust
      dens2=2.65 
      absor2=absor1/dens2
      
      densi=2.65/0.917
      ALOAD=1.8*densi*polut/absor1
      
c     5- retrieved impurity load (ppmw- ppm weight):      
      aload1=1.e+6*aload

c     6- retrieved mass absorption coefficient (MAC) of dust in snow at 1000nm(m**2/g)      
      absor1000=absor2*1.e-3
c     7-retrieved mass absorption coefficient (MAC) of dust in snow at 660nm(m**2/g) 
       absef660=absor1000*(660./1000.)**(-powe)

c     no retrieval for too low impurity load (below 2ppm):
       if (aload1.le.2.) deff=0.
       if (aload1.le.2.) absor1=0.
       if (aload1.le.2.) absef660=0.
       if (aload1.le.2.) absor1000=0.
          if (aload1.le.2.)  powe  =0.
          if (aload1.le.2.)   polut =0.
          if (aload1.le.2.) aload1=0.
        
                   go to 1942
c***********************************************      
c     special case of soot  pollution:    
 1941              continue
      NPOLE=2            
      aload1=0.
      deff=0.
      absor1=0.
      absef660=0.
      absor1000=0.
      if (powe.lt.0.0) isnow=1
      if (powe.lt.0.0) polut=0
      if (powe.lt.0.0) powe=0
      if (powe.lt.0.0) go to 1942
      
      aco=2.06e+3
      aload1=polut/aco
       if (factor.lt.0.99) isnow=3
       
       if (aload1.le.2.) powe=0.
       if (aload1.le.2.) polut=0.
       if (aload1.le.2.) aload1=0.
781    continue
1942   continue      
7754   continue

c     FINISH:
c     PROPERTIES OF SNOW POLLUTANTS (IF ANY) ARE
c     RETRIEVED

c         *****************Step 3.3******************
c*********************************************************
c       calculation of snow spectral albedo( plane/spherical) and
c       bottom of atmosphere snow reflectance
       
        DO 765 IK=1,21 
         ALFA = 4.*pi*bai(ik)/ (wls(ik)*1.e-6)
         SDU = polut*(wls(ik)/1000.)**(-powe)
         albs(ik)=factor*exp(-sqrt(ALFA+SDU)*al))
         albp(ik)=albs(ik)**u1
         botka(ik)=factor*r0*exp(-sqrt(ALFA+SDU)*al))** (u1*u2/r0)  
765     continue

c     *************Step 3.4***********************        
c                   BBA calculation      
                if (isnow.eq.1) go to 8998
    
c                   3.4.1. Polluted snow           
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       call bbalb(2,rvis,rnir,rsw)
       call bbalb(1,rviss,rnirs,rsws)
c*************end of calculation of BBA for polluted snow

       go to 877   
8998   continue      
c                  3.4.2: clean snow
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      rsw=0.5271+0.3612* exp(-u1*sqrt(0.02350*al))
      rvis=exp(-u1*sqrt (7.86e-5*al))
      rnir=0.2335+0.56*exp(-u1*sqrt(0.0327*al))
      rsws=0.5271+0.3612* exp(-sqrt(0.02350*al))
      rviss=exp(-sqrt (7.86e-5*al))
      rnirs=0.2335+0.56*exp(-sqrt(0.0327*al))
      if (toa_cal_cor(1).lt.THV0)rsw=rsw*factor
      if (toa_cal_cor(1).lt.THV0)rvis=rvis*factor
      if (toa_cal_cor(1).lt.THV0)rnir=rnir*factor
      if (toa_cal_cor(1).lt.THV0)rsws=rsws*factor
      if (toa_cal_cor(1).lt.THV0)rviss=rviss*factor
      if (toa_cal_cor(1).lt.THV0)rnirs=rnirs*factor
877        continue

c    ***** end of calculation for BBA for clean snow *****
c                         END of BBA calculation
c                         END OF RETRIEVAL
c                         FINISH: all parameters are retrieved

c     **********************STEP 4*******************
c         *****  OLCI SPECTRAL TOA modelling  *****
c        calculation of gaseous transmittance at 620, 940, and 761nm:             
         tt620 = toa_cal(7)  /botka(7)
         tt940 = toa_cal(20) /botka(20)
         tt761 = toa_cal(13)/botka(13)

c        calculation of gaseous vertical optical depth:         
          vodox =-alog(tt761)/inv_cos_za
          vodvod=-alog(tt940)/inv_cos_za
          vodka=-alog(tt620)/inv_cos_za

c     calculation of TOA reflectance
                    tauaer=aot*(0.62/0.5)**(-anna)
                    taumol=0.0053/0.62**(4.0932)
                    tau=tauaer+taumol
                        
c          asymmetry parameter
                       g0=0.5263
                       g1=0.4627
                       wave0=0.4685                       
                    
                    refatm=func1(tau)
                    tatm=  func2(tau)
                    albatm=func3(tau)
                    rspher=answ2(7)
                    rsnow=r0*rspher**(u1*u2/r0)
                    ckk=refatm+tatm*rsnow/(1.-rspher*albatm)
                    tt620=toa_cal_cor(7)/ckk

         vodoz=-alog(tt620)/inv_cos_za
         tocos=vodka*9349.3
         
             do 1972 jt=1,21
                alam=wls(jt)/1000.
                rspher=answ2(jt)
                welt=factor
                TOAR=fox(rspher)
                 tauaer=aot*(alam/0.5)**(-anna)
                     taumol=0.0053/alam**(4.0932)
                     tau=tauaer+taumol
                        
c          asymmetry parameter
                       gaer=g0+g1*exp(-alam/wave0)
                       g=tauaer*gaer/tau

                    refatm=func1(tau)
                    tatm=  func2(tau)
                    albatm=func3(tau)
                    rsnow=r0*rspher**andre
         Tozone=t620**f(jt)
         
         TOX=1.
         TVODA=1.
         
         if (jt.eq.13) TOX= tt761**1.
         if (jt.eq.14) TOX= tt761**0.532
         if (jt.eq.15) TOX= tt761**0.074  
         if (jt.eq.19) TVODA=tt940**0.25
         if (jt.eq.20) TVODA=tt940**1.
         
         TOK(jt)=TOAR*TVODA*TOX*tozone
         
         if (toa_cal(1).lt.THV0) TOK(jt)=TOAR*TVODA*TOX*tozone*factor
1972   continue
1973   continue
91     continue
191    continue

c     spectral TOA and BOA refletances with account for gaseous absorption
c     *******************has been derived ****************************  
c************************end of step 4************************************

c                                         STEP 5
c                                       FINAL STEP
c                                        study of
c          quality of retrievals and print of retrieved snow characteristics
c                                     postprocessing         

         sum1=0
         sum2=0
         
         do 2022 MNK=1,21
            sum2=sum2+toa_cal_cor(mnk)
         sum1=sum1+ (toa_cal_cor(mnk)-tok(mnk))* (toa_cal_cor(mnk)-tok(mnk))
 2022    continue
         sum1=sqrt(sum1/21.)
         sum2=sum2/21.
         cv1=100.*sum1/sum2
         um1=0.
         um2=0.
         do 202 MNK=1,21
            if (MNK.eq.13) go to 202
            if (MNK.eq.14) go to 202
            if (MNK.eq.15) go to 202
            if (MNK.eq.19) go to 202
            if (MNK.eq.20) go to 202
            um2=um2+toa_cal_cor(mnk)
         um1=um1+ (toa_cal_cor(mnk)-tok(mnk))* (toa_cal_cor(mnk)-tok(mnk))
 202    continue
        um1=sqrt(um1/16.)
        um2=um2/16.
        cv2=100.*um1/um2
        difoz=0.0
        icloud=0
        if (tocos.gt.0.)difoz=100.*(tocos-akozon)/tocos
        if (abs(difoz).gt.10.) icloud=1
        DIFKA=abs(difoz)
        ccv2=abs(cv2)
        asinka=inka
        
c                     **OUTPUT**
 3092          format(i5,2x,2e12.4,i3,f8.4,20(2x,f8.4))
          if (powe.lt.0.9)  aload1 =0.0
          if (powe.gt.10.)  powe=   0.0
          if (powe.gt.10.)  aload1= 0.0
          if (powe.lt.0.9)  powe=   0.0
22022    continue
        DIFKA=abs(difoz)
		
        NBARE=0
        NSNOW=0
c       dark ice        
        if (ANDBI.lt.0.65.and. toa(1).lt.0.75)  NBARE=2
c       clean        
        if (ANDSI.gt.0.33.and. NBARE.NE.2)      NBARE=1
c          SNOW INDEX            
        if (toa(1).gt.0.75.and.andsi.lt.0.1)   NSNOW=1
        polut=1.e+3*polut
		
c     *****    output: no postprocessing *********   
        write(1004,*)
c       pixel number, latitude, longitude,
     c               j,alat,alon,
        
c     class ( 1- clean, 2 - polluted, 3-partially snow -covered),
c     cloud fraction, diameter of grains, SSA, EAL,R0,        
     c       isnow,factor,diam,ssa,al,rv,
        
c     impurity load (ppm_weight), Angstrom Absorption Exponent,
c     normalized volumetric absorption coefficient at 1000nm,        
     c       aload1,powe,polut,
c     effective radius of dust grains, volumetric absorption coefficients at 1000nm,
c     dust mass absorption coefficient at 660 and 1000nm,        
     c       deff,absor1,absef660,absor1000,
c      plane visible and NIR BBA        
     c  rsw, rvis,rnir,
c      spherical  visible and NIR BBA          
     c  rsws, rviss,rnirs,
c      OLCI spectral indices (NDBI,NDSI,OSI)
     c  andbi,andsi,ratka,
c      type of pollutants  ( 1-dust, 2- soot)     
     c      NPOLE,
c     bare ice index (0-no bare ice, 1-bare ice - clean, 2-bare ice - polluted)
     c      NBARE,
c     snow index     (0-no snow, 1- snow)
     c       NSNOW,
c     SZA,VZA,RAA, reflectance at channel 1, reflectance at channel 2   
     c sza,vza,raa,toa(1),toa(21),
c    retrieved TOC,ECMWF, TOC difference(%),cv1,cv2
     c       tocos,akozon,difka,cv1,cv2
	 
        if (j.gt.jend) go to 2403
        write(1001,*) j,alat,alon,   
     c      isnow,factor,(albs(ir),ir=1,21)
        
        write(1002,*) j,alat,alon,   
     c       isnow,factor,(albp(ir),ir=1,21)
          
        write(1003,*) j,alat,alon,   
     c       isnow,factor,(botka(ir),ir=1,21)
 2403   continue
c**************************************************
c     if (DIFKA.gt.12.0.or. cv2.gt.5.0.or.diam.ge.0.1) go to 932
        if (ccv2.gt.THV1.or.diam.lt.THV2.or.DIFKA.gt.THV3) go to 932
        
c     *****    output: postprocessing *****
      write(7004,*) j,alat,alon,   
     c isnow,factor,diam,ssa,al,rv,
     c aload1,powe,polut,
     c deff,absor1,absef660,absor1000,
     c rsw,rvis,rnir,rsws,rviss,rnirs,ansi,andbi,
     c       ratka,npole,nbare,nsnow,
     c sza,vza,raa,toa(1),toa(21),
     c tocos,akozon,difka,cv1,cv2 
         if (j.gt.jend) go to 2404
        write(7001,*) j,alat,alon,   
     c      isnow,factor,(albs(ir),ir=1,21)  
        write(7002,*) j,alat,alon,   
     c       isnow,factor,(albp(ir),ir=1,21)     
        write(7003,*) j,alat,alon,   
     c       isnow,factor,(botka(ir),ir=1,21)
2404    continue      
932       continue
87        continue
      STOP
      END

c****************************************************
c               MINIMIZATION FUNCTION:      
      FUNCTION  fun(x)
       common sza,vza,raa, cos_sza,cos_vza,u1,
     c  u2,co,alam,reflec,height,aot,anna,pi,rv,tauaer,
     c taumol,gaer,foto
                        r0=rv                                    
                    RFINAL=refatm+tatm*r0*x**z/(1.-albatm*x)
                    fun=reflec-RFINAL
           write(*, *) x,fun,reflec,rfinal,refatm,tatm,r0,x,z,albatm
c**********************************                    
                                             return  
                                             end
											
             function func1(tau)
        common sza,vza,raa, cos_sza,cos_vza,u1,u2,co,
     c  alam,reflec,height,aot,anna,pi,r0,tauaer,
     c     taumol,gaer,foto
        
       wav=alam
c****************************************************
c                  HG phase function for aerosol
                       g11= 0.80
                       g22=-0.45
                       pa1=(1-g11*g11)/(1.-2.*g11*co+g11*g11)**1.5
                       pa2=(1-g22*g22)/(1.-2.*g22*co+g22*g22)**1.5
                       cp=(gaer-g11)/(g11-g22)
                       pa=cp*pa1+(1.-cp)*pa2
                        pr=0.75*(1.+co*co)
                        p=(taumol*pr+tauaer*pa)/tau
                        g=tauaer*gaer/tau
                       
c******************************************************************
c                                  SOBOLEV
                       ASTRA=(1.-exp(-tau*inv_cos_za))/(cos_sza+cos_vza)/4.
                       OSKAR=4.+3.*(1.-g)*tau
                       
                       b1=1.+1.5*cos_sza+(1.-1.5*cos_sza)*exp(-tau/cos_sza)
                       b2=1.+1.5*cos_vza+(1.-1.5*cos_vza)*exp(-tau/cos_vza)
                       
               RSS=p*ASTRA
               RMS=1.-b1*b2/OSKAR+(3.*(1.+g)*cos_sza*cos_vza-2.*(cos_sza+cos_vza))*ASTRA
               
               Ratm=RSS+RMS
c******************************************************************
                                              func1=ratm
                                              return  
                                              end

      function func2(tau)
       common sza,vza,raa, cos_sza,cos_vza,u1,u2,co,
     c alam,reflec,height,aot,anna,pi,rv,tauaer,
     c taumol,gaer,foto
               TZ=1.+gaer*gaer+(1.-gaer*gaer)*sqrt(1.+gaer*gaer)              
               Baer=0.5 +gaer*(gaer*gaer-3.)/TZ/2.             
               if (gaer.lt.1.e-3) go to 89 
               GT=(1.+gaer)/  sqrt ( 1.+gaer*gaer)  -1.
               Baer=(1.-gaer)*GT/2./gaer
 89            continue
               B=(0.5*taumol+Baer*tauaer)/tau             
               t1=exp(-B*tau/cos_sza)
               t2=exp(-B*tau/cos_vza)       
               func2=t1*t2 
         RETURN
         END

       function func3(tau)
       common sza,vza,raa, cos_sza,cos_vza,u1,u2,co,
     c alam,reflec,height,aot,anna,pi,rv,tauaer,
     c      taumol,gaer ,foto
         gasa=0.5772157
         y=(1.+tau)*tau*exp(-tau)/4.
         Z3=tau*tau*(-alog(tau)-gasa)
         Z4=tau*tau*(tau-tau*tau/4.+tau*tau*tau/18.)
         Z=z3+z4
         f=(1.+0.5*tau)*Z/2.-y
         W1=1.+f
         W2=1.+0.75*tau*(1.-g)
         func3=1.-W1/W2
       return
       end

        FUNCTION zbrent(fun,x1,x2,tol)
      INTEGER ITMAX
      REAL zbrent,tol,x1,x2,func,EPS
      EXTERNAL fun,func1,func2,func3
      PARAMETER (ITMAX=100,EPS=3.e-8)
      INTEGER iter
      REAL a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      fa=fun(a)
      fb=fun(b)
      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))then
	    zbrent=-999
		RETURN
		endif
c    *'root must be bracketed for zbrent'
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=fun(b)
11    continue
c      pause 'zbrent exceeding maximum iterations'
      zbrent=b
      
      return
      END


       FUNCTION RINFF(cos_sza,cos_vza)
       a=1.247
       b=1.186
       c=5.157
       p=0.0
       rinff=( a+b*(cos_sza+cos_vza)+c*cos_sza*cos_vza+p)/4./ (cos_sza+cos_vza)
       return
       end

      
           FUNCTION  fox(x)
       common sza,vza,raa, cos_sza,cos_vza,u1,u2,co,
     c alam,reflec,height,aot,anna,pi,rv,tauaer,
     c taumol,gaer,foto
          common /CD/ WELT                         
c****************************************************

                  r0=rv
                    tauaer=aot*(alam/0.5)**(-anna)
                     taumol=0.0053/alam**(4.0932)
                     tau=tauaer+taumol
                        
c          asymmetry parameter
                       g0=0.5263
                       g1=0.4627
                       wave0=0.4685                       
                       gaer=g0+g1*exp(-alam/wave0)
                       g=tauaer*gaer/tau
                    
                    refatm=func1(tau)
                    tatm=  func2(tau)
                    albatm=func3(tau)
                    z=u1*u2/r0
                    FOX=1.
                    vv=x/welt
          if (x.gt.0.05)FOX=refatm+tatm*welt*r0*vv**z/(1.-albatm*x*welt)
c          write(*,*) fox,refatm,tatm,welt,r0,vv,z
                                             return  
                                             end

      
      subroutine BBALB (jdom,rp1,rp2,rp3)
      real answ2(21),panswer(21),albs(21),albp(21)
      common /ans/answ2,pa
      common /bns/albs,albp

      alam2=0.4
      alam3=0.56
      alam5=0.709
      alam6=0.753
      alam7=0.865
      alam8=1.02
      
c**************************************************      
c     sphrerical BBA      
      if (jdom.eq.1) r2=answ2(1)
       if (jdom.eq.1)r3=answ2(6)
       if (jdom.eq.1)r5=answ2(11)
       if (jdom.eq.1)r6=answ2(12)
      if (jdom.eq.1) r7=answ2(17)
       if (jdom.eq.1)r8=answ2(21)
	   
c         plane bba      
      if (jdom.eq.2) r2=panswer(1)
       if (jdom.eq.2)r3=panswer(6)
       if (jdom.eq.2)r5=panswer(11)
       if (jdom.eq.2)r6=panswer(12)
       if (jdom.eq.2)r7=panswer(17)
       if (jdom.eq.2)r8=panswer(21)

c     spherical BBA    - clean  
       if (jdom.eq.3)r2=albs(1)
       if (jdom.eq.3)r3=albs(6)
       if (jdom.eq.3)r5=albs(11)
       if (jdom.eq.3)r6=albs(12)
       if (jdom.eq.3)r7=albs(17)
      if (jdom.eq.3)r8=albs(21)
	  
c         plane bba     -clean 
       if (jdom.eq.4)r2=albp(1)
       if (jdom.eq.4)r3=albp(6)
       if (jdom.eq.4) r5=albp(11)
       if (jdom.eq.4)r6=albp(12)
       if (jdom.eq.4)r7=albp(17)
       if (jdom.eq.4)r8=albp(21)
c      ------------------------------
c     QUADRATIC POLYNOMIAL for the range 400-709nm
      x0=alam2
      x1=alam3
      x2=alam5
      y0=r2
      y1=r3
      y2=r5
                          d1=(x0-x1)*(x0-x2)
                          d2=(x1-x0)*(x1-x2)
                          d3=(x2-x0)*(x2-x1)
                           a1=x1*x2*y0/d1+x0*x2*y1/d2+x0*x1*y2/d3
                           b1=-(x1+x2)*y0/d1-(x0+x2)*y1/d2-(x0+x1)*y2/d3
                           c1=y0/d1+y1/d2+y2/d3
						   
c     QUADRATIC POLYNOMIAL for the range 709-865nm
      x0=alam5
      x1=alam6
      x2=alam7
      y0=r5
      y1=r6
      y2=r7
                          d1=(x0-x1)*(x0-x2)
                          d2=(x1-x0)*(x1-x2)
                          d3=(x2-x0)*(x2-x1)
                           a2=x1*x2*y0/d1+x0*x2*y1/d2+x0*x1*y2/d3
                           b2=-(x1+x2)*y0/d1-(x0+x2)*y1/d2-(x0+x1)*y2/d3
                           c2=y0/d1+y1/d2+y2/d3

c     exponential approximation for the range above 0.865nm
                           rati=r7/r8
                           alasta=(alam8-alam7)/alog(rati)
                           an=1./alasta
                           p=r7*exp(alam7/alasta)
                       
c     END of approximations for 3 intervals

c     approximation for the solar flux
c     f=f0+f1*exp(-bet*lambda/star1)+f2*exp(-gam*lambda)
                             f0=32.38
                             f1=-160140.33
                             f2=7959.53
                             bet= 1./0.08534
                             gam=1./0.40179

c     ANALYTICAL INTEGRATION OF SOLAR FLUX (DENOMINATOR)
                             z1=0.3
                             z2=0.7
                           sol1a=f0*(z2-z1)
                           sol1b=-f1*(  exp(-bet*z2) -exp(-bet*z1)) /bet
                           sol1c=-f2*(  exp(-gam*z2) -exp(-gam*z1))/gam
                           sol1=sol1a+sol1b+sol1c
						   
                              z1=0.7
                              z2=2.4
                           sol1a=f0*(z2-z1)
                           sol1b=-f1*(  exp(-bet*z2) -exp(-bet*z1)) /bet
                           sol1c=-f2*(  exp(-gam*z2) -exp(-gam*z1)) /gam
                     sol2=sol1a+sol1b+sol1c
                     sol3=sol1+sol2
c*************************************************                     
                              z1=0.7
                              z2=0.865
                                   asol1a=f0*(z2-z1)
                     asol1b=-f1*(  exp(-bet*z2) -exp(-bet*z1)) /bet
                    asol1c=-f2*(  exp(-gam*z2)-exp(-gam*z1))/gam
                     asol=asol1a+asol1b+asol1c
c*******************************************************      
c     ANALYTICAL EQUATION FOR THE NOMINATOR
c     integration over 3 segments
c     segment 1                     
                     z1=0.3
                     z2=0.7
                     ajx1=a1*sol1
c     first coef                 
                     ak1=(z2**2.-z1**2.)/2.
                     ak2=(z2/bet+1./bet/bet)*exp(-bet*z2)-
     c                    (z1/bet+1./bet/bet)*exp(-bet*z1)
                     ak3=(z2/gam+1./gam/gam)*exp(-gam*z2)-
     c                    (z1/gam+1./gam/gam)*exp(-gam*z1)
                     ajx2=b1*(f0*ak1  -f1*ak2  -f2*ak3 )
c     second coef                 
                     am1=(z2**3.-z1**3.)/3.
                     
      am2=(z2**2./bet+2.*z2/bet/bet+2./bet/bet/bet)*exp(-bet*z2)-
     c         (z1**2./bet+2.*z1/bet/bet+2./bet/bet/bet)*exp(-bet*z1)
            
        am3=(z2**2./gam+2.*z2/gam/gam+2./gam**3.)*exp(-gam*z2)-
     c           (z1**2./gam+2.*z1/gam/gam+2./gam**3.)*exp(-gam*z1)
    
             ajx3=c1*(f0*am1 -f1*am1 -f2*am2)
             aj1=ajx1+ajx2+ajx3

c     segment 2                   
                      z1=0.7
                      z2=0.865
                         cajx1=a2*asol
c     first coef                 
                     ak1=(z2**2.-z1**2.)/2.
                     ak2=(z2/bet+1./bet/bet)*exp(-bet*z2)-
     c                    (z1/bet+1./bet/bet)*exp(-bet*z1)
                     ak3=(z2/gam+1./gam/gam)*exp(-gam*z2)-
     c                    (z1/gam+1./gam/gam)*exp(-gam*z1)
                     cajx2=b2*(f0*ak1  -f1*ak2  -f2*ak3 )
c     second coef                 
                     am1=(z2**3.-z1**3.)/3.
      am2=(z2**2./bet+2.*z2/bet/bet+2./bet/bet/bet)*exp(-bet*z2)-
     c         (z1**2./bet+2.*z1/bet/bet+2./bet/bet/bet)*exp(-bet*z1)
            
        am3=(z2**2./gam+2.*z2/gam/gam+2./gam**3.)*exp(-gam*z2)-
     c           (z1**2./gam+2.*z1/gam/gam+2./gam**3.)*exp(-gam*z1)
             cajx3=c2*(f0*am1 -f1*am2 -f2*am3)
                  caj2=cajx1+cajx2+cajx3
c     segment 3
                      z1=0.865
                      z2=2.4

             aj31=(1./an)*(exp(-an*z2)-exp(-an*z1))
             aj32=(1./(bet+an))*(exp(-(bet+an)*z2)-exp(-(an+bet)*z1))
             aj33=(1./(gam+an))*(exp(-(gam+an)*z2)-exp(-(an+gam)*z1))
               caj3=(-f0*aj31-f1*aj32-f2*aj33)*p
c*************************************************************************
                     ajto=aj1+caj2+caj3
                     bbavis=aj1/sol1
                     bbanir=(caj2+caj3)/sol2
                     bbat=ajto/sol3
                     
                    rp1=bbavis
                    rp2=bbanir
                    rp3=bbat
 892                 CONTINUE
                
                     return
                     end
