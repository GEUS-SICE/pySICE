           program SICE
      
c                      VERSION 3.1 or 2.0.2
      
c                   AUGUST 30, 2019
      
c     This code retrieves snow/ice  albedo
c               and related snow products
c     for clean Arctic atmosphere



c***********************************************      

c     The errors increase with the load of
c     pollutants in air
c***********************************************
 
c***************************************************     
c     1rst version:        21.05.2019
c     correction: Alex: 09.06.2019
c     correction: Alex: 14.06.2019
c     modification: Alex: 28.06.2019 - ozone retrieval is performed
c     modification: Alex: 07.08.2019-if OLCI reflectance at channel 21 is below 0.5, then the clean snow reflectance is not retrieved but assumed.
c***************************************************

c     modification
c     AUGUST 29 2019
c               new code for BBA calculation> bare ice/dark ice      
c     modification
c     AUGUST 30 2019
c               reading new input file: 2 new column added: total ozone, water vapour

      
c                      Alexander  KOKHANOVSKY
c     a.kokhanovsky@vitrocisetbelgium.com

c     input:      spectral OLCI TOA reflectance (R=pi*I_reflec/cos(SZA)/E_0),
c                    satellite and solar angles, lat, lon, height_of_ground_in_m

c     output:    albedo, snow grain size, snow specific surface area,
c                    concentration and type/properties  of pollutants in snow, etc.
      
           real           answer(21),   tozon(21), voda(21),
     c                      w(21),            toa(21),
     c                      refl(21),rp(21),bai(21),xa(168),ya(168)
               integer ndate(6) 
                external fun,sobthv,psi,funp,funs
                
                common sza,vza,raa, g, wave,rtoa ,r0, height,aot
                
                common /quatro/x0,x1,x2,y0,y1,y2,ak1
                
                common /qua/ xa,ya,am1,AL,NSOLO
                
c             OLCI channels
                
         DATA w/  0.4000E+00,      
     c                  0.4125E+00,      
     c                  0.4425E+00,     
     c                  0.4900E+00,      
     c                  0.5100E+00,      
     c                  0.5600E+00,      
     c                  0.6200E+00,      
     c                  0.6650E+00,      
     c                  0.6737E+00,      
     c                  0.6812E+00,      
     c                  0.7088E+00,      
     c                  0.7538E+00,      
     c                  0.7613E+00,      
     c                  0.7644E+00,      
     c                  0.7675E+00,      
     c                  0.7788E+00,      
     c                  0.8650E+00,      
     c                  0.8850E+00,      
     c                  0.9000E+00,      
     c                  0.9400E+00,      
     c                  0.1020E+01     /



         
c     Imaginary part of ice refractive index at OLCI channels
         
                     DATA bai/         
     c                  2.365E-11,      
     c                  2.7E-11,     
     c                  7.0E-11,      
     c                  4.17E-10,      
     c                  8.04E-10,      
     c                  2.84E-09,      
     c                  8.58E-09,      
     c                  1.78E-08,      
     c                  1.95E-08,      
     c                  2.1E-08,      
     c                  3.3E-08,      
     c                  6.23E-08,      
     c                  7.1E-08,      
     c                  7.68E-08,      
     c                  8.13E-08,      
     c                  9.88E-08,      
     c                  2.4E-07,      
     c                  3.64E-07,      
     c                  4.2E-07,
     c                  5.53e-07,                
     c                  2.25E-06     /
                     
    
c     input
c**************************************************************************         
c                     ns,alat,alon,sza,vza,raa,height,(toa(iks),iks=1,21)
c     ns-pixel number ( could be substituted bz measurement time, etc.)
c     latitude, longitude, solar zenith angle, solar azimuthal angle,
c     viewing zenith angle, viewing azimuthal angle, height of underlying surface(meters),
c     OLCI TOA reflectance at 21 channels
                     
c                     open(1,file='olci_toa.dat')
c                     open(1,file='/home/vtcb/snow/sice1mirror/proc_input
c     c/v2_2_test/output/max_egp_olci.dat')
                     
c     number of lines to read
                     
                     open(1,file='olci_toa_newformat.dat')

                     
c     number of lines to read                    
                     open(101,file='nlines.dat')

c                     assumed date
                     ndate(3)=010120
                     
c     assumed AOT at 500nm for THV clean and polluted snow
                            AOT=thv                   
                            open(102,file='thv.dat')
                            
c     assumed THVs for the limiting value of the diameter of grains
c      and for the TOA reflectance at 1020nm                            
                            open( 1984, file='limits.dat')
                            read(1984,*) ALR21,ALDI
                            
c     imaginary part of ice refractive index
                     open(5000,file='ice_index.dat')
                     open(2099,file='interm.dat')
                     open(1212, file='interm_alb.dat')
                     
c***************************************************************************
c                  Alex 09.06.2019
c     ozone vertical absorption thickness
                     
                     open(1975,file='tg_vod.dat')
                         open(1985,file='tg_water_vod.dat')
c     ozone concentration for a given place in kg/m/m
c               it is given now in the file olci_toa_newformat.dat                     
c                     open(1973,file='ozone.dat')
                
c                retrieved ozone
                       open(1914,file='retrieved_O3.dat')


                     
                     read(1975,*)
                     read(1975,*)
                     do 1976 np=1,21
                        read(1975,*)dlin,tozon(np)
                        read(1985,*)dlin,voda(np)
1976            continue

                     
c     output
                 
                 open(157,  file= 'spherical_albedo.dat' )
                 open(257,  file= 'planar_albedo.dat'     )
                 open(357,  file= 'boar.dat'                    )
                 open(1001,file= 'size.dat'                     )
                 open(457,file=   'impurity.dat'              )
                 open(557,file=   'bba.dat'                     )
                 open(5570,file=   'bba_alex_reduced.dat')
                         open(55701,file=   'notsnow.dat')
c               number of lines to be processed                                
                 read( 101,*) nlines

                 
c                pre-assumed number: AOT-aerosol optical thickness  at 500nm       
                 read(102,*) AOT

c                reading ice refractive index
                  do 8867 jg=1,168  
             read(5000,*) xa(jg),an,ya(jg)
 8867     continue


          icloud=0
          iice=0
          scale=acos(-1.)/180.





          
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
c              -START-           OF MAIN ROUTINE:
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
                 do 87 J=1,nlines
c                                           READING OLCI DATA:           
          read(1,*,end=87) ns,alat,alon,sza,vza,saa,vaa,height,
     c                   (toa(iks),iks=1,21),OZON,WATER

             
       
          
c     ecmwf ozone in DOBSON UNITS from OLCI file:
          
     
              totadu=46729.*ozon
              AKOEF=totadu/404.59
c                   kg/m**2. transfer to mol/cm**2         
           roznov=2.99236e-22          
          vap=water/roznov
          
      
          
          AKOWAT=vap/3.847e+22
          
c!!!!Alex 28.06.2019         
     
     
        
c     AMF:  Alex
          
          amf=1./cos(sza*scale)+1./cos(vza*scale)
          
          eps=1.55
          
                 BX=(toa(21)**(1.-eps))  * (toa(17)**eps) / toa(7)
                 BXX=alog(BX)
                 BXXX=BXX/1.11e-4/amf
                 if (BXXX.GT.500.) BXXX=999.
                 deltak=100.*(1.-BXXX/totadu)
               write(1914,*) ns,alat,alon,BXXX,totadu,deltak,sza,vza,amf
                 
               do 1941 nv=1,21
          
                  tvoda= exp(amf*voda(nv)*AKOWAT)
                 
             toa(nv)=toa(nv)*tvoda*exp(amf*tozon(nv)*AKOEF)

                  
1941        continue
         


            
c the case of ice and open water - no retrieval
             if(toa(21).lt.ALR21) iice=1         
             if(toa(21).lt.ALR21) go to 8797
             
c     transfer of OLCI relative azimuthal angle to the definition used in
c     radiative transfer code                    
                    raa=180.-(vaa-saa)
                    
c****************************************************                    
                               pi=acos(-1.)
c****************************************************
                    
                    am1=cos(sza*pi/180.)
                    am2=cos(vza*pi/180.)
                    
                    ak1=3.*(1.+2.*am1)/7.
                    ak2=3.*(1.+2.*am2)/7.
                    
c     STEP 1
c     calculation of NDSI, NDBI and flags
                        rr1=toa(17)   
                        rr2=toa(21)
                        
c****************************************************************************
                        arr1=toa(1)
c****************************************************************************
c     derivation of indices
                        
                        indexs=0
                        indexi=0
                        indexd=0
                       
                        
                        andsi=(rr1-rr2)/(rr1+rr2)
                        andbi=(arr1-rr2)/(arr1+rr2)

                          if (andsi.gt.0.03.and.arr1.gt.0.5)indexs=1
                          if (andbi.gt.0.33)    indexi=1
                          if (andbi.gt.0.66)    indexd=1

                          
c STEP 2                         
c     retrieval of snow properties ( R_0, size of grains from OLCI channels 865(17) and 1020nm(21)
c     assumed not influenced by atmospheric scattering and absorption processes)                       
                        akap1=2.40e-7
                        akap2=2.25e-6
                        
                        
                        alpha1=4.*pi*akap1/0.865
                        alpha2=4.*pi*akap2/1.020
                        
                         b=sqrt(alpha1/alpha2)
                         eps=1./(1.-b)
                    
 
c********************************************************
c                      reflectivity of nonabsorbing snow layer                       
                        r0=(rr1**eps)*(rr2**(1.-eps))
                        
                              xxx=ak1*ak2/r0                
c*********************************************************                    
                       
c                  effective absorption length(mm)
                    zz=(rr2/r0)
                    bal=alog(zz)*alog(zz)/xxx/xxx/alpha2
                    al=bal/1000.
c                  effective grain size(mm):diameter
                    D=al/16.36
                    
                    if (d.lt.ALDI) icloud=1
                    if (d.lt.ALDI) go to 8797
c                  snow specific area ( dimension: m*m/kg)
                    densit=0.917
                    area=   6./D/densit
                    

c*********************************************************

c            STEP 3
c     checking for type of snow: clean or polluted
c     polluted snow - small TOA(1)
c      atmosphere with AOT=0.8 at 400nm is assumed
c     AOT=0.6 ( see thv.dat)
c                    Alex 09.06.2019                    
c                    thv=0.0
                    THV=SOBTHV(AOT)
                 
                           isnow=1
                        write(2099,*) ndate(3),toa(1),thv
                       if ( toa(1).lt.THV) go to 1961
c**********************************************************
 9393                       continue
c     STEP 4a
                       
c     1. clean snow retrieval
                       isnow=0
                       bf=0.
                       bm=0.
                       ntype=0
                       conc=0.
                   
c**********************************************************
                 do 7777 isk=1,21
                          
                       wave=w(isk)           
                       alpha=4.*pi*bai(isk)/wave
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c                                          NB                       
c     the spherical albedo derivation: answer
                       ABSOR=1000.*alpha*AL
                      if (absor.gt. 1.e-6) answer(isk)=exp(-sqrt(ABSOR))
                      if (absor.le. 1.e-6) answer(isk)=1.0
                      
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      
c                 derivation of plane albedo                  
                      rp(isk)=answer(isk)**ak1
                      
c                  derivation of snow reflectance function                      
                      refl(isk)=r0*answer(isk)**xxx
                   
7777       continue


                       
c*************end of clean snow retrieval*******************
                                     GO TO 1964
c**********************************************************

                       
 1961                                CONTINUE

c                            STEP 4b
                                     
c                       2. polluted snow retrieval                      

c                       it is assumed that albedo is in the range 0.1-1.0           
                    
                       x1=0.1
                       x2=1.2
                       
c                     convergence parameter for albedo search:           
                       tol=1.e-6
           
           
                  do 577 isk=1,21
                          
                       wave=w(isk)           
                       rtoa=toa(isk)
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c                                          NB                       
c     the solution of transcendent equation to find spherical albedo: answer

                       
c************************************************************************                       
c     Alex AUGUST 7, 2019
                       az=1.247
                       bz=1.186
                       cz=5.157
                       am11=sqrt(1.-am1**2.)
                       am12=sqrt(1.-am2**2.)
                       tz=acos(-am1*am2+am11*am12*cos(raa*3.14159/180.))
                       tz=tz*180./acos(-1.)
                       pz=11.1*exp(-0.087*tz)+1.1*exp(-0.014*tz)
                       rclean=az+bz*(am1+am2)+cz*am1*am2+pz
                       rclean=rclean/4./(am1+am2)
                       step=rclean/ak1/ak2

                       if ( toa(21).lt.0.5) r0=rclean
c               END OF CHANGE: AUGUST 7, 2019
c************************************************************************                       
                       answer(isk)=zbrent(fun,x1,x2,tol)
                 

                  
                    
c                       answer(isk)=(rtoa/r0)**(1./xxx)
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 577                continue
                  
c*****************************************************
c                    write(1212,1213)  (answer(ik),ik=1,21),
c     c                   (toa(ik),ik=1,21),rtoa,r0,al,d,xxx
 1213               format(21e12.4,5x,21e12.4,5x,5f6.2)
c             INTERNAL CHECK FOR CLEAN PIXELS
                    
                           if (answer(1).gt.0.98) go to 9393
                           if (answer(2).gt.0.98) go to 9393
c****************************************************

                           
                    do 1577 isk=1,21
                      
c***********************************************************************
c                    the analysis of snow impurities
c                   ( the concentrations below 0.0001 are not reliable )
                      
                      
c     the Angstroem parameter
                       bm=0.0
                       bf=0.0
                     
                         p1=alog(answer(1))*alog(answer(1))
                         p2=alog(answer(2))*alog(answer(2))
                         bm=alog( p1/p2)/alog(w(2)/w(1))
                        
                         if (bm.lt.0.1) go to 1511
                         soda=(w(1))**bm
                         bf=soda*p1/bal
 1511                    continue
                         
c     normalized absorption coefficient of pollutants at
c                     the wavelength  1000nm
                         
                         alaska=w(1)
                         bff=p1/bal
c                    bal   -effective absorption length in microns


                         
c     type of pollutants
                         
c                         soot
                         if (bm.le.1.2)   ntype=1

c                        dust
                           if (bm.gt.1.2) ntype=2
                          
                           
c       enhancement factors for soot and ice grains respectively                   
                           BBBB=1.6
                           FFFF= 0.9
                           
c     bulk soot absorption coefficient at 1000nm
                           alfa=4.*pi*0.47/w(1)
                           
c     volumetric absorption coefficient of dust
                           DUST=0.01
                           
                 if (ntype.eq.1) conc=BBBB*bff/FFFF/alfa
                 if (ntype.eq.2) conc=BBBB*bff/DUST
                 if (bm.le.0.5)    ntype=3
                 if (bm.ge.10.)   ntype=4
                 
 1577          continue




c     Alex   09.06.2019


                       
c                   reprocessing of albedo
c     to remove gaseous absorption
c     using LPA ( linear polynomial approximation in the range 753-778nm)
              
c      comment:       in the range >865nm the physical interpolation is used
             
              x1=0.753
              x2=0.778
             
              y1=answer(12)
              y2=answer(16)
              
               afirn=(y2-y1)/(x2-x1)
               bfirn=y2-afirn*x2
               answer(13)=bfirn+afirn*0.761
               answer(14)=bfirn+afirn*0.764
               answer(15)=bfirn+afirn*0.767

               
c                 first channels:              
                 do 1989 m=1,17
                 
                       rp(m)=answer(m)**ak1
                       refl(m)=r0*answer(m)**xxx
                       
 1989               continue

                    
                if (toa(21).lt.0.5) go to 19644
                    
c:          remaining channels                    
                     do 1890 m=18,21                                
                        answer(m)=exp(-sqrt(4.*1000.*AL*pi*bai(m)/w(m)))
c     derivation of plane albedo                  
                       rp(m)=answer(m)**ak1
c     derivation of snow reflectance function                      
                       refl(m)=r0*answer(m)**xxx
 1890               continue
                    go to 1964


                    
19644           continue
                 do 1891 m=18,21                                
                  
c     derivation of plane albedo                  
                       rp(m)=answer(m)**ak1
c     derivation of snow reflectance function                      
                       refl(m)=r0*answer(m)**xxx
 1891                  continue
                    
1964       CONTINUE
            
c                                    output:            
c***************************************************            
             write(157,*) ns,ndate(3),alat,alon,(answer(i),i=1,21),isnow
             write(257,*) ns,ndate(3),alat,alon,(rp(i),i=1,21),isnow
             write(357,*) ns,ndate(3),alat,alon,(refl(i),i=1,21),isnow
             write(457,*) ns,ndate(3),alat,alon,ntype,conc,bf,bm,thv,
     c       toa(1),isnow
             write(1001,*) ns,ndate(3),alat,alon,D,area,al,r0,
     c       andsi,andbi,indexs,indexi,indexd,isnow
c***************************************************             




             
c                               FINAL STEP


             
c                                STEP  5

             
c                   CALCULATION OF BBA

c     solar flux calculation
             
             wave1=0.3
             wave2=0.7
             wave3=2.4
             wosk=0.4
             
                pp0=32.38
                pp1=-160140.33
                pp2=7959.53

        tt1=85.34*1.e-3
        tt2=401.79*1.e-3
        sol0=(pp0+pp1*exp(-0.4/tt1)+pp2*exp(-0.4/tt2))*0.1
          
                sol1= sol (wave2)      -     sol(wosk)+sol0
                sol2= sol (wave3)       -    sol(wave2)
                
                sol3= sol1  +  sol2
                
                if (isnow.eq.1) go to 1924
                
c                clean snow
                
c     plane albedo
                NSOLO=0
                    call qsimp(funp,wave1,wave2,p1)
                    rp1=p1/sol1
                   
                    call qsimp(funp,wave2,wave3,p2)
                    rp2=p2/sol2
                    rp3=(p1+p2)/sol3
                    
c     spherical albedo
                    NSOLO=1
                    call qsimp(funp,wave1,wave2,s1)
                    rs1=s1/sol1
                    call qsimp(funp,wave2,wave3,s2)
                    rs2=s2/sol2
                    rs3=(s1+s2)/sol3
                    
                    go to 1945
 1924           continue

                
c     polluted snow
                

c              NEW CODE FOR BBA OF BARE ICE
c     ALEX 29.08.2019

c     this code calculates bba analytically
      
c     input


c     wavelengths
      alam1=0.3
      alam2=0.4
      alam3=0.56
      alam4=0.7
      alam5=0.709
      alam6=0.753
      alam7=0.865
      alam8=1.02
      alam9=2.4
      
c     spherical albedo
c 1)      
c     400nm    
c      r2=0.5735
c     560
c      r3=0.5321
c    709
c      r5=0.5226
      
c     2)
      
c     754
c      r6=0.5072

c     865      
c      r7=0.4472
      
c    1020
c      r8=0.2561
c-------------------------
      r2=toa(1)
      r3=toa(6)
      r5=toa(11)
      r6=toa(12)
      r7=toa(17)
      r8=toa(21)
c      ------------------------------
c     QUADRATIC POLYNOMIAL for the range 709-865nm


      x0=alam5
      x1=alam6
      x2=alam7
      y0=r5
      y1=r6
      y3=r7
                          d1=(x0-x1)*(x0-x2)
                          d2=(x1-x0)*(x1-x2)
                          d3=(x2-x0)*(x2-x1)

                           a2=x1*x2*y0/d1+x0*x2*y1/d2+x0*x1*y2/d3
                           b2=-(x1+x2)*y0/d1-(x0+x2)*y1/d2-(x0+x1)*y2/d3
                           c2=y0/d1+y1/d2+y2/d3
                           x=x1
                           sa1=a2+b2*x+c2*x*x
                       
c     QUADRATIC POLYNOMIAL for the range 400-709nm


      x0=alam2
      x1=alam3
      x2=alam5
      y0=r2
      y1=r3
      y3=r5
                          d1=(x0-x1)*(x0-x2)
                          d2=(x1-x0)*(x1-x2)
                          d3=(x2-x0)*(x2-x1)

                           a1=x1*x2*y0/d1+x0*x2*y1/d2+x0*x1*y2/d3
                           b1=-(x1+x2)*y0/d1-(x0+x2)*y1/d2-(x0+x1)*y2/d3
                           c1=y0/d1+y1/d2+y2/d3
                           x=x1
                           sa1=a1+b1*x+c1*x*x
                        

                           
c     exponential approximation for the range above 0.865nm
                           x0=    alam7
                           x1=    alam8
                           rati=r7/r8
                           alasta=(x1-x0)/alog(rati)
                           an=1./alasta
                           p=r7*exp(x0/alasta)
                         
                           x=alam8
                           sa1=p*exp(-x/alasta)
                       
c     END of approximations for 3 intervals

c     approximation for the solar flux
c     f=f0+f1*exp(-bet*lambda/star1)+f2*exp(-gam*lambda)
                             f0=32.38
                             f1=-160140.33
                             f2=7959.53
                             bet= 1./0.08534
                             gam=1./0.40179

c     ANALYTICAL INTEGRATION OF SOLAR FLUX (DOMINATOR)

                           
                             z1=0.3
                             z2=0.7
                             
                             sol1a=f0*(z2-z1)
                             adx1= exp(-bet*z2)
                           
                             adx2=exp(-bet*z1)
                          
                             addx=dx1-dx2
                           sol1b=-f1*(  exp(-bet*z2) -exp(-bet*z1)) /bet
                           
                           sol1c=-f2*(  exp(-gam*z2)-exp(-gam*z1))/gam
                           
                           sol1=sol1a+sol1b+sol1c
c                       write(*,*) sol1a,sol1b,sol1c,bet,z1,z2,gam,f1,f2
                     
                           
                              z1=0.7
                              z2=2.4
                                   sol1a=f0*(z2-z1)
                             adx1= exp(-bet*z2)
                            
                             adx2=exp(-bet*z1)
                          
                             addx=dx1-dx2
                           sol1b=-f1*(  exp(-bet*z2) -exp(-bet*z1)) /bet
                           
                           sol1c=-f2*(  exp(-gam*z2)-exp(-gam*z1))/gam
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
              

                     
c                     write(*,*) sol1,sol2,sol3
c                     write(*,*)
c                     write(*,*)
c     END


                     
c     ANALYTICAL EQUATION FOR THE NOMINATOR
c     integration over 3 segments
                     
c     segment 1                     
                     z1=0.3
                     z2=0.7


                     
                     ajx1=a1*sol1
c*************************************                   
                     
                     ak1=(z2**2.-z1**2.)/2.
                     
                     ak2=(z2/bet+1./bet/bet)*exp(-bet*z2)-
     c                    (z1/bet+1./bet/bet)*exp(-bet*z1)
                     
                     ak3=(z2/gam+1./gam/gam)*exp(-gam*z2)-
     c                    (z1/gam+1./gam/gam)*exp(-gam*z1)
                     
                     ajx2=b1*(f0*ak1  -f1*ak2  -f2*ak3 )
c************************************************************



                     
                     am1=(z2**3.-z1**3.)/3.
                     
      am2=(z2**2./bet+2.*z2/bet/bet+2./bet/bet/bet)*exp(-bet*z2)-
     c         (z1**2./bet+2.*z1/bet/bet+2./bet/bet/bet)*exp(-bet*z1)
            
        am3=(z2**2./gam+2.*z2/gam/gam+2./gam**3.)*exp(-gam*z2)-
     c           (z1**2./gam+2.*z1/gam/gam+2./gam**3.)*exp(-gam*z1)
    
             ajx3=c1*(f0*am1 -f1*am2 -f2*am3)
c**************************************************************


             
             aj1=ajx1+ajx2+ajx3
             
                
                     
                      z1=0.7
                      z2=0.865
                 
                      
c             segment 2                      
c                     write (*,*) 'segment2'
                    

                         cajx1=a2*asol
c*************************************                   
                     
                     ak1=(z2**2.-z1**2.)/2.
                     
                     ak2=(z2/bet+1./bet/bet)*exp(-bet*z2)-
     c                    (z1/bet+1./bet/bet)*exp(-bet*z1)
                     
                     ak3=(z2/gam+1./gam/gam)*exp(-gam*z2)-
     c                    (z1/gam+1./gam/gam)*exp(-gam*z1)
                     
                     
                     cajx2=b2*(f0*ak1  -f1*ak2  -f2*ak3 )
c************************************************************
                  


                     
                     am1=(z2**3.-z1**3.)/3.
                     
      am2=(z2**2./bet+2.*z2/bet/bet+2./bet/bet/bet)*exp(-bet*z2)-
     c         (z1**2./bet+2.*z1/bet/bet+2./bet/bet/bet)*exp(-bet*z1)
            
        am3=(z2**2./gam+2.*z2/gam/gam+2./gam**3.)*exp(-gam*z2)-
     c           (z1**2./gam+2.*z1/gam/gam+2./gam**3.)*exp(-gam*z1)
c********************************************************************


        
             cajx3=c2*(f0*am1 -f1*am2 -f2*am3)
c**************************************************************

             
              
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
                 
                     bbat=ajto/sol3
                     bbavis=aj1/sol1
                     bbanir=(caj2+caj3)/sol2
c                     write(*,*) aj1,caj2,caj3,ajto
                     
c                    write(*,*) bbat, bbavis, bbanir
                     
                     rp3=bbat
                     rp1=bbavis
                     rp2=bbanir



                
c     end of modification  AUGUST 29, 2019
                     

                    
1945       CONTINUE
           write(557,33) ns,ndate(3),alat,alon,rp3,rp1,rp2,rs3,rs1,
     c      rs2,isnow
           write(5570,333)ns,ndate(3),rp3,isnow
           
           go to 87
 8797      continue
           
          
           write(55701,*)ns,ndate(3),icloud,iice
           iice=0
           icloud=0
           
 87     continue
 33     format(i15,i19,3x,8f10.4,i4)
 333    format(i15,i19,2x,f10.4,i4)
        stop
        end
      
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF RETRIEVALS!!!!!!!!!!!!!!!!

c     USED FUNCTIONS:



c                          POLLUTION MASK

                       FUNCTION sobthv(AOTs)
      
               real as(4), bs(4), cs(4), als(4), bets(4)
            
            data as/  .18016,  -0.18229,  0.15535,     -0.14223  /
            data bs/  .58331,  -0.50662,  -0.09012,        0.0207/
            data cs/  0.21475,   -0.1,  0.13639,            -0.21948/
            data als/  0.16775, -0.06969,  0.08093,     -0.08903/
            data bets/  1.09188,  0.08994,  0.49647,   -0.75218/
            
           common sza,vza,raa, g, wave,rtoa ,r0,height,aot
               wawa=0.4
                    pi=acos(-1.)

                    
                    am1=cos(sza*pi/180.)
                    am2=cos(vza*pi/180.)
                    
            
                    
                    as1=sin(sza*pi/180.)
                    as2=sin(vza*pi/180.)
                       
                    cofi=cos(raa*pi/180.)
                    amf=1./am1+1./am2
                    co=-am1*am2+as1*as2*cofi
                    
                    theta=acos(co)*180./pi

                    
c     AOT
c                    tauaer=aot
c                Alex  09.06.2019                   
                    tauaer=AOTs*(0.4/0.5)**(-1.3)
                    ad=height/7400.
                    AK=1.
                    if (ad.gt.1.e-6) AK=exp(-ad)
                    taumol=AK*0.00877/wawa**(4.05)
                    
                        tau=tauaer+taumol
                        
c          asymmetry parameter
                       g0=0.5263
                       g1=0.4627
                       wave0=0.4685                       
                       gaer=g0+g1*exp(-wawa/wave0)
                       g=tauaer*gaer/tau
                    
c                  HG phase function for aerosol
                    
                        pa=(1-g*g)/(1.-2.*g*co+g*g)**1.5
                        pr=0.75*(1.+co*co)
                        p=(taumol*pr+tauaer*pa)/tau
                   
c******************************************************************                       
c                                  SOBOLEV
                          
                       ASTRA=(1.-exp(-tau*amf))/(am1+am2)/4.
                       OSKAR=4.+3.*(1.-g)*tau
                       
                       b1=1.+1.5*am1+(1.-1.5*am1)*exp(-tau/am1)
                       b2=1.+1.5*am2+(1.-1.5*am2)*exp(-tau/am2)
                       
               RSS=p*ASTRA
               RMS=1.-b1*b2/OSKAR+(3.*(1.+g)*am1*am2-2.*(am1+am2))*ASTRA
               
               R=RSS+RMS
               
               t1=exp(-(1.-g)*tau/am1/2.)
               t2=exp(-(1.-g)*tau/am2/2.)




            
                a=0.
                b=0.
                c=0.
                
                al=0.
                bet=0.
            
                do 1    i=1,4
                   
                   if (i.eq.1) aks=1.
                   if (i.gt.1) aks=g**(i-1)
                   
                     a=a       +as(i)*aks
                     b=b       +bs(i)*aks
                     c=c        +cs(i)*aks
                     al=al      +als(i)*aks
                     bet=bet +bets(i)*aks
          
 1              continue
               
                  
                     ratm=tau*(a*exp(-tau/al)+b*exp(-tau/bet)+c)
               


               surf=t1*t2/(1.-ratm)
               
               RS=R+SURF
            
c******************************************************************
               sobthv=rs
              
               

                          RETURN
                          END

      
c     MAIN FUNCTION
c     DETERMINATION OF THE TOA reflectance as function of snow spherical albedo a




      
                         FUNCTION  fun(a)
            common sza,vza,raa, g, wave,rtoa ,r0,height,aot
                                        pi=acos(-1.)
c****************************************************
                    
                    am1=cos(sza*pi/180.)
                    am2=cos(vza*pi/180.)
                    
                    ak1=3.*(1.+2.*am1)/7.
                    ak2=3.*(1.+2.*am2)/7.
                    
                    as1=sin(sza*pi/180.)
                    as2=sin(vza*pi/180.)
                       
                    cofi=cos(raa*pi/180.)
                    amf=1./am1+1./am2
                    co=-am1*am2+as1*as2*cofi
                    
                    theta=acos(co)*180./pi

                    
c     AOT
c     tauaer=0.07*(wave/0.5)**(-1.3)
c                 ALEX 09.06.2019                    
                                 tauaer=AOT*(wave/0.5)**(-1.3)
                    ad=height/7400.
                    AK=1.
                    if (ad.gt.1.e-6) AK=exp(-ad)
                    taumol=AK*0.00877/wave**(4.05)
                    
                        tau=tauaer+taumol
                        
c          asymmetry parameter
                       g0=0.5263
                       g1=0.4627
                       wave0=0.4685                       
                       gaer=g0+g1*exp(-wave/wave0)
                       g=tauaer*gaer/tau
                    
c                  HG phase function for aerosol
                    
                        pa=(1-g*g)/(1.-2.*g*co+g*g)**1.5
                        pr=0.75*(1.+co*co)
                        p=(taumol*pr+tauaer*pa)/tau


                       

                              xxx=ak1*ak2/r0
                    
                       
                       
                      
                       
c******************************************************************                       
c                                  SOBOLEV
                          
                       ASTRA=(1.-exp(-tau*amf))/(am1+am2)/4.
                       OSKAR=4.+3.*(1.-g)*tau
                       
                       b1=1.+1.5*am1+(1.-1.5*am1)*exp(-tau/am1)
                       b2=1.+1.5*am2+(1.-1.5*am2)*exp(-tau/am2)
                       
               RSS=p*ASTRA
               RMS=1.-b1*b2/OSKAR+(3.*(1.+g)*am1*am2-2.*(am1+am2))*ASTRA
               
               R=RSS+RMS
               
               t1=exp(-(1.-g)*tau/am1/2.)
               t2=exp(-(1.-g)*tau/am2/2.)
               ratm=salbed(tau)


               surf=t1*t2*r0*a**xxx/(1-a*ratm)
               
               RS=R+SURF
            
c******************************************************************
             fun=rtoa-rs
                                             return  
                                             end

c                SPHERICAL ALBEDO OF TERRESTRIAL ATMOSPHERE:      

             FUNCTION SALBED(tau)
            common sza,vza,raa, g, wave,rtoa ,r0,height,aot
             
            real as(4), bs(4), cs(4), als(4), bets(4)
            
            data as/  .18016,  -0.18229,  0.15535,     -0.14223  /
            data bs/  .58331,  -0.50662,  -0.09012,        0.0207/
            data cs/  0.21475,   -0.1,  0.13639,            -0.21948/
            data als/  0.16775, -0.06969,  0.08093,     -0.08903/
            data bets/  1.09188,  0.08994,  0.49647,   -0.75218/
            
                a=0.
                b=0.
                c=0.
                
                al=0.
                bet=0.
            
                do 1    i=1,4
                   
                   if (i.eq.1) aks=1.
                   if (i.gt.1) aks=g**(i-1)
                   
                     a=a       +as(i)*aks
                     b=b       +bs(i)*aks
                     c=c        +cs(i)*aks
                     al=al      +als(i)*aks
                     bet=bet +bets(i)*aks
          
 1              continue
               
                  
                     salbed=tau*(a*exp(-tau/al)+b*exp(-tau/bet)+c)
                     
                  return
                  end
      
c     SOLAR SPECTRUM at GROUND level
      
                      function sol(x)
      
                           a=   32.38
                           b=   -160140.33
                           c=    7959.53

                            bet=   1./(85.34*1.e-3)
                            gam= 1./(401.79*1.e-3)
                            
                            sol=a*x-b*exp(-bet*x)/bet-c*exp(-gam*x)/gam
                            
                            if (x.lt.0.4)
     c  sol=a*0.4-b*exp(-bet*0.4)/bet-c*exp(-gam*0.4)/gam
                            
                         return
                         end


c                        analytical calculation of partial BBA integrals
      
                                  function psi(x)
      
                                  common /quatro/ x0,x1,x2,y0,y1,y2,ak1
                           p=32.38
                           q=-160140.33
                           r=7959.53
                           

                          d1=(x0-x1)*(x0-x2)
                          d2=(x1-x0)*(x1-x2)
                          d3=(x2-x0)*(x2-x1)

                           a=x1*x2*y0/d1+x0*x2*y1/d2+x0*x1*y2/d3
                           b=(x1+x2)*y0/d1+(x0+x2)*y1/d2+(x0+x1)*y2/d3
                           c=y0/d1+y1/d2+y2/d3
                           b=-b
                            bet= 1./(85.34*1.e-3)
                            gam=1./(401.79*1.e-3)
                            
                    psi1=p*(a*x+b*x*x/2.+c*x*x*x/3.)
                            
                    si2=-q*exp(-bet*x)/bet
              psi2=si2*(a+b*(x+1./bet)+c*(x*x+2.*x/bet+2./bet/bet))
                    
                    si3=-r*exp(-gam*x)/gam
             psi3=si3*(a+b*(x+1./gam)+c*(x*x+2.*x/gam+2./gam/gam))
                    
                    psi=psi1+psi2+psi3
                    
                         return
                         end

      
c     quadratic interpolation routine
      
                                           function funs(x)
      
                 real xa(168),ya(168)
                 common /quatro/ x0,x1,x2,y0,y1,y2,ak1
                 common /qua/ xa,ya,am1,AL,NSOLO              
                           p=32.38
                           q=-160140.33
                           r=7959.53
                           

                          d1=(x0-x1)*(x0-x2)
                          d2=(x1-x0)*(x1-x2)
                          d3=(x2-x0)*(x2-x1)

                           a=x1*x2*y0/d1+x0*x2*y1/d2+x0*x1*y2/d3
                           b=(x1+x2)*y0/d1+(x0+x2)*y1/d2+(x0+x1)*y2/d3
                           c=y0/d1+y1/d2+y2/d3
                           b=-b
                           
                            bet= 1./(85.34*1.e-3)
                            gam=1./(401.79*1.e-3)
                            
                  T1=p+q*exp(-bet*x)+r*exp(-gam*x)
                            
                  if (x.lt.0.4)T1=p+q*exp(-bet*0.4)+r*exp(-gam*0.4)
                            
                  T2new=a+b*x+c*x*x
                  if (T2new.LT.0.05)         T2=0.0
                  
                  
                   if(T2new.ge.0.05)         T2=T2new**ak1
                   if (NSOLO.EQ.1)           T2=T2new
                   
                  funs=T1*T2
c                  write(*,*) x,t1,t2,funs,ak1,x0,x1,x2,y0,y1,y2,ak1
                         return
                         end
      
c     BBA calculations : integrand
      
                        real function funp(x)

                 real a(6),xa(168),ya(168),dif(168)
                 integer kss(2)
                 common /qua/ xa,ya,am1,AL,NSOLO

                       ak1=3.*(1.+2.*am1)/7.
         
                   pi=acos(-1.)
  
                   ks=1
                   
          do 68 k=1,168
               
               dif(k)=abs(x-xa(k))              
               
               if ( dif(k).le.1.e-6) ks=k
               if ( dif(k).le.1.e-6) go to 69
68       continue
         
                 kss(1)=minloc(dif,1)
                 l=kss(1)
         
                 delta=x-xa(l)
                 if (delta.le.0) LS1=L-1
                 if (delta.le.0) LS2=L
                       if (delta.ge.0) LS1=L
                       if (delta.ge.0) LS2=L+1
            x0=xa(ls1)
            x1=xa(ls2)
            y0=ya(ls1)
            y1=ya(ls2)
        y=y0+(x-x0)*(y1-y0)/(x1-x0)
             go to 18
 69   continue
            y=ya(ks)
 18         continue

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
c                  imaginary part of ice refractive index
c                   at the wavelength x
                 ASTRA=Y
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
         dega=1000.*al*4.*pi*astra/x
         pow=sqrt(dega)
         
         if (pow.gt.1.e-6)rsd=exp(-pow)
         if (pow.le.1.e-6) rsd=1.
         
         if (NSOLO.eq.0) f1=rsd**ak1
         if (NSOLO.eq.1) f1=rsd
         
          

        p0=32.38
        p1=-160140.33
        p2=7959.53

        t1=85.34*1.e-3
        t2=401.79*1.e-3
        funcs=p0+p1*exp(-x/t1)+p2*exp(-x/t2)
        if (x.lt.0.4)funcs=p0+p1*exp(-0.4/t1)+p2*exp(-0.4/t2)
     
        funp=f1*funcs
                      
            return
            END


      

c                       integration routine

             SUBROUTINE qsimp(func,a,b,s)
              INTEGER JMAX
               REAL a,b,func,s,EPS
                 EXTERNAL func
                    PARAMETER (EPS=1.e-3, JMAX=20)

                       INTEGER j
                       REAL os,ost,st
               
                        ost=-1.e30
                        os= -1.e30
      
                       do 11 j=1,JMAX
          
                            call trapzd(func,a,b,st,j)
                                s=(4.*st-ost)/3.
       
                                  if (j.gt.5) then
          if (abs(s-os).lt.EPS*abs(os).or.(s.eq.0..and.os.eq.0.)) return
          endif
          os=s
          ost=st
11        continue
          
          END
      
      
c     integration routine
      
                      SUBROUTINE trapzd(func,a,b,s,n)
                       INTEGER n
                          REAL a,b,s,func
                          EXTERNAL func
                            INTEGER it,j
                            REAL del,sum,tnm,x
                 
                                if (n.eq.1) then
                            s=0.5*(b-a)*(func(a)+func(b))
       
                          else
                          it=2**(n-2)
                         tnm=it
                          del=(b-a)/tnm
                          x=a+0.5*del
                              sum=0.
                       do 11 j=1,it
                        sum=sum+func(x)
                       x=x+del
11                         continue
                              s=0.5*(s+(b-a)*sum/tnm)
        
                     endif

                     END


      
c                                  Function used to solve the equation
c                                  f(a)=0,
c                                  a - snow spherical albedo
      
            FUNCTION zbrent(func,x1,x2,tol)
      INTEGER ITMAX
      REAL zbrent,tol,x1,x2,func,EPS
      EXTERNAL func
      PARAMETER (ITMAX=100,EPS=3.e-8)
      INTEGER iter
      REAL a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      fa=func(a)
      fb=func(b)
c      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))pause
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
        fb=func(b)
11    continue
c      pause 'zbrent exceeding maximum iterations'
      zbrent=b
      
      return
      END
