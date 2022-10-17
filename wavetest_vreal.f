C****************************************************************************
C WAVETEST: Example Fortran program for WAVELET, using ADCP dataset
C (modified from NINO3 SST dataset program)
C
C COMPILE:  f77 -fno-second-underscore -o wave.o  -L/usr/local/lib
C           -I/usr/include wavetest.f  -lnetcdf 
C
C INPUT: change n: no. of ensembles, m: depth, infile: input file name
C        dt: sampling interval in hours, jtot: check formula in wavelet.f       
C See "http://paos.colorado.edu/research/wavelets/"
C
C  Copyright (C) 1998, Christopher Torrence and Gilbert P. Compo
C This software may be used, copied, or redistributed as long as it is not
C sold and this copyright notice is reproduced on each copy made.  This
C routine is provided as is without any express or implied warranties
C whatsoever.
C
C Modified: November 1999 by Arjan van Dijk to include IMPLICIT NONE and
C           to convert all routines to DOUBLE precision.
c Modified: May 2009 by Amol Prakash to include input moored adcp netcdf file
C           and returns output netcdf file
C****************************************************************************
c to run
c f77 -o exe  -fno-second-underscore -L/usr/local/netcdf/liblcs -I/usr/local/netcdf/includelcs -O wavetest_real.f -lnetcdf -lnetcdff 
      PROGRAM wavetest

      include '/usr/include/netcdf.inc'
c      IMPLICIT none


      INTEGER n,m,subscale,jtot
      DOUBLE PRECISION dt,s0,dj

      PARAMETER (n=8895,m=1,dt=0.166667,s0=dt)                     
C          xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx change
      INTEGER ncidin,alowid,olowid,depthid,tid,ilatid,ilonid
      INTEGER idepthid,depthdim,dime2(2)  
      INTEGER depthlen,tlen 
  
      CHARACTER*60 outfile,infile

            
      
      DOUBLE PRECISION time(n),depth(m),coidat,vrange(2),miss      
      REAL  alow(n,m),olow(n,m),alow1,olow1
C these parameters depend on the particular time series
c      PARAMETER (n=5367,dt=1D0,s0=dt)
      PARAMETER (subscale=12)
      PARAMETER (dj=1.D0/subscale,jtot=10*subscale)



      INTEGER status,ncid        
      INTEGER periodid, timeid , waveid,waverealid, coiid,kmax
      INTEGER perioddim, timedim, dime1(4),iindex(4)
      INTEGER londim, latdim, imax ,jmax ,lonid,latid
      DOUBLE PRECISION lon, lat ,ncoi(m,n),wavedat,waverealdat     
        
C Note: for accurate reconstruction and wavelet-derived variance
C     do not pad with zeroes, set s0=dt (for Paul set s0=dt/4), and use
C     a large "jtot" (even though the extra scales will be within
C     the cone of influence).
C     For plotting purposes, it is only necessary to use
C     s0=2dt (for Morlet) and "jtot" from Eqn(10) Torrence&Compo(1998).

      INTEGER mother,ibase2,npad
      DOUBLE PRECISION sst(n),recon_sst(n),param,pi !,date(n)
      DOUBLE PRECISION scale(jtot),period(jtot),coi(n),nperiod(jtot)
      DOUBLE COMPLEX wave(n,jtot),nwave(m,n,jtot),nwavereal(m,n,jtot)

      INTEGER i,j,isigtest,javg1,javg2
      DOUBLE PRECISION lag1,siglvl,dof(jtot)
      DOUBLE PRECISION fft_theor(jtot),signif(jtot),ymean,variance
      DOUBLE PRECISION recon_mean,recon_vari
      DOUBLE PRECISION Cdelta,psi0
      DOUBLE PRECISION global_ws(jtot),global_signif(jtot)
      DOUBLE PRECISION savg_dof(jtot),savg_signif(jtot),sstENSO(n)
      
      DATA vrange /-999.0D0,999.0D0/
      DATA miss /99999.0D0/    
C*************Read Input Netcdf file**********************************              
C        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx change
      infile  = 'residuals_detided_l1_1.nc'  
      outfile = 'l1_1_vreal.nc'
   
      status = nf_open(infile, 0, ncidin)
C        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx change
      status = nf_inq_dimid(ncidin,'LON',ilonid)      
      status = nf_inq_dimid(ncidin,'LAT',ilatid)
      status = nf_inq_dimid(ncidin,'DEPTH',idepthid)
      status = nf_inq_dimid(ncidin,'TIME',tid)
     
      status = nf_inq_dimlen(ncidin,idepthid,depthlen)          
      status = nf_inq_dimlen(ncidin,tid,tlen)

        

C        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx change
      status = nf_inq_varid(ncidin,'URES1', alowid)          
      status = nf_inq_varid(ncidin,'VRES1', olowid)

      status = nf_get_var_double (ncidin, tid, time)
      status = nf_get_var_double (ncidin, idepthid, depth)
      status = nf_get_var_double (ncidin, ilonid, lon)
      status = nf_get_var_double (ncidin, ilatid, lat)

        imax = 1
        jmax = 1


c            status = nf_get_var_real(ncidin, alowid,alow1)

c            status = nf_get_var_real(ncidin, olowid,olow1)
c            write(*,*)alow1    
      DO 500 l=1,n 
         iindex(4)=l
        DO 520 k=1,m
         iindex(3)=k
         DO 540 j=1,jmax
          iindex(2)=j
          DO 560 i= 1,imax
            iindex(1)=i    
            status = nf_get_var1_real(ncidin, alowid,iindex,alow1)
            status = nf_get_var1_real(ncidin, olowid,iindex,olow1)
            alow(l,k) = alow1
            olow(l,k) = olow1
560       CONTINUE  
540      CONTINUE  
520     CONTINUE  
500   CONTINUE           
                
  
       pi = 4.D0*ATAN(1.D0)
       ibase2 = NINT(LOG(DBLE(n))/LOG(2.D0))+1
       npad = INT(2.D0**ibase2)
C      npad = n  ! this is for no padding with zeroes

C*************************************************** Wavelet transform

C** let the WAVELET subroutine choose the defaults for these:
       mother = 0
       param = 6.D0

C** read in the NINO3 SST data
c      OPEN(UNIT=11,FILE='input.dat',STATUS='old')
c      READ(11,*) sst
c      CLOSE(11)
c      PRINT'(/,"sst(1)=",F6.2,"  sst(n) = ",F6.2,/)',sst(1),sst(n)

C        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx change
      DO 400 jj = 1,m  
         DO 420 ii = 1, n
            sst(ii)=olow(ii,jj)
420      CONTINUE    
C** get the wavelet transform
      CALL WAVELET(n,sst,dt,mother,param,s0,dj,jtot,npad,
     &             wave,scale,period,coi)


C*************************************************** Significance testing

C** local significance test
      isigtest = 0
      lag1 = 0.72D0
      siglvl = 0.05D0
      CALL WAVE_SIGNIF (isigtest,n,sst,dt,mother,param,dj,jtot,
     &       scale,period,lag1,siglvl,dof,fft_theor,signif,
     &       ymean,variance,Cdelta,psi0)


C** global wavelet spectrum & significance test
      isigtest = 1
      lag1 = 0.72D0
      siglvl = 0.05D0
      DO 10 j=1,jtot
        DO 20 i=1,n
          global_ws(j) = global_ws(j) + ABS(wave(i,j))**2
20      CONTINUE
        global_ws(j) = global_ws(j)/n
        dof(j) = n - scale(j)
10    CONTINUE

      CALL WAVE_SIGNIF (isigtest,n,sst,dt,mother,param,dj,jtot,
     &       scale,period,lag1,siglvl,dof,fft_theor,global_signif,
     &       ymean,variance,Cdelta,psi0)


C** scale-average time series & significance test
      isigtest = 2
      lag1 = 0.72D0
      siglvl = 0.05D0
C    scale average between 2 and 7.9 years
      savg_dof(1) = 2.0D0
      savg_dof(2) = 7.9D0
C    find the "j"-values that correspond to savg_dof(1) & savg_dof(2)
      javg1 = 0
      javg2 = 0
      DO 30 j=1,jtot
        IF ((scale(j).GE.savg_dof(1)).AND.(javg1.EQ.0)) javg1 = j
        IF (scale(j).LE.savg_dof(2)) javg2 = j
30    CONTINUE
C   call wave_signif first, to get the value of "Cdelta"
      CALL WAVE_SIGNIF (isigtest,n,sst,dt,mother,param,dj,jtot,
     &     scale,period,lag1,siglvl,savg_dof,fft_theor,savg_signif,
     &     ymean,variance,Cdelta,psi0)
C   construct the scale-averaged time series [Eqn(24)]
      DO 50 i=1,n
        sstENSO(i) = 0.D0
        DO 60 j=javg1,javg2
          sstENSO(i) = sstENSO(i) + (ABS(wave(i,j))**2)/scale(j)
60      CONTINUE
        sstENSO(i) = dj*dt*sstENSO(i)/Cdelta
50    CONTINUE

C** construct the wavelet derived variance (Parseval's theorem)  [Eqn(14)]
C   Cdelta & psi0 are returned from WAVE_SIGNIF
      recon_vari = 0.D0
      DO 900 i=1,n
        DO 1000 j=1,jtot
          recon_vari = recon_vari + (ABS(wave(i,j))**2)/scale(j)
1000    CONTINUE
900   CONTINUE
      recon_vari = dj*dt*recon_vari/(Cdelta*n)
c      PRINT'(A,F14.5)',' Reconstructed variance=',recon_vari
c      PRINT'(A,F14.5)',' Original variance   =',variance
c      PRINT'(A,F14.5,A,/)',' Ratio = ',recon_vari/variance,
c     &     ' (this is low due to padding with zeroes)'

C** reconstruct the time series [Eqn(11)]
C   check mean and RMS difference of reconstructed time series
      recon_mean=0.D0
      recon_vari = 0.D0
      DO 1100 i=1,n
        recon_sst(i)=0.D0
        DO 1200 j=1,jtot
          recon_sst(i) = recon_sst(i)+(DBLE(wave(i,j)))/SQRT(scale(j))
1200    CONTINUE
        recon_sst(i) = dj*SQRT(dt)*recon_sst(i)/(Cdelta*psi0)
        recon_vari = recon_vari+(sst(i)-ymean-recon_sst(i))**2
        recon_mean = recon_mean + recon_sst(i)
1100  CONTINUE
      recon_mean = recon_mean/n
      recon_vari = SQRT(recon_vari/n)

c      PRINT'(A,F14.6)',' Reconstructed mean=',recon_mean
c      PRINT'(A,F14.6)',' Original mean   =',ymean
c      PRINT'(A,F14.6,/)',' Root-mean-square difference of time series=',
c     &      recon_vari

     
   
        DO 440 j = 1,jtot
          nperiod(j) = period(j)
          DO 460 i=1,n
           IF(i.EQ.1.or.i.EQ.n)THEN
           ELSEIF(j.EQ.1)THEN
         nwavereal(jj,i,j) = realpart(wave(i,j))
           ELSE 
          nwavereal(jj,i,j) = realpart(wave(i,j))
           ENDIF     
460       CONTINUE
440     CONTINUE      
        

 

       
400   CONTINUE
        
      do i=1,jtot
         nperiod(i)=1*(nperiod(i))/24
!        nperiod(i) = log10(nperiod(i))/(log10(2.))
      enddo
        
!     date(1) = 1 
!     do i=2,n
!       date(i)=date(i-1)+dt
!     enddo
    
C**************************************************************
C                   creation of netcdf file

        
C create output file      
      status = nf_create(outfile,nf_clobber,ncid)
      if (status.ne.nf_noerr) stop 'error with nf_open(outfil)'

c                   #########################                           
C define dimensions

      status = nf_def_dim(ncid,'depth',m,depthdim)
      if (status.ne.nf_noerr) stop 'error with depthdim'      
  
      status = nf_def_dim(ncid,'period',jtot,perioddim)
      if (status.ne.nf_noerr) stop 'error with perioddim'


      status = nf_def_dim(ncid,'Time',n,timedim)
      if (status.ne.nf_noerr) stop 'error with time'

c                   #########################                
C define variables
      dime1(1) = perioddim
      dime1(2) = depthdim         
      dime1(3) = timedim

      dime2(1)=  timedim   
      dime2(2)=  depthdim

        
C#######  depth 
      status=nf_def_var(ncid,'depth',nf_float,1,depthdim,depthid)
      if (status.ne.nf_noerr) stop 'error with depthid'

      status=nf_copy_att(ncidin,idepthid,'long_name',ncid,depthid)
      if (status.ne.nf_noerr) stop 'error with depthatt1'
      status=nf_copy_att(ncidin,idepthid,'units',ncid,depthid)
      if (status.ne.nf_noerr) stop 'error with depthatt2'
      status=nf_copy_att(ncidin,idepthid,'positive',ncid,depthid)
      if (status.ne.nf_noerr) stop 'error with depthatt3'
      status=nf_copy_att(ncidin,idepthid,'point_spacing',ncid,depthid)
      if (status.ne.nf_noerr) stop 'error with depthatt4'
      status=nf_copy_att(ncidin,idepthid,'axis',ncid,depthid)
  
C####### period
      status = nf_def_var(ncid,'period',nf_double,1,perioddim,
     &    periodid)
      if (status.ne.nf_noerr) stop 'error with period1'

      status=nf_put_att_text(ncid,periodid,'axis',1,'Y')
      if (status.ne.nf_noerr) stop 'error with periodatt1' 
 

C######## Time 
      status=nf_def_var(ncid,'Time',nf_float,1,timedim,timeid)
      if (status.ne.nf_noerr) stop 'error with time1'

!      status=nf_copy_att(ncidin,tid,'long_name',ncid,timeid)  
!      if (status.ne.nf_noerr) stop 'error with timeatt1'
      

!      status=nf_put_att_text(ncid,timeid,'units',31,
!      &'hours since 2008-03-04 09:39:00')  
      status=nf_copy_att(ncidin,tid,'units',ncid,timeid)  
      if (status.ne.nf_noerr) stop 'error with timeatt2'
      status=nf_copy_att(ncidin,tid,'axis',ncid,timeid)  
      if (status.ne.nf_noerr) stop 'error with timeatt3'
!      status=nf_put_att_text(ncid,timeid,'time_origin',20,
!     &'04-MAR-2008 09:39:00')  
      status=nf_copy_att(ncidin,tid,'time_origin',ncid,timeid)  
!      if (status.ne.nf_noerr) stop 'error with timeatt4'


C wavereal 

      status=nf_def_var(ncid,'wavereal',nf_float,3,dime1,waverealid)
      if (status.ne.nf_noerr) stop 'error with wave1'

   
      
       status=nf_enddef(ncid)
C                   ######################### 
C write data

      status= nf_put_var_double(ncid,depthid,depth)
      if (status.ne.nf_noerr) stop 'error with depth'


       status = nf_put_var_double(ncid,periodid,nperiod)
       if (status.ne.nf_noerr) stop 'error with period'

                
      status = nf_put_var_double(ncid,timeid,time)
      if (status.ne.nf_noerr) stop 'error with time'
        

      DO 600  l=1, n
         iindex(3)=l
        DO 610  k=1,jtot
         iindex(1)=k
         DO 620 j=1,m
          iindex(2)=j      
           waverealdat = nwavereal(j,l,k)     
         status = nf_put_var1_double(ncid,waverealid,iindex,
     &     waverealdat)
             if (status.ne.nf_noerr) stop 'error with wavereal'

620      CONTINUE
610     CONTINUE
600   CONTINUE
  
      END



      include 'wavelet.f'
      include 'cfftpack.f'
      include 'chisqr.f'  
