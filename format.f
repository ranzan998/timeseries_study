        program format

!Program to format input file for task (tide.in).
!Input file for the program is variable listed from netcdf file
!without formatting and no header (tide.txt) 
        real time1,time2,d1,d2,d3,d4,d5
        integer day,year,sno,flag,jday,n,nl
        character*3 month
        character*14 nc
        character*100 st1, st 
        character*1 check
               
        data d2,d3,d4,d5/4*0.00/
        
        open(11,file='tide.txt',status='old')
        open(12,file='tide.in',status='unknown')
        
        n = 0
        do i = 1,100000
         read(11,*,end=80)
         n= n+1
        enddo 
80      write(nc,*) n
        nl = len_trim(adjustl(nc))
        write(nc,*)nl
        nc = adjustl(nc)
        st ='(I3,1x,A3,1x,I4,1x,F2.0,1x,F2.0,A1)'
        rewind 11
        read(11,st) day, month, year, time1, time2, check
        if (check.eq.':')then
          st1 ='(I3,1x,A3,1x,I4,1x,F2.0,1x,F2.0,1x,F2.0,3x,I'//trim(nc)
     &        //',2x,F7.2)'  
        else
          st1 ='(I3,1x,A3,1x,I4,1x,F2.0,1x,F2.0,3x,I'//trim(nc)
     &        //',2x,F7.2)'  

        endif
        

        do i = 1,20
                write(12,*)
        enddo

        rewind 11
        do i =1,n
           if (check.eq.':')then
              read(11,st1)day,month,year,time1,time2,time3,sno,d1
              time1 = time1 + (time2/60) + (time3/3600) 
              
           else 
               read(11,st1)day,month,year,time1,time2,sno,d1
               time1 = time1 + (time2/60) 
           endif                    
               call doy(day,month,year,jday) 
               if(d1.lt.-100)then
                  flag = 1
                  write(12,52)sno,flag,year,jday,time1,d1,d2,d3,d4,d5
                              
               else
                  flag = 0 
                  write(12,52)sno,flag,year,jday,time1,d1
     &                       ,d2,d3,d4,d5
               endif   
                
        enddo       
        
52      format(I6,1X,I1,1X,I4,1X,I3,F7.3,5F8.2)
        stop
        end 


        subroutine doy(day,cmonth,year,dayofyear)
!       subroutine to calculate day of the year        

        integer day,year,dayofyear,imonth
        character*3 cmonth

        dayofyear = day
 
        if(cmonth.eq.'JAN')then
                imonth=1
        elseif(cmonth.eq.'FEB')then
                imonth=2
        elseif(cmonth.eq.'MAR')then
                imonth=3
        elseif(cmonth.eq.'APR')then
                imonth=4
        elseif(cmonth.eq.'MAY')then
                imonth=5
        elseif(cmonth.eq.'JUN')then
                imonth=6
        elseif(cmonth.eq.'JUL')then
                imonth=7
        elseif(cmonth.eq.'AUG')then
                imonth=8
        elseif(cmonth.eq.'SEP')then
                imonth=9
        elseif(cmonth.eq.'OCT')then
                imonth=10
        elseif(cmonth.eq.'NOV')then
                imonth=11
        elseif(cmonth.eq.'DEC')then
                imonth=12
        else
               write(*,*)'ERROR:SUBROUTINE(doy)-month wrong format'
        endif

        
        if(mod(year,400).eq.0)then
                leapday = 1
        elseif(mod(year,100).eq.0)then
                leapday = 0
        elseif(mod(year,4).eq.0)then
                leapday = 1
        else
                leapday = 0
        endif

        do i=1,imonth-1
           select case(i)
                case(1,3,5,7,8,10,12)
                        dayofyear = dayofyear + 31
                case(4,6,9,11)
                        dayofyear = dayofyear + 30
                case(2)
                        dayofyear = dayofyear + 28 + leapday
           endselect             
        enddo

        return

        end

