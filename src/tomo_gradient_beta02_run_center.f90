PROGRAM TOMOGRADIENT
!========================================================================================================
! Calculate the tomography gradient for each layer
! Requirement:
!            1. tomography models are stored in different layer files 
!            2. each layer is for one depth range
!            3. first coloumn is latitude ,second coloumn is the longitude,then dvs
!              3.1. first line is the minimum depth range and maximum depth of this layer
!            4. model is stored along latitude value line (different longitude)
!            5. latitude  and longitude order does not matter, as long as the first point
! is either minimum(increasing) or maximu (decreasing) value. otherwise  my location 
! searching algorithm would fail.
! Algorithm:
!           1. reorganize the location to 2-D array grids
!           2. for each grid location, calculate great circle sample points for every 1 degree
! distances for every 10 degree azimuth direction.
!           3. search for grid locations surrounding these great circle sample points
!           4. calculate an average dvs using these surrounding grid locations for this sample point
!           5. so in azimuth direction, you have many sample locations, with different dvs and distance
! from the grid location
!           6. Using linear_square to fit a line to get the great for this azimuth directions
!           7. Get the maximum or minimum ( whichever their absolute value is bigger) of these gradients
! for different azimuth directions.
!           8. Then this is the gradient for this grid location, and store the direction of this gradient

!OUTPUT lat, lon, gradient, gradient_az, av_standard_deviation (of all directions used)
!=========================================================================================================

IMPLICIT NONE
INTEGER, PARAMETER :: maxp=100000
INTEGER, PARAMETER :: maxg=1000
INTEGER, PARAMETER :: layer=500
! tomography model position and Vs velocity array
REAL, DIMENSION(layer,maxp) :: gridlat,gridlon,griddvs
REAL, DIMENSION(maxg,maxg) :: gradient_absmax,gradientmax,gradientmin,gradient,gradientmax_az,gradientmin_az,gradient_az,ave_std
REAL, DIMENSION(layer)      :: layerdep,layerdepmin,layerdepmax
REAL, DIMENSION(maxg,maxg)  :: paz_lat,paz_lon,paz_dvs,dist,grid,std_maxaz,std_minaz
! Error
REAL, DIMENSION(maxp) :: az,azgradient,std
INTEGER, DIMENSION(layer) :: layerpnum
CHARACTER(len=60), DIMENSION(layer) :: tomofilename
CHARACTER(len=60) :: tomofilelist,outputfilename,form
CHARACTER(len=5) :: char_radius
INTEGER :: i,j,k,l,m,n,status_write,status_read,num_layer,laz,naz,njunk,nump_bin,con,n_az
INTEGER :: npts_lon,npts_lat,j_search,k_search,nump_lat_radius,nump_lon_radius,loc,num_tmp
INTEGER :: j_limit,k_limit,lat_up,lat_low,lon_up,lon_low,num_dist,num_az,con_lon_sys
REAL :: lat,lon,snr,slat,slon,sdv,depth,radius,baz,az_con,slope,ypoint,lat_max,lat_min,dist_tmp
REAL :: lat_step,lon_step,az_step,dist_step,err,err_lat,err_lon,diff_lon,diff_lat
REAL :: dvs_tmp1,dvs_tmp2,gridlattmp,gridlontmp,lon_max,lon_min,az_tmp
!*************************************************************************************


n=IARGC()
IF(n<1) THEN
WRITE(*,*)"USAGE: TOMOGRAD tomofilelist searching_radius"
ELSE 
CALL GETARG(1,tomofilelist)
CALL GETARG(2,char_radius)
READ(char_radius,*)radius
!WRITE(*,*)tomofilelist,radius
ENDIF

form="(F5.1 F7.1 F10.6 F7.1 F9.6)" !Output file format

!-------------------------------------------------------------------------
! Read in tomo file
!-------------------------------------------------------------------------
OPEN(UNIT=11, FILE=tomofilelist, STATUS='OLD', ACTION='read', IOSTAT=status_read)
DO i=1,layer
read(11,*,IOSTAT=status_read)tomofilename(i)
IF(status_read /= 0 ) EXIT
END DO
num_layer=i-1
CLOSE(11)

DO i=1,num_layer
OPEN(11,FILE=tomofilename(i),STATUS='OLD',ACTION='read',IOSTAT=status_read)
  READ(11,*,IOSTAT=status_read)layerdepmin(i),layerdepmax(i)
  DO j=1,maxp
   READ(11,*,IOSTAT=status_read)gridlat(i,j),gridlon(i,j),griddvs(i,j)
   IF(griddvs(i,j)>900)griddvs(i,j)=0.0
!   IF(gridlon(i,j)<0)gridlon(i,j)=360.0+gridlon(i,j)
   IF(status_read /= 0 ) EXIT
  END DO
layerpnum(i)=j-1
!write(*,*)'layer:',i,'num points:',layerpnum(i)
CLOSE(11)
END DO
!-------------------------------------------------------------------------
! Rearrange grid location to make it 2-demension
!-------------------------------------------------------------------------
DO i=1,num_layer
 ! Find out the demension,assuming gridlat is the first repeating value
write(*,*)'Processing layer:',tomofilename(i)
lat_max=MAXVAL(gridlat(i,1:layerpnum(i)))
lat_min=MINVAL(gridlat(i,1:layerpnum(i)))
lon_max=MAXVAL(gridlon(i,1:layerpnum(i)))
lon_min=MINVAL(gridlon(i,1:layerpnum(i)))
write(*,*) "lon_min",lon_min
write(*,*) "lon_max",lon_max
write(*,*) "lat_min",lat_min
write(*,*) "lat_max",lat_max
! Determine with longitude system:
!  0-360 or -180 - 180
con_lon_sys=0
if(lon_max>185)con_lon_sys=1 ! This means it is 360 system
!write(*,*)'Latitude extremes:',lat_max,lat_min,'Longitude extremes:',lon_max,lon_min
   npts_lon=0
   DO j=1,maxp
     IF(gridlat(i,j)==gridlat(i,j+1)) THEN
       npts_lon=npts_lon+1
     ELSE
       EXIT
     ENDIF
   END DO
   npts_lon=npts_lon+1
   npts_lat=layerpnum(i)/npts_lon
   IF(npts_lat*npts_lon < layerpnum(i))WRITE(*,*)'Your matrix is not complete!',npts_lat,npts_lon,layerpnum(i)
 ! Transfer to 2D array 
   DO j=1,npts_lat
    DO k=1,npts_lon
     m=k+(j-1)*npts_lon
     grid(j,k)=griddvs(i,m)
    ! write(*,*)grid(j,k),griddvs(i,m)
    END DO ! k loop
   END DO ! j loop
! So now latitude step is:
  lat_step = gridlat(i,1+npts_lon)-gridlat(i,1)
! longitude step is: 
  lon_step = gridlon(i,2)-gridlon(i,1)
! Initial value for the 2D layer is: gridlat(i,1),gridlon(i,1)
! So any grid(j,k)'s locations is: gridlat(i,1)+(j-1)*lat_step, gridlon(i,1)+(k-1)*lon_step
!write(*,*)'2D formed:',npts_lat,npts_lon,lat_step,lon_step

!-------------------------------------------------------------------------
! Calculate the gradient for each layer
!-------------------------------------------------------------------------
outputfilename='GRADIENT.'//tomofilename(i)
write(*,*)'Output file:',outputfilename
OPEN(UNIT=11,FILE=outputfilename,ACTION='WRITE',STATUS='UNKNOWN',IOSTAT=status_write)
WRITE(11,*)layerdepmin(i),layerdepmax(i),radius
   DO j=1,npts_lat
     DO k=1,npts_lon
! Search around this point grid(j,k)
     gridlattmp=gridlat(i,1)+lat_step*(j-1)
     gridlontmp=gridlon(i,1)+lon_step*(k-1)
     
! Calculate positions in every azimth directions 
      az_step = 10.0 ! in degree
      dist_step = 1 ! in degree
      num_dist=int(radius/dist_step)+1
      num_az  = int(360/az_step)
      DO m=1,num_az
       az(m)=az_step*(m-1)
! initialize
       num_tmp=num_dist*2-1
       paz_lat(m,1:num_tmp)=0.0
       paz_lon(m,1:num_tmp)=0.0
       paz_dvs(m,1:num_tmp)=0.0

! The first point is always the to be searched grid point of tomo model
       n=1
       dist(m,1)=radius*(-1)
       paz_lat(m,num_dist+1)=gridlattmp
       paz_lon(m,num_dist+1)=gridlontmp
       paz_dvs(m,num_dist+1)=grid(j,k)

         DO n=1,num_dist*2-1
           dist(m,n)=dist(m,1)+dist_step*(n-1)
           if(dist(m,n)<0) then
            dist_tmp=dist(m,n)*(-1)
            az_tmp=az(m)+180.0
           else
            dist_tmp=dist(m,n)
            az_tmp=az(m)
           endif
           if(az_tmp>360)az_tmp=az_tmp-360
           CALL azloc(gridlattmp,gridlontmp,paz_lat(m,n),paz_lon(m,n),az_tmp,dist_tmp)
         END DO ! End of n loop

         num_tmp=num_dist*2-1
         DO n=1,num_tmp
! the azloc spit out longitude in 0-360 system
           if(con_lon_sys==0 .AND. paz_lon(m,n) > 180 )paz_lon(m,n)=paz_lon(m,n)-360

! Find the surrounding tomogrid points to this position
          IF(paz_lat(m,n) < lat_max .AND. paz_lat(m,n) > lat_min) THEN
          lat_low=int((paz_lat(m,n)-gridlat(i,1))/lat_step)+1
          lat_up=lat_low+1
          lon_low=int((paz_lon(m,n)-gridlon(i,1))/lon_step)+1
          lon_up=lon_low+1
          diff_lat=paz_lat(m,n)-(gridlat(i,1)+lat_step*(lat_low-1))
          diff_lon=paz_lon(m,n)-(gridlon(i,1)+lon_step*(lon_low-1))
! Calculate the average dvs using these four closest grid          
          dvs_tmp1=grid(lat_low,lon_low)+(grid(lat_up,lon_low)-grid(lat_low,lon_low))*diff_lat/lat_step
          dvs_tmp2=grid(lat_low,lon_up)+(grid(lat_up,lon_up)-grid(lat_low,lon_up))*diff_lat/lat_step
          paz_dvs(m,n)=dvs_tmp1+(dvs_tmp2-dvs_tmp1)*diff_lon/lon_step
          ELSE
          IF(paz_lat(m,n)==90)paz_lon(m,n)=gridlontmp
          paz_dvs(m,n)=0.0
          ENDIF
         END DO ! n loop
       num_tmp = 2*num_dist-1
       CALL LINEARFIT(dist(m,1:num_tmp),paz_dvs(m,1:num_tmp),num_tmp,slope,ypoint,std(m))
!       azgradient(m)=abs(slope)
       azgradient(m)=slope
!       write(*,*) "slope",slope
      END DO ! m loop
! Find the biggest gradient
     gradientmax(j,k)=MAXVAL(azgradient(1:num_az))
     gradientmin(j,k)=MINVAL(azgradient(1:num_az))

!     IF(gradientmax(j,k)>(-1*(gradientmin(j,k))))THEN
!     gradient_absmax(j,k)=(gradientmax(j,k))
!     ELSE
!     gradient_absmax(j,k)=(gradientmin(j,k))    
!     ENDIF
!DEBUG     write(*,*) "max",gradientmax(j,k),"min",gradientmin(j,k),"absmax",gradient_absmax(j,k)

     loc=MAXLOC(azgradient(1:num_az),DIM=1)
     gradientmax_az(j,k)=az(loc)
     std_maxaz(j,k)=std(loc)
     loc=MINLOC(azgradient(1:num_az),DIM=1)
     gradientmin_az(j,k)=az(loc)
     std_minaz(j,k)=std(loc)
! Find the best gradient
     IF(abs(gradientmax(j,k))>abs(gradientmin(j,k))) THEN
       gradient(j,k)=gradientmax(j,k)
       gradient_az(j,k)=gradientmax_az(j,k)
       ave_std(j,k)=std_maxaz(j,k)
     ELSE
       gradient(j,k)=gradientmin(j,k)
       gradient_az(j,k)=gradientmin_az(j,k)
       ave_std(j,k)=std_minaz(j,k)
     ENDIF
     IF(gradient_az(j,k)>180)gradient_az(j,k)=(gradient_az(j,k)-360)
     gradient(j,k)=abs(gradient(j,k)) !Only write positive gradients
!     gradient(j,k)=(gradient(j,k)) !Write positive and negative gradients
     WRITE(11,form)gridlattmp,gridlontmp,gradient(j,k),gradient_az(j,k),ave_std(j,k)
     !WRITE(*,*)gridlattmp,gridlontmp,gradient(j,k),gradient_az(j,k)
    END DO ! k loop
   END DO ! j loop
CLOSE(UNIT=11)
END DO ! i loop
END PROGRAM TOMOGRADIENT

