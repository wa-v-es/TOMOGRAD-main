          subroutine azloc(elat_input,elon_input,stlat,stlon,az,dist)
! Using user input source location ( lat lon  ), azimuth and distance
! To calculate the location from the source.
! NOTE: THIS CODE WILL GIVE OUT LONGITUDE: -180~180
          parameter (PI=3.1415926)
          real elat,elon,plon,elat_input,elon_input
          real stlat,stlon
          real az,dist,tmp
! Note the distance here is always the great circle distance

! convert the elon into uniform elon I want
           if(elon_input.LT.0) then
           p=1
           elon=360+elon_input
           else
           p=0
           elon=elon_input
           endif

           darc=3.1415926/180
           arc=1/darc
           plon=elon*darc
           a=(90-elat_input)*darc
           C=az*darc
           b=dist*darc
! Solve the great circle triangle.
!                        /\
!                       /  \
!                      /    \
!                   a /      \c
!                    /        \
!                   /)C________\
!                        b
! we know length of the side a and b and their angle C, to calculate c

          side_c=acos(cos(a)*cos(b)+sin(a)*sin(b)*cos(C))
! NOTE: side_c is only from 0 to PI, but stlat=90-side_c*arc will move that to correct domain

          tmp=sin(b)/sin(side_c)*sin(C)
          if(tmp.GT.1.0)tmp=1.0
          if(tmp.LT.-1.0)tmp=-1.0
          angleB1=asin(tmp)
! NOTE: angleB is only from -PI/2 to PI/2, we want to extend this to -PI to PI
! SO for every angleB, then we should have two result within -PI to PI, which is: angleB and sign(angleB)*PI-angleB
          angleB2=sign(PI,angleB1)-angleB1

! Determine which value should we use,using cosine law for the side b which is coresponded to angleB
          condition=(cos(b)-cos(a)*cos(side_c))/(sin(a)*sin(side_c))
          con=abs(sign(condition,cos(angleB1))+condition)
  
          if(con.GT.0.0000001) then
           angleB=angleB1
          else
           angleB=angleB2
          endif

          stlat=90-side_c*arc 

          stlon=(plon+angleB)*arc

          !if(stlon.GT.180)stlon=stlon-360

!          write(*,*)'STATION LOCATION FOR AZ:',az,dist,'is:',stlat,stlon

          end
