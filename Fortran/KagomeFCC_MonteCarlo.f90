!stack_module
!   Module for creating and using the Stack datatype

module stack_module

!Define types for each element and the stack (last element in linked list)
type ELEMENT_TYPE
 integer(4) :: VAL
 type(ELEMENT_TYPE), Pointer :: PREV
end type ELEMENT_TYPE

type STACK_TYPE
 integer(4) :: SIZE=0
 type(ELEMENT_TYPE), Pointer :: LASTIN=>NULL()
end type STACK_TYPE

!Define subroutines and functions for push, pop and isEmpty
contains
 !push takes a value and adds it to the end of the linked list/top of the stack
 subroutine push(VAL_,STACK)
  implicit none
  !data types for push
  integer(4), intent(in) :: VAL_
  type(STACK_TYPE), intent(inout) :: STACK
  type(ELEMENT_TYPE), pointer :: CURRENT
  !Initialize Current
  allocate(CURRENT)
  CURRENT%VAL= VAL_
  CURRENT%PREV => STACK%LASTIN
  !Add Current to Stack
  STACK%LASTIN => CURRENT
  STACK%SIZE = STACK%SIZE+1
  RETURN
 end subroutine push

 !pop takes the value from the end of the linked list/top of the stack and returns and removes it from the stack
 subroutine pop(VAL_,STACK)
  implicit none
  !data types for pop
  type(STACK_TYPE), intent(inout) :: STACK
  integer(4), intent(out) :: VAL_
  type(ELEMENT_TYPE), pointer :: B4LASTIN
  !Write to VAL_
  if (associated(STACK%LASTIN)) then
   VAL_ = STACK%LASTIN%VAL
   !Take out the last in element
   B4LASTIN => STACK%LASTIN%PREV
   deallocate(STACK%LASTIN)
   STACK%LASTIN => B4LASTIN
   STACK%SIZE=STACK%SIZE-1
  else
   if(STACK%SIZE.ne.0) then
    print*, STACK%SIZE
    stop 'MISMATCH BETWEEN STACKSIZE AND POINTER: BAD BOOK KEEPING'
   endif
  endif
  RETURN
 end subroutine pop

 !isEmpty returns true if and only if the stack is empty (accounts for mismatches)
 logical function isEmpty(STACK)
  implicit none
  !data types for isEmpty
  type(STACK_TYPE), intent(in) :: STACK
  if(.not. (associated(STACK%LASTIN))) then
   if(STACK%SIZE.eq.0) then
    isEmpty=.true.
   else
    stop 'MISTMATCH BETWEEN STACKSIZE AND POINTER: BAD BOOK KEEPING'
   endif
  else
   isEmpty=.false.
  endif
 end function isEmpty
end module stack_module

module input_module_3d
! define parameters
 integer,parameter :: L=6
         
  double precision, parameter :: pi=3.14159265358979

  double precision :: Tmax1,Tmin1
  
  integer,parameter    :: ntemp1=0, ntrans=0, nmeas = 4000,nsite=L*L*L

  integer, parameter :: TIME_SERIES_STEPS=100
  integer, parameter :: BACKUP_STEPS=4000

  double precision, parameter :: P_CLUSTER=1.0d0

  logical, parameter :: USE_CLUSTER=.TRUE.
       
       
! declare types
double precision :: E,mav,mav2,mx,my,mz,T,Eav2,Eav,r,ct,phi,st,m1,m2,dh,x,mzav,mxav,myav,eav4,UE,jex,di
INTEGER ::TABLE(L,L,L),ITABLE(L*L*L,3),sm(L*L*L),isub(L*L*L)
integer :: nmax,IS,m,it,ix,iyy,iz,ns
double precision :: r250_
double precision :: R2,cp,sp,s1,s2,rr,CVT,dum,thet,u1,u2,u3,rnscx,rnscy,rnscz,sdn,rnd
double precision :: hpx1,hpy1,hpz1,hmp,hm,hmi,hx,hy,hz,sth,sthi,cph,sph,y,m0z,m0zav,&
    & m0x,m0xav,m0y,m0yav, e_Prev
double precision ::  sx(L*L*L),sy(L*L*L),sz(L*L*L),D(3,3,0:L-1,0:L-1,0:L-1),s(L*L*L,3),nx(L*L*L),ny(L*L*L),nz(L*L*L)
double precision :: smx(0:3),smy(0:3),smz(0:3),sm2av(0:3),smav(0:3),s0x(L*L*L),s0y(L*L*L),s0z(L*L*L),smxalt(0:3),smyalt(0:3)
double precision :: rmx,rmy,rmz,rnum,rmav,rmav2,rmxav,rmyav,rmzav,rm2root,rmav4,smzalt(0:3),sm2avalt(0:3),smavalt(0:3)
double precision :: UB,nxsum,nysum,nzsum,ent,entav,beta1,beta2,t1,t2
integer, dimension(0:3,0:(nsite/4)) :: inv_isub    !Inverse sublattice lookup
integer :: imeas,irand0,iseed,ran,rmeas,I2count,tmeas,I2check
integer :: flush
integer :: i,j,k,dx,dy,dz,ax,ay,az,bx,by,bz,cx,cy,cz
integer :: a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3
Logical :: backUp_file, tran_file, series_file
external r250_, flush

common /table1/ Table,itable,D,sm,isub,smx,smy,smz,mx,my,mz,inv_isub
common /spins/ sx,sy,sz,nx,ny,nz,beta1,beta2    
!
   end module input_module_3d

module cluster_module
use input_module_3d
use stack_module

!Declare cluster parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Test values
double precision, parameter :: kb=1.d0
double precision, parameter :: Js=1.d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Declare cluster types
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer, dimension(0:3) :: s0_Start !   Starting point for sublattice cumulative sums
integer, dimension(0:3,1:3) :: s0_Pos, s0_Rel  ! Starting and Relative position in cumulative sums
integer, dimension(1:3) :: curr_Pos,wc_Pos   !Position of the relative and seed spins of the selected spins

integer :: sub_I,sub_J,sub_K    !   Subroutine iterators
integer :: spin_Curr, sub_Curr    ! Spin currently being considered and it's sublattice
integer :: wc,zc, wc_Prime    ! Holds lattice sites for cluster growing as defined by Luijten, Blote and Baek
integer :: jseed, counterrr, cc1, cc2, cs1, cs2, cs3, cs4 ! Utility

double precision, dimension(1:nsite) :: sx_New, sy_New, sz_New  ! New vector components of each spin added.
double precision, dimension(1:3) :: ref_Norm    ! The normal vector to the reflection plane

double precision, dimension(0:3,1:nsite) :: s0_Cum ! Cumulative distribution as defined by Luijten, Blote and Baek
double precision, dimension(1:nsite,1:nsite) :: Rij, Jij ! Distance and interaction strength between two lattice points.
double precision, dimension(1:nsite,1:nsite) :: Rx, Ry, Rz   ! Distance in each direction between any two sites

double precision :: ref_Norm_L, s_New_L   ! The length of the vectors of the normal to the reflection plane, and the new spins.
double precision :: Rx_Curr, Ry_Curr, Rz_Curr   ! Holds position elements of current spin pair
double precision :: expo1_1x,expo1_1y,expo1_1z,expo1_2x,expo1_2y,expo1_2z   ! Terms and directions for P_Add numerator exponent, labelled expo1
double precision :: expo1, expo2    ! Numerator and denominator exponent of P_Add, respectively.
double precision :: P_Add, P_Acc    ! Probabilities to add to cluster, and accept flip
double precision :: proj    ! The projection of the spin vector to the norm of the plane
double precision :: dE_Bulk, dE_Surf, surf_Curr, bulk_Curr  ! The change in energy in both bulk and surface and the current term for each sum
double precision :: rand, uc    !  Holds random numbers

logical, dimension(1:nsite) :: cluster_Spin ! Returns true if spin is within cluster

Type(STACK_TYPE) :: cluster_Stack   ! Stack to hold spins while cluster is being built


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module cluster_module



program kag_dipole_3d
use input_module_3d

call init_random_seed() !initialize random number generator

!Open(UNIT=1,FILE='temperature.input',STATUS='OLD')
!read(1,*) Tmax1, Tmin1
Tmax1=.10
Tmin1=.10
close(1)

! ns counts number of sites
!nz0 counts sublattice D  (0,0,0)
!nz1 counts sublattice A  (0,1,1)
!nz2 counts sublattice B  (1,1,0)
!nz3 counts sublattice C  (1,0,1)

ns=0
nz0=0
nz1=0
nz2=0
nz3=0
!ix,iyy,iz are integer components along non-orthogonal basis vectors
do ix=1,L
do iyy=1,L
do iz=1,L
Ns=Ns+1
!set all spin lengths to unity
sm(ns)=1
!set spin length to zero on D sites for fcckagome and to unity for fcc
if(((ix+iyy-2*((ix+iyy)/2)).eq.0).and.((iyy+iz-2*((iyy+iz)/2)).eq.0).and.((ix+iz-2*((ix+iz)/2)).eq.0)) then
sm(ns)=0
nz0=nz0+1
isub(ns)=0
inv_isub(0,nz0)=Ns
endif

!label sublattice given site number  

if(((ix+iyy-2*((ix+iyy)/2)).eq.0).and.((iyy+iz-2*((iyy+iz)/2)).ne.0).and.((ix+iz-2*((ix+iz)/2)).ne.0)) then 
isub(ns)=1
nz1=nz1+1
inv_isub(1,nz1)=Ns
endif
if(((ix+iyy-2*((ix+iyy)/2)).ne.0).and.((iyy+iz-2*((iyy+iz)/2)).eq.0).and.((ix+iz-2*((ix+iz)/2)).ne.0)) then
isub(ns)=2
nz2=nz2+1
inv_isub(2,nz2)=Ns
endif
if(((ix+iyy-2*((ix+iyy)/2)).ne.0).and.((iyy+iz-2*((iyy+iz)/2)).ne.0).and.((ix+iz-2*((ix+iz)/2)).eq.0)) then
isub(ns)=3
nz3=nz3+1
inv_isub(3,nz3)=Ns
endif
!record site number for each space vector
table(ix,iyy,iz)=ns
!create inverse table
itable(ns,1)=ix
itable(ns,2)=iyy
itable(ns,3)=iz
enddo
enddo
enddo

!set all spin components to zero
sx=0
sy=0
sz=0
nocc=0
do i=1,nsite
nocc=nocc+sm(i)
enddo
!construct the dipole matrix for the given choice of basis vectors(in the subroutine)
call gamma(L,D)
!set up a ground state when D has sm=0
do i=1,nsite
m=itable(i,1)+itable(i,2)+itable(i,3)
s0x(i)=0
s0y(i)=0
s0z(i)=0
if(isub(i).eq.1) then
s0x(i)=1/dsqrt(2.d0)*(-1)**m
s0y(i)=-1/dsqrt(2.d0)*(-1)**m
s0z(i)=0
endif
if(isub(i).eq.2) then
s0x(i)=0
s0y(i)=1/dsqrt(2.d0)*(-1)**m
s0z(i)=-1/dsqrt(2.d0)*(-1)**m
endif
if(isub(i).eq.3) then
s0x(i)=-1/dsqrt(2.d0)*(-1)**m
s0y(i)=0
s0z(i)=1/dsqrt(2.d0)*(-1)**m
endif
enddo
sx=s0x
sy=s0y
sz=s0z
call measure
write(6,*) e/nocc

!set up a rotation of the ground state about a given unit vector uu1,uu2,uu3 to check energy
uu1=1/dsqrt(2.d0)
uu2=-uu1
uu3=0
do i5=0,200
thet=i5*pi/100.
do i=1,nsite
if(sm(i).ne.0) then
sdn=s0x(i)*uu1+s0y(i)*uu2+s0z(i)*uu3
rnscx=uu2*s0z(i)-uu3*s0y(i)
rnscy=uu3*s0x(i)-uu1*s0z(i)
rnscz=uu1*s0y(i)-uu2*s0x(i)
sxx=sdn*uu1*(1-dcos(thet))+s0x(i)*dcos(thet)+rnscx*dsin(thet)
syy=sdn*uu2*(1-dcos(thet))+s0y(i)*dcos(thet)+rnscy*dsin(thet)
szz=sdn*uu3*(1-dcos(thet))+s0z(i)*dcos(thet)+rnscz*dsin(thet)
sx(i)=sxx
sy(i)=syy
sz(i)=szz

endif
enddo
call measure
write(98,*) i5,thet,e/nocc
enddo
do i=1,nsite
do j=0,3
if(isub(i).eq.j) then
write(93+j,990) itable(i,1)+itable(i,2)-2,itable(i,1)+itable(i,3)-2,itable(i,2)+itable(i,3)-2,sx(i),sy(i),sz(i)
endif
enddo
enddo

jex=0
di=1
e_Prev=0.d0

!Inquire and if added to check if input file exists, if not create random spins (April 21 2016)

!Change file name to appropriate file name, or automate this perhaps?
Inquire(file='backUpSpins.dat', exist=backUp_file)

if(backUp_file) then
call spinitPre

else
Call random_number(rnd)
ran=INT(rnd*100)

Do i=1,ran
call spinit
EndDo

call spinit
rmeas=1
I2count=0
endif

DO I2=I2Count,ntemp1
if(rmeas.eq.1) then !!Only needed if using backup
sm2av=0.
entav=0.
mcav2=0.
mcav=0.
msav2=0.
msav=0.
Eav=0.
Eav2=0.
Eav4=0.
mav=0.
mav2=0.
mzav=0.
myav=0.
mxav=0.
rmav2=0.0
rmav=0.0
rmxav=0.0
rmyav=0.0
rmzav=0.0
rmav4=0.0
cmxav=0.
cmyav=0.
cmzav=0.
cm2av=0.
rm2root=0.
endif

T=Tmax1
!T=Tmax1-I2*(Tmax1-Tmin1)/(ntemp1)  !Cooling
!T=Tmin1+I2*(Tmax1-Tmin1)/ntemp1  !Heating
! run for ntrans iterations to eliminate transients
if(rmeas.eq.1) then

Inquire(file='tranSpins.dat', exist=tran_file)
tmeas=0

if(tran_file) then
    Open(unit=3,file='tranSpins.dat', Status="OLD")   !Open the spin configuration data.

    read(3,fmt=*) I2check, tmeas          !Read the position of the run.

    do i=1,nsite                    !Read in spins in the correct positions.
    do j=1,3                            !Determine correct order (Necessary to skip D sublattice.)
    if(isub(i).eq.j) then                   !IF sublattice matches
        read(3,fmt=*) sx(i),sy(i),sz(i)               !Read spins.
endif
    enddo
    enddo

flush(3)
Close(unit=3)
endif

Do jmeas=tmeas,ntrans
!call sweep
call metro
!call cluster

if(Modulo(jmeas,BACKUP_STEPS).eq.0.or.jmeas.eq.tmeas) then
open(unit=2, file='tranSpins.dat', status='REPLACE', action='WRITE')
write(2, fmt=*) I2, jmeas
do i=1,nsite
do j=1,3
if(isub(i).eq.j) then
write(2, fmt=*) sx(i),sy(i),sz(i)
endif
enddo
enddo
endif
close(unit=2)

!Output current state of the system for time series every so many iterations

if(Modulo(jmeas,TIME_SERIES_STEPS).eq.0) then         !Every x steps, output the data.

open(unit=2, file="timeSeries.dat",position="append",action='WRITE',status='unknown')    !Open time series data file.

write(2, fmt=*) I2, jmeas, E, mx*mx+my*my+mz*mz, smx(0)*smx(0)+smy(0)*smy(0)+smz(0)*smz(0), &
 smx(1)*smx(1)+smy(1)*smy(1)+smz(1)*smz(1), smx(2)*smx(2)+smy(2)*smy(2)+smz(2)*smz(2), &
  smx(3)*smx(3)+smy(3)*smy(3)+smz(3)*smz(3),rmx*rmx+rmy*rmy+rmz*rmz         !Write out all relevant data line by line.

flush(2)        !Flush and close the output.
close(unit=2)
endif

 enddo
endif

!Set up cluster
if(USE_CLUSTER) then
    call cluster_setup
endif

! perform nmeas measurements and average
Do imeas=rmeas,nmeas
!call sweep
if(USE_CLUSTER) then
    call RANDOM_NUMBER(rnd)
    if(rnd.le.P_CLUSTER) then
        call cluster
    else
        call metro
    endif
else
    call metro
endif

call measure

if(E.gt.e_Prev) then
    write(12480,*) imeas, E, e_Prev
endif

e_Prev=E

do i=0,3
sm2av(i)=sm2av(i)+smx(i)*smx(i)+smy(i)*smy(i)+smz(i)*smz(i)
sm2avalt=sm2avalt(i)+smxalt(i)*smxalt(i)+smyalt(i)*smyalt(i)+smzalt(i)*smzalt(i)
enddo

 mzav=mzav+mz
 mxav=mxav+mx
 myav=myav+my
 mav2=mav2+mx*mx+my*my+mz*mz
 EAV=Eav+e
 EAV2=EAV2+(E*E)
 Eav4=eav4+(E*E*E*E)

rmav2=rmav2+rmx*rmx+rmy*rmy+rmz*rmz
rmav4=rmav4+(rmx*rmx+rmy*rmy+rmz*rmz)**2
rmxav=rmxav+(rmx)
rmyav=rmyav+(rmy)
rmzav=rmzav+(rmz)
rm2root=rm2root+dsqrt(rmx*rmx+rmy*rmy+rmz*rmz)

!Output current state of the system for time series every so many iterations

if(Modulo(imeas,TIME_SERIES_STEPS).eq.0) then         !Every x steps, output the data.

open(unit=2, file="timeSeries.dat",position="append",action='WRITE',status='unknown')    !Open time series data file.

write(2, fmt=*) I2, imeas+ntrans, E, mx*mx+my*my+mz*mz, smx(0)*smx(0)+smy(0)*smy(0)+smz(0)*smz(0), &
 smx(1)*smx(1)+smy(1)*smy(1)+smz(1)*smz(1), smx(2)*smx(2)+smy(2)*smy(2)+smz(2)*smz(2), &
  smx(3)*smx(3)+smy(3)*smy(3)+smz(3)*smz(3),rmx*rmx+rmy*rmy+rmz*rmz         !Write out all relevant data line by line.

flush(2)        !Flush and close the output.
close(unit=2)
endif

!Output current temp and spin configuration every so many iterations in case a restart is needed (Apr 21 2016) !Keep backup now works for cooling and heating up to 499? I2 iterations
if(Modulo(imeas,BACKUP_STEPS).eq.0.or.imeas.eq.rmeas) then
open(unit=2, file="backUpSpins.dat", status='REPLACE', action='WRITE')
write(2, fmt=*) I2,imeas
do i=1,nsite
do j=1,3
if(isub(i).eq.j) then
write(2, fmt=*) sx(i),sy(i),sz(i)
endif
enddo
enddo

open(unit=15, file="backUpData.dat", status='REPLACE', action='WRITE')
write(15, fmt=*) mxav,myav,mzav,mav2
write(15, fmt=*) eav,eav2,eav4,entav
write(15, fmt=*) rmxav,rmyav,rmzav,rm2root,rmav2,rmav4
write(15, fmt=*) sm2av
write(15, fmt=*) sm2avalt
write(15, fmt=*) msav,msav2
write(15, fmt=*) mcav,mcav2
write(15, fmt=*) cmxav, cmyav, cmzav, cm2av
flush(2)
flush(15)
close(unit=2)
close(unit=15)
endif

enddo
sm2av=sm2av/nmeas
smav=dsqrt(sm2av)/nz1
sm2avalt=sm2avalt/nmeas
smavalt=dsqrt(sm2avalt)/nz1

rmav2=rmav2/nmeas
rmav4=rmav4/nmeas
rmxav=rmxav/nmeas
rmyav=rmyav/nmeas
rmzav=rmzav/nmeas
rm2root=rm2root/nmeas
rmav=dsqrt(rmav2)/(nocc)


chi=(rmav2-rm2root**2)/T
chi=chi/nocc

ub=1-rmav4/(3*rmav2*rmav2)
ub=rmav4/(rmav2*rmav2)


mzav=mzav/nmeas
myav=myav/nmeas
mxav=mxav/nmeas
 mav2=mav2/nmeas
 eav2=eav2/nmeas
 eav4=eav4/nmeas
 mcav2=mcav2/nmeas
 msav2=msav2/nmeas
 eav=eav/nmeas
 mzav=mzav/(nocc)
mxav=mxav/(nocc)
myav=myav/(nocc)
UE=1-Eav4/(3*Eav2*eav2)
 mav=sqrt(mav2)/(nocc)
 CVT=(EAV2-EAV*EAV)/(T*T)/(nocc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write to each file. Inquire and Open in case of restart.
!
! When restarted, carriage of each file is placed at beginning and must be set to append.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Open(unit=5,file='fort.77',access='append',action='write')
write(unit=5, fmt=999) t, ub,UE
close(unit=5)

Open(unit=5,file='fort.69',access='append',action='write')
write(unit=5, fmt=999) t, rmav
Close(unit=5)

Open(unit=5,file='fort.71',access='append',action='write')
write(unit=5, fmt=999) t, chi
close(unit=5)

Open(unit=5,file='fort.25',access='append',action='write')
write(unit=5,fmt=999) T,smav
close(unit=5)

Open(unit=5,file='fort.26',access='append',action='write')
 write(unit=5,fmt=999)T,EAV/(nocc),Cvt,mav
close(unit=5)

Open(unit=5,file='fort.27',access='append',action='write')
write(unit=5,fmt=999) T,smavalt
close(unit=5)

Open(unit=5,file='fort.28',access='append',action='write')
 write(unit=5,fmt=999) T,mxav,myav,mzav
close(unit=5)

Open(unit=5,file='fort.30',access='append',action='write')
 write(unit=5,fmt=999) T,rmxav
close(unit=5)

Open(unit=5,file='fort.31',access='append',action='write')
 write(unit=5,fmt=999) T,rmyav
close(unit=5)

Open(unit=5,file='fort.32',access='append',action='write')
 write(unit=5,fmt=999) T,rmzav
close(unit=5)

Open(unit=5,file='fort.40',access='append',action='write')
 write(unit=5,fmt=999) t1,t2,entav
close(unit=5)

Open(unit=5,file='fort.29',access='append',action='write')
 write(unit=5,fmt=999) T,smx
close(unit=5)

Open(unit=5,file='fort.33',access='append',action='write')
 write(unit=5,fmt=999) T,smy
close(unit=5)

Open(unit=5,file='fort.34',access='append',action='write')
 write(unit=5,fmt=999) T,smz
close(unit=5)

 999 format(1x,8e16.6)
 998 format(1x,i10,1x,8f16.6)

do i=1,nsite
do j=1,3
if(isub(i).eq.j) then
write(100+i2,*) sx(i),sy(i),sz(i)
endif
enddo
enddo

!Delete the tranSpins backUp when no longer necessary. Allows for new tranSpins to be created.
open(unit=1234, file='tranSpins.dat', status='old')
close(1234, status='delete')

rmeas=1

end do

991 format(1x,3F10.3)
990 format(1x,3I10,3F10.3)
99 format(1x,4I10,3F10.3)
97 format(1x,2I10,3F10.3)
end program kag_dipole_3d
!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!
subroutine spinit
use input_module_3d

! initialize all spins in the particle in random directions

            DO IT=1,nsite
               
                   
            CT=-1.+2.*r250_(iseed)
             st=sqrt(1.-ct*ct)
             PHI=2*pi*r250_(iseed)
sx(it)=sm(it)*st*cos(phi)
sy(it)=sm(it)*st*sin(phi)
sz(it)=sm(it)*ct
end do
end subroutine spinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read input file and use it to populate spins in specific directions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine spinitPre
use input_module_3d

!initialze all spins in the lattice according to an input file
Open(unit=3,file='backUpSpins.dat', Status="OLD")   !Open the spin configuration data.

read(3,fmt=*) I2count, rmeas          !Read the position of the run.

do i=1,nsite                    !Read in spins in the correct positions.
do j=1,3                            !Determine correct order (Necessary to skip D sublattice.)
if(isub(i).eq.j) then                   !IF sublattice matches
read(3,fmt=*) sx(i),sy(i),sz(i)               !Read spins.
endif
enddo
enddo

flush(3)
Close(unit=3)

Open(unit=4,file='backUpData.dat', Status="OLD")    !Open the simulation data

read(4,fmt=*) mxav,myav,mzav,mav2
read(4,fmt=*) eav,eav2,eav4,entav
read(4,fmt=*) rmxav,rmyav,rmzav,rm2root,rmav2,rmav4
read(4,fmt=*) sm2av
read(4,fmt=*) sm2avalt
read(4,fmt=*) msav,msav2
read(4,fmt=*) mcav,mcav2
read(4,fmt=*) cmxav, cmyav, cmzav, cm2av
flush(4)
close(unit=4)

end subroutine spinitPre

!r250_
!Produces random number

  double precision function r250_(idum)
      implicit double precision (a-h,o-z)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ntab,ndiv
 real*8 ran2_new,am,eps,rnmx
      parameter (im1=2147483563,im2=2147483399,am=1.d0/im1,imm1=im1-1,&
     &    ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211, &
     &     ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2d-14,rnmx=1.d0-eps)
!     
      integer idum2,j,k,iv(ntab),iy
      save    iv,iy,idum2
      data idum2/123456789/,iv/ntab*0/,iy/0/
      
      if(idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do j=ntab+8,1,-1
            k=idum/iq1
            idum=ia1*(idum-k*iq1)-k*ir1
            if(idum.lt.0) idum=idum+im1
            if(j.le.ntab) iv(j) = idum
         end do
         iy = iv(1)
      end if
!     
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
                if(idum.lt.0) idum=idum+im1
!     
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2=idum2+im2
      j = 1+iy/ndiv
      iy = iv(j)-idum2
      iv(j) = idum
      if(iy.lt.1) iy=iy+imm1 
      r250_ = min(am*iy,rnmx)
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!sweep
! Performs Metropolis algorithm iteratively over all spins

subroutine sweep
use input_module_3d

Do it=1,nsite
if(sm(it).eq.0) go to 55
hpx1=0.
 hpy1=0.
 hpz1=0.
do j=1,nsite
ix=(itable(j,1)-itable(it,1))
iyy=(itable(j,2)-itable(it,2))
iz=(itable(j,3)-itable(it,3))
if(ix<0) then 
ix=ix+L 
endif
if(iyy<0) then 
iyy=iyy+L 
endif
if(iz<0) then 
iz=iz+L 
endif
s(j,1)=sx(j)
s(j,2)=sy(j)
s(j,3)=sz(j)

do k=1,3
hpx1=hpx1-D(1,k,ix,iyy,iz)*s(j,k)
hpy1=hpy1-D(2,k,ix,iyy,iz)*s(j,k)
hpz1=hpz1-D(3,k,ix,iyy,iz)*s(j,k)
enddo 
 enddo

 hpx1=hpx1/t
 hpy1=hpy1/t
 hpz1=hpz1/t

 hmp=hpx1**2+hpy1**2+hpz1**2
 if(hmp.lt.1.d-06) go to 55
 hm=sqrt(hmp)
 hmi=1./hm

 hx=hpx1*hmi
 hy=hpy1*hmi
 hz=hpz1*hmi
 sth=sqrt(hx**2+hy**2)
 sthi=1./sth
 cph=hx*sthi
 sph=hy*sthi

 phi=2*pi*r250_(iseed)
 cp=cos(phi)
 sp=sin(phi)
 y=r250_(iseed)
 ct=1.0+log((1.0-y)*exp(-2.0*hm)+y)*hmi
 IF (ct.gt.1)ct=1.
 st=sqrt(abs(1.0-ct*ct))  
 sx(it)=(hz*cph*st*cp-sph*st*sp+hx*ct)
 sy(it)=(hz*sph*st*cp+cph*st*sp+hy*ct)
 sz(it)=(-sth*st*cp+hz*ct)
 55 continue
 
  end do
  return
  end subroutine sweep
  !!!!!!!!!!!!!!!!!!!!!!!
  !Measure
  !Measures observables of system

  subroutine measure
use input_module_3d
 
  E=0
  mx=0.
  my=0.
  mz=0.
  mcx=0.
  mcy=0.
  mcz=0.
  msx=0.
  msy=0.
  msz=0.
  smx=0.
smy=0.
smz=0.
smxalt=0.
smyalt=0.
smzalt=0.
rmx=0.0
rmy=0.0
rmz=0.0
  Do IT=1,Nsite
  do i=0,3
  if(isub(it).eq.i) then
  smx(i)=smx(i)+sx(it)*(-1)**(itable(it,1)+itable(it,2)+itable(it,3))
  smy(i)=smy(i)+sy(it)*(-1)**(itable(it,1)+itable(it,2)+itable(it,3))
  smz(i)=smz(i)+sz(it)*(-1)**(itable(it,1)+itable(it,2)+itable(it,3))
  smxalt(i)=smxalt(i)+sx(it)
  smyalt(i)=smyalt(i)+sy(it)
  smzalt(i)=smzalt(i)+sz(it)
  endif
enddo
rmx=rmx+sx(it)
rmy=rmy+sy(it)
rmz=rmz+sz(it)

  mx=mx+sx(it)
  my=my+sy(it)
  mz=mz+sz(it)
 

hpx1=0.d0
hpy1=0.d0
hpz1=0.


do j=1,nsite
ix=(itable(j,1)-itable(it,1))
iyy=(itable(j,2)-itable(it,2))
iz=(itable(j,3)-itable(it,3))
if(ix<0) then 
ix=ix+L 
endif
if(iyy<0) then 
iyy=iyy+L 
endif
if(iz<0) then 
iz=iz+L 
endif
s(j,1)=sx(j)
s(j,2)=sy(j)
s(j,3)=sz(j)

do k=1,3
hpx1=hpx1-D(1,k,ix,iyy,iz)*s(j,k)
hpy1=hpy1-D(2,k,ix,iyy,iz)*s(j,k)
hpz1=hpz1-D(3,k,ix,iyy,iz)*s(j,k)
enddo 


 enddo
 E=E-hpx1*sx(it)-hpy1*sy(it)-hpz1*sz(it)
 enddo
E=E/2
 return
 end subroutine measure


!metro
! Perform metropolis algorithm by randomly selecting spins

  subroutine metro
  use input_module_3d

real eold,enew,sx1,sy1,sz1,sdot1,de,stnew,seed

Do iii=1,nsite
55 continue

is=INT(nsite*r250_(iseed))
if(is.lt.1.or.is.gt.nsite) go to 55
if(sm(is).eq.0) go to 55
it=is

hpx1=0.
 hpy1=0.
 hpz1=0.

do j=1,nsite
ix=(itable(j,1)-itable(it,1))
iyy=(itable(j,2)-itable(it,2))
iz=(itable(j,3)-itable(it,3))
if(ix<0) then 
ix=ix+L 
endif
if(iyy<0) then 
iyy=iyy+L 
endif
if(iz<0) then 
iz=iz+L 
endif
s(j,1)=sx(j)
s(j,2)=sy(j)
s(j,3)=sz(j)

do k=1,3
hpx1=hpx1-D(1,k,ix,iyy,iz)*s(j,k)
hpy1=hpy1-D(2,k,ix,iyy,iz)*s(j,k)
hpz1=hpz1-D(3,k,ix,iyy,iz)*s(j,k)
enddo 
 enddo
  
 eold=-(hpx1*sx(it)+hpy1*sy(it)+hpz1*sz(it))
 
  phi=2*pi*r250_(iseed)
  cp=cos(phi)
 sp=sin(phi)
 ct=-1.+2*r250_(iseed)
 st=sqrt(abs(1.0-ct*ct))
 sx1=st*cp
 sy1=st*sp
 sz1=ct

stnew=sqrt(sx1*sx1+sy1*sy1+sz1*sz1)

sx1=sx1/stnew
sy1=sy1/stnew
sz1=sz1/stnew


 enew= -(hpx1*sx1+hpy1*sy1+hpz1*sz1)

DE=Enew-Eold
                 y=r250_(iseed)
                 If(y.le.(exp(-de/t)/(1.+exp(-de/t)))) then
                  sx(it)=sx1
                 sy(it)=sy1
                 sz(it)=sz1
                                  endif 
end do
  return
  end subroutine metro

subroutine cluster_setup
use input_module_3d
use cluster_module

!Build J0 and R0 arrays, calling sub_I=1 as origin spin
do sub_J=1,nsite

    !Set self interaction to 0
    Jij(sub_J,sub_J)=0.d0

    do sub_K=sub_J+1,nsite

        !Periodic boundary condition must be satisfied. (If multiple equal shortest paths Jij should be 0, how do?)
        Rx(sub_J,sub_K)=(ITABLE(sub_K,1)-ITABLE(sub_J,1))
        if(Rx(sub_J,sub_K).lt.(-1)*L/2) Rx(sub_J,sub_K)=Rx(sub_J,sub_K)+L
        if(Rx(sub_J,sub_K).ge.L/2) Rx(sub_J,sub_K)=Rx(sub_J,sub_K)-L

        Ry(sub_J,sub_K)=(ITABLE(sub_K,2)-ITABLE(sub_J,2))
        if(Ry(sub_J,sub_K).lt.(-1)*L/2) Ry(sub_J,sub_K)=Ry(sub_J,sub_K)+L
        if(Ry(sub_J,sub_K).ge.L/2) Ry(sub_J,sub_K)=Ry(sub_J,sub_K)-L

        Rz(sub_J,sub_K)=(ITABLE(sub_K,3)-ITABLE(sub_J,3))
        if(Rz(sub_J,sub_K).lt.(-1)*L/2) Rz(sub_J,sub_K)=Rz(sub_J,sub_K)+L
        if(Rz(sub_J,sub_K).ge.L/2) Rz(sub_J,sub_K)=Rz(sub_J,sub_K)-L

        Rx(sub_K,sub_J)=(-1)*Rx(sub_J,sub_K)
        Ry(sub_K,sub_J)=(-1)*Ry(sub_J,sub_K)
        Rz(sub_K,sub_J)=(-1)*Rz(sub_J,sub_K)

        Rij(sub_J,sub_K)=(Rx(sub_J,sub_K)*Rx(sub_J,sub_K))+(Ry(sub_J,sub_K)*Ry(sub_J,sub_K))+(Rz(sub_J,sub_K)*Rz(sub_J,sub_K))
        Rij(sub_J,sub_K)=dsqrt(Rij(sub_J,sub_K))

        Rij(sub_K,sub_J)=Rij(sub_J,sub_K)

        !Interaction strength between any two
        if(Rx(sub_J,sub_K).eq.(L/2) .OR. Ry(sub_J,sub_K).eq.(L/2) .OR. Rz(sub_J,sub_K).eq.(L/2)) then   !Set interaction to zero if two equal displacements using PBC

        Jij(sub_J,sub_K)=0.d0 
        Jij(sub_K,sub_J)=0.d0

        else if(isub(sub_J).eq.0 .OR. isub(sub_K).eq.0) then

        Jij(sub_J,sub_K)=0.d0 
        Jij(sub_K,sub_J)=0.d0

        else

        Jij(sub_J,sub_K)=Js/(Rij(sub_J,sub_K)*Rij(sub_J,sub_K)*Rij(sub_J,sub_K))
        Jij(sub_K,sub_J)=Js/(Rij(sub_J,sub_K)*Rij(sub_J,sub_K)*Rij(sub_J,sub_K))

        endif

        if(imeas.eq.1) then
            write(12321,*) sub_K,sub_J,Rij(sub_J,sub_K),Jij(sub_J,sub_K)
        endif

    enddo ! sub_K
enddo ! sub_J


!Build an array of partial sums

!Sum of every point up to the current point (See Baek)
s0_Cum(0:3,1)=(/ 0.d0, 0.d0, 0.d0, 0.d0 /)    !Set first spin to 0 for each sublattice
s0_Start(0:3)=(/ inv_isub(0,1), inv_isub(1,1), inv_isub(2,1), inv_isub(3,1) /)

write(11122,*) 0,0,0, s0_Start(1), s0_Start(2), s0_Start(3), s0_Start(0)

do sub_I=0,3

    s0_Pos(sub_I,1:3)=(/ &
        ITABLE(s0_Start(sub_I),1),ITABLE(s0_Start(sub_I),2),ITABLE(s0_Start(sub_I),3) &
    /)

enddo ! sub_I

write(21131,*) sub_I,sub_J,sub_K,s0_Pos(1,1:3)
write(21132,*) sub_I,sub_J,sub_K,s0_Pos(2,1:3)
write(21133,*) sub_I,sub_J,sub_K,s0_Pos(3,1:3)
write(21134,*) sub_I,sub_J,sub_K,s0_Pos(0,1:3)

cc1=2

do sub_I=0,L-1
    do sub_J=0,L-1
        do sub_K=0,L-1

            if(sub_I.eq.0 .AND. sub_J.eq.0 .AND. sub_K.eq.0) go to 111

            do cc2=0,3

                s0_Rel(cc2,1:3)=(/ &
merge(s0_Pos(cc2,1)+sub_I,s0_Pos(cc2,1)+sub_I-L,s0_Pos(cc2,1)+sub_I.le.L), &
merge(s0_Pos(cc2,2)+sub_J,s0_Pos(cc2,2)+sub_J-L,s0_Pos(cc2,2)+sub_J.le.L), &
merge(s0_Pos(cc2,3)+sub_K,s0_Pos(cc2,3)+sub_K-L,s0_Pos(cc2,3)+sub_K.le.L) &
                /)
            
            enddo ! cc2

            write(11131,*) sub_I,sub_J,sub_K,s0_Rel(1,1:3)
            write(11132,*) sub_I,sub_J,sub_K,s0_Rel(2,1:3)
            write(11133,*) sub_I,sub_J,sub_K,s0_Rel(3,1:3)
            write(11134,*) sub_I,sub_J,sub_K,s0_Rel(0,1:3)

            cs1=TABLE(s0_Rel(1,1),s0_Rel(1,2),s0_Rel(1,3))
            cs2=TABLE(s0_Rel(2,1),s0_Rel(2,2),s0_Rel(2,3))
            cs3=TABLE(s0_Rel(3,1),s0_Rel(3,2),s0_Rel(3,3))
            cs4=TABLE(s0_Rel(0,1),s0_Rel(0,2),s0_Rel(0,3))
   
            s0_Cum(1,cc1)=s0_Cum(1,cc1-1)+(Jij(s0_Start(1),cs1)/(kb*T)) ! Make a cumulative sum for each sublattice
            s0_Cum(2,cc1)=s0_Cum(2,cc1-1)+(Jij(s0_Start(2),cs2)/(kb*T))
            s0_Cum(3,cc1)=s0_Cum(3,cc1-1)+(Jij(s0_Start(3),cs3)/(kb*T))
            s0_Cum(0,cc1)=0.d0    ! Non magnetic always zero

            write(11122,*) sub_I,sub_J,sub_K,cs1,cs2,cs3,cs4

            cc1=cc1+1 !   Increment the counter

            111 continue

        enddo ! sub_K
    enddo ! sub_J
enddo ! sub_I

end subroutine cluster_setup

!cluster
!Test cluster algorithm as defined by Baek.

subroutine cluster
use stack_module
use input_module_3d
use cluster_module

!cluster_Spin all sites as false. Change cluster_Spin to true when thespin_Curry are added to the stack/cluster
!cluster_Spin sites where isub is 0 to remove unoccupied sites from being randomly picked
do sub_I=1,nsite
    if(isub(sub_I).eq.0) then
    cluster_Spin(sub_I)=.true.
    else 
    cluster_Spin(sub_I)=.false.
    endif
enddo ! sub_I

!Pick randomly a spin and add it's index to the stack.
123 continue

!spin_Curr=int(r250_(jseed)*nsite)
call random_number(rand)
spin_Curr=int(rand*nsite)

!Make sure it is a acceptable site (within 1:nsite and not isub=0)
if(spin_Curr.gt.nsite.or.spin_Curr.lt.1) go to 123
if(cluster_Spin(spin_Curr)) go to 123

call push(spin_Curr,cluster_Stack)
cluster_Spin(spin_Curr)=.true.

!Randomly pick a direction vector ref_Norm to define a plane of reflection
ref_Norm_L=0.0

do sub_K=1,3

    !ref_Norm(sub_K)=r250_(jseed)
    call random_number(rand)
    ref_Norm(sub_K)=rand
    ref_Norm_L=ref_Norm_L+ref_Norm(sub_K)*ref_Norm(sub_K)

enddo ! sub_K

ref_Norm_L=dsqrt(ref_Norm_L)

!Normalize the direction vector ref_Norm
do sub_K=1,3

    ref_Norm(sub_K)=ref_Norm(sub_K)/ref_Norm_L

enddo   ! sub_K

counterrr=0

!Continue while the stack is not empty
do while(.not. isEmpty(cluster_Stack))

    call pop(spin_Curr,cluster_Stack)
    counterrr=counterrr+1

    sub_Curr=isub(spin_Curr)    !   First spin sublattice

    curr_Pos(1:3)=(/ &
        ITABLE(spin_Curr,1),ITABLE(spin_Curr,2),ITABLE(spin_Curr,3) &
    /)

    !Find reflection of spin spin_Curr
    proj=(sx(spin_Curr)*ref_Norm(1)+sy(spin_Curr)*ref_Norm(2)+sz(spin_Curr)*ref_Norm(3))

    sx_New(spin_Curr)=sx(spin_Curr)-2*((proj)*ref_Norm(1))/ref_Norm_L
    sy_New(spin_Curr)=sy(spin_Curr)-2*((proj)*ref_Norm(2))/ref_Norm_L
    sz_New(spin_Curr)=sz(spin_Curr)-2*((proj)*ref_Norm(3))/ref_Norm_L

    s_New_L=sx_New(spin_Curr)*sx_New(spin_Curr)+sy_New(spin_Curr)*sy_New(spin_Curr)+sz_New(spin_Curr)*sz_New(spin_Curr)
    s_New_L=dsqrt(s_New_L)

    sx_New(spin_Curr)=sx_New(spin_Curr)/s_New_L
    sy_New(spin_Curr)=sy_New(spin_Curr)/s_New_L
    sz_New(spin_Curr)=sz_New(spin_Curr)/s_New_L

    zc=1    ! Set zc to next spin position (i.e. 0 ahead)
    wc=2
    
    122 continue

    !Draw a uniform random number uc
    !uc=r250_(jseed)
    call random_number(uc)

    if(uc.ge.1) go to 122

    !Find an index that satisfies the inequality
    do while(wc.le.nsite)

        wc_Pos(1:3)=(/ &    !Find the actual position of the new spin
merge(curr_Pos(1)+((wc-1)/(L**2)),curr_Pos(1)+((wc-1)/(L**2))-L,curr_Pos(1)+((wc-1)/(L**2)).le.L), &
merge(curr_Pos(2)+(MODULO(wc-1,L**2)/L),curr_Pos(2)+(MODULO(wc-1,L**2)/L)-L,curr_Pos(2) &
    +(MODULO(wc-1,L**2)/L).le.L), &
merge(curr_Pos(3)+MODULO((MODULO(wc-1,L**2)),L),curr_Pos(3)+MODULO((MODULO(wc-1,L**2)),L)-L, &
    curr_Pos(3)+MODULO((MODULO(wc-1,L**2)),L).le.L) &
        /)

        wc_Prime=TABLE(wc_Pos(1),wc_Pos(2),wc_Pos(3))    !Get the spin number corresponding with that position

        if(imeas.eq.1) then
            write(18882,*) spin_Curr, wc, wc_Prime, curr_Pos(1:3), wc_Pos(1:3)
        endif

        if(cluster_Spin(wc_Prime) .AND. isub(wc_Prime).ne.0) go to 125

        if(imeas.eq.1) then
            write(12348,*) spin_Curr,zc,wc-1,wc,s0_Cum(sub_Curr,wc-1),(-1)*dlog(1.0-uc)/(2)+s0_Cum(sub_Curr,zc),s0_Cum(sub_Curr,wc)
        endif

        if(s0_Cum(sub_Curr,wc-1).le.(-1)*(dlog(1.0-uc)/(2))+s0_Cum(sub_Curr,zc).and. &
    s0_Cum(sub_Curr,wc).gt.(-1)*(dlog(1.0-uc)/(2))+s0_Cum(sub_Curr,zc)) then
        zc=wc

        !Find probability to add new spin

        expo1_1x=Jij(spin_Curr,wc_Prime)*(sx_New(spin_Curr)-sx(spin_Curr))
        expo1_1y=Jij(spin_Curr,wc_Prime)*(sy_New(spin_Curr)-sy(spin_Curr))
        expo1_1z=Jij(spin_Curr,wc_Prime)*(sz_New(spin_Curr)-sz(spin_Curr))

        expo1_2x=sx(wc_Prime)/(kb*T)
        expo1_2y=sy(wc_Prime)/(kb*T)
        expo1_2z=sz(wc_Prime)/(kb*T)

        expo1=expo1_1x*expo1_2x+expo1_1y*expo1_2y+expo1_1z*expo1_2z
        expo2=(-2)*(Jij(spin_Curr,wc_Prime))/(kb*T)

        P_Add=max(0.d0,(1-dexp(expo1))/(1-dexp(expo2)))

        !if(r250_(jseed).le.P_Add) then
        call random_number(rand)
        if(rand.le.P_Add) then

        cluster_Spin(wc_Prime)=.true.
        call push(wc_Prime,cluster_Stack)

        endif

        124 continue

        !Draw a uniform random number uc
        !uc=r250_(jseed)
        call random_number(uc)

        if(uc.ge.1) go to 124

        endif

        125 continue
        wc=wc+1
    enddo ! wc
enddo ! Stack empty

if(imeas.lt.10000) then
    write(12346,*) imeas,counterrr
endif

!Calculate bulk_Curr and surface energy change

!bulk_Curr
dE_Bulk=0.d0
do sub_I=1,nsite-1

    if(cluster_Spin(sub_I).and.isub(sub_I).ne.0) then
    do sub_K=sub_I+1,nsite

        if(cluster_Spin(sub_K).and.isub(sub_K).ne.0)then

        bulk_Curr=0.d0

        Rx_Curr=Rx(sub_I,sub_K)
        Ry_Curr=Ry(sub_I,sub_K)
        Rz_Curr=Rz(sub_I,sub_K)

        !Use relative position vector Rij
        bulk_Curr=3*(sx(sub_I)*Rx_Curr+sy(sub_I)*Ry_Curr+sz(sub_I)*Rz_Curr)*(sx(sub_K)*Rx_Curr+sy(sub_K)*Ry_Curr+sz(sub_K)*Rz_Curr)
        bulk_Curr=bulk_Curr-3*(sx_New(sub_I)*Rx_Curr+sy_New(sub_I)*Ry_Curr+sz_New(sub_I)*Rz_Curr)*(sx_New(sub_K)*Rx_Curr &
    +sy_New(sub_K)*Ry_Curr+sz_New(sub_K)*Rz_Curr)
        bulk_Curr=bulk_Curr/(Rij(sub_I,sub_K)**3)

        dE_Bulk=dE_Bulk+bulk_Curr
        endif
    enddo ! sub_K
    endif
enddo ! sub_I

dE_Bulk=Js*dE_Bulk

!surface
dE_Surf=0.0
do sub_I=1,nsite

    if(cluster_Spin(sub_I).and.isub(sub_I).ne.0) then
    do sub_K=1,nsite

        if(.not. cluster_Spin(sub_K).and.isub(sub_I).ne.0) then

        surf_Curr=0.d0

        Rx_Curr=Rx(sub_I,sub_K)
        Ry_Curr=Ry(sub_I,sub_K)
        Rz_Curr=Rz(sub_I,sub_K)

        surf_Curr=3*(sx(sub_I)*Rx_Curr+sy(sub_I)*Ry_Curr+sz(sub_I)*Rz_Curr)*(sx(sub_K)*Rx_Curr+sy(sub_K)*Ry_Curr+sz(sub_K)*Rz_Curr)
        surf_Curr=surf_Curr-3*(sx_New(sub_I)*Rx_Curr+sy_New(sub_I)*Ry_Curr+sz_New(sub_I)*Rz_Curr)*(sx(sub_K)*Rx_Curr &
    +sy(sub_K)*Ry_Curr+sz(sub_K)*Rz_Curr)
        surf_Curr=surf_Curr/(Rij(sub_I,sub_K)**3)

        dE_Surf=dE_Surf+surf_Curr

        endif
    enddo   ! sub_K
    endif
enddo   ! sub_I

dE_Surf=Js*dE_Surf

!Flip the cluster with probability P_Acc

P_Acc=min(1.0,dexp((-dE_Surf-dE_Bulk)/(kb*T)))

counterrr=0
!if(r250_(jseed).le.P_Acc) then
call random_number(rand)
if(rand.le.P_Acc) then
write(12399,*) imeas,dE_Surf,dE_Bulk,dE_Surf+dE_Bulk,rand,P_Acc,dexp(-(dE_Surf+dE_Bulk))

do sub_I=1,nsite

    if(cluster_Spin(sub_I).and.isub(sub_I).ne.0) then
    sx(sub_I)=sx_New(sub_I)
    sy(sub_I)=sy_New(sub_I)
    sz(sub_I)=sz_New(sub_I)
    counterrr=counterrr+1
    endif

enddo   ! sub_I

!Output cluster size flipped
write(12345,*) imeas,counterrr

endif

end subroutine cluster

!Gamma
!Build dipole interaction, etc.

Subroutine Gamma(L,D)
double precision :: a(3),b(3),c(3),u1(3),u2(3),u3(3),u2xu3(3),u3xu1(3),u1xu2(3),b1(3),b2(3),b3(3),vu
double precision :: Q(3),r(3),R1(3),Q2,Q1,rdotQ,rRab,rR2,rR1,rR4,rR5,F,G,D(3,3,0:L-1,0:L-1,0:L-1)
double precision, parameter :: pi=3.14159265358979323, sqrt_pi=1.772453850905516027,tau=1.d-3,eta=0.5d0  ,eta2=eta*eta
integer, parameter :: nsize=15
integer :: R3(3)
a(1)=1
a(2)=0
a(3)=0
b(1)=0.5d0
b(2)=dsqrt(3.d0)/2
b(3)=0
c(1)=0.5d0
c(2)=dsqrt(3.d0)/6
c(3)=dsqrt(2.d0)/dsqrt(3.d0)
a(1)=0
a(2)=1/dsqrt(2.d0)
a(3)=1/dsqrt(2.d0)
b(1)=1/dsqrt(2.d0)
b(2)=0
b(3)=1/dsqrt(2.d0)
c(1)=1/dsqrt(2.d0)
c(2)=1/dsqrt(2.d0)
c(3)=0

do i=1,3
u1(i)=L*a(i)
u2(i)=L*b(i)
u3(i)=L*c(i)
enddo
u2xu3(1)=u2(2)*u3(3)-u2(3)*u3(2)
u2xu3(2)=u2(3)*u3(1)-u2(1)*u3(3)
u2xu3(3)=u2(1)*u3(2)-u2(2)*u3(1)
vu=0
do i=1,3
vu=vu+u1(i)*u2xu3(i)
enddo
vu=dabs(vu)
u3xu1(1)=u3(2)*u1(3)-u3(3)*u1(2)
u3xu1(2)=u3(3)*u1(1)-u3(1)*u1(3)
u3xu1(3)=u3(1)*u1(2)-u3(2)*u1(1)
u1xu2(1)=u1(2)*u2(3)-u1(3)*u2(2)
u1xu2(2)=u1(3)*u2(1)-u1(1)*u2(3)
u1xu2(3)=u1(1)*u2(2)-u1(2)*u2(1)
do i=1,3
b1(i)=2*pi*u2xu3(i)/vu
b2(i)=2*pi*u3xu1(i)/vu
b3(i)=2*pi*u1xu2(i)/vu
enddo
do ial=1,3
do ibet=ial,3

do kk=0,L-1
r3(3)=kk
do jj=0,L-1
r3(2)=jj
do ii=0,L-1
r3(1)=ii
F=0
G=0

do i=1,3
r(i)=r3(1)*a(i)+r3(2)*b(i)+r3(3)*c(i)
enddo

do n=-nsize,nsize
do m=-nsize,nsize
do k=-nsize,nsize
do i=1,3
Q(i)=n*b1(i)+m*b2(i)+k*b3(i)
R1(i)=n*u1(i)+m*u2(i)+k*u3(i)+r(i)
enddo
Q2=Q(1)*Q(1)+Q(2)*Q(2)+Q(3)*Q(3)
Q1=dsqrt(Q2)
rdotQ=r(1)*Q(1)+r(2)*Q(2)+r(3)*Q(3)
if (Q1.lt.1.d-05) go to 50
if(kk.eq.0.and.jj.eq.0.and.ii.eq.0) then
F=F-4*pi*dexp(-Q2/(4*eta2))*(Q2*ndelta(ial,ibet)-3*Q(ial)*Q(ibet))/(3*vu*Q2)
else
F=F+4*pi*Q(ial)*Q(ibet)*dcos(rdotQ)*dexp(-Q2/(4*eta2))/(Q2*vu)
endif
50 continue 
rRab=R1(ial)*R1(ibet)
rR2=R1(1)*R1(1)+R1(2)*R1(2)+R1(3)*R1(3)
rR1=dsqrt(rR2)
rR4=rR2*rR2
rR5=rR4*rR1
if(rR1.lt.1.d-5) go to 51
if(kk.eq.0.and.jj.eq.0.and.ii.eq.0) then
G=G+(rR2*ndelta(ial,ibet)-3*R1(ial)*R1(ibet))*(2*eta*dexp(-rR2*eta2)*(3+2*rR2*eta2)/(3*rR4*sqrt_pi)+erfc(rR1*eta)/rR5)
else
G=G+2*eta*dexp(-eta2*rR2)*(rR2*ndelta(ial,ibet)-3*rRab-2*eta2*rRab*rR2)/(sqrt_pi*rR4)+(rR2*ndelta(ial,ibet)&
&-3*rRab)*derfc(eta*rR1)/rR5
endif
51 continue
enddo
enddo
enddo
69 continue
write(30+10*ial+ibet,*) r3,F+G
D(ial,ibet,ii,jj,kk)=F+G
D(ibet,ial,ii,jj,kk)=F+G
enddo
enddo
enddo

enddo
enddo

end subroutine Gamma
integer function ndelta(i,j)
integer i,j
ndelta=0
if(i.eq.j) ndelta=1
end

!init_random_seed
! initiate a random number seed.

Subroutine init_random_seed()
         Integer :: i, n, clock
         Integer, Dimension(:), Allocatable :: seed

         Call random_seed(size = n)
         Allocate(seed(n))

         Call system_clock(Count=clock)

         seed = clock + 37 * (/ (i - 1, i = 1, n) /)
         Call random_seed(PUT = seed)

         Deallocate(seed)
End Subroutine