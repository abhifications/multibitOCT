Program HTST

Implicit none

integer,parameter         :: r8b=kind(1.D200)! selected_real_kind(15,307)!kind(1.D200)
integer,parameter         :: i4b= kind(2147483647)


real(r8b),parameter::Mub = 0.057883818060d0 !meV/T
real(r8b),parameter::gama = 176.085963023d0 !1/(ns*T)
real(r8b),parameter::Pi = 3.14159265358979323846264338327950d0


real(r8b)::K_hard, K_easy
real(r8b)::Anis_vec_H(3),Anis_vec_E(3), Mag_stable(3), Mag_saddle(3), unit_tangent(3)
real(r8b)::Hessian(3,3),Phys_Hess_E(2,2),V(3,2),eig(2), eig_stable(2), eig_saddle(2)
real(r8b)::lambda, energy_barrier, thermal_energy, tau, magnetic_moment



! Input related variables
character(len=100) :: buffer, label,check,Images,Initial
integer(i4b) :: pos 
integer(i4b), parameter :: fh = 15
integer(i4b) :: ios = 0
integer(i4b) :: line = 0
logical::dir_e

real(r8b),external         :: rlfsr113
external                      lfsrinit
call lfsrinit(10)




call read_input_file




call Hessian_Analy(Hessian)
call Hessian_Phys_E(Mag_stable,Hessian,V,Phys_Hess_E,eig)
eig_stable(:) = eig(:)



call Hessian_Analy(Hessian)
call Hessian_Phys_E(Mag_saddle,Hessian,V,Phys_Hess_E,eig)
eig_saddle(:) = eig(:)

!print*,eig_stable
!print*,eig_saddle

call dynamical_factor(Phys_Hess_E,V,lambda)

call lifetime(eig_stable, eig_saddle, lambda, tau)

print*, "Mean lifetime in second(s) is: ", tau*1d-9






contains




subroutine read_input_file

implicit none
real(r8b):: norm,tmp

!-------check if file exist!------------------------

  inquire( file = "Input_File.in" , exist=dir_e )

  if ( dir_e ) then
     goto 7
  else
     print*,"File is Missing, please check Input_File.in"
     stop
  end if
!---------------------------------------------------

 7 open(fh, file='Input_File.in')

  ! ios is negative if an end of record condition is encountered or if
  ! an endfile condition was detected.  It is positive if an error was
  ! detected.  ios is zero otherwise.

  do while (ios == 0)
  
     read(fh, '(A)', iostat=ios) buffer
     if (ios == 0) then
        line = line + 1

        ! Find the first instance of whitespace.  Split label and data.
        pos = scan(buffer, '         ')
        label = buffer(1:pos)
        buffer = buffer(pos+1:)

        select case (label)

          case ('Mag_stable')
          read(buffer, *, iostat=ios)  Mag_stable 
          
          case ('Mag_saddle')
          read(buffer, *, iostat=ios)  Mag_saddle 
          
          
          case ('unit_tangent')
          read(buffer, *, iostat=ios)  unit_tangent
          
          case ('Anis_vec_E')
          read(buffer, *, iostat=ios)  Anis_vec_E   
          
          case ('Anis_vec_H')
          read(buffer, *, iostat=ios)  Anis_vec_H
          
          case ('K_easy')
          read(buffer, *, iostat=ios)  K_easy  
          
          case ('K_hard')
          read(buffer, *, iostat=ios)  K_hard 
          
          case ('magnetic_moment')
          read(buffer, *, iostat=ios)  magnetic_moment  
                       
           case ('energy_barrier')
          read(buffer, *, iostat=ios)  energy_barrier  
          
          case ('thermal_energy')
          read(buffer, *, iostat=ios)  thermal_energy  
          
                 
        end select
     end if



  end do


  norm = dsqrt( Anis_vec_E(1)*Anis_vec_E(1) + Anis_vec_E(2)*Anis_vec_E(2) + Anis_vec_E(3)*Anis_vec_E(3) )
  Anis_vec_E(1) = Anis_vec_E(1)/norm
  Anis_vec_E(2) = Anis_vec_E(2)/norm
  Anis_vec_E(3) = Anis_vec_E(3)/norm


  norm = dsqrt( Anis_vec_H(1)*Anis_vec_H(1) + Anis_vec_H(2)*Anis_vec_H(2) + Anis_vec_H(3)*Anis_vec_H(3) )
  Anis_vec_H(1) = Anis_vec_H(1)/norm
  Anis_vec_H(2) = Anis_vec_H(2)/norm
  Anis_vec_H(3) = Anis_vec_H(3)/norm


  !--if |M|>1 then we renormalize--
  norm = dsqrt (dot_product(Mag_stable,Mag_stable))
  if(norm .gt. 1) then
    tmp = normv(3,Mag_stable)
    Mag_stable(:) = Mag_stable(:)/tmp
  endif  

  !--if |M|>1 then we renormalize--
  norm = dsqrt (dot_product(Mag_saddle,Mag_saddle))
  if(norm .gt. 1) then
    tmp = normv(3,Mag_saddle)
    Mag_saddle(:) = Mag_saddle(:)/tmp
  endif 

end subroutine



!**************************************************************************
!                  Subroutine to Calculate Force / Effective field
!**************************************************************************


subroutine dF(Mag,F)

 implicit none

 real(r8b),intent(in)::Mag(3)
 real(r8b),intent(out)::F(3)
 real(r8b) :: ani1,ani2

 ani1 = dot_product(Mag,Anis_vec_E)
 ani2 = dot_product(Mag,Anis_vec_H)

 F(:) =  2.0d0*( K_easy*ani1*Anis_vec_E(:) +  K_hard*ani2*Anis_vec_H(:) )

end subroutine



!**************************************************************************
!                     Function for Normalization
!**************************************************************************
real(r8b) function normv(N,v)
implicit none
integer(i4b), intent(in) :: N
real(r8b), intent(in) :: v(N)
integer(i4b) :: i
real(r8b) :: s
   s = 0.0d0
   do i = 1,N
      s = s + v(i)*v(i)
   end do
   normv = dsqrt(s)
   
end function normv
!------------------------------------------------------------------------------



!**************************************************************************
!                    subroutine  for cross product
!**************************************************************************
subroutine cross_product(N,A,B,C)
implicit none
integer(i4b), intent(in) :: N
real(r8b), intent(in) :: A(N),B(N)
real(r8b),intent(out)::c(3)

  C(1) = A(2)*B(3) - A(3)*B(2)
  C(2) = -( A(1)*B(3) - A(3)*B(1) )
  C(3) = A(1)*B(2) - A(2)*B(1)
      

end subroutine
!------------------------------------------------------------------------------




!**************************************************************************
! Hessian using Analytical Calculation/second derivative of the energy
!**************************************************************************

 ! Energy of the system is defined as follows: E = -Sum_k(S.e_k)^2; and 
 ! to have easy axis in z-direction K should take +ve values, whereas 
 ! Kh should take -ve values
  
subroutine Hessian_Analy(Hessian)

 implicit none

 real(r8b),intent(out)::Hessian(3,3)
 integer(i4b)::alpha,beta


    18 format(1000(es16.8E3,x))          ! 1000 means # of columns
    open(20,file='Hessian.txt',action='write')
    101 format(1000(A,x))   
    write(20,101)'==============================================================='



    Hessian(:,:)  = 0.0d0


    do alpha = 1 , 3
  
       do beta = 1 , 3

          !------case for Anistropy Term-----------
          Hessian(alpha,beta) = -2.0d0*(K_hard*Anis_vec_H(alpha)*Anis_vec_H(beta) + K_easy*Anis_vec_E(alpha)*Anis_vec_E(beta))

       enddo     

       write(20,18) Hessian(alpha,:)
    enddo


    close(20)

end subroutine





!**************************************************************************
!                         Physical Hessian_Energy 
!**************************************************************************

subroutine Hessian_Phys_E(Mag,Hessian,V,Phys_Hess_E,eig)

 implicit none

 real(r8b),intent(in)::Mag(3)
 real(r8b),intent(inout)::Hessian(3,3)
 real(r8b),intent(out):: V(3,2),Phys_Hess_E(2,2),eig(2)
 real(r8b)::F(3),H_mult_V(3,2),eta(3),Xi(3)
 real(r8b)::B_dot_M,Lagrange_multiplier,tmp
 integer(i4b)::alpha,beta,info


12 format(1000(es16.8E3,x))          ! 1000 means # of columns
open(10,file='Physical_Hessian.txt',action='write')
102 format(1000(A,x))   
write(10,102)'======================================================='




!-----Lagrange Multipliers calculation-------

call dF(Mag,F)

Lagrange_multiplier = dot_product(F,Mag)
 


!-------Physical Hessian (3N*3N matrix)------
! here we add lagrange multipliers to the diagonal 
! elements of the hessian matrix

do alpha = 1 , 3
   Hessian(alpha,alpha) = Hessian(alpha,alpha) + Lagrange_multiplier
enddo




!------------Create random vector-------------

eta(1) =  rlfsr113()
eta(2) =  rlfsr113()
eta(3) =  rlfsr113()


!------------orthonormalization of the random vector---------

eta(:) = eta(:) - Mag(:)*(dot_product(eta,Mag))

tmp = normv(3,eta)

eta(:) = eta(:)/tmp

call cross_product(3,eta,Mag,Xi)



!---------Define Transformation matrix  (3*2 matrix)---------------
! U is a projection matrix

do alpha = 1 , 3

   do beta = 1 , 2

      if( beta .eq. 1) then

          V(alpha,beta) =  Xi(alpha)

      else

          V(alpha,beta) =  eta(alpha)

      endif

   enddo

enddo


!----------The projected Hessian Matrix  (2*2 matrix)--------------------

 do alpha = 1 , 3
    
    do beta = 1 , 2

       H_mult_V(alpha,beta) = dot_product( Hessian(alpha,:),V(:,beta) )

    enddo

 enddo

 do alpha = 1 , 2
    
    do beta = 1 , 2

       Phys_Hess_E(alpha,beta) = dot_product(  V(:,alpha),H_mult_V(:,beta) )

    enddo
    
       write(10,12) Phys_Hess_E(alpha,:)
       
 enddo

 close (10)
 
 
 
 
 
 call calc_hess_eig(1,Phys_Hess_E,eig,info)

 !print*,"eigenvalues are:",eig(:)
 

end subroutine




!**************************************************************************
!                          Dynamical factor 
!**************************************************************************

subroutine dynamical_factor(Phys_Hess_E,V,lambda)

 implicit none

 real(r8b),intent(in)::Phys_Hess_E(2,2),V(3,2)
 real(r8b),intent(out):: lambda
 real(r8b)::A(2,2),Hsp_dot_A(2,2),At_Hsp_A(2,2), S(2), At_Hsp_A_Dot_s(2)
 integer(i4b)::i,j


!------Pauli matrix-------
 
 A(1,1) = 0
 A(1,2) = 1
 A(2,1) = -1
 A(2,2) = 0
 
 
 do i = 1 , 2
    do j = 1 , 2
       Hsp_dot_A(i,j)= Dot_product(Phys_Hess_E(i,:) , A(:,j))
    enddo
 enddo
 
 
 do i = 1 , 2
    do j = 1 , 2
       At_Hsp_A(i,j)= Dot_product(A(:,i) , Hsp_dot_A(:,j) )
    enddo
 enddo
 
 

 
 !---eigenvector of Hsp corresponding to the negative eigenvalue is defined in the input file as unit_tangent---
 ! it is defined as 3D vector and we need to convert it to 2D vector then we use the projection matrix V
 
 do i = 1 , 2
    s(i) = Dot_product(V(:,i) , unit_tangent)
 enddo
 
 do i = 1 , 2
    At_Hsp_A_Dot_s(i) = Dot_product(At_Hsp_A(i,:) , s)
 enddo
 
 lambda = ( gama/(magnetic_moment*Mub) )*( dsqrt( Dot_product( s , At_Hsp_A_Dot_s ) ) )

end subroutine





!**************************************************************************
!                          lifetime calculations 
!**************************************************************************

subroutine lifetime(eig_stable, eig_saddle, lambda, tau)

 implicit none

 real(r8b),intent(in)::eig_stable(2), eig_saddle(2)
 real(r8b),intent(in)::lambda
 real(r8b),intent(out):: tau
 real(r8b):: pre_expofactor


 pre_expofactor = (2.0d0*Pi)/(lambda)*( dsqrt( eig_saddle(2)/(eig_stable(1)*eig_stable(2)) ) )
 
 tau = pre_expofactor*dexp(energy_barrier/thermal_energy)


 end subroutine
 
 
 
!**************************************************************************
!                            EigenValues 
!**************************************************************************

   !> Calculates eigenvalues of the projected Hessian
   subroutine calc_hess_eig(Natom,hess,eig,info)
  
      implicit none
      integer, intent(in) :: Natom
      
      real(r8b), dimension(2*Natom,2*Natom), intent(in) :: hess!< projected Hessian
      real(r8b), dimension(2*Natom), intent(out) :: eig !< eigenvalues
      integer, intent(out) :: info
      integer :: lda,lwork,i,j
      integer(8) :: lwmax
      real(r8b), allocatable :: work (:),hwork(:,:)
      
      lwmax = 5*Natom*Natom
      lda = 2*Natom
      allocate(work(lwmax))
      allocate(hwork(2*Natom,2*Natom))
      
      !print *,'allocatad:',allocated(work),allocated(hwork)
      
      
      do i=1,2*Natom
         do j=1,2*Natom
            hwork(i,j) = hess(i,j)
         end do
      end do
      
      
     lwork = -1
      
      
      call dsyev('N','U',2*Natom,hwork,lda,eig,work,lwork,info)
      
      lwork = min(lwmax,int(work(1)))
     ! print *,lwork,lwmax,info
      
      call dsyev('N','U',2*Natom,hwork,lda,eig,work,lwork,info)

      deallocate(work,hwork)     
   end subroutine calc_hess_eig


end program




!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! Random number generator RLFSR113  - real version
!
! Following a suggestion of Pierre L'Ecuyer 1997
! "Tables of maximally equidistributed combined LFSR generators"
! see http://www.iro.umontreal.ca/~lecuyer/papers.html
!
! A call to rlfsr113() gives one random real in the open
! interval (0,1).
!
! Before using rlfsr113 call lfsrinit(seed) to initialize
! the generator by random integers produced by Park/Millers
! minimal standard LCG.
! Seed should be any positive integer.
! 
! FORTRAN version by Thomas Vojta, vojta@physik.tu-chemnitz.de
! 
! History:
!  04 Feb 1998    v0.9    first FORTRAN implementation
!  05 Feb 1998    v0.91   corrected multiplicator am to 1/(2^31)
!  15 Apr 1999    v0.92   added real*8 rlfs113 in lfsrinit
!  10 Feb 2000    v0.93   changed lfsrinit to allow for CONSTANT arguments
! 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      FUNCTION rlfsr113 ()
      implicit none
      integer,parameter         :: i4b= kind(2147483647)
      integer,parameter         :: r8b= kind(1.D200)

      real(r8b),parameter       :: am=4.656612873077d-10

      real(r8b)    rlfsr113             
      integer(i4b) b,z1,z2,z3,z4
      common /lfsrcom/z1,z2,z3,z4

      b  = ishft(ieor(ishft(z1,6),z1),-13)
      z1 = ieor(ishft(iand(z1,-2),18),b)

      b  = ishft(ieor(ishft(z2,2),z2),-27)
      z2 = ieor(ishft(iand(z2,-8),2),b)

      b  = ishft(ieor(ishft(z3,13),z3),-21)
      z3 = ieor(ishft(iand(z3,-16),7),b)

      b  = ishft(ieor(ishft(z4,3),z4),-12)
      z4 = ieor(ishft(iand(z4,-128),13),b)

      rlfsr113=ishft( ieor(ieor(ieor(z1,z2),z3),z4) , -1)*am
!      print *, 'RLFSR113 ', rlfsr113

      return
      end


      SUBROUTINE lfsrinit(iinit)
      implicit none
      integer,parameter         :: r8b= kind(1.D200)
      integer,parameter         :: i4b= kind(2147483647)

      integer(i4b) idum,ia,im,iq,ir,iinit
      integer(i4b) k,z1,z2,z3,z4,c1,c2,c3,c4
      real(r8b)    rlfsr113,rdum
      parameter (ia=16807,im=2147483647,iq=127773,ir=2836)
      common /lfsrcom/z1,z2,z3,z4

      data c1 /B'11111111111111111111111111111110'/
      data c2 /B'11111111111111111111111111111000'/
      data c3 /B'11111111111111111111111111110000'/
      data c4 /B'11111111111111111111111110000000'/
      if ((c1.ne.-2).or.(c2.ne.-8).or.(c3.ne.-16).or.(c4.ne.-128)) then
         print *,'Nonstandard integer representation. Stoped.'
         stop
      endif

      idum=iinit
      if (idum.le.0) idum=1
      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.2) then
         z1=idum+2 
      else 
         z1=idum
      endif
      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.8) then 
         z2=idum+8 
      else 
         z2=idum
      endif
      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.16) then
         z3=idum+16 
      else 
         z3=idum
      endif
      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.128) then
         z4=idum+128 
      else 
         z4=idum
      endif

      rdum=rlfsr113()
      
      return
      end 
