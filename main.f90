

program main


use lib


implicit none
integer :: i,j,iiseed,count1,count2,tr_n,number_copies,Equilib_steps,Switch_steps,counter,info

real(q2), allocatable:: ene(:),traj1(:,:,:),traj2(:,:,:),zero_field(:,:),Bopt(:,:),Bopt1(:,:)
real(q2), allocatable:: Hess_Analy(:,:),Phys_Hess_E(:,:)
real(q2):: Te,dt,alpha,u,mc,s,c,d,tol,cutoff,magnetic_moment
real(q2)::Anis_Vec1(3),Anis_Vec2(3),Initail_state(3),time(2),eig(2) 
real(q2)::K_ani1,psi,Equilib_time,ST,alph_pertr,tmp1,tmp2,success
character(len=1000)::statuss

logical :: exi
	character(15)::fname_ini,fname_fin,mode_elli,mode_dw,fname_part
	character(50)::x
	character(25)::fname_work
	character(1) :: mode_ud !direction of dw motion
	character(3) :: num1,num2,tmp_ch,mode_ts
	character(9):: num

   iiseed = 10000*seed()

   
   !***************************************************************************************   
   !***************************************************************************************
   !VERY IMPORRTANT NOTES ARE WRITTEN IN THE END OF THIS CODE, PLEASE READ THEM CARFULLY
   !BEFORE SUBMITING ANY CALCULATIONS. 
   !***************************************************************************************
   !***************************************************************************************
   
   
      
   
   !**************************************TESTING**************************************** 
   !---here we test the subroutine that generates s = sn(u|m), c=cn(u|m), d=dn(u|m)******
   !u = 4.73d0
   !mc = 11d0
   !call gscdv(u,mc,s,c,d)
   !print *,'s:',s,'c:',c,'d:',d
   !*************************************************************************************
   

   !****READING DATA FROM INPUT_FILE****
   call read_input_file   

   !**********ALLOCATING MEMORY FOR VARIABLES***********
   !******NOTES REALATED TO THE ALLOCATED VARIABLES*****
   !"ene" is the energy of the system     
   !"traj1" is used for the equilibration process
   !"traj2" is used for the switching process
   !"zero_field" has zero entries because it is used in the LLG subroutine during the equilibration process, i.e., we do not want to apply pulse during the process
   !"Bopt" is the deterministic pulse that triggers the reversal motion
   !"Hess_Analy" produces a 3*3 matrix, which is a 2nd order derivative of the energy in Cartesian coordinates
   !"Phys_Hess_E" projects the Hess_Analy matrix on the tangent space, and produces a 2*2 matrix
   !***************************************************   
   allocate(ene(Equilib_steps+Switch_steps+Equilib_steps),traj1(3,Equilib_steps,number_copies),traj2(3,Switch_steps,number_copies)) 
   allocate(Bopt(4,Switch_steps-1),Bopt1(4,Switch_steps-1),zero_field(4,Equilib_steps-1),Hess_Analy(3,3),Phys_Hess_E(2,2))
   

   !***********INITIALIZATION*********
   zero_field=0d0
   

   print*,"----starting ASD calculations----"

   
   !******DEFINING THE INITIAL State******
   do i = 1 , number_copies
      traj1(1,1,i) = Initail_state(1)
      traj1(2,1,i) = Initail_state(2)
      traj1(3,1,i) = Initail_state(3)
      call normalize_vec(3,traj1(:,1,i))
   end do
   
   !print *, 'Initial energy:',cene(K_ani1,Anis_Vec1,psi,Anis_Vec2,traj1(:,1,1))


   !********************************TESTING*************************************
   !***************Caculating the eigenvalues of the initial state**************
   !---here "Hessian_Analy" produces a 3*3 matrix, which is used as input for the subroutine "Hessian_Phys_E". 
   !---here "Hessian_Phys_E" produces a 2*2 matrix, and from this matrix we compute its eigenvalues
   !---If the two eigenvalues are +ve, the initial state is a minimum
   !---If the two eigenvalues are -ve, the initial state is a maximum
   !---If one of the eigenvalues is +ve and the other one is -ve, the initial state is a saddle point
   
   !call Hessian_Analy(K_ani1,Anis_Vec1,psi,Anis_Vec2,Hess_Analy)
   !call Hessian_Phys_E(traj1(:,1,1),K_ani1,Anis_Vec1,psi,Anis_Vec2,Hess_Analy,Phys_Hess_E)
   !pause
   !****************************************************************************   
   


   call read_pulse(Switch_steps,magnetic_moment,K_ani1,Bopt) !---here we read the pre-calculated optimal pulse using our OCT code
   
   !******INITIAL EQUILIBRATION PROCESS******
   !Note: here LLG subroutine takes zero_field(1:3,:) array as input since we are just equilibrating the system
   !*****************************************
   call LLG(Equilib_steps,number_copies,K_ani1,Anis_Vec1,psi,Anis_Vec2,alph_pertr,Te,Equilib_time,zero_field(1:3,:),traj1) 
   

   !*********************************************************************************
   !**********Testing if spontaneous switching occurs during the equilibration*******
   !*********************************************************************************     
   !do i = 1 , Equilib_steps
   !   do j = 1 , number_copies
   !      IF (traj1(3,i,j)<-dabs(cutoff)) THEN
   !      call nrerror('WARNING! spontaneous switching occurs during the thermalization process!')
   !      EndIF
   !   end do
   !end do

!pause


   !*********************************************************************************
   !*************CALCULATING AVERAGE ENERGY AND WRITING DATA INTO FILE**************
   !*********************************************************************************     
   do i = 1 , Equilib_steps
      ene(i) = 0d0
      do j = 1 , number_copies
         ene(i) = ene(i) + cene(K_ani1,Anis_Vec1,psi,Anis_Vec2,traj1(:,i,j))
      end do
      ene(i) = ene(i)/real(number_copies,q2)
   end do
   
   dt = Equilib_time/real(Equilib_steps-1,q2)
   
   open(99, file='traj.dat', access = 'sequential',action = 'write', status = 'replace')
   do j = 1 , Equilib_steps
      write(99,'(es16.8E3,a,es16.8E3,a,es16.8E3,a,es16.8E3,a,es16.8E3,a,es16.8E3,a,es16.8E3)',advance = 'yes') dt*real(j-1,q2),'   ',traj1(1,j,tr_n),'   ',traj1(2,j,tr_n),'   ',traj1(3,j,tr_n),'   ',dotprod(3,traj1(:,j,tr_n),traj1(:,j,tr_n)),'   ',cene(K_ani1,Anis_Vec1,psi,Anis_Vec2,traj1(:,j,tr_n)),'   ',ene(j)
   end do
      
   
   
   !********************************************
   !After finishing the initial equilibration, the fial state we obtain, i.e., traj1(:,Equilib_steps,:)
   !becomes the initial state for the path, i.e., traj2(:,1,:) we use during the reversal process
   !It is very important to normalize the vector just to make sure that |vec{s}|=1, as described below
   !********************************************
   do i = 1 , number_copies
      traj2(1,1,i) = traj1(1,Equilib_steps,i)
      traj2(2,1,i) = traj1(2,Equilib_steps,i)
      traj2(3,1,i) = traj1(3,Equilib_steps,i)
      call normalize_vec(3,traj2(:,1,i))
   end do
   
 
 
   !****************SWITCHING PROCESS***********
   !Note: here LLG subroutine takes Bopt(1:3,:) array as input which will be used for inducing the magnetization reversal process
   !********************************************
   call LLG(Switch_steps,number_copies,K_ani1,Anis_Vec1,psi,Anis_Vec2,alph_pertr,Te,ST,Bopt(1:3,:),traj2) 
      
      
      
   !********************************************    
   !After applying the pulse, the fial state we obtain, i.e., traj2(:,Switch_steps,:)
   !becomes the initial state for the path, i.e., traj1(:,1,:) we use for the thermalization process
   !********************************************   
   do i = 1 , number_copies
      traj1(1,1,i) = traj2(1,Switch_steps,i)
      traj1(2,1,i) = traj2(2,Switch_steps,i)
      traj1(3,1,i) = traj2(3,Switch_steps,i)
      call normalize_vec(3,traj1(:,1,i))
   end do

   !******SECOND EQUILIBRATION PROCESS****** 
   !this step is important as we need to equilibrate the system again after the pulse is terminated    
   !****************************************
   call LLG(Equilib_steps,number_copies,K_ani1,Anis_Vec1,psi,Anis_Vec2,alph_pertr,Te,Equilib_time,zero_field(1:3,:),traj1) 
   


   !*********************************************************************************
   !*************CALCULATING AVERAGE ENERGY AND WRITING DATA INTO FILES**************
   !*********************************************************************************     
   do i = Equilib_steps+1,Equilib_steps+Switch_steps
      ene(i) = 0d0
      do j = 1 , number_copies
         ene(i) = ene(i) + cene(K_ani1,Anis_Vec1,psi,Anis_Vec2,traj2(:,i-Equilib_steps,j))
      end do
      ene(i) = ene(i)/real(number_copies,q2)
   end do
   
   do i = Equilib_steps+Switch_steps+1,Equilib_steps+Switch_steps+Equilib_steps
      ene(i) = 0d0
      do j = 1 , number_copies
         ene(i) = ene(i) + cene(K_ani1,Anis_Vec1,psi,Anis_Vec2,traj1(:,i-Equilib_steps-Switch_steps,j))
      end do
      ene(i) = ene(i)/real(number_copies,q2)
   end do
   
   
      
   dt = ST/real(Switch_steps-1,q2)
      
   do j = 2 , Switch_steps
      write(99,'(es16.8E3,a,es16.8E3,a,es16.8E3,a,es16.8E3,a,es16.8E3,a,es16.8E3,a,es16.8E3)',advance = 'yes') Equilib_time+dt*real(j-1,q2),'   ',traj2(1,j,tr_n),'   ',traj2(2,j,tr_n),'   ',traj2(3,j,tr_n),'   ',dotprod(3,traj2(:,j,tr_n),traj2(:,j,tr_n)),'   ',cene(K_ani1,Anis_Vec1,psi,Anis_Vec2,traj2(:,j,tr_n)),'   ',ene(j+Equilib_steps)
   end do
   
   dt = Equilib_time/real(Equilib_steps-1,q2)
      
   do j = 2 , Equilib_steps
      write(99,'(es16.8E3,a,es16.8E3,a,es16.8E3,a,es16.8E3,a,es16.8E3,a,es16.8E3,a,es16.8E3)',advance = 'yes') Equilib_time+ST+dt*real(j-1,q2),'   ',traj1(1,j,tr_n),'   ',traj1(2,j,tr_n),'   ',traj1(3,j,tr_n),'   ',dotprod(3,traj1(:,j,tr_n),traj1(:,j,tr_n)),'   ',cene(K_ani1,Anis_Vec1,psi,Anis_Vec2,traj1(:,j,tr_n)),'   ',ene(j+Equilib_steps+Switch_steps)
   end do
      
   close(99)

   dt = ST/real(Switch_steps-1,q2)

   open(99, file='bopt.dat', access = 'sequential',action = 'write', status = 'replace')
   do j = 1 , Switch_steps-1
      write(99,'(es16.8E3,a,es16.8E3,a,es16.8E3,a,es16.8E3,a,es16.8E3)',advance = 'yes') dt*(real(j,q2)-0.5d0),'   ',Bopt(1,j),'   ',Bopt(2,j),'   ',Bopt(3,j),'   ',Bopt(4,j)
   end do
   close(99)

     
      
   count1 = 0
   count2 = 0
      
   !do i = 1 , number_copies
   !   if (traj1(3,Equilib_steps,i)<-dabs(cutoff)) then
   !      count1 = count1+1
   !   elseif (traj1(3,Equilib_steps,i)>cutoff) then
   !      count2 = count2+1
   !      print *,'unsuccessful copy:',i
   !   else
   !      print *,'undefined copy:',i
   !      print *,'mz:', traj1(3,Equilib_steps,i)
   !   end if
   !end do

   do i = 1 , number_copies
      if (traj1(1,Equilib_steps,i)<-0.5d0 .and. traj1(2,Equilib_steps,i)<0.5d0) then
         count1 = count1+1
      else
         count2 = count2+1
         print *,'unsuccessful copy:',i
      end if
   end do

      
   print *,'success:',count1
   print *,'fail:',count2
   print *,'undefined:',number_copies-count1-count2
   print *,'success rate,%:',real(count1,q2)/real(number_copies,q2)*100d0
   print*, 'error:',1.96d0*dsqrt( real(count1,q2)/real(number_copies,q2)*(1.0d0 - real(count1,q2)/real(number_copies,q2)) / real(number_copies,q2) )   
   !In the error calculation, we have multiply it by 1.96 which is for 95% confidence level
   !see: https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval    
 



   !*****************************************************************
   !*************Writing data into the file "results.out"************
   !*****************************************************************
   
   
   open(10,file="results.out",action="write",status="replace")
   11 format(1000(A,x))   
   write(10,11)'===================================================================='

   count1 = 0
   count2 = 0
      
!   do i = 1 , number_copies
!      if (traj1(1,Equilib_steps,i)<-0.5d0 .and. traj1(2,Equilib_steps,i)<0.5d0) then
!         count1 = count1+1
!      elseif (traj1(1,Equilib_steps,i)>-0.5d0 .and. traj1(2,Equilib_steps,i)>0.5d0) then
!         count2 = count2+1
!         write(10,'(a,I6)') 'unsuccessful copy:',i
!      else
!         write(10,'(a,I6,a,a,es16.8E3)') 'undefined copy:',i,'   ','mz:', traj1(3,Equilib_steps,i)
!      end if
!   end do

   do i = 1 , number_copies
      if (traj1(1,Equilib_steps,i)<-0.5d0 .and. traj1(2,Equilib_steps,i)<0.5d0) then
         count1 = count1+1
      else 
         count2 = count2+1
         write(10,'(a,I6)') 'unsuccessful copy:',i
      end if
   end do

   write(10,*)
   write(10,*) '*********Summary********'
   write(10,'(a,I6)') 'number of copies:',number_copies
   write(10,'(a,I6)') 'success:',count1
   write(10,'(a,I6)') 'fail:',count2
   write(10,'(a,I6)') 'undefined:',number_copies-count1-count2
   write(10,'(a,es16.8E3)') 'success rate,%:',real(count1,q2)/real(number_copies,q2)*100d0
   write(10,'(a,es16.8E3)') 'error:',1.96d0*dsqrt( real(count1,q2)/real(number_copies,q2)*(1.0d0 - real(count1,q2)/real(number_copies,q2)) / real(number_copies,q2) )  



 
   

   
   
   !**********DEALLOCATING ARRAYS************
   deallocate(ene,traj1,traj2,zero_field,Bopt,Hess_Analy,Phys_Hess_E)

contains

!**********************************************************
!               Reading data from Input File
!**********************************************************
subroutine read_input_file
  implicit none
  character(len=100) :: buffer, label,check,Images,Initial      
  integer:: pos  
  integer, parameter:: fh = 15 
  integer:: ios = 0
  integer:: line = 0
  logical:: dir_e  



  !*********checking if file exists!************
  inquire( file = "Input_File.in" , exist=dir_e )
         if ( dir_e ) then
            goto 7
         else
            call nrerror('File is missing: please check Input_File.in!!')
         end if
  !*********************************************

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
 
          case ('Initail_state')
            read(buffer, *, iostat=ios) Initail_state         
          case ('Anis_vec1')
            read(buffer, *, iostat=ios) Anis_vec1
          case ('Anis_vec2')
            read(buffer, *, iostat=ios) Anis_vec2
          case ('K_ani1')
            read(buffer, *, iostat=ios) K_ani1
          case ('psi')
            read(buffer, *, iostat=ios) psi
          
          case ('number_copies')
            read(buffer, *, iostat=ios) number_copies
          case ('Equilib_time')
            read(buffer, *, iostat=ios) Equilib_time 
          case ('Equilib_steps')
            read(buffer, *, iostat=ios) Equilib_steps  
          case ('Te')
            read(buffer, *, iostat=ios) Te  
          case ('ST')
            read(buffer, *, iostat=ios) ST 
          case ('Switch_steps')
            read(buffer, *, iostat=ios) Switch_steps
          case ('magnetic_moment')
            read(buffer, *, iostat=ios) magnetic_moment                        
          case ('alpha')
            read(buffer, *, iostat=ios) alpha  
          case ('alph_pertr')
            read(buffer, *, iostat=ios) alph_pertr              
          case ('tol')
            read(buffer, *, iostat=ios) tol                 
          case ('cutoff')
            read(buffer, *, iostat=ios) cutoff  
          case ('tr_n')
            read(buffer, *, iostat=ios) tr_n                                  
                   
                                                          
        end select
        
     end if

  end do

  call normalize_vec(3,Anis_vec1)
  call normalize_vec(3,Anis_vec2)
  

end subroutine read_input_file




!**********************************************************
!       Reading the components of the optimal pulse
!**********************************************************
subroutine read_pulse(Switch_steps,magnetic_moment,K_ani1,pulse)
      implicit none
      integer, intent(in) :: Switch_steps
      real(q2), intent(in) :: magnetic_moment,K_ani1
      real(q2), intent(out) :: pulse(4,Switch_steps-1)
      real(q2), parameter::mu_b = 0.057883818060d0
      real(q2) :: read_file(6)
      real(q2) :: mu
      integer :: i
      character(len=1000)::tmp_line
      logical:: dir_e  



      !*********checking if file exists!************
      inquire( file = "pulse.in" , exist=dir_e )
             if ( dir_e ) then
                goto 12
             else
                call nrerror('File is missing: please check pulse.in!!')
             end if
      !*********************************************

      12 open(1,file="pulse.in",action="read")

         read(1,*) tmp_line
         read(1,*) tmp_line
         
         mu = magnetic_moment*mu_b
         do i = 1 , Switch_steps-1
                read(1,*) read_file
!                pulse(1,i) = read_file(4)*mu/(2.0d0*K_ani1)   !X-component of the pulse in units of [K/mu]
!                pulse(2,i) = read_file(5)*mu/(2.0d0*K_ani1)   !Y-component of the pulse in units of [K/mu]
!                pulse(3,i) = read_file(6)*mu/(2.0d0*K_ani1)    !Z-component of the pulse in units of [K/mu]
!                pulse(4,i) = dsqrt( pulse(1,i)*pulse(1,i) + pulse(2,i)*pulse(2,i) + pulse(3,i)*pulse(3,i) )  !---magnitude of the pulse

                pulse(1,i) = read_file(4)*mu/(2.0d0*K_ani1)   !X-component of the pulse in units of [K/mu]
                pulse(2,i) = read_file(5)*mu/(2.0d0*K_ani1)   !Y-component of the pulse in units of [K/mu]
                pulse(3,i) = read_file(6)*mu/(2.0d0*K_ani1)    !Z-component of the pulse in units of [K/mu]
                pulse(4,i) = dsqrt( pulse(1,i)*pulse(1,i) + pulse(2,i)*pulse(2,i) + pulse(3,i)*pulse(3,i) )  !---magnitude of the pulse
         enddo
         close(1)      

end subroutine read_pulse

   
end program main

