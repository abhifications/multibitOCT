module lib
     
   USE MKL_VSL_TYPE
   USE MKL_VSL
   
  
   INTEGER,PARAMETER :: q2=SELECTED_REAL_KIND(15,305)
   
  
   real(q2), parameter :: pi = 3.14159265358979323d0, pid = .318309886183790671d0, &
                          othree = 1d0/3d0
  
   real(q2),external :: dznrm2,dnrm2
   
   contains
   
   subroutine gscdv(u,m,s,c,d)
      implicit none
      real(q2), intent(in) :: u,m
      real(q2), intent(out) :: s,c,d
      real(q2) :: saux,caux,daux
      real(q2) :: uaux
      call gscd(u*dsqrt(1d0+m),1d0/(m+1d0),saux,caux,daux)
      
      c = caux/daux
      s = saux/(daux*dsqrt(1d0+m))
      d = dsqrt(1d0 + m*s*s)

   
   end subroutine gscdv
   
   subroutine gscd(u,mc,s,c,d)
!
!	Double precision subroutine to compute three Jacobian elliptic functions simultaneously
!
!   For general argument: -infty < u < infty
!
!     Reference: T. Fukushima, (2012) Numer. Math. DOI 10.1007/s00211-012-0498-0
!       "Precise and Fast Computation of Jacobian Elliptic Functions by
!        Conditional Duplication"
!
!     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!
!     Used subprograms: scd2, elk
!
!     Inputs: u = argument, mc = 1-m, 0 < mc <= 1
!
!     Output: s = sn(u|m), c=cn(u|m), d=dn(u|m)
!
real(q2) u,mc,s,c,d
real(q2) m,kc,ux,k,kh,kh3,kh5,kh7,k2,k3,k4,sx,cx,dx
!real(q2) elk
!
m=1.d0-mc
kc=sqrt(mc)
ux=abs(u)
if(ux.lt.0.785d0) then
    call scd2(ux,mc,s,c,d)
else
    k=elk(mc)
    kh=k*0.5d0; kh3=k*1.5d0; kh5=k*2.5d0; kh7=k*3.5d0;
    k2=k*2.d0; k3=k*3.d0; k4=k*4.d0
    ux=ux-k4*dble(int(ux/k4))
    if(ux.lt.kh) then
        call scd2(ux,mc,s,c,d)
    elseif(ux.lt.k) then
        ux=k-ux
        call scd2(ux,mc,s,c,d)
        sx=c/d; c=kc*s/d; s=sx; d=kc/d
    elseif(ux.lt.kh3) then
        ux=ux-k
        call scd2(ux,mc,s,c,d)
        sx=c/d; c=-kc*s/d; s=sx; d=kc/d
    elseif(ux.lt.k2) then
        ux=k2-ux
        call scd2(ux,mc,s,c,d)
        c=-c
    elseif(ux.lt.kh5) then
        ux=ux-k2
        call scd2(ux,mc,s,c,d)
        s=-s; c=-c
    elseif(ux.lt.k3) then
        ux=k3-ux
        call scd2(ux,mc,s,c,d)
        sx=-c/d; c=-kc*s/d; s=sx; d=kc/d
    elseif(ux.lt.kh7) then
        ux=ux-k3
        call scd2(ux,mc,s,c,d)
        sx=-c/d; c=kc*s/d; s=sx; d=kc/d
    else
        ux=k4-ux
        call scd2(ux,mc,s,c,d)
        s=-s
    endif
endif
if(u.lt.0.d0) s=-s
return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine scd2(u,mc,s,c,d)
!
!	Double precision subroutine to compute three Jacobian elliptic functions simultaneously
!
!   For limited argument: 0 <= u < K/2
!
!     Reference: T. Fukushima, (2012) Numer. Math. DOI 10.1007/s00211-012-0498-0
!       "Precise and Fast Computation of Jacobian Elliptic Functions by
!        Conditional Duplication"
!
!     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!
!     Inputs: u = argument, mc = 1-m, 0 < mc <= 1
!
!     Output: s = sn(u|m), c=cn(u|m), d=dn(u|m)
!
    real(q2) u,mc,s,c,d
    real(q2) B10,B11,B20,B21,B22,m,uA,uT,u0,v,a,b,y,z,my,mc2,m2,xz,w
    integer n,j,i
    parameter (B10=1.d0/24.d0,B11=1.d0/6.d0,B20=1.d0/720.d0)
    parameter (B21=11.d0/180.d0,B22=1.d0/45.d0)
    m=1.d0-mc; uA=1.76269d0+mc*1.16357d0; uT=5.217d-3-m*2.143d-3; u0=u
    do n=0,20
        if(u0.lt.uT) goto 1
        u0=u0*0.5d0
    enddo
    write(*,*) "(scd2) Too large input argument: u=", u
1 continue
    v=u0*u0; a=1.d0; b=v*(0.5d0-v*(B10+m*B11-v*(B20+m*(B21+m*B22))))
    if(u.lt.uA) then
        do j=1,n
            y=b*(a*2.d0-b); z=a*a; my=m*y; b=(y*2.d0)*(z-my); a=z*z-my*y
        enddo
    else
        do j=1,n
            y=b*(a*2.d0-b); z=a*a; my=m*y
            if(z.lt.my*2.d0) goto 2
            b=(y*2.d0)*(z-my); a=z*z-my*y
        enddo
    endif
    b=b/a; y=b*(2.d0-b); c=1.d0-b; s=sqrt(y); d=sqrt(1.d0-m*y)
     return
2 continue
    c=a-b; mc2=mc*2.d0; m2=m*2.d0
    do i=j,n
        x=c*c; z=a*a; w=m*x*x-mc*z*z; xz=x*z; c=mc2*xz+w; a=m2*xz-w
    enddo
    c=c/a; x=c*c; s=sqrt(1.d0-x); d=sqrt(mc+m*x)
    return; end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(q2) function elk(mc)
!c
!c	Double precision complete elliptic integral of the first kind
!c
!c     Reference: T. Fukushima, (2009) Celest. Mech. Dyn. Astron. 105, 305-328
!c        "Fast Computation of Complete Elliptic Integrlals and Jacobian
!c         Elliptic Functions"
!c
!c     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!c
!c     Inputs: mc   = complementary parameter 0 <= mc   <= 1
!c
!c     Output: elk
!c
	real(q2) mc
	real(q2) mcold,PIHALF,PIINV,elkold,TINY,m,mx,P,Q
	real(q2) kkc,nome
!c
	real(q2) D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,D12,D13,D14
	parameter (D1=1.d0/16.d0,D2=1.d0/32.d0,D3=21.d0/1024.d0)
	parameter (D4=31.d0/2048.d0,D5=6257.d0/524288.d0)
	parameter (D6=10293.d0/1048576.d0,D7=279025.d0/33554432.d0)
	parameter (D8=483127.d0/67108864.d0)
	parameter (D9=435506703.d0/68719476736.d0)
	parameter (D10=776957575.d0/137438953472.d0)
	parameter (D11=22417045555.d0/4398046511104.d0)
	parameter (D12=40784671953.d0/8796093022208.d0)
	parameter (D13=9569130097211.d0/2251799813685248.d0)
	parameter (D14=17652604545791.d0/4503599627370496.d0)
!c
	logical first/.TRUE./
	
	save first,mcold,PIHALF,PIINV,elkold,TINY
!c
   
	if(first) then
		first=.FALSE.
		mcold=1.d0
		PIHALF=atan(1.d0)*2.d0
		PIINV=0.5d0/PIHALF
		elkold=PIHALF
		TINY=1.d-99
	endif
	m=1.d0-mc
	if(abs(m).lt.1.d-16) then
		elk=PIHALF
	elseif(abs(mc-mcold).lt.1.11d-16*mc) then
		elk=elkold
	elseif(mc.lt.TINY) then
	  elk=1.3862943611198906d0-0.5d0*log(TINY)
	elseif(mc.lt.1.11d-16) then
	  elk=1.3862943611198906d0-0.5d0*log(mc)
	elseif(mc.lt.0.1d0) then
		nome=mc*(D1+mc*(D2+mc*(D3+mc*(D4+mc*(D5+mc*(D6 &
     		+mc*(D7+mc*(D8+mc*(D9+mc*(D10+mc*(D11+mc*(D12 &
     		+mc*(D13+mc*D14)))))))))))))
		mx=mc-0.05d0
!c
!c	K'
!c
		kkc=1.591003453790792180d0+mx*( &
     		0.416000743991786912d0+mx*( &
     		0.245791514264103415d0+mx*( &
     		0.179481482914906162d0+mx*( &
     		0.144556057087555150d0+mx*( &
     		0.123200993312427711d0+mx*( &
     		0.108938811574293531d0+mx*( &
     		0.098853409871592910d0+mx*( &
     		0.091439629201749751d0+mx*( &
     		0.085842591595413900d0+mx*( &
     		0.081541118718303215d0))))))))))
!c
		elk=-kkc*PIINV*log(nome)
	elseif(m.le.0.1d0) then
		mx=m-0.05d0
		elk=1.591003453790792180d0+mx*( &
     		0.416000743991786912d0+mx*( &
     		0.245791514264103415d0+mx*( &
     		0.179481482914906162d0+mx*( &
     		0.144556057087555150d0+mx*( &
     		0.123200993312427711d0+mx*( &
     		0.108938811574293531d0+mx*( &
     		0.098853409871592910d0+mx*( &
     		0.091439629201749751d0+mx*( &
     		0.085842591595413900d0+mx*( &
     		0.081541118718303215d0))))))))))
	elseif(m.le.0.2d0) then
		mx=m-0.15d0
		elk=1.635256732264579992d0+mx*( &
     		0.471190626148732291d0+mx*( &
     		0.309728410831499587d0+mx*( &
     		0.252208311773135699d0+mx*( &
     		0.226725623219684650d0+mx*( &
     		0.215774446729585976d0+mx*( &
     		0.213108771877348910d0+mx*( &
     		0.216029124605188282d0+mx*( &
     		0.223255831633057896d0+mx*( &
     		0.234180501294209925d0+mx*( &
     		0.248557682972264071d0+mx*( &
     		0.266363809892617521d0+mx*( &
     		0.287728452156114668d0))))))))))))
	elseif(m.le.0.3d0) then
		mx=m-0.25d0
		elk=1.685750354812596043d0+mx*( &
     		0.541731848613280329d0+mx*( &
     		0.401524438390690257d0+mx*( &
     		0.369642473420889090d0+mx*( &
     		0.376060715354583645d0+mx*( &
     		0.405235887085125919d0+mx*( &
     		0.453294381753999079d0+mx*( &
     		0.520518947651184205d0+mx*( &
     		0.609426039204995055d0+mx*( &
     		0.724263522282908870d0+mx*( &
     		0.871013847709812357d0+mx*( &
     		1.057652872753547036d0)))))))))))
	elseif(m.le.0.4d0) then
		mx=m-0.35d0
		elk=1.744350597225613243d0+mx*( &
     		0.634864275371935304d0+mx*( &
     		0.539842564164445538d0+mx*( &
     		0.571892705193787391d0+mx*( &
     		0.670295136265406100d0+mx*( &
     		0.832586590010977199d0+mx*( &
     		1.073857448247933265d0+mx*( &
     		1.422091460675497751d0+mx*( &
     		1.920387183402304829d0+mx*( &
     		2.632552548331654201d0+mx*( &
     		3.652109747319039160d0+mx*( &
     		5.115867135558865806d0+mx*( &
     		7.224080007363877411d0))))))))))))
	elseif(m.le.0.5d0) then
		mx=m-0.45d0
		elk=1.813883936816982644d0+mx*( &
     		0.763163245700557246d0+mx*( &
     		0.761928605321595831d0+mx*( &
     		0.951074653668427927d0+mx*( &
     		1.315180671703161215d0+mx*( &
     		1.928560693477410941d0+mx*( &
     		2.937509342531378755d0+mx*( &
     		4.594894405442878062d0+mx*( &
     		7.330071221881720772d0+mx*( &
     		11.87151259742530180d0+mx*( &
     		19.45851374822937738d0+mx*( &
     		32.20638657246426863d0+mx*( &
     		53.73749198700554656d0+mx*( &
     		90.27388602940998849d0)))))))))))))
	elseif(m.le.0.6d0) then
		mx=m-0.55d0
		elk=1.898924910271553526d0+mx*( &
     		0.950521794618244435d0+mx*( &
     		1.151077589959015808d0+mx*( &
     		1.750239106986300540d0+mx*( &
     		2.952676812636875180d0+mx*( &
     		5.285800396121450889d0+mx*( &
     		9.832485716659979747d0+mx*( &
     		18.78714868327559562d0+mx*( &
     		36.61468615273698145d0+mx*( &
     		72.45292395127771801d0+mx*( &
     		145.1079577347069102d0+mx*( &
     		293.4786396308497026d0+mx*( &
     		598.3851815055010179d0+mx*( &
     		1228.420013075863451d0+mx*( &
     		2536.529755382764488d0))))))))))))))
	elseif(m.le.0.7d0) then
		mx=m-0.65d0
		elk=2.007598398424376302d0+mx*( &
     		1.248457231212347337d0+mx*( &
     		1.926234657076479729d0+mx*( &
     		3.751289640087587680d0+mx*( &
     		8.119944554932045802d0+mx*( &
     		18.66572130873555361d0+mx*( &
     		44.60392484291437063d0+mx*( &
     		109.5092054309498377d0+mx*( &
     		274.2779548232413480d0+mx*( &
     		697.5598008606326163d0+mx*( &
     		1795.716014500247129d0+mx*( &
     		4668.381716790389910d0+mx*( &
     		12235.76246813664335d0+mx*( &
     		32290.17809718320818d0+mx*( &
     		85713.07608195964685d0+mx*( &
     		228672.1890493117096d0+mx*( &
     		612757.2711915852774d0))))))))))))))))
	elseif(m.le.0.8d0) then
		mx=m-0.75d0
		elk=2.156515647499643235d0+mx*( &
     		1.791805641849463243d0+mx*( &
     		3.826751287465713147d0+mx*( &
     		10.38672468363797208d0+mx*( &
     		31.40331405468070290d0+mx*( &
     		100.9237039498695416d0+mx*( &
     		337.3268282632272897d0+mx*( &
     		1158.707930567827917d0+mx*( &
     		4060.990742193632092d0+mx*( &
     		14454.00184034344795d0+mx*( &
     		52076.66107599404803d0+mx*( &
     		189493.6591462156887d0+mx*( &
     		695184.5762413896145d0+mx*( &
     		2.567994048255284686d6+mx*( &
     		9.541921966748386322d6+mx*( &
     		3.563492744218076174d7+mx*( &
     		1.336692984612040871d8+mx*( &
     		5.033521866866284541d8+mx*( &
     		1.901975729538660119d9+mx*( &
     		7.208915015330103756d9)))))))))))))))))))
	elseif(m.le.0.85d0) then
		mx=m-0.825d0
		elk=2.318122621712510589d0+mx*( &
     		2.616920150291232841d0+mx*( &
     		7.897935075731355823d0+mx*( &
     		30.50239715446672327d0+mx*( &
     		131.4869365523528456d0+mx*( &
     		602.9847637356491617d0+mx*( &
     		2877.024617809972641d0+mx*( &
     		14110.51991915180325d0+mx*( &
     		70621.44088156540229d0+mx*( &
     		358977.2665825309926d0+mx*( &
     		1.847238263723971684d6+mx*( &
     		9.600515416049214109d6+mx*( &
     		5.030767708502366879d7+mx*( &
     		2.654441886527127967d8+mx*( &
     		1.408862325028702687d9+mx*( &
     		7.515687935373774627d9)))))))))))))))
      else
		mx=m-0.875d0
		elk=2.473596173751343912d0+mx*( &
     		3.727624244118099310d0+mx*( &
     		15.60739303554930496d0+mx*( &
     		84.12850842805887747d0+mx*( &
     		506.9818197040613935d0+mx*( &
     		3252.277058145123644d0+mx*( &
     		21713.24241957434256d0+mx*( &
     		149037.0451890932766d0+mx*( &
     		1.043999331089990839d6+mx*( &
     		7.427974817042038995d6+mx*( &
     		5.350383967558661151d7+mx*( &
     		3.892498869948708474d8+mx*( &
     		2.855288351100810619d9+mx*( &
     		2.109007703876684053d10+mx*( &
     		1.566998339477902014d11+mx*( &
     		1.170222242422439893d12+mx*( &
     		8.777948323668937971d12+mx*( &
     		6.610124275248495041d13+mx*( &
     		4.994880537133887989d14+mx*( &
     		3.785974339724029920d15)))))))))))))))))))
	endif
!c
	mcold=mc
	elkold=elk
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
   
   
   SUBROUTINE nrerror(string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	write (*,*) 'nrerror: ',string
	STOP 'program terminated by nrerror'
	END SUBROUTINE nrerror
   
   FUNCTION gaussian(idum,rand,mean,std)
    IMPLICIT NONE

    INTEGER,PARAMETER :: q1=SELECTED_REAL_KIND(6,30)
    INTEGER,PARAMETER :: q2=SELECTED_REAL_KIND(15,305)

    INTEGER,INTENT(IN) :: idum
    REAL(q2),INTENT(IN) :: mean,std
    INTEGER,SAVE :: iset=0
    REAL(q2),SAVE :: gset
    REAL(q2) :: fac,rsq,gaussian,v1,v2
    REAL(q2),EXTERNAL :: rand

    IF ( iset == 0 ) THEN
      DO
        v1=2._q2*rand(idum)-1._q2
        v2=2._q2*rand(idum)-1._q2
        rsq=v1*v1+v2*v2
        IF (rsq < 1.0_q2 .and. rsq /= 0.0_q2) EXIT
      END DO
      fac=SQRT(-2.0_q2*LOG(rsq)/rsq)
      gset=v1*fac
      iset=1
      gaussian=v2*fac
    ELSE
      iset=0
      gaussian=gset
    END IF
    gaussian=mean+std*gaussian           ! std*std=var
  RETURN
  END FUNCTION gaussian



  FUNCTION seed()
    IMPLICIT NONE

    INTEGER,DIMENSION(8) :: iseed
    INTEGER :: seed

    CALL DATE_AND_TIME(values=iseed)
    seed=MOD(1000*iseed(7)+iseed(8),1024)

  RETURN
  END FUNCTION seed

  FUNCTION ran2(idum)
    IMPLICIT NONE

    INTEGER,PARAMETER :: q1=SELECTED_REAL_KIND(6,30)
    INTEGER,PARAMETER :: q2=SELECTED_REAL_KIND(15,305)

    INTEGER,INTENT(INOUT) :: idum
    INTEGER,PARAMETER :: im1=2147483563,im2=2147483399,imm1=im1-1,ia1=40014,        &
                         ia2=40692,iq1=53668,iq2=52774,ir1=12211,ir2=3791,ntab=32,  &
                         ndiv=1+imm1/ntab
    REAL(q1),PARAMETER :: am=1._q1/im1,eps=1.2e-7_q1,rnmx=1._q1-eps
    INTEGER :: j,k
    INTEGER,SAVE :: iy=0,idum2=123456789
    INTEGER,DIMENSION(ntab),SAVE :: iv=0
    REAL(q2) :: ran2

    IF (idum <= 0) THEN
      idum=MAX(-idum,1)
      idum2=idum
      DO j=ntab+8,1,-1
        k=idum/iq1
        idum=ia1*(idum-k*iq1)-k*ir1
        IF (idum < 0) idum=idum+im1
        IF (j <= ntab) iv(j)=idum
      END DO
      iy=iv(1)
    END IF
    k=idum/iq1
    idum=ia1*(idum-k*iq1)-k*ir1
    IF (idum < 0) idum=idum+im1
    k=idum2/iq1
    idum2=ia2*(idum2-k*iq2)-k*ir2
    IF (idum2 < 0) idum2=idum2+im2
    j=1+iy/ndiv
    iy=iv(j)-idum2
    iv(j)=idum
    IF (iy < 1) iy=iy+imm1
    ran2=MIN(am*iy,rnmx)
  RETURN
  END FUNCTION ran2

  real(q2) function elp(T,alpha,p)
   implicit none
   real(q2), intent(in) :: T,alpha,p
   
   
   !elp = 4d0*p*elk(-alpha*alpha*p*p)-T
   elp = 4d0*(1d0+alpha*alpha)*p*elkv(-alpha*alpha*p*p)-T
   !elp = 4d0*p*elkv(-alpha*alpha*p*p)-T
   
  end function elp
  
  real(q2) function elkv(x)
  implicit none
  real(q2), intent(in) :: x 
  
   elkv = elk(1d0/(-x+1d0))/dsqrt(-x+1d0)
    
  end function elkv
   
!   function wich finds the root of function 'func' between x1 and x2 with accuracy 'tol'
   FUNCTION bin_root(T,alpha,x1,x2,tol)
	
	IMPLICIT NONE

	REAL(q2), INTENT(IN) :: T,alpha,x1,x2,tol
	REAL(q2) :: bin_root
	
	INTEGER, PARAMETER :: ITMAX=100
	REAL(q2), PARAMETER :: EPS2=epsilon(x1)
	INTEGER :: iter
	REAL(q2) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
	a=x1
	b=x2
	fa=elp(T,alpha,a)
	fb=elp(T,alpha,b)

	if ((fa > 0.0d0 .and. fb > 0.0d0) .or. (fa < 0.0d0 .and. fb < 0.0d0)) then

		call nrerror('root must be bracketed for bin_root')
	end if	
	c=b
	fc=fb
	do iter=1,ITMAX
		if ((fb > 0.0d0 .and. fc > 0.0d0) .or. (fb < 0.0d0 .and. fc < 0.0d0)) then
			c=a
			fc=fa
			d=b-a
			e=d
		end if
		if (abs(fc) < abs(fb)) then
			a=b
			b=c
			c=a
			fa=fb
			fb=fc
			fc=fa
		end if
		tol1=2.0d0*EPS2*abs(b)+0.5d0*tol
		xm=0.5d0*(c-b)
		if (abs(xm) <= tol1 .or. fb == 0.0) then
			bin_root=b
			RETURN
		end if
		if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
			s=fb/fa
			if (a == c) then
				p=2.0d0*xm*s
				q=1.0d0-s
			else
				q=fa/fc
				r=fb/fc
				p=s*(2.0d0*xm*q*(q-r)-(b-a)*(r-1.0d0))
				q=(q-1.0d0)*(r-1.0d0)*(s-1.0d0)
			end if
			if (p > 0.0) q=-q
			p=abs(p)
			if (2.0d0*p  <  min(3.0d0*xm*q-abs(tol1*q),abs(e*q))) then
				e=d
				d=p/q
			else
				d=xm
				e=d
			end if
		else
			d=xm
			e=d
		end if
		a=b
		fa=fb
		b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
		fb=elp(T,alpha,b)
	end do
	call nrerror('bin_root: exceeded maximum iterations')
	bin_root=b
	END FUNCTION bin_root

   
   real(q2) function find_p(T,alpha,tol)
      implicit none
      real(q2), intent(in) :: T, alpha,tol
      real(q2) :: p2
      p2 = 1d0
      do while (elp(T,alpha,p2).le.0d0)
         p2 = p2*2d0
      end do
      
      find_p = bin_root(T,alpha,0d0,p2,tol)
   
   end function find_p
   
   
   !********************************************
   !    Optimal pulse of the external field
   !********************************************
   subroutine calc_bopt(N,T,alpha,tol,bopt)
      implicit none
      integer, intent(in) :: N
      real(q2), intent(in) :: T,alpha,tol
      real(q2), intent(out) :: bopt(4,N-1)
      real(q2) :: dt,b,tt,p,s,c,d,phi, sig, dphi
      integer(q2) :: i
      
      !---------------------------------------------------------
      !The optimal pulse is reconstructed from OCP,
      !where we have the solution of OCP in spherical
      !coordinates with the following form:
      !theta(t) = 0.5*am[t/(tau0*p*(1+al^2)), -al^2p^2]
      !where am[.,.] is the Jacobi amplitude
      !Also, we have phi(t) explicitly depends on theta(t)
      !where phi(t) = Int_{0,t'} [ Cos(theta(tau))dtau/(1+al^2) ] + phi0
      !where phi0 is a phase shift at t=0
      !Therefore, in this subroutine we need to calculate dphi
        
      
      dt = T/(1d0+alpha*alpha)/real(N-1,q2)          !---time step
      p = find_p(T,alpha,tol)                        !---T = 4τ0(1+al^2)pK(-al2p2) ---check our PRL paper!
      !print *, 'p value:',p
      
      phi = 0d0                                      !---Initialization
      dphi = calc_dphi(T,dt,0d0,0.5d0*dt,alpha,p)    !---calcualte dphi using Simpson's rule
      
      do i = 1 , N-1                                 !---counter that is running over number of images
          
         tt = (real(i,q2)-0.5d0)*dt                  !---time slice at the midpoint
         if (tt.le.(0.5d0*T/(1d0+alpha*alpha))) then !---before the magnetic moment crosses the energy barrier
            sig = 1d0
         else
            sig =-1d0
         end if
         !sig=1.0d0        
         call gscdv(tt/p,alpha*alpha*p*p,s,c,d)         !---calculate sn[.,.], cn[.,.] 
         
         b = 0.5d0/(p*(1d0+alpha*alpha))*(d+alpha*p*s)  !---amplitude of the pulse defined by Eq.(7) in our PRL paper
                                                        !---here the pulse is calculated at the midpoint as indicated by the variable "tt"
                  
         phi = phi + dphi                               !---update phi
         
         dphi = calc_dphi(T,dt,(real(i,q2)-0.5d0)*dt,(real(i,q2)+0.5d0)*dt,alpha,p)
         
         !print *,'phi:',phi

         bopt(1,i) = b*(sig*alpha*dcos(phi)*dsqrt(0.5d0*(1d0+c))-dsin(phi))
         bopt(2,i) = b*(sig*alpha*dsin(phi)*dsqrt(0.5d0*(1d0+c))+dcos(phi))
         bopt(3,i) = -b*(alpha*dsqrt(0.5d0*(1d0-c)))
         bopt(4,i) = b


                
      end do

      !---------------------------------------------------------------------------------------
      !-------------------------------------------Notes---------------------------------------
      !The optimal switching field is defined according to Eq.(6) in our PRL paper as follows:
      !vec{bm} = b/(sqrt(1+al^2))*[al Hat{e_theta} + Hat{e_phi}]
      !Hat{e_theta} = cos(theta)cos(phi)Hat{i} + cos(theta)sin(phi)Hat{j} - sin(theta)Hat{k}
      !Hat{e_phi} = -sin(phi)Hat{i} + cos(phi)Hat{j}
      !Hat{e_r} = sin(theta)cos(phi)Hat{i} + sin(theta)sin(phi)Hat{j} + cos(theta)Hat{k}
      !In the above equations of "bopt", we have dsqrt(0.5d0*(1d0+c)) which is equivilant
      !to cos(theta) where theta=0.5*am[t/(tau0*p*(1+al^2)), -al^2p^2]
      !Also, we have dsqrt(0.5d0*(1d0-c)) which is equivilant to sin(theta)
      !---------------------------------------------------------------------------------------
      
   end subroutine calc_bopt
   
   
   !********************************************
   !    phi calculations usig Simpson's rule
   !********************************************   
   real(q2) function calc_dphi(T,dt,a,b,alpha,p)
   implicit none
      real(q2), intent(in) :: T,dt,a,b,alpha,p
      real(q2) :: s,d,ca,cb,cc, sig1,sig2,sig3, k1,k2,k3
      
      !--------------------------------------------------
      !phi(t) = Int_{0,t'} [ Cos(theta(tau))dtau/(1+al^2) ] + phi0
      !t' belongs to the time interval [0,T]
      !phi0 is considered to be zero
      !In general, Simpson's rule is given by the following formula:
      !(b-a)/6*[f(a) + 4*f( (a+b)/2) ) + f(b)]
      
      
      if (a.le.(0.5d0*T/(1d0+alpha*alpha))) then
            sig1 = 1d0
         else
            sig1 =-1d0
      end if
      if (b.le.(0.5d0*T/(1d0+alpha*alpha))) then
            sig3 = 1d0
         else
            sig3 =-1d0
      end if
      
      if ((0.5d0*(a+b)).le.(0.5d0*T/(1d0+alpha*alpha))) then
            sig2 = 1d0
         else
            sig2 =-1d0
      end if
      
      call gscdv(a/p,alpha*alpha*p*p,s,ca,d)              !---calculate sn[.,.], cn[.,.] at the end point a
      call gscdv(b/p,alpha*alpha*p*p,s,cb,d)              !---calculate sn[.,.], cn[.,.] at the end point b
      call gscdv(0.5d0*(a+b)/p,alpha*alpha*p*p,s,cc,d)    !---calculate sn[.,.], cn[.,.] at the midpoint (a+b)/2
      
      k1 = sig1*dsqrt(0.5d0*(1d0+ca))                     !---Cos(theta) at the end point a, i.e., f(a)
      k2 = sig2*dsqrt(0.5d0*(1d0+cc))                     !---Cos(theta) at the end point a, i.e., f(b)
      k3 = sig3*dsqrt(0.5d0*(1d0+cb))                     !---Cos(theta) at the midpoint, i.e., f[(a+b)/2]
      
      calc_dphi = dt/6d0*(k1+4d0*k2+k3)                   !---approximate solution for phi(t) = Int_{0,t'} [ Cos(theta(tau))dtau/(1+al^2) ] + phi0
        
   end function calc_dphi



   !********************************************
   !            Equation of motion
   !********************************************     
   subroutine LLG(N,Nc,k1,e1,psi,e2,alpha,Te,T,bopt,traj)

   implicit none
   integer, intent(in) :: N,Nc
   real(q2), intent(in) :: k1,e1(3),psi,e2(3),alpha,Te,T
   real(q2), intent(in) :: bopt(3,N-1)   
   real(q2), intent(inout) :: traj(3,N,Nc)
   integer :: i,j
   real(q2) :: dt,adt,D,func,theta,phi,xi(3)
   TYPE (VSL_STREAM_STATE) :: stream
 
   integer(kind=4) :: errcode
   integer :: brng,method,seed

   !---------------------------------------
   !In the following:"brng", "method", and
   !"seed" are needed because we will use function from
   !the mkl-library, vdrnggaussian() function, for generating
   !random number with a given mean and standard deviation  
   brng=VSL_BRNG_MT19937
   method=VSL_RNG_METHOD_GAUSSIAN_ICDF
   seed=777
   !---------------------------------------- 
  
   !***** Initializing *****
   errcode = vslnewstream(stream,brng,seed)   
   
   dt = T/(1.0d0+alpha*alpha)/real(N-1,q2)  !---time step

   adt = dsqrt(2d0*dabs(dlog(dt)))        !---check A_h given by Eq(16) in Mentink's paper 
   
   !D = alpha*Te/(1d0+alpha*alpha)         !---"D" is the strength of the fluctuation, check Eq(6) in Mentink's paper 
   D = 0.5d0*alpha*Te/(1d0+alpha*alpha)         !---"D" is the strength of the fluctuation, check Eq(6) in Mentink's paper 
   do i = 1 , Nc
      do j = 2 , N
         errcode = vdrnggaussian(method, stream, 3, xi, 0d0, 1d0) !---here number 3 means xi is a 3D vector, 0 means zero mean and 1 means unit variance
         
         call gen_xi(adt,xi) !---after generating xi using vdrnggaussian function, we choose its value depending on adt (see Eq. (16) in Mentink's paper)
                                                 
         call SIB(k1,e1,psi,e2,xi,alpha,D,dt,bopt(:,j-1),traj(:,j-1,i),traj(:,j,i)) !---Semi-implicit method (see Mentink's paper) for solving sLLG equation 
         
      end do
   
   end do
   
   errcode=vsldeletestream( stream )
   end subroutine LLG

  
   
   !********************************************
   !        Random variable distribution
   !********************************************  
   subroutine gen_xi(adt,xi)
      implicit none 
      real(q2), intent(in) :: adt
      real(q2), intent(inout) :: xi(3)
      integer :: i
      real(q2) :: aux(3)
      !------------------------------------------
      !This subroutine is based on xi(3) which is   
      !generated using vdrnggaussian function from 
      !mkl-library. However, if xi(3) is solely chosen
      !from a random number xi~N(0,1) where N(0,1) is
      !a Gaussian random ariable with zero mean and 
      !unit variance the implicit midpoint rule can
      !in general diverge!. Hence, we decide on the 
      !value of xi(3) using Eq.(16) in Mentink's paper!
      !This subroutine is an implementation of Eq(16)
      !given by Mentink's paper
         
      aux = xi 
      do i = 1 , 3
         if (aux(i) > adt) then
            xi(i) = adt
         elseif (aux(i) < -adt) then
            xi(i) = -adt
         else
            xi(i) = aux(i)
         end if
      end do
   end subroutine gen_xi

 
   
   !********************************************
   !          Semi-Implicit B solver
   !********************************************   
   subroutine SIB(k1,e1,psi,e2,xi,alpha,D,dt,bopt,momi,momf)
      implicit none
      real(q2), intent(in) :: k1,psi,e1(3),e2(3),alpha,D,dt,momi(3),xi(3),bopt(3)
      real(q2), intent(out) :: momf(3)
      real(q2) :: moma(3),beff(3),xbv(3),xbxi(3),vec(3),xbov(3)
      integer :: i
      !----------------------------------------------------------
      !This solver is based on Mentink's paper:
      !Journal of Physics: Condensed Matter 22.17 (2010): 176001.
      !In order to understand this subroutine, see the attached notes to Mentink's paper


      !first step
      call calc_beff(k1,e1,psi,e2,momi,beff)    
      call vecprod(momi,beff,xbv)
      call vecprod(momi,xi,xbxi)
      call vecprod(momi,bopt,xbov)
      
      !xbov = 0d0
      
      do i=1,3
         vec(i) = -0.5d0*dt*(beff(i)+bopt(i)+alpha*(xbv(i)+xbov(i)))+ dsqrt(dt*D*0.5d0)*(xi(i)+alpha*xbxi(i))   !---This is G vector in my notes (see Mentink's paper)
      end do
      
      call solve_imp(momi,vec,moma) !---here we find the predicted step by solving an implicit equation (Eq. 1 in my notes about SIB; attached to Mentink's paper)
      
      do i=1,3
         moma(i) = 0.5d0*(moma(i)+momi(i)) !---here we find the midpoint between Vec{s}_i and and Vec{s}_pred 
      end do
      
      
      
      !second step
      call calc_beff(k1,e1,psi,e2,moma,beff)
      
      call vecprod(moma,beff,xbv)
      call vecprod(moma,xi,xbxi)
      call vecprod(moma,bopt,xbov)
      !xbov = 0d0
      
      do i=1,3
         vec(i) = -0.5d0*dt*(beff(i)+bopt(i)+alpha*(xbv(i)+xbov(i)))+ dsqrt(dt*D*0.5d0)*(xi(i)+alpha*xbxi(i))  !---This is F vector in my notes (see Mentink's paper)
      end do
      
      call solve_imp(momi,vec,momf) !---here we find Vec{s}_{i+1} by solving an implicit equation (Eq. 2 in my notes about SIB; attached to Mentink's paper)
      
   
   end subroutine SIB
   
   
   !********************************************
   !       Solving an implicit equation
   !********************************************    
   subroutine solve_imp(x,b,y) !solving equation: y = x + a*(x+y) \times b
      implicit none
      real(q2), intent(in) :: x(3),b(3)
      real(q2), intent(out) :: y(3)
      real(q2) :: bm,bxd,bxv(3),den,num
      !-----------------------------------------------------
      !In this subroutine, we find the value of the vector y
      !which is defined by the following implicit equation:
      !y = x + a*(x+y) \times b, where x, y, and b are 3D 
      !vector. Here, the parameter a is set to unity
      !The derivation of the above implicit equation
      !is described in my notes (see Mentink's paper!)
           
      bm = dotprod(3,b,b)
      num = 1d0 - bm
      den = 1d0 + bm
      bxd = 2d0*dotprod(3,x,b)
      call vecprod(x,b,bxv)
      y(1) = (num*x(1)+bxd*b(1)+2d0*bxv(1))/den
      y(2) = (num*x(2)+bxd*b(2)+2d0*bxv(2))/den
      y(3) = (num*x(3)+bxd*b(3)+2d0*bxv(3))/den
   
   end subroutine solve_imp

   !********************************************
   !      Function defines the dot product
   !********************************************    
   real(q2) function dotprod(N,a,b) ! dotprod = a*b
      implicit none
      integer, intent(in) :: N
      real(q2), intent(in) :: a(N),b(N)
      integer :: i
      real(q2) :: s
   
      s = 0d0
      do i=1,N
         s = s + a(i)*b(i)
      end do
   
      dotprod = s
   
   end function dotprod

   !********************************************
   !      Function defines the cross product
   !********************************************     
   subroutine vecprod(a,b,c) ! c = a \times b
   implicit none
      real(q2), intent(in) :: a(3),b(3)
      real(q2), intent(out) :: c(3)
      
      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) =-a(1)*b(3) + a(3)*b(1)
      c(3) = a(1)*b(2) - a(2)*b(1)
   
   end subroutine vecprod


!********************************************
!            Energy of the system (Modified)
!********************************************
   real(q2) function cene(k1, e1, psi, e2, mom)
   implicit none
   ! Inputs:
   !   k1  = K_ani1, psi = K_ani2
   !   e1, e2 are provided for backwards compatibility (not used here)
   !   mom(3) = magnetization vector (m1, m2, m3)
   real(q2), intent(in) :: k1, psi
   real(q2), intent(in) :: e1(3), e2(3)
   real(q2), intent(in) :: mom(3)
   
   cene = k1 * ( mom(1)**2 * mom(2)**2 + mom(1)**2 * mom(3)**2 + mom(2)**2 * mom(3)**2 ) + &
          psi * ( mom(1)**2 * mom(2)**2 * mom(3)**2 )
 end function cene
 
 
 !********************************************
 !             Effective field (Modified)
 !********************************************
 subroutine calc_beff(k1, e1, psi, e2, mom, beff)
  implicit none
  ! Inputs:
  !   k1  = K_ani1, psi = K_ani2
  !   e1, e2 are provided for compatibility (not used in the field calculation)
  !   mom(3) = magnetization vector (m1, m2, m3)
  ! Output:
  !   beff(3) = effective field = -dE/dm
  real(q2), intent(in) :: k1, psi
  real(q2), intent(in) :: e1(3), e2(3)
  real(q2), intent(in) :: mom(3)
  real(q2), intent(out) :: beff(3)
  
  ! The effective field is given by the negative derivative of cene with respect to each m_i
  beff(1) = - 2.0d0 * mom(1) * ( k1 * ( mom(2)**2 + mom(3)**2 ) + psi * mom(2)**2 * mom(3)**2 )
  beff(2) = - 2.0d0 * mom(2) * ( k1 * ( mom(1)**2 + mom(3)**2 ) + psi * mom(1)**2 * mom(3)**2 )
  beff(3) = - 2.0d0 * mom(3) * ( k1 * ( mom(1)**2 + mom(2)**2 ) + psi * mom(1)**2 * mom(2)**2 )
end subroutine calc_beff

 




   
   real(q2) function norm_vec(N,vec)
   implicit none
   integer, intent(in) :: N
   real(q2), intent(in) :: vec(N)
   real(q2) :: tmp
   integer :: i
      tmp = 0d0
      do i=1,N
         tmp = tmp + vec(i)*vec(i)
      end do
   
      norm_vec = dsqrt(tmp)
   
   end function norm_vec
   
   subroutine normalize_vec(PP,vec)
   implicit none
   integer, intent(in) :: PP
   real(q2), intent(inout) :: vec(PP)
   real(q2) :: tmp
   integer :: i
   
      tmp = norm_vec(PP,vec)
   
      do i=1,PP
         vec(i) = vec(i)/tmp
      end do
   end subroutine normalize_vec
   
 
   
   subroutine extract_tp(n,t,p)
   implicit none
   real(q2), intent(in) :: n(3)
   real(q2), intent(out) :: t,p
   
   !print *,'ok:'
   !print *,n
   t = datan2(dsqrt(n(1)*n(1)+n(2)*n(2)),n(3))
   p = datan2(n(2),n(1))
   
   if (p<0d0) then
      p = p+2d0*pi
   end if
   
   end subroutine extract_tp
   
   

   
   function calc_ang(n1,n2)
   implicit none
   real(q2), intent(in) :: n1(3), n2(3) !n1 and n2 have to be normalized
   real(q2) :: calc_ang
   real(q2) :: n(3),prod,tmp
   
      n(1) = n1(2)*n2(3)-n1(3)*n2(2)
      n(2) =-n1(1)*n2(3)+n1(3)*n2(1)
      n(3) = n1(1)*n2(2)-n1(2)*n2(1)
      
      tmp = dsqrt(n(1)*n(1)+n(2)*n(2)+n(3)*n(3))
      prod = n1(1)*n2(1)+n1(2)*n2(2)+n1(3)*n2(3)
      
      calc_ang = datan2(tmp,prod)
      
   end function calc_ang
   

subroutine spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, c )


  implicit none
!
  integer, intent(in) :: n, ibcbeg, ibcend
  real(q2), intent(in) :: t(n), y(n), ybcbeg, ybcend
  real(q2), intent(out) :: c(4,n)
!
  real(q2) :: diag(n), sub(2:n), sup(1:n-1),ypp(n),h
  integer i
  
 
  
!
!  Check.
!
  if ( n <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The number of knots must be at least 2.'
    write ( *, '(a,i6)' ) '  The input value of N = ', n
    stop
  end if

  do i = 1, n-1
    if ( t(i) >= t(i+1) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
      write ( *, '(a)' ) '  The knots must be strictly increasing, but'
      write ( *, '(a,i6,a,g14.6)' ) '  T(',  i,') = ', t(i)
      write ( *, '(a,i6,a,g14.6)' ) '  T(',i+1,') = ', t(i+1)
      stop
    end if
  end do
!
!  Set the first equation.
!
  if ( ibcbeg == 0 ) then
    ypp(1) = 0d0
    diag(1) = 1d0
    sup(1) = -1d0
  else if ( ibcbeg == 1 ) then
    ypp(1) = ( y(2) - y(1) ) / ( t(2) - t(1) ) - ybcbeg
    diag(1) = ( t(2) - t(1) ) / 3d0 
    sup(1) = ( t(2) - t(1) ) / 6d0
  else if ( ibcbeg == 2 ) then
    ypp(1) = ybcbeg
    diag(1) = 1d0
    sup(1) = 0d0
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The boundary flag IBCBEG must be 0, 1 or 2.'
    write ( *, '(a,i6)' ) '  The input value is IBCBEG = ', ibcbeg
    stop
  end if
!
!  Set the intermediate equations.
!
  do i = 2, n-1
    ypp(i) = ( y(i+1) - y(i) ) / ( t(i+1) - t(i) ) &
           - ( y(i) - y(i-1) ) / ( t(i) - t(i-1) )
    sub(i) = ( t(i) - t(i-1) ) / 6d0
    diag(i) = ( t(i+1) - t(i-1) ) / 3d0
    sup(i) = ( t(i+1) - t(i) ) / 6d0
  end do
!
!  Set the last equation.
!
  if ( ibcend == 0 ) then
    ypp(n) = 0d0
    sub(n) = -1d0
    diag(n) = 1d0
  else if ( ibcend == 1 ) then
    ypp(n) = ybcend - ( y(n) - y(n-1) ) / ( t(n) - t(n-1) )
    sub(n) = ( t(n) - t(n-1) ) / 6d0
    diag(n) = ( t(n) - t(n-1) ) / 3d0
  else if ( ibcend == 2 ) then
    ypp(n) = ybcend
    sub(n) = 0d0
    diag(n) = 1d0
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The boundary flag IBCEND must be 0, 1 or 2.'
    write ( *, '(a,i6)' ) '  The input value is IBCEND = ', ibcend
    stop
  end if
!
!  Special case:
!    N = 2, IBCBEG = IBCEND = 0.
!
  if ( n == 2 .and. ibcbeg == 0 .and. ibcend == 0 ) then

    ypp(1) = 0d0
    ypp(2) = 0d0
!
!  Solve the linear system.
!
  else

    call s3_fs ( sub, diag, sup, n, ypp, ypp )
    

  end if
	
	do i=1,n-1
		h = t(i+1) - t(i)
		c(1,i) = y(i)
		c(2,i) = (y(i+1) - y(i))/h - othree*(0.5d0*ypp(i+1) + ypp(i))*h
		c(3,i) = 0.5d0*ypp(i)
		c(4,i) = (ypp(i+1) - ypp(i) )/(6d0*h)
	end do
  		c(1,n) = y(n)
  		c(2,n) = 0d0
  		c(3,n) = 0d0
  		c(4,n) = 0d0
  
end subroutine spline_cubic_set

subroutine spline_cubic_valED(n, t, c, tval, yval, ypval)
	implicit none
	
	integer, intent(in) :: n
  real(q2), intent(in) :: t(n), c(4,n), tval
  real(q2), intent(inout) :: yval
  real(q2), intent(out) :: ypval
  

  real(q2) :: dt, h
  
  integer :: left, right

!
!  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
!  Values below T(1) or above T(N) use extrapolation.
!
  call rvec_bracket ( n, t, tval, left, right )
!
!  Evaluate the polynomial.
!
  dt = tval - t(left)
  
	
	yval = c(1,left) + dt*(c(2,left) + dt*(c(3,left) + dt*c(4,left)))
	
	ypval = c(2,left) + dt*(2d0*c(3,left) + 3d0*c(4,left)*dt)
  
  
	

end subroutine spline_cubic_valED

subroutine spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, yppval )


  implicit none
!
  integer, intent(in) :: n
  real(q2), intent(in) :: t(n), y(n), ypp(n), tval
  real(q2), intent(out) :: yppval, ypval, yval
  

  real(q2) :: dt, h
  
  integer :: left, right

!
!  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
!  Values below T(1) or above T(N) use extrapolation.
!
  call rvec_bracket ( n, t, tval, left, right )
!
!  Evaluate the polynomial.
!
  dt = tval - t(left)
  h = t(right) - t(left)

  yval = y(left) &
       + dt * ( ( y(right) - y(left) ) / h &
              - ( ypp(right) / 6d0 + ypp(left) / 3d0 ) * h &
       + dt * ( 0.5d0 * ypp(left) &
       + dt * ( ( ypp(right) - ypp(left) ) / ( 6d0 * h ) ) ) )

  ypval = ( y(right) - y(left) ) / h &
       - ( ypp(right) / 6d0 + ypp(left) / 3d0 ) * h &
       + dt * ( ypp(left) &
       + dt * ( 0.5d0 * ( ypp(right) - ypp(left) ) / h ) )

  yppval = ypp(left) + dt * ( ypp(right) - ypp(left) ) / h 
  
  
!---------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------  
  

end subroutine spline_cubic_val

subroutine coeffs_set(n,c,p)
	!n - number of degrees of freedom
	implicit none
	integer, intent(in) :: n
	real(q2), intent(in) :: c(4,n)
	real(q2), intent(out) :: p(5)
	integer :: i
	p = 0d0
	
	do i=1,n
		p(1) = p(1) + c(2,i)*c(2,i)
		p(2) = p(2) + c(2,i)*c(3,i)
		p(3) = p(3) + 2d0*c(3,i)*c(3,i) + 3d0*c(2,i)*c(4,i)
		p(4) = p(4) + c(3,i)*c(4,i)
		p(5) = p(5) + c(4,i)*c(4,i)
	end do
	
	p(2) = p(2)*4d0
	p(3) = p(3)*2d0
	p(4) = p(4)*12d0
	p(5) = p(5)*9d0
	

end subroutine coeffs_set

function curve_val(n, t, p, tval)
	implicit none
	
	integer, intent(in) :: n
	real(q2), intent(in) :: t(n), p(5,n), tval
	real(q2) :: dt,curve_val
	integer :: left,right
	
	call rvec_bracket ( n, t, tval, left, right )
	
	dt = tval - t(left)
	
	curve_val = dsqrt(p(1,left) + dt*(p(2,left) + &
					dt*(p(3,left) + dt*(p(4,left) + dt*p(5,left)))))
	

end function curve_val

subroutine path_len_val(n,t,p,pathl,prec)
	implicit none
	integer, intent(in) :: n
	real(q2), intent(in) :: t(n),p(5,n),prec
	real(q2), intent(out) :: pathl(n)
	integer :: i
	
	pathl(1) = 0d0
	do i=2,n
		!pathl(i) = pathl(i-1) + simpson_len(n,t,p,t(i-1),t(i),2,0.0000001d0)
		pathl(i) = pathl(i-1) + simpson_len(n,t,p,t(i-1),t(i),2,prec)
	end do
	

end subroutine path_len_val

function simpson_len(n,t,p,a,b,NN,prec)
		implicit none
		integer, intent(in) :: n,NN
		real(q2), intent(in) :: a,b,prec,t(n),p(5,n)
		real(q2) :: s1,s2,h,res,cri,tmp,simpson_len
		integer :: i,nknot
		
		nknot = NN
		h = (b-a)/(2*nknot)
		s1 = 0d0
		s2 = 0d0
!		curve_val(, tval)
		
		do i = 1,nknot-1
			s1 = s1 + curve_val(n,t,p,a+(2*i-1)*h)
			s2 = s2 + curve_val(n,t,p,a+2*i*h)
		end do
		
		res = othree*h*(curve_val(n,t,p,a) + curve_val(n,t,p,b) + &
				4d0*(s1+curve_val(n,t,p,a+(2*nknot-1)*h)) + 2d0*s2)
		
		cri = 1d0
		
		do while (cri>prec)
			nknot = 2d0*nknot
			h = (b-a)/(2*nknot)
			s1 = 0d0
			s2 = 0d0
			do i = 1,nknot-1
				s1 = s1 + curve_val(n,t,p,a+(2*i-1)*h)
				s2 = s2 + curve_val(n,t,p,a+2*i*h)
			end do
			
!			print *,'ok'
			tmp = othree*h*(curve_val(n,t,p,a) + curve_val(n,t,p,b) + &
					4d0*(s1+curve_val(n,t,p,a+(2*nknot-1)*h)) + 2d0*s2)
			cri = dabs(tmp-res)
			res = tmp
		end do
		
		simpson_len = res
	end function simpson_len


subroutine s3_fs ( a1, a2, a3, n, b, x )

  implicit none
!
  integer, intent(in) :: n
!
  real(q2), intent(inout) :: a1(2:n), a2(1:n), a3(1:n-1), b(n)
 
  integer i
  real(q2), intent(out):: x(n)
  real(q2) :: xmult
!
!  The diagonal entries can't be zero.
!
  do i = 1, n
    if ( a2(i) == 0d0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'S3_FS - Fatal error!'
      write ( *, '(a,i6,a)' ) '  A2(', i, ') = 0.'
      return
    end if
  end do

  do i = 2, n-1

    xmult = a1(i) / a2(i-1)
    a2(i) = a2(i) - xmult * a3(i-1)

    b(i) = b(i) - xmult * b(i-1)

  end do

  xmult = a1(n) / a2(n-1)
  a2(n) = a2(n) - xmult * a3(n-1)

  x(n) = ( b(n) - xmult * b(n-1) ) / a2(n)
  do i = n-1, 1, -1
    x(i) = ( b(i) - a3(i) * x(i+1) ) / a2(i)
  end do

  
end subroutine s3_fs



subroutine rvec_bracket ( n, x, xval, left, right )

  implicit none
!
  integer, intent(in) :: n
!
  integer i
  integer, intent(out) :: left,right
  real(q2), intent(in) :: x(n),xval
  
!
  do i = 2, n - 1

    if ( xval < x(i) ) then
      left = i - 1
      right = i
      return
    end if

   end do

  left = n - 1
  right = n

end subroutine rvec_bracket

   subroutine rvec_bracket_int( n, x, xval, left, right )
   implicit none
!
   integer, intent(in) :: n
   integer, intent(in) :: x(n)
   real(q2), intent(in) :: xval
   integer :: i
   integer, intent(out) :: left,right
  
!
      do i = 2, n - 1
         if ( xval < real(x(i),q2) ) then
            left = i - 1
            right = i
            return
         end if
      end do

      left = n - 1
      right = n

   end subroutine rvec_bracket_int

      
end module lib      
