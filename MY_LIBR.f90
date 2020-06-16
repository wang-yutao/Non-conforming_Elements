module my_library
contains   
subroutine wilson_der(wder,points,i)! 计算原始wilson形函数导数
! N = [1-xi^2 1-eta^2]

! dN= [dN5/dxi 0       ]
!     [0       dN6/deta]
 implicit none
 integer,intent(in)::i; real,intent(in):: points(:,:)
 real,intent(out)::wder(:,:)
  real::xi,eta
  wder=0.
  xi=points(i,1)
  eta=points(i,2)
  wder(1,1)=-2*xi
  wder(2,2)=-2*eta
return
end subroutine wilson_der

subroutine wilson_beemat(wbee,wderiv)! 计算原始wilson应变矩阵
! N_bar = [1-xi^2 1-eta^2 0	     0      ]
!		  [0      0       1-xi^2 1-eta^2]

! wBee =  [dN5/dx 0      0      0     ]
!         [0      0      0      dN6/dy]
!         [0      dN6/dy dN5/dx 0     ]
 implicit none
 real,intent(in)::wderiv(:,:)
 real,intent(out)::wbee(:,:)
  real::x,y
  wbee=0.0
  x=wderiv(1,1)! dN5/dx
  y=wderiv(2,2)! dN6/dy
  wbee(1,1)=x
  wbee(3,3)=x
  wbee(2,4)=y
  wbee(3,2)=y
 return
end subroutine wilson_beemat

subroutine smoothmat(smooth)! 计算应力磨平矩阵
 implicit none
 real,intent(in out)::smooth(:,:)
  real::a,b,c
  a=1+sqrt(3.)/2.
  b=-0.5
  c=1-sqrt(3.)/2.
  smooth(1,1)=b;smooth(1,2)=c;smooth(1,3)=a;smooth(1,4)=b;
  smooth(2,1)=a;smooth(2,2)=b;smooth(2,3)=b;smooth(2,4)=c;
  smooth(3,1)=b;smooth(3,2)=a;smooth(3,3)=c;smooth(3,4)=b;
  smooth(4,1)=c;smooth(4,2)=b;smooth(4,3)=b;smooth(4,4)=a;
 return
end subroutine smoothmat
 
subroutine ptcfit_der(wder,points,i,coord)! 计算满足小片检验的形函数导数
! N = [1-xi^2-** 1-eta^2-**]

! dN= [dN5/dxi  dN6/dxi ]
!     [dN5/deta dN6/deta]
 implicit none
 integer,intent(in)::i; real,intent(in):: points(:,:),coord(:,:)
 real,intent(out)::wder(:,:)
  real::xi,eta,p_lambda(2,2),p_star(2,2),abmat(4,2),onemat(4,4),   &
        a1,a2,a3,b1,b2,b3,p(2,2),detp,j1,j2,j0
        
  onemat=1.0
  onemat(1,1)=-1.0; onemat(1,2)=-1.0
  onemat(2,2)=-1.0; onemat(2,4)=-1.0
  onemat(3,1)=-1.0; onemat(3,4)=-1.0
  
  abmat=0.25*matmul(onemat,coord)
  a1=abmat(1,1); a2=abmat(2,1); a3=abmat(3,1)
  b1=abmat(1,2); b2=abmat(2,2); b3=abmat(3,2)
  
  !j0=a1*b3-a3*b1
  !j1=a1*b2-a2*b1
  !j2=a2*b3-a3*b2
  !
  !wder=0.0
  !xi=points(i,1)
  !eta=points(i,2)
  !wder(1,1)=-2.*xi+2./3.*j1/j0; wder(1,2)=-2./3.*j1/j0
  !wder(2,1)=-2./3.*j2/j0; wder(2,2)=-2.*eta+2./3.*j2/j0
  
  p_lambda(1,1)=-8.0/3.0*b2; p_lambda(1,2)=8.0/3.0*b2
  p_lambda(2,1)=8.0/3.0*a2; p_lambda(2,2)=-8.0/3.0*a2
  
  detp=a1*b3-a3*b1
  p_star(1,1)=0.25*a1/detp; p_star(1,2)=0.25*b1/detp
  p_star(2,1)=0.25*a3/detp; p_star(2,2)=0.25*b3/detp
  
  p=matmul(p_star,p_lambda)
  
  wder=0.
  xi=points(i,1)
  eta=points(i,2)
  wder(1,1)=-2*xi-p(1,1); wder(1,2)=-p(1,2)
  wder(2,1)=-p(2,1); wder(2,2)=-2*eta-p(2,2)
return
end subroutine ptcfit_der

subroutine ptcfit_beemat(wbee,wderiv)! 计算满足小片检验的应变矩阵
! N_bar = [1-xi^2-** 1-eta^2-** 0	      0         ]
!		  [0         0          1-xi^2-** 1-eta^2-**]

! wBee =  [dN5/dx dN6/dx 0   	0     ]
!         [0      0      dN5/dy dN6/dy]
!         [dN5/dy dN6/dy dN5/dx dN6/dx]
 implicit none
 real,intent(in)::wderiv(:,:)
 real,intent(out)::wbee(:,:)
  real::x5,y5,x6,y6
  wbee=0.0
  x5=wderiv(1,1)! dN5/dx
  x6=wderiv(1,2)! dN6/dx
  y5=wderiv(2,1)! dN5/dy
  y6=wderiv(2,2)! dN6/dy
  wbee(1,1)=x5
  wbee(1,2)=x6
  wbee(3,1)=y5
  wbee(3,2)=y6
  wbee(2,3)=y5
  wbee(2,4)=y6
  wbee(3,3)=x5
  wbee(3,4)=x6
 return
end subroutine ptcfit_beemat
end module my_library