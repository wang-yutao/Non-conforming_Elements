 program ptc_fit_wilson       
!------------------------------------------------------------------------
!      program 5.2 plane strain of an elastic solid using uniform 4-node
!      quadrilateral elements numbered in the y direction.
!      Analytical forms of km and bee matrices
!
!   整数：
!   NELS            单元数
!   NCE             x方向列数
!   NEQ             网格中的自由度数 方程个数
!   NBAND           网格半带宽
!   NN              网格中节点数
!   NIP             积分点数目
!   LOADED_NODES    施加载荷的节点数
!   NDOF            每个单元的自由度数（8）
!   NST             应力-应变矩阵的大小（3）
!   NOD             每个单元的节点数（4）
!   NODOF           每个节点的自由度数（2）
!   NDIM            问题的维数（2）
!   实数：
!   E               杨氏模量
!   V               泊松比
!   DET             雅可比矩阵的行列式值
!   实型数组：
!   DEE             应力-应变矩阵
!   POINTS          积分点
!   COORD           单元节点坐标
!   JAC             雅可比矩阵
!   DER             形函数在局部坐标系下的偏导数
!   DERIV           形函数在整体坐标系下的偏导数
!   BEE             应变-位移矩阵
!   KM              单元刚度矩阵
!   ELD             单元位移矢量
!   SIGMA           单元应力矢量
!   FUN             局部坐标系中的单元形函数
!   G_COORD         COORD的整体坐标形式
!   KV              整体刚度矩阵
!   LOADS           整体载荷（位移）矢量
!   WEIGHTS         积分点的加权系数
!	KM              单刚矩阵
!   KV              总刚矩阵
!   整型数组：
!   G               单元定位矢量
!   NUM             单元节点号矢量
!   G_NUM           全体NUM
!   G_G             全体G
!   NF              节点自由度数组
!   NR              约束节点数组
!------------------------------------------------------------------------
 use new_library  ;  use  geometry_lib ;  use vlib  ;  use my_library
 implicit none
 integer::nels,neq,nband,nn,nr,nip,nodof=2,nod=4,nst=3,ndof,loaded_nodes,  &
          i,k,iel,ndim=2
 real:: e,v,det;  character(len=15) :: element='quadrilateral'  
!----------------------------- dynamic arrays---------------------------------
 real  , allocatable  :: kb(:,:),loads(:),points(:,:),dee(:,:),coord(:,:), &
                         jac(:,:),der(:,:),deriv(:,:),weights(:),          &
                         bee(:,:),km(:,:),eld(:),sigma(:),g_coord(:,:),	   &
                         kuu(:,:),kua(:,:),kaa(:,:),                       &
                         wder(:,:),wbee(:,:),wderiv(:,:),                  &
                         smooth(:,:),sigmat(:,:),signod(:),				   &
                         ptc(:,:)
 integer, allocatable :: nf(:,:),g(:),num(:),g_num(:,:),g_g(:,:)
!---------------------------input and initialisation---------------------------
   open (10,file='Data\R25P.dat',status=    'old',action='read')
   open (11,file='Results\R25P.res',status='replace',action='write')                   
  read (10,*) nels,nn,nip,e,v 
  ndof=nod*nodof 
  allocate (nf(nodof,nn),points(nip,ndim),g(ndof),g_coord(ndim,nn),        & 
            dee(nst,nst),coord(nod,ndim),jac(ndim,ndim),g_g(ndof,nels),    &
            weights(nip),der(ndim,nod),deriv(ndim,nod),bee(nst,ndof),      &
            km(ndof,ndof),eld(ndof),sigma(nst),num(nod),g_num(nod,nels),   &
            kuu(ndof,ndof),kua(ndof,4),kaa(4,4),                           &
            wder(2,2),wbee(nst,4),wderiv(ndim,2),						   &
            smooth(4,4),sigmat(3,4),signod(4),                             &
            ptc(3,4))
  nf=1; read(10,*) nr ; if(nr>0)read(10,*)(k,nf(:,k),i=1,nr)
  call formnf(nf); neq=maxval(nf)    ;  nband = 0
  call deemat(dee,e,v)   ;  call sample(element,points,weights)               
!-------- loop the elements to set up global geometry and kdiag -------------
  ! 读取单元节点坐标
  read(10,*)(g_coord(:,i),i=1,nn)
  ! 读取单元节点号矢量
  read(10,*)(g_num(:,i),i=1,nels)
  elements_1  : do iel =1,nels
                  num=g_num(:,iel)
                  call num_to_g( num , nf, g)
                  g_g(:,iel)=g
                  if(nband<bandwidth(g))nband=bandwidth(g)
  end do elements_1               
    write(11,'(a)') "Global G "
    do k=1,nels;write(11,'(a,i5,a,8i5)')"Element",k,"       ",g_g(:,k);end do
        
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,4i5)')                                  &
                              "Element ",k,"        ",g_num(:,k); end do  
    write(11,'(2(a,i5))')                                                      &
                "There are ",neq ,"  equations and the half-bandwidth is", nband
                allocate(kb(neq,nband+1),loads(0:neq)); kb=.0   
!---------------- element stiffness integration and assembly-------------------
 ptc=0.0
 elements_2: do iel = 1 , nels
             num = g_num(: , iel);  g = g_g( : , iel)
             coord = transpose(g_coord(:,num))
             ! 数组归零
             km=0.0; kuu=0.0; kua=0.0; kaa=0.0
          gauss_pts_1: do i = 1 , nip
               ! 计算 jac
               call shape_der (der,points,i)
               call ptcfit_der (wder,points,i,coord)
               jac = matmul(der,coord) 
               det = determinant(jac)
               call invert(jac)
               deriv = matmul(jac,der)
               wderiv = matmul(jac,wder)
               ! 计算 kuu
               call beemat (bee,deriv) 
               kuu = kuu + matmul(matmul(transpose(bee),dee),bee) *det* weights(i)
               ! 计算 kua
               call ptcfit_beemat (wbee,wderiv)
               kua = kua + matmul(matmul(transpose(bee),dee),wbee) *det* weights(i)
               ! 记录小片检验矩阵
               if(iel==1) ptc=ptc+wbee*det
               ! 计算 kaa
               kaa = kaa + matmul(matmul(transpose(wbee),dee),wbee) *det* weights(i)
          end do gauss_pts_1
   ! 计算单刚矩阵 km
   call invert(kaa)
   km = kuu - matmul(matmul(kua,kaa),transpose(kua))
   call formkb (kb,km,g)
 end do elements_2    
 loads=0.0 ; read (10,*) loaded_nodes,(k,loads(nf(:,k)), i=1,loaded_nodes)
!---------------------------equation solution--------------------------------
    call cholin(kb) ;call chobac(kb,loads) 
    write(11,'(a)') "The nodal displacements Are :"
    write(11,'(a)') "Node         Displacement"
    do k=1,nn; write(11,'(i5,a,2e12.4)') k,"   ",loads(nf(:,k)); end do

!-------------------recover stresses at element Gauss-points-----------------
sigmat=0.0
elements_3:do iel = 1 , nels
               num = g_num(:,iel);  coord =transpose( g_coord(: ,num)) 
               g = g_g(:,iel)     ;    eld=loads(g)
              write(11,'(a,i5,a)')                                             &
                       "The Gauss Point stresses for element",iel,"   are :" 
            integrating_pts_2: do i = 1 , nip
                 call bee4(coord,points,i,det,bee)
                 sigma = matmul (dee,matmul(bee,eld)) 
                 write(11,'(a,i5)') "Point",i    ;   write(11,'(3e12.4)') sigma
                 ! 将单元 1 的应力存入 sigmat
                 if(iel==1) sigmat(:,i)=sigma
            end do integrating_pts_2
end do elements_3
!-------------------output stress of i and displaceemt of j-----------------
 ! 使用应力差值外推获得 i 点应力值
 call smoothmat(smooth)
 stress: do i = 1 , 3
   signod=matmul(smooth,sigmat(i,:))
   sigma(i)=signod(2)
   !.25*(1-1./sqrt(3.))*(1-1./sqrt(3.))*signod(1)&
   !        +.25*(1-1./sqrt(3.))*(1+1./sqrt(3.))*signod(2)&
   !        +.25*(1+1./sqrt(3.))*(1+1./sqrt(3.))*signod(3)&
   !        +.25*(1+1./sqrt(3.))*(1-1./sqrt(3.))*signod(4)
 end do stress
 ! 输出 i 点应力
 write(11,'(a,4e12.4)')                                                         &
         "The stress of i is ",sqrt(sigma(1)**2+sigma(2)**2+sigma(3)**2),       &
                                sigma(1),sigma(2),sigma(3)
 ! 输出 j 点位移
 write(11,'(a,e12.4)')												           &
        "The displaceemt of j is ",sqrt(loads(nf(1,nn))**2+loads(nf(2,nn))**2)
 ! 输出小片检验矩阵	全0通过
 write(11,'(a)') "Patch test condition:"
 do k=1,3; write(11,'(4e12.4)') ptc(k,:); end do
end  program ptc_fit_wilson  
