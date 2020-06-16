 program ptc_fit_wilson       
!------------------------------------------------------------------------
!      program 5.2 plane strain of an elastic solid using uniform 4-node
!      quadrilateral elements numbered in the y direction.
!      Analytical forms of km and bee matrices
!
!   ������
!   NELS            ��Ԫ��
!   NCE             x��������
!   NEQ             �����е����ɶ��� ���̸���
!   NBAND           ��������
!   NN              �����нڵ���
!   NIP             ���ֵ���Ŀ
!   LOADED_NODES    ʩ���غɵĽڵ���
!   NDOF            ÿ����Ԫ�����ɶ�����8��
!   NST             Ӧ��-Ӧ�����Ĵ�С��3��
!   NOD             ÿ����Ԫ�Ľڵ�����4��
!   NODOF           ÿ���ڵ�����ɶ�����2��
!   NDIM            �����ά����2��
!   ʵ����
!   E               ����ģ��
!   V               ���ɱ�
!   DET             �ſɱȾ��������ʽֵ
!   ʵ�����飺
!   DEE             Ӧ��-Ӧ�����
!   POINTS          ���ֵ�
!   COORD           ��Ԫ�ڵ�����
!   JAC             �ſɱȾ���
!   DER             �κ����ھֲ�����ϵ�µ�ƫ����
!   DERIV           �κ�������������ϵ�µ�ƫ����
!   BEE             Ӧ��-λ�ƾ���
!   KM              ��Ԫ�նȾ���
!   ELD             ��Ԫλ��ʸ��
!   SIGMA           ��ԪӦ��ʸ��
!   FUN             �ֲ�����ϵ�еĵ�Ԫ�κ���
!   G_COORD         COORD������������ʽ
!   KV              ����նȾ���
!   LOADS           �����غɣ�λ�ƣ�ʸ��
!   WEIGHTS         ���ֵ�ļ�Ȩϵ��
!	KM              ���վ���
!   KV              �ܸվ���
!   �������飺
!   G               ��Ԫ��λʸ��
!   NUM             ��Ԫ�ڵ��ʸ��
!   G_NUM           ȫ��NUM
!   G_G             ȫ��G
!   NF              �ڵ����ɶ�����
!   NR              Լ���ڵ�����
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
  ! ��ȡ��Ԫ�ڵ�����
  read(10,*)(g_coord(:,i),i=1,nn)
  ! ��ȡ��Ԫ�ڵ��ʸ��
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
             ! �������
             km=0.0; kuu=0.0; kua=0.0; kaa=0.0
          gauss_pts_1: do i = 1 , nip
               ! ���� jac
               call shape_der (der,points,i)
               call ptcfit_der (wder,points,i,coord)
               jac = matmul(der,coord) 
               det = determinant(jac)
               call invert(jac)
               deriv = matmul(jac,der)
               wderiv = matmul(jac,wder)
               ! ���� kuu
               call beemat (bee,deriv) 
               kuu = kuu + matmul(matmul(transpose(bee),dee),bee) *det* weights(i)
               ! ���� kua
               call ptcfit_beemat (wbee,wderiv)
               kua = kua + matmul(matmul(transpose(bee),dee),wbee) *det* weights(i)
               ! ��¼СƬ�������
               if(iel==1) ptc=ptc+wbee*det
               ! ���� kaa
               kaa = kaa + matmul(matmul(transpose(wbee),dee),wbee) *det* weights(i)
          end do gauss_pts_1
   ! ���㵥�վ��� km
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
                 ! ����Ԫ 1 ��Ӧ������ sigmat
                 if(iel==1) sigmat(:,i)=sigma
            end do integrating_pts_2
end do elements_3
!-------------------output stress of i and displaceemt of j-----------------
 ! ʹ��Ӧ����ֵ���ƻ�� i ��Ӧ��ֵ
 call smoothmat(smooth)
 stress: do i = 1 , 3
   signod=matmul(smooth,sigmat(i,:))
   sigma(i)=signod(2)
   !.25*(1-1./sqrt(3.))*(1-1./sqrt(3.))*signod(1)&
   !        +.25*(1-1./sqrt(3.))*(1+1./sqrt(3.))*signod(2)&
   !        +.25*(1+1./sqrt(3.))*(1+1./sqrt(3.))*signod(3)&
   !        +.25*(1+1./sqrt(3.))*(1-1./sqrt(3.))*signod(4)
 end do stress
 ! ��� i ��Ӧ��
 write(11,'(a,4e12.4)')                                                         &
         "The stress of i is ",sqrt(sigma(1)**2+sigma(2)**2+sigma(3)**2),       &
                                sigma(1),sigma(2),sigma(3)
 ! ��� j ��λ��
 write(11,'(a,e12.4)')												           &
        "The displaceemt of j is ",sqrt(loads(nf(1,nn))**2+loads(nf(2,nn))**2)
 ! ���СƬ�������	ȫ0ͨ��
 write(11,'(a)') "Patch test condition:"
 do k=1,3; write(11,'(4e12.4)') ptc(k,:); end do
end  program ptc_fit_wilson  
