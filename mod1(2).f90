Module module
  Implicit none

  ! Déclaration des coeff physiques, et des conditions aux limites
  Integer, parameter:: PR=4
  real(kind=PR)::D_O2,D_Zn,dt,dx,dy,C1_0, C2_0,Tf,pi
  integer::N1,N2
  Real(kind=PR), dimension(:),allocatable::x, y

  ! On discrétise l'équation d_t c = div(D*grad(c)) dans un domaine
  ! rectangulaire découpé en (N1+1)*(N2+2) mailles notées
  ! K_{ij}. L'inconnue discrète est la fonction constante par morceaux
  ! sur le maillage. Elle est donnée par un vecteur C(t) =
  ! (c_ij(t))_{i=0..N1, j=0..N2} de R^{(N1+1)(N2+1)}, qui est solution
  ! de l'équation
  !
  ! |K_ij|dc_ij/dt = sum_{kl, voisin de ij} F(c_ij,c_kl) |arète_{ij|kl}|
  !
  ! Avec F(c_ij,c_kl) = D*(c_kl-c_ij)/h_{ij|kl}.
  !
  ! |K_{ij}| dc_ij/dt = (F(c_ij,c_{i+1j}) - F(c_{i-1j},c_ij)) h_y + (F(c_ij,c_{ij+1}) - F(c_{ij-1},c_ij)) h_x.
  !
  ! On programme une fonction G(D,cg,cd) = D*(cd-cg), de telle sorte que
  ! F(c_ij,c_{i+1j}) = G(D,c_ij,c_{i+1j})/h_x et F(c_ij,c_{ij+1}) =
  ! G(D,c_ij,c_{ij+1})/h_y.
  !
  ! Ensuite on programme un fonction qui -div(D*grad(.)) [C'est une
  ! fonction affine de c, elle dépend aussi des conditions aux limites],
  ! qui calcule...
Contains

  ! Calcul du flux numérique int_e (D*grad(c)).n_e, où cg etcd
  ! sont les valeurs à gauche et à droite (au sens de la normale).
  function G(D,cg,cd)
    real(kind=PR), intent(in)::D,cd,cg
    real(kind=PR)::G
    G=D*(cd-cg)
  end function G


  ! Application c \in \R^{N1*N2} --> Ac \in \R^{n1*N2} où Ac_{ij} \simeq 1/{dxdy} \int_{K_{ij}} -div(D*grad(c)).
  ! CL_type est un vecteur de taille 4, haut, droite, bas, gauche contenant le type de CL: 1 pour Neumann, 0 Pour Dirichlet
  ! CL est un vect de taille 4 contenant les CL en haut, droite, bas, gauche

  ! L'application (C,CL) \in \R^{(N1+1)(N2+1)} x R^4 --> A*c + B \in R^^{(N1+1)(N2+1)}.
  !
  ! Pour CL=0 (conditions aux limites homogènes), c'est une application
  ! linéaire A*c. C'est ce qu'on utilise pour résoudre le système
  ! linaire par uneméthode itérative.
  !
  ! Pour C=0, on retrouve le vecteur B, qui contient la contribution des
  ! conditions aux limites.
  subroutine moins_div_D_grad(D,C,AC_out,CL_type,CL_h,CL_d,CL_b,CL_g)
    real(kind=PR),intent(in)::D
    real(kind=PR), dimension(0:N1,0:N2), intent(out)::AC_out
    real(kind=PR), dimension(0:N1,0:N2), intent(in)::C
    character(len=*), dimension(4), intent (in):: CL_type
    real(kind=PR), dimension(0:N1), intent (in):: CL_h, CL_b
    real(kind=PR), dimension(0:N2), intent (in):: CL_d, CL_g
    integer::i,j
    real(kind=PR)::F
  !  real(kind=PR), dimension(N1):: Flux_h, Flux_b
  !  real(kind=PR), dimension(N2):: Flux_d, Flux_g




    AC_out=0
    ! Boucle sur les interfaces verticales (ij)|(i+1j)
    Do j = 0,N2
       i = 0  ! Bord bas
       if (CL_type(4)=='Dir') then
         AC_out(i,j)=AC_out(i,j)+CL_b(j)
       else
         AC_out(i,j)=AC_out(i,j)
       end if

       Do i =1,N1-1
         F=G(D,C(i-1,j),C(i,j))*dy/dx
        AC_out(i-1,j)=AC_out(i-1,j)+F
        AC_out(i,j)=AC_out(i,j)-F
       end Do
       i = N1 ! Bord haut
       F=G(D,C(i-1,j),C(i,j))*dy/dx
       if (CL_type(2)=='Dir') then
         AC_out(i,j)=AC_out(i,j)+CL_h(j)-F
         AC_out(i-1,j)=AC_out(i-1,j)+F
       else
         AC_out(i,j)=AC_out(i,j)
       end if
    end Do

    ! Boucle sur les interfaces horizontales (ij)|(ij+1)
    Do i = 0,N1
       j = 0 ! Bord gauche
       if (CL_type(3)=='Dir') then
         AC_out(i,j)=AC_out(i,j)+CL_g(i)
       else
         AC_out(i,j)=AC_out(i,j)
       end if

       Do j = 1,N2-1
        F=G(D,C(i,j-1),C(i,j))*dx/dy
        AC_out(i,j-1)=AC_out(i,j-1)+F
        AC_out(i,j)=AC_out(i,j)-F
       end Do
       j = N2 ! Bord droit
       F=G(D,C(i,j-1),C(i,j))*dx/dy
       if (CL_type(1)=='Dir') then
         AC_out(i,j)=AC_out(i,j)+CL_d(i)-F
         AC_out(i,j-1)=AC_out(i,j-1)+F
       else
         AC_out(i,j)=AC_out(i,j)
       end if
    end Do
    AC_out=AC_out*1./(dx*dy)
  end subroutine moins_div_D_grad



  subroutine Euler_Explicite(C,C_out,AC)
    real(kind=PR), dimension(0:N1,0:N2), intent(out)::C_out
    real(kind=PR), dimension(0:N1,0:N2), intent(in)::C,AC
    C_out(:,:)=C(:,:)+dt*AC
  end subroutine Euler_Explicite



  Subroutine Lecture_fichier(name_file)
    character(len=8), intent(in)::name_file
    OPEN(10,file=name_file,status='OLD')
    READ(10,*)
    read(10,*)C1_0
    READ(10,*)
    read(10,*)C2_0
    READ(10,*)
    read(10,*)Tf
    READ(10,*)
    read(10,*)N1 !Nombre de points selon x
    READ(10,*)
    read(10,*)N2 !Nombre de points selon y
    READ(10,*)
    read(10,*)D_O2
    READ(10,*)
    read(10,*)D_Zn
  End Subroutine Lecture_fichier

  Subroutine Remplissage_Vect_CL(CL_h,CL_d,CL_b,CL_g)
    real(kind=PR), dimension(0:N1), intent(out)::CL_h, CL_b
    real(kind=PR), dimension(0:N2), intent(out)::CL_d, CL_g
    integer:: i
    Do i=0,N1
      CL_h(i)=0.00001
      CL_b(i)=0.00001
    End do
    Do i=0,N2
      CL_d(i)=0.00001
      CL_g(i)=0.00001
    End do
  End Subroutine Remplissage_Vect_CL


  Subroutine Ecriture_fichier(t, name_exacte, name_exp,k, sol)
    real(kind=PR), intent(in)::t
    real(kind=PR), dimension(0:N1,0:N2), intent(in)::sol
    character(len=*), intent (in)::name_exp
    character(len=*), intent (in)::name_exacte
    character(len=50) :: mystr
    integer::nbre, i, j,nf
    integer, intent(in) ::k


    print*,'temps=',t
    nf=1
    nbre = k/nf
    write(mystr,'(i4.4)')nbre


    open(21,file=name_exp//trim(mystr)//'.vtk')
    open(22,file=name_exacte//trim(mystr)//'.vtk')
    open(23,file=name_exp//trim(mystr)//'.dat')
    open(24,file=name_exacte//trim(mystr)//'.dat')

    do   i=21,22
      write(i,'(a)')"# vtk DataFile Version 3.0"
      write(i,'(a)')"cell"
      write(i,'(a)')"ASCII"
      write(i,'(a)')"DATASET STRUCTURED_POINTS"
      write(i,'(a)')"DIMENSIONS 50 50 1 "!, N2 ," ", N1 ," 1"
      write(i,'(a)')"ORIGIN 0 0 0"!, 0 , " " , 0 , " " , 0
      write(i,'(a)')"SPACING 1 1 0"! , 1 , " " , 1 ," " , 1
      write(i,'(a)')"POINT_DATA 2500" !, N2*N1
      write(i,'(a)')"SCALARS cell float"
      write(i,'(a)')"LOOKUP_TABLE default"
    End do



    do i=0, N1
          write(21,*)sol(i,0:N2)
          write(22,*)cos(2*pi*x(i))*cos(2*pi*y(0:N2))*exp(-2*(2*pi)**2*D_O2*t)
          do j=0,N2
            write(23,*)x(i),y(j),sol(i,j)
            write(24,*)x(i),y(j),cos(2*pi*x(i))*cos(2*pi*y(j))*exp(-2*(2*pi)**2*D_O2*t)
        enddo
    enddo
    close(21)
    close(22)
    close(23)
    close(24)

  End Subroutine Ecriture_fichier








!  Subroutine Cm(C, C_1, D)
!    real(kind=PR), dimension(0:N1,0:N2), intent(out)::C_1
!    real(kind=PR), dimension(0:N1,0:N2), intent(in)::C
!    real(kind=PR), intent(in)::dt,dx,dy,D
!    integer::i,j
!    integer, intent(in)::N1,N2
!    C_1(0,:)=C(0,:)
!    C_1(:,0)=!!
!    C_1(N1,:)=Cinit1
!    C_1(:,N2)=!!
!    Do i=1,N1-1
!      Do j=1,N2-1
!        C_1(i,j)=C(i,j)+D*dt/(dx**2)*(C(i+1,j)-2*C(i,j)+C(i-1,j))+D*dt/(dy**2)*(C(i,j+1)-2*C(i,j)+C(i,j-1))
!       !C1(i)=C(i)*exp(-Ro(i)*t)
!      End Do
!    End do
!  End Subroutine



End module
