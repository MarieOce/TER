Module module
  Implicit none

  ! Déclaration des coeff physiques, et des conditions aux limites
  Integer, parameter:: PR=8
  real(kind=PR)::D_O2,D_Zn,dt,dx,dy,C1_0, C2_0,Tf,pi,k_Zn,k_Fe,k_O2,F,F_O2,Lx,Ly
  integer::N1,N2
  character(len=3), dimension(4):: CL_type
  Real(kind=PR), dimension(:),allocatable::x, y

  integer, parameter :: i_bord_h = 1
  integer, parameter :: i_bord_d = 2
  integer, parameter :: i_bord_b = 3
  integer, parameter :: i_bord_g = 4

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
    real(kind=PR), dimension(0:N2), intent (in):: CL_h, CL_b
    real(kind=PR), dimension(0:N1), intent (in):: CL_d, CL_g
    integer::i,j
    real(kind=PR)::F
  !  real(kind=PR), dimension(N1):: Flux_h, Flux_b
  !  real(kind=PR), dimension(N2):: Flux_d, Flux_g




    ! |K_{ij}| dc_ij/dt = (F(c_ij,c_{i+1j}) - F(c_{i-1j},c_ij)) h_y + (F(c_ij,c_{ij+1}) - F(c_{ij-1},c_ij)) h_x.
    !
    ! On programme une fonction G(D,cg,cd) = D*(cd-cg), de telle sorte que
    ! F(c_ij,c_{i+1j}) = G(D,c_ij,c_{i+1j})/h_x et F(c_ij,c_{ij+1}) =
    ! G(D,c_ij,c_{ij+1})/h_y.
    AC_out=0
    ! Boucle sur les interfaces horizontales (i-1,j)|(i,j)
    Do j = 0,N2

       i = 0  ! Bord bas
       if (CL_type(i_bord_b)=='Dir') then
          F = 2.0_pr*G(D,CL_b(j),C(i,j))/dy*dx
       else
          F = -CL_b(j)*dx
       end if
       AC_out(i,j) = AC_out(i,j) - F

       Do i =1,N1
          F = G(D,C(i-1,j),C(i,j))/dy*dx
          AC_out(i-1,j) = AC_out(i-1,j) + F
          AC_out(i,j)   = AC_out(i,j)   - F
       end Do

       i = N1+1 ! Bord haut
       if (CL_type(i_bord_h)=='Dir') then
          F = 2.0_pr*G(D,C(i-1,j),CL_h(j))/dy*dx
       else
          F = CL_h(j)*dx
       end if
       AC_out(i-1,j) = AC_out(i-1,j) + F

    end Do

    ! Boucle sur les interfaces verticales (i,j-1)|(i,j)
    Do i = 0,N1

       j = 0  ! Bord gauche
       if (CL_type(i_bord_g)=='Dir') then
          F = 2.0_pr*G(D,CL_g(j),C(i,j))/dx*dy
       else
          F = -CL_g(i)*dy
       end if
       AC_out(i,j) = AC_out(i,j) - F

       Do j =1,N2
          F = G(D,C(i,j-1),C(i,j))/dx*dy
          AC_out(i,j-1) = AC_out(i,j-1) + F
          AC_out(i,j)   = AC_out(i,j)   - F
       end Do

       j = N2+1 ! Bord droit
       if (CL_type(i_bord_d)=='Dir') then
          F = 2.0_pr*G(D,C(i,j-1),CL_d(i))/dx*dy
       else
          F = CL_d(i)*dy
       end if
       AC_out(i,j-1) = AC_out(i,j-1) + F

    end Do

    AC_out = AC_out / (dx*dy)

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
    read(10,*)dt
    READ(10,*)
    read(10,*)Lx
    READ(10,*)
    read(10,*)Ly
    READ(10,*)
    read(10,*)N1 !Nombre de points selon x
    READ(10,*)
    read(10,*)N2 !Nombre de points selon y
    READ(10,*)
    read(10,*)D_O2
    READ(10,*)
    read(10,*)D_Zn
    READ(10,*)
    READ(10,*)k_Zn
    read(10,*)
    READ(10,*)k_Fe
    read(10,*)
    READ(10,*)k_O2
    read(10,*)
    READ(10,*)F
    read(10,*)
    READ(10,*)F_O2
    read(10,*)
    read(10,*)CL_type(1)
    READ(10,*)
    read(10,*)CL_type(2)
    READ(10,*)
    read(10,*)CL_type(3)
    READ(10,*)
    read(10,*)CL_type(4)


  End Subroutine Lecture_fichier

  Subroutine Remplissage_Vect_CL(CL_h,CL_d,CL_b,CL_g,C_out)
    integer::a
    real(kind=PR), dimension(0:N1,0:N2), intent(in)::C_out
    real(kind=PR), dimension(0:N2), intent(out)::CL_h, CL_b
    real(kind=PR), dimension(0:N1), intent(out)::CL_d, CL_g
    integer:: i,j

    a=int(N2/4)

    Do i=0,N1
      CL_g(i)=0._pr
      CL_d(i)=0._pr
    End do
    Do j=0,N2
      CL_h(j)=F_O2*(1._PR-C_out(N1,j)/0.26_PR)!(j+0.5_pr)*dx


    End do
    Do j=0,a
      CL_b(j)=0._pr!(j+0.5_pr)*dx
    End do
    Do j=a+1,3*a
      CL_b(j)=-k_O2*exp(0.16)*C_out(0,j)!(j+0.5_pr)*dx
    End do
    Do j=3*a+1,N2
      CL_b(j)=0._pr!(j+0.5_pr)*dx
    End do
  End Subroutine Remplissage_Vect_CL


  Subroutine Ecriture_fichier(t, name_exacte, name_exp,k, sol)
    real(kind=PR), intent(in)::t
    real(kind=PR), dimension(0:N1,0:N2,2), intent(in)::sol
    character(len=*), intent (in)::name_exp
    character(len=*), intent (in)::name_exacte
    character(len=50) :: m, Nx, Ny, dim
    integer::nbre, i, j
    integer, intent(in) ::k


    print*,'temps=',t
    !nf=1
    nbre = k !/nf
    write(m,'(i4.4)')nbre
    write(Nx,'(i4)')N1+1
    write(Ny,'(i4)')N2+1
    write(dim,'(i7)')(N2+1)*(N1+1)


    open(21,file=name_exp//trim(m)//'.vtk')
  !  open(22,file=name_exacte//trim(m)//'.vtk')
    !open(23,file=name_exp//trim(m)//'.dat')
  !  open(24,file=name_exacte//trim(m)//'.dat')

    do   i=21,21
      write(i,'(a)')"# vtk DataFile Version 3.0"
      write(i,'(a)')"cell"
      write(i,'(a)')"ASCII"
      write(i,'(a)')"DATASET STRUCTURED_POINTS"
      write(i,'(a)')"DIMENSIONS"//trim(Nx)//trim(Ny)//" 1"! 50 50 1" !, N2 ," ", N1 ," 1"
      write(i,'(a)')"ORIGIN 0 0 0"!, 0 , " " , 0 , " " , 0
      write(i,'(a)')"SPACING 1 1 0"! , 1 , " " , 1 ," " , 1
      write(i,'(a)')"POINT_DATA"//trim(dim)! 900" !, N2*N1
      write(i,'(a)')"SCALARS cell float"
      write(i,'(a)')"LOOKUP_TABLE default"
    End do



    do i=0, N1
          write(21,*)sol(i,0:N2,1)

          !write(22,*)cos(2*pi*x(i))*cos(2*pi*y(0:N2))*exp(-2*(2*pi)**2*D_O2*t)
      !    do j=0,N2
        !    write(22,*)cos(2*pi*(0.5_pr+j)*dx)*cos(2*pi*(0.5_pr+i)*dy)*exp(-2*(2*pi)**2*D_O2*t)
            !write(23,*)x(i),y(j),sol(i,j,1)
            !write(24,*)x(i),y(j),cos(2*pi*x(i))*cos(2*pi*y(j))*exp(-2*(2*pi)**2*D_O2*t)
      !  enddo
    enddo
   close(21)
  !  close(22)
  !  close(23)
  !  close(24)

  End Subroutine Ecriture_fichier








End module
