Program TER
  Use module


  integer, parameter :: nb_especes = 2
  integer, parameter :: i_Zn = 1
  integer, parameter :: i_O2 = 2

  !Déclarations
  Real(kind=PR), dimension(:,:,:),allocatable::c,c_aux
  Real(kind=PR), dimension(:),allocatable:: CL_h,CL_d,CL_b,CL_g
  Real(kind=PR):: cfl, t
  integer::nb_iter, i, j, k, nfreq, nsort

  !Lecture des constantes
  Call Lecture_fichier('C_IN.txt')



  !Initialisation
  cfl=0.4_PR
  dx=Lx/real(N2+1,pr)
  dy=Ly/real(N1+1,pr)
  !cfl*dx*dy !! à revoir
  pi=4*atan(1._PR)
  nb_iter=floor(Tf/dt)
  nsort=int(Tf)
  nfreq=1 !nb_iter /nsort

  ! On a plusieurs entités chimiques qui sont régies par des équations identiques
  allocate(c(0:N1,0:N2,nb_especes), c_aux(0:N1,0:N2,nb_especes))
  Allocate(CL_h(0:N2),CL_b(0:N2),CL_d(0:N1),CL_g(0:N1))

  !Definition de la condition initiale
  do j = 0,N2
     do i = 0,N1
        !c(i,j,i_O2) = 1._pr
        !c(i,j,i_O2) = (j+0.5_pr)*dx
        c(i,j,i_O2) = 0._PR!cos(2*pi*(0.5_pr+j)*dx)*cos(2*pi*(0.5_pr+i)*dy)
     end do
  end do

  k=0
  t=0._PR



  Do While (t<Tf)
!!$     if (mod(k,nfreq).eq.0) then
!!$        ! Ecriture de la solution exacte et experimentale dans des fichiers .txt et .dat
        Call Ecriture_fichier(t,'Sol_exacte','Sol_exp',k,c(:,:,i_O2))
!!$     endif
     ! Schéma explicite, on doit faire c <- c - dt*(-div(G*grad(c)))
     !print *,c(:,:,i_O2)
     !read *
     Call Remplissage_Vect_CL(CL_h,CL_d,CL_b,CL_g,c(:,:,i_O2))
     Call moins_div_D_grad(D_O2,c(:,:,i_O2),c_aux(:,:,i_O2),CL_type,CL_h,CL_d,CL_b,CL_g)
     c(:,:,i_O2) = c(:,:,i_O2) + dt*c_aux(:,:,i_O2)
     k = k+1
     t = k*dt
  End do





End program TER
