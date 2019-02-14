Program TER
  Use module


  !Déclarations
  Real(kind=PR), dimension(:,:),allocatable::C1,C2, C1_1,C2_1,AC_out
  Real(kind=PR), dimension(:),allocatable:: tn, CL_h,CL_d,CL_b,CL_g
  Real(kind=PR):: cfl, t
  integer::nb_iter, i, j, k, nfreq, nsort

  !C1 := concentration de O2
  !C2 := concentration de Zn

  !Lecture des constantes
  Call Lecture_fichier('C_IN.txt')



  !Initialisation
  cfl=0.4_PR
  dx=1._PR/N1
  dy=1._PR/N2
     !cfl*dx*dy !! à revoir
  pi=4*atan(1._PR)

  nb_iter=floor(Tf/dt)
  nsort=int(Tf)
  nfreq=1 !nb_iter /nsort


  !Allocations des tailles des différents vecteurs
  Allocate(x(0:N1), y(0:N2), tn(0:nb_iter), C1(0:N1,0:N2), C2(0:N1,0:N2), C1_1(0:N1,0:N2) , C2_1(0:N1,0:N2),AC_out(0:N1,0:N2))
  Allocate(CL_h(0:N1),CL_b(0:N1),CL_d(0:N2),CL_g(0:N2))

  !Definition des vecteurs en espace et en temps
   x=(/(i*dx, i=0,N1)/)
   y=(/(i*dy, i=0,N2)/)
   tn=(/(i*dt, i=0,nb_iter)/)

  !Definition de la condition initiale
  Do i=0,N1
    Do j=0,N2
      C1(i,j)=cos(2*pi*x(i))*cos(2*pi*y(j))
    End do
  end do


  k=0
  t=0._PR
  Call Remplissage_Vect_CL(CL_h,CL_d,CL_b,CL_g)
  Do While (t<Tf)
       if (mod(k,nfreq).eq.0) then
         ! Ecriture de la solution exacte et experimentale dans des fichiers .txt et .dat
       Call Ecriture_fichier(t,'Sol_exacte','Sol_exp',k,C1)
       endif
       Call moins_div_D_grad(D_O2,C1,AC_out,CL_type,CL_h,CL_d,CL_b,CL_g)
       call Euler_Explicite(C1,C1_1,AC_out)
       C1=C1_1
       k=k+1
       t=t+dt
    End do




End program TER
