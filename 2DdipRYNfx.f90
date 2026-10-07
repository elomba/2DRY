Module datatrans
  implicit none
  integer, parameter :: nsp=2
  integer, parameter :: mxNr=2500, next=2500,iterMax=30000
  real (kind=8), parameter ::  pi=3.141592653589793D0
  real (kind=8) :: W(mxNr,mxNr), BJ1sq(mxNr), root(mxNr)
End Module datatrans


Program twoDdipRY
  !     N-component  HDs+1/r**3 in 2D, from 2008's F. Lado code.
  !     Included excess two-particle energy calculation Jan. 2021. 
  !     Modified for parallel dipole-dipole interaction March 2019. 
  !     adapted and modified to F90 syntax by E. Lomba, Jan. 2018.
  !     Last FL change:  FL   20 Nov 2008   11:49 am
  !     Uses setupW, E1, Bpc2D, BpyLR, Hankel

  !     This program solves the Ornstein-Zernike equation with RHNC closure
  !     for the pair distribution functions and thermodynamics of a
  !     two-dimensional N-component charged hard disks with ln(r) potentials.
  !     Variations of the RHNC closure included are:
  !       nBref = 1:  Bref=0 (HNC closure)
  !       nBref = 2:  Bref=B(hard disk) using PC closure for Restricted Primitive Model
  !       nBref = 3:  Bref=B(short range) for arbitrary electrolyte using PY for SR part

  !     NOTE: Let FTf(q;j,k) be the 2D Fourier transform of a function f(r;j,k).
  !           Then, in the program, the computed transform is
  !                         xFTf(q;j,k) = (qMax/2.0*pi*rMax)*FTf(q;j,k).
  !           Also,           Tf(q;j,k) = sqrt(rhoHat(j)*rhoHat(k))*xFTf(q;j,k),
  !                                     = sqrt(rho(j)*rho(k))*FTf(q;j,k),
  !           where rhoHat = (2.0*pi*rMax/qMax)*rho.
  use datatrans
  implicit none
  integer :: ipiv(nsp)
  integer :: lwork
  integer :: Ncore(nsp,nsp)
  character :: FTtables*60, INfile*60, OUTfile*60, fname*20, fnameo*21
  real(kind=8) :: s11, s12, s22, ssum11, ssum12, ssum22, blend,blend0,&
       & sexsum, chempot(nsp), csr0(nsp,nsp), lambda, s11_0, s12_0,&
       & s11_q, s12_q, aconst, bconst, denom, rmsMax, c112,c122,&
       & v(nsp), work(2*nsp),sjk(nsp,nsp),ssumjk(nsp,nsp), ac=1.0d0,&
       & shrink, ulr(nsp,nsp), virlr(nsp,nsp), chempot0(nsp), deltarho=0.001
  real(kind=8) :: rho(nsp), rhoHat(nsp), z(nsp), sigma(nsp,nsp),  &
       ICinv(nsp,nsp),   TcSq(nsp,nsp), sum22(nsp,nsp), sum11(nsp,nsp), &
       sum12(nsp,nsp),  sum00(nsp,nsp), sum01(nsp,nsp), sum02(nsp,nsp), &
       r(0:mxNr), q(0:mxNr), dr(0:mxNr), dq(0:mxNr), rhoi(nsp), A10,&
       & A20, pres(-1:1), uint(-1:1), xc(-1:1), rhot(-1:1),P10, P20 
  real(kind=8) :: phi(0:mxNr,nsp,nsp), phiSR(0:mxNr,nsp,nsp),  phiLR(0:mxNr,nsp,nsp), &
       TphiLR(0:mxNr,nsp,nsp),   sSR(0:mxNr,nsp,nsp),    cSR(0:mxNr,nsp,nsp), &
       xFTcSR(0:mxNr,nsp,nsp),  TcSR(0:mxNr,nsp,nsp), xFTsSR(0:mxNr,nsp,nsp), &
       sSRnew(0:mxNr,nsp,nsp),     g(0:mxNr,nsp,nsp),   Bref(0:mxNr,nsp,nsp), &
       s0(0:mxNr,nsp,nsp),    s1(0:mxNr,nsp,nsp),     s2(0:mxNr,nsp,nsp), &
       d0(0:mxNr,nsp,nsp),    d1(0:mxNr,nsp,nsp),     d2(0:mxNr,nsp,nsp), &
       sc(0:mxNr), scint(0:mxNr), sg(nsp,nsp), sgi(0:mxNr,nsp,nsp),&
       & flr(0:mxNr), tflr(0:mxNr), dphi(0:mxNr,nsp,nsp),fint(1:mxNr)

  integer :: Nr, nBref, iStart, newW, n, i, j, k, l, NrIn, i2, i3, i4&
       &, iter, iext, info, lcount, irho, its, nxf, ixf, nrt, inr 
  real(kind=8) :: rMax, qMax, rhoTotal, rhoTotal0,rholist(100),Gamma, a0, a1, a2, a3, a4, a5&
       &, sumMat, sum0, rms, a01, a02, a12, a22, a11, c1, c2, sSRng, sumtsx,&
       & sumst, ssumx, s2xc, sumP1, sumP2, sumU, sumUb, sumP0, sumU0,&
       & sumUb0, sumX0, P1, P2, U, X, chemp, sumc, sumint,&
       & sumA1, sumA10, sumA20, TH11, TH12, TH22, ssum, sumX, eta, U0&
       &, X0, dPr, fopt, y1, y2, y3, sqmax, sq0, dsig, ersig, tol,&
       & eti, eto, fopto, fp, rmscut, rtmin, rtmax, xf(100), R01, R02, sq110, sq220, sums2, s2ex
  real(kind=8), external :: E1, mk
  logical :: unfound, osig, solve
  lwork = 2*nsp
  !     Read input parameters.
  write(fname,"('2DdipHD',i1,'c_list.dat')") nsp
  write(fnameo,"('2DdipdHD',i1,'c_out.dat')") nsp
  open (2,file=fname,status='old')
  read (2,*) Nr, nBref, iStart, newW, &
       Gamma, blend0, rmsMax, rmscut
  if (Nr > mxNr) Then
     print *, ' ** Eror Nr > ',mxNr
     stop
  Endif
  read (2,*) nrt
  read(2,*) rholist(1:nrt)
  read(2,*) nxf
  read(2,*) xf(1:nxf)
  read(2,*) z(1:nsp)
  do i=1,nsp
     read (2,*) Ncore(i,i:nsp)
  Enddo
  read (2,*) FTtables
  read (2,*) INfile
  read (2,*) OUTfile
  read (2,*) lambda, shrink, r01, R02
  r02 = r02*r01
  read (2, *) eta, dsig, tol, osig
  close (2,status='keep')
  !     Echo input parameters.
  open  (3,file=fnameo,status='unknown')
  write (3,'(4i10/11f15.7)') Nr,  nBref, iStart, newW, &
       Gamma, z(1:nsp)
  do i=1,nsp
     write (3,'(3i10)')Ncore(i,i:nsp)
  Enddo
  write (3,30) FTtables
  write (3,30) INfile
  write (3,30) OUTfile
  write (*,'(4i10/11f15.7)') Nr,  nBref, iStart, newW, &
       Gamma, z(1:nsp)
  do i=1,nsp
     write (*,'(5i10)')Ncore(i,i:nsp)
  Enddo
  write (*,30) FTtables
  write (*,30) INfile
  write (*,30) OUTfile
10 format (6i5/5f10.5)
20 format (/' PROGRAM 2DdipRY' &
       /' Iterative solution of the RY equation for 2D N-component&
       & charged 1/r^3 particles' &
       /' Input data:'/6i5/5f10.5)
30 format (a60)
  open(88,file="thermol.dat")
  write(88,"(5x,'# rho     xf         u/Nkt    p/(kT rho)      xc(0)    dP*/d rho    eta&
       &          sqmax       sq0       sq0/sqmax     packing      S2ex')")     

  !     Set up tables for Hankel transforms.
  open (15,file=FTtables,status='unknown')
  if (newW .eq. 1) then ! calculate needed FT tables; save for later reuse.
     write (*,*) ' Please wait. Creating W table...'
     call setupW(Nr)
     write (15,40) ((W(i,j), i = 1,Nr-1), j = 1,Nr-1)
     write (15,40) (BJ1sq(i), i = 1,Nr)
     write (15,40) (root(i), i = 1,Nr)
  else ! read in stored tables.
     read (15,40) ((W(i,j), i = 1,Nr-1), j = 1,Nr-1)
     read (15,40) (BJ1sq(i), i = 1,Nr)
     read (15,40) (root(i), i = 1,Nr)
  end if
  close (15,status='keep')
40 format (6f12.6)

  !     Initialize constants and arrays.
  n = Nr-1
  do i=1,nsp
     do j=i,nsp
        Ncore(j,i) = Ncore(i,j)
     Enddo
  Enddo
  rMax = root(Nr)/root(Ncore(1,1))
  qMax = root(Nr)/rMax
  r(0) = 0.0d0
  q(0) = 0.0d0
  do i = 1,Nr
     r(i) = root(i)/qMax
     q(i) = root(i)/rMax
     dr(i) = (2.0d0/qMax)/(root(i)*BJ1sq(i))
     dq(i) = (2.0d0/rMax)/(root(i)*BJ1sq(i))
  end do
  ! shrink core when needed
  Ncore(:,:) = shrink*Ncore(:,:)
  do j = 1,nsp
     do k = 1,nsp
        sigma(j,k) = root(Ncore(j,k))/qMax
        do i = 0,Ncore(j,k)-1
           g(i,j,k) = 0.0d0
        end do
     end do
  end do
  do i=1,nsp
     write (3,50) sigma(i,i:nsp)
     write (*,50) sigma(i,i:nsp)
     write (*,"(' gamma(i,j) =', 5f15.10)") (gamma*z(i)*z(j),j=1,nsp)
  Enddo

50 format (/5x, 'sigma(i,j) =',5f15.10)

  !     Construct potentials.
  flr(:) = 0
  tflr(:) = 0
  !  ac = 0
  if (istart .eq. 1) then ! from existing solution.
     open (16,file=INfile,status='old')
     read (16,70) NrIn, i2, i3, i4, &
          a1, a2, a3, a4, a5
     do j = 1,nsp
        do k = j,nsp
           read (16,40) (sSR(i,j,k), i = 0,NrIn-1)
        end do
     end do
     close (16,status='keep')
     write (3,60)
     write (*,60)
     write (3,70) NrIn, i2, i3, i4, &
          a1, a2, a3, a4, a5
     write (*,70) NrIn, i2, i3, i4, &
          a1, a2, a3, a4, a5
  else ! from SR ideal gas.
     do j = 1,nsp
        do k = j,nsp
           do i = 0,Nr
              sSR(i,j,k) = 0.0d0
           end do
        end do
     end do
  end if
60 format (/' Starting from diskfile:')
70 format (4i5/5f10.5)
  !
  ! Calculate long range functions
  !
  call lrfuncs(flr(0:Nr),tflr(0:Nr),r(0:Nr),q(0:Nr),Nr,ac)
  !
  ! Loop over densities to compute derivatives
  !
  do j = 1,nsp
     do k = 1,nsp
        if (k.eq.j) then 
           a0 = z(j)*z(k)*Gamma
        else
           a0 = lambda*z(j)*z(k)*Gamma
        endif
        ulr(j,k) = a0/r(Nr)
        virlr(j,k) = 3*a0/r(Nr)
        do i = 1,Nr
           phi(i,j,k) = a0/r(i)**3
           dphi(i,j,k) = 3*a0/r(i)**3 ! - r*du/dr 
           phiLR(i,j,k) = a0*flr(i)
           phiSR(i,j,k) = phi(i,j,k)-phiLR(i,j,k)
        end do
     end do
  end do
  do inr = 1, nrt
     rhoTotal0=rholist(inr)
     do ixf=1,nxf
        ersig = 100.0
        its = 0
        solve = .true.
        rhoi(1) = (1-xf(ixf))*rhoTotal0
        rhoi(2) = xf(ixf)*rhoTotal0
        do while (ersig > tol .and. solve)
           if (.not. osig) solve = .false.

           its= its+1
           do irho=-1,1
              rho(1:nsp)=rhoi(1:nsp)*(1+irho*deltarho/sum(rhoi(1:nsp)))
              rhoTotal = sum(rho(1:nsp))
              rhot(irho) = rhoTotal
              print *, rhot(irho)
              print *, rho(1:nsp)
              rhoHat(1:nsp) = (2.0d0*pi*rMax/qMax)*rho(1:nsp)
              do j = 1,nsp
                 do k = 1,nsp
                    a1 = sqrt(rho(j)*rho(k))*a0
                    TphiLR(0:nr,j,k) = a1*tflr(0:nr)
                 end do
              end do
              !     Start iterations ...
              write (3,80)(i,i,i=1,nsp)
              write (*,80)(i,i,i=1,nsp)
80            format (/13x, 'iter', 7x, 'rms', 8x, 5('sSR(',i1,i1,')',4x:))
              iext = next
              iter = 0
              !     Begin Picard iteration.
              rms = 1000.0
              fint(1:nr) = 1-exp(-eta*r(1:nr))
              do while(rms > rmsMax)
                 iter = iter+1
                 !     Begin iteration core 
                 do j = 1,nsp
                    do k = j,nsp
                       do i = Ncore(j,k),n
                          g(i,j,k) = exp(-phiSR(i,j,k))*(1+(exp((sSR(i,j,k))*fint(i))-1)/fint(i))
                       end do
                       g(Ncore(j,k),j,k) = 0.5d0*g(Ncore(j,k),j,k)
                       do i = 0,n
                          cSR(i,j,k) = g(i,j,k)-1.0d0-sSR(i,j,k)
                          cSR(i,k,j) = cSR(i,j,k)
                          g(i,k,j) = g(i,j,k)
                       end do
                    end do
                 end do
                 !     Transform cSR(r).
                 do j = 1,nsp
                    do k = j,nsp
                       call Hankel(n,r,dr,rMax,cSR(0,j,k),qMax,xFTcSR(0,j,k))
                       !!        rewind(300+j+k)
                       do i = 0,n
                          TcSR(i,j,k) = sqrt(rhoHat(j)*rhoHat(k))*xFTcSR(i,j,k)
                          TcSR(i,k,j) = TcSR(i,j,k)
                          xFTcSR(i,k,j) = xFTcSR(i,j,k)
                          !!           write(300+j+k,'(3f15.7)')q(i),TcSR(i,j,k),TphiLR(i,j,k)
                       end do
                    end do
                 end do
                 !!  stop
                 !     Ornstein-Zernike equation for xFTsSR(k).
                 do i = 0 ,n
                    do j = 1,nsp
                       do k = j,nsp
                          sumMat = 0.0d0
                          do l = 1,nsp
                             sumMat = sumMat+(TcSR(i,j,l)-TphiLR(i,j,l)) &
                                  *(TcSR(i,l,k)-TphiLR(i,l,k))
                          end do
                          TcSq(j,k) = sumMat
                          TcSq(k,j) = sumMat
                       end do
                    end do
                    do j=1,nsp-1
                       ICinv(j,j) = 1-(TcSR(i,j,j)-TphiLR(i,j,j))
                       do k=j+1,nsp
                          ICinv(j,k) = -(TcSR(i,j,k)-TphiLR(i,j,k))
                          ICinv(k,j) = ICinv(j,k)
                       Enddo
                    Enddo
                    ICinv(nsp,nsp) = 1-(TcSR(i,nsp,nsp)-TphiLR(i,nsp,nsp))
                    call dsytrf('U',nsp,ICinv,nsp,ipiv,work,lwork,info)
                    call dsytri('U',nsp,ICinv,nsp,ipiv,work,info)
                    do k=1,nsp
                       do l=k,nsp
                          ICinv(l,k) = ICinv(k,l)
                       Enddo
                    Enddo
                    do j = 1,nsp
                       do k = j,nsp
                          sumMat = 0.0d0
                          do l = 1,nsp
                             sumMat = sumMat+ICinv(j,l)*TcSq(l,k)
                          end do
                          xFTsSR(i,j,k) = (sumMat-TphiLR(i,j,k))/sqrt(rhoHat(j)*rhoHat(k))
                          xFTsSR(i,k,j) = xFTsSR(i,j,k)
                       end do
                    end do
                 end do
                 !     Inverse transform of xFTsSR(k).
                 do j = 1,nsp
                    do k = j,nsp
                       call Hankel(n,q,dq,qMax,xFTsSR(0,j,k),rMax,sSRnew(0,j,k))
                       sSRnew(0:n,k,j) = sSRnew(0:n,j,k)
                    end do
                 end do
                 !     End iteration core.

                 !     Calculate next s(i) iterate using Ng acceleration. See J. Chem. Phys. 61, 2680 (1974).
                 do j = 1,nsp
                    do k = j,nsp
                       sum22(j,k) = sum11(j,k)
                       sum11(j,k) = sum00(j,k)
                       sum12(j,k) = sum01(j,k)
                       sum00(j,k) = 0.0d0
                       sum01(j,k) = 0.0d0
                       sum02(j,k) = 0.0d0
                       sum0 = 0.0d0
                       do i = 0,n
                          s2(i,j,k) = s1(i,j,k)
                          s1(i,j,k) = s0(i,j,k)
                          s0(i,j,k) = sSRnew(i,j,k)
                          d2(i,j,k) = d1(i,j,k)
                          d1(i,j,k) = d0(i,j,k)
                          d0(i,j,k) = sSRnew(i,j,k)-sSR(i,j,k)
                          sum00 = sum00+d0(i,j,k)*d0(i,j,k)
                          sum01 = sum01+d0(i,j,k)*d1(i,j,k)
                          sum02 = sum02+d0(i,j,k)*d2(i,j,k)
                          sum0 = sum0+(r(i)*d0(i,j,k))**2
                       end do
                       rms = sqrt(dr(Nr)*sum0)
                       if (rms .gt. rmsCut) then
                          blend = blend0
                       else
                          blend = 1.0d0-(1.0d0-blend0)*(rms/rmsCut)**2
                       end if
                       if (iter .eq. 1) then
                          do i = 0,n
                             sSR(i,j,k) = sSR(i,j,k)+blend*(sSRnew(i,j,k)-sSR(i,j,k))
                             sSR(i,k,j) = sSR(i,j,k)
                          end do
                       else if (iter .eq. 2) then
                          a01 = sum00(j,k)-sum01(j,k)
                          a11 = sum00(j,k)-2.0d0*sum01(j,k)+sum11(j,k)
                          c1 = a01/a11
                          do i = 0,n
                             sSRng = (1.0d0-c1)*s0(i,j,k)+c1*s1(i,j,k)
                             sSR(i,j,k) = sSR(i,j,k)+blend*(sSRng-sSR(i,j,k))
                             sSR(i,k,j) = sSR(i,j,k)
                          end do
                       else if (iter .le. iterMax) then
                          a01 = sum00(j,k)-sum01(j,k)
                          a02 = sum00(j,k)-sum02(j,k)
                          a11 = sum00(j,k)-2.0d0*sum01(j,k)+sum11(j,k)
                          a22 = sum00(j,k)-2.0d0*sum02(j,k)+sum22(j,k)
                          a12 = sum00(j,k)-sum01(j,k)-sum02(j,k)+sum12(j,k)
                          c1 = (a01*a22-a02*a12)/(a11*a22-a12**2)
                          c2 = (a02*a11-a01*a12)/(a11*a22-a12**2)
                          do i = 0,n
                             sSRng = (1.0d0-c1-c2)*s0(i,j,k)+c1*s1(i,j,k)+c2*s2(i,j,k)
                             sSR(i,j,k) = sSR(i,j,k)+blend*(sSRng-sSR(i,j,k))
                             sSR(i,k,j) = sSR(i,j,k)
                          end do
                       else
                          write (*,*) '*** Too many iterations. Quitting.'
                          stop
                       end if
                    end do
                 end do
                 if (mod(iter,10) == 0) then
                    write (3,90) iter, rms, (sSR(0,i,i),i=1,nsp)
                    write (*,90) iter, rms, (sSR(0,i,i),i=1,nsp)
                 endif
90               format (10x, i6, 1pe13.2, 0pf12.4, 5f12.4)
                 !     End Ng acceleration.

                 !     Check for convergence.
                 ! check for possible extrapolation attempt.
                 if (iter .eq. iext) then
                    iext = iext+next
                    do j = 1,nsp
                       do k = j,nsp
                          do i = 0,n
                             s2(i,j,k) = s0(i,j,k)-(s0(i,j,k)-s1(i,j,k))**2 &
                                  /(s0(i,j,k)-2.0d0*s1(i,j,k)+s2(i,j,k))
                             if (r(i)*dabs(s2(i,j,k)-sSR(i,j,k)) .gt. 1.0d0) cycle
                          end do
                          do i = 0,n
                             sSR(i,j,k) = s2(i,j,k)
                          end do
                       end do
                    end do
                    write (3,100) (sSR(0,i,i), i=1,nsp)
                    write (*,100) (sSR(0,i,i), i=1,nsp)
                 end if
                 !       End extrapolation attempt.
              enddo
100           format (13x, 'ext', 13x, 5f12.4)
              !     End Picard iteration.

              !     Done! Calculate g(r;j,k). Print heading.
              do j = 1,nsp
                 do k = j,nsp
                    do i = Ncore(j,k),n
                       g(i,j,k) = exp(-phi(i,j,k))*(1+(exp((sSR(i,j,k)+phiLR(i,j,k))*fint(i))-1)/fint(i))
                       g(i,k,j) = g(i,j,k)
                    end do
                 end do
              end do

!!$           scint(:) = 0.d0
!!$           sumtsx = 0.0d0
!!$           sumst = 0.0d0
!!$           sg(:,:) =0.d0
!!$           sgi(:,:,:) = 0.0d0
!!$           do i=1, n
!!$              ssumx = 0
!!$              do j = 1,nsp
!!$                 do k = 1,nsp
!!$                    if (i == Ncore(j,k)) Then
!!$                       ssumx = ssumx + 0.5*pi*rho(j)*rho(k)*(g(i,j,k)-2)&
!!$                            &+0.5*pi*rho(j)*rho(k)*(g(i,j,k)*log(g(i,j&
!!$                            &,k)))  
!!$                       sg(j,k) = sg(j,k) + 0.5*dr(i)*r(i)*((g(i,j,k)&
!!$                            &-1.0d0)**2+1)
!!$                    else
!!$                       ssumx = ssumx + pi*rho(j)*rho(k)*(g(i,j,k)-1)
!!$                       sg(j,k) = sg(j,k) + dr(i)*r(i)*(g(i,j,k)-1.0d0)**2
!!$                    Endif
!!$                    sgi(i,j,k) = 2*pi*sg(j,k)
!!$                    if (g(i,j,k) > 0) Then
!!$                       ssumx = ssumx-pi*rho(j)*rho(k)*(g(i,j,k)*log(g(i,j&
!!$                            &,k)))
!!$                    Endif
!!$                 Enddo
!!$              Enddo
!!$              sumst = sumst -(g(i,1,1)-1)*r(i)*dr(i)
!!$              if (g(i,1,1)>0) Then
!!$                 sumst = sumst +(g(i,1,1)*log(g(i,1,1)))*r(i)*dr(i)
!!$              Endif
!!$              sc(i) = ssumx*r(i)/rhoTotal
!!$              sumtsx = sumtsx + sc(i)*dr(i)
!!$              scint(i) = sumtsx
!!$           Enddo
!!$           s2exc = -pi*rhoTotal*sumst

              !     Calculate thermodynamics.
              sumP1 = 0.0d0
              sumP2 = 0.0d0
              sumU = 0.0d0
              sumUb = 0.0d0
              sumX = 0.0d0
              lcount = 1
              do j = 1,nsp
                 do k = 1,nsp
                    sumP0 = 0.0d0
                    sumU0 = 0.0d0
                    sumUb0 = 0.0d0
                    sumX0 = 0.0d0
                    g(Ncore(j,k),j,k) = 0.5d0*g(Ncore(j,k),j,k)
                    sexsum = 0.0d0
                    !
                    ! The -1 in the virial and energy accounts for the effect
                    ! of electroneutrality
                    !
                    do i = 1,n
                       sumP0 = sumP0+dr(i)*r(i)*g(i,j,k)*dphi(i,j,k)
                       sumU0 = sumU0+dr(i)*r(i)*g(i,j,k)*phi(i,j,k)
                       sumX0 = sumX0+dr(i)*r(i)*cSR(i,j,k)
                    end do
                    g(Ncore(j,k),j,k) = 2.0d0*g(Ncore(j,k),j,k)
                    sumP2 = sumP2+rho(j)*rho(k)*(sumP0+virlr(j,k))
                    sumP1 = sumP1+rho(j)*rho(k)*sigma(j,k)**2*g(Ncore(j,k),j,k)
                    sumU = sumU+rho(j)*rho(k)*(sumU0+ulr(j,k))
                    sumX = sumX+rho(j)*rho(k)*(2*pi*sumX0-z(j)*z(k)*Gamma*tflr(0))
                 end do
              end do
              P1 = 1.0d0+(pi/(2.0d0*rhoTotal))*sumP1
              P2 = (pi/(2.0d0*rhoTotal))*sumP2
              U = (pi/rhoTotal)*(sumU-sumUb)
              X = (1.0d0-(sumX/rhoTotal))
              pres(irho) = P1+P2
              uint(irho) = U
              xc(irho) = X
              if(irho==0) Then
                 P10 = P1
                 P20 = P2
                 U0 = U
                 X0 = X
              Endif


              !     if (nBref .eq. 1) then ! calculate Helmholtz free energy in HNC approximation.
              ! Chemical potential in the HNC approx
              do j=1,nsp
                 chemp = 0.0d0
                 do k=1,nsp
                    sumc =0
                    do I=1,n
                       sumc = sumc + cSR(i,j,k)*dr(i)*r(i)
                    Enddo
                    if (j.eq.k) Then
                       csr0(j,k) = 2*pi*(sumc-z(j)*z(k)*Gamma*sqrt(ac*pi)/2)
                    else
                       csr0(j,k) = 2*pi*(sumc-z(j)*z(k)*Gamma*lambda*sqrt(ac*pi)/2)
                    endif

                    sums2 = 0
                    sumint = 0
                    do i=1, n
                       if (i.eq.Ncore(j,k)) Then
                          sumint = sumint + 0.5*(g(i,j,k)-2)*(sSR(i,j,k)+phiLR(i,j,k))*dr(i)*r(i)
                          sums2 = sums2 +  0.5*(g(i,j,k)*log(g(i,j,k))-g(i,j,k)+1)*dr(i)*r(i)
                       else
                          sumint = sumint + (g(i,j,k)-1)*(sSR(i,j,k)+phiLR(i,j,k))*dr(i)*r(i)
                          if (i > Ncore(j,k)) sums2 = sums2 + (g(i,j,k)*log(g(i,j,k))-g(i,j,k)+1)*dr(i)*r(i)
                       Endif
                    Enddo
                    chemp = chemp - rho(k)*csr0(j,k) + pi*rho(k)*sumint
                 Enddo
                 chempot(j) = chemp
                 if (irho == 0) Then
                    chempot0(j) = chemp
                    s2ex = -pi*rho(j)*rho(k)*sums2/rhoTotal0
                 endif
              Enddo
              print *, 'sums=',sums2 
              sumA1 = 0.0d0
              do j = 1,nsp
                 do k = 1,nsp
                    sumA10 = 0.0d0
                    do i = 1,n
                       sumA10 = sumA10+dr(i)*r(i)*(cSR(i,j,k) &
                            +((cSR(i,j,k)-phiLR(i,j,k))**2-(sSR(i,j,k)+phiLR(i,j,k))**2)/2.0d0)
                    end do
                    sumA1 = sumA1+rho(j)*rho(k)*sumA10
                 end do
              end do
              sumc = 0
              do j=1,nsp
                 do k=1,nsp
                    if (j.eq.k) then
                       sumc = sumc -rho(j)*rho(k)*z(j)*z(k)*sqrt(ac*pi)/2
                    else
                       sumc = sumc -rho(j)*rho(k)*lambda*z(j)*z(k)*sqrt(ac*pi)/2
                    Endif
                 Enddo
              Enddo
              sumc = -pi*sumc*Gamma/rhoTotal
              A1 = -(pi/rhoTotal)*sumA1+sumc
              sumA20 = 0.0d0
              do i = 1,n
                 TH11 = rhoHat(1)*(xFTcSR(i,1,1)+xFTsSR(i,1,1))
                 TH22 = rhoHat(2)*(xFTcSR(i,2,2)+xFTsSR(i,2,2))
                 TH12 = sqrt(rhoHat(1)*rhoHat(2))*(xFTcSR(i,1,2)+xFTsSR(i,1,2))
                 sumA20 = sumA20+dq(i)*q(i)*(log((1.0d0+TH11)*(1.0d0+TH22)-TH12**2)-(TH11+TH22))
              end do
              A2 = -sumA20/(4.0d0*pi*rhoTotal)
              if (irho==0) Then
                 A20= A2
                 A10 = A1
              Endif
           Enddo
           dPr = (rhot(1)*pres(1)-rhot(-1)*pres(-1))/(2*deltarho)
           fopt = (xc(0)-dPr)
           ersig = abs(fopt/xc(0))
           if(osig) then
              if(its.lt.2)then
                 eti = eta+dsig
              else
                 !
                 !    discrete NR
                 !
                 fp = (fopt-fopto)/(eta-eto)
                 !           if (abs(fopt/fp) > abs(2*dsig)) then
                 !              eti = eta - sign(2*dsig,fopt/fp)
                 !           else
                 eti = eta - fopt/fp
                 !           endif
              endif
              eto = eta
              fopto = fopt
              eta = eti
           endif
           Write(*,'(" iter=", i4, " Ersig =", f10.6, " eta=", f10.6," fopt=&
                &",f10.6)') iter, ersig, eta, fopt
           write(*, "(' xc0=',f10.6,' dP/drho =',f10.6)")xc(0), dPr
        enddo
        write(*, "(' drho/dP =',f12.7,' Fopt',f12.7)")1.0/dPr,fopt
        write (3,101) Gamma, (i,rhoi(i),i=1,nsp)
        write (*,101) Gamma, (i,rhoi(i),i=1,nsp)
        write (3,"(5(' z(',i1,') =',f8.4,',':))") (i,z(i),i=1,nsp)
        write (*,"(5(' z(',i1,') =',f8.4,',':))") (i,z(i),i=1,nsp)
101     format (/' Thermodynamics of a 2D 1/r^3 using the RY equation' &
             /' with Bref = 0 (HNC)' &
             /' Gamma =', f8.4, ',', 5('rho(',i1,') =', f8.4,',':))
        write (3,110) pres(0), uint(0), xc(0), P10, P20
        write (*,110) P10+P20, U0, X0, P10, P20
110     format (5x, 'pA/NkT =', f15.4, ', U/NkT =', f15.4, ', NkTX/A =', f8.4/ &
             12x, '=', f15.4, '  (HD)'/12x, '=', f15.4, '  (QQ)')
        !  write (3,120) A10, A20, A10+A20
        !  write (*,120) A10, A20, A10+A20
        ssum =0
        do j=1,nsp
           ssum = rhoi(j)*chempot0(j)+ssum
           !    write(*,"(' mu^ex(',i1,')/kT=',f10.5)")j,chempot0(j)
           !    write(3,"(' mu^ex(',i1,')/kT=',f10.5)")j,chempot0(j)
        Enddo
        ! write(*,'(" A^ex/NkT(from mu)=",f10.5)')ssum/sum(rhoi(1:nsp))-(pres(0)-1)
        ! write(3,'(" A^ex/NkT(from mu)=",f10.5)')ssum/sum(rhoi(1:nsp))-(pres(0)-1)

        open(95,file='srq.dat')

120     format (5x, 'HNC free energy: A1/NkT =', f8.5, ', A2/NkT =', f8.5, &
             ', Aex/NkT =', f8.5)

        !     Save solution.
        open  (17,file=OUTfile,status='unknown')
        write (17,70) Nr, Ncore(1,1), Ncore(2,2), nBref, &
             Gamma, rho(1), rho(2), z(1), z(2)
        do j = 1,nsp
           do k = j,nsp
              write (17,40) (sSR(i,j,k), i = 0,n)
           end do
        end do
        close (17,status='keep')
        close (3,status='keep')
        open(22,file='gr.dat')
        open(23,file='sq.dat')
        do j=1,nsp
           do k=1, nsp
              y1 = rho(j)*rho(k)/rhoTotal*(2.0d0*pi*rMax/qMax)&
                   &*(xFTcSR(1,j,k)+xFTsSR(1,j,k))
              y2 = rho(j)*rho(k)/rhoTotal*(2.0d0*pi*rMax/qMax)&
                   &*(xFTcSR(2,j,k)+xFTsSR(2,j,k))
              y3 = rho(j)*rho(k)/rhoTotal*(2.0d0*pi*rMax/qMax)&
                   &*(xFTcSR(3,j,k)+xFTsSR(3,j,k))
              if (j==k) then
                 y1 = y1+rho(j)/rhoTotal
                 y2 = y2+rho(j)/rhoTotal
                 y3 = y3+rho(j)/rhoTotal
              endif
!!$        a= -((-(r(2)*y1) + r(3)*y1 + r(1)*y2 - r(3)*y2 - r(1)*y3 + r(2)*y3)&
!!$             &/(r(1)**2*r(2) - r(1)*r(2)**2 - r(1)**2*r(3) + r(2)**2*r(3) +&
!!$             & r(1)*r(3)**2 - r(2)*r(3)**2))) 
!!$        b= -((r(2)**2*y1 - r(3)**2*y1 - r(1)**2*y2 + r(3)**2*y2 + r(1)&
!!$             &**2*y3 - r(2)**2*y3)/(r(1)**2*r(2) - r(1)*r(2)**2 - r(1)&
!!$             &**2*r(3) + r(2)**2*r(3) + r(1)*r(3)**2 - r(2)*r(3)**2))) 
              sjk(j,k)=-((-(r(2)**2*r(3)*y1) + r(2)*r(3)**2*y1 + r(1)**2*r(3)*y2 -&
                   & r(1)*r(3)**2*y2 - r(1)**2*r(2)*y3 +  r(1)*r(2)**2*y3)&
                   &/(r(1)**2*r(2) - r(1)*r(2)**2 - r(1)**2*r(3) + r(2)**2&
                   &*r(3) + r(1)*r(3)**2 - r(2)*r(3)**2))
           enddo
        enddo
        write(23,'(12f15.7)')q(0), q(0)/sqrt(rhoTotal0),((sjk(j,k),j=1,nsp),k=1,nsp),&
             & sum(sjk(:,:)),((rho(1)*mk(0.0d0,R01))**2*sjk(1,1)+(rho(2)&
             &*mk(0.0d0,R02))**2*sjk(2,2)+2*rho(1)*rho(2)*mk(0.0d0,R01)&
             &*mk(0.0d0,R02)*sjk(1,2))/rhoTotal,mk(0.0,R02)
        sqmax = 0.0
        sq0 = sum(sjk(:,:))
        sq110 = sjk(1,1)
        sq220 = sjk(2,2)
        unfound = .true.
        write(95,'(12f15.7)')q(0),(sjk(j,j)*rhoTotal/rho(j),j=1,nsp)
        do i=1, n
           do j=1,nsp
              sjk(j,j) = rho(j)/rhoTotal+rho(j)*rho(j)/rhoTotal*(2.0d0*pi&
                   &*rMax/qMax)*(xFTcSR(i,j,j)+xFTsSR(i,j,j)) 
              do k=1,nsp
                 if (j.ne.k) Then
                    sjk(j,k) = rho(j)*rho(k)/rhoTotal*(2.0d0*pi*rMax/qMax)&
                         &*(xFTcSR(i,j,k)+xFTsSR(i,j,k))
                 Endif
              Enddo
           Enddo
           write(22,'(18f15.7)')r(i),((g(i,j,k),j=1,nsp),k=1,nsp)
           write(23,'(12f15.7)')q(i), q(i)/sqrt(rhoTotal0),((sjk(j,k),j=1,nsp),k=1,nsp),&
                & sum(sjk(:,:)),((rho(1)*mk(q(i),R01))**2*sjk(1,1)+(rho(2)&
                &*mk(q(i),R02))**2*sjk(2,2)+2*rho(1)*rho(2)*mk(q(i),R01)&
                &*mk(q(i),R02)*sjk(1,2))/rhoTotal,mk(q(i),R02)
           if (sum(sjk(:,:)) > 1) then
              if (sum(sjk(:,:)) > sqmax .and. unfound) then
                 sqmax = sum(sjk(:,:))
              else
                 unfound = .false.
              endif
           endif
           write(95,'(12f15.7)')q(i),(sjk(j,j)*rhoTotal/rho(j),j=1,nsp)
        Enddo
        close(22)
        close(23)
        close(95)
        write(88,"(16f12.6)")rhoTotal0,xf(ixf),uint(0), pres(0), xc(0), dPr, eta&
             &, sqmax, sq0, sq0/sqmax,(sq220/(sq110))**0.25, s2ex 
        R02 = R01/((sq220)/(sq110))**0.25
        print *, ' Optimum R01/R02=', R01/R02
        open(230,file='sqopt.dat')
        do j=1,nsp
           do k=1, nsp
              y1 = rho(j)*rho(k)/rhoTotal*(2.0d0*pi*rMax/qMax)&
                   &*(xFTcSR(1,j,k)+xFTsSR(1,j,k))
              y2 = rho(j)*rho(k)/rhoTotal*(2.0d0*pi*rMax/qMax)&
                   &*(xFTcSR(2,j,k)+xFTsSR(2,j,k))
              y3 = rho(j)*rho(k)/rhoTotal*(2.0d0*pi*rMax/qMax)&
                   &*(xFTcSR(3,j,k)+xFTsSR(3,j,k))
              if (j==k) then
                 y1 = y1+rho(j)/rhoTotal
                 y2 = y2+rho(j)/rhoTotal
                 y3 = y3+rho(j)/rhoTotal
              endif
!!$        a= -((-(r(2)*y1) + r(3)*y1 + r(1)*y2 - r(3)*y2 - r(1)*y3 + r(2)*y3)&
!!$             &/(r(1)**2*r(2) - r(1)*r(2)**2 - r(1)**2*r(3) + r(2)**2*r(3) +&
!!$             & r(1)*r(3)**2 - r(2)*r(3)**2))) 
!!$        b= -((r(2)**2*y1 - r(3)**2*y1 - r(1)**2*y2 + r(3)**2*y2 + r(1)&
!!$             &**2*y3 - r(2)**2*y3)/(r(1)**2*r(2) - r(1)*r(2)**2 - r(1)&
!!$             &**2*r(3) + r(2)**2*r(3) + r(1)*r(3)**2 - r(2)*r(3)**2))) 
              sjk(j,k)=-((-(r(2)**2*r(3)*y1) + r(2)*r(3)**2*y1 + r(1)**2*r(3)*y2 -&
                   & r(1)*r(3)**2*y2 - r(1)**2*r(2)*y3 +  r(1)*r(2)**2*y3)&
                   &/(r(1)**2*r(2) - r(1)*r(2)**2 - r(1)**2*r(3) + r(2)**2&
                   &*r(3) + r(1)*r(3)**2 - r(2)*r(3)**2))
           enddo
        enddo
        write(230,'(12f15.7)')q(0), q(0)/sqrt(rhoTotal0), ((sjk(j,k),j=1,nsp),k=1,nsp),&
             & sum(sjk(:,:)),((rho(1)*mk(0.0d0,R01))**2*sjk(1,1)+(rho(2)&
             &*mk(0.0d0,R02))**2*sjk(2,2)+2*rho(1)*rho(2)*mk(0.0d0,R01)&
             &*mk(0.0d0,R02)*sjk(1,2))/rhoTotal,mk(0.0,R02)
        sqmax = 0.0
        sq0 = sum(sjk(:,:))
        sq110 = sjk(1,1)
        sq220 = sjk(2,2)
        unfound = .true.
        do i=1, n
           do j=1,nsp
              sjk(j,j) = rho(j)/rhoTotal+rho(j)*rho(j)/rhoTotal*(2.0d0*pi&
                   &*rMax/qMax)*(xFTcSR(i,j,j)+xFTsSR(i,j,j)) 
              do k=1,nsp
                 if (j.ne.k) Then
                    sjk(j,k) = rho(j)*rho(k)/rhoTotal*(2.0d0*pi*rMax/qMax)&
                         &*(xFTcSR(i,j,k)+xFTsSR(i,j,k))
                 Endif
              Enddo
           Enddo
           write(230,'(12f15.7)')q(i), q(i)/sqrt(rhoTotal0),((sjk(j,k),j=1,nsp),k=1,nsp),&
                & sum(sjk(:,:)),((rho(1)*mk(q(i),R01))**2*sjk(1,1)+(rho(2)&
                &*mk(q(i),R02))**2*sjk(2,2)+2*rho(1)*rho(2)*mk(q(i),R01)&
                &*mk(q(i),R02)*sjk(1,2))/rhoTotal,mk(q(i),R02)
           if (sum(sjk(:,:)) > 1) then
              if (sum(sjk(:,:)) > sqmax .and. unfound) then
                 sqmax = sum(sjk(:,:))
              else
                 unfound = .false.
              endif
           endif
        Enddo
        close(230)


     enddo
  Enddo
end program twoDdipRY




SUBROUTINE BesChb(x,gam1,gam2,gampl,gammi)
  !     Uses Chebev
  !     Taken from "Numerical Recipes" by W.H. Press et al.
  implicit none
  real(kind=8), save :: c1(7)=(/-1.142022680371172d0,6.516511267076d-3&
       &,3.08709017308d-4,-3.470626964d-6,6.943764d-9,3.6780d-11,&
       &-1.36d-13/) , c2(8)=(/1.843740587300906d0,-.076852840844786d0,1.271927136655d-3, &
       &-4.971736704d-6,-3.3126120d-8,2.42310d-10,-1.70d-13,-1.d-15/)
  integer, parameter :: nuse1=5,nuse2=5
  real(kind=8) :: x, xx, gam1, gam2, gammi, gampl
  real(kind=8), external :: Chebev
  xx = 8.0d0*x*x-1.0d0
  gam1 = Chebev(-1.0d0,1.0d0,c1,nuse1,xx)
  gam2 = Chebev(-1.0d0,1.0d0,c2,nuse2,xx)
  gampl = gam2-x*gam1
  gammi = gam2+x*gam1
end subroutine BesChb



   SUBROUTINE BessJY(x,xnu,rj,ry,rjp,ryp)
     !     Uses BesChb
     !     Taken from "Numerical Recipes" by W.H. Press et al.
     !     Modified to calculate J0(x) and J0'(x) by direct series summation
     !     for x<2 to get double precision in this range.
     implicit real*8(a-h,o-z), integer(i-n)
     parameter (eps=1.0d-16,fpmin=1.0d-30,maxit=10000,xmin=2.0d0)

     pi = dacos(-1.0d0)
     if (x.le.0.0d0 .or. xnu.lt.0.0d0) then
        write (*,*) '*** Bad arguments in BessJY'
        stop
     end if
     if(x .lt. xmin) then
        nl = idint(xnu+0.5d0)
     else
        nl = max0(0,idint(xnu-x+1.5d0))
     end if
     xmu = xnu-dble(nl)
     xmu2 = xmu*xmu
     xi = 1.0d0/x
     xi2 = 2.0d0*xi
     w = xi2/pi
     isign = 1
     h = xnu*xi
     if (h.lt.fpmin) h = fpmin
     b = xi2*xnu
     d = 0.0d0
     c = h
     do i = 1,maxit
        b = b+xi2
        d = b-d
        if (dabs(d) .lt. fpmin) d=fpmin
        c = b-1.0d0/c
        if (dabs(c) .lt. fpmin) c=fpmin
        d = 1.0d0/d
        del = c*d
        h = del*h
        if (d .lt. 0.0d0) isign=-isign
        if (dabs(del-1.0d0) .lt. eps) go to 1
     end do
     write (*,*) '*** x too large in BessJY; try asymptotic expansion.'
     stop
1    continue
     rjl = isign*fpmin
     rjpl = h*rjl
     rjl1 = rjl
     rjp1 = rjpl
     fact = xnu*xi
     do l = nl,1,-1
        rjtemp = fact*rjl+rjpl
        fact = fact-xi
        rjpl = fact*rjtemp-rjl
        rjl = rjtemp
     end do
     if (rjl .eq. 0.0d0) rjl=eps
     f = rjpl/rjl
     if (x .lt. xmin) then
        x2 = 0.5d0*x
        pimu = pi*xmu
        if (dabs(pimu) .lt. eps) then
           fact = 1.0d0
        else
           fact = pimu/sin(pimu)
        end if
        d = -log(x2)
        e = xmu*d
        if (dabs(e) .lt. eps) then
           fact2 = 1.0d0
        else
           fact2 = dsinh(e)/e
        end if
        call BesChb(xmu,gam1,gam2,gampl,gammi)
        ff = 2.0d0/pi*fact*(gam1*cosh(e)+gam2*fact2*d)
        e = dexp(e)
        p = e/(gampl*pi)
        q = 1.0d0/(e*pi*gammi)
        pimu2 = 0.5d0*pimu
        if (dabs(pimu2) .lt. eps) then
           fact3 = 1.0d0
        else
           fact3 = sin(pimu2)/pimu2
        end if
        r = pi*pimu2*fact3*fact3
        c = 1.0d0
        d = -x2*x2
        sum0 = ff+r*q
        sum1 = p
        do i = 1,maxit
           ff = (dble(i)*ff+p+q)/(dble(i*i)-xmu2)
           c = c*d/dble(i)
           p = p/(dble(i)-xmu)
           q = q/(dble(i)+xmu)
           del = c*(ff+r*q)
           sum0 = sum0+del
           del1 = c*p-dble(i)*del
           sum1 = sum1+del1
           if (dabs(del) .lt. (1.0d0+dabs(sum0))*eps) go to 2
        end do
        write (*,*) '*** BessY series failed to converge in Subroutine BessJY'
        stop
2       continue
        rymu = -sum0
        ry1 = -sum1*xi2
        rymup = xmu*xi*rymu-ry1
        rjmu = w/(rymup-f*rymu)
     else
        a = 0.25d0-xmu2
        p = -0.5d0*xi
        q = 1.0d0
        br = 2.0d0*x
        bi = 2.0d0
        fact = a*xi/(p*p+q*q)
        cr = br+q*fact
        ci = bi+p*fact
        den = br*br+bi*bi
        dr = br/den
        di = -bi/den
        dlr = cr*dr-ci*di
        dli = cr*di+ci*dr
        temp = p*dlr-q*dli
        q = p*dli+q*dlr
        p = temp
        do i = 2,maxit
           a = a+dble(2*(i-1))
           bi = bi+2.0d0
           dr = a*dr+br
           di = a*di+bi
           if (dabs(dr)+dabs(di) .lt. fpmin) dr = fpmin
           fact = a/(cr*cr+ci*ci)
           cr = br+cr*fact
           ci = bi-ci*fact
           if (dabs(cr)+dabs(ci) .lt. fpmin) cr = fpmin
           den = dr*dr+di*di
           dr = dr/den
           di = -di/den
           dlr = cr*dr-ci*di
           dli = cr*di+ci*dr
           temp = p*dlr-q*dli
           q = p*dli+q*dlr
           p = temp
           if (dabs(dlr-1.0d0)+dabs(dli) .lt. eps) go to 3
        end do
        write (*,*) '*** cf2 failed in Subroutine BessJY'
        stop
3       continue
        gam = (p-f)/q
        rjmu = sqrt(w/((p-f)*gam+q))
        rjmu = dsign(rjmu,rjl)
        rymu = rjmu*gam
        rymup = rymu*(p+q/gam)
        ry1 = xmu*xi*rymu-rymup
     end if
     fact = rjmu/rjl
     rj = rjl1*fact
     rjp = rjp1*fact
     do i = 1,nl
        rytemp = (xmu+dble(i))*xi2*ry1-rymu
        rymu = ry1
        ry1 = rytemp
     end do
     ry = rymu
     ryp = xnu*xi*rymu-ry1

     !     Recalculate J0(x) and J0'(x) for small x by direct series summation.
     if (xnu .eq. 0.0d0 .and. x .lt. 2.0d0) then
        xsq = -0.25d0*x*x
        k = 0
        add0 = 1.0d0
        add1 = 1.0d0
        rj  = 1.0d0
        rjp = 1.0d0
4       k = k+1
        add0 = add0*xsq/dble(k*k)
        add1 = add1*xsq/dble(k*(k+1))
        rj  = rj +add0
        rjp = rjp+add1
        if (dabs(add0) .gt. 1.0d-15) go to 4
        rjp = -0.5d0*x*rjp
     end if
   end subroutine BessJY



   SUBROUTINE Bpc2D(Nr,Ncore,rho,r,q,dr,dq,Bpc)
     !     Uses Hankel
     !     Obtain pressure-consistent (PC) bridge function for hard disks.
     !     See J. Chem. Phys. 49, 3092 (1968).
     use datatrans, only : mxNR, pi
     implicit none
     integer :: n, Nr, Ncore, i, k, iter
     real(kind=8) :: mu, eta, rMax, qMax, rhoHat, a, bigs, biga,&
          & alpha, beta, c0, c1, c2, c0p, x, d11, d22, Ccore, rho, d00,&
          & d12, d01, d02, sum0, rms, dd01, dd02, dd11, dd22, dd12,&
          & Sng, Pr
     real(kind=8) :: rmsMax=1.0d-5,blend=0.8d0
     real(kind=8) :: r(0:mxNr), q(0:mxNr), dr(0:mxNr), dq(0:mxNr), Bpc(0:mxNr)
     real(kind=8) :: s(0:mxNr),    c(0:mxNr),    P(0:mxNr), xFTc(0:mxNr), &
          xFTp(0:mxNr), xFTs(0:mxNr), snew(0:mxNr), &
          f0(0:mxNr),   f1(0:mxNr),   f2(0:mxNr), &
          d0(0:mxNr),   d1(0:mxNr),   d2(0:mxNr)

     !     Initialize parameters and arrays.
     eta = pi*rho/4.0d0
     n = Nr-1
     rMax = r(Nr)
     qMax = q(Nr)
     rhoHat = (2.0d0*pi*rMax/qMax)*rho
     write (3,10) rho, Ncore, rMax, qMax
     write (*,10) rho, Ncore, rMax, qMax
10   format (/3x, 'SUBROUTINE Bpc2D with rho =', f8.5/ &
          3x, 'Ncore =', i4, ', rMax =', f10.5, ', qMax =', f10.5)

     !     Start from modeled PY c(r). See Leutheuser, J. Chem. Phys. 84, 1050 (1986).
     a = 1.0d0+eta
     bigs = 0.5d0*sqrt(1.0d0-a*a/4.0d0)+dasin(a/2.0d0)/a
     biga = (bigs-(1.0d0-a*a/4.0d0)*sqrt(1.0d0-a*a/4.0d0))*16.0d0/(pi*a*a)
     alpha = 2.0d0*eta**2*biga
     beta = 8.0d0*eta*bigs/pi
     c1 = ((4.0d0*eta-1.0d0)+sqrt((4.0d0*eta-1.0d0)**2 &
          -4.0d0*(alpha-beta)))/(2.0d0*(alpha-beta))
     c0 = c1-beta*c1**2
     c0p = 8.0d0*eta*c1**2/pi
     do i = 0,Ncore
        x = a*r(i)/2.0d0
        c(i) = c0+c0p*(dasin(x)+x*sqrt(1.0d0-x*x))/a
     end do
     c(Ncore) = 0.5d0*c(Ncore)
     do i = Ncore+1,n
        c(i) = 0.0d0
     end do
     !     Compute Hankel transform of c(r).
     call Hankel(n,r,dr,rMax,c,qMax,xFTc)
     !     OZ equation for xFTs(k).
     do k = 1,n
        xFTs(k) = rhoHat*xFTc(k)**2/(1.0d0-rhoHat*xFTc(k))
     end do
     !     Compute inverse Hankel transform of xFTs(k) to get modeled s(r).
     s(0) = -1.0d0-c(0)
     call Hankel(n,q,dq,qMax,xFTs,rMax,s)

     iter = 0
     mu = 0.0920d0+0.1222d0*rho+0.1642d0*rho**2+0.1100d0*rho**3
     !     Begin Picard iteration for S(r) = g(r)exp[beta*phi(r)]-1.
     !     Begin iteration core.
     rms = 100.0
     do while (rms > rmsMax)
        iter = iter+1
        do i = 0,Ncore
           P(i) = S(i)-log(1.0d0+S(i))
           C(i) = -(1.0+S(i))+mu*P(i)
        end do
        Ccore = C(Ncore)
        do i = Ncore,n
           P(i) = S(i)-log(1.0d0+S(i))
           C(i) = mu*P(i)
        end do
        C(Ncore) = 0.5d0*(Ccore+C(Ncore))
        !     Compute Hankel transform of C(r) and P(r).
        call Hankel(n,r,dr,rMax,C,qMax,xFTc)
        call Hankel(n,r,dr,rMax,P,qMax,xFTp)
        !     OZ equation for xFTs(k).
        do k = 0,n
           xFTs(k) = rhoHat*xFTc(k)**2/(1.0d0-rhoHat*xFTc(k))+mu*xFTp(k)
        end do
        !     Compute inverse Hankel transform of xFTs(k).
        call Hankel(n,q,dq,qMax,xFTs,rMax,Snew)
        !     End iteration core.

        !     Check convergence and compute next S iterate using Ng acceleration.
        d22 = d11
        d11 = d00
        d12 = d01
        d00 = 0.0d0
        d01 = 0.0d0
        d02 = 0.0d0
        sum0 = 0.0d0
        do i = 0,n
           f2(i) = f1(i)
           f1(i) = f0(i)
           f0(i) = Snew(i)
           d2(i) = d1(i)
           d1(i) = d0(i)
           d0(i) = Snew(i)-S(i)
           d00 = d00+d0(i)*d0(i)
           d01 = d01+d0(i)*d1(i)
           d02 = d02+d0(i)*d2(i)
           sum0 = sum0+d0(i)*d0(i)
        end do
        rms = sqrt(sum0/dble(n))
        !     write (3,30) iter, rms, Snew(1)
        !     write (*,30) iter, rms, Snew(1)
        !  30 format (i6, e13.2, f12.5)
        dd01 = d00-d01
        dd02 = d00-d02
        dd11 = d00-2.0d0*d01+d11
        dd22 = d00-2.0d0*d02+d22
        dd12 = d00-d01-d02+d12
        if (iter .ge. 3) then
           c1 = (dd01*dd22-dd02*dd12)/(dd11*dd22-dd12**2)
           c2 = (dd02*dd11-dd01*dd12)/(dd11*dd22-dd12**2)
        else
           c1 = dd01/dd11
           c2 = 0.0d0
        end if
        if (iter .eq. 1) then
           do i = 0,n
              S(i) = S(i)+blend*(Snew(i)-S(i))
           end do
        else 
           do i = 0,n
              Sng = f0(i)+c1*(f1(i)-f0(i))+c2*(f2(i)-f0(i))
              S(i) = S(i)+blend*(Sng-S(i))
           end do
        end if
     Enddo
     !     End Picard iteration

     !     Done! Get hard disk thermodynamics.
     Pr = 1.0d0+pi*rho*(1.0d0+S(Ncore))/2.0d0
     sum0 = 0.0d0
     do i = 1,n
        sum0 = sum0+dr(i)*r(i)*C(i)
     end do
     X = 1.0d0/(1.0d0 -2.0d0*pi*rho*sum0)
     write (3,40) rho, eta, mu
     write (*,40) rho, eta, mu
     write (3,50) iter, Pr, X
     write (*,50) iter, Pr, X
40   format (3x, 'rho =', f7.4, ', eta =', f7.4, ', mu =', f7.4/ &
          3x, 'Hard disk thermodynamics from the PC equation.')
50   format (3x, i3, ' iterations for S(r): pA/NkT =', f8.4, ', NkTX/A =', f8.4)

     !     Calculate PC bridge function of hard disks.
     do i = 0,n
        Bpc(i) = (1.0d0-mu)*(log(1.0d0+S(i))-S(i))
     end do
   end subroutine Bpc2D



   SUBROUTINE BpySR(Nr,Ncore,Gamma,sigma,z,rho,lambda, r,q,dr,dq,phiSR0,BpySR0)
     !     Uses Hankel
     !     Obtain the Percus-Yevick bridge functions for a "short-range electrolyte"
     !     of hard disks with short-range electrostatic potentials phiSR0(r;j,k).
     use datatrans, only : pi, mxNr, nsp
     implicit none
     real(kind=8), parameter :: rmsMax=1.0d-5,rmsCut=1.0d-2,blend0=0.5d0
     integer, parameter :: next=25,iterMax=1500
     integer ::  Ncore(nsp,nsp)
     integer :: i, j, k, l, info, iter, iext, Nr, n
     real(kind=8) :: rMax, qMax, rhoTotal, sumMat, sum0, rms, blend,&
          & a01, a11, a02, a22, a12, c1, c2, sSRng, Gamma, sumP1,&
          & sumP2, sumU, sumX, sumP0, sumU0, sumX0, lambda, P1, P2, U&
          &, X
     real(kind=8) :: sigma(nsp,nsp), z(nsp), rho(nsp), r(0:mxNr), q(0:mxNr),   &
          dr(0:mxNr), dq(0:mxNr), phiSR0(0:mxNr,nsp,nsp), BpySR0(0:mxNr,nsp,nsp)
     real(kind=8) :: rhoHat(nsp), ICinv(nsp,nsp),  TcSq(nsp,nsp), sum22(nsp,nsp), sum11(nsp,nsp), &
          sum12(nsp,nsp), sum00(nsp,nsp), sum01(nsp,nsp), sum02(nsp,nsp)
     real(kind=8) :: eBond(0:mxNr,nsp,nsp),    sSR(0:mxNr,nsp,nsp),    g(0:mxNr,nsp,nsp), &
          cSR(0:mxNr,nsp,nsp), xFTcSR(0:mxNr,nsp,nsp), TcSR(0:mxNr,nsp,nsp), &
          xFTsSR(0:mxNr,nsp,nsp), sSRnew(0:mxNr,nsp,nsp),                   &
          s0(0:mxNr,nsp,nsp),     s1(0:mxNr,nsp,nsp),   s2(0:mxNr,nsp,nsp), &
          d0(0:mxNr,nsp,nsp),     d1(0:mxNr,nsp,nsp),   d2(0:mxNr,nsp,nsp)
     integer :: ipiv(nsp)
     integer :: lwork
     real (kind=8) ::  work(2*nsp)
     lwork = 2*nsp
     rMax = r(Nr)
     qMax = q(Nr)
     n = Nr-1
     do j = 1,nsp
        do k = 1,nsp
           do i = 0,Ncore(j,k)-1
              eBond(i,j,k) = 0.0d0
              sSR(i,j,k) = 0.0d0
              g(i,j,k) = 0.0d0
           end do
           do i = Ncore(j,k),n
              eBond(i,j,k) = dexp(-phiSR0(i,j,k))
              sSR(i,j,k) = 0.0d0
           end do
        end do
     end do
     rhoTotal = sum(rho(1:nsp))
     rhoHat(1:nsp) = (2.0d0*pi*rMax/qMax)*rho(1:nsp)
     iext = next
     iter = 0

     !     Begin Picard iteration.
     rms = 100.0
     do while (rms > rmsMax)
        iter = iter+1
        !     Begin iteration core.
        do j = 1,nsp
           do k = j,nsp
              do i = Ncore(j,k),n
                 g(i,j,k) = eBond(i,j,k)*(1.0d0+sSR(i,j,k))
              end do
              g(Ncore(j,k),j,k) = 0.5d0*g(Ncore(j,k),j,k)
              do i = 0,n
                 cSR(i,j,k) = g(i,j,k)-1.0d0-sSR(i,j,k)
                 cSR(i,k,j) = cSR(i,j,k)
                 g(i,k,j) = g(i,j,k)
              end do
           end do
        end do
        !     Transform cSR(r).
        do j = 1,nsp
           do k = j,nsp
              call Hankel(n,r,dr,rMax,cSR(0,j,k),qMax,xFTcSR(0,j,k))
              do i = 0,n
                 TcSR(i,j,k) = sqrt(rhoHat(j)*rhoHat(k))*xFTcSR(i,j,k)
                 TcSR(i,k,j) = TcSR(i,j,k)
                 xFTcSR(i,k,j) = xFTcSR(i,j,k)
              end do
           end do
        end do
        !     Ornstein-Zernike equation for xFTsSR(k).
        do i = 1,n
           do j = 1,nsp
              do k = j,nsp
                 sumMat = 0.0d0
                 do l = 1,nsp
                    sumMat = sumMat+TcSR(i,j,l)*TcSR(i,l,k)
                 end do
                 TcSq(j,k) = sumMat
                 TcSq(k,j) = sumMat
              end do
           end do
           do j=1,nsp-1
              ICinv(j,j) = 1-TcSR(i,j,j)
              do k=j+1,nsp
                 ICinv(j,k) = -TcSR(i,j,k)
                 ICinv(k,j) = ICinv(j,k)
              Enddo
           Enddo
           ICinv(nsp,nsp) = 1-TcSR(i,nsp,nsp)
           call dsytrf('U',nsp,ICinv,nsp,ipiv,work,lwork,info)
           call dsytri('U',nsp,ICinv,nsp,ipiv,work,info)
           do k=1,nsp
              do l=k,nsp
                 ICinv(l,k) = ICinv(k,l)
              Enddo
           Enddo
           do j = 1,nsp
              do k = j,nsp
                 sumMat = 0.0d0
                 do l = 1,nsp
                    sumMat = sumMat+ICinv(j,l)*TcSq(l,k)
                 end do
                 xFTsSR(i,j,k) = sumMat/sqrt(rhoHat(j)*rhoHat(k))
                 xFTsSR(i,k,j) = xFTsSR(i,j,k)
              end do
           end do
        end do
        !     Inverse transform of xFTsSR(k).
        do j = 1,nsp
           do k = j,nsp
              call Hankel(n,q,dq,qMax,xFTsSR(0,j,k),rMax,sSRnew(0,j,k))
              sSRnew(0:n,k,j) = sSRnew(0:n,j,k)
           end do
        end do
        !     End iteration core.

        !     Calculate next s(i) iterate using Ng acceleration. See J. Chem. Phys. 61, 2680 (1974).
        do j = 1,nsp
           do k = j,nsp
              sum22(j,k) = sum11(j,k)
              sum11(j,k) = sum00(j,k)
              sum12(j,k) = sum01(j,k)
              sum00(j,k) = 0.0d0
              sum01(j,k) = 0.0d0
              sum02(j,k) = 0.0d0
              sum0 = 0.0d0
              do i = 0,n
                 s2(i,j,k) = s1(i,j,k)
                 s1(i,j,k) = s0(i,j,k)
                 s0(i,j,k) = sSRnew(i,j,k)
                 d2(i,j,k) = d1(i,j,k)
                 d1(i,j,k) = d0(i,j,k)
                 d0(i,j,k) = sSRnew(i,j,k)-sSR(i,j,k)
                 sum00 = sum00+d0(i,j,k)*d0(i,j,k)
                 sum01 = sum01+d0(i,j,k)*d1(i,j,k)
                 sum02 = sum02+d0(i,j,k)*d2(i,j,k)
                 sum0 = sum0+(r(i)*d0(i,j,k))**2
              end do
              rms = sqrt(dr(Nr)*sum0)
!!$           if (rms .gt. rmsCut) then
!!$              blend = blend0
!!$           else
!!$              blend = 1.0d0-(1.0d0-blend0)*(rms/rmsCut)**2
!!$           end if
              blend=blend0
              if (iter .eq. 1) then
                 do i = 0,n
                    sSR(i,j,k) = sSR(i,j,k)+blend*(sSRnew(i,j,k)-sSR(i,j,k))
                    sSR(i,k,j) = sSR(i,j,k)
                 end do
              else if (iter .eq. 2) then
                 a01 = sum00(j,k)-sum01(j,k)
                 a11 = sum00(j,k)-2.0d0*sum01(j,k)+sum11(j,k)
                 c1 = a01/a11
                 do i = 0,n
                    sSRng = (1.0d0-c1)*s0(i,j,k)+c1*s1(i,j,k)
                    sSR(i,j,k) = sSR(i,j,k)+blend*(sSRng-sSR(i,j,k))
                    sSR(i,k,j) = sSR(i,j,k)
                 end do
              else if (iter .le. iterMax) then
                 a01 = sum00(j,k)-sum01(j,k)
                 a02 = sum00(j,k)-sum02(j,k)
                 a11 = sum00(j,k)-2.0d0*sum01(j,k)+sum11(j,k)
                 a22 = sum00(j,k)-2.0d0*sum02(j,k)+sum22(j,k)
                 a12 = sum00(j,k)-sum01(j,k)-sum02(j,k)+sum12(j,k)
                 c1 = (a01*a22-a02*a12)/(a11*a22-a12**2)
                 c2 = (a02*a11-a01*a12)/(a11*a22-a12**2)
                 do i = 0,n
                    sSRng = (1.0d0-c1-c2)*s0(i,j,k)+c1*s1(i,j,k)+c2*s2(i,j,k)
                    sSR(i,j,k) = sSR(i,j,k)+blend*(sSRng-sSR(i,j,k))
                    sSR(i,k,j) = sSR(i,j,k)
                 end do
              else
                 write (3,*) '*** Too many iterations in BpySR. Quitting.'
                 write (*,*) '*** Too many iterations in BpySR. Quitting.'
                 stop
              end if
           end do
        end do
        !     write (3,90) iter, rms, sSR(0,1,1), sSR(0,1,2), sSR(0,2,2)
        write (*,90) iter, rms, (sSR(0,i,i),i=1,nsp)
90      format (10x, i6, 1pe13.2, 0pf12.4, 5f12.4)
        !     End Ng acceleration.

        !     Check for convergence.
        !       check for possible extrapolation attempt.
        if (iter .eq. iext) then
           iext = iext+next
           do j = 1,nsp
              do k = j,nsp
                 do i = 0,n
                    s2(i,j,k) = s0(i,j,k)-(s0(i,j,k)-s1(i,j,k))**2 &
                         /(s0(i,j,k)-2.0d0*s1(i,j,k) +s2(i,j,k))
                    if (r(i)*dabs(s2(i,j,k)-sSR(i,j,k)) .gt. 2.0d0) cycle
                 end do
                 do i = 0,n
                    sSR(i,j,k) = s2(i,j,k)
                 end do
              end do
           end do
        end if
        !       End extrapolation attempt.
     Enddo
     !     End Picard iteration.

     !     Done! Calculate g(r;j,k) and BpySR0(r;j,k). Print heading.
     do j = 1,nsp
        do k = 1,nsp
           do i = 1,n
              g(i,j,k) = eBond(i,j,k)*(1.0d0+sSR(i,j,k))
              BpySR0(i,j,k) = log(1.0d0+sSR(i,j,k))-sSR(i,j,k)
           end do
        end do
     end do
     write (3,103) iter, Gamma, rho(1), rho(2)
     write (*,103) iter, Gamma, rho(1), rho(2)
103  format (/3x, 'SUBROUTINE BpySR ' &
          /3x, 'Thermodynamics of a "short-ranged 2D electrolyte" using the' &
          ' PY equation' &
          /3x, 'iter =', i3, ': Gamma =', f8.4, ', rho1 =', f8.4, ', rho2 =', f8.4)

     !     Calculate PY thermodynamics of SR potentials.
     sumP1 = 0.0d0
     sumP2 = 0.0d0
     sumU = 0.0d0
     sumX = 0.0d0
     do j = 1,nsp
        do k = 1,nsp
           sumP0 = 0.0d0
           sumU0 = 0.0d0
           sumX0 = 0.0d0
           g(Ncore(j,k),j,k) = 0.5d0*g(Ncore(j,k),j,k)
           do i = 1,n
              sumP0 = sumP0+dr(i)*r(i)*g(i,j,k)*dexp(-r(i)**2)
              sumU0 = sumU0+dr(i)*r(i)*g(i,j,k)*phiSR0(i,j,k)
              sumX0 = sumX0+dr(i)*r(i)*(g(i,j,k)-1.0d0)
           end do
           g(Ncore(j,k),j,k) = 2.0d0*g(Ncore(j,k),j,k)
           sumP1 = sumP1+rho(j)*rho(k)*sigma(j,k)**2*g(Ncore(j,k),j,k)
           if (j.eq.k) Then
              sumP2 = sumP2+Gamma*rho(j)*rho(k)*z(j)*z(k)*sumP0
           else
              sumP2 = sumP2+Gamma*rho(j)*rho(k)*lambda*z(j)*z(k)*sumP0
           Endif

           sumU = sumU+rho(j)*rho(k)*sumU0
           sumX = sumX+rho(j)*rho(k)*sumX0
        end do
     end do
     P1 = 1.0d0+(pi/(2.0d0*rhoTotal))*sumP1
     P2 = (pi/(2.0d0*rhoTotal))*sumP2
     U = (pi/rhoTotal)*sumU
     X = 1.0d0+(2.0d0*pi/rhoTotal)*sumX
     write (3,110) P1+P2, U, X, P1, P2
     write (*,110) P1+P2, U, X, P1, P2
110  format (8x, 'pA/NkT =', f8.4, ', U/NkT =', f8.4, ', NkTX/A =', f8.4/ &
          15x, '=', f8.4, '  (HD)'/15x, '=', f8.4, '  (QQ-SR)')
   end subroutine BpySR



   FUNCTION Chebev(a,b,c,m,x)
     !     Taken from "Numerical Recipes" by W.H. Press et al.
     implicit none
     integer :: m, j
     real(kind=8) :: c(m), x, a, b, Chebev
     real(kind=8) :: d,dd, sv, y2, y
     if ((x-a)*(x-b) .gt. 0.0d0) then
        write (*,*) '*** x not in range in function Chebev'
        stop
     end if
     d = 0.0d0
     dd = 0.0d0
     y = (2.0d0*x-a-b)/(b-a)
     y2 = 2.0d0*y
     do j = m,2,-1
        sv = d
        d = y2*d-dd+c(j)
        dd = sv
     end do
     Chebev = y*d-dd+0.5d0*c(1)
   end function Chebev



   FUNCTION E1(x)
     !     Returns the exponential integral E_1(x) for positive real x.
     !     See AMS55, p.231.
     implicit none
     real(kind=8) :: x, E1, y
     real(kind=8), save :: p0,p1,p2,p3,p4,p5
     data p0,p1,p2,p3,p4,p5/-0.57721566d0, 0.99999193d0,-0.24991055d0, &
          0.05519968d0,-0.00976004d0, 0.00107857d0/
     real(kind=8), save :: q1,q2,q3,q4
     data q1,q2,q3,q4/8.5733287401d0,18.0590169730d0, &
          8.6347608925d0, 0.2677737343d0/
     real(kind=8), save :: r1,r2,r3,r4
     data r1,r2,r3,r4/ 9.5733223454d0,25.6329561486d0, &
          21.0996530827d0, 3.9584969228d0/
     if (x .lt. 1.0d0) then
        E1 = -log(x)+(p0+x*(p1+x*(p2+x*(p3+x*(p4+x*p5)))))
     else
        y = 1.0d0/x
        E1 = (1.0d0+y*(q1+y*(q2+y*(q3+y*q4))))/(1.0d0+y*(r1+y*(r2+y*(r3+y*r4)))) &
             *dexp(-x)/x
     end if
   end function E1



   SUBROUTINE Hankel(n,x,dx,Cin,Fin,Cout,Fout)
     use datatrans, only : mxNr, W, BJ1sq, root
     implicit none
     real(kind=8) ::  x(0:mxNr), dx(0:mxNr), Fin(0:mxNr), Fout(0:mxNr)
     real(kind=8) :: Cout, Cin, sumFt
     integer :: n, i, k
     sumFT = 0.0d0
     do i = 1,n
        sumFT = sumFT+dx(i)*x(i)*Fin(i)
     end do
     Fout(0) = (Cout/Cin)*sumFT
     do k = 1,n
        sumFT = 0.0d0
        do i = 1,n
           sumFT = sumFT+W(k,i)*Fin(i)
        end do
        Fout(k) = sumFT
     end do
   end subroutine Hankel



   FUNCTION SerAdd(x,n,a)
     !     SerAdd = a(1)+x(a(2)+x(a(3)+x( ...+x(a(n))...)))
     implicit none
     integer :: n
     real(kind=8) :: a(n), x, SerAdd
     integer :: i
     SerAdd = a(n)
     do i = n-1,1,-1
        SerAdd = a(i)+x*SerAdd
     end do
   end function SerAdd



   SUBROUTINE setupW(Nr)
     !     Uses SerAdd, BessJY
     use datatrans, only : mxNr, W, BJ1sq, root, pi 
     implicit none
     integer :: i, j, k, Nr 
     real (kind=8), parameter :: xerrMx=1.0d-12
     real(kind=8), dimension(4) :: croot = (/ 1.0d0, -41.33333333333333d0,&
          & 8061.86666666667d0, -3826125.40952d0 /)  
     real(kind=8), dimension(20) :: root1 = (/ 2.4048255577d0,&
          & 5.5200781103d0, 8.6537279129d0,11.7915344391d0,14.9309177086d0,18.0710639679d0,&
          21.2116366299d0,24.3524715308d0,27.4934791320d0&
          &,30.6346064684d0, 33.7758202136d0,36.9170983537d0&
          &,40.0584257646d0,43.1997917132d0, 46.3411883717d0&
          &,49.4826098974d0,52.6240518411d0,55.7655107550d0,&
          & 58.9069839261d0,62.0484691902d0 /) 
     real (kind=8) :: x, BessJ,BessY,BessJp,BessYp,BesJ0,BesY0,BesJ0p&
          &,BesY0p, xerror, xnew, a1, a2, b, xsq
     real (kind=8), external :: SerAdd
     
     !     Calculate roots of zero-order Bessel function J0(x).
     !     First pass. Use explicit roots or construct polynomial approximations.
     do i = 1,20
        root(i) = root1(i)
     end do
     do i = 21,Nr
        b = (dble(i)-0.25d0)*pi
        x = 1.0d0/(8.0d0*b)
        xsq = x**2
        root(i) = b+x*seradd(xsq,4,croot)
     end do
     !     Second pass. Refine first-pass roots with Newton's method.
     do i = 1,Nr
        x = root(i)
        k = 0
1       continue
        k = k+1
        call BessJY(x,0.0d0,BesJ0,BesY0,BesJ0p,BesY0p)
        xnew = x-BesJ0/BesJ0p
        xerror = dabs((xnew-x)/xnew)
        x = xnew
        if (k .gt. 200) then
           write (*,*) '*** Count exceeds 200 in subroutine SETUPW'
           stop
        end if
        if (xerror .gt. xerrMx) go to 1
        root(i) = x
        BJ1sq(i) = BesJ0p**2
     end do

     !     Set up W table.
     a1 = 2.0d0/root(Nr)
     do i = 1,Nr-1
        a2 = root(i)/root(Nr)
        do j = 1,Nr-1
           x = root(j)*a2
           call BessJY(x,0.0d0,BessJ,BessY,BessJp,BessYp)
           W(i,j) = a1*BessJ/BJ1sq(j)
        end do
     end do
   end subroutine setupW

   subroutine lrfuncs(flr,tflr,r,q,Nr,a)
     implicit none
     real(kind=8) :: flr(0:Nr), tflr(0:Nr), r(0:Nr), q(0:Nr)
     real(kind=8), parameter :: pi=3.141592653589793d0 
     real(kind=8) :: i0lr, i1lr, q2a, q2, I0, I1, a
     real(kind=8), external :: dbsI0E, dbsI1E
     integer :: i, Nr
     flr(0) = 0
     tflr(0) = pi*sqrt(a*pi)
     do i=1, Nr
        flr(i) = (1-(1+a*r(i)**2)*exp(-a*r(i)**2))/r(i)**3
        q2 = q(i)**2
        q2a= q2/a
        I0 = dbsI0E(q2a/8)
        I1 = dbsi1E(q2a/8)
        tflr(i) = 2*pi*(-Q(i) + sqrt(pi/a)*((2*a+q2)&
             &*I0+q2*I1)/4)
!        write(99,'(4f15.10)')q(i),tflr(i),r(i),flr(i)
     enddo
end subroutine lrfuncs
function mk(k,R)
  use datatrans
  implicit none
  real(kind=8), external :: dbsj1
  real(kind=8) :: k, r, mk
  if (k*R > 1.0d-6) then
     mk = 2*pi*r*dbsj1(k*r)/k
  else
     mk = pi*r**2
  endif
end function mk
