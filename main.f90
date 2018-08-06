program main
  use params
  call alloc
  UB_shape=NS+4+(NC-8)+4+ND
! allocating shape for UB
  allocate(UB(UB_shape))

! creating UB array of size NS + 4 + (NC-8) + 4 + ND
  UB(1:NS+1)=0.0d0
  UB(NS+1:NS+1+4)=0.40d0
  UB(NS+1+4:NS+1+4+NC-8)=0.0d0
  UB(NS+1+4+NC-8:NS+1+4+NC-8+4)=0.40d0
  UB(NS+1+4+NC-8+4:)=0.0d0
  open(1,file='pot.dat', action='write')
  do i=1,UB_shape
    write(1,*) dble(i)*0.01, UB(i)
  end do
  close(1)
! setting up T matrix
  T=0.0d0
  T(1,1)=2.0d0*t0
  do i=2,Np
    T(i,i)=2.0d0*t0+UB(i)  ! I did T=T+diag(UB) in here since diag(UB) is just
    T(i,i-1)=-t0       ! the main diagonal matrix
    T(i-1,i)=-t0
  end do

! Bias
! Defining VV matrix
  dx=0.50d0/(NV-1)
  do i=1,NV
    VV(i)=0.0d0+dble(i-1)*dx
  end do

! Creating U1_basis vector for U1 later required
  U1_basis(1:NS+1)=0.50d0
  dx_U1=(-0.50d0-0.50d0)/(NC-1)
  do i=1,NC
    U1_basis(NS+i)=0.50d0+(i-1)*dx_U1
  end do
  U1_basis(NS+NC+1:)=-0.50d0

! linspacing E as vector
  do i=1, NE
    E(i)=-0.20d0+dble(i-1)*(0.80d0+0.20d0)/dble(NE-1)
  end do

  do i=1, Np
    ceye(i,i)=(1.0d0, 0.0d0) ! pure real number
  end do

  do iv=1, NV
    diagU1=0.0d0
    V=VV(iv)
    mu1=Ef+(V/2.0d0)
    mu2=Ef-(V/2.0d0)
    U1=V*U1_basis
    do i=1, Np
      diagU1(i,i)=U1(i)
    end do
! Energy grid for Green's function method
    zplus=(0,1.0d-12)
    dE=E(2)-E(1)
    f1=1.0d0/(1+exp((E-mu1)/Kt))
    f2=1.0d0/(1+exp((E-mu2)/Kt))
! For infinite 2D cross-section
! f1=(2*m*Kt*q)/(2*pi*hbar**2)*log(1+exp((mu1-E)/Kt))
! f2=(2*m*Kt*q)/(2*pi*hbar**2)*log(1+exp((mu2-E)/Kt))

! Transmission

    It=0.0d0
    do k=1, NE
      sig1=(0.0d0,0.0d0)
      sig2=(0.0d0,0.0d0)
      sig3=(0.0d0,0.0d0)
      gam1=sig1
      gam2=sig2
      gam3=sig3
      ck=1.0d0-((E(k)+zplus-U1(1)-UB(1))/(2*t0))
      ka=acos(ck)
      sig1(1,1)=-t0*exp(j*ka)
      gam1(1,1)=j*2.0d0*imag(sig1(1,1))*j ! gam1=i*(sig1-sig1')
      ck=1.0d0-((E(k)+zplus-U1(Np)-UB(Np))/(2.0d0*t0))
      ka=acos(ck)
      sig2(Np,Np)=-t0*exp(j*ka)
      gam2(Np,Np)=j*2.0d0*imag(sig2(Np,Np))*j ! gam2=i*(sig2-sig2')
      sig3(Np/2,Np/2)=-1.0d0*j*0.000250d0
      gam3(Np/2,Np/2)=j*2.0d0*imag(sig3(Np/2,Np/2))*j
! gam3=j*(sig3-sig3') no need to do as (Buttiker probe)
! sig3=0 matrix by default and sig3' as well
      Gbinv=(E(k)+zplus)*ceye-T-diagU1-sig1-sig2-sig3
      call ZGETRF(Np,Np,Gbinv,LDA,IPIV,INFO)
      call ZGETRI(Np,Gbinv,LDA,IPIV,WORKK,LWORK,INFO)
      G=Gbinv
      TransG=transpose(G)
      TransG=conjg(G)
      gaG1=matmul(gam1,G)
      gaGT2=matmul(gam2,TransG)
      gaGT3=matmul(gam3,TransG)
      gaG2=matmul(gam2,G)
      T12_mat=matmul(gaG1,gaGT2)
      T13_mat=matmul(gaG1,gaGT3)
      T23_mat=matmul(gaG2,gaGT3)
      do i=1,Np
        T12_mat_arr(i)=real(T12_mat(i,i))
        T13_mat_arr(i)=real(T13_mat(i,i))
        T23_mat_arr(i)=real(T23_mat(i,i))
      end do
      T12=sum(T12_mat_arr)
      T13=sum(T13_mat_arr)
      T23=sum(T23_mat_arr)
      TM(k)=T12+(T13*T23/(T12+T23))
      It=It+(dE*IE*TM(k)*(f1(k)-f2(k)))
    end do
    II(iv)=It
  end do

  open(2,file='iv_tunnel.dat',action='write')
  do i=1, NV
    write(2,*) VV(i), II(i)
  end do
  close(2)
  open(3,file='trans.dat', action='write')
  do i=i, NE
    write(3,*) TM(i), E(i)
  end do
  close(3)
  call dealloc

end program main
