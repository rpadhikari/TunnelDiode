  implicit none
  complex(8),parameter :: j=(0,1)

  real(8),parameter :: pi=3.14159265358979 ! Constant Pi
  real(8),parameter :: hbar=1.05457180013d-34 ! Joule -Sec
  real(8),parameter :: q=1.602176620898d-19 ! Coulomb
  real(8),parameter :: m=0.250d0*9.109383560d-31  ! Kg
  real(8),parameter :: ie=3.87404586591771d-05 ! (q*q)/(2*pi*hbar);
  real(8),parameter :: ef=0.10d0 ! eV
  real(8),parameter :: Kt=0.0250d0 ! eV
!inputs
  real(8),parameter :: a=3.0e-10
  real(8),parameter :: t0=1.69332536885689d0 !(hbar^2)/(2*m*(a^2)*q)
  real(8) dx, V, mu1, mu2, dx_U1, dE, It

  real(8) T12,T13,T23

  integer(8) INFO, r
  integer(4),parameter :: NS=15
  integer(4),parameter :: NC=16
  integer(4),parameter :: ND=15
  integer(4),parameter :: Np=NS+NC+ND, NV=52, NE=201
  integer(4),parameter :: LDA=Np,LWORK=Np
  integer(4) UB_shape, i, k, iv
  complex(8) zplus, ck, ka
  real(8),allocatable :: TM(:)
  real(8),allocatable :: II(:)
  integer(8),allocatable :: IPIV(:)
  real(8), allocatable :: UB(:),T(:,:),VV(:)
  real(8),allocatable :: U1(:), U1_basis(:), E(:),f1(:),f2(:),diagU1(:,:)
  complex(8),allocatable :: sig1(:,:),sig2(:,:),sig3(:,:)
  complex(8),allocatable :: gam1(:,:),gam2(:,:),gam3(:,:)
  complex(8),allocatable :: Gbinv(:,:), G(:,:),ceye(:,:),WORKK(:)
  complex(8),allocatable :: gaG1(:,:), gaG2(:,:), gaGT2(:,:), gaGT3(:,:)
  complex(8),allocatable :: T12_mat(:,:),T13_mat(:,:),T23_mat(:,:),TransG(:,:)
  real(8),allocatable :: T12_mat_arr(:),T13_mat_arr(:),T23_mat_arr(:)

end module params
