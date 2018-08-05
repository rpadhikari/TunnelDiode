! allocate the variarrays

subroutine alloc
  use params, only : Np, NV, NE, NS, ND, LWORK, &
                     T12_mat_arr, T13_mat_arr, T23_mat_arr,&
                     II, TM, IPIV, &
                     T, VV, U1, U1_basis, E, f1, f2, diagU1, &
                     sig1, sig2, sig3, gam1, gam2, gam3, &
                     Gbinv, G, WORKK, ceye, gaG1, gaG2, gaGT2, gaGT3, &
                     T12_mat, T13_mat, T23_mat, TransG
  allocate(T12_mat_arr(Np),T13_mat_arr(Np),T23_mat_arr(Np))
  allocate(II(NV))
  allocate(TM(NE))
  allocate(IPIV(Np))
  allocate(T(Np,Np),VV(NV))
  allocate(U1(Np), U1_basis(Np),E(NE),f1(NE),f2(NE),diagU1(Np,Np))
  allocate(sig1(Np,Np),sig2(Np,Np),sig3(Np,Np))
  allocate(gam1(Np,Np),gam2(Np,Np),gam3(Np,Np))
  allocate(Gbinv(Np,Np),G(Np,Np),WORKK(LWORK))
  allocate(ceye(Np,Np))
  allocate(gaG1(Np,Np),gaG2(Np,Np),gaGT2(Np,Np),gaGT3(Np,Np))
  allocate(T12_mat(Np,Np),T13_mat(Np,Np),T23_mat(Np,Np))
  allocate(TransG(Np,Np))
end subroutine alloc
