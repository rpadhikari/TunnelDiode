subroutine dealloc
  use params, only : TransG, T12_mat_arr, T13_mat_arr, T23_mat_arr, &
                     II, TM, IPIV, UB, T, VV, &
                     U1, U1_basis, E, f1, f2, diagU1, &
                     sig1, sig2, sig3, gam1, gam2, gam3, &
                     Gbinv, G, WORKK, gaG1, gaG2, gaGT2, gaGT3, &
                     T12_mat, T13_mat, T23_mat
  deallocate(TransG)
  deallocate(T12_mat_arr, T13_mat_arr, T23_mat_arr)
  deallocate(II)
  deallocate(TM)
  deallocate(IPIV)
  deallocate(UB, T, VV)
  deallocate(U1, U1_basis, E, f1, f2, diagU1)
  deallocate(sig1, sig2, sig3)
  deallocate(gam1, gam2, gam3)
  deallocate(Gbinv, G, WORKK)
  deallocate(gaG1, gaG2, gaGT2, gaGT3)
  deallocate(T12_mat, T13_mat, T23_mat)
end subroutine dealloc
