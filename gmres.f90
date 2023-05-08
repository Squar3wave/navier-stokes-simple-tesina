module gmres_solver

  use matdef
  use common
  use sparsealg
  implicit none
  
  public :: arnoldi, apply_givens_rotation, gmres
  
  contains
  
  subroutine arnoldi(A_in,Q_in,k_in, OUT)
    
    real(dp), dimension(:,:), intent(out)  :: OUT
    !real(dp), dimension(:,:), intent(in)  :: A_in, Q_in
    real(dp), dimension(:,:), intent(in)   :: Q_in
    type(rCSR), intent(in)                 :: A_in
    integer, intent(in)                    :: k_in
    real(dp), dimension(:), allocatable    :: q, h
    integer                                :: i
    integer, dimension(0:1)                :: Qsize
    
    Qsize = shape(Q_in)
    
    allocate(q(Qsize(0)))
    allocate(h(Qsize(0)))
    
    ! Krylov vector
    !q = matmul(A_in, Q_in(:,k_in))
    
    call matvec(A_in, Q_in(:,k_in), q)
    
    ! Modified Gram-Schmidt, keeping the Hessenberg matrix
    do i = 1, k_in
      
      h(i) = dot_product(q,Q_in(:,i))
      q    = q - h(i)*Q_in(:,i)
    end do
    
    h(k_in) = norm2(q)
    q = q/norm2(q)
    
    OUT(:,1) = h
    OUT(:,2) = q
    
    deallocate(q)
    deallocate(h)
    
  end subroutine

  subroutine givens_rotation(v1,v2,cs_out,sn_out)
  
    real(dp)              :: t
    real(dp), intent(out) :: cs_out,  sn_out
    real(dp), intent(in)  :: v1,v2
    
    t=sqrt(v1**2 + v2**2)
    
    cs_out = v1/t
    sn_out = v2/t
    
  end subroutine
  
  subroutine apply_givens_rotation(h_io, cs_in, sn_in, k_in, cs_k, sn_k)
    
    real(dp), dimension(:), intent(out) :: h_io
    real(dp), intent(out)               :: cs_k,  sn_k
    real(dp), dimension(:), intent(in)  :: cs_in, sn_in
    integer, intent(in)                 :: k_in
    real(dp)                            :: temp
    integer                             :: i
    
    !applied to i-th column
    
    do i = 1, k_in-1
      
      temp      =  cs_in(i)*h_io(i) + sn_in(i)*h_io(i + 1)
      h_io(i+1) = -sn_in(i)*h_io(i) + cs_in(i)*h_io(i + 1)
      h_io(i)   = temp
      
    end do
    
    ! update the next sin cos values for rotation
    call givens_rotation(h_io(k_in-1),h_io(k_in), cs_k, sn_k)
    
    ! eliminate H(i,i-1)
    h_io(k_in) = cs_k*h_io(k_in) + sn_k*h_io(k_in+1)
    h_io(k_in+1)   = 0.0_dp
    
  end subroutine
  
  subroutine gmres(A_in, b_in, x_out, e, Iter_in, threshold_in)
    
    real(dp), dimension(:), intent(inout)  :: x_out, e
    real(dp), dimension(:), intent(in)     :: b_in
    type(rCSR), intent(in)                 :: A_in
    integer, intent(in)                    :: Iter_in
    real(dp), intent(in)                   :: threshold_in
    real(dp), dimension(:), allocatable    :: r, cs, sn, e1, beta, error, y, tmp_v
    integer                                :: n, k, i, IPIV, LDB, INFO, xsize
    real(dp)                               :: b_norm, r_norm
    real(dp), dimension(:,:), allocatable  :: mQ, mH, temp_var, tmp_m
    integer, dimension(2)                  :: tmp_m_shape
    
    xsize = size(x_out)
    
    allocate(       r(xsize))
    allocate(   tmp_v(xsize))
    allocate(      cs(Iter_in))
    allocate(      sn(Iter_in))
    allocate(   error(Iter_in))
    allocate(       y(Iter_in))
    allocate(      e1(Iter_in+1))
    allocate(    beta(Iter_in+1))
    allocate(temp_var(Iter_in+1,2))
    allocate(   tmp_m(Iter_in,Iter_in))
    allocate(      mQ(xsize,Iter_in))
    allocate(      mH(xsize,Iter_in))
    
    cs = 0.0_dp
    e1 = 0.0_dp
    
    e1(1) = 1.0_dp
    
    call matvec(A_in,x_out,tmp_v)
    
    r = b_in - tmp_v
    
    b_norm = norm2(b_in)
    
    error = norm2(r) / b_norm
    
    r_norm = norm2(r)
    mQ(:,0) = r/r_norm
    beta   = r_norm * e1
    
    mH = 0.0_dp
    mQ = 0.0_dp
    
    print*, "start iterations and error minimization"
    
    do k = 1, Iter_in
      
      call arnoldi(A_in, mQ, k, temp_var)
      
      mH(1:k+1,k)   = temp_var(:,1)
      mQ(    :,k+1) = temp_var(:,2)
      
      
      ! eliminate the last element in H ith row and update the rotation matrix
      call apply_givens_rotation(mH(1:k+1,k), cs, sn, k , cs(k), sn(k))
      
      
      ! Residual vector update
      beta(k+1) = -sn(k)*beta(k)
      beta(k)   =  cs(k)*beta(k)
      error(k)  = abs(beta(k+1))/b_norm
      
      if (error(k) <= threshold_in) then
        exit
      end if
      
    end do
    
    e = error
    
    print *, "done"
    print *,""

    print *, "solving system with DGBVS"
    
    tmp_m = mH(1:Iter_in,1:Iter_in)
    !call dgbsv(size(tmp_m,1), size(tmp_m,1)-1, size(tmp_m,1)-1, 1,tmp_m, size(tmp_m,1), IPIV, y, size(y), INFO)
    
    tmp_m_shape = shape(tmp_m)
    
    call dgbsv(Iter_in, &        ! N: The number of linear equations
               1, &              ! KL: The number of subdiagonals
               1, &              ! KU: The number of superdiagonals
               1, &              ! NRHS: The number of column for B
               tmp_m, &          ! AB: The matrix A in band storage
               tmp_m_shape(1), & ! LDAB: AB's leading dimension
               IPIV, &           ! IPIV: The pivot indices that define the permutation matrix P
               y, &              ! B: out vector
               size(y), &        ! LDB: The leading dimension of the array B
               INFO)             ! INFO: if 0 all good, if <0 failed, if > 0 no solution
    
    print *, "done"
    print *, ""

    print *, "outputting result"
    
    
    !$OMP PARALLEL, public(Q,y,x)
    !$OMP DO
    
    do i=1,size(x_out)
      
      !print *, dot_product(Q(i,1:Iter_in),y)
      
      x_out(i)= x_out(i) + dot_product(mQ(i,1:Iter_in),y)
      
    end do
    
    !$OMP END DO
    !$OMP END PARALLEL
    print *, "done!"
    print *, ""
    !out = out + matmul(Q(:,0:Iter_in-1), &
    !matmul(inv(H(0:Iter_in-1,0:Iter_in-1)),beta(0:Iter_in-1)))
    
    print *, "deallocating arrays"
    
    deallocate(r)
    deallocate(tmp_v)
    deallocate(cs)
    deallocate(sn)
    deallocate(error)
    deallocate(y)
    deallocate(e1)
    deallocate(beta)
    deallocate(temp_var)
    deallocate(tmp_m)
    deallocate(mQ)
    deallocate(mH)
    
    print *, "exiting subroutine"
    print *,""
    end subroutine
 
  
end module
