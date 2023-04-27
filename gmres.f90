module gmres_solver

  use common
  implicit none
  
  public :: arnoldi, apply_givens_rotation, gmres
  
  contains
  
  subroutine arnoldi(A_in,Q_in,k_in, OUT)
    
    real(dp), dimension(:,:), intent(out) :: OUT(0:,0:)
    real(dp), dimension(:,:), intent(in)  :: A_in, Q_in
    integer, intent(in)                   :: k_in
    real(dp), dimension(:), allocatable   :: q, h
    integer                               :: i
    integer, dimension(0:1)               :: Qsize
    
    Qsize = shape(Q_in)
    
    allocate(q(0:Qsize(0)-1))
    allocate(h(0:Qsize(0)-1))
    
    ! Krylov vector
    q = matmul(A_in, Q_in(:,k_in))
    
    ! Modified Gram-Schmidt
    do i = 0, k_in-1
      
      h(i) = matmul(q,Q_in(:,i))
      q    = q - matmul(h(i),Q_in(:,i))
      
    end do
    
    h(k_in) = norm2(q)
    q = q/norm2(q)
    
    OUT(:,0) = h
    OUT(:,1) = q
    
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
    
    real(dp), dimension(:), intent(out) :: h_io(0:)
    real(dp), intent(out)               :: cs_k,  sn_k
    real(dp), dimension(:), intent(in)  :: cs_in(0:), sn_in(0:)
    integer, intent(in)                 :: k_in
    real(dp)                            :: temp
    integer                             :: i
    
    !applied to i-th column
    
    do i = 0, k_in-2
      
      temp      =  cs_in(i)*h_io(i) + sn_in(i)*h_io(i + 1)
      h_io(i+1) = -sn_in(i)*h_io(i) + cs_in(i)*h_io(i + 1)
      h_io(i)   = temp
      
    end do
    
    ! update the next sin cos values for rotation
    call givens_rotation(h_io(k_in-1),h_io(k_in), cs_k, sn_k)
    
    ! eliminate H(i,i-1)
    h_io(k_in-1) = cs_k*h_io(k_in) + sn_k*h_io(k_in)
    h_io(k_in)   = 0.0_dp
    
  end subroutine
  
  subroutine gmres(A_in, b_in, x_in, Iter_in, threshold_in, out, e)
    
    real(dp), dimension(:), intent(out)   :: out, e
    real(dp), dimension(:), intent(in)    :: b_in,x_in
    real(dp), dimension(:,:), intent(in)  :: A_in
    integer, intent(in)                   :: Iter_in
    real(dp), intent(in)                  :: threshold_in
    real(dp), dimension(:), allocatable   :: r, cs, sn, e1, beta, error, y
    integer                               :: n, i, j, k
    integer, dimension(0:1)               :: in_shape
    real(dp)                              :: b_norm, r_norm
    real(dp), dimension(:,:), allocatable :: Q, H, temp_var
    
    in_shape = shape(A_in)
    n = max(in_shape(0), in_shape(1))
    
    allocate(r(0:size(x_in)-1))
    allocate(cs(0:Iter_in-1))
    allocate(sn(0:Iter_in-1))
    allocate(error(0:Iter_in-1))
    allocate(e1(0:Iter_in))
    allocate(beta(0:Iter_in))
    allocate(temp_var(0:Iter_in,0:1))
    allocate(Q(0:size(x_in)-1,0:Iter_in-1))
    allocate(H(0:size(x_in)-1,0:Iter_in-1))
    
    cs = 0.0_dp
    e1 = 0.0_dp
    
    e1(0) = 1.0_dp
    
    r = b_in - matmul(A_in,x_in)
    
    b_norm = norm2(b_in)
    
    error = norm2(r) / b_norm
    
    r_norm = norm2(r)
    Q(:,0) = r/r_norm
    beta   = r_norm * e1
    
    do k = 0, Iter_in-1
      
      call arnoldi(A_in, Q, k, temp_var)
      
      H(:,k)   = temp_var(:,0)
      Q(:,k+1) = temp_var(:,1)
      
      
      ! eliminate the last element in H ith row and update the rotation matrix
      call apply_givens_rotation(H(:,k), cs, sn, k , cs(k), sn(k))
      
      
      ! Residual vector update
      beta(k+1) = -sn(k)*beta(k)
      beta(k)   =  cs(k)*beta(k)
      error(k)  = abs(beta(k+1))/b_norm
      
      if (error <= threshold_in) exit
      
    end do
    
    e = error
    
    out = out + matmul(Q(:,0:Iter_in-1), matmul(inv(H(0:Iter_in-1,0:Iter_in-1)),beta(0:Iter_in-1)))
    
    deallocate(r)
    deallocate(cs)
    deallocate(sn)
    deallocate(Q)
    
    end subroutine
    
    
! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
! found on https://fortranwiki.org/fortran/show/Matrix+inversion
function inv(A) result(Ainv)
  real(dp), dimension(:,:), intent(in) :: A
  real(dp), dimension(size(A,1),size(A,2)) :: Ainv

  real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end function inv
  
end module
