module gmres_solver

  use matdef
  use common
  use sparsealg
  implicit none
  
  public :: gmres, gmres_preliminary
  
  type(rCSR), public                              :: A_csr
  real(dp), dimension(:), allocatable, public     :: x_gmres, b_gmres, error 
  integer, public                                 :: max_iter, vecsize, Nrow, flag=0
  real(dp), public                                :: threshold
  real(dp), dimension(:), allocatable, public     :: r, cs, sn, e1, beta_gmres, error_gmres, y_gmres, tmp_v
  real(dp), dimension(:,:), allocatable, public   :: mQ, mH, temp_var, tmp_m
  
  contains

  subroutine allocate_gmres_stuff()
    
    allocate(          r(vecsize))
    allocate(      tmp_v(vecsize))
    allocate(         cs(max_iter))
    allocate(         sn(max_iter))
    allocate(error_gmres(max_iter))
    allocate(    y_gmres(max_iter))
    allocate(         e1(max_iter+1))
    allocate( beta_gmres(max_iter+1))
    allocate(   temp_var(max_iter+1,2))
    allocate(      tmp_m(max_iter,max_iter))
    allocate(         mQ(vecsize,max_iter))
    allocate(         mH(vecsize,max_iter))
    
    flag=1
    
  end subroutine
  
  subroutine apply_constraint_gmres(nf)
    
    integer :: i, j
    integer, intent(in) :: nf 
    
    Sp(1:NXmax+1,:,nf)       = F(1:NXmax+1,:,nf)
    Sp(1:NXmax+1,NYmax+1,nf) = F(1:NXmax+1,NYmax+1,nf)
    Sp(:,NYmax+1,nf)         = F(:,NYmax+1,nf)
    Sp(NXmax+1,1:NYmax+1,nf) = F(NXmax+1,1:NYmax+1,nf)
    
  end subroutine

  subroutine gmres_preliminary()
    
    print *, "======================="
    print *, "preliminary procedures"
    print *, "======================="
    print *,""
    
    print *, "checking if CSR is already allocated"
    print *,""
    
    if(.not. allocated(A_csr%nzval)) then
      print *, "create CSR"
      call create(A_csr, Nrow, 5*Nrow)
      print *, "done"
      print *,""
    else
      print *,"A_csr already allocated"
      print *,""
    end if
    
    print *, "apply costraints"
    call apply_constraint_gmres(1)
    print *, "done"
    print *,""
    
    print *, "convert matrix to CRS"
    call convert2csr(Ap,Ae,An,As,Aw, A_csr)
    print *, "done"
    print *,""
    
    print *, "checking if vectors are already allocated"
    if(.not. allocated(x_gmres) .and. .not. allocated(b_gmres) .and. .not. allocated(error)) then
      print *, "allocating x"
      allocate(x_gmres(Nrow))
      print *, "allocating b"
      allocate(b_gmres(Nrow))
      print *, "allocating error"
      allocate(error(size(x_gmres)))
    else
      print *, "vectors already allocated"
    end if
    
    print *, "done"
    print *,""
    
    vecsize = size(x_gmres)
    
    print *, "checking if gmres auxiliary vectors are already allocated"
    if(flag==0) then
      print *,"allocating auxiliary GMRES vectors"
      call allocate_gmres_stuff()
    else
      print *, "gmres auxiliary vectors are already allocated"
      print *, ""
    end if
    print *, "done"
    print *,""
    
    print *, "reshape arrays"
    x_gmres(:) = reshape( F(:,:,1), shape=[size( F(:,:,1))])
    b_gmres(:) = reshape(Sp(:,:,1), shape=[size(Sp(:,:,1))])
    print *, "done"
    print *,""
    
    error= 0.0_dp
    
  end subroutine  
  
  subroutine gmres(A_in, b_in, x_out, e, Iter_in, threshold_in)
    
    real(dp), dimension(:), intent(out) :: x_out, e
    real(dp), dimension(:), intent(in)  :: b_in
    type(rCSR), intent(in)              :: A_in
    integer, intent(in)                 :: Iter_in
    real(dp), intent(in)                :: threshold_in
    integer                             :: n, k, i, IPIV, LDB, INFO, xsize
    real(dp)                            :: b_norm, r_norm
    integer, dimension(2)               :: tmp_m_shape
    
    xsize = size(x_out)
    
    
    print *, "======================="
    print *, "starting solver"
    print *, "======================="
    print *,""
    
    cs = 0.0_dp
    e1 = 0.0_dp
    
    e1(1) = 1.0_dp
    
    call matvec(A_in,x_out,tmp_v)
    
    r = b_in - tmp_v
    
    b_norm = norm2(b_in)
    
    error_gmres = norm2(r) / b_norm
    
    r_norm       = norm2(r)
    mQ(:,0)      = r/r_norm
    beta_gmres   = r_norm * e1
    
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
      beta_gmres(k+1) = -sn(k)*beta_gmres(k)
      beta_gmres(k)   =  cs(k)*beta_gmres(k)
      error_gmres(k)  = abs(beta_gmres(k+1))/b_norm
      
      if (error_gmres(k) <= threshold_in) then
        exit
      end if
      
    end do
    
    e = error_gmres
    
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
               y_gmres, &        ! B: out vector
               size(y_gmres), &  ! LDB: The leading dimension of the array B
               INFO)             ! INFO: if 0 all good, if <0 failed, if > 0 no solution
    
    print *, "done"
    print *, ""

    print *, "outputting result"
    
    
    !$OMP PARALLEL, public(mQ,y_gmres,x_out, Iter_in)
    !$OMP DO
    
    do i=1,size(x_out)
      
      !print *, dot_product(mQ(i,1:Iter_in),y_gmres)
      
      x_out(i)= x_out(i) + dot_product(mQ(i,1:Iter_in),y_gmres)
      
    end do
    
    !$OMP END DO
    !$OMP END PARALLEL
    print *, "done!"
    print *, ""
    
    
    
    !out = out + matmul(Q(:,0:Iter_in-1), &
    !matmul(inv(H(0:Iter_in-1,0:Iter_in-1)),beta(0:Iter_in-1)))
    
    
    print *, "exiting subroutine"
    print *,""
    
  end subroutine
  
  subroutine arnoldi(A_in,Q_in,k_in, OUT)
    
    real(dp), dimension(:,:), intent(out)  :: OUT
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
 

 
end module
