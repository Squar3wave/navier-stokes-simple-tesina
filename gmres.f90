module gmres_solver

  use common
  implicit none
  public :: arnoldi, apply_givens_rotation
  
  subroutine arnoldi(A_in,Q_in,k_in, OUT)
    
    real(dp), dimension(:,:), intent(out) :: OUT
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
      
      h(i) = matmul(conjg(q),Q_in(:,i))
      q    = q - matmul(h(i),Q_in(:,i))
      
    end do
    
    h(k_in) = norm2(q)
    q = q/norm2(q)
     
    OUT(:,0) = q
    OUT(:,1) = h
     
    deallocate(q)
    deallocate(h)
    
  end subroutine

  subroutine givens_rotation(v1,v2,cs_out,sn_out)
  
    real(dp)                            :: t
    real(dp), dimension(:), intent(out) :: cs_k,  sn_k
    real(dp), intent(in)                :: v1,v2
    
    t=sqrt(v1^2 + v2^2)
    
    cs_out = v1/t
    sn_out = v2/t
    
  end subroutine
  
  subroutine apply_givens_rotation(h_io, cs_in, sn_in, k_in, cs_k, sn_k)
    
    real(dp), dimension(:), intent(out) :: cs_k,  sn_k, h_io
    real(dp), dimension(:), intent(in)  :: cs_in, sn_in
    integer, intent(in)                 :: k_in
    real(dp)                            :: temp
    
    !applied to i-th column
    
    do i = 0, k_in-2
      
      temp      =  cs_in(i)*h_io(i) + sn_in(i)*h_io(i + 1)
      h_io(i+1) = -sn_in(i)*h_io(i) + cs_in(i)*h_io(i + 1)
      h_io(i)   = temp
      
    end do
    
    ! update the next sin cos values for rotation
    call givens_rotation(h(k_in-1),h(k_in), cs_k, sn_k)
    
    ! eliminate H(i,i-1)
    h_io(k_in-1) = cs_k*h_io(k_in) + sn_k*h(k_in)
    h_io(k_in)   = 0.0_dp
    
  end subroutine
  
  subroutine gmres(A_in, b_in, x_in, Iter_in, threshold_in, out, e)
    
    
    
  end subroutine
  
end module
