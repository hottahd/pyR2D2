!===================================================
subroutine d_x(qq, x, ix, jx, kx, qqd) bind(c)
!===================================================
  use iso_c_binding
  use omp_lib
  implicit none

  integer :: i, j, k
  real(8) :: a, b, c, d
  integer, parameter :: margin = 2
  integer(c_int), intent(in) :: ix, jx, kx
  real(c_double), dimension(ix, jx, kx), intent(in) :: qq
  real(c_double), dimension(ix), intent(in) :: x
  real(c_double), dimension(ix, jx, kx), intent(out) :: qqd

  real(4), dimension(ix) :: rk1, rk2, rk3, rk4

!---------------------------------------------------

  do i = 1 + margin, ix - margin
    a = x(i + 2) - x(i)
    b = x(i + 1) - x(i)
    c = x(i) - x(i - 1)
    d = x(i) - x(i - 2)
    rk1(i) = (c*d - b*(c + d))/(a - b)/(a + c)/(a + d)
    rk2(i) = (a*(c + d) - c*d)/(a - b)/(b + c)/(b + d)
    rk3(i) = (b*d + a*(d - b))/(a + c)/(b + c)/(c - d)
    rk4(i) = (a*(b - c) - b*c)/(c - d)/(a + d)/(b + d)
  end do

  !$omp parallel do private(i,j,k)
  do k = 1, kx
  do j = 1, jx
  do i = 1 + margin, ix - margin
    qqd(i, j, k) = &
         & +rk1(i)*qq(i + 2, j, k) &
         & + rk2(i)*qq(i + 1, j, k) &
         & + rk3(i)*qq(i - 1, j, k) &
         & + rk4(i)*qq(i - 2, j, k)
  end do
  end do
  end do
  !$omp end parallel do

  return
end subroutine d_x

!---------------------------------------------------
!===================================================
subroutine d_y(qq, y, ix, jx, kx, qqd) bind(c)
!===================================================
  use iso_c_binding
  use omp_lib
  implicit none

  integer :: i, j, k
  real(8) :: a, b, c, d
  integer, parameter :: margin = 2
  integer(c_int), intent(in) :: ix, jx, kx
  real(c_double), dimension(ix, jx, kx), intent(in) :: qq
  real(c_double), dimension(jx), intent(in) :: y
  real(c_double), dimension(ix, jx, kx), intent(out) :: qqd

  real(4), dimension(jx) :: rk1, rk2, rk3, rk4

!---------------------------------------------------

  do j = 1 + margin, jx - margin
    a = y(j + 2) - y(j)
    b = y(j + 1) - y(j)
    c = y(j) - y(j - 1)
    d = y(j) - y(j - 2)
    rk1(j) = (c*d - b*(c + d))/(a - b)/(a + c)/(a + d)
    rk2(j) = (a*(c + d) - c*d)/(a - b)/(b + c)/(b + d)
    rk3(j) = (b*d + a*(d - b))/(a + c)/(b + c)/(c - d)
    rk4(j) = (a*(b - c) - b*c)/(c - d)/(a + d)/(b + d)
  end do

  !$omp parallel do private(i,j,k)
  do k = 1, kx
  do j = 1 + margin, jx - margin
  do i = 1, ix
    qqd(i, j, k) = &
         & +rk1(j)*qq(i, j + 2, k) &
         & + rk2(j)*qq(i, j + 1, k) &
         & + rk3(j)*qq(i, j - 1, k) &
         & + rk4(j)*qq(i, j - 2, k)
  end do
  end do
  end do
  !$omp end parallel do

  return
end subroutine d_y

!---------------------------------------------------
!===================================================
subroutine d_z(qq, z, ix, jx, kx, qqd) bind(c)
!===================================================
  use iso_c_binding
  use omp_lib
  implicit none

  integer :: i, j, k
  real(8) :: a, b, c, d
  integer, parameter :: margin = 2

  integer(c_int), intent(in) :: ix, jx, kx
  real(c_double), dimension(ix, jx, kx), intent(in) :: qq
  real(c_double), dimension(kx), intent(in) :: z
  real(c_double), dimension(ix, jx, kx), intent(out) :: qqd

  real(4), dimension(kx) :: rk1, rk2, rk3, rk4

!---------------------------------------------------

  do k = 1 + margin, kx - margin
    a = z(k + 2) - z(k)
    b = z(k + 1) - z(k)
    c = z(k) - z(k - 1)
    d = z(k) - z(k - 2)
    rk1(k) = (c*d - b*(c + d))/(a - b)/(a + c)/(a + d)
    rk2(k) = (a*(c + d) - c*d)/(a - b)/(b + c)/(b + d)
    rk3(k) = (b*d + a*(d - b))/(a + c)/(b + c)/(c - d)
    rk4(k) = (a*(b - c) - b*c)/(c - d)/(a + d)/(b + d)
  end do

  !$omp parallel do private(i,j,k)
  do j = 1, jx
  do k = 1 + margin, kx - margin
  do i = 1, ix
    qqd(i, j, k) = &
         & +rk1(k)*qq(i, j, k + 2) &
         & + rk2(k)*qq(i, j, k + 1) &
         & + rk3(k)*qq(i, j, k - 1) &
         & + rk4(k)*qq(i, j, k - 2)
  end do
  end do
  end do
  !$omp end parallel do

  return
end subroutine d_z

!---------------------------------------------------