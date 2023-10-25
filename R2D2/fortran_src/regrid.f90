subroutine interp(x,y,z,xu,yu,zu,qq,ix,jx,kx,ixu,jxu,kxu,qu) bind(C)
   use iso_c_binding
   implicit none

  integer :: i,j,k
  integer(c_int), intent(in) :: ix,jx,kx,ixu,jxu,kxu
  real(c_double), dimension(ix), intent(in) :: x
  real(c_double), dimension(jx), intent(in) :: y
  real(c_double), dimension(kx), intent(in) :: z
  real(c_double), dimension(ixu), intent(in) :: xu
  real(c_double), dimension(jxu), intent(in) :: yu
  real(c_double), dimension(kxu), intent(in) :: zu
  
  real(c_double), dimension(ix ,jx ,kx ), intent(in)  :: qq
  real(c_double), dimension(ixu,jxu,kxu), intent(out) :: qu
  
  !real(8), dimension(kx ,jx ,ix ) :: qq
  !real(8), dimension(kxu,jxu,ixu) :: qu
  !intent(in)  :: qq
  !intent(out) :: qu

  integer :: imin,imax
  integer :: jmin,jmax
  integer :: kmin,kmax
  integer, dimension(ixu) :: iloc
  integer, dimension(jxu) :: jloc
  integer, dimension(kxu) :: kloc
  real(8), dimension(ixu) :: dx0,dx1
  real(8), dimension(jxu) :: dy0,dy1
  real(8), dimension(kxu) :: dz0,dz1
  
  call id_loc(x,xu,ix,ixu,imin,imax,iloc,dx0,dx1)
  call id_loc(y,yu,jx,jxu,jmin,jmax,jloc,dy0,dy1)
  call id_loc(z,zu,kx,kxu,kmin,kmax,kloc,dz0,dz1)
  
  qu = 0.d0

  do i = imin,imax
  do j = jmin,jmax
  do k = kmin,kmax
!     qu(k,j,i) = ( &
!          & + qq(kloc(k)+0,jloc(j)+0,iloc(i)+0)*dz1(k)*dy1(j)*dx1(i) &
!          & + qq(kloc(k)+0,jloc(j)+0,iloc(i)+1)*dz1(k)*dy1(j)*dx0(i) &
!          & + qq(kloc(k)+0,jloc(j)+1,iloc(i)+0)*dz1(k)*dy0(j)*dx1(i) &
!          & + qq(kloc(k)+0,jloc(j)+1,iloc(i)+1)*dz1(k)*dy0(j)*dx0(i) &
!          & + qq(kloc(k)+1,jloc(j)+0,iloc(i)+0)*dz0(k)*dy1(j)*dx1(i) &
!          & + qq(kloc(k)+1,jloc(j)+0,iloc(i)+1)*dz0(k)*dy1(j)*dx0(i) &
!          & + qq(kloc(k)+1,jloc(j)+1,iloc(i)+0)*dz0(k)*dy0(j)*dx1(i) &
!          & + qq(kloc(k)+1,jloc(j)+1,iloc(i)+1)*dz0(k)*dy0(j)*dx0(i) &
!          & )/(dx0(i) + dx1(i))/(dy0(j) + dy1(j))/(dz0(k) + dz1(k))

       qu(i,j,k) = ( &
          & + qq(iloc(i)+0,jloc(j)+0,kloc(k)+0)*dz1(k)*dy1(j)*dx1(i) &
          & + qq(iloc(i)+0,jloc(j)+0,kloc(k)+1)*dz0(k)*dy1(j)*dx1(i) &
          & + qq(iloc(i)+0,jloc(j)+1,kloc(k)+0)*dz1(k)*dy0(j)*dx1(i) &
          & + qq(iloc(i)+0,jloc(j)+1,kloc(k)+1)*dz0(k)*dy0(j)*dx1(i) &
          & + qq(iloc(i)+1,jloc(j)+0,kloc(k)+0)*dz1(k)*dy1(j)*dx0(i) &
          & + qq(iloc(i)+1,jloc(j)+0,kloc(k)+1)*dz0(k)*dy1(j)*dx0(i) &
          & + qq(iloc(i)+1,jloc(j)+1,kloc(k)+0)*dz1(k)*dy0(j)*dx0(i) &
          & + qq(iloc(i)+1,jloc(j)+1,kloc(k)+1)*dz0(k)*dy0(j)*dx0(i) &
          & )/(dx0(i) + dx1(i))/(dy0(j) + dy1(j))/(dz0(k) + dz1(k))
  enddo
  enddo
  enddo

     
  return
end subroutine interp

subroutine id_loc(x,xu,ix,ixu,imin,imax,iloc,dx0,dx1) 
   use iso_c_binding
   implicit none
  
  integer :: i,ii
  integer(c_int), intent(in) :: ix,ixu
  real(c_double), dimension(ix), intent(in) :: x
  real(c_double), dimension(ixu), intent(in) :: xu

  integer, intent(out) :: imin,imax
  integer, dimension(ixu), intent(out) :: iloc
  real(8), dimension(ixu), intent(out) :: dx0,dx1

  ! serch minmum
  do i = 1,ixu
     if(xu(i) > x(1)) then
        imin = i
        goto 2000
     endif
  enddo
2000 continue

  ! serch maximum
  do i = ixu,1,-1
     if(xu(i) < x(ix)) then
        imax = i
        goto 3000
     endif
  enddo
3000 continue
    
  do i = 1,ixu
     do ii = ix,1,-1
        if(xu(i) >= x(ii)) then
           iloc(i) = ii
           dx0(i) = xu(i) - x(ii)
           dx1(i) = x(ii+1) - xu(i)
           goto 1000
        endif
     enddo
1000 continue
  enddo
  
  return
end subroutine id_loc
