!============================================================
subroutine spherical2cartesian(qqs,rr,th,ph,ixs,jxs,kxs,ixc,jxc,kxc,qqc,xc,yc,zc) bind(C)
!============================================================
    ! Currently applicable only to uniform grid
    use iso_c_binding
    implicit none

    integer :: i,j,k
    integer :: ic,jc,kc
    integer :: ii0,jj0,kk0,ii1,jj1,kk1
    integer(c_int), intent(in) :: ixs,jxs,kxs
    integer(c_int), intent(in) :: ixc,jxc,kxc

    ! original spherical geometry
    real(c_double), dimension(ixs,jxs,kxs), intent(in) :: qqs
    real(c_double), dimension(ixs), intent(in) :: rr ! radius
    real(c_double), dimension(jxs), intent(in) :: th ! co-latitude
    real(c_Double), dimension(kxs), intent(in) :: ph ! longitude

    real(8) :: drr,dth,dph
    real(8) :: drr0,dth0,dph0,drr1,dth1,dph1

    ! converted Cartesian geometry
    real(c_double), dimension(ixc,jxc,kxc), intent(out) :: qqc
    real(c_double), dimension(ixc), intent(out) :: xc
    real(c_double), dimension(jxc), intent(out) :: yc
    real(c_double), dimension(kxc), intent(out) :: zc

    real(8) :: xmax,xmin,ymax,ymin,zmax,zmin
    real(8) :: dxc,dyc,dzc
    real(8) :: rrc,thc,phc
    
    ! Grid spacing in spherical geometry
    drr = rr(2) - rr(1)
    dth = th(2) - th(1)
    dph = ph(2) - ph(1)

    ! set Cartesian geometry
    xmax = maxval(rr)
    xmin = -xmax
    ymax = xmax
    ymin = xmin
    zmax = xmax
    zmin = xmin

    dxc = (xmax - xmin)/real(ixc-1)
    xc(1) = xmin
    do i = 2,ixc
        xc(i) = xc(i-1) + dxc
    enddo

    dyc = (ymax - ymin)/real(jxc-1)
    yc(1) = ymin
    do j = 2,jxc
        yc(j) = yc(j-1) + dyc
    enddo

    dzc = (zmax - zmin)/real(kxc-1)
    zc(1) = zmin
    do k = 2,kxc
        zc(k) = zc(k-1) + dzc
    enddo

    ! Convert Spherical to Cartesian
    do k = 1,kxc
    do j = 1,jxc
    do i = 1,ixc
        ! Position of Cartesian grid in spherical geometry
        rrc = sqrt(xc(i)**2 + yc(j)**2 + zc(k)**2)
        thc = atan2(sqrt(xc(i)**2 + yc(j)**2),zc(k))
        phc = atan2(yc(j),xc(i))

    if(     rrc < rr(ixs-1) .and. rrc > rr(1) .and. &
        &   thc < th(jxs-1) .and. thc > th(1) ) then !.and. &
!        &   phc < ph(kxs-1) .and. phc > ph(1) ) then
        ic = int((rrc - rr(1))/drr) + 1
        jc = int((thc - th(1))/dth) + 1
        kc = int((phc - ph(1))/dph) + 1
        
        drr0 = rrc - rr(ic)
        dth0 = thc - th(jc)
        dph0 = phc - ph(kc)
        drr1 = drr - drr0
        dth1 = dth - dth0
        dph1 = dph - dph0

        ii0 = ic
        ii1 = ic + 1
        jj0 = jc
        jj1 = jc + 1
        kk0 = kc
        kk1 = kc + 1

        if(kc == kxc) then
            kk1 = 1
        endif

        qqc(i,j,k) = &
             (  qqs(ii0,jj0,kk0)*drr1*dth1*dph1 &
             & +qqs(ii1,jj0,kk0)*drr0*dth1*dph1 &
             & +qqs(ii0,jj1,kk0)*drr1*dth0*dph1 &
             & +qqs(ii1,jj1,kk0)*drr0*dth0*dph1 &
             & +qqs(ii0,jj0,kk1)*drr1*dth1*dph0 &
             & +qqs(ii1,jj0,kk1)*drr0*dth1*dph0 &
             & +qqs(ii0,jj1,kk1)*drr1*dth0*dph0 &
             & +qqs(ii1,jj1,kk1)*drr0*dth0*dph0 )/(drr*dth*dph)
        else
            qqc(i,j,k) = 0.d0
        endif
    enddo
    enddo
    enddo

    return    
end subroutine spherical2cartesian

