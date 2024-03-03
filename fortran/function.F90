! This file is part of GTC version 4.5 
! GTC version 4.5 is released under the 3-Clause BSD license:

! Copyright (c) 2002,2010,2016, GTC Team (team leader: Zhihong Lin, zhihongl@uci.edu)
! All rights reserved.

! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the following conditions are met:

! 1. Redistributions of source code must retain the above copyright notice, 
!    this list of conditions and the following disclaimer.

! 2. Redistributions in binary form must reproduce the above copyright notice, 
!    this list of conditions and the following disclaimer in the documentation 
!    and/or other materials provided with the distribution.

! 3. Neither the name of the GTC Team nor the names of its contributors may be 
!    used to endorse or promote products derived from this software without 
!    specific prior written permission.
! ==============================================================================

!> \file
!> numerical functions related to generic smoothing and spline interpolation
!> \author Z. Lin, L. Shi
!> \date 01/24/2017


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODIFIED BY Lei Shi to have unified interface for 1D, 2D, and 3D quadratic spline
! Date: Nov. 19, 2016
! Modified by Lei Shi to have generic 1D smooth function
! Date: Jan 24, 2017
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!##############################################################################
!####  MODULE SPLINE FUNCTIONS
!##############################################################################
!> \brief quadratic spline interpolation functions
!> \details generic 1D, 2D, and 3D quadratic spline interpolation functions are
!> provided. Specific evaluators for equilibrium quantities are also provided. 
module spline_function
  implicit none

!  public
!  private:: spline0, dspline0, spline1, dspline1, spline_periodic, dspline_periodic, &
!       spline2d_dep, construct_spline2d_dep, construct_spline0, construct_spline1, &
!       construct_spline_periodic, construct_spline2d_periodic_dep
  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Generic spline evaluation functions !!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !> \brief Generic 1D spline evaluations
  !> \details Evaluate f(x) at given x
  !> By default, spline starts from x=0,y=0, and get up to x=(nx-1)*delx,
  !> out of bound x values are not checked, and the function is evaluated using
  !> the last cell spline coefficients.

  !> Boundary condition needs to be provided by programmer
  !> BC 0: use spline0
  !> BC 1: use spline1
  !> BC 2: use spline_periodic
  !> \sa construct_spline1d

  !> @param[in] x
  !> evaluation location

  !> @param[in] deriv
  !> derivative order: 0 for f(x), 1 for f'(x)

  !> @param[in] nx
  !> total number of grids in spline mesh

  !> @param[in] delx
  !> spline grid step size

  !> @param[in] y
  !> spline coefficients array

  !> @param[in] bcx
  !> Boundary condition 
  !> - 0: normal spline with calculated 1st derivative
  !> - 1: special 1st cell square root spline
  !> - 2: periodic boundary condition
  real(8) function spline1d(x, deriv, nx, delx, y, bcx)
    integer, intent(in) :: nx, deriv, bcx
    real(8), intent(in) :: x, delx, y(3,nx)
    integer i
    real(8) dx
    
    if (bcx == 0) then
       if (deriv==0) spline1d = spline0(x, nx, delx, y)
       if (deriv==1) spline1d = dspline0(x, nx, delx, y)
    else if (bcx == 1) then
       if (deriv==0) spline1d = spline1(x, nx, delx, y)
       if (deriv==1) spline1d = dspline1(x, nx, delx, y)
    else if (bcx == 2) then
       if (deriv==0) spline1d = spline_periodic(x, nx, delx, y)
       if (deriv==1) spline1d = dspline_periodic(x, nx, delx, y)
    else
       write(*,*) 'Spline1d Error: Invalid boundary condition bcx=', bcx
       stop
    endif

  end function spline1d


  !> normal 1D spline evaluator 
  real(8) function spline0(pdum,nsp,delx,y)
    integer nsp,i
    real(8) pdum,y(3,nsp),delx,dpx

    i=max(1,min(nsp-1,ceiling(pdum/delx)))
    dpx=pdum-delx*real(i-1,8)
    spline0=y(1,i)+dpx*y(2,i)+dpx*dpx*y(3,i)

  end function spline0

  !> derivative of normal 1D spline function
  real(8) function dspline0(pdum,nsp,delx,y)
    integer nsp,i
    real(8) pdum,y(3,nsp),delx,dpx

    i=max(1,min(nsp-1,ceiling(pdum/delx)))
    dpx=pdum-delx*real(i-1,8)
    dspline0=y(2,i)+2.0*dpx*y(3,i)

  end function dspline0

  !> 1D spline with first point being linear function y=sqrt(x)
  real(8) function spline1(pdum,nsp,delx,y)
    integer nsp,i
    real(8) pdum,y(3,nsp),delx,dpx,dp2

    i=max(1,min(nsp-1,ceiling(pdum/delx)))
    dpx=pdum-delx*real(i-1,8)
  ! expand y(x) using sprt(x) near x=0
    if(i==1)dpx=sqrt(dpx)
    dp2=dpx*dpx
    spline1=y(1,i)+dpx*y(2,i)+dp2*y(3,i)

  end function spline1

  !> derivative of 1D spline with 1st cell x being \f$ \sqrt{\psi} \f$
  real(8) function dspline1(pdum,nsp,delx,y)
    integer nsp, i
    real(8) pdum,y(3,nsp),delx,dpx

    i=max(1,min(nsp-1,ceiling(pdum/delx)))
    dpx=pdum-delx*real(i-1,8)
    if(i==1)dpx=sqrt(dpx)

    if(i==1)then
       dspline1=0.5*y(2,i)/dpx+y(3,i) ! y(x)=y1+y2*sprt(x)+y3*x near x=0
    else
       dspline1=y(2,i)+2.0*dpx*y(3,i) ! y(x)=y1+y2*x+y3*x*x otherwise
    endif

  end function dspline1

  !> \brief periodic spline evaluation
  !> \details Same as spline0 but out of range values are wrapped back in to its periodic
  !> image location.
  real(8) function spline_periodic(x, nx, delx, y)
    integer, intent(in) :: nx
    real(8), intent(in) :: x, delx, y(3,nx)
    integer i
    real(8) dx, x_ev
    x_ev = modulo(x, (nx-1)*delx)
    i=max(1,min(nx-1,ceiling(x/delx)))
    dx=x-delx*real(i-1, 8)
    spline_periodic = y(1,i) + dx*y(2,i) + dx*dx*y(3,i)
  end function spline_periodic

  !> \brief derivative of 1D spline defined on periodic domain
  real(8) function dspline_periodic(x, nx, delx, y)
    integer, intent(in) :: nx
    real(8), intent(in) :: x, delx, y(3,nx)
    integer i
    real(8) dx, x_ev
    x_ev = modulo(x, (nx-1)*delx)
    i=max(1,min(nx-1,ceiling(x/delx)))
    dx=x-delx*real(i-1, 8)
    dspline_periodic = y(2,i) + 2*dx*y(3,i)
  end function dspline_periodic

  !#############################################################################


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> \brief Generic 2D spline function
  !> \details  Evaluate f(x,y) or its partial derivatives using the created spline
  !> coefficients for f
  !> By default, spline starts from x=0,y=0, and get up to x=(nx-1)*delx,
  !> and y=(ny-1)*dely. Out of range (x,y) are not checked

  !> \author Lei Shi

  !> @param[in] x
  !> the first coordinate where the value of f is evaluated

  !> @param[in] y
  !> the second coordinate where the value of f is evaluated

  !> @param[in] deriv 
  !> derivative flag.
  !> value | meaning
  !> ----- | -------
  !> 0     | \f$ f(x,y) \f$
  !> 1     | \f$ \partial f/ \partial x \f$
  !> 2     | \f$ \partial f/ \partial y \f$

  !> @param[in] nx
  !> number of spline nodes in x

  !> @param[in] ny
  !> number of spline nodes in y

  !> @param[in] delx
  !> stepsize for each spline cell in x

  !> @param[in] dely
  !> stepsize for each spline cell in y

  !> @param[in] f
  !> 2D spline coefficient array constructed by construct_spline2d subroutine

  !> @param[in] bcx
  !> Boundary condition in x used when f is created.

  !> @param[in] bcy
  !> Boundary condition in y used when f is created.
  !> \warning Boundary Condition consistency is not checked! It is the programer's
  !> responsibility to make sure bcx,bcy passed in is the same as those
  !> used for construct_spline2d.
  real(8) function spline2d(x, y, deriv, nx, ny, delx, dely, f, bcx, bcy)
    integer, intent(in):: deriv, nx, ny, bcx, bcy
    real(8), intent(in):: x, y, delx, dely
    real(8), dimension(:,:,:), intent(in):: f
    integer i,j
    real(8) dx, dx2, dy, dy2
    
  ! locate the indexes for the spline cell, out of range values are not checked!
    i = max(1, min(nx-1, ceiling(x/delx)))
    j = max(1, min(ny-1, ceiling(y/dely)))
    dx = x - delx*real(i-1, 8)
    dy = y - dely*real(j-1, 8)
  ! Boundary Condition 1 requires first cell in form f = f1+f2*sqrt(x)+f3*x
    if (i==1 .and. bcx==1) dx=sqrt(dx)
    if (j==1 .and. bcy==1) dy=sqrt(dy)
    dx2 = dx*dx
    dy2 = dy*dy
    
    spline2d=0.0_8
    if (deriv==0) then
  ! evaluate f(x,y)
       spline2d = f(1,i,j)+f(2,i,j)*dx + f(3,i,j)*dx2 +&
            ( f(4,i,j)+f(5,i,j)*dx + f(6,i,j)*dx2)*dy +&
            (f(7,i,j) + f(8,i,j)*dx + f(9,i,j)*dx2)*dy2
    else if (deriv == 1) then
  ! evaluate df/dx
       if (i==1 .and. bcx==1) then
               spline2d = 0.5_8*(f(2,i,j) + f(5,i,j)*dy + f(8,i,j)*dy2)/max(dx,1e-6) +&
               (f(3,i,j) + f(6,i,j)*dy + f(9,i,j)*dy2)
       else
          spline2d = (f(2,i,j) + f(5,i,j)*dy + f(8,i,j)*dy2) +&
               2.0_8*(f(3,i,j) + f(6,i,j)*dy + f(9,i,j)*dy2)*dx
       endif
    else if (deriv == 2) then
       if (j==1 .and. bcy==1) then
          spline2d = 0.5_8*(f(4,i,j) + f(5,i,j)*dx + f(6,i,j)*dx2)/dy +&
               (f(7,i,j) + f(8,i,j)*dx + f(9,i,j)*dx2)
       else
          spline2d = (f(4,i,j) + f(5,i,j)*dx + f(6,i,j)*dx2) +&
               2.0_8*(f(7,i,j) + f(8,i,j)*dx + f(9,i,j)*dx2)*dy
       endif
    else
       write(*, *) "Spline2d Error: wrong derivative flag deriv=",deriv
       stop
    endif
  end function spline2d

  !> \deprecated
  !> Old 2D/3D spline evaluator. Should be removed when all changes have been made to other parts of the code.
  real(8) function spline2d_dep(iflag,sd,pdum,tdum,nsp,nst,delp,delt,f)
    integer iflag,nsp,nst,i,j,is,sd
    real(8) pdum,tdum,f(sd,nsp,nst),dpx,dp2,dtx,dt2,delp,delt,dx(sd),tempres

    i=max(1,min(nsp-1,ceiling(pdum/delp)))
    dpx=pdum-delp*real(i-1,8)
    ! Boundary Condition 1 corresponds to construct_spline1, first cell uses sqrt(x)
    if(i==1)dpx=sqrt(dpx)
    dp2=dpx*dpx

    j=max(1,min(nst-1,ceiling(tdum/delt)))
    dtx=tdum-delt*real(j-1,8)
    dt2=dtx*dtx

    dx=0.0_8
    dx(1)=1.0_8
    dx(2)=dpx
    dx(3)=dp2
    dx(4:6)=dx(1:3)*dtx
    dx(7:9)=dx(1:3)*dt2

    if(iflag==0)then !2D spline value
       spline2d_dep=f(1,i,j)    +f(2,i,j)*dpx    +f(3,i,j)*dp2 &
            +f(4,i,j)*dtx+f(5,i,j)*dtx*dpx+f(6,i,j)*dtx*dp2 &
            +f(7,i,j)*dt2+f(8,i,j)*dt2*dpx+f(9,i,j)*dt2*dp2
    elseif(iflag==1)then !derivative with respect to x
       if(i==1)then
          spline2d_dep=0.5_8*(f(2,i,j)+f(5,i,j)*dtx+f(8,i,j)*dt2)/dpx+f(3,i,j)+f(6,i,j)*dtx+f(9,i,j)*dt2
       else
          spline2d_dep=f(2,i,j)+f(5,i,j)*dtx+f(8,i,j)*dt2+2.0_8*dpx*(f(3,i,j)+f(6,i,j)*dtx+f(9,i,j)*dt2)
       endif
    elseif(iflag==2) then !derivative with respect to y
       spline2d_dep=f(4,i,j)+f(5,i,j)*dpx+f(6,i,j)*dp2+2.0_8*dtx*(f(7,i,j)+f(8,i,j)*dpx+f(9,i,j)*dp2)
    elseif(iflag==3) then !derivative with respect to zeta
       if(sd==27)then
         dx(10:18)=dx(1:9)
         tempres=0
         do is = 10, 18
            tempres=tempres+f(is,i,j)*dx(is)
         enddo
         spline2d_dep = tempres
       else
         spline2d_dep=0.0_8
       endif
    endif
  end function spline2d_dep
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> \brief Generic 3D spline function
  !> \details Evaluate f(x,y,z) or its partial derivatives using the created spline 
  !> coefficients for f
  !> By default, spline starts from x=0,y=0,z=0 and get up to x=(nx-1)*delx,
  !> and y=(ny-1)*dely, z=(nz-1)*delz. Out of range (x,y,z) are not checked

  !> @author Lei Shi

  !> @param[in] x
  !> the first coordinate where the value of f is evaluated

  !> @param[in] y
  !> the second coordinate where the value of f is evaluated

  !> @param[in] z
  !> the third coordinate where the value of f is evaluated

  !> @param[in] deriv 
  !> derivative flag.
  !> value | meaning
  !> ----- | -------
  !> 0     | \f$ f(x,y) \f$
  !> 1     | \f$ \partial f/ \partial x \f$
  !> 2     | \f$ \partial f/ \partial y \f$
  !> 3     | \f$ \partial f/ \partial z \f$

  !> @param[in] nx
  !> number of spline nodes in x

  !> @param[in] ny
  !> number of spline nodes in y

  !> @param[in] nz
  !> number of spline nodes in z

  !> @param[in] delx
  !> stepsize for each spline cell in x

  !> @param[in] dely
  !> stepsize for each spline cell in y

  !> @param[in] delz
  !> stepsize for each spline cell in z

  !> @param[in] f
  !> 3D spline coefficient array constructed by construct_spline3d subroutine

  !> @param[in] bcx
  !> Boundary condition in x used when f is created.

  !> @param[in] bcy
  !> Boundary condition in y used when f is created.

  !> @param[in] bcz
  !> Boundary condition in z used when f is created.

  !> \warning Boundary Condition consistency is not checked! It is the programer's
  !> responsibility to make sure bcx,bcy passed in is the same as those
  !> used for construct_spline3d.
  real(8) function spline3d(x, y, z, deriv, nx, ny, nz, delx, dely, delz,  f, bcx, bcy, bcz)
    integer, intent(in):: deriv, nx, ny, nz, bcx, bcy, bcz
    real(8), intent(in):: x, y, z, delx, dely, delz, f(27, nx, ny, nz)
    integer i, j, k
    real(8) dx, dxinv, dy, dyinv, dz, dzinv, dvec(27)
    
  ! locate the indexes for the spline cell, out of range values are not checked!
    i = max(1, min(nx-1, ceiling(x/delx)))
    j = max(1, min(ny-1, ceiling(y/dely)))
    k = max(1, min(nz-1, ceiling(z/delz)))
    dx = x - delx*real(i-1, 8)
    dy = y - dely*real(j-1, 8)
    dz = z - delz*real(k-1, 8)
    dxinv = 1.0_8/dx
    dyinv = 1.0_8/dy
    dzinv = 1.0_8/dz

  ! Boundary Condition 1 requires first cell in form f = f1+f2*sqrt(x)+f3*x
    if (i==1 .and. bcx==1) dx=sqrt(dx)
    if (j==1 .and. bcy==1) dy=sqrt(dy)
    if (k==1 .and. bcz==1) dz=sqrt(dz)
  ! Construct the spline vector
  ! take care of the boundary condition
    if (deriv==0) then
       dvec(1) = 1.0_8
       dvec(2) = dx
       dvec(3) = dx*dx
       dvec(4:6) = dvec(1:3)*dy
       dvec(7:9) = dvec(4:6)*dy
       dvec(10:18) = dvec(1:9)*dz
       dvec(19:27) = dvec(10:18)*dz
    else if (deriv==1) then
  ! partial derivative in x
       if (bcx==1 .and. i==1) then
          dvec(1) = 0.0_8
          dvec(2) = 0.5_8*dxinv
          dvec(3) = 1.0_8
          dvec(4:6) = dvec(1:3)*dy
          dvec(7:9) = dvec(4:6)*dy
          dvec(10:18) = dvec(1:9)*dz
          dvec(19:27) = dvec(10:18)*dz
       else
          dvec(1) = 0.0_8
          dvec(2) = 1.0_8
          dvec(3) = 2.0_8*dx
          dvec(4:6) = dvec(1:3)*dy
          dvec(7:9) = dvec(4:6)*dy
          dvec(10:18) = dvec(1:9)*dz
          dvec(19:27) = dvec(10:18)*dz
       endif
    else if (deriv == 2) then
  ! partial derivative in y
       if (bcy==1 .and. j==1) then
          dvec(1:3) = 0.0_8
          dvec(7) = 1.0_8
          dvec(8) = dx
          dvec(9) = dx*dx
          dvec(4:6) = dvec(7:9)*0.5_8*dyinv
          dvec(10:18) = dvec(1:9)*dz
          dvec(19:27) = dvec(10:18)*dz
       else
          dvec(1:3) = 0.0_8
          dvec(4) = 1.0_8
          dvec(5) = dx
          dvec(6) = dx*dx
          dvec(7:9) = dvec(4:6)*2.0_8*dy
          dvec(10:18) = dvec(1:9)*dz
          dvec(19:27) = dvec(10:18)*dz
       endif
    else if (deriv == 3) then
       if (bcz==1 .and. k==1) then
          dvec(1:9)=0.0_8
          dvec(19) = 1.0_8
          dvec(20) = dx
          dvec(21) = dx*dx
          dvec(22:24) = dvec(19:21)*dy
          dvec(25:27) = dvec(22:24)*dy
          dvec(10:18) = dvec(19:27)*dzinv*0.5_8
       else
          dvec(1:9)=0.0_8
          dvec(10) = 1.0_8
          dvec(11) = dx
          dvec(12) = dx*dx
          dvec(13:15) = dvec(10:12)*dy
          dvec(16:18) = dvec(13:15)*dy
          dvec(19:27) = dvec(10:18)*dz*2.0_8
       endif
    else
       write(*, *) "Spline2d Error: wrong derivative flag deriv=",deriv
       stop
    endif

    spline3d = sum(f(:,i,j,k)*dvec(:))

  end function spline3d

  !#######################################################################
  !> \defgroup construct_spline ''Generic spline construction''
  !> @{

  !> Generic construct 1D spline

  !> @param[in] nx
  !> total spline grid number. The cell number is nx-1

  !> @param[in] delx
  !> spline step size

  !> @param[in,out] f
  !> spline coefficients array. f(1,1:nx) is passed in as the given function values on grids.
  !> f(1:3, 1:nx) are then calculated to represent the spline function

  !> @param[in] bcx
  !> Boundary condition.

  !> Three kinds of boundary conditions are available for each dimension
  !>   - BC 0: Use 341 formula to obtain the estimated derivative at x=0
  !>   - BC 1: Use 341 formula to obtain the estimated derivative at x=xmax, and use
  !>   \f[ 
  !>   f(x) = a + b\sqrt{x} + cx 
  !>   \f]
  !>   formula for the first cell near x=0.
  !>   - BC 2: Periodic condition, df/dx(0)=df/dx(xmax)
  subroutine construct_spline1d(nx, delx, f, bcx)
    integer, intent(in):: nx, bcx
    real(8), intent(in):: delx
!f2py intent(in,out)  f(3,nx)
    real(8), intent(inout):: f(3,nx)

    if (bcx == 0) then
       call construct_spline0(0,nx,delx,f)
    else if (bcx == 1) then
       call construct_spline1(nx,delx,f)
    else if (bcx == 2) then
       call construct_spline_periodic(nx, delx, f)
    else
       print *,'Error: Wrong spline boundary condition choice:', bcx
       stop
    endif
  end subroutine construct_spline1d

  !> Generic construct 2D spline

  !> @param[in] nx
  !> total spline grid number in x. The cell number is nx-1

  !> @param[in] ny
  !> total spline grid number in y. The cell number is ny-1

  !> @param[in] delx
  !> spline step size in x

  !> @param[in] dely
  !> spline step size in y

  !> @param[in,out] f
  !> spline coefficients array. f(1,1:nx,1:ny) is passed in as the given function values on grids.
  !> f(1:9,1:nx,1:ny) are then calculated to represent the spline function

  !> @param[in] bcx
  !> Boundary condition in x. 

  !> @param[in] bcy
  !> Boundary condition in y. 

  !> \sa construct_spline1d
  subroutine construct_spline2d(nx,ny,delx,dely,f, bcx, bcy)
    integer, intent(in):: nx, ny, bcx, bcy
    real(8), intent(in):: delx, dely
!f2py intent(in,out) f(9,nx,ny)
    real(8), intent(inout):: f(9,nx,ny)
    integer i,j,s
    real(8) ddum1(3,nx),ddum2(3,ny)

    ddum1=0.0
    ddum2=0.0
    
  ! periodic condition needs to be enforced if bc is set to 2
    if (bcx==2) then
       f(1,1,:) = 0.5_8*(f(1,1,:)+f(1,nx,:))
       f(1,nx,:) = f(1,1,:)
    endif
    if (bcy==2) then
       do j = 1, ny-1
          ddum1(1,:)=f(1,:,j)
          call construct_spline1d(nx, delx, ddum1, bcx)
          f(1,:,j)=ddum1(1,:) ! the function value may be altered due to smoothing
          f(2,:,j)=ddum1(2,:)
          f(3,:,j)=ddum1(3,:)
       enddo
       f(1,:,ny) = f(1,:,1)
       f(2,:,ny) = f(2,:,1)
       f(3,:,ny) = f(3,:,1)
    else
       do j = 1, ny
          ddum1(1,:)=f(1,:,j)
          call construct_spline1d(nx, delx, ddum1, bcx)
          f(1,:,j)=ddum1(1,:) ! the function value may be altered due to smoothing
          f(2,:,j)=ddum1(2,:)
          f(3,:,j)=ddum1(3,:)
       enddo
    endif

    do i = 1, nx
       do s = 1, 3
          ddum2(1,:)=f(s,i,:)
          call construct_spline1d(ny, dely, ddum2, bcy)
          f(s,i,:)=ddum2(1,:)
          f(s+3,i,:)=ddum2(2,:)
          f(s+6,i,:)=ddum2(3,:)
       enddo
    enddo
  end subroutine construct_spline2d

  !> Generic construct 3D spline

  !> @param[in] nx
  !> total spline grid number in x. The cell number is nx-1

  !> @param[in] ny
  !> total spline grid number in y. The cell number is ny-1

  !> @param[in] nz
  !> total spline grid number in z. The cell number is nz-1

  !> @param[in] delx
  !> spline step size in x

  !> @param[in] dely
  !> spline step size in y

  !> @param[in] delz
  !> spline step size in z

  !> @param[in,out] f
  !> spline coefficients array. f(1,:,:,:) is passed in as the given function values on grids.
  !> f(1:9,:,:,:) are then calculated to represent the spline function

  !> @param[in] bcx
  !> Boundary condition in x. 

  !> @param[in] bcy
  !> Boundary condition in y. 

  !> @param[in] bcz
  !> Boundary condition in z. 

  !> \sa construct_spline1d
  subroutine construct_spline3d(nx,ny,nz,delx,dely,delz,f,bcx,bcy,bcz)
    integer, intent(in):: nx, ny, nz, bcx, bcy, bcz
    real(8), intent(in):: delx, dely, delz
    real(8), intent(inout):: f(27,nx,ny,nz)

    integer i,j,k,s
    real(8) temp1d(3,nz), temp2d(9,nx,ny)

  ! construct 2D spline on each x-y plane
  ! periodic in z will be enforced if bcz==2

    if (bcz==2) then
       do i=1,nz-1
          temp2d(1,:,:) = f(1,:,:,i)
          call construct_spline2d(nx,ny,delx,dely,temp2d,bcx,bcy)
          f(1:9,:,:,i) = temp2d(:,:,:)
       enddo
       f(1:9,:,:,nz)=f(1:9,:,:,1)
    else
       do i=1,nz
          temp2d(1,:,:) = f(1,:,:,i)
          call construct_spline2d(nx,ny,delx,dely,temp2d,bcx,bcy)
          f(1:9,:,:,i) = temp2d(:,:,:)
       enddo
    endif

  ! generate spline coefficient on z direction
    do i=1,nx
       do j=1,ny
          do s=1,9
             temp1d(1,:) = f(s,i,j,:)
             call construct_spline1d(nz, delz, temp1d, bcz)
             f(s,i,j,:) = temp1d(1,:)
             f(s+9,i,j,:) = temp1d(2,:)
             f(s+18,i,j,:) = temp1d(3,:)
          enddo
       enddo
    enddo
    
  end subroutine construct_spline3d

  !> construct 1D spline with on x=[0,xmax], grid x_i = (i-1)*delx, delx=xmax/(nsp-1)
  !> in domain i, y = y(1,i) + y(2,i)*delx + y(3,i)*delx**2
  subroutine construct_spline0(iflag,nsp,delx,y)
    integer i,nsp,ipp,iflag
    real(8) delx,y(3,nsp)

  ! first point
    if(iflag==0)then
  ! iflag=0: first point being y=y1+y2*x+y3*x*x
       y(2,1)=(4.0_8*y(1,2)-y(1,3)-3.0_8*y(1,1))/(2.0_8*delx)
       y(3,1)=(y(1,2)-y(1,1)-y(2,1)*delx)/(delx*delx)

    elseif(iflag==1)then
  ! iflag=1: first point being linear function y=y_1+y2*x
       y(2,1)=(y(1,2)-y(1,1))/delx
       y(3,1)=0.0_8

    elseif(iflag==2)then
  ! iflag=2: first point being quadratic function y=y1+y3*x*x
       y(2,1)=0.0_8
       y(3,1)=(y(1,2)-y(1,1))/(delx*delx)
    endif

    do i=2,nsp-2
       ipp=min(i+2,nsp)
       y(2,i)=-y(2,i-1)+2.0_8*(y(1,i)-y(1,i-1))/delx

  ! smooth f1
       y(1,i+1)=0.5_8*delx*y(2,i)+0.25_8*y(1,ipp)+0.75_8*y(1,i)
    enddo

    y(2,nsp-1)=-y(2,nsp-2)+2.0_8*(y(1,nsp-1)-y(1,nsp-2))/delx
    y(2,nsp)=-y(2,nsp-1)+2.0_8*(y(1,nsp)-y(1,nsp-1))/delx

    do i=2,nsp-1
       y(3,i)=(y(2,i+1)-y(2,i))/(2.0_8*delx)
    enddo

  ! last point is not used;
    y(3,nsp)=0.0_8

  end subroutine construct_spline0

  !> spline1 for cases with first point being function y=y1+y2*sqrt(x)+y3*x
  !> Detailed documentation can be found in Doc/Developer_Notes/periodic_spline.pdf Appendix C
  subroutine construct_spline1(nsp,delx,f)
    integer i,nsp,ipp
    real(8) delx,f(3,nsp)

    ! first point
    f(2,1)=(2.0_8*f(1,2)-f(1,3)-f(1,1))/((2.0_8-sqrt(2.0_8))*sqrt(delx))
    f(3,1)=(f(1,2)-f(1,1)-f(2,1)*sqrt(delx))/delx

    ! second point
    f(2,2)=0.5_8*f(2,1)/sqrt(delx)+f(3,1)
    f(3,2)=(f(1,3)-f(1,2)-delx*f(2,2))/(delx*delx)

    do i=3,nsp-2
       ipp=min(i+2,nsp)
       f(2,i)=-f(2,i-1)+2.0_8*(f(1,i)-f(1,i-1))/delx

       ! smooth f1
       f(1,i+1)=0.5_8*delx*f(2,i)+0.25_8*f(1,ipp)+0.75_8*f(1,i)
    enddo

    f(2,nsp-1)=-f(2,nsp-2)+2.0_8*(f(1,nsp-1)-f(1,nsp-2))/delx
    f(2,nsp)=-f(2,nsp-1)+2.0_8*(f(1,nsp)-f(1,nsp-1))/delx

    do i=3,nsp-1  !!change start from 1 to 3
       f(3,i)=(f(2,i+1)-f(2,i))/(2.0_8*delx)
    enddo
  ! last point is not used;
    f(3,nsp)=0.0_8
  end subroutine construct_spline1

  !> Construct 1D Quadratic Spline using periodic boundary condition
  !>
  !> Check the documentation in Doc/Developers_Notes/periodic_spline.pdf for a detailed discussion
  !> \author Lei Shi
  !> \date Oct 16, 2016

  !> @param[in] nsp
  !> spline grid number

  !> @param[in] delx
  !> spline cell size

  !> @param[in,out] y
  !> spline coefficients. y(1,:) is passed in as given values.
  subroutine construct_spline_periodic(nsp,delx,y)
    integer, intent(in) :: nsp
    real(8), intent(in) :: delx
    real(8), intent(inout):: y(3, nsp)
  ! local variables
    integer :: i
    real(8) :: dxinv
    real(8) :: delta

    ! Periodic boundary condition requires the first order derivative to be continuous at the end points
    ! Check the documentation for a detailed discussion about the spline algorithm

    dxinv = 1.0_8/delx

    ! First, let's check if the input data satisfies the requirements
    ! Check grid point number
    if (mod(nsp, 2) /= 0) then
       write(*,*) 'Periodic Spline Error: Total number of data must be even. nsp=',nsp,' doesn''t work. Use normal spline or request a new input data set.'
       return
    endif

    ! check periodicity in raw data
    if (abs(y(1,1)-y(1,nsp)) > SPACING(y(1,1))+SPACING(y(1,nsp)) ) then
       write(*,*) 'Periodic Spline Error: function values at end points must be the same. y(1)=',y(1,1),', y(n)=',y(1,nsp),' doesn''t work. Use normal spline or request a new input data set.'
       write(*,*) 'SPACING(y(1,1)) = ', SPACING(y(1,1)), 'SPACING(y(1,nsp)) = ', SPACING(y(1,nsp))
       return
    endif
    ! reset all the higher order spline coefficients
    y(2:3,:) = 0.0_8
    ! Now we let y = a + b*x + c*x**2 on each section, first step is to calculate b1

    do i=1,nsp/2-1
       y(2,1) = y(2,1) + y(1,2*i) - y(1, 2*i+1)
    enddo
    y(2,1) = y(2,1)*2.0_8*dxinv
    ! we use y(2,nsp) as a buffer to save the periodic b1
    y(2, nsp)=y(2,1)
    ! Then, we can obtain bi one by one
    do i=2, nsp-1
       y(2, i) = 2.0_8*(y(1,i) - y(1,i-1))*dxinv - y(2,i-1)
    enddo
!    delta = 0.0_8
!    do i=1,nsp/2-1
!       delta = delta - y(2,i*2+1) + y(2,i*2)
!    enddo
!    do i=1,(nsp+1)/2-1
!       delta = delta + y(2,i*2) - y(2,i*2-1)
!    enddo
!    delta = delta / (nsp-2) / 2.0_8
!    do i=0,nsp/2-1
!       y(2,i*2+1) = y(2,i*2+1) + delta
!    enddo
!    do i=1,(nsp+1)/2-1
!       y(2,i*2) = y(2,i*2) - delta
!    enddo
!    y(2,nsp) = y(2,1)
   

    ! Then the ci's
    do i=1, nsp-1
       y(3, i) = (y(2, i+1)-y(2, i))*0.5_8*dxinv
    enddo
    y(3,nsp) = y(3,1)
  end subroutine construct_spline_periodic

  !======================================================

  !> inversion of spline0 y(x) to x(y)
  subroutine invert_spline0(iflag,nsp,delx,dely,y,x)
    integer i,nsp,j,iflag
    real(8) delx,dely,y(3,nsp),x(3,nsp),ydum,y0,y1,y2

  ! first point given by input
    x(1,1)=0.0_8
  ! other points
    do i=2,nsp-1
  ! y grid
       ydum=dely*real(i-1,8)
  ! search x grid for ydum
       j=1
       do while (ydum>y(1,j+1))
          j=j+1
       enddo

  ! x(1,i)=j grid location + distance from j grid
       y0=y(1,j)
       y1=y(2,j)
       y2=y(3,j)
       if (abs(y2)>0.000001_8) then
         x(1,i)=delx*real(j-1,8)+(sqrt(y1*y1+4.0_8*y2*(ydum-y0))-y1)/(2.0_8*y2)
       else
         x(1,i)=delx*real(j-1,8)+(ydum-y0)/y1
       endif

    enddo

  ! last point
    x(1,nsp)=x(1,1)+delx*real(nsp-1,8)

  ! spline fit x function
    if(iflag==0)call construct_spline0(0,nsp,dely,x)

  end subroutine invert_spline0

  !> inversion of spline1 y(x) to x(y)
  subroutine invert_spline1(nsp,delx,dely,y,x)
    integer i,nsp,j
    real(8) delx,dely,y(3,nsp),x(3,nsp),ydum,y0,y1,y2

  ! first point is given by inputs
    x(1,1)=0.0_8
  ! other points
    do i=2,nsp-1
  ! y grid
       ydum=dely*real(i-1,8)
  ! search x grid for ydum
       j=1
       do while (ydum>y(1,j+1))
          j=j+1
       enddo

  ! x(1,i)=j grid location + distance from j grid
       y0=y(1,j)
       y1=y(2,j)
       y2=y(3,j)

       if (abs(y2)>0.000001_8) then
         x(1,i)=delx*real(j-1,8)+(sqrt(y1*y1+4.0_8*y2*(ydum-y0))-y1)/(2.0_8*y2)
       else
         x(1,i)=delx*real(j-1,8)+(ydum-y0)/y1
       endif

       if (j==1) x(1,i)=x(1,i)*x(1,i)

    enddo
  ! last point
    x(1,nsp)=delx*real(nsp-1,8)

  ! spline fit x function
  ! call spline with first point being ~ sqrt(r)
    call construct_spline1(nsp,dely,x)

  end subroutine invert_spline1

  !> \deprecated
  !> legacy construction of 2D spline coefficients. Should be removed after other parts are changed.
  subroutine construct_spline2d_dep(nx,ny,delx,dely,y)
    integer i,j,s,nx,ny,ipp
    real(8) delx,dely,y(9,nx,ny),ddum1(3,nx),ddum2(3,ny)

    ddum1=0.0_8
    ddum2=0.0_8

    do j = 1, ny
       ddum1(1,:)=y(1,:,j)     
       call construct_spline1(nx,delx,ddum1)
       y(1,:,j)=ddum1(1,:)
       y(2,:,j)=ddum1(2,:)
       y(3,:,j)=ddum1(3,:)
    enddo
    do i = 2, nx
       do s = 1, 3
          ddum2(1,:)=y(s,i,:)
          call construct_spline0(0,ny,dely,ddum2)
          y(s,i,:)=ddum2(1,:)
          y(s+3,i,:)=ddum2(2,:)
          y(s+6,i,:)=ddum2(3,:)
       enddo
    enddo
  end subroutine construct_spline2d_dep
  !!!!!!!!!!!!!

  !> \deprecated
  !> should be removed after other parts are changed.
  !> Construct 2D spline with the second dimension periodic

  !> \author Lei Shi
  !> \date Oct 16, 2016
  !> \sa construct_spline_periodic
  subroutine construct_spline2d_periodic_dep(nx,ny,delx,dely,z)
    integer i,j,s,nx,ny,ipp
    real(8) delx,dely,z(9,nx,ny),ddum1(3,nx),ddum2(3,ny)

    ddum1=0.0_8
    ddum2=0.0_8

    do j = 1, ny-1
       ddum1(1,:)=z(1,:,j)
       call construct_spline1(nx,delx,ddum1)
       z(1,:,j)=ddum1(1,:)
       z(2,:,j)=ddum1(2,:)
       z(3,:,j)=ddum1(3,:)
    enddo
    ! Values at periodic boundary are forced to be equal
    z(1,:,ny)=z(1,:,1)
    z(2,:,ny)=z(2,:,1)
    z(3,:,ny)=z(3,:,1)

    do i = 2, nx
       do s = 1, 3
          ddum2(1,:)=z(s,i,:)
          call construct_spline_periodic(ny,dely,ddum2)
          z(s,i,:)=ddum2(1,:)
          z(s+3,i,:)=ddum2(2,:)
          z(s+6,i,:)=ddum2(3,:)
       enddo
    enddo
  end subroutine construct_spline2d_periodic_dep

  !====================================================================
  ! 3D Data creation and spline subroutines for VMEC
  ! Written by Lei Shi, Nov 19, 2016
  !==================================================================== 
  ! NOT IMPLEMENTED YET

  !==================================================================== 


  !> Cross product
  !> \f$ \vec{v}_3 = \vec{v}_1 \times \vec{v}_2 \f$
  subroutine cross(vec1,vec2,vec3)
    real(8):: vec1(3),vec2(3),vec3(3)
    
    vec3(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2)
    vec3(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3)
    vec3(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1)
    
  end subroutine cross
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> dot product
  !>\f$ f=\vec{v}_1 \cdot \vec{v}_2 \f$

  !> @param[in] n
  !> length of the vectors
  real(8) function dots(vec1,vec2,n)
    integer i,n
    real(8) ::temp,vec1(n),vec2(n)
    temp=0.0_8
    do i=1,n
      temp=temp+vec1(i)*vec2(i)
    enddo
    dots=temp
  end function dots

end module spline_function
