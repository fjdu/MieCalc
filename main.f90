subroutine miev0_simple(xx, indreal, indimg, &
                        Qabs, Qsca, g)

  implicit none
  !
  real, intent(in) :: xx, indreal, indimg
  real, intent(out) :: Qabs, Qsca, g
  !
  integer, parameter :: maxang=11, momdim=1
  complex crefin
  real mimcut
  integer nmom, IPOLZN, numang
  real :: xmu(maxang)
  logical perfect, anyang, PRT(2)
  !
  real Qext, GQsc
  complex :: S1(maxang), S2(maxang)
  complex sforw, sback
  complex :: tforw(2), tback(2)
  real spike
  real :: pmom(0:momdim, 4)
  !
  integer i
  real dmu, Qback
  !
  external miev0, adt_sphere, bhmie
  !
  crefin = cmplx(indreal, indimg)
  !
  perfect = .false.
  anyang = .true.
  PRT = .false.
  mimcut = 1e-8
  nmom = momdim
  IPOLZN = 0
  numang = maxang
  !
  dmu = 2e0 / real(numang-1)
  xmu(1) = 1e0
  do i=2, numang
    xmu(i) = xmu(i-1) - dmu
  end do
  xmu(numang) = -1e0
  if (mod(numang, 2) .eq. 1) then
    xmu(numang/2+1) = 0e0
  end if
  !
  if (xx .gt. 1e4) then
    call adt_sphere(xx, crefin, Qext, Qabs, Qsca, g)
  else if (((real(crefin) .lt. 0e0) .and. (.not. perfect)) .or. &
           (real(crefin) .gt. 2e2) .or. &
           (abs(aimag(crefin)) .gt. 2e2)) then
    Qext = -1.0
    Qabs = -1.0
    Qsca = -1.0
    GQsc = -1.0
    g    = sqrt(Qext-1.0)
  else
    call miev0(xx, crefin, perfect, mimcut, anyang, numang, xmu, &
               nmom, IPOLZN, momdim, PRT, &
               Qext, Qsca, GQsc, pmom, &
               sforw, sback, S1, S2, tforw, tback, spike)
    Qabs = Qext - Qsca
    g = GQsc / (Qsca + 1e-30)
  end if
  !
  if ((Qabs .le. 0.0) .and. (xx .le. 1e5)) then
    call bhmie(xx, crefin, 2, S1, S2, Qext, Qsca, Qback, g)
    Qabs = Qext - Qsca
  end if
  !
  Qabs = max(Qabs, 0.0)
  Qsca = max(Qsca, 0.0)
end subroutine miev0_simple



subroutine adt_sphere(xx, crefin, Qext, Qabs, Qsca, g)
  implicit none
  real, intent(in) :: xx
  complex, intent(in) :: crefin
  real, intent(out) :: Qext, Qabs, Qsca, g
  complex rho
  real rho1, rho2, beta, rr
  !
  rho = 2.0 * xx * (crefin - 1.0)
  rr = abs(rho)
  rho1 = real(rho)
  rho2 = aimag(rho)
  !
  if (rr .le. 1e-3) then
    Qext = 2.0
    Qabs = 4.0/3.0 * rho2
    Qsca = Qext - Qabs
    g = 1.0
  else
    !
    beta = atan2(rho2, rho1)
    !
    Qext = 2.0 + 4.0/(rr*rr) * &
      (cos(2.0*beta) - &
       exp(-rho2) * &
       (cos(rho1-2.0*beta) - rr*sin(rho1-beta)))
    Qabs = 1.0 + (exp(-2.0*rho2)*(2.0*rho2+1.0) - 1.0) / (2.0*rho2*rho2)
    Qsca = Qext - Qabs
    g = 1.0
  end if
end subroutine adt_sphere
