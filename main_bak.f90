program main

  implicit none
  !
  real xx
  complex crefin
  real mimcut
  integer nmom, IPOLZN, numang
  integer, parameter :: maxang=51, momdim=1
  real :: xmu(maxang)
  logical perfect, anyang, PRT(2)
  !
  real Qext, Qsca, GQsc, Qabs, g
  complex :: S1(maxang), S2(maxang)
  complex sforw, sback
  complex :: tforw(2), tback(2)
  real spike
  real :: pmom(0:momdim, 4)
  !
  integer i
  real dmu
  !
  character(len=128) filename
  real lam, indr, indi
  !
  external miev0


  xx = 0.2
  crefin = (1.0D0, -2.3D0)
  perfect = .false.
  mimcut = 1e-6
  anyang = .true.
  nmom = 1
  IPOLZN = 0
  PRT(1) = .false.
  PRT(2) = .false.
  numang = maxang

  dmu = 2D0 / dble(numang-1)
  xmu(1) = 1D0
  do i=2, numang
    xmu(i) = xmu(i-1) - dmu
  end do
  xmu(numang) = -1D0
  xmu(numang/2+1) = 0D0
  !
  call miev0(xx, crefin, perfect, mimcut, anyang, numang, xmu, &
             nmom, IPOLZN, momdim, PRT, &
             Qext, Qsca, GQsc, pmom, &
             sforw, sback, S1, S2, tforw, tback, spike)
  !
  Qabs = Qext - Qsca
  g = GQsc / Qsca
  write(*, *) Qext, Qsca, Qabs, g, spike
  do i=1, numang
    write(*, *) i, acos(xmu(i))*180/3.14, abs(S1(i)), abs(S2(i))
  end do

end program main
