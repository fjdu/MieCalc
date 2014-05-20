      program main
      
        implicit none
        !
        real xx
        complex crefin
        real mimcut
        integer nmom, IPOLZN, numang
        integer, parameter :: maxang=7, momdim=200
        real :: xmu(maxang)
        logical perfect, anyang, PRT(2)
        !
        real Qext, Qsca, GQsc, Qabs, g
        complex :: S1(maxang), S2(maxang)
        complex sforw, sback
        complex :: tforw(2), tback(2)
        real spike, ang
        real :: pmom(0:momdim, 4)
        !
        integer i
        real dmu
        !
        external miev0
      
        xx = 0.001
        crefin = (1.2D0, 0.3D0)
        perfect = .true.
        mimcut = 1e-6
        anyang = .true.
        nmom = 1
        IPOLZN = +1234
        PRT(1) = .true.
        PRT(2) = .true.
        numang = maxang
      
        !dmu = 2D0 / dble(numang-1)
        !xmu(1) = 1D0
        !do i=2, numang
        !  xmu(i) = xmu(i-1) - dmu
        !end do
        !xmu(numang) = -1D0
        !xmu(numang/2+1) = 0D0
        do i=1, numang
          ang = (i-1) * 180.0 / (numang-1)
          xmu(i) = cos(3.1415926/180.0*ang)
        end do
        !do i=1, numang
        !  write(*,*) i, xmu(i)
        !end do
        !
        call miev0(xx, crefin, perfect, mimcut, anyang, numang, xmu,
     &             nmom, IPOLZN, momdim, PRT,
     &             Qext, Qsca, GQsc, pmom,
     &             sforw, sback, S1, S2, tforw, tback, spike)
        !
        Qabs = Qext - Qsca
        g = GQsc / Qsca
        write(*, *) Qext, Qsca, Qabs, g
      
      end program main
