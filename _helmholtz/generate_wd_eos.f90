! Generate white dwarf EOS table using Helmholtz EOS
! Output format for RNS: density[g/cm^3]  pressure[dyn/cm^2]  spec_enthalpy [cm^2 /s^2], baryon density [cm^-3]

      integer          ionmax, imax, i
      parameter (ionmax = 2)
      double precision xmass(ionmax),aion(ionmax),zion(ionmax),temp,dens,abar,zbar
      double precision pres, eint, cs, enthalpy, cp, gammac, gammae
      double precision rho_min, rho_max, log_rho_min, log_rho_max, log_rho, baryondens, amu
      
! Set composition for white dwarf
! Pure He white dwarf (adjust as needed)
!      xmass(1) = 1.0d0;    aion(1)  = 4.0d0;  zion(1) = 2.0d0

! Alternative compositions:
! Pure C white dwarf: xmass(1) = 1.0; aion(1) = 12.0; zion(1) = 6.0
! Pure O white dwarf: xmass(1) = 1.0; aion(1) = 16.0; zion(1) = 8.0
! 50/50 C/O : set ionmax = 2
      xmass(1) = 0.5; aion(1) = 12.0; zion(1) = 6.0
      xmass(2) = 0.5; aion(2) = 16.0; zion(2) = 8.0

! Calculate average atomic weight and charge
      abar   = 1.0d0/sum(xmass(1:ionmax)/aion(1:ionmax))
      zbar   = abar * sum(xmass(1:ionmax) * zion(1:ionmax)/aion(1:ionmax))

! Set temperature (white dwarfs are typically cooler)
      temp = 1.0d4  ! 10^4 K (adjust as needed for your WD model)

! Density range for white dwarf
      rho_min = 1.0d0   ! g/cm^3 (surface)
      rho_max = 1.0d10  ! g/cm^3 (center of massive WD)
      imax = 1000        ! Number of density points

! Logarithmic spacing in density
      log_rho_min = log10(rho_min)
      log_rho_max = log10(rho_max)

! Read the Helmholtz data table
      call read_helm_table

! Open output file
      open(unit=10, file='wd_eos_helmholtz.dat', status='replace')
      write(10,'(a)') '# White dwarf EOS from Helmholtz code'
      write(10,'(a,f6.1,a,f6.1)') '# Composition: abar = ', abar, ', zbar = ', zbar
      write(10,'(a,es10.3,a)') '# Temperature = ', temp, ' K'
      write(10,'(a)') '# density[g/cm^3]  pressure[dyn/cm^2]  internal energy density [erg/cm^3]'
      write(10,'(a)') '# Note internal energy density = rho e_int excludes rest mass energy contribution'

! Generate EOS table
      do i = 1, imax
        ! Logarithmic spacing in density
        log_rho = log_rho_min + (real(i-1) / real(imax-1)) * (log_rho_max - log_rho_min)
        dens = 10.0d0**log_rho
        
        ! Call Helmholtz EOS
        call helm_wrapper(dens, temp, abar, zbar, pres, eint, cs, enthalpy, cp, gammac, gammae)
      
        ! Calculate baryon density (cm^-3)
!        amu = 1.66053906660d-24  ! Atomic mass unit in g
!        baryondens = dens / (abar * amu)
        ! Write to file: density, pressure, internal energy density
        write(10,'(4es15.6)') dens, pres, dens * eint 
      end do

      close(10)

      print *, 'EOS table generated: wd_eos_helmholtz.dat'
      print *, 'Density range: ', rho_min, ' to ', rho_max, ' g/cm^3'
      print *, 'Number of points: ', imax
      print *, 'Composition: abar =', abar, ', zbar =', zbar
      print *, 'Temperature =', temp, ' K'

      end
