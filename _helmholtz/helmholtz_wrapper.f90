      subroutine helm_wrapper (dens, temp, abar, zbar, pres, eint, cs, enthalpy, cp, gammac, gammae)

! This is a simplified wrapper for the Helmholtz EOS. It takes only
! density, temperature, and composition (through abar and zbar) for a single
! cell, and computes the resulting pressures, internal energies, entropies,
! and specific heats, and returns these through the common blocks
! in the include decks as well as in the subroutine call itself.
!
! Notes on enthalpy calculation:
! - Returns relativistic specific enthalpy potential in CGS units [erg/g]
! - h_relativistic = e_specific + p/ρ + c²
! - Includes rest mass energy contribution (c² ≈ 9.0 × 10²⁰ erg/g) for use in XNS

      include 'implno.dek'
      include 'vector_eos.dek'

      double precision dens, temp, abar, zbar ! inputs
      double precision pres, eint, cs, enthalpy, h, gammac, gammae, cp
      double precision, parameter :: c = 2.99792458d10 ! speed of light in CGS units [cm/s]

      temp_row (1) = temp
      den_row (1) = dens
      abar_row (1) = abar
      zbar_row (1) = zbar

! here is the tabular helmholtz free energy eos:
!
! routine read_helm_table reads an electron helm free energy table
! routine read_helm_iontable reads an ion free energy table
! routine helmeos computes the pressure, energy and entropy via tables
! routine helmeos3 uses fermionic ions
! routine helmeos2 adds third derivatives
! routine helmeos_orig as helmeos with the table read inside the routine

      jlo_eos = 1 ; jhi_eos = 1 ! 1 entry only

! Call helmeos

      call helmeos

      cs = cs_row (1)
      eint = etot_row (1)
      pres = ptot_row (1)
!      enthalpy = etot_row (1) + ptot_row (1) / den_row (1)
      h = etot_row (1) + ptot_row (1) / den_row (1) + c**2 ! Add rest mass contribution for relativistic enthalpy
      enthalpy = c**2 * log (h / c**2) ! Potential energy contribution for enthalpy
      cp   = cp_row (1)
      gammac = gam1_gas_row(1)
      gammae = gam2_gas_row(1)

      end
