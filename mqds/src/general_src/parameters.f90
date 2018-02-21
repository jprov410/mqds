MODULE parameters
  USE kinds
  IMPLICIT NONE

  REAL(dp), PARAMETER :: pi=4.0_dp*DATAN(1.0_dp),    & !4.0d0*atan(1.d0),
       kb=1.3806488e-23_dp,                      & ! Boltzmann constant J/K
       planck=6.626068e-34_dp,                   & ! Planck's constant in Js
       hbar=planck/2.e0_dp/pi,                   & ! hbar in Js
       avog = 6.0221415e23_dp,                   & ! Avogadro's number
       light = 2.99792458e10_dp,                 & ! Speed of light (SI) in cm/s
       ! *************** ENERGY *************** !
       joules2wvnbr = 1.e0_dp/planck/light,      & ! Joules (SI) to wavenumber
       wvnbr2joules = 1.e0_dp/joules2wvnbr,      & ! wavenumber to Joules
       ! **************** LENGTH ************** !
       m2ang = 1.e10_dp,                         & ! meter to Angstrom
       ang2m = 1.e0_dp/m2ang,                    & ! Angstrom to meter
       ! ***************** RMASS *************** !
       amu2kg = 1.66053892e-27_dp,               & ! atomic rmass unit to kg
       ! ************** TIME ****************** !
       s2ps = 1.0e12_dp,                          & ! seconds to picoseconds
       ps2s = 1.0_dp/s2ps,                      & ! picoseconds to seconds
       s2fs = 1.0e15_dp,                          & ! seconds to femtoseconds
       fs2s = 1.0_dp/s2fs,                      & ! femtoseconds to seconds
       ! ************** FREQUENCY ************* !
       wvnbr2Hz = light,                      & ! wavenumbers to Hertz
       Hz2wvnbr = 1.e0_dp/wvnbr2Hz,              & ! Hertz to wavenumbers
       ! ***************** 2 a.u. ************* !
       au2joules = 4.35974417e-18_dp,            & ! a.u. energy to Joules
       joules2au = 1.0_dp/au2joules,            & ! Joules to a.u. energy
       ev2au=0.03674932_dp,                    & ! eV to a.u. energy
       au2ev=1.0_dp/ev2au,                      & ! a.u. energy to eV
       autime2s = 2.418884326505e-17_dp,         & ! a.u. time to seconds
       s2autime = 1.0_dp/autime2s,              & ! seconds to a.u. time
       autemp2K = 3.1577464e5_dp,                & ! a.u. temperature to Kelvin
       K2autemp = 1.0_dp/autemp2K,              & ! Kelvin to a.u. temperature
       au2ang = 0.52917725_dp,                 & ! atomic length to Angstrom
       ang2au = 1.0_dp/au2ang                    ! Angstrom to atomic length
  
  COMPLEX(dp), PARAMETER :: eye=(0.0_dp, 1.0_dp)
  
END MODULE parameters
