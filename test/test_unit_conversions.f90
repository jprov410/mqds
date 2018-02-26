PROGRAM test_unit_conversions
    USE unit_conversions
    USE kinds
    IMPLICIT NONE
    REAL(dp) :: time=1.0_dp, temp=300.0_dp, freq=1.0_dp, energy=100.0_dp
    REAL(dp) :: beta=0.0_dp
    REAL(dp) :: eps=1.0E-10

    time = time * convert('fs','au_time')
    beta = 1.0_dp / ( temp * convert('kelvin','au_energy') )
    freq = freq * convert('wvnbr','au_ang_freq')
    energy = energy * convert('wvnbr','au_energy')

    IF ( ABS(time - 41.341373336561368) >= 41.341373336561368/1.0E6_dp ) THEN
        WRITE(*,*) 'failed unit conversion for time, outside tolerance'
        STOP 1
    ELSE IF ( ABS(beta - 1052.5834351695135) >= 1052.5834351695135/1.0E6_dp ) THEN
        WRITE(*,*) 'failed unit conversion for beta/temp, outside tolerance'
        STOP 1
    ELSE IF ( ABS(freq - 4.5563352527599625E-006) >= 4.5563352527599625E-006/1.0E6_dp ) THEN
        WRITE(*,*) 'failed unit conversion for frequency, outside tolerance'
        STOP 1
    ELSE IF ( ABS(energy - 4.5563343515983042E-004) >= 4.5563343515983042E-004/1.0E6_dp ) THEN
        WRITE(*,*) 'failed unit conversion for energy, outside tolerance'
        STOP 1
    END IF

END PROGRAM test_unit_conversions