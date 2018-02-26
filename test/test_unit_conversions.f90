PROGRAM test_unit_conversions
    USE unit_conversions
    USE kinds
    IMPLICIT NONE
    REAL(dp) :: time=1.0_dp, temp=300.0_dp, freq=1.0_dp, energy=100.0_dp
    REAL(dp) :: beta=0.0_dp

    time = time * convert('fs','au_time')
    beta = 1.0_dp / ( temp * convert('kelvin','au_energy') )
    freq = freq * convert('wvnbr','au_ang_freq')
    energy = energy * convert('wvnbr','au_energy')

END PROGRAM test_unit_conversions