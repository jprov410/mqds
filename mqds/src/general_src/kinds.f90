!> Module that defines the double precision and single
!! precision kinds.
MODULE kinds
  SAVE

  !single-precision kind
  INTEGER, PARAMETER :: sp = selected_real_kind( 6, 37 )

  !double-precision kind
  INTEGER, PARAMETER :: dp = selected_real_kind( 15, 307 )

END MODULE kinds

