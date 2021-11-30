module mod_aux_precision

   implicit none

   integer, parameter :: SP = selected_real_kind( p = 6, r = 37 )
   integer, parameter :: DP = selected_real_kind( p = 15, r = 307 )
   integer, parameter :: QP = selected_real_kind( p = 33, r = 4931 )
   integer, parameter :: WP = DP

end module mod_aux_precision
