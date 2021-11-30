module mod_aux_csr_matrix
!***************************************************************************************************
!> This module contains a type declaration for csr matrices.
!> This type has be adopted from Ayato's code.
!
!***************************************************************************************************

   use mod_aux_precision, only:WP_ => WP

   implicit none
   private

   type, public :: csr_matrix_t
      real(WP_), public, allocatable :: values(:)
      integer,   public, allocatable :: columns(:)
      integer,   public, allocatable :: row_index(:)
      integer,   public              :: nonzero_count
   end type csr_matrix_t


end module mod_aux_csr_matrix
