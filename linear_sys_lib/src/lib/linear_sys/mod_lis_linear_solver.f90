module mod_lis_linear_solver
!***************************************************************************************************
!> This module contains the declaration of a type for LIS solver.
!> This type is a realization of the abstract type for linear system solvers (base_linear_solver_t)
!
!***************************************************************************************************
#include "lisf.h"

   use mod_base_linear_solver
   use mod_aux_precision
   use mod_aux_csr_matrix

   implicit none
   private
   public :: lis_linear_solver_t

   !> Module parameters
   integer,  parameter :: CLENGTH     = 30
   integer,  parameter :: IZERO       = 0
   integer,  parameter :: ITWO        = 2
   integer,  parameter :: MAXITER_DIF = 1000



   !> Options that will be aggregated into the LIS solver
   type :: lis_options_t
      real(WP) :: tol            = epsilon (0.0_WP)
      integer  :: maxiter        = MAXITER_DIF
      integer  :: mat_block_size = ITWO
      integer  :: print_option   = 0
      logical  :: is_zero_init   = .true.
      logical  :: is_verbose     = .false.
      logical  :: add_Schwarz    = .false.
      character(len=:), allocatable :: solver_name
      character(len=:), allocatable :: solver_aux_option
      character(len=:), allocatable :: convergence_cond
      character(len=:), allocatable :: mat_form
      character(len=:), allocatable :: pcond_name
      character(len=:), allocatable :: pcond_aux_option
      character(len=:), allocatable :: scale
      character(len=:), allocatable :: add_Schw_opttion
   end type lis_options_t



   type, extends(base_linear_solver_t) :: lis_linear_solver_t

      type(lis_options_t), public :: options
      LIS_MATRIX :: mat_a
      LIS_VECTOR :: vec_x
      LIS_VECTOR :: vec_b
      LIS_SOLVER :: solver
      logical    :: has_matrices = .false.

   contains

      ! All overridden procedures are private unless abstract type make them public
!      private
      ! Overridden procedures
      procedure :: Set_options_from_file         => Set_lis_opt_from_file
      procedure :: Build_lsys_solver_wth_mat     => Build_lis_solver_wth_mat
      procedure :: Build_lsys_solver_wth_arr     => Build_lis_solver_wth_arr
      procedure :: Build_lsys_solver_wthout_ops  => Build_lis_solver_wo_ops
      procedure :: Keep_lsys_ops_and_solve       => Keep_lis_ops_and_solve
      procedure :: Update_lsys_wth_mat_and_solve => Update_lis_wth_mat_and_solve
      procedure :: Update_lsys_wth_arr_and_solve => Update_lis_wth_arr_and_solve
      ! Local procedures
      ! All local procedures are private unless otherwise is set
      procedure :: Make_string_of_lis_ops
      procedure :: Make_lis_matrix
      procedure :: Make_lis_solver
      procedure :: Destroy_all_operators
      final :: Finalize_lis_solver
   end type lis_linear_solver_t



   !> Interfaces of the sub-module procedure for LIS library
   interface

      module subroutine Build_lis_solver_wth_arr( this, nvar, ptr, indx, val, file_name )
         class(lis_linear_solver_t), intent(inout) :: this
         character(len=*), optional, intent(in)    :: file_name
         integer,                    intent(in)    :: nvar
         integer,  contiguous,       intent(in)    :: indx(:)
         integer,  contiguous,       intent(in)    :: ptr(:)
         real(WP), contiguous,       intent(in)    :: val(:)
      end subroutine Build_lis_solver_wth_arr


      module subroutine Build_lis_solver_wth_mat( this, nvar, mat, file_name )
         class(lis_linear_solver_t), intent(inout) :: this
         class(csr_matrix_t),        intent(in)    :: mat
         character(len=*), optional, intent(in)    :: file_name
         integer,                    intent(in)    :: nvar
      end subroutine Build_lis_solver_wth_mat


      module subroutine Set_lis_opt_from_file( this, file_name )
         class(lis_linear_solver_t), intent(inout) :: this
         character(len=*),           intent(in)    :: file_name
      end subroutine  Set_lis_opt_from_file


      module subroutine Keep_lis_ops_and_solve( this, xi, rhs, xo )
         class(lis_linear_solver_t), intent(in)  :: this
         real(WP),                   intent(in)  :: xi(:)
         real(WP),                   intent(in)  :: rhs(:)
         real(WP),                   intent(out) :: xo(:)
      end subroutine Keep_lis_ops_and_solve


      module subroutine Update_lis_wth_mat_and_solve( this, mat, xi, rhs, xo )
         class(lis_linear_solver_t), intent(inout) :: this
         class(csr_matrix_t),        intent(in)    :: mat
         real(WP),                   intent(in)    :: xi(:)
         real(WP),                   intent(in)    :: rhs(:)
         real(WP),                   intent(out)   :: xo(:)
      end subroutine Update_lis_wth_mat_and_solve


      module subroutine Update_lis_wth_arr_and_solve( this, nvar, ptr, indx, val, xi, rhs, xo)
         class(lis_linear_solver_t), intent(inout) :: this
         integer,                    intent(in   ) :: nvar
         integer,                    intent(in   ) :: indx(:)
         integer,                    intent(in   ) :: ptr(:)
         real(WP),                   intent(in   ) :: val(:)
         real(WP),                   intent(in   ) :: rhs(:)
         real(WP),                   intent(in   ) :: xi(:)
         real(WP),                   intent(  out) :: xo(:)
      end subroutine Update_lis_wth_arr_and_solve


      module subroutine Make_string_of_lis_ops( this, str_options )
         class(lis_linear_solver_t),    intent(inout) :: this
         character(len=:), allocatable, intent(inout) :: str_options
      end subroutine Make_string_of_lis_ops


      module subroutine Make_lis_matrix( this, nvar, ptr, indx, val )
         class(lis_linear_solver_t), intent(inout) :: this
         integer,              intent(in) :: nvar
         integer,  contiguous, intent(in) :: indx(:)
         integer,  contiguous, intent(in) :: ptr(:)
         real(WP), contiguous, intent(in) :: val(:)
      end subroutine Make_lis_matrix


      module subroutine Make_lis_solver( this, fname )
         class(lis_linear_solver_t), intent(inout) :: this
         character(len=*),           intent(in)    :: fname
      end subroutine Make_lis_solver


      module subroutine Finalize_lis_solver( this )
         type(lis_linear_solver_t), intent (in) :: this
      end subroutine Finalize_lis_solver


      module subroutine Build_lis_solver_wo_ops( this, file_name )
         class(lis_linear_solver_t), intent(inout) :: this
         character(len=*), optional, intent(in)    :: file_name
      end subroutine Build_lis_solver_wo_ops


      module subroutine Destroy_all_operators( this )
         class(lis_linear_solver_t), intent(inout) :: this
      end subroutine Destroy_all_operators


   end interface

!contains
!
!
!   subroutine Build_lis_solver_wth_arr( this, nvar, ptr, indx, val, file_name )
!   !************************************************************************************************
!   !> Creates the linear system solver from the given arrays
!   !
!   !************************************************************************************************
!      use mod_aux_string_functions, only: CBLANK
!
!      class(lis_linear_solver_t), intent(inout) :: this
!      character(len=*), optional, intent(in)    :: file_name
!      integer,                    intent(in)    :: nvar
!      integer,  contiguous,       intent(in)    :: indx(:)
!      integer,  contiguous,       intent(in)    :: ptr(:)
!      real(WP), contiguous,       intent(in)    :: val(:)
!
!      !> Local
!      integer                                   :: ierr
!      character(len=:), allocatable             :: fname
!
!      !> lis starts
!      call LIS_INITIALIZE( ierr )
!
!!      !> Matrix (operator)
!      call this%Make_lis_matrix( nvar, ptr, indx, val )
!
!      !> The  vector for the solution
!      call LIS_VECTOR_CREATE( izero, this%vec_x, ierr )
!      call LIS_VECTOR_SET_SIZE( this%vec_x, izero, nvar, ierr )
!!
!      !> The vector for the right hand side
!      call LIS_VECTOR_CREATE( izero, this%vec_b, ierr )
!      call LIS_VECTOR_SET_SIZE( this%vec_b, izero, nvar, ierr )
!!
!      !> Should the solver be made from the options given by a file or by its defaults?
!      if ( present( file_name ) ) then
!         fname = trim( adjustl( file_name ) )
!      else
!         fname = CBLANK
!      end if
!      !> Make solver with the decided options
!      call this%Make_lis_solver( fname )
!
!      this%has_matrices = .true.
!   end subroutine Build_lis_solver_wth_arr
!   !************************************************************************************************

end module mod_lis_linear_solver
