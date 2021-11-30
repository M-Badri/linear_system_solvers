module mod_base_linear_solver
!***************************************************************************************************
!> This module contains the declaration of a base (abstract) type for all linear system solvers.
!> This type could be extend for both iterative and direct solvers.
!
!***************************************************************************************************

   use mod_aux_precision, only: WP
   use mod_aux_csr_matrix, only: csr_matrix_t

   implicit none
   private
   public :: base_linear_solver_t



   !> Abstract type for sparse linear system solver
   !> It will be implemented by other explicit particular solvers (either with direct or iterative)
   type, abstract :: base_linear_solver_t

   contains

      ! Only the generic interface of the procedures are public
      procedure(Set_options_from_file_interface) ,       deferred :: Set_options_from_file
      procedure(Build_lsys_solver_wth_mat_interface),    deferred :: Build_lsys_solver_wth_mat
      procedure(Build_lsys_solver_wth_array_interface),  deferred :: Build_lsys_solver_wth_arr
      procedure(Build_lsys_solver_wthout_ops_interface), deferred :: Build_lsys_solver_wthout_ops
      generic, public :: Build =>                       &
                         Build_lsys_solver_wth_mat,     &
                         Build_lsys_solver_wth_arr,     &
                         Build_lsys_solver_wthout_ops
      procedure(Keep_lsys_ops_and_solve_interface),      deferred :: Keep_lsys_ops_and_solve
      procedure(Update_lsys_wth_mat_and_solve_interface),deferred :: Update_lsys_wth_mat_and_solve
      procedure(Update_lsys_wth_arr_and_solve_interface),deferred :: Update_lsys_wth_arr_and_solve
      generic, public :: Solve =>                       &
                         Keep_lsys_ops_and_solve,       &
                         Update_lsys_wth_arr_and_solve, &
                         Update_lsys_wth_mat_and_solve
   end type base_linear_solver_t



   !> Interfaces of the deferred procedures used in the abstract type
   abstract interface

      subroutine Set_options_from_file_interface (this, file_name)
         import
         class(base_linear_solver_t), intent(inout)  :: this
         character(len=*),            intent(in)     :: file_name
      end subroutine Set_options_from_file_interface


      subroutine Build_lsys_solver_wth_mat_interface( this, nvar, mat, file_name )
         import
         class(base_linear_solver_t), intent(inout) :: this
         class(csr_matrix_t),         intent(in   ) :: mat
         character(len=*), optional,  intent(in   ) :: file_name
         integer,                     intent(in   ) :: nvar
      end subroutine Build_lsys_solver_wth_mat_interface


      subroutine Build_lsys_solver_wth_array_interface( this, nvar, ptr, indx, val, file_name )
         import
         class(base_linear_solver_t), intent(inout) :: this
         character(len=*), optional,  intent(in   ) :: file_name
         integer,                     intent(in   ) :: nvar
         integer,  contiguous,        intent(in   ) :: indx(:)
         integer,  contiguous,        intent(in   ) :: ptr(:)
         real(WP), contiguous,        intent(in   ) :: val(:)
      end subroutine Build_lsys_solver_wth_array_interface


      subroutine Build_lsys_solver_wthout_ops_interface( this, file_name )
         import
         class(base_linear_solver_t), intent(inout) :: this
         character(len=*), optional,  intent(in)    :: file_name
      end subroutine Build_lsys_solver_wthout_ops_interface


      subroutine Keep_lsys_ops_and_solve_interface( this, xi, rhs, xo )
         import
         class(base_linear_solver_t), intent(in   ) :: this
         real(WP),                    intent(in   ) :: xi(:)
         real(WP),                    intent(in   ) :: rhs(:)
         real(WP),                    intent(  out) :: xo(:)
      end subroutine Keep_lsys_ops_and_solve_interface


      subroutine Update_lsys_wth_mat_and_solve_interface( this, mat, xi, rhs, xo )
         import
         class(base_linear_solver_t), intent(inout) :: this
         class(csr_matrix_t),         intent(in   ) :: mat
         real(WP),                    intent(in   ) :: xi(:)
         real(WP),                    intent(in   ) :: rhs(:)
         real(WP),                    intent(  out) :: xo(:)
      end subroutine Update_lsys_wth_mat_and_solve_interface


      subroutine Update_lsys_wth_arr_and_solve_interface( this, nvar, ptr, indx, val, &
         &                                                xi, rhs, xo )
         import
         class(base_linear_solver_t), intent(inout) :: this
         integer,                     intent(in   ) :: nvar
         integer,                     intent(in   ) :: indx(:)
         integer,                     intent(in   ) :: ptr(:)
         real(WP),                    intent(in   ) :: val(:)
         real(WP),                    intent(in   ) :: rhs(:)
         real(WP),                    intent(in   ) :: xi(:)
         real(WP),                    intent(  out) :: xo(:)
      end subroutine Update_lsys_wth_arr_and_solve_interface

   end interface


end module mod_base_linear_solver
