
program lSysTestWithExLib

   use mod_aux_precision, only: WP
   use mod_lis_linear_solver, only: lis_linear_solver_t
   use mod_aux_string_functions, only: itoa, clength
   use TOAST, only: TESTCASE, PRINTSUMMARY
   use mod_aux_csr_matrix, only: csr_matrix_t

   implicit none

   integer,  parameter    :: nnz  = 7
   integer,  parameter    :: nvar = 5
   real(WP), parameter    :: tol_for_real = 1.0E-12
   integer                :: indx_in(1:nnz)
   integer                :: ptr_in(1:nvar+1)
   integer                :: i
   integer                :: ierr
   integer                :: arg_i
   real(WP)               :: val_in(1:nnz)
   real(WP)               :: val_0(1:nvar*nvar)
   real(WP)               :: uin(1:nvar)
   real(WP)               :: rhs(1:nvar)
   real(WP)               :: uout(1:nvar)
   character(len=clength) :: arg_str
   !
   type(lis_linear_solver_t)    :: lsolver
   type(TESTCASE)               :: TEST

   !> Initial arrays for testing the solver
   val_0(:) = [ 2.0_WP, -2.0_WP, 0.0_WP, 0.0_WP, 3.0_WP,  &
                0.0_WP,  2.0_WP, 0.0_WP, 0.0_WP, 0.0_WP,  &
                0.0_WP,  0.0_WP, 2.0_WP, 0.0_WP, 0.0_WP,  &
                0.0_WP,  0.0_WP, 0.0_WP, 2.0_WP, 0.0_WP,  &
                0.0_WP,  0.0_WP, 0.0_WP, 0.0_WP, 2.0_WP     ]

   val_in(:) = pack( val_0, val_0 /= 0.0_WP )
   indx_in(:) = [0, 1, 4, 1, 2, 3, 4 ]
   ptr_in(:)  = [0, 3, 4, 5, 6, nnz]

   !>  Build a LIS solver for a system of linear equations by using the given arrays for CSR format
   call lsolver%Build( nvar, ptr_in, indx_in, val_in, file_name="example_2_input.cfg" )
!   call this%Make_lis_matrix( nvar, ptr, indx, val )
   !> Initial condition and right-hand-side for testing the solver
   rhs(:) = [3.0_WP, 2.0_WP, 2.0_WP, 2.0_WP, 2.0_WP]
   uin(:) = [0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP]

   !> Solve the system with given right-hand-side and the initial condition arrays.
   call lsolver%Solve( uin(:), rhs(:), uout(:) )

   do i = 1, size( uout )
      call TEST%ASSERTEQUAL( uout(i), 1.0_wp, abs_tol=tol_for_real,   &
      &                      message=itoa(i) // "-th element did not match" )
   end do

   !> The summary of the test function
   call PRINTSUMMARY( test )

end program lSysTestWithExLib

