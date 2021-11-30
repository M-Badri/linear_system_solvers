submodule (mod_lis_linear_solver) smod_lis_linear_solver
!***************************************************************************************************
!> This sub-module contains the module procedures for lis_lsolver_t type.
!
!***************************************************************************************************
#include "lisf.h"

   use mod_aux_string_functions, only: CBLANK, itoa, rtoa
   use ZOO, only: FILE_INI
   implicit none



contains



   module subroutine Build_lis_solver_wth_arr( this, nvar, ptr, indx, val, file_name )
   !************************************************************************************************
   !> Creates the linear system solver from the given arrays
   !
   !************************************************************************************************
      class(lis_linear_solver_t), intent(inout) :: this
      character(len=*), optional, intent(in   ) :: file_name
      integer,                    intent(in   ) :: nvar
      integer,  contiguous,       intent(in   ) :: indx(:)
      integer,  contiguous,       intent(in   ) :: ptr(:)
      real(WP), contiguous,       intent(in   ) :: val(:)

      !> Local
      integer                                   :: ierr
      character(len=:), allocatable             :: fname

      !> lis starts
      call LIS_INITIALIZE( ierr )

      !> Matrix (operator)
      call this%Make_lis_matrix( nvar, ptr, indx, val )

      !> The  vector for the solution
      call LIS_VECTOR_CREATE( izero, this%vec_x, ierr )
      call LIS_VECTOR_SET_SIZE( this%vec_x, izero, nvar, ierr )

      !> The vector for the right hand side
      call LIS_VECTOR_CREATE( izero, this%vec_b, ierr )
      call LIS_VECTOR_SET_SIZE( this%vec_b, izero, nvar, ierr )

      !> Should the solver be made from the options given by a file or by its defaults?
      if ( present( file_name ) ) then
         fname = trim( adjustl( file_name ) )
      else
         fname = cblank
      end if
      !> Make solver with the decided options
      call this%Make_lis_solver( fname )

      this%has_matrices = .true.
   end subroutine Build_lis_solver_wth_arr
   !************************************************************************************************



   module subroutine Build_lis_solver_wth_mat( this, nvar, mat, file_name )
   !************************************************************************************************
   !> Constructor of the linear system with a given CSR matrix
   !
   !************************************************************************************************
      class(lis_linear_solver_t),   intent(inout)  :: this
      class(csr_matrix_t),          intent(in   )  :: mat
      character(len=*), optional,   intent(in   )  :: file_name
      integer,                      intent(in   )  :: nvar

      ! Local
      real(Wp),         allocatable                :: val(:)
      integer,          allocatable                :: indx(:)
      integer,          allocatable                :: ptr(:)
      integer                                      :: ierr
      integer                                      :: nnz
      character(len=:), allocatable                :: str_options
      character(len=:), allocatable                :: fname

      !> LIS starts
      call LIS_INITIALIZE( ierr )

      !> Matrix (operator)
      nnz = mat%nonzero_count
      allocate( val,  source = mat%values )
      allocate( indx, source = mat%row_index )
      allocate( ptr,  source = mat%columns )
      call this%Make_lis_matrix( nvar, ptr, indx, val )

      !> The vector for the solution
      call LIS_VECTOR_CREATE( izero, this%vec_x, ierr )
      call LIS_VECTOR_SET_SIZE( this%vec_x, izero, nvar, ierr )

      !> The vector for the right hand side
      call LIS_VECTOR_CREATE( izero, this%vec_b, ierr )
      call LIS_VECTOR_SET_SIZE( this%vec_b, izero, nvar, ierr )

      !> Should the solver be made from the options given by a file or by its defaults?
      if ( present( file_name ) ) then
         fname = trim( adjustl( file_name ) )
      else
         fname = cblank
      end if
      !> Make solver with the decided options
      call this%Make_lis_solver( fname )

      !> Free the memory
      deallocate( val, indx, ptr )

      this%has_matrices = .true.
   end subroutine Build_lis_solver_wth_mat
   !************************************************************************************************



   module subroutine Build_lis_solver_wo_ops( this, file_name )
   !************************************************************************************************
   !> Builds a solver without operator.
   !> Therefor, it does not have matrices and the subroutine "solve" must always provide the info
   !>    for the matrices.
   !
   !************************************************************************************************
      class(lis_linear_solver_t),  intent(inout) :: this
      character(len=*), optional,  intent(in   ) :: file_name

      ! Local
      integer                                    :: ierr
      character(len=:), allocatable              :: str_options
      character(len=:), allocatable              :: fname

      !> lis starts
      call LIS_INITIALIZE( ierr )

      !> Should the solver be made from the options given by a file or by its defaults?
      if ( present( file_name ) ) then
         fname = trim( adjustl( file_name ) )
      else
         fname = cblank
      end if
      !> Make solver with the decided options
      call this%Make_lis_solver( fname )

      this%has_matrices = .false.
   end subroutine Build_lis_solver_wo_ops
   !************************************************************************************************



   module subroutine Set_lis_opt_from_file( this, file_name )
   !************************************************************************************************
   !> Routine for setting the options of the solver from a configuration file (put in "inputs" dir)
   !
   !************************************************************************************************
      class(lis_linear_solver_t), intent(inout) :: this
      character(len=*),           intent(in   ) :: file_name

      ! Local
      type(FILE_INI)                            :: FINI
      integer                                   :: ierr
      character(len=CLENGTH)                    :: temp_string

      call FINI%LOAD( filename='inputs/' // trim( file_name ), error=ierr )

      call FINI%GET( section_name = 'solver',                         &
                     option_name  = 'name',                           &
                     val          = temp_string)
      this%options%solver_name = trim( adjustl( temp_string ) )


      call FINI%GET(section_name = 'solver',                          &
                     option_name  = 'auxiliary_option',               &
                     val          = temp_string)
      this%options%solver_aux_option = trim( adjustl( temp_string ) )


      call FINI%GET( section_name = 'solver',                         &
                     option_name  = 'stop_residual',                  &
                     val          = this%options%tol )

      call FINI%GET( section_name = 'solver',                         &
                     option_name  = 'maximum_iteration',              &
                     val          = this%options%maxiter )

      call FINI%GET( section_name = 'solver',                         &
                     option_name  = 'matrix_format',                  &
                     val          = this%options%mat_form )

      call FINI%GET( section_name = 'solver',                         &
                     option_name  = 'block_matrix_size',              &
                     val          = this%options%mat_block_size )

      call FINI%GET( section_name = 'solver',                         &
                     option_name  = 'is_initiated_with_zero',         &
                     val          = this%options%is_zero_init )

      call FINI%GET( section_name = 'solver',                         &
                     option_name  = 'is_verbose',                     &
                     val          = this%options%is_verbose )

      call FINI%GET( section_name = 'solver',                         &
                     option_name  = 'convergence_condition',          &
                     val          = temp_string)
      this%options%convergence_cond = trim( adjustl( temp_string ) )


      call FINI%GET( section_name = 'solver',                         &
                     option_name  = 'print_option',                   &
                     val          = this%options%print_option )

      call FINI%GET( section_name = 'preconditioner',                 &
                     option_name  = 'name',                           &
                     val          = temp_string )
      this%options%pcond_name = trim( adjustl( temp_string ) )


      call FINI%GET( section_name = 'preconditioner',                 &
                     option_name  = 'auxiliary_option',               &
                     val          = temp_string )
      this%options%pcond_aux_option = trim( adjustl( temp_string ) )


      call FINI%GET( section_name = 'additive_Schwarz',               &
                     option_name  = 'is_added',                       &
                     val          = this%options%add_Schwarz )

      call FINI%GET( section_name = 'additive_Schwarz',               &
                     option_name  = 'auxiliary_option',               &
                     val          = temp_string )
      this%options%add_Schw_opttion = trim( adjustl( temp_string ) )


      if (ierr == 0) then
         if ( this%options%is_verbose ) then
            write (*, fmt ="(/, a, /, 4x ,a)")                &
               "The options have been read from the file:",   &
               "inputs/" // trim( file_name )
         end if
      else
         write (*, fmt ="(/, a, /, 4x, a)")                   &
            "The options cannot be read from this file:",     &
            "inputs/" // trim( file_name )
         stop
      end if
   end subroutine Set_lis_opt_from_file
   !************************************************************************************************



   module subroutine Keep_lis_ops_and_solve( this, xi, rhs, xo )
   !************************************************************************************************
   !> Solve the system of the equation with the operators (matrices) already set for the solver.
   !
   !************************************************************************************************
      class(lis_linear_solver_t), intent(in   ) :: this
      real(WP),                   intent(in   ) :: xi(:)
      real(WP),                   intent(in   ) :: rhs(:)
      real(WP),                   intent(  out) :: xo(:)

      ! Local
      integer                                   :: i
      integer                                   :: ierr
      real(WP)                                  :: gval

      !> Check if the operator has already been set
      if (.not. this%has_matrices) then
         error stop "The system is being used without setting the operator properly!"
      end if

      !> Check if the input has enough data
      call Assert_lisvec_and_arr_len( this%vec_x, xi )
      call Assert_lisvec_and_arr_len( this%vec_x, xo )

      !> Fill out the initial values of the solution vector
      do i = 1, size(xi, 1)
         call LIS_VECTOR_SET_VALUE( LIS_INS_VALUE, i, xi(i), this%vec_x, ierr )
         call LIS_VECTOR_SET_VALUE( LIS_INS_VALUE, i, rhs(i),this%vec_b, ierr )
      end do

      !> Solve the system
       call LIS_SOLVE( this%mat_a, this%vec_b, this%vec_x, this%solver, ierr)

       !> Return the output in an array
       do i =1, size( xo )
          call LIS_VECTOR_GET_VALUE( this%vec_x, i, gval, ierr )
          xo(i) = gval
       end do
   end subroutine Keep_lis_ops_and_solve
   !************************************************************************************************



   module subroutine Update_lis_wth_mat_and_solve( this, mat, xi, rhs, xo )
   !************************************************************************************************
   !> Set the matrix, initial condition and RHS, then solves the system
   !> Before using this subroutine, the object must be constructed without operators.
   !
   !************************************************************************************************
      class(lis_linear_solver_t), intent(inout) :: this
      class(csr_matrix_t),        intent(in   ) :: mat
      real(WP),                   intent(in   ) :: xi(:)
      real(WP),                   intent(in   ) :: rhs(:)
      real(WP),                   intent(  out) :: xo(:)

      ! Local
      integer                                   :: nvar
      real(Wp),         allocatable             :: val(:)
      real(WP)                                  :: gval
      integer,          allocatable             :: indx(:)
      integer,          allocatable             :: ptr(:)
      integer                                   :: ierr
      integer                                   :: nnz
      integer                                   :: i
      character(len=:), allocatable             :: str_options
      character(len=:), allocatable             :: fname

      ! The size of output, initial, and RHS arrays
      nvar = size( xi )

      !> Matrix (operator)
      nnz = mat%nonzero_count
      allocate( val,  source = mat%values )
      allocate( indx, source = mat%row_index )
      allocate( ptr,  source = mat%columns )
      call this%Make_lis_matrix( nvar, ptr, indx, val )

      !> The vector for the solution
      call LIS_VECTOR_CREATE( izero, this%vec_x, ierr )
      call LIS_VECTOR_SET_SIZE( this%vec_x, izero, nvar, ierr )

      !> The vector for the right hand side
      call LIS_VECTOR_CREATE( izero, this%vec_b, ierr )
      call LIS_VECTOR_SET_SIZE( this%vec_b, izero, nvar, ierr )

      !> Fill out the initial values of the solution vector
      do i = 1, size(xi, 1)
         call LIS_VECTOR_SET_VALUE( LIS_INS_VALUE, i, xi(i), this%vec_x, ierr )
         call LIS_VECTOR_SET_VALUE( LIS_INS_VALUE, i, rhs(i),this%vec_b, ierr )
      end do

      !> Solve the system
       call LIS_SOLVE( this%mat_a, this%vec_b, this%vec_x, this%solver, ierr)

       !> Return the output in an array
       do i =1, size( xo )
          call LIS_VECTOR_GET_VALUE( this%vec_x, i, gval, ierr )
          xo(i) = gval
       end do

      ! Free the memory taken by LIS operators
      call this%Destroy_all_operators()

   end subroutine Update_lis_wth_mat_and_solve
   !************************************************************************************************



   module subroutine Update_lis_wth_arr_and_solve( this, nvar, ptr, indx, val, xi, rhs, xo )
   !************************************************************************************************
   !>
   !
   !************************************************************************************************
      class(lis_linear_solver_t), intent(inout) :: this
      integer,                    intent(in   ) :: nvar
      integer,                    intent(in   ) :: indx(:)
      integer,                    intent(in   ) :: ptr(:)
      real(WP),                   intent(in   ) :: val(:)
      real(WP),                   intent(in   ) :: rhs(:)
      real(WP),                   intent(in   ) :: xi(:)
      real(WP),                   intent(  out) :: xo(:)

      ! local
      integer                                   :: i
      integer                                   :: ierr
      real(WP)                                  :: gval

      !> The vector for the solution
      call LIS_VECTOR_CREATE( izero, this%vec_x, ierr )
      call LIS_VECTOR_SET_SIZE( this%vec_x, izero, nvar, ierr )

      !> The vector for the right hand side
      call LIS_VECTOR_CREATE( izero, this%vec_b, ierr )
      call LIS_VECTOR_SET_SIZE( this%vec_b, izero, nvar, ierr )

      !> Fill out the initial values of the solution vector
      do i = 1, nvar
         call LIS_VECTOR_SET_VALUE( LIS_INS_VALUE, i, xi(i), this%vec_x, ierr )
         call LIS_VECTOR_SET_VALUE( LIS_INS_VALUE, i, rhs(i),this%vec_b, ierr )
      end do

      call this%Make_lis_matrix( nvar, ptr, indx, val )

      call lis_vector_print (this%vec_b, ierr)

      !> Solve the system
!       call LIS_SOLVE( this%mat_a, this%vec_b, this%vec_x, this%solver, ierr)

!       !> Return the output in an array
!       do i =1, size( xo )
!          call LIS_VECTOR_GET_VALUE( this%vec_x, i, gval, ierr )
!          xo(i) = gval
!       end do
!
!      ! Free the memory taken by LIS operators
!      call this%Destroy_all_operators()
   end subroutine Update_lis_wth_arr_and_solve
   !************************************************************************************************



   module subroutine Destroy_all_operators( this )
   !************************************************************************************************
   !> Frees the memory taken by the LIS operators
   !
   !************************************************************************************************
      class(lis_linear_solver_t), intent(inout) :: this

      ! Local
      integer                                   :: ierr

       call LIS_VECTOR_DESTROY( this%vec_b, ierr )
       call LIS_VECTOR_DESTROY( this%vec_x, ierr )
       call LIS_MATRIX_DESTROY( this%mat_a, ierr )
   end subroutine Destroy_all_operators
   !************************************************************************************************



   module subroutine Make_string_of_lis_ops( this, str_options )
   !************************************************************************************************
   !> From the options given in the input file, makes a string which will be passed to LIS solver.
   !
   !************************************************************************************************
      class (lis_linear_solver_t),   intent(inout) :: this
      character(len=:), allocatable, intent(inout) :: str_options

      ! Local
      logical                                      :: p_opt

      str_options = " "

      if ( this%options%solver_name == "empty" ) then
         str_options = str_options
      else
         str_options = "-i" // cblank // trim( this%options%solver_name )
      end if

      if ( this%options%solver_aux_option == "empty" ) then
         str_options = str_options
      else if ( this%options%solver_aux_option /= "empty" &
               & .and. &
               this%options%solver_name       /= "empty" ) then
         str_options = str_options // cblank // "-" // trim( this%options%solver_aux_option )
      end if

      if ( this%options%pcond_name == "empty" ) then
         str_options = str_options
      else
         str_options = str_options // cblank // "-p" // cblank // trim( this%options%pcond_name )
      end if

      if ( this%options%pcond_aux_option == "empty" ) then
         str_options = str_options
      else if ( this%options%pcond_aux_option /= "empty" &
               & .and. &
               this%options%pcond_name       /= "empty" )  then
         str_options = str_options // cblank // "-" // trim( this%options%pcond_aux_option )
      end if

      if (this%options%maxiter == maxiter_dif)  then
         str_options = str_options
      else
         str_options = &
               str_options // cblank // "-maxiter" // cblank // trim( itoa( this%options%maxiter ) )
      end if

      if (abs( this%options%tol - epsilon (0.0_WP) ) < tiny( 0.0_WP ) ) then
         str_options = str_options
      else
         str_options = str_options // cblank // "-tol" // cblank // trim( rtoa( this%options%tol ) )
      end if

      if ( this%options%is_zero_init ) then
         str_options = str_options
      else
         str_options = str_options // cblank // "-initx_zeros"// cblank // "false"
      end if

      if ( this%options%print_option == 0 .or. &
           this%options%print_option == 1 .or. &
           this%options%print_option == 2 .or. &
           this%options%print_option == 3        ) then

           p_opt = .true.
      end if

      if ( p_opt ) then
         str_options = &
            str_options // cblank // "-print" // cblank // itoa( this%options%print_option )
      end if
   end subroutine Make_string_of_lis_ops
   !************************************************************************************************



   module subroutine Make_lis_matrix( this, nvar, ptr, indx, val )
   !************************************************************************************************
   !> Creates a (operator) matrix for the system
   !
   !************************************************************************************************
      class(lis_linear_solver_t), intent(inout) :: this
      integer,                    intent(in   ) :: nvar
      integer,  contiguous,       intent(in   ) :: indx(:)
      integer,  contiguous,       intent(in   ) :: ptr(:)
      real(WP), contiguous,       intent(in   ) :: val(:)

      !> Local
      integer                                   :: ierr

      call LIS_MATRIX_CREATE( IZERO, this%mat_a, ierr )
      call LIS_MATRIX_SET_SIZE( this%mat_a, IZERO, nvar, ierr )
      call LIS_MATRIX_SET_TYPE( this%mat_a, LIS_MATRIX_CSR, ierr )
      call LIS_MATRIX_SET_CSR( size( val, 1 ), ptr, indx, val, this%mat_a, ierr)
      call LIS_MATRIX_ASSEMBLE( this%mat_a, ierr )
   end subroutine Make_lis_matrix
   !************************************************************************************************



   module subroutine Make_lis_solver( this, fname )
   !************************************************************************************************
   !> Creates the solver using  LIS library methods
   !> It sets the option for the solver as well.
   !
   !************************************************************************************************
      class(lis_linear_solver_t), intent(inout) :: this
      character(len=*),           intent(in   ) :: fname

      ! Local
      integer                                   :: ierr
      character(len=:), allocatable             :: str_options

      !> Solver
      call LIS_SOLVER_CREATE( this%SOLVER, ierr )

      !> Make needed options either from an input or use the default ones
      if (fname /= " ") then
         call this%Set_options_from_file( fname )
         call this%Make_string_of_lis_ops( str_options )
      else
         str_options = " "
         write (*, fmt ="(/, a)") "The default / manually set options are used"
      end if

      !> Set options for the solver in LIS
      call LIS_SOLVER_SET_OPTION( str_options, this%SOLVER, ierr )
   end subroutine Make_lis_solver
   !************************************************************************************************




   module subroutine Finalize_lis_solver( this )
   !************************************************************************************************
   !>
   !
   !************************************************************************************************
      type(lis_linear_solver_t), intent(in) :: this

      ! Local
      integer                               :: ierr

      call LIS_VECTOR_DESTROY( this%vec_b, ierr )
      call LIS_VECTOR_DESTROY( this%vec_x, ierr )
      call LIS_MATRIX_DESTROY( this%mat_a, ierr )
      call LIS_SOLVER_DESTROY( this%solver, ierr )
      call LIS_FINALIZE( ierr )
   end subroutine Finalize_lis_solver
   !************************************************************************************************



   subroutine Assert_lisvec_and_arr_len( l_vec, f_array )
   !************************************************************************************************
   !> Checks if the length of an array is equal with a LIS vector
   !
   !************************************************************************************************
      real(WP), intent(in) :: f_array(:)
      LIS_VECTOR           :: l_vec

      ! Local
      integer              :: local_n, ierr, global_n

      call lis_vector_get_size( l_vec, local_n, global_n, ierr )

      if ( local_n /= size(f_array, 1) ) then
         error stop "The length of the given vector is not matched with that of the array!"
      end if
   end subroutine Assert_lisvec_and_arr_len
   !************************************************************************************************


end submodule smod_lis_linear_solver
