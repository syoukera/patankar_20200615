!  patankar_20200615.f90 
!
!  FUNCTIONS:
!  patankar_20200615 - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: patankar_20200615
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program patankar_20200615
    
    use variables

    implicit none
    
    integer i
    
    call initialize_variables()
    ! call calc_cond()
    call calc_coeff()
    call calc_tdma()    

    do i = 1, num_grid
        print *, Positions(i), Temperatures_K(i)
    end do

    end program patankar_20200615
    
    ! ----------
    
    subroutine initialize_variables()
    
    use variables
    
    implicit none
    
    integer i
    integer :: pnt_bound = int(num_grid/5)
    
    ! print *, 'Pointer of Boundary is ', pnt_bound
    Temperatures_K(:pnt_bound) = 800.0d0
    Temperatures_K(pnt_bound:) = 300.0d0
    
    ! print *, num_grid, len_grid, len_cell
    
    do i = 1, num_grid
        Positions(i) = len_cell*dble(i-1)
    end do
    
    do i = 1, num_grid
        print *, Positions(i), Temperatures_K(i)
    end do
    
    end subroutine initialize_variables
    
    ! ----------
    
    subroutine calc_cond()
    
    use variables
    
    implicit none
    
    k_heat_conductivity = 1.0d0
    
    end subroutine calc_cond
    
    ! ----------

    subroutine calc_coeff()
    
    use variables
    
    implicit none
    
    integer i
    real(8) delta_x, delta_t
    
    do i = 2, num_grid-1
        a_E(i) = k_heat_conductivity(i-1)/(Positions(i) - Positions(i-1))
        a_W(i) = k_heat_conductivity(i+1)/(Positions(i+1) - Positions(i))
        delta_x = (Positions(i+1) - Positions(i-1))/2.0d0
        delta_t = 1.0d-10
        a_P0(i) = rho_density*c_heat_capacity*delta_x/delta_t
        a_P(i) = f_weight*a_E(i) + f_weight*a_W(i) + a_P0(i)
    end do
    
    !do i = 2, num_grid-1
    !    a_E(i) = k_heat_conductivity(i+1)/(Positions(i+1) - Positions(i))
    !    a_W(i) = k_heat_conductivity(i-1)/(Positions(i) - Positions(i-1))
    !    delta_x = (Positions(i+1) - Positions(i-1))/2.0d0
    !    delta_t = 1.0d0-5
    !    a_P0(i) = rho_density*c_heat_capacity*delta_x/delta_t
    !    a_P(i) = a_E(i) + a_W(i) + a_P0(i)
    !end do
    
    print *, 'Print from calc_coeff'    
    
    end subroutine calc_coeff
    
    
    ! ----------
    
    subroutine calc_tdma()
    
    use variables
    
    implicit none
    
    integer i
    
    real(8) :: a(num_grid)
    real(8) :: b(num_grid)
    real(8) :: c(num_grid)
    real(8) :: d(num_grid)
    
    real(8) :: P(num_grid)
    real(8) :: Q(num_grid)
    
    ! Upper boundary
    ! Temperature = d(1)/a(1) = 800.0d0
    a(1) = 1.0d0
    b(1) = 0.0d0
    c(1) = 0.0d0
    d(1) = 800.0d0
    
    ! Inner cells    
    a(2:num_grid-1) = a_P(2:num_grid-1)
    b(2:num_grid-1) = a_E(2:num_grid-1)*f_weight
    c(2:num_grid-1) = a_W(2:num_grid-1)*f_weight
    do i = 2, num_grid-1
        d(i) = a_E(i)*(1 - f_weight)*Temperatures_K(i+1)    &
             + a_W(i)*(1 - f_weight)*Temperatures_K(i-1)    &
             + (a_P0(i) - (1 - f_weight)*a_E(i) - (1 - f_weight)*a_W(i))*Temperatures_K(i)
    end do
    
    !! Inner cells    
    !a(2:num_grid-1) = a_P(2:num_grid-1)
    !b(2:num_grid-1) = a_E(2:num_grid-1)
    !c(2:num_grid-1) = a_W(2:num_grid-1)
    !d(2:num_grid-1) = a_P0(2:num_grid-1)*Temperatures_K(2:num_grid-1)
    
    ! Lower Boundery
    ! Temperature = d(1)/a(1) = 300.0d0
    a(num_grid) = 1.0d0
    b(num_grid) = 0.0d0
    c(num_grid) = 0.0d0
    d(num_grid) = 300.0d0    
    
    P(1) = b(1)/a(1)
    Q(1) = d(1)/a(1)
    
    do i = 2, num_grid
        P(i) = b(i)/(a(1) - c(i)*P(i-1))
        Q(i) = (d(i) + c(i)*Q(i-1))/(a(i) - c(i)*P(i-1))
    end do
    
    Temperatures_K(num_grid) = Q(num_grid)
    
    do i = num_grid-1, 1, -1
        Temperatures_K(i) = P(i)*Temperatures_K(i+1) + Q(i)
    end do
    
    print *, 'Print from calc_tdma'    
    
    end subroutine calc_tdma
