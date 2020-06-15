    module variables
    
    implicit none
    
    integer, parameter :: num_grid = 101
    real(8), parameter :: len_grid = 1.0d0
    real(8), parameter :: len_cell = len_grid/dble(num_grid-1)
     
    real(8), parameter :: rho_density = 1.204d0 ! kg/m3
    real(8), parameter :: c_heat_capacity = 1.006d4 ! J/kg/K
    
    real(8) :: Temperatures_K(num_grid) ! K
    real(8) :: Positions(num_grid)      ! m
    real(8) :: k_heat_conductivity(num_grid) = 2.559d-2 ! J/s/m/K
    
    
    real(8) :: a_E(num_grid)
    real(8) :: a_W(num_grid)
    real(8) :: a_P0(num_grid)
    real(8) :: a_P(num_grid)
    real(8) :: a_B(num_grid)
    
    real(8), parameter :: f_weight = 1.0d0
    
    end module variables