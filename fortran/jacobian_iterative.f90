program jacobi
implicit none

!these will hold our final answer, if we converge
double precision 					:: current_x1
double precision 					:: current_x2
double precision 					:: current_x3
double precision 					:: current_x4
double precision 					:: current_x5

!used in calculation for Xn 
double precision 					:: previous_x1
double precision 					:: previous_x2
double precision 					:: previous_x3
double precision 					:: previous_x4
double precision 					:: previous_x5

!safeguard exit condition
integer								:: iteration;
integer, parameter					:: max_iterations = 40;

!exit condition if x converges to 10^-5
double precision, parameter			:: tolerance = .00001 

! ||k||∞
double precision					:: infinity_norm_k

! || k - k-1 ||∞
double precision					:: absolute_error_infinity_norm

! This is the resulting vector of k - k-1
double precision, dimension(5)		:: k_minus_k_minus_one_vector

!initial guess at 0,0,0
current_x1 = 2.0
current_x2 = 0.0
current_x3 = 0.66667
current_x4 = 1.83
current_x5 = -2.5

infinity_norm_k = 0.0
absolute_error_infinity_norm = 0.0
k_minus_k_minus_one_vector = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)

!initialize potential exit condition
iteration = 0

do
	if (iteration >= max_iterations) then
		write(*,*) "Reached max iterations and did not converge to tolerance.", iteration
		write(*,*) "x1: ", current_x1
		write(*,*) "x2: ", current_x2
		write(*,*) "x3: ", current_x3
		write(*,*) "x4: ", current_x4
		write(*,*) "x5: ", current_x5
		EXIT
	end if
	
	previous_x1 = current_x1
	previous_x2 = current_x2
	previous_x3 = current_x3
	previous_x4 = current_x4
	previous_x5 = current_x5
	
	current_x1 = (-previous_x2 + previous_x3 - previous_x4 + previous_x5 + 2.0d0)
	current_x2 = (4.0d0 - (2.0d0 * previous_x1) - previous_x3 + previous_x4 - previous_x5)/(2.0d0)
	current_x3 = ((-3.0d0 * previous_x1) - (previous_x2) + (2.0d0 * previous_x4) - (3.0d0 * previous_x5) + 8.0d0)/(-3.0d0)
	current_x4 = ((previous_x1 * -4.0d0) - previous_x2 + previous_x3 + (5.0d0 * previous_x5) + 16.0d0)/(4.0d0)
	current_x5 = ((-16.0d0 * previous_x1) + previous_x2 - previous_x3 + previous_x4 + 32.0d0)/(-1.0d0)

	!calculate k - k(n-1)
	k_minus_k_minus_one_vector(1) = current_x1 - previous_x1
	k_minus_k_minus_one_vector(2) = current_x2 - previous_x2
	k_minus_k_minus_one_vector(3) = current_x3 - previous_x3
	k_minus_k_minus_one_vector(4) = current_x4 - previous_x4
	k_minus_k_minus_one_vector(5) = current_x5 - previous_x5

	!calculate ||k||∞
	infinity_norm_k = max(abs(current_x1), abs(current_x2), abs(current_x3), abs(current_x4), abs(current_x5))

	!calculate ||k-k(n-1)||∞
	absolute_error_infinity_norm = max(abs(k_minus_k_minus_one_vector(1)), abs(k_minus_k_minus_one_vector(2)), &
			abs(k_minus_k_minus_one_vector(3)), abs(k_minus_k_minus_one_vector(4)), abs(k_minus_k_minus_one_vector(5)))

	!check if relative error is within tolerance
	if (absolute_error_infinity_norm <= tolerance) then
		write (*,*) "Converged in : ", iteration, " iterations"
		write(*,*) "x1: ", current_x1
		write(*,*) "x2: ", current_x2
		write(*,*) "x3: ", current_x3
		write(*,*) "x4: ", current_x4
		write(*,*) "x5: ", current_x5
		EXIT
	end if

	!don't forget to update a possible exit condition
	iteration = iteration + 1
end do


end program jacobi