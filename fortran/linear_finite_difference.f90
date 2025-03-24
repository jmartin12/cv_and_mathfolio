program linear_finite_difference

! a <= x <= b
real, parameter			::	a		= 1.0
real, parameter			::	b		= 2.0

! y(a) = 1 and y(b) = 2
real, parameter			::	alpha	= 1.0
real, parameter			::	beta	= 2.0

!used to control our loops 
integer					::	i

!configurable
integer, parameter		::	n		= 9.0

!step size
real					::	h

!which step we are at
real					::	x

!used throughout algorithm 11.3 to hold various intermediate approximation values
real, dimension(1:n)	::	a_vector
real, dimension(1:n)	::	b_vector
real, dimension(1:n)	::	c_vector
real, dimension(1:n)	::	d_vector

!vectors that represent the tridiagonal linear system (lower, upper, middle)
real, dimension(1:n)	::	l_diagonal
real, dimension(1:n)	::	u_diagonal
real, dimension(1:n)	::	z_diagonal

!final values that we will approximate
real, dimension(0:n+1)	::	w_approximations

!step 1
h = (b - a)/(n + 1)
write(*,*) "h size:  ", h
x = a + h
write(*,*) "x is:  ", x
a_vector(1) = 2 + h**2 * q_x(x)
b_vector(1) = -1 + (h/2) * p_x(x)
d_vector(1) = -h**2 * r_x(x) + (1 + (h/2)* p_x(x)) * alpha

!step 2
do i = 2, n-1
	x = a + i * h
	a_vector(i) = 2 + h**2 * q_x(x)
	b_vector(i) = -1 + (h/2) * p_x(x)
	c_vector(i) = -1 - (h/2) * p_x(x)
	d_vector(i) = -h**2 * r_x(x)
end do

!step 3
x = b - h
a_vector(n) = 2 + h**2 * q_x(x)
c_vector(n) = -1 - (h/2) * p_x(x)
d_vector(n) = -h**2 * r_x(x) + (1 - (h/2) * p_x(x)) * beta

!steps 4-8
!call crout_factorization(n, a_vector, b_vector, c_vector, d_vector, l_diagonal, u_diagonal, z_diagonal, w_approximations)
call jacobian(n, a_vector, b_vector, c_vector, d_vector, w_approximations, alpha, beta, h)

!step 9
do i = 0, n+1
	x = a + i * h
	write(*,*) "x:", x, &
		 "w(i):", w_approximations(i), &
		 "y(i):", actual_solution(x), &
	 	 "error:", actual_solution(x) - w_approximations(i)
end do

write(*,*) "L^2 error is: ", l_2_error(n, a, h, w_approximations)

contains
	subroutine jacobian(n, a_vector, b_vector, c_vector, d_vector, w_approx, alpha, beta, h)
		integer,			  	intent(in)		::	n
		real, dimension(1:n), 	intent(inout)	::	a_vector
		real, dimension(1:n), 	intent(inout)	::	b_vector
		real, dimension(1:n), 	intent(inout)	::	c_vector
		real, dimension(1:n), 	intent(inout)	::	d_vector
		real, dimension(0:n+1),	intent(inout)	::	w_approx
		real,					intent(in)		::	alpha
		real,					intent(in)		::	beta
		real,					intent(in)		::	h
		
		real, dimension(1:n)					::	previous_w_approx

		!failsafe stopping condition
		integer									::	current_iteration
		integer,	parameter					::	max_iterations 					= 1000
		
		!convergence stopping condition
		real									::	potential_max					!This is used to detect max values in list of values
		real									::	absolute_error_infinity_norm
		real,		parameter					::	tolerance 						= .00001

		current_iteration = 1

		!The boundaries are already given
		w_approx(0) = alpha
		w_approx(n+1) = beta

		!all initial guesses set to 0
		do i = 1, n
			w_approx(i) = 0
		end do

		do
			if (current_iteration > max_iterations) then
				write(*,*) "Reached max iterations in Jacobian and did not converge."
				exit
			end if

			do i = 1, n
				previous_w_approx(i) = w_approx(i)
			end do

			do i = 1, n
				if (i == 1) then
					w_approx(i) = (-b_vector(i) * previous_w_approx(i + 1) + d_vector(i)) / a_vector(i)
				else if (i == n) then
					w_approx(i) = (-c_vector(i) * previous_w_approx(i - 1) + d_vector(i)) / a_vector(i)
				else
					w_approx(i) = (-c_vector(i) * previous_w_approx(i - 1) - b_vector(i) * previous_w_approx(i+1) + d_vector(i)) / a_vector(i)
				end if
			end do

			!find infinity norm by taking max of abs values in the set
			potential_max = -1.0
			absolute_error_infinity_norm = -1.0
			do i = 1, n
				potential_max = abs(w_approx(i) - previous_w_approx(i))

				if (potential_max > absolute_error_infinity_norm) then
					absolute_error_infinity_norm = potential_max
				end if
			end do

			if (absolute_error_infinity_norm <= tolerance) then
				write(*,*) "Converged in ", current_iteration + 1, " iterations."
				exit
			end if
			
			current_iteration = current_iteration + 1
		end do
	end subroutine jacobian

	subroutine crout_factorization(n, a_vector, b_vector, c_vector, d_vector, l_diagonal, u_diagonal, z_diagonal, w_approx)
		integer,			  	intent(in)		::  n
		real, dimension(1:n), 	intent(in)		::	a_vector
		real, dimension(1:n), 	intent(in)		::	b_vector
		real, dimension(1:n), 	intent(in)		::	c_vector
		real, dimension(1:n), 	intent(in)		::	d_vector
		real, dimension(1:n),	intent(inout)	::	l_diagonal
		real, dimension(1:n),	intent(inout) 	::	u_diagonal
		real, dimension(1:n),	intent(inout)	::	z_diagonal
		real, dimension(0:n+1),	intent(inout)	::	w_approx

		!step 4
		l_diagonal(1) = a_vector(1)
		u_diagonal(1) = b_vector(1)/a_vector(1)
		z_diagonal(1) = d_vector(1)/l_diagonal(1)

		!step 5
		do i = 2, n-1
			l_diagonal(i) = a_vector(i) - (c_vector(i) * u_diagonal(i-1))
			u_diagonal(i) = b_vector(i)/l_diagonal(i)
			z_diagonal(i) = (d_vector(i) - (c_vector(i) * z_diagonal(i-1)))/l_diagonal(i)
		end do

		!step 6
		l_diagonal(n) = a_vector(n) - (c_vector(n) * u_diagonal(n-1))
		z_diagonal(n) = (d_vector(n) - (c_vector(n) * z_diagonal(n-1)))/l_diagonal(n)

		!step 7
		w_approx(0) = alpha
		w_approx(n+1) = beta
		w_approx(n) = z_diagonal(n)

		!step 8
		do i = n-1, 1, -1
			w_approx(i) = z_diagonal(i) - (u_diagonal(i) * w_approx(i+1))
		end do
	end subroutine crout_factorization

	!p(x)y'
	real function p_x(x)
		implicit none 
		real, intent(in)	::	x
		
		p_x = -2/x
	end function p_x


	!q(x)y
	real function q_x(x)
		implicit none
		real, intent(in)	::	x

		q_x = 2/(x**2)
	end function q_x

	!r(x)
	real function r_x(x)
		implicit none
		real, intent(in)	::	x

		r_x = (sin(log(x)))/(x**2)	
	end function r_x

	! c1 * x + c2/x^2 - 3/10 sin(log(x)) - 1/10 cos(log(x)) 
	real function actual_solution(x)
		implicit none
		real, intent(in)	::	x
		
		actual_solution = (c_1(x) * x) + (c_2(x)/x**2) - (3.0/10.0) * sin(log(x)) - (1.0/10.0) * cos(log(x))
	end function actual_solution

	! c1 = 11/10 - c2 
	real function c_1(x)
		implicit none
		real, intent(in)	::	x

		c_1 = (11.0/10.0) - c_2(x) 
	end function c_1

	! c2 = 1/70(8-12sin(log2))-4cos(log(2))
	real function c_2(x)
		implicit none
		real, intent(in)	::	x

		c_2 = (1.0/70.0) * (8.0 - 12.0 * sin(log(2.0)) - 4.0 * cos(log(2.0)))
	end function c_2

	real function l_2_error(n, a, h, y_approx)
		implicit none
		integer, 	intent(in)	::	n
		real, 		dimension(0:n+1),	intent(in)	::	y_approx
		real, 		intent(in)	::	a
		real, 		intent(in)	::	h


		real, 		dimension(1:n)	::	y_actual
		real						::	summation

		!calculate actual values
		do i = 1, n
			x = a + i * h
			y_actual(i) = actual_solution(x)
		end do

		summation = 0.0
		!The magnitude summation of errors
		! E| y_exact(xi) - y_approx(xi) |^2 for i = 1 , n (where n is end of all interior points)
		do i = 1, n
			summation = summation + (abs(y_actual(i) - y_approx(i))**2)
		end do
		
		! sqrt( 1 / n- 1 * summation)
		l_2_error = sqrt((1.0 / (n)) * summation)
	end function l_2_error

end program linear_finite_difference


