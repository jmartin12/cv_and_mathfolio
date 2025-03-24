program jacobian_parallel

! required MPI include file
include 'mpif.h'

! a <= x <= b
real, parameter			::	a		= 1.0
real, parameter			::	b		= 2.0

! y(a) = 1 and y(b) = 2
real, parameter			::	alpha	= 1.0
real, parameter			::	beta	= 2.0

!configurable
integer, parameter		::	n		= 9999

!step size
real					::	h

!which step we are at
real					::	x

!used throughout algorithm 11.3 to hold various intermediate approximation values
real, dimension(1:n)	::	a_vector
real, dimension(1:n)	::	b_vector
real, dimension(1:n)	::	c_vector
real, dimension(1:n)	::	d_vector

!final values that we will approximate
real, dimension(0:n+1)	::	w_approximations

!step 1
h = (b - a)/(n + 1)
write(*,*) "h size:  ", h
x = a + h
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
call jacobian_parallelized(n, a_vector, b_vector, c_vector, d_vector, w_approximations, alpha, beta, h)

contains
	subroutine jacobian_parallelized(n, a_vector, b_vector, c_vector, d_vector, w_approx, alpha, beta, h)
		implicit none
		integer,			  	intent(in)		::	n
		real, dimension(1:n), 	intent(inout)	::	a_vector
		real, dimension(1:n), 	intent(inout)	::	b_vector
		real, dimension(1:n), 	intent(inout)	::	c_vector
		real, dimension(1:n), 	intent(inout)	::	d_vector
		real, dimension(0:n+1),	intent(inout)	::	w_approx
		real,					intent(in)		::	alpha
		real,					intent(in)		::	beta
		real,					intent(in)		::	h

		!failsafe stopping conditions
		integer,	parameter					::	max_iterations 					= 1000000
		integer									::	current_iteration
		
		!convergence stopping condition
		real,		parameter					::	tolerance 						= .000001
		real									:: 	absolute_error_infinity_norm
		real									:: 	potential_max

		!used when gathering previous approx's from child processes
		real, dimension(1:n)					:: prev_w_approx

		!mpi related variables
		integer numtasks, rank, ierr  
		integer stat(MPI_STATUS_SIZE)   ! required variable for receive routines
		real	outmsg, inmsg			! for passing ghost values to each other 
		logical quit					! when master process lets the child process know to stop iterating

		!our starting left boundry, extending to the right side with an additional ghost value
		!Ex: n = 9. 
		!This will hold:  | 1 | 2 | 3 | 4 | 5 | 6 | 
		!6 is ghost value
		real, dimension(1: ((n+1)/2)+1)			::	w_approx_rank_one 
		real, dimension(1: ((n+1)/2)+1)			::	prev_rank_one_approx
		
		!an additional ghost value on the left side, extending to our ending right boundry
		!Ex: n = 9. 
		!This will hold: |5 | 6 | 7 | 8 | 9 |  
		!5 is ghost value
		real, dimension(((n+1)/2): n)			::	w_approx_rank_two
		real, dimension(((n+1)/2): n)			::	prev_rank_two_approx
		
		!loop controls,rank 0
		integer									::	l
		!loop controls,rank 1
		integer									::	j
		!loop controls,rank 2
		integer									::	k

		! initialize MPI
		call MPI_INIT(ierr)

		! get number of tasks
		call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)

		! get my rank
		call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

		if (rank == 1) then			
			!all initial guesses to 0, including ghost value
			do j = 1, ((n+1)/2)+1
				w_approx_rank_one(j) = 0 
			end do

			do 				
				do j = 1, ((n+1)/2)+1
					prev_rank_one_approx(j) = w_approx_rank_one(j)
				end do

				!be sure not to calculate the value for the ghost value
				do j = 1, (n+1)/2
					if (j == 1) then
						w_approx_rank_one(j) = (-b_vector(j) * prev_rank_one_approx(j+1) + d_vector(j)) / a_vector(j)
					else
						w_approx_rank_one(j) = (-c_vector(j) * prev_rank_one_approx(j-1) - b_vector(j) * prev_rank_one_approx(j+1) & 
						+ d_vector(j)) &
						/ a_vector(j)
					end if
				end do

				call MPI_SEND(w_approx_rank_one, ((n+1)/2)+1, MPI_REAL, 0, 1, MPI_COMM_WORLD, ierr)
				call MPI_SEND(prev_rank_one_approx, ((n+1)/2)+1, MPI_REAL, 0, 1, MPI_COMM_WORLD, ierr)
				call MPI_RECV(quit, 1, MPI_LOGICAL, 0, 1, MPI_COMM_WORLD, stat, ierr)

				if (quit) then
					write(*,*) "Exiting rank 1."
					exit
				end if

				!send MPI ghost value to other process, specifically in the (n+1)/2 index to the buffer
				outmsg = w_approx_rank_one((n+1)/2)
				call MPI_SEND(outmsg, 1, MPI_REAL, 2, 1, MPI_COMM_WORLD, ierr)
				!write(*,*) "Sent a value of ", w_approx_rank_one((n+1)/2), "from rank 1 to rank 2"
				!receive the ghost value this process will use, specifically in the ((n+1)/2 + 1) index
    			call MPI_RECV(inmsg, 1, MPI_REAL, 2, 1, MPI_COMM_WORLD, stat, ierr)
				w_approx_rank_one(((n+1)/2) + 1) = inmsg
				!write(*,*) "Received a value of ", w_approx_rank_one(((n+1)/2) + 1), "from rank 2 to rank 1" 
			end do
		end if

		if (rank == 2) then		
			!all initial guesses to 0, including ghost value
			do k = (n+1)/2, n
				w_approx_rank_two(k) = 0 
			end do

			do 
				do k = (n+1)/2, n
					prev_rank_two_approx(k) = w_approx_rank_two(k)
				end do

				!be sure not to calculate the value for the ghost value
				do k = ((n+1)/2)+1, n
					if (k == n) then
						w_approx_rank_two(k) = (-c_vector(k) * prev_rank_two_approx(k-1) + d_vector(k)) / a_vector(k)
					else
						w_approx_rank_two(k) = (-c_vector(k) * prev_rank_two_approx(k-1) - b_vector(k) * prev_rank_two_approx(k+1) &
						+ d_vector(k)) &
						/ a_vector(k)
					end if
				end do


				call MPI_SEND(w_approx_rank_two, (n+1)/2, MPI_REAL, 0, 1, MPI_COMM_WORLD, ierr)
				call MPI_SEND(prev_rank_two_approx, (n+1)/2, MPI_REAL, 0, 1, MPI_COMM_WORLD, ierr)
				call MPI_RECV(quit, 1, MPI_LOGICAL, 0, 1, MPI_COMM_WORLD, stat, ierr)

				if (quit) then
					write(*,*) "Exiting rank 2."
					exit
				end if

				!receive the ghost value this process will use, specifically in the (n+1)/2 index
    			call MPI_RECV(inmsg, 1, MPI_REAL, 1, 1, MPI_COMM_WORLD, stat, ierr)
				w_approx_rank_two((n+1)/2) = inmsg
				!write(*,*) "Received a value of ", w_approx_rank_two((n+1)/2), "from rank 1 to rank 2" 
				!send MPI values to other process, specifically in the  ((n+1)/2)+1 index to the buffer
				outmsg = w_approx_rank_two(((n+1)/2)+1)
				call MPI_SEND(outmsg, 1, MPI_REAL, 1, 1, MPI_COMM_WORLD, ierr)
				!write(*,*) "Sent a value of ", w_approx_rank_two(((n+1)/2)+1), "from rank 2 to rank 1"
			end do
		end if

		!rank 0 will be in charge of detecting if we pass tolerance or not.
		!If we do pass our tolerance, simply print the output and calculate l2 error
		!If we dont pass our tolerance, send back a message to the other processes and unblock them
		if (rank == 0) then
			current_iteration = 1
			do l = 1, n
				w_approx(l) = 0
			end do

			do
				if (current_iteration > max_iterations) then
					!Tell other processes to exit
					quit = .true.
					call MPI_SEND(quit, 1, MPI_LOGICAL, 1, 1, MPI_COMM_WORLD, ierr)
					call MPI_SEND(quit, 1, MPI_LOGICAL, 2, 1, MPI_COMM_WORLD, ierr) 
					write(*,*) "Failed to converge in jacobian iterative method."
					exit
				end if

				do l = 1, n
					prev_w_approx(l) = w_approx(l)
				end do

				!get all the new values for things
				call MPI_RECV(w_approx_rank_one, ((n+1)/2)+1, MPI_REAL, 1, 1, MPI_COMM_WORLD, stat, ierr)
				call MPI_RECV(prev_rank_one_approx, ((n+1)/2)+1, MPI_REAL, 1, 1, MPI_COMM_WORLD, stat, ierr)
				call MPI_RECV(w_approx_rank_two, (n+1)/2, MPI_REAL, 2, 1, MPI_COMM_WORLD, stat, ierr)
				call MPI_RECV(prev_rank_two_approx, (n+1)/2, MPI_REAL, 2, 1, MPI_COMM_WORLD, stat, ierr)

				!gather results
				do l = 1, n
					if (l <= (n+1)/2) then
						w_approx(l) = w_approx_rank_one(l)
					else
						w_approx(l) = w_approx_rank_two(l)
					end if
				end do

				potential_max = -1.0
				absolute_error_infinity_norm = -1.0			

				do l = 1, n
					potential_max = abs(w_approx(l) - prev_w_approx(l))

					if (potential_max > absolute_error_infinity_norm) then
						absolute_error_infinity_norm = potential_max
					end if
				end do

				if (absolute_error_infinity_norm <= tolerance) then
					write(*,*) "Converged. Number of iterations: ", current_iteration + 1
					!Don't forget to set our IC's
					w_approx(0) = alpha
					w_approx(n+1) = beta

					do l = 0, n+1
						x = a + l * h
						if (mod(l, 10) == 0) then
									write(*,*) "x:", x, &
										"w(i):", w_approximations(l), &
										"y(i):", actual_solution(x), &
										"error:", actual_solution(x) - w_approximations(l)
							end if
					end do

					write(*,*) "L^2 error is: ", l_2_error(n, a, h, w_approx)

					!Tell other processes to exit
					quit = .true.
					call MPI_SEND(quit, 1, MPI_LOGICAL, 1, 1, MPI_COMM_WORLD, ierr)
					call MPI_SEND(quit, 1, MPI_LOGICAL, 2, 1, MPI_COMM_WORLD, ierr) 

					exit
				else
					!Unblock the other processes so that we can continue iterating
					quit = .false.
					call MPI_SEND(quit, 1, MPI_LOGICAL, 1, 1, MPI_COMM_WORLD, ierr)
					call MPI_SEND(quit, 1, MPI_LOGICAL, 2, 1, MPI_COMM_WORLD, ierr) 
				end if

				current_iteration = current_iteration + 1
			end do
		end if

		call MPI_FINALIZE(ierr)
	end subroutine jacobian_parallelized

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
		
		actual_solution = (c_1() * x) + (c_2()/x**2) - (3.0/10.0) * sin(log(x)) - (1.0/10.0) * cos(log(x))
	end function actual_solution

	! c1 = 11/10 - c2 
	real function c_1()
		implicit none

		c_1 = (11.0/10.0) - c_2() 
	end function c_1

	! c2 = 1/70(8-12sin(log2))-4cos(log(2))
	real function c_2()
		implicit none

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

end program jacobian_parallel


