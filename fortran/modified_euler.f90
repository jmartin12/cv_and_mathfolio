program midpoint_and_euler
implicit none

!This will be used to hold the output for the midpoint method, also initialize an initial guess
double precision						:: y1
!This will be used to hold the output for the modified euler method, also initialize an initial guess
double precision						:: y2

! Used to hold the current iteration.
integer									:: i

!Boundaries for the I.V.P. [ 0 <= t <= 2]
!The lower_bounds will be used as t(i), and be updated as the algorithm progresses
double precision						:: lower_bounds
double precision,		parameter		:: upper_bounds = 2.0d0

!Actual values to the IVP so that we can compare the error between methods
double precision, dimension(11)			:: actual_ivp_values

!These will hold the calculated error between the exact answers and approximate answers
double precision, dimension(11)			:: meuler_error_margin
double precision, dimension(11)			:: midpoint_error_margin

!These will hold the calculated values, mainly used for printing a table at the end
double precision, dimension(11)			:: midpoint_values
double precision, dimension(11)			:: meuler_values

!Step size of 0.2.
!This is calculated via (upper_bounds - lower_bounds) / N. 
!We could make this configurable by the user, but for this problem we have
! (2 - 0) / 10 = 0.2
double precision,		parameter		:: step_size 	= 0.2d0

!initialize non constant variables
lower_bounds 		= 0.0d0
y1 					= 0.5d0
y2 					= 0.5d0
i 					= 1

!values taken from the book
actual_ivp_values 	= 	&
		(/ 0.5, 		&
		0.8292986, 		&
		1.2140877,		&
		1.64894060,		&
		2.1272295, 		&
		2.6408591,		&
		3.1799415,		&
		3.732400,		&
		4.2834838,		&
		4.8151763,		&
		5.3054720 /)

do	
	if (lower_bounds > upper_bounds) then
		write(*,*) "Finished midpoint and euler methods."
		exit
	end if

	call midpt(lower_bounds, y1)
	midpoint_error_margin(i) = abs(y1) - abs(actual_ivp_values(i))
	midpoint_values(i) = y1

	call meuler(lower_bounds, y2)
	meuler_error_margin(i) = abs(y2) - abs(actual_ivp_values(i))
	meuler_values(i) = y2
	
	!don't forget to update our boundaries
	lower_bounds = lower_bounds + step_size
	i = i + 1
end do


!reset i and lower bounds so we can use it again to print out output
i = 1
lower_bounds = 0.0d0

print 100
100 format ('t(i)', 7x, 'Yn (actual)', 7x, 'Yn (Midpoint)', 7x, 'Midpoint Error',7x, 'Yn (MEuler)', 7x, 'MEuler Error' )

do 
	if (i > 11) then
		write(*,*) "Finished."
		exit 
	end if

	write(*, '(f7.4)', advance="no") lower_bounds
	write(*, '(11xf7.4)', advance="no") actual_ivp_values(i)
	write(*, '(11xf7.4)', advance="no") midpoint_values(i)
	write(*, '(11xf7.4)', advance="no") midpoint_error_margin(i)
	write(*, '(11xf7.4)', advance="no") meuler_values(i)
	write(*, '(11xf7.4)')  meuler_error_margin(i)

	i = i + 1
	lower_bounds = lower_bounds + step_size
end do

contains
	!unsure what to do with dt at this point.
	subroutine midpt(tnow, y1)
		implicit none
		
		double precision, 	intent(in)		::	tnow
		double precision, 	intent(inout) 	::	y1
		
		! The O.D.E. y' = y - t^2 + 1 
		! with range [0 <= t <= 2]
		! and initial value of y(0) = .5 
		! using the midpoint formula, this reduces to the following equation
		! y(i+1) = 1.22 * y(i) - 0.088 * i^2 - 0.008 * i + 0.218
		y1 = 1.22d0 * y1 - 0.0088d0 * (tnow ** 2) - 0.008d0 * tnow + 0.218d0
	end subroutine midpt

	subroutine meuler(tnow, y2)
		implicit none

		double precision, 	intent(in)		:: tnow
		double precision, 	intent(inout) 	:: y2 

		y2 = 1.22d0 * y2 - 0.0088d0 * (tnow ** 2) - 0.008d0 * tnow + 0.216d0
	end subroutine meuler

end program midpoint_and_euler