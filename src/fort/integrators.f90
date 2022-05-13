function rectangular(f, a, b, n)
!f2py threadsafe
	real*8 :: rectangular
	interface
		function f(x)
			real*8 :: f
			real*8, intent(in) :: x
		end function
	end interface
	real*8, intent(in) :: a, b
	integer*2, intent(in) :: n
	integer :: i
	real*8 :: x, step

	rectangular = 0.0
	step = (b-a)/(real(n)-1.0)

	!$OMP PARALLEL DO
	do i = 1, n
		x = a + (real(i)-0.5)*step
		rectangular = rectangular + step*f(x)
	end do
	!$OMP END PARALLEL DO
end function

function trapezoidal(f, a, b, n)
!f2py threadsafe
	real*8 :: trapezoidal
	interface
		function f(x)
			real*8 :: f
			real*8, intent(in) :: x
		end function
	end interface
	real*8, intent(in) :: a, b
	integer*2, intent(in) :: n
	integer :: i
	real*8 :: x1, x2, step

	trapezoidal = 0.0
	step = (b-a)/(real(n)-1.0)

	!$OMP PARALLEL DO
	do i = 1, n
		x1 = a + (real(i)-1.0)*step
		x2 = a + (real(i))*step
		trapezoidal = trapezoidal + step*(f(x1)+f(x2))/2.0
	end do
	!$OMP END PARALLEL DO
end function


function simpsons(f, a, b, n)
!f2py threadsafe
	real*8 :: simpsons
	interface
		function f(x)
			real*8 :: f
			real*8, intent(in) :: x
		end function
	end interface
	real*8, intent(in) :: a, b
	integer*2, intent(in) :: n
	integer :: i
	real*8 :: x1, x2, step

	simpsons = 0.0
	step = (b-a)/(real(n)-1.0)

	!$OMP PARALLEL DO
	do i = 1, n
		x1 = a + (real(i)-1.0)*step
		x2 = a + (real(i))*step
		simpsons = simpsons + step*(f(x1)+4.0*f((x1+x2)/2.0)+f(x2))/6.0
	end do
	!$OMP END PARALLEL DO
end function