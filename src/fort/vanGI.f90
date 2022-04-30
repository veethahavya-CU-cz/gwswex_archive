function vanGI_simps(d, simpsons_n)
    real*8 :: vanGI_simps
    real*8, intent(inout) :: d
    integer, intent(in) :: simpsons_n
    d = -d/100
    vanGI_simps = simpsons(theta, d, dble(0), simpsons_n)*100
    contains
        function theta(h_c)
            real*8, intent(in) :: h_c
            real*8 :: theta, theta_r, theta_s, alpha, n, m
            theta_r = parms_vanGI(1)
            theta_s = parms_vanGI(2)
            alpha = parms_vanGI(3)
            n = parms_vanGI(4)
            m = (1-(1/n))
            theta = theta_r + ((theta_s - theta_r)/((1+(alpha*(abs(h_c)))**n))**m)
        end function theta
end function vanGI_simps

function vanGI_quad(d)
    USE quadpack
    real*8 :: vanGI_quad
    real*8, intent(inout) :: d
    integer ( kind = 4 ), parameter :: limit = 500
    integer ( kind = 4 ), parameter :: lenw = 4 * limit
    real ( kind = 8 ), parameter :: a = 0.0D+00
    real ( kind = 8 ) abserr
    real ( kind = 8 ) b
    real ( kind = 8 ), parameter :: epsabs = 0.0D+00
    real ( kind = 8 ), parameter :: epsrel = 0.001D+00
    real ( kind = 8 ), external :: f02
    integer ( kind = 4 ) ier
    integer ( kind = 4 ) iwork(limit)
    integer ( kind = 4 ), parameter :: key = 6
    integer ( kind = 4 ) last
    integer ( kind = 4 ) neval
    real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
    real ( kind = 8 ) result
    real ( kind = 8 ), parameter :: true = 0.06278740D+00
    real ( kind = 8 ) work(lenw)

    d = -d/100
    b = dble(0)

    CALL dqag(theta, d, b, epsabs, epsrel, key, vanGI_quad, abserr, neval, ier, limit, lenw, last, iwork, work)
    vanGI_quad = vanGI_quad*100

    contains
        function theta(h_c)
            real*8, intent(in) :: h_c
            real*8 :: theta, theta_r, theta_s, alpha, n, m
            theta_r = parms_vanGI(1)
            theta_s = parms_vanGI(2)
            alpha = parms_vanGI(3)
            n = parms_vanGI(4)
            m = (1-(1/n))
            theta = theta_r + ((theta_s - theta_r)/((1+(alpha*(abs(h_c)))**n))**m)
        end function theta
end function vanGI_quad