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


function vanGI_rect(d, rect_n)
    real*8 :: vanGI_rect
    real*8, intent(inout) :: d
    integer, intent(in) :: rect_n
    d = -d/100
    vanGI_rect = simpsons(theta, d, dble(0), rect_n)*100
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
end function vanGI_rect


function vanGI_trap(d, trap_n)
    real*8 :: vanGI_trap
    real*8, intent(inout) :: d
    integer, intent(in) :: trap_n
    d = -d/100
    vanGI_trap = simpsons(theta, d, dble(0), trap_n)*100
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
end function vanGI_trap


function vanGI_fgsl(d)
    USE fgsl
    real*8 :: vanGI_trap
    real*8, intent(inout) :: d
    integer, intent(in) :: trap_n
    d = -d/100
    vanGI_trap = simpsons(theta, d, dble(0), trap_n)*100
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
end function vanGI_trap