subroutine init(chd_l, p_l, et_l)
    real*8, intent(in) :: p_l(:,:), et_l(:,:)
    logical, intent(in) :: chd_l(:)
    chd = chd_l
    p = p_l
    et = et_l
    open(unit=42, file="fort.log", status="old")
    write(42,*) "initialised"
end subroutine