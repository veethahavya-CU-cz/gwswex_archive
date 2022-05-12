subroutine init()
    character(*) :: chd_file, p_file, et_file

    chd_file = trim(input_path)//"chd.ip"
    p_file = trim(input_path)//"p.ip"
    et_file = trim(input_path)//"et.ip"

    open(tu, file=chd_file, form='unformatted')
    read(tu) chd
    close(tu, status='keep')
    open(tu, file=p_file, form='unformatted')
	read(tu) p
    close(tu, status='keep')
    open(tu, file=et_file, form='unformatted')
	read(tu) et
    close(tu, status='keep')

    write(lu,*) "initialised"
end subroutine