subroutine init()
	character(255) :: chd_file, p_file, et_file

	chd_file = trim(input_path)//"chd.ip"
	p_file = trim(input_path)//"p.ip"
	et_file = trim(input_path)//"et.ip"

	open(tu, file=chd_file, form='unformatted', action='read')
	read(tu) chd
	close(tu, status='keep')
	open(tu, file=p_file, form='unformatted', action='read')
	read(tu) p
	close(tu, status='keep')
	open(tu, file=et_file, form='unformatted', action='read')
	read(tu) et
	close(tu, status='keep')

	write(lu,*) "initialised"
end subroutine