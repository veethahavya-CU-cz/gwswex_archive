SUBROUTINE init()
	CHARACTER(255) :: chd_file, p_file, et_file

	chd_file = TRIM(input_path)//"chd.ip"
	p_file = TRIM(input_path)//"p.ip"
	et_file = TRIM(input_path)//"et.ip"

	OPEN(tu, file=chd_file, form='unformatted', action='READ')
	READ(tu) chd
	CLOSE(tu, status='keep')
	OPEN(tu, file=p_file, form='unformatted', action='READ')
	READ(tu) p
	CLOSE(tu, status='keep')
	OPEN(tu, file=et_file, form='unformatted', action='READ')
	READ(tu) et
	CLOSE(tu, status='keep')

	CALL logger%log("initialised")
END SUBROUTINE