SUBROUTINE init(self)
    CLASS(logger_type), INTENT(INOUT) :: self

    self%lu = lu
	self%info = 0
	self%moreinfo = 1
	self%trace = 2
	self%debug = 3
	self%warn = -1
	self%error = -2
	self%fatal = -3

    OPEN(UNIT=self%lu, FILE=self%fname, STATUS='REPLACE', ACTION='WRITE')
END SUBROUTINE init


FUNCTION lv_name(lv)
    INTEGER(1), INTENT(IN) :: lv
    CHARACTER(len=10) :: lv_name

    SELECT CASE(lv)
    CASE(0)
        lv_name = 'INFO'
    CASE(1)
        lv_name = 'INFO'
    CASE(2)
        lv_name = 'TRACE'
    CASE(3)
        lv_name = 'DEBUG'
    CASE(-1)
        lv_name = 'WARN'
    CASE(-2)
        lv_name = 'ERROR'
    CASE(-3)
        lv_name = 'FATAL'
    END SELECT 
END FUNCTION lv_name


SUBROUTINE log_real(self, lv, msg, val, addnl_val)
    CLASS(logger_type), INTENT(INOUT) :: self
    INTEGER(1), INTENT(IN) :: lv
    CHARACTER(len=*), INTENT(IN) :: msg
    CHARACTER(len=20) :: buffer
    REAL*8, INTENT(IN) :: val
    REAL*8, OPTIONAL :: addnl_val

    timenow = timefetch()-tglobal_start
    WRITE(buffer, '(f6.2)') timenow
    buffer = "["//TRIM(lv_name(lv))//"] " // "["//TRIM(buffer)//" s]"

    IF (lv<self%level .OR. lv==self%level) THEN
        IF (PRESENT(addnl_val)) THEN
            WRITE(self%lu,*) buffer, ": ", lv, msg, val, addnl_val
        ELSE
            WRITE(self%lu,*) buffer, ": ", msg, val
        ENDIF
    END IF
END SUBROUTINE log_real


SUBROUTINE log_int(self, lv, msg, val, addnl_val)
    CLASS(logger_type), INTENT(INOUT) :: self
    INTEGER(1), INTENT(IN) :: lv
    CHARACTER(len=*), INTENT(IN) :: msg
    CHARACTER(len=20) :: buffer
    INTEGER, INTENT(IN) :: val
    INTEGER, OPTIONAL :: addnl_val

    timenow = timefetch()-tglobal_start
    WRITE(buffer, '(f6.2)') timenow
    buffer = "["//TRIM(lv_name(lv))//"] " // "["//TRIM(buffer)//" s]"

    IF (lv<self%level .OR. lv==self%level) THEN
        IF (PRESENT(addnl_val)) THEN
            WRITE(self%lu,*) buffer, ": ", lv, msg, val, addnl_val
        ELSE
            WRITE(self%lu,*) buffer, ": ", msg, val
        ENDIF
    END IF
END SUBROUTINE log_int


SUBROUTINE log_str(self, lv, msg)
    CLASS(logger_type), INTENT(INOUT) :: self
    INTEGER(1), INTENT(IN) :: lv
    CHARACTER(len=*), INTENT(IN) :: msg
    CHARACTER(len=20) :: buffer

    timenow = timefetch()-tglobal_start
    WRITE(buffer, '(f6.2)') timenow
    buffer = "["//TRIM(lv_name(lv))//"] " // "["//TRIM(buffer)// " s]"

    IF (lv<self%level .OR. lv==self%level) THEN
        WRITE(self%lu,*) buffer, ": ", msg
    END IF
END SUBROUTINE log_str