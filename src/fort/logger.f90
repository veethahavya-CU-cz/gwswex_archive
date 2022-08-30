SUBROUTINE init(self)
    CLASS(logger_type), INTENT(INOUT) :: self

    IF (self%switch) THEN
        OPEN(UNIT=self%lu, FILE=self%fname, STATUS='REPLACE', ACTION='WRITE')
    END IF

END SUBROUTINE init

SUBROUTINE log_real(self, msg, val, addnl_val)
    CLASS(logger_type), INTENT(INOUT) :: self
    CHARACTER(len=*), INTENT(IN) :: msg
    REAL*8, INTENT(IN) :: val
    REAL*8, OPTIONAL :: addnl_val

    IF (self%switch) THEN
        IF (PRESENT(addnl_val)) THEN
            WRITE(self%lu,*) msg, val, addnl_val
        ELSE
            WRITE(self%lu,*) msg, val
        ENDIF
    END IF
END SUBROUTINE log_real

SUBROUTINE log_int(self, msg, val, addnl_val)
    CLASS(logger_type), INTENT(INOUT) :: self
    CHARACTER(len=*), INTENT(IN) :: msg
    INTEGER, INTENT(IN) :: val
    INTEGER, OPTIONAL :: addnl_val

    IF (self%switch) THEN
        IF (PRESENT(addnl_val)) THEN
            WRITE(self%lu,*) msg, val, addnl_val
        ELSE
            WRITE(self%lu,*) msg, val
        ENDIF
    END IF
END SUBROUTINE log_int

SUBROUTINE log_str(self, msg)
    CLASS(logger_type), INTENT(INOUT) :: self
    CHARACTER(len=*), INTENT(IN) :: msg

    IF (self%switch) THEN
        WRITE(self%lu,*) msg
    END IF
END SUBROUTINE log_str