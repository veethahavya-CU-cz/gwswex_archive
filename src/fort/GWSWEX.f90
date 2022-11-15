MODULE GWSWEX
	USE core, only: build, init, run
	USE timing

	IMPLICIT NONE

    CONTAINS
        FUNCTION initialize() RESULT(status)
            INTEGER :: status
            status = 0
            CALL build()
            CALL init()
            status = 1
        END FUNCTION initialize

        FUNCTION finalize() RESULT(status)
            INTEGER :: status
            status = 0
            CALL run()
            status = 1
        END FUNCTION finalize
END MODULE GWSWEX