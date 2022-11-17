MODULE GWSWEX
	USE core, only: build, init, run
	USE timing
    USE logger_mod, only: logger_type

	IMPLICIT NONE

    INTEGER(4), PUBLIC :: ne, nt

    CONTAINS
        SUBROUTINE initialize(attribs, gok_l, bot_l, n_l, k_l, macropore_inf_degree_l, vanG_pars_l, chd_l, p_l, et_l)
            IMPLICIT NONE
            INTEGER(4), INTENT(IN) :: attribs(4)
            LOGICAL, INTENT(IN)  :: chd_l(:)
	        REAL(8), INTENT(IN) :: gok_l(:), bot_l(:), n_l(:), k_l(:), macropore_inf_degree_l(:), vanG_pars_l(4), p_l(:,:), et_l(:,:)

            ne = attribs(1)
            nt = attribs(2)

            CALL build(attribs, gok_l, bot_l, n_l, k_l, macropore_inf_degree_l, vanG_pars_l)
            CALL init(chd_l, p_l, et_l)
        END SUBROUTINE initialize


        SUBROUTINE finalize(gws_l, sws_l, sm_l, epv_l, gw_sm_interconnectivity_l)
            IMPLICIT NONE
            !custom type usage (not def) possible here
            !TYPE(logger_type), DIMENSION(7) :: logger_array
            REAL(8), INTENT(INOUT) :: gws_l(:,:), sws_l(:,:), sm_l(:,:), epv_l(:,:), gw_sm_interconnectivity_l(:)

            CALL run(gws_l, sws_l, sm_l, epv_l, gw_sm_interconnectivity_l(:))
        END SUBROUTINE finalize

        SUBROUTINE fetch_1d(attr, res)
            USE core, only: gw_sm_interconnectivity

            IMPLICIT NONE
            CHARACTER(*), INTENT(IN) :: attr
            REAL(8), INTENT(INOUT) :: res(:)
            
            SELECT CASE(attr)
                CASE('gw_sm_interconnectivity')
                    res = gw_sm_interconnectivity
                CASE DEFAULT
                    ERROR STOP 'unknown attribute'
            END SELECT
        END SUBROUTINE fetch_1d

        SUBROUTINE fetch_2d(attr, res)
            USE core, only: gw_dis, sw_dis, sm_dis, Qin, Qout, Qdiff
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN) :: attr
            REAL(8), INTENT(INOUT) :: res(:,:)
            
            SELECT CASE(attr)
                CASE('gw_dis')
                    res = gw_dis
                CASE('sw_dis')
                    res = sw_dis
                CASE('sm_dis')
                    res = sm_dis
                CASE('Qin')
                    res = Qin
                CASE('Qout')
                    res = Qout
                CASE('Qdiff')
                    res = Qdiff
                CASE DEFAULT
                    ERROR STOP 'ERROR: unknown attribute'
            END SELECT

        END SUBROUTINE fetch_2d
END MODULE GWSWEX