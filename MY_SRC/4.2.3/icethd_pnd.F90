MODULE icethd_pnd
   !!======================================================================
   !!                     ***  MODULE  icethd_pnd   ***
   !!   sea-ice: Melt ponds on top of sea ice
   !!======================================================================
   !! history :       !  2012     (O. Lecomte)       Adaptation from Flocco and Turner
   !!                 !  2017     (M. Vancoppenolle, O. Lecomte, C. Rousset) Implementation
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3' :                                     SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_thd_pnd_init : some initialization and namelist read
   !!   ice_thd_pnd      : main calling routine
   !!----------------------------------------------------------------------
   USE phycst         ! physical constants
   USE dom_oce        ! ocean space and time domain
   USE ice            ! sea-ice: variables
   USE ice1D          ! sea-ice: thermodynamics variables
   USE icetab         ! sea-ice: 1D <==> 2D transformation
   USE sbc_ice        ! surface energy budget
   USE sbc_oce , ONLY : tprecip, sprecip
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_thd_pnd_init    ! routine called by icestp.F90
   PUBLIC   ice_thd_pnd         ! routine called by icestp.F90

   INTEGER ::              nice_pnd    ! choice of the type of pond scheme
   !                                   ! associated indices:
   INTEGER, PARAMETER ::   np_pndNO   = 0   ! No pond scheme
   INTEGER, PARAMETER ::   np_pndCST  = 1   ! Constant ice pond scheme
   INTEGER, PARAMETER ::   np_pndLEV  = 2   ! Level ice pond scheme
   INTEGER, PARAMETER ::   np_pndTOPO = 3   ! Level ice pond scheme

   !--------------------------------------------------------------------------
   ! Diagnostics for pond volume per area
   !
   ! dV/dt = mlt + drn + lid + rnf
   ! mlt   = input from surface melting
   ! drn   = drainage through brine network
   ! lid   = lid growth & melt
   ! rnf   = runoff (water directly removed out of surface melting + overflow)
   !
   ! In topo mode, the pond water lost because it is in the snow is not included in the budget
   ! In level mode, all terms are incorporated
   !
   LOGICAL ::   ll_diag_pnd
   !
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   diag_dvpn_mlt       ! meltwater pond volume input      [kg/m2/s]
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   diag_dvpn_drn       ! pond volume lost by drainage     [-]
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   diag_dvpn_lid       ! exchange with lid / refreezing   [-]
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   diag_dvpn_rnf       ! meltwater pond lost to runoff    [-]
   REAL(wp), ALLOCATABLE, DIMENSION(:)   ::   diag_dvpn_mlt_1d    ! meltwater pond volume input      [-]
   REAL(wp), ALLOCATABLE, DIMENSION(:)   ::   diag_dvpn_drn_1d    ! pond volume lost by drainage     [-]
   REAL(wp), ALLOCATABLE, DIMENSION(:)   ::   diag_dvpn_lid_1d    ! exchange with lid / refreezing   [-]
   REAL(wp), ALLOCATABLE, DIMENSION(:)   ::   diag_dvpn_rnf_1d    ! meltwater pond lost to runoff    [-]

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION ice_thd_pnd_alloc()
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_thd_pnd_alloc ***
      !!-------------------------------------------------------------------
      ALLOCATE( diag_dvpn_mlt   (jpi,jpj), diag_dvpn_lid   (jpi,jpj), diag_dvpn_drn   (jpi,jpj), diag_dvpn_rnf(jpi,jpj) , &
         &      diag_dvpn_mlt_1d(jpij)  , diag_dvpn_lid_1d(jpij)  , diag_dvpn_drn_1d(jpij)  , diag_dvpn_rnf_1d(jpij), &
         &      STAT=ice_thd_pnd_alloc )
      
      CALL mpp_sum ( 'icethd_pnd', ice_thd_pnd_alloc )
      IF( ice_thd_pnd_alloc /= 0 )   CALL ctl_stop( 'STOP',  'ice_thd_pnd_alloc: failed to allocate arrays'  )
      !
   END FUNCTION ice_thd_pnd_alloc

   
   SUBROUTINE ice_thd_pnd

      !!-------------------------------------------------------------------
      !!               ***  ROUTINE ice_thd_pnd   ***
      !!
      !! ** Purpose :   change melt pond fraction and thickness
      !!
      !! ** Note    : Melt ponds affect only radiative transfer for now
      !!              No heat, no salt.
      !!              The current diagnostics lacks a contribution from drainage
      !!-------------------------------------------------------------------
      INTEGER ::   ji, jj, jl        ! loop indices
      !!-------------------------------------------------------------------

      ALLOCATE( diag_dvpn_mlt(jpi,jpj), diag_dvpn_lid(jpi,jpj), diag_dvpn_drn(jpi,jpj), diag_dvpn_rnf(jpi,jpj) )
      ALLOCATE( diag_dvpn_mlt_1d(jpij), diag_dvpn_lid_1d(jpij), diag_dvpn_drn_1d(jpij), diag_dvpn_rnf_1d(jpij) )
      
         diag_dvpn_mlt (:,:) = 0._wp   ;   diag_dvpn_drn (:,:) = 0._wp
         diag_dvpn_lid (:,:) = 0._wp   ;   diag_dvpn_rnf (:,:) = 0._wp
         diag_dvpn_mlt_1d(:) = 0._wp   ;   diag_dvpn_drn_1d(:) = 0._wp
         diag_dvpn_lid_1d(:) = 0._wp   ;   diag_dvpn_rnf_1d(:) = 0._wp
      
      !-------------------------------------
      !  Remove ponds where ice has vanished
      !-------------------------------------
      at_i(:,:) = SUM( a_i, dim=3 )
      !
      DO jl = 1, jpl
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            IF( v_i(ji,jj,jl) < epsi10 .OR. at_i(ji,jj) < epsi10 ) THEN
               wfx_pnd  (ji,jj)    = wfx_pnd(ji,jj) + ( v_ip(ji,jj,jl) + v_il(ji,jj,jl) ) * rhow * r1_Dt_ice
               a_ip     (ji,jj,jl) = 0._wp
               v_ip     (ji,jj,jl) = 0._wp
               v_il     (ji,jj,jl) = 0._wp
               h_ip     (ji,jj,jl) = 0._wp
               h_il     (ji,jj,jl) = 0._wp
               a_ip_frac(ji,jj,jl) = 0._wp
            ENDIF
         END_2D
      END DO

      !------------------------------
      !  Identify grid cells with ice
      !------------------------------
      npti = 0   ;   nptidx(:) = 0
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         IF( at_i(ji,jj) >= epsi10 ) THEN
            npti = npti + 1
            nptidx( npti ) = (jj - 1) * jpi + ji
         ENDIF
      END_2D

      !------------------------------------
      !  Select melt pond scheme to be used
      !------------------------------------
      IF( npti > 0 ) THEN
         SELECT CASE ( nice_pnd )
            !         CASE (np_pndCST)   ;   CALL pnd_CST     !== Constant melt ponds  ==!
            !
         CASE (np_pndLEV)   ;   CALL pnd_LEV                  !== Level ice melt ponds  ==!
            !
         CASE (np_pndTOPO)  ;   CALL pnd_TOPO                 !== Topographic melt ponds  ==!
         END SELECT
      ENDIF
      
      ! the following fields need to be updated in the halos (done in icethd): a_ip, v_ip, v_il

      !------------------------------------
      !  Diagnostics
      !------------------------------------
       CALL iom_put( 'dvpn_mlt', diag_dvpn_mlt ) ! input from melting
       CALL iom_put( 'dvpn_lid', diag_dvpn_lid ) ! exchanges with lid
       CALL iom_put( 'dvpn_drn', diag_dvpn_drn ) ! vertical drainage
       CALL iom_put( 'dvpn_rnf', diag_dvpn_rnf ) ! runoff + overflow
   
      DEALLOCATE( diag_dvpn_mlt   , diag_dvpn_lid   , diag_dvpn_drn   , diag_dvpn_rnf    )
      DEALLOCATE( diag_dvpn_mlt_1d, diag_dvpn_lid_1d, diag_dvpn_drn_1d, diag_dvpn_rnf_1d )
      
   END SUBROUTINE ice_thd_pnd


   SUBROUTINE pnd_CST
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE pnd_CST  ***
      !!
      !! ** Purpose :   Compute melt pond evolution
      !!
      !! ** Method  :   Melt pond fraction and thickness are prescribed
      !!                to non-zero values when t_su = 0C
      !!
      !! ** Tunable parameters : pond fraction (rn_apnd), pond depth (rn_hpnd)
      !!
      !! ** Note   : Coupling with such melt ponds is only radiative
      !!             Advection, ridging, rafting... are bypassed
      !!
      !! ** References : Bush, G.W., and Trump, D.J. (2017)
      !!-------------------------------------------------------------------
      INTEGER  ::   ii, jl    ! loop indices
      REAL(wp) ::   zdv_pnd   ! Amount of water going into the ponds & lids
      !!-------------------------------------------------------------------
      DO jl = 1, jpl

         CALL tab_2d_1d( npti, nptidx(1:npti), a_i_1d    (1:npti), a_i    (:,:,jl) )
         CALL tab_2d_1d( npti, nptidx(1:npti), t_su_1d   (1:npti), t_su   (:,:,jl) )
         CALL tab_2d_1d( npti, nptidx(1:npti), a_ip_1d   (1:npti), a_ip   (:,:,jl) )
         CALL tab_2d_1d( npti, nptidx(1:npti), h_ip_1d   (1:npti), h_ip   (:,:,jl) )
         CALL tab_2d_1d( npti, nptidx(1:npti), h_il_1d   (1:npti), h_il   (:,:,jl) )
         CALL tab_2d_1d( npti, nptidx(1:npti), wfx_pnd_1d(1:npti), wfx_pnd(:,:)    )

         DO ii = 1, npti
            !
            zdv_pnd = ( h_ip_1d(ii) + h_il_1d(ii) ) * a_ip_1d(ii)
            !
            IF( a_i_1d(ii) >= 0.01_wp .AND. t_su_1d(ii) >= rt0 ) THEN
               h_ip_1d(ii)      = rn_hpnd
               a_ip_1d(ii)      = rn_apnd * a_i_1d(ii)
               h_il_1d(ii)      = 0._wp    ! no pond lids whatsoever
            ELSE
               h_ip_1d(ii)      = 0._wp
               a_ip_1d(ii)      = 0._wp
               h_il_1d(ii)      = 0._wp
            ENDIF
            !
            v_ip_1d(ii) = h_ip_1d(ii) * a_ip_1d(ii)
            v_il_1d(ii) = h_il_1d(ii) * a_ip_1d(ii)
            !
            zdv_pnd = ( h_ip_1d(ii) + h_il_1d(ii) ) * a_ip_1d(ii) - zdv_pnd
            wfx_pnd_1d(ii) = wfx_pnd_1d(ii) - zdv_pnd * rhow * r1_Dt_ice
            !
         END DO

         CALL tab_1d_2d( npti, nptidx(1:npti), a_ip_1d   (1:npti), a_ip   (:,:,jl) )
         CALL tab_1d_2d( npti, nptidx(1:npti), h_ip_1d   (1:npti), h_ip   (:,:,jl) )
         CALL tab_1d_2d( npti, nptidx(1:npti), h_il_1d   (1:npti), h_il   (:,:,jl) )
         CALL tab_1d_2d( npti, nptidx(1:npti), v_ip_1d   (1:npti), v_ip   (:,:,jl) )
         CALL tab_1d_2d( npti, nptidx(1:npti), v_il_1d   (1:npti), v_il   (:,:,jl) )
         CALL tab_1d_2d( npti, nptidx(1:npti), wfx_pnd_1d(1:npti), wfx_pnd(:,:)    )

      END DO
      !
   END SUBROUTINE pnd_CST


   SUBROUTINE pnd_LEV
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE pnd_LEV  ***
      !!
      !! ** Purpose : Compute melt pond evolution
      !!
      !! ** Method  : A fraction of meltwater is accumulated in ponds and sent to ocean when surface is freezing
      !!              We  work with volumes and then redistribute changes into thickness and concentration
      !!              assuming linear relationship between the two.
      !!
      !! ** Action  : - pond growth:      Vp = Vp + dVmelt                                          --- from Holland et al 2012 ---
      !!                                     dVmelt = (1-r)/rhow * ( rhoi*dh_i + rhos*dh_s ) * a_i
      !!                                        dh_i  = meltwater from ice surface melt
      !!                                        dh_s  = meltwater from snow melt
      !!                                        (1-r) = fraction of melt water that is not flushed
      !!
      !!              - limtations:       a_ip must not exceed (1-r)*a_i
      !!                                  h_ip must not exceed 0.5*h_i
      !!
      !!              - pond shrinking:
      !!                       if lids:   Vp = Vp -dH * a_ip
      !!                                     dH = lid thickness change. Retrieved from this eq.:    --- from Flocco et al 2010 ---
      !!
      !!                                                                   rhoi * Lf * dH/dt = ki * MAX(Tp-Tsu,0) / H
      !!                                                                      H = lid thickness
      !!                                                                      Lf = latent heat of fusion
      !!                                                                      Tp = -0.15C
      !!
      !!                                                                And solved implicitely as:
      !!                                                                   H(t+dt)**2 -H(t) * H(t+dt) -ki * (Tp-Tsu) * dt / (rhoi*Lf) = 0
      !!
      !!                    if no lids:   Vp = Vp * exp(0.01*MAX(Tp-Tsu,0)/Tp)                      --- from Holland et al 2012 ---
      !!
      !!              - Flushing:         w = -perm/visc * rho_oce * grav * Hp / Hi * flush         --- from Flocco et al 2007 ---
      !!                                     perm = permability of sea-ice                              + correction from Hunke et al 2012 (flush)
      !!                                     visc = water viscosity
      !!                                     Hp   = height of top of the pond above sea-level
      !!                                     Hi   = ice thickness thru which there is flushing
      !!                                     flush= correction otherwise flushing is excessive
      !!
      !!              - Corrections:      remove melt ponds when lid thickness is 10 times the pond thickness
      !!
      !!              - pond thickness and area is retrieved from pond volume assuming a linear relationship between h_ip and a_ip:
      !!                                  a_ip/a_i = a_ip_frac = h_ip / zaspect
      !!
      !! ** Tunable parameters : rn_apnd_max, rn_apnd_min, rn_pnd_flush
      !!
      !! ** Note       :   Mostly stolen from CICE but not only. These are between level-ice ponds and CESM ponds.
      !!
      !! ** References :   Flocco and Feltham (JGR, 2007)
      !!                   Flocco et al       (JGR, 2010)
      !!                   Holland et al      (J. Clim, 2012)
      !!                   Hunke et al        (OM 2012)
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(nlay_i) ::   ztmp           ! temporary array
      !!
      REAL(wp), PARAMETER ::   zaspect =  0.8_wp      ! pond aspect ratio
      REAL(wp), PARAMETER ::   zTp     = -0.15_wp     ! reference temperature
      REAL(wp), PARAMETER ::   zvisc   =  1.79e-3_wp  ! water viscosity
      !!
      REAL(wp) ::   zfr_mlt, zdv_mlt, zdv_avail       ! fraction and volume of available meltwater retained for melt ponding
      REAL(wp) ::   zdv_frz, zdv_flush                ! Amount of melt pond that freezes, flushes
      REAL(wp) ::   zdv_pnd                           ! Amount of water going into the ponds & lids
      REAL(wp) ::   zhp                               ! heigh of top of pond lid wrt ssh
      REAL(wp) ::   zv_ip_max                         ! max pond volume allowed
      REAL(wp) ::   zdT                               ! zTp-t_su
      REAL(wp) ::   zsbr, ztmelts                     ! Brine salinity
      REAL(wp) ::   zdum, zperm                       ! permeability of sea ice
      REAL(wp) ::   z1_rhow, z1_aspect, z1_Tp         ! inverse
      !!
      INTEGER  ::   ii, jk, jl                        ! loop indices
      !!-------------------------------------------------------------------
      z1_rhow   = 1._wp / rhow
      z1_aspect = 1._wp / zaspect
      z1_Tp     = 1._wp / zTp

      CALL tab_2d_1d( npti, nptidx(1:npti), at_i_1d          (1:npti), at_i          )
      CALL tab_2d_1d( npti, nptidx(1:npti), wfx_pnd_1d       (1:npti), wfx_pnd       )
      IF( ll_diag_pnd ) THEN
         CALL tab_2d_1d( npti, nptidx(1:npti), diag_dvpn_mlt_1d (1:npti), diag_dvpn_mlt )
         CALL tab_2d_1d( npti, nptidx(1:npti), diag_dvpn_drn_1d (1:npti), diag_dvpn_drn )
         CALL tab_2d_1d( npti, nptidx(1:npti), diag_dvpn_lid_1d (1:npti), diag_dvpn_lid )
         CALL tab_2d_1d( npti, nptidx(1:npti), diag_dvpn_rnf_1d (1:npti), diag_dvpn_rnf )
      ENDIF
      DO jl = 1, jpl

         CALL tab_2d_1d( npti, nptidx(1:npti), a_i_1d  (1:npti), a_i (:,:,jl) )
         CALL tab_2d_1d( npti, nptidx(1:npti), h_i_1d  (1:npti), h_i (:,:,jl) )
         CALL tab_2d_1d( npti, nptidx(1:npti), t_su_1d (1:npti), t_su(:,:,jl) )
         CALL tab_2d_1d( npti, nptidx(1:npti), a_ip_1d (1:npti), a_ip(:,:,jl) )
         CALL tab_2d_1d( npti, nptidx(1:npti), h_ip_1d (1:npti), h_ip(:,:,jl) )
         CALL tab_2d_1d( npti, nptidx(1:npti), h_il_1d (1:npti), h_il(:,:,jl) )

         CALL tab_2d_1d( npti, nptidx(1:npti), dh_i_sum(1:npti), dh_i_sum_2d(:,:,jl) )
         CALL tab_2d_1d( npti, nptidx(1:npti), dh_s_mlt(1:npti), dh_s_mlt_2d(:,:,jl) )

         DO jk = 1, nlay_i
            CALL tab_2d_1d( npti, nptidx(1:npti), sz_i_1d(1:npti,jk), sz_i(:,:,jk,jl) )
            CALL tab_2d_1d( npti, nptidx(1:npti), t_i_1d (1:npti,jk), t_i (:,:,jk,jl) )
         END DO

         !-----------------------
         ! Melt pond calculations
         !-----------------------
         DO ii = 1, npti
            !
            zdv_pnd = ( h_ip_1d(ii) + h_il_1d(ii) ) * a_ip_1d(ii)
            !                                                            !----------------------------------------------------!
            IF( h_i_1d(ii) < rn_himin .OR. a_i_1d(ii) < 0.01_wp ) THEN   ! Case ice thickness < rn_himin or tiny ice fraction !
               !                                                         !----------------------------------------------------!
               !--- Remove ponds on thin ice or tiny ice fractions
               a_ip_1d(ii) = 0._wp
               h_ip_1d(ii) = 0._wp
               h_il_1d(ii) = 0._wp
               !                                                         !--------------------------------!
            ELSE                                                         ! Case ice thickness >= rn_himin !
               !                                                         !--------------------------------!
               v_ip_1d(ii) = h_ip_1d(ii) * a_ip_1d(ii)   ! retrieve volume from thickness
               v_il_1d(ii) = h_il_1d(ii) * a_ip_1d(ii)
               !
               !------------------!
               ! case ice melting !
               !------------------!
               !
               !--- available meltwater for melt ponding (zdv_avail) ---!
               zdv_avail = -( dh_i_sum(ii)*rhoi + dh_s_mlt(ii)*rhos ) * z1_rhow * a_i_1d(ii) ! > 0
               zfr_mlt   = rn_apnd_min + ( rn_apnd_max - rn_apnd_min ) * at_i_1d(ii) !  = ( 1 - r ) = fraction of melt water that is not flushed
               zdv_mlt   = MAX( 0._wp, zfr_mlt * zdv_avail ) ! max for roundoff errors?
               !
               !--- overflow ---!
               !
               ! area driven overflow
               !    If pond area exceeds zfr_mlt * a_i_1d(ii) then reduce the pond volume
               !       a_ip_max = zfr_mlt * a_i
               !       => from zaspect = h_ip / (a_ip / a_i), set v_ip_max as:
               zv_ip_max = zfr_mlt * zfr_mlt * a_i_1d(ii) * zaspect
               zdv_mlt   = MAX( 0._wp, MIN( zdv_mlt, zv_ip_max - v_ip_1d(ii) ) )

               ! depth driven overflow
               !    If pond depth exceeds half the ice thickness then reduce the pond volume
               !       h_ip_max = 0.5 * h_i
               !       => from zaspect = h_ip / (a_ip / a_i), set v_ip_max as:
               zv_ip_max = z1_aspect * a_i_1d(ii) * 0.25 * h_i_1d(ii) * h_i_1d(ii)
               zdv_mlt   = MAX( 0._wp, MIN( zdv_mlt, zv_ip_max - v_ip_1d(ii) ) )

               !--- Pond growing ---!
               v_ip_1d(ii) = v_ip_1d(ii) + zdv_mlt
               !
               !--- Lid melting ---!
               IF( ln_pnd_lids )   v_il_1d(ii) = MAX( 0._wp, v_il_1d(ii) - zdv_mlt ) ! must be bounded by 0
               !
               !-------------------!
               ! case ice freezing ! i.e. t_su_1d(ii) < (zTp+rt0)
               !-------------------!
               !
               zdT = MAX( zTp+rt0 - t_su_1d(ii), 0._wp )
               !
               !--- Pond contraction (due to refreezing) ---!
               IF( ln_pnd_lids ) THEN
                  !
                  !--- Lid growing and subsequent pond shrinking ---!
                  zdv_frz = - 0.5_wp * MAX( 0._wp, -v_il_1d(ii) + & ! Flocco 2010 (eq. 5) solved implicitly as aH**2 + bH + c = 0
                     & SQRT( v_il_1d(ii)*v_il_1d(ii) + a_ip_1d(ii)*a_ip_1d(ii) * 4._wp * rcnd_i * zdT * rDt_ice / (rLfus * rhow) ) ) ! max for roundoff errors

                  ! Lid growing
                  v_il_1d(ii) = MAX( 0._wp, v_il_1d(ii) - zdv_frz )

                  ! Pond shrinking
                  v_ip_1d(ii) = MAX( 0._wp, v_ip_1d(ii) + zdv_frz )

               ELSE
                  zdv_frz = v_ip_1d(ii) * ( EXP( 0.01_wp * zdT * z1_Tp ) - 1._wp )  ! Holland 2012 (eq. 6)
                  ! Pond shrinking
                  v_ip_1d(ii) = MAX( 0._wp, v_ip_1d(ii) + zdv_frz )
               ENDIF
               !
               !--- Set new pond area and depth ---! assuming linear relation between h_ip and a_ip_frac
               ! v_ip     = h_ip * a_ip
               ! a_ip/a_i = a_ip_frac = h_ip / zaspect (cf Holland 2012, fitting SHEBA so that knowing v_ip we can distribute it to a_ip and h_ip)
               a_ip_1d(ii) = MIN( a_i_1d(ii), SQRT( v_ip_1d(ii) * z1_aspect * a_i_1d(ii) ) ) ! make sure a_ip < a_i
               h_ip_1d(ii) = zaspect * a_ip_1d(ii) / a_i_1d(ii)
               !

               !------------------------------------------------!
               ! Pond drainage through brine network (flushing) !
               !------------------------------------------------!
               ! height of top of the pond above sea-level
               zhp = ( h_i_1d(ii) * ( rho0 - rhoi ) + h_ip_1d(ii) * ( rho0 - rhow * a_ip_1d(ii) / a_i_1d(ii) ) ) * r1_rho0

               ! Calculate the permeability of the ice (Assur 1958, see Flocco
               ! 2010)
               DO jk = 1, nlay_i
                  ! MV Assur is inconsistent with SI3
                  !!zsbr = - 1.2_wp                                  &
                  !!   &   - 21.8_wp    * ( t_i_1d(ji,jk) - rt0 )    &
                  !!   &   - 0.919_wp   * ( t_i_1d(ji,jk) - rt0 )**2 &
                  !!   &   - 0.0178_wp  * ( t_i_1d(ji,jk) - rt0 )**3
                  !!ztmp(jk) = sz_i_1d(ji,jk) / zsbr
                  ! MV linear expression more consistent & simpler: zsbr = - (
                  ! t_i_1d(ji,jk) - rt0 ) / rTmlt
                  ztmelts  = -rTmlt * sz_i_1d(ii,jk)
                  ztmp(jk) = ztmelts / MIN( ztmelts, t_i_1d(ii,jk) - rt0 )
               END DO
               zperm = MAX( 0._wp, 3.e-08_wp * MINVAL(ztmp)**3 )

               ! Do the drainage using Darcy's law
               zdv_flush   = -zperm * rho0 * grav * zhp * rDt_ice / (zvisc * h_i_1d(ii)) * a_ip_1d(ii) * rn_pnd_flush 
               ! zflush comes from Hunke et al. (2012)
               zdv_flush   = MAX( zdv_flush, -v_ip_1d(ii) ) ! < 0
               v_ip_1d(ii) = v_ip_1d(ii) + zdv_flush


               !--- Set new pond area and depth ---! assuming linear relation between h_ip and a_ip_frac
               a_ip_1d(ii) = MIN( a_i_1d(ii), SQRT( v_ip_1d(ii) * z1_aspect * a_i_1d(ii) ) ) ! make sure a_ip < a_i
               h_ip_1d(ii) = zaspect * a_ip_1d(ii) / a_i_1d(ii)

               !--- Corrections and lid thickness ---!
               IF( ln_pnd_lids ) THEN
                  !--- retrieve lid thickness from volume ---!
                  IF( a_ip_1d(ii) > 0.01_wp ) THEN   ;   h_il_1d(ii) = v_il_1d(ii) / a_ip_1d(ii)
                  ELSE                               ;   h_il_1d(ii) = 0._wp
                  ENDIF
                  !--- remove ponds if lids are much larger than ponds ---!
                  IF ( h_il_1d(ii) > h_ip_1d(ii) * 10._wp ) THEN
                     a_ip_1d(ii) = 0._wp
                     h_ip_1d(ii) = 0._wp
                     h_il_1d(ii) = 0._wp
                  ENDIF
               ENDIF

               ! diagnostics: dvpnd = mlt+rnf+lid+drn
                  diag_dvpn_mlt_1d(ii) = diag_dvpn_mlt_1d(ii) + rhow *   zdv_avail             * r1_Dt_ice   ! > 0, surface melt input
                  diag_dvpn_rnf_1d(ii) = diag_dvpn_rnf_1d(ii) + rhow * ( zdv_mlt - zdv_avail ) * r1_Dt_ice   ! < 0, runoff
                  diag_dvpn_lid_1d(ii) = diag_dvpn_lid_1d(ii) + rhow *   zdv_frz               * r1_Dt_ice   ! < 0, shrinking
                  diag_dvpn_drn_1d(ii) = diag_dvpn_drn_1d(ii) + rhow *   zdv_flush             * r1_Dt_ice   ! < 0, drainage
               !
            ENDIF
            !
            v_ip_1d(ii) = h_ip_1d(ii) * a_ip_1d(ii)
            v_il_1d(ii) = h_il_1d(ii) * a_ip_1d(ii)
            !
            zdv_pnd = ( h_ip_1d(ii) + h_il_1d(ii) ) * a_ip_1d(ii) - zdv_pnd
            wfx_pnd_1d(ii) = wfx_pnd_1d(ii) - zdv_pnd * rhow * r1_Dt_ice
            !
         END DO

         !--------------------------------------------------------------------
         ! Retrieve 2D arrays
         !--------------------------------------------------------------------
         CALL tab_1d_2d( npti, nptidx(1:npti), a_ip_1d(1:npti), a_ip(:,:,jl) )
         CALL tab_1d_2d( npti, nptidx(1:npti), h_ip_1d(1:npti), h_ip(:,:,jl) )
         CALL tab_1d_2d( npti, nptidx(1:npti), h_il_1d(1:npti), h_il(:,:,jl) )
         CALL tab_1d_2d( npti, nptidx(1:npti), v_ip_1d(1:npti), v_ip(:,:,jl) )
         CALL tab_1d_2d( npti, nptidx(1:npti), v_il_1d(1:npti), v_il(:,:,jl) )
         !
      END DO
      !
      CALL tab_1d_2d( npti, nptidx(1:npti), wfx_pnd_1d(1:npti), wfx_pnd )
      !
         CALL tab_1d_2d( npti, nptidx(1:npti), diag_dvpn_mlt_1d (1:npti), diag_dvpn_mlt )
         CALL tab_1d_2d( npti, nptidx(1:npti), diag_dvpn_drn_1d (1:npti), diag_dvpn_drn )
         CALL tab_1d_2d( npti, nptidx(1:npti), diag_dvpn_lid_1d (1:npti), diag_dvpn_lid )
         CALL tab_1d_2d( npti, nptidx(1:npti), diag_dvpn_rnf_1d (1:npti), diag_dvpn_rnf )
      !
   END SUBROUTINE pnd_LEV



   SUBROUTINE pnd_TOPO

      !!-------------------------------------------------------------------
      !!                ***  ROUTINE pnd_TOPO  ***
      !!
      !! ** Purpose :   Compute melt pond evolution based on the ice
      !!                topography inferred from the ice thickness distribution
      !!
      !! ** Method  :   This code is initially based on Flocco and Feltham
      !!                (2007) and Flocco et al. (2010).
      !!
      !!                - Calculate available pond water base on surface meltwater
      !!                - Redistribute water as a function of topography, drain water
      !!                - Exchange water with the lid
      !!
      !! ** Tunable parameters :
      !!
      !! ** Note :
      !!
      !! ** References
      !!
      !!    Flocco, D. and D. L. Feltham, 2007.  A continuum model of melt pond
      !!    evolution on Arctic sea ice.  J. Geophys. Res. 112, C08016, doi:
      !!    10.1029/2006JC003836.
      !!
      !!    Flocco, D., D. L. Feltham and A. K. Turner, 2010.  Incorporation of
      !!    a physically based melt pond scheme into the sea ice component of a
      !!    climate model.  J. Geophys. Res. 115, C08012,
      !!    doi: 10.1029/2009JC005568.
      !!
      !!-------------------------------------------------------------------
      REAL(wp), PARAMETER :: &   ! shared parameters for topographic melt ponds
         zTd     = 273._wp       , & ! temperature difference for freeze-up (K)
         zvp_min = 1.e-4_wp          ! minimum pond volume (m)


      ! local variables
      REAL(wp) :: &
         zdHui,   &      ! change in thickness of ice lid (m)
         zomega,  &      ! conduction
         zdTice,  &      ! temperature difference across ice lid (C)
         zdvice,  &      ! change in ice volume (m)
         zfsurf,  &      ! net heat flux, excluding conduction and transmitted radiation (W/m2)
         zTp,     &      ! pond freezing temperature (C)
         zrhoi_L, &      ! volumetric latent heat of sea ice (J/m^3)
         zfr_mlt, &      ! fraction and volume of available meltwater retained for melt ponding
         z1_rhow, &      ! inverse water density
         zv_pnd  , &     ! volume of meltwater contributing to ponds
         zv_mlt          ! total amount of meltwater produced

      REAL(wp), DIMENSION(jpi,jpj) ::   zvolp_ini , &   !! total melt pond water available before redistribution and drainage
                                       zvolp           !! total melt pond water volume

      REAL(wp), DIMENSION(jpi,jpj) ::   zdv_pnd1, zdv_pnd2

      INTEGER  ::   ji, jj, jk, jl                    ! loop indices

      ! Note
      ! equivalent for CICE translation
      ! a_ip      -> apond
      ! a_ip_frac -> apnd

      !---------------------------------------------------------------
      ! Initialise
      !---------------------------------------------------------------

      ! Parameters & constants (move to parameters)
      zrhoi_L   = rhoi * rLfus      ! volumetric latent heat (J/m^3)
      zTp       = rt0 - 0.15_wp          ! pond freezing point, slightly below 0C (ponds are bid saline)
      z1_rhow   = 1._wp / rhow

      at_i (:,:) = SUM( a_i (:,:,:), dim=3 ) ! ice fraction
      vt_i (:,:) = SUM( v_i (:,:,:), dim=3 ) ! volume per grid area

      !---------------------------------------------------------------
      ! Change 2D to 1D
      !---------------------------------------------------------------
      ! MV
      ! a less computing-intensive version would have 2D-1D passage here
      ! use what we have in iceitd.F90 (incremental remapping)

      !--------------------------------------------------------------
      ! Collect total available pond water volume
      !--------------------------------------------------------------
      ! Assuming that meltwater (+rain in principle) runsoff the surface
      ! Holland et al (2012) suggest that the fraction of runoff decreases with total ice fraction
      ! I cite her words, they are very talkative
      ! "grid cells with very little ice cover (and hence more open water area)
      ! have a higher runoff fraction to rep- resent the greater proximity of ice to open water."
      ! "This results in the same runoff fraction r for each ice category within a grid cell"

      !clem: conservation (quick fix)
      zdv_pnd1(:,:) = 0._wp
      DO jl = 1, jpl
         DO_2D( 0, 0, 0, 0 )
           zdv_pnd1(ji,jj) = zdv_pnd1(ji,jj) + v_ip(ji,jj,jl) + v_il(ji,jj,jl)
         END_2D
      ENDDO

      zvolp(:,:) = 0._wp

      DO jl = 1, jpl
         DO_2D( 0, 0, 0, 0 )
            !                                                                                !----------------------------------------------------!
            IF( v_i(ji,jj,jl) < rn_himin*a_i(ji,jj,jl) .OR. a_i(ji,jj,jl) < 0.01_wp ) THEN   ! Case ice thickness < rn_himin or tiny ice fraction !
               !                                                                             !----------------------------------------------------!
               !--- Remove ponds on thin ice or tiny ice fractions
               a_ip(ji,jj,jl) = 0._wp
               v_ip(ji,jj,jl) = 0._wp
               h_ip(ji,jj,jl) = 0._wp
               !                                                                             !--------------------------------!
            ELSE                                                                             ! Case ice thickness >= rn_himin !
               !                                                                             !--------------------------------!
               !--- Available and contributing meltwater for melt ponding
               zv_mlt  = - ( dh_i_sum_2d(ji,jj,jl) * rhoi + dh_s_mlt_2d(ji,jj,jl) * rhos ) * z1_rhow * a_i(ji,jj,jl) ! available volume of surface melt water per grid area
               ! MV -> could move this directly in ice_thd_dh and get an array (ji,jj,jl) for surface melt water volume per grid area
               IF( ln_pnd_rain )    zv_mlt = zv_mlt + ( tprecip(ji,jj) - sprecip(ji,jj) ) * z1_rhow * a_i(ji,jj,jl) * rDt_ice

               zfr_mlt = rn_apnd_min + ( rn_apnd_max - rn_apnd_min ) * at_i(ji,jj)                  ! fraction of surface meltwater going to ponds
               zv_pnd  = zfr_mlt * zv_mlt                                                           ! contributing meltwater volume for category jl
               
                  diag_dvpn_mlt(ji,jj) = diag_dvpn_mlt(ji,jj) + rhow * zv_mlt * r1_Dt_ice              ! diags
                  diag_dvpn_rnf(ji,jj) = diag_dvpn_rnf(ji,jj) - rhow * ( 1._wp - zfr_mlt ) * zv_mlt * r1_Dt_ice
               
               !--- Create possible new ponds
               ! if pond does not exist, create new pond over full ice area
               !IF ( a_ip_frac(ji,jj,jl) < epsi10 ) THEN
               IF( a_ip(ji,jj,jl) < epsi10 ) THEN
                  a_ip     (ji,jj,jl) = a_i(ji,jj,jl)
                  a_ip_frac(ji,jj,jl) = 1.0_wp    ! pond fraction of sea ice (apnd for CICE)
               ENDIF
               
               !--- Deepen existing ponds with no change in pond fraction, before redistribution and drainage
               v_ip(ji,jj,jl) = v_ip(ji,jj,jl) + zv_pnd                                            ! use pond water to increase thickness
               h_ip(ji,jj,jl) = v_ip(ji,jj,jl) / a_ip(ji,jj,jl)

            ENDIF
            
            !--- Total available pond water volume (pre-existing + newly produced)
            zvolp(ji,jj) = zvolp(ji,jj) + v_ip(ji,jj,jl)

         END_2D
      END DO

      zvolp_ini(:,:) = zvolp(:,:)

      !--------------------------------------------------------------
      ! Redistribute and drain water from ponds
      !--------------------------------------------------------------
      CALL ice_thd_pnd_area( zvolp, a_ip, v_ip, h_ip )

      !--------------------------------------------------------------
      ! Melt pond lid growth and melt
      !--------------------------------------------------------------

      IF( ln_pnd_lids ) THEN

         DO_2D( 0, 0, 0, 0 )

            !-----------------------------------------------------------------------------!
            ! Case ice thickness < rn_himin or tiny ice fraction or pond volume too small !
            !-----------------------------------------------------------------------------!
            IF ( at_i(ji,jj) < 0.01 .OR. vt_i(ji,jj) < (rn_himin*at_i(ji,jj)) .OR. zvolp_ini(ji,jj) < (zvp_min*at_i(ji,jj)) ) THEN

                  diag_dvpn_lid(ji,jj) = diag_dvpn_lid(ji,jj) + rhow * SUM( v_il(ji,jj,:) ) * r1_Dt_ice  ! diags
                  diag_dvpn_rnf(ji,jj) = diag_dvpn_rnf(ji,jj) - rhow * SUM( v_ip(ji,jj,:) + v_il(ji,jj,:) ) * r1_Dt_ice
               a_ip(ji,jj,:) = 0._wp
               v_ip(ji,jj,:) = 0._wp
               v_il(ji,jj,:) = 0._wp
               h_ip(ji,jj,:) = 0._wp

               !------------------------------------------------------------!
               ! Case ice thickness > rn_himin and pond volume large enough !
               !------------------------------------------------------------!
            ELSE

               !--------------------------
               ! Pond lid growth and melt
               !--------------------------
               DO jl = 1, jpl-1

                  IF ( v_il(ji,jj,jl) > epsi10 ) THEN

                     zdvice = MIN( -dh_i_sum_2d(ji,jj,jl)*a_ip(ji,jj,jl), v_il(ji,jj,jl) )

                     !----------------------------------------------------------------
                     ! Lid melting: floating upper ice layer melts in whole or part
                     !----------------------------------------------------------------
                     ! Use Tsfc for each category
                     IF ( t_su(ji,jj,jl) > zTp ) THEN


                        IF ( zdvice > epsi10 ) THEN

                           v_il (ji,jj,jl) = v_il (ji,jj,jl)   - zdvice
                           v_ip(ji,jj,jl)  = v_ip(ji,jj,jl)    + zdvice ! MV: not sure i understand dh_i_sum seems counted twice -
                           ! as it is already counted in surface melt

                           IF ( v_il(ji,jj,jl) < epsi10 .AND. v_ip(ji,jj,jl) > epsi10) THEN
                              ! ice lid melted and category is pond covered
                              v_ip(ji,jj,jl)  = v_ip(ji,jj,jl)  + v_il(ji,jj,jl)
                              diag_dvpn_lid(ji,jj) = diag_dvpn_lid(ji,jj) + rhow * v_il(ji,jj,jl) * r1_Dt_ice  ! diag
                              v_il(ji,jj,jl)   = 0._wp
                           ENDIF
                           h_ip(ji,jj,jl) = v_ip(ji,jj,jl) / a_ip(ji,jj,jl) !!! could get a division by zero here

                           diag_dvpn_lid(ji,jj) = diag_dvpn_lid(ji,jj) + rhow * zdvice * r1_Dt_ice ! diag

                        ENDIF

                        !----------------------------------------------------------------
                        ! Freeze pre-existing lid
                        !----------------------------------------------------------------

                     ELSE IF ( v_ip(ji,jj,jl) > epsi10 ) THEN ! Tsfcn(i,j,n) <= Tp

                        ! differential growth of base of surface floating ice layer
                        zdTice = MAX( - ( t_su(ji,jj,jl) - zTd ) , 0._wp ) ! > 0
                        zomega = rcnd_i * zdTice / zrhoi_L
                        ! DML: The official code has a_i here where it should
                        ! have a_ip (two instances).
                        ! The following corrects that.
                        zdHui  = SQRT( 2._wp * zomega * rDt_ice + ( v_il(ji,jj,jl) / a_ip(ji,jj,jl) )**2 ) &
                           - v_il(ji,jj,jl) / a_ip(ji,jj,jl)
                        zdvice = min(zdHui*a_ip(ji,jj,jl) , v_ip(ji,jj,jl) )

                        IF ( zdvice > epsi10 ) THEN
                           v_il (ji,jj,jl)  = v_il(ji,jj,jl)   + zdvice
                           v_ip(ji,jj,jl)   = v_ip(ji,jj,jl)   - zdvice
                           h_ip(ji,jj,jl)   = v_ip(ji,jj,jl) / a_ip(ji,jj,jl)

                           diag_dvpn_lid(ji,jj) = diag_dvpn_lid(ji,jj) - rhow * zdvice * r1_Dt_ice  ! diag
                        ENDIF

                     ENDIF ! Tsfcn(i,j,n)

                     !----------------------------------------------------------------
                     ! Freeze new lids
                     !----------------------------------------------------------------
                     !  upper ice layer begins to form
                     ! note: albedo does not change

                  ELSE ! v_il < epsi10

                     ! thickness of newly formed ice
                     ! the surface temperature of a meltpond is the same as that
                     ! of the ice underneath (0C), and the thermodynamic surface
                     ! flux is the same

                     ! we need net surface energy flux, excluding conduction
                     ! fsurf is summed over categories in CICE
                     ! we have the category-dependent flux, let us use it ?
                     zfsurf = qns_ice(ji,jj,jl) + qsr_ice(ji,jj,jl)
                     zdHui  = MAX ( -zfsurf * rDt_ice/zrhoi_L , 0._wp )
                     zdvice = MIN ( zdHui * a_ip(ji,jj,jl) , v_ip(ji,jj,jl) )
                     IF ( zdvice > epsi10 ) THEN
                        v_il (ji,jj,jl)  = v_il(ji,jj,jl)   + zdvice
                        v_ip(ji,jj,jl)   = v_ip(ji,jj,jl)   - zdvice

                        diag_dvpn_lid(ji,jj) = diag_dvpn_lid(ji,jj) - rhow * zdvice * r1_Dt_ice ! diag
                        h_ip(ji,jj,jl)   = v_ip(ji,jj,jl) / a_ip(ji,jj,jl) ! MV - in principle, this is useless as h_ip is computed in icevar
                     ENDIF

                  ENDIF  ! v_il

               END DO ! jl

               
            ENDIF

         END_2D

      ENDIF ! ln_pnd_lids

      !clem: conservation (quick fix)
      zdv_pnd2(:,:) = 0._wp
      DO jl = 1, jpl 
         DO_2D( 0, 0, 0, 0 )
            zdv_pnd2(ji,jj) = zdv_pnd2(ji,jj) + v_ip(ji,jj,jl) + v_il(ji,jj,jl)
         END_2D
      ENDDO

      DO_2D( 0, 0, 0, 0 )
         wfx_pnd(ji,jj) = wfx_pnd(ji,jj) - ( zdv_pnd2(ji,jj) - zdv_pnd1(ji,jj) ) * rhow * r1_Dt_ice
      END_2D

   END SUBROUTINE pnd_TOPO


   SUBROUTINE ice_thd_pnd_area( zvolp, a_ip, v_ip, h_ip )

       !!-------------------------------------------------------------------
       !!                ***  ROUTINE ice_thd_pnd_area ***
       !!
       !! ** Purpose : Given the total volume of available pond water,
       !!              redistribute and drain water
       !!
       !! ** Method
       !!
       !-----------|
       !           |
       !           |-----------|
       !___________|___________|______________________________________sea-level
       !           |           |
       !           |           |---^--------|
       !           |           |   |        |
       !           |           |   |        |-----------|              |-------
       !           |           |   | alfan  |           |              |
       !           |           |   |        |           |--------------|
       !           |           |   |        |           |              |
       !---------------------------v-------------------------------------------
       !           |           |   ^        |           |              |
       !           |           |   |        |           |--------------|
       !           |           |   | betan  |           |              |
       !           |           |   |        |-----------|              |-------
       !           |           |   |        |
       !           |           |---v------- |
       !           |           |
       !           |-----------|
       !           |
       !-----------|
       !
       !!
       !!------------------------------------------------------------------

       REAL (wp), DIMENSION(jpi,jpj), INTENT(INOUT) :: zvolp ! total available pond water
       REAL (wp), DIMENSION(:,:,:) , INTENT(INOUT) :: a_ip, v_ip, h_ip

       INTEGER ::  &
          ns,      &
          m_index, &
          permflag

       REAL (wp), DIMENSION(jpl) :: &
          hicen, &
          hsnon, &
          asnon, &
          alfan, &
          betan, &
          cum_max_vol, &
          reduced_aicen

       REAL (wp), DIMENSION(0:jpl) :: &
          cum_max_vol_tmp

       REAL (wp) :: &
          hpond, &
          drain, &
          floe_weight, &
          pressure_head, &
          hsl_rel, &
          deltah, &
          perm, &
          msno

       REAL (wp), parameter :: &
          viscosity = 1.79e-3_wp     ! kinematic water viscosity in kg/m/s

      REAL(wp), PARAMETER :: &   ! shared parameters for topographic melt ponds
         zvp_min = 1.e-4_wp          ! minimum pond volume (m)

      INTEGER  ::   ji, jj, jk, jl                    ! loop indices


       DO_2D( 0, 0, 0, 0 )


          !-----------------------------------------------------------------------------!
          ! Case ice thickness < rn_himin or tiny ice fraction or pond volume too small !
          !-----------------------------------------------------------------------------!
          IF ( at_i(ji,jj) < 0.01 .OR. vt_i(ji,jj) < (rn_himin*at_i(ji,jj)) .OR. zvolp(ji,jj) < (zvp_min*at_i(ji,jj)) ) THEN
             a_ip(ji,jj,:) = 0._wp
             h_ip(ji,jj,:) = 0._wp
             v_ip(ji,jj,:) = 0._wp
          ELSE
             
             !-------------------------------------------------------------------
             ! initialize
             !-------------------------------------------------------------------
             
             DO jl = 1, jpl

                !----------------------------------------
                ! compute the effective snow fraction
                !----------------------------------------

                IF (a_i(ji,jj,jl) < epsi10)  THEN
                   hicen(jl) =  0._wp
                   hsnon(jl) =  0._wp
                   reduced_aicen(jl) = 0._wp
                   asnon(jl) = 0._wp         !js: in CICE 5.1.2: make sense as the compiler may not initiate the variables
                ELSE
                   hicen(jl) = v_i(ji,jj,jl) / a_i(ji,jj,jl)
                   hsnon(jl) = v_s(ji,jj,jl) / a_i(ji,jj,jl)
                   reduced_aicen(jl) = 1._wp ! n=jpl

                   !js: initial code in NEMO_DEV
                   !IF (n < jpl) reduced_aicen(jl) = aicen(jl) &
                   !                     * (-0.024_wp*hicen(jl) + 0.832_wp)

                   !js: from CICE 5.1.2: this limit reduced_aicen to 0.2 when hicen is too large
                   IF( jl < jpl )   reduced_aicen(jl) = a_i(ji,jj,jl) * MAX(0.2_wp,(-0.024_wp*hicen(jl) + 0.832_wp))

                   asnon(jl) = reduced_aicen(jl)  ! effective snow fraction (empirical)
                   ! MV should check whether this makes sense to have the same effective snow fraction in here
                   ! OLI: it probably doesn't
                END IF

                ! This choice for alfa and beta ignores hydrostatic equilibium of categories.
                ! Hydrostatic equilibium of the entire ITD is accounted for below, assuming
                ! a surface topography implied by alfa=0.6 and beta=0.4, and rigidity across all
                ! categories.  alfa and beta partition the ITD - they are areas not thicknesses!
                ! Multiplying by hicen, alfan and betan (below) are thus volumes per unit area.
                ! Here, alfa = 60% of the ice area (and since hice is constant in a category,
                ! alfan = 60% of the ice volume) in each category lies above the reference line,
                ! and 40% below. Note: p6 is an arbitrary choice, but alfa+beta=1 is required.

                ! MV:
                ! Note that this choice is not in the original FF07 paper and has been adopted in CICE
                ! No reason why is explained in the doc, but I guess there is a reason. I'll try to investigate, maybe

                ! Where does that choice come from ? => OLI : Coz' Chuck Norris said so...

                alfan(jl) = 0.6 * hicen(jl)
                betan(jl) = 0.4 * hicen(jl)

                cum_max_vol(jl)     = 0._wp
                cum_max_vol_tmp(jl) = 0._wp

             END DO ! jpl

             cum_max_vol_tmp(0) = 0._wp
             drain = 0._wp

             !----------------------------------------------------------
             ! Drain overflow water, update pond fraction and volume
             !----------------------------------------------------------

             !--------------------------------------------------------------------------
             ! the maximum amount of water that can be contained up to each ice category
             !--------------------------------------------------------------------------
             ! If melt ponds are too deep to be sustainable given the ITD (OVERFLOW)
             ! Then the excess volume cum_max_vol(jl) drains out of the system
             ! It should be added to wfx_pnd_out

             DO jl = 1, jpl-1 ! last category can not hold any volume

                IF (alfan(jl+1) >= alfan(jl) .AND. alfan(jl+1) > 0._wp ) THEN

                   ! total volume in level including snow
                   cum_max_vol_tmp(jl) = cum_max_vol_tmp(jl-1) + (alfan(jl+1) - alfan(jl)) * SUM(reduced_aicen(1:jl))

                   ! subtract snow solid volumes from lower categories in current level
                   DO ns = 1, jl
                      cum_max_vol_tmp(jl) = cum_max_vol_tmp(jl) &
                         &                  - rhos/rhow * &     ! free air fraction that can be filled by water
                         &                    asnon(ns) * &     ! effective areal fraction of snow in that category
                         &                    MAX(MIN(hsnon(ns)+alfan(ns)-alfan(jl), alfan(jl+1)-alfan(jl)), 0._wp)
                   END DO

                ELSE ! assume higher categories unoccupied
                   cum_max_vol_tmp(jl) = cum_max_vol_tmp(jl-1)
                END IF
                !IF (cum_max_vol_tmp(jl) < z0) THEN
                !   CALL abort_ice('negative melt pond volume')
                !END IF
             END DO
             cum_max_vol_tmp(jpl) = cum_max_vol_tmp(jpl-1)  ! last category holds no volume
             cum_max_vol  (1:jpl) = cum_max_vol_tmp(1:jpl)

             !----------------------------------------------------------------
             ! is there more meltwater than can be held in the floe?
             !----------------------------------------------------------------
             IF (zvolp(ji,jj) >= cum_max_vol(jpl)) THEN
                drain = zvolp(ji,jj) - cum_max_vol(jpl) + epsi10
                zvolp(ji,jj) = zvolp(ji,jj) - drain ! update meltwater volume available

                diag_dvpn_rnf(ji,jj) = - diag_dvpn_rnf(ji,jj) - rhow * drain * r1_Dt_ice  
                ! diag - overflow counted in the runoff part (arbitrary choice)

                IF (zvolp(ji,jj) < epsi10) THEN
                   zvolp (ji,jj) = 0._wp
                END IF
             END IF

             ! height and area corresponding to the remaining volume
             ! routine leaves zvolp unchanged
             CALL ice_thd_pnd_depth(reduced_aicen, asnon, hsnon, alfan, zvolp(ji,jj), cum_max_vol, hpond, m_index)

             DO jl = 1, m_index
                !h_ip(jl) = hpond - alfan(jl) + alfan(1) ! here oui choulde update
                !                                         !  volume instead, no ?
                h_ip(ji,jj,jl) = MAX((hpond - alfan(jl) + alfan(1)), 0._wp)      !js: from CICE 5.1.2
                a_ip(ji,jj,jl) = reduced_aicen(jl)
                ! in practise, pond fraction depends on the empirical snow fraction
                ! so in turn on ice thickness
             END DO
             !zapond = sum(a_ip(1:m_index))     !js: from CICE 5.1.2; not in Icepack1.1.0-6-gac6195d

        !------------------------------------------------------------------------
        ! Drainage through brine network (permeability)
        !------------------------------------------------------------------------
        !!! drainage due to ice permeability - Darcy's law

        ! sea water level
        msno = 0._wp
        DO jl = 1 , jpl
          msno = msno + v_s(ji,jj,jl) * rhos
        END DO
        floe_weight = ( msno + rhoi*vt_i(ji,jj) + rho0*zvolp(ji,jj) ) / at_i(ji,jj)
        hsl_rel = floe_weight / rho0 &
                - ( ( sum(betan(:)*a_i(ji,jj,:)) / at_i(ji,jj) ) + alfan(1) )

        deltah = hpond - hsl_rel
        pressure_head = grav * rho0 * max(deltah, 0._wp)

        ! drain if ice is permeable
        permflag = 0

             IF (pressure_head > 0._wp) THEN
                DO jl = 1, jpl-1
                  IF ( hicen(jl) /= 0._wp ) THEN

                    perm = 0._wp ! MV ugly dummy patch
                    CALL ice_thd_pnd_perm(t_i(ji,jj,:,jl),  sz_i(ji,jj,:,jl), perm) ! bof
                    IF (perm > 0._wp) permflag = 1

                    drain = perm*a_ip(ji,jj,jl)*pressure_head*rDt_ice / &
                                          (viscosity*hicen(jl))
                    diag_dvpn_drn(ji,jj) = diag_dvpn_drn(ji,jj) - rhow * min(drain, zvolp(ji,jj)) * r1_Dt_ice  ! diag
                    zvolp(ji,jj) = max(zvolp(ji,jj) - drain, 0._wp)

                    IF (zvolp(ji,jj) < epsi10) THEN
                      diag_dvpn_drn(ji,jj) = diag_dvpn_drn(ji,jj) - rhow * zvolp(ji,jj) * r1_Dt_ice ! diag
                      zvolp(ji,jj) = 0._wp
                    END IF
                  END IF
                END DO

                ! adjust melt pond dimensions
                IF (permflag > 0) THEN
                   ! recompute pond depth
                   CALL ice_thd_pnd_depth(reduced_aicen, asnon, hsnon, alfan, zvolp(ji,jj), cum_max_vol, hpond, m_index)
                   DO jl = 1, m_index
                      h_ip(ji,jj,jl) = hpond - alfan(jl) + alfan(1)
                      a_ip(ji,jj,jl) = reduced_aicen(jl)
                   END DO
                   !zapond = sum(a_ip(1:m_index))       !js: from CICE 5.1.2; not in Icepack1.1.0-6-gac6195d
                END IF
             END IF ! pressure_head

             !-------------------------------
             ! remove water from the snow
             !-------------------------------
             !------------------------------------------------------------------------
             ! total melt pond volume in category does not include snow volume
             ! snow in melt ponds is not melted
             !------------------------------------------------------------------------
             
             ! MV here, it seems that we remove some meltwater from the ponds, but I can't really tell
             ! how much, so I did not diagnose it
             ! so if there is a problem here, nobody is going to see it...
             
             
             ! Calculate pond volume for lower categories
             DO jl = 1,m_index-1
                v_ip(ji,jj,jl) = a_ip(ji,jj,jl) * h_ip(ji,jj,jl) & ! what is not in the snow
                   &             - (rhos/rhow) * asnon(jl) * MIN(hsnon(jl), h_ip(ji,jj,jl))
             END DO

             ! Calculate pond volume for highest category = remaining pond volume

             ! The following is completely unclear to Martin at least
             ! Could we redefine properly and recode in a more readable way ?

             ! m_index = last category with melt pond

             IF (m_index == 1) v_ip(ji,jj,m_index) = zvolp(ji,jj) ! volume of mw in 1st category is the total volume of melt water

             IF (m_index > 1) THEN
                IF (zvolp(ji,jj) > SUM( v_ip(ji,jj,1:m_index-1))) THEN
                   v_ip(ji,jj,m_index) = zvolp(ji,jj) - SUM(v_ip(ji,jj,1:m_index-1))
                ELSE
                   v_ip(ji,jj,m_index) = 0._wp
                   h_ip(ji,jj,m_index) = 0._wp
                   a_ip(ji,jj,m_index) = 0._wp
                   ! If remaining pond volume is negative reduce pond volume of
                   ! lower category
                   IF ( zvolp(ji,jj) + epsi10 < SUM(v_ip(ji,jj,1:m_index-1))) &
                      v_ip(ji,jj,m_index-1) = v_ip(ji,jj,m_index-1) - SUM(v_ip(ji,jj,1:m_index-1)) + zvolp(ji,jj)
                END IF
             END IF

             DO jl = 1,m_index
                IF (a_ip(ji,jj,jl) > epsi10) THEN
                   h_ip(ji,jj,jl) = v_ip(ji,jj,jl) / a_ip(ji,jj,jl)
                ELSE
                   diag_dvpn_rnf(ji,jj) = diag_dvpn_rnf(ji,jj) - rhow * v_ip(ji,jj,jl) * r1_Dt_ice  ! diag
                   h_ip(ji,jj,jl) = 0._wp
                   v_ip(ji,jj,jl)  = 0._wp
                   a_ip(ji,jj,jl) = 0._wp
                END IF
             END DO
             DO jl = m_index+1, jpl
                h_ip(ji,jj,jl) = 0._wp
                a_ip(ji,jj,jl) = 0._wp
                v_ip(ji,jj,jl) = 0._wp
             END DO

          ENDIF

       END_2D

    END SUBROUTINE ice_thd_pnd_area


    SUBROUTINE ice_thd_pnd_depth(aicen, asnon, hsnon, alfan, zvolp, cum_max_vol, hpond, m_index)
       !!-------------------------------------------------------------------
       !!                ***  ROUTINE ice_thd_pnd_depth  ***
       !!
       !! ** Purpose :   Compute melt pond depth
       !!-------------------------------------------------------------------

       REAL (wp), DIMENSION(jpl), INTENT(IN) :: &
          aicen, &
          asnon, &
          hsnon, &
          alfan, &
          cum_max_vol

       REAL (wp), INTENT(IN) :: &
          zvolp

       REAL (wp), INTENT(OUT) :: &
          hpond

       INTEGER, INTENT(OUT) :: &
          m_index

       INTEGER :: n, ns

       REAL (wp), DIMENSION(0:jpl+1) :: &
          hitl, &
          aicetl

       REAL (wp) :: &
          rem_vol, &
          area, &
          vol, &
          tmp, &
          z0   = 0.0_wp

       !----------------------------------------------------------------
       ! hpond is zero if zvolp is zero - have we fully drained?
       !----------------------------------------------------------------
       
       IF (zvolp < epsi10) THEN
          hpond = z0
          m_index = 0
       ELSE

          !----------------------------------------------------------------
          ! Calculate the category where water fills up to
          !----------------------------------------------------------------

          !----------|
          !          |
          !          |
          !          |----------|                                     -- --
          !__________|__________|_________________________________________ ^
          !          |          |             rem_vol     ^                | Semi-filled
          !          |          |----------|-- -- -- - ---|-- ---- -- -- --v layer
          !          |          |          |              |
          !          |          |          |              |hpond
          !          |          |          |----------|   |     |-------
          !          |          |          |          |   |     |
          !          |          |          |          |---v-----|
          !          |          | m_index  |          |         |
          !-------------------------------------------------------------

          m_index = 0  ! 1:m_index categories have water in them
          DO n = 1, jpl
             IF (zvolp <= cum_max_vol(n)) THEN
                m_index = n
                IF (n == 1) THEN
                   rem_vol = zvolp
                ELSE
                   rem_vol = zvolp - cum_max_vol(n-1)
                END IF
                exit ! to break out of the loop
             END IF
          END DO
          m_index = min(jpl-1, m_index)

          !----------------------------------------------------------------
          ! semi-filled layer may have m_index different snow in it
          !----------------------------------------------------------------

          !-----------------------------------------------------------  ^
          !                                                             |  alfan(m_index+1)
          !                                                             |
          !hitl(3)-->                             |----------|          |
          !hitl(2)-->                |------------| * * * * *|          |
          !hitl(1)-->     |----------|* * * * * * |* * * * * |          |
          !hitl(0)-->-------------------------------------------------  |  ^
          !                various snow from lower categories          |  |alfa(m_index)

          ! hitl - heights of the snow layers from thinner and current categories
          ! aicetl - area of each snow depth in this layer

          hitl(:) = z0
          aicetl(:) = z0
          DO n = 1, m_index
             hitl(n)   = MAX(MIN(hsnon(n) + alfan(n) - alfan(m_index), alfan(m_index+1) - alfan(m_index)), z0)
             aicetl(n) = asnon(n)

             aicetl(0) = aicetl(0) + (aicen(n) - asnon(n))
          END DO

          hitl(m_index+1) = alfan(m_index+1) - alfan(m_index)
          aicetl(m_index+1) = z0

          !----------------------------------------------------------------
          ! reorder array according to hitl
          ! snow heights not necessarily in height order
          !----------------------------------------------------------------

          DO ns = 1, m_index+1
             DO n = 0, m_index - ns + 1
                IF (hitl(n) > hitl(n+1)) THEN ! swap order
                   tmp = hitl(n)
                   hitl(n) = hitl(n+1)
                   hitl(n+1) = tmp
                   tmp = aicetl(n)
                   aicetl(n) = aicetl(n+1)
                   aicetl(n+1) = tmp
                END IF
             END DO
          END DO

          !----------------------------------------------------------------
          ! divide semi-filled layer into set of sublayers each vertically homogenous
          !----------------------------------------------------------------

          !hitl(3)----------------------------------------------------------------
          !                                                       | * * * * * * * *
          !                                                       |* * * * * * * * *
          !hitl(2)----------------------------------------------------------------
          !                                    | * * * * * * * *  | * * * * * * * *
          !                                    |* * * * * * * * * |* * * * * * * * *
          !hitl(1)----------------------------------------------------------------
          !                 | * * * * * * * *  | * * * * * * * *  | * * * * * * * *
          !                 |* * * * * * * * * |* * * * * * * * * |* * * * * * * * *
          !hitl(0)----------------------------------------------------------------
          !    aicetl(0)         aicetl(1)           aicetl(2)          aicetl(3)

          ! move up over layers incrementing volume
          DO n = 1, m_index+1

             area = SUM(aicetl(:)) - &                 ! total area of sub-layer
                &   (rhos/rho0) * SUM(aicetl(n:jpl+1)) ! area of sub-layer occupied by snow

             vol = (hitl(n) - hitl(n-1)) * area      ! thickness of sub-layer times area

             IF (vol >= rem_vol) THEN  ! have reached the sub-layer with the depth within
                hpond = rem_vol / area + hitl(n-1) + alfan(m_index) - alfan(1)

                EXIT
             ELSE  ! still in sub-layer below the sub-layer with the depth
                rem_vol = rem_vol - vol
             END IF

          END DO

       END IF

    END SUBROUTINE ice_thd_pnd_depth

    SUBROUTINE ice_thd_pnd_perm(ticen, salin, perm)
       !!-------------------------------------------------------------------
       !!                ***  ROUTINE ice_thd_pnd_perm ***
       !!
       !! ** Purpose :   Determine the liquid fraction of brine in the ice
       !!                and its permeability
       !!-------------------------------------------------------------------

       REAL (wp), DIMENSION(nlay_i), INTENT(IN) :: &
          ticen,  &     ! internal ice temperature (K)
          salin         ! salinity (ppt)     !js: ppt according to cice

       REAL (wp), INTENT(OUT) :: &
          perm      ! permeability

       REAL (wp) ::   &
          Sbr       ! brine salinity

       REAL (wp), DIMENSION(nlay_i) ::   &
          Tin, &    ! ice temperature
          phi       ! liquid fraction

       INTEGER :: k

       !-----------------------------------------------------------------
       ! Compute ice temperatures from enthalpies using quadratic formula
       !-----------------------------------------------------------------

       DO k = 1,nlay_i
          Tin(k) = ticen(k) - rt0   !js: from K to degC
       END DO

       !-----------------------------------------------------------------
       ! brine salinity and liquid fraction
       !-----------------------------------------------------------------

       DO k = 1, nlay_i

          Sbr    = - Tin(k) / rTmlt ! Consistent expression with SI3 (linear liquidus)
          ! Best expression to date is that one (Vancoppenolle et al JGR 2019)
          ! Sbr  = - 18.7 * Tin(k) - 0.519 * Tin(k)**2 - 0.00535 * Tin(k) **3
          phi(k) = salin(k) / Sbr

       END DO
       !-----------------------------------------------------------------
       ! permeability
       !-----------------------------------------------------------------

       perm = 3.0e-08_wp * (minval(phi))**3 ! Golden et al. (2007)

   END SUBROUTINE ice_thd_pnd_perm

   SUBROUTINE ice_thd_pnd_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_thd_pnd_init   ***
      !!
      !! ** Purpose : Physical constants and parameters linked to melt ponds
      !!              over sea ice
      !!
      !! ** Method  :  Read the namthd_pnd  namelist and check the melt pond
      !!               parameter values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namthd_pnd
      !!-------------------------------------------------------------------
      INTEGER  ::   ios, ioptio   ! Local integer
      !!
      NAMELIST/namthd_pnd/  ln_pnd, ln_pnd_LEV , rn_apnd_min, rn_apnd_max, rn_pnd_flush, &
         &                          ln_pnd_CST , rn_apnd, rn_hpnd,         &
         &                          ln_pnd_TOPO,                           &
         &                          ln_pnd_lids, ln_pnd_rain, ln_pnd_alb,  &
         &                          rn_pnd_hl_min, rn_pnd_hl_max
      !!-------------------------------------------------------------------
      !
      READ  ( numnam_ice_ref, namthd_pnd, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namthd_pnd  in reference namelist' )
      READ  ( numnam_ice_cfg, namthd_pnd, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namthd_pnd in configuration namelist' )
      IF(lwm) WRITE ( numoni, namthd_pnd )
      !
      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_thd_pnd_init: ice parameters for melt ponds'
         WRITE(numout,*) '~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namicethd_pnd:'
         WRITE(numout,*) '      Melt ponds activated or not ln_pnd         = ', ln_pnd
         WRITE(numout,*) '         Topographic melt pond scheme ln_pnd_TOPO    = ', ln_pnd_TOPO
         WRITE(numout,*) '         Level ice melt pond scheme ln_pnd_LEV     = ', ln_pnd_LEV
         WRITE(numout,*) '            Minimum ice fraction that contributes to melt ponds   rn_apnd_min    = ', rn_apnd_min
         WRITE(numout,*) '            Maximum ice fraction that contributes to melt ponds   rn_apnd_max    = ', rn_apnd_max
         WRITE(numout,*) '            Pond flushing efficiency rn_pnd_flush   = ', rn_pnd_flush
         WRITE(numout,*) '         Constant ice melt pond scheme ln_pnd_CST     = ', ln_pnd_CST
         WRITE(numout,*) '            Prescribed pond fraction rn_apnd        = ', rn_apnd
         WRITE(numout,*) '            Prescribed pond depth rn_hpnd        = ', rn_hpnd
         WRITE(numout,*) '         Frozen lids on top of melt ponds ln_pnd_lids    = ', ln_pnd_lids
         WRITE(numout,*) '         Rain added to melt ponds ln_pnd_rain    = ', ln_pnd_rain
         WRITE(numout,*) '         Melt ponds affect albedo or not ln_pnd_alb     = ', ln_pnd_alb
         WRITE(numout,*) '         Lid thickness below which full pond area used in albedo  rn_pnd_hl_min  = ', rn_pnd_hl_min
         WRITE(numout,*) '         Lid thickness above which ponds disappear from albedo    rn_pnd_hl_max  = ', rn_pnd_hl_max
      ENDIF
      !
      !                             !== set the choice of ice pond scheme ==!
      ioptio = 0
      IF( .NOT.ln_pnd ) THEN   ;   ioptio = ioptio + 1   ;   nice_pnd = np_pndNO     ;   ENDIF
      IF( ln_pnd_CST  ) THEN   ;   ioptio = ioptio + 1   ;   nice_pnd = np_pndCST    ;   ENDIF
      IF( ln_pnd_LEV  ) THEN   ;   ioptio = ioptio + 1   ;   nice_pnd = np_pndLEV    ;   ENDIF
      IF( ln_pnd_TOPO ) THEN   ;   ioptio = ioptio + 1   ;   nice_pnd = np_pndTOPO   ;   ENDIF
      IF( ioptio /= 1 )   &
         & CALL ctl_stop( 'ice_thd_pnd_init: choose either none (ln_pnd=F) or only one pond scheme (ln_pnd_LEV, ln_pnd_CST or ln_pnd_TOPO)' )
      !
      SELECT CASE( nice_pnd )
      CASE( np_pndNO )
         IF( ln_pnd_alb  ) THEN ; ln_pnd_alb  = .FALSE. ; CALL ctl_warn( 'ln_pnd_alb=false when no ponds' )  ; ENDIF
         IF( ln_pnd_lids ) THEN ; ln_pnd_lids = .FALSE. ; CALL ctl_warn( 'ln_pnd_lids=false when no ponds' ) ; ENDIF
      CASE( np_pndCST )
         IF( ln_pnd_lids ) THEN ; ln_pnd_lids = .FALSE. ; CALL ctl_warn( 'ln_pnd_lids=false when constant ponds' ) ; ENDIF
      END SELECT
      !
   END SUBROUTINE ice_thd_pnd_init

#else
   !!----------------------------------------------------------------------
   !!   Default option          Empty module          NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icethd_pnd
