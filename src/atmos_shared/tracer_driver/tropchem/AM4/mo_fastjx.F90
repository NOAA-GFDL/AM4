module MO_FASTJX_MOD
!----------------------------------------------------------------------
!       FAST-JX Photolysis rate 
!       Revised by Jingyi Li
!       
!       Aug 12, 2015
!       Jingyi.Li@noaa.gov
!----------------------------------------------------------------------
!         Below is the original information from FASTJX code
! ---------- JX71_notes.f90

!  same as JX70b except that test for solar flux = 0 so that can run 200 nm data set
!  the 2200 nm data set cuts all J's at exactly 200 nm - used to combien with WACCM tables
!  to allow J's above 60 km, including Lyman-alpha
!
! ----subroutines and calls:
!       main standalone
!   >>>   call INIT_FJX (TITLJXX,NJX_,NJXX)
!         call RD_JS_JX(9,'FJX_j2j.dat', TITLJXX,NJXX)
!         call SOLAR_JX(GMTAU,IDAY,YGRD,XGRD, SZA,U0,SOLF)
!         call ACLIM_FJX (YLAT, MONTH, PPP,TTT,ZZZ,DDD,OOO, L1_)
!         call JP_ATM0(PPP,TTT,DDD,OOO,ZZZ, L_)
!   >>>   call PHOTO_JX(U0,SZA,REFLB,SOLF, LPRTJ,  PPP,ZZZ,TTT,DDD,RRR,OOO,
!                 LWP,IWP,REFFL,REFFI, AERSP,NDXAER,L1U,ANU,  VALJXX,NJXU)
!
! ----notes:   >>> = only two essential fast-JX calls are denoted above with >>>
!
!   >>> INIT_FJX (TITLJXX,NJX_,NJXX)
!          called once to read in data files, returns thru calling sequence:
!             TITLJXX(1:NJX_) the char*6 title for each fast-JX J-value
!             NJXX the actual number of J-values that fast-JX will compute.
!
!       RD_JS_JX(9,'FJX_j2j.dat', TITLJXX,NJXX)
!           called once after INIT_FJX to map the CTM J's onto the fast-JX J's
!           this is an example, see the data faile 'FJX_j2j.dat'
!
!   >>> PHOTO_JX(U0,SZA,REFLB,SOLF, LPRTJ,  PPP,ZZZ,TTT,DDD,RRR,OOO,
!                  CLDWP,AERSP,NDXCLD,NDXAER,L1_,AN_,    VALJXX,NJX_)
!          called every time for each Indep. Colm. Atmos. (ICA) to compute J's
!          all the information is passed through the calling sequence:
!     ----------------------------------------------------------------------------
!      L1_    CTM has layers 1:L_, for JX must add top-of-atmos layer, L1_ = L_+1
!      AN_    dimension of max number of aerosol types given to JX
!      U0     cosine of the solar zenith angle
!      SZA    solar zenith angle (degrees)
!      REFLB  Lambertian reflectance of lower boundary
!        (NB JX could be adapted to add lower layers for canopy/snow scatt/abs)
!      SOLF   Solar radiation factor, correcting for sun-earth distance
!        (NB JX could be modified for solar cycle, but needs link to internal data)
!      LPRTJ  logical to produce internal printout (stdout) of fast-JX J's, fluxes,
!         absorption, etc.  This information is now only internal to PHOTO_JX
!      PPP(1:L1_+1)   edge press (hPa)
!      ZZZ(1:L1_+1)   edge altitude (cm)
!      TTT(1:L1_)     mid-layer temp (K)
!      DDD(1:L1_)     layer dens-path (# molec /cm2)
!      RRR(1:L1_)     mid-layer relative humidity (0.0 to 1.0)
!      OOO(1:L1_)     layer O3 path (# O3 /cm2)
!      CLDWP(1:L1_)   layer cloud water path (kg/m2), liquid and ice
!      AERSP(1:L1_,1:AN_)  aerosol path (g/m2)
!      NDXCLD(1:L1_)  layer cloud index (type)
!          only a single cloud type is allowed for optical properties, pick
!          the dominant one in terms of optical depth,
!          see notes on cloud types allowed in fast-JX: 'FJX_scatt-cld.dat'
!      NDXAER(1:L1_,1:AN_) aerosol index (type)
!          sample aerosol types are in 'FJX_scat-aer.dat' and 'FJX_scat-UMa.dat'
!          the UMa data allows for relative humidity to be included
!          other aerosol types can be added.
!     ----------------------------------------------------------------------------
!      VALJXX(1:,NJX_,1:L) & NJX_ (first dimension of VALJXX) are returned
!          VALJXX is the array of fast-JX J's, the second dimension is not given
!             but is large enough to accommodate the CTM layers 1:L1_
!          the main code must use the information calcualted by RD_JS_JX to
!             re-map the VALJXX onto the CTM J's.  A useful example is given.
!     ----------------------------------------------------------------------------
!
!       SOLAR_JX calculates solar zenith angle & distance correction (if needed)
!
!       ACLIM_FJX fills in T & O3 from a climatology
!             may be needed for the layer above the CTM to account for O3 & O2
!
!       JP_ATM0 does a simple printout (stdout) of the atmosphere
!
!
! ---------- fjx70sub.f  fast-JX core routines ver 7.0+ (10/2012, mjp)
!
! ----subroutines and calls:  >>> only subroutines called from outside >>>
!      one include 'cmn_FJX.f' is common to several and has parameters, etc.
!      only other connection with CTM code is in call sequence and noted above.
!
!   >>> subroutine INIT_FJX (TITLEJXX,NJXU,NJXX)
!         call RD_XXX(JXUNIT,'FJX_spec.dat')
!         call RD_MIE(JXUNIT,'FJX_scat.dat')
!         call RD_UM (JXUNIT,'FJX_UMaer.dat')
!         call RD_PROF(JXUNIT,'atmos_std.dat')
!       subroutine RD_XXX(NJ1,NAMFIL)
!       subroutine RD_MIE(NJ1,NAMFIL)
!       subroutine RD_UM(NJ1,NAMFIL)
!       subroutine RD_PROF(NJ2,NAMFIL)
!       subroutine EXITC(T_EXIT)
!   >>> subroutine SOLAR_JX(GMTIME,NDAY,YGRDJ,XGRDI, SZA,COSSZA,SOLFX)
!   >>> subroutine ACLIM_FJX (YLATD, MONTH, PPP,TTT,ZZZ,DDD,OOO, L1U)
!   >>> subroutine PHOTO_JX(U0,SZA,REFLB,SOLF,LPRTJ, PPP,ZZZ,TTT,DDD,RRR,OOO,
!                          CLDWP,AERSP,NDXCLD,NDXAER,L1U,ANU,  VALJXX,NJXU)
!         call SPHERE2 (U0,RAD,ZZJ,ZZHT,AMF2, L1U,JXL1_)
!         call OPTICL (OPTX,SSAX,SLEGX,  ODCLD,NDCLD)
!         call OPTICA (OPTX,SSAX,SLEGX,  PATH,RH, NAER)
!         call OPTICM (OPTX,SSAX,SLEGX,  PATH,RH,-NAER)
!         call EXTRAL(OD600,L1U,L2U,N_,JTAUMX,ATAU,ATAU0, JXTRA)
!         call X_interp (TTTX,XQO2, TQQ(1,1),QO2(K,1), TQQ(2,1),QO2(K,2),.,,)
!         call X_interp (TTTX,XQO3, TQQ(1,2),QO3(K,1), TQQ(2,2),QO3(K,2),,,,)
!         call OPMIE (DTAUX,POMEGAX,U0,RFL,AMF2,JXTRA,
!                          FJACT,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0, LU)
!         call JRATET(PPJ,TTJ,FFF, VALJXX, LU,NJXU)
!         call JP_ATM(PPJ,TTJ,DDJ,OOJ,ZZJ,DTAU600,POMG600,JXTRA, LU)
!       subroutine OPMIE (DTAUX,POMEGAX,U0,RFL,AMF2,JXTRA,
!                          FJACT,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0, LU)
!         call MIESCT(FJ,FJT,FJB,POMEGA,FZ,ZTAU,ZFLUX,RFL,U0,ND)
!       subroutine MIESCT(FJ,FJT,FJB, POMEGA,FZ,ZTAU,ZFLUX,RFL,U0,ND)
!         call LEGND0 (EMU(I),PM0,M2_)
!         call LEGND0 (-U0,PM0,M2_)
!         call BLKSLV(FJ,POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,FJT,FJB, ND)
!       subroutine LEGND0 (X,PL,N)
!       subroutine BLKSLV(FJ,POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,FJTOP,FJBOT,ND)
!         call GEN_ID (POMEGA(1,1,K),FZ(1,K),ZTAU(1,K),ZFLUX(K),RFL(K),PM,PM0,
!              B(1,1,1,K),CC(1,1,1,K),AA(1,1,1,K),A(1,1,K),H(1,1,K),C(1,1,K),ND)
!       subroutine GEN_ID(POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,B,CC,AA,A,H,C,ND)
!       subroutine OPTICL (OPTD,SSALB,SLEG, ODCLD,NDCLD)
!       subroutine OPTICA (OPTD,SSALB,SLEG, PATH,RELH,L)
!       subroutine OPTICM (OPTD,SSALB,SLEG, PATH,RELH,L)
!       subroutine JRATET(PPJ,TTJ,FFF, VALJL,LU,NJXU)
!         call X_interp (TT,QO2TOT, TQQ(1,1),QO2(K,1),...)
!         call X_interp (TT,QO3TOT, TQQ(1,2),QO3(K,1),...)
!         call X_interp (TT,QO31DY, TQQ(1,3),Q1D(K,1),...)
!         call X_interp (PP,QQQT, TQQ(1,J),QQQ(K,1,J),...)
!         call X_interp (TT,QQQT, TQQ(1,J),QQQ(K,1,J),...)
!       subroutine X_interp (TINT,XINT, T1,X1, T2,X2, T3,X3, L123)
!       subroutine JP_ATM(PPJ,TTJ,DDJ,OOJ,ZZJ,DTAU6,POMEG6,JXTRA,LU)
!   >>> subroutine JP_ATM0(PPJ,TTJ,DDJ,OOJ,ZZJ, LU)
!       subroutine SPHERE2(U0,RAD,ZHL,ZZHT,AMF2, L1U,LJX1U)
!       subroutine EXTRAL(DTAUX,L1X,L2X,NX,JTAUMX,ATAU,ATAU0, JXTRA)
! note - the calls to EXITC are not listed here
!       e.g.,  call EXITC(' INIT_JX: invalid no. wavelengths')
!
!
!
!  >>>>>>>>>>>>>>>>current code revised to JX ver 7.0+ (10/12)<<<<<<<<<<<<
!
!  fastJX version 7.0+ (f90) - Prather notes (Jan 2013)
!
!---calculation of cloud optical depth in FJX-70b !!!
!---    assumes that clouds are 100% if in layer
!
!   IWP = ice water path (in layer, in cloud) in g/m**2
!   LWP = liquid water path (in layer, in cloud) in g/m**2
!   REFFI = effective radius of ice cloud (microns)
!   REFFL = effective radius of liquid cloud (microns)
!
!>>>>method for calculating cloud OD (in layer) used by FJX core or outside
!>>>>FJX core needs only the _WP and the REFF_
!>>>> note that FJX can use correct Q's and density updates if need be.
!   ODL = LWP * 0.75 * 2.10 / REFFL
!   ODI = IWP * 0.75 * 2.00 / (REFFI * 0.917)
!
!>>>R-effective determined by main code, not FJX
!   REFF determined by user - some recommendations below (from Neu & Prather)
!       REFFI is a simple function of ice water content IWC (g/m3, 0.0001 to 0.1)
!          IWC = IWP / delta-Z (of layer in m, approx OK)
!          REFFI = 50. * (1. + 8.333 * IWC)
!   prefer Heymsfield++ 2003 JAM, log-log fit ext(/m) vs. IWC, Fig B1a, p.1389
!              EXT (/m) = 1.7e-3 * (IWC/0.1)**0.77
!          REFFI = 164. * IWC**0.23     (33 microns at 0.001 --- 164 at 1.0)
!
!          REFFL is a simple function of pressure (PCLD):
!            FACTOR = (PCLD - 610.) / 200.
!            FACTOR = min(1.0, max(0.0, FACTOR))
!          REFFL = 9.60*FACTOR + 12.68*(1.-FACTOR)
!
!>>>indices for cloud scattering determined by FJX core, not main code.
!   NDX = cloud scattering index based on REFFL or TCLD(for ice cloud)
!          NDXI = 6  ! ice hexag (cold)
!             if (TCLD .ge. 233.15) then
!          NDXI = 7  ! ice irreg
!             endif
!          NDXC = 1
!             do I=2,4
!             if (REFFL .gt. 0.5*(RCC(I-1)+RCC(I))) then
!          NDXC = I
!             endif
!             enddo
!
!
! version 7.0
!       New modular structure:
!           Designed for .f90 and CAM5
!           All set up routines grouped into one super-call
!           Separate interface with CTM and ICAs calls PHOT_J
!           PHOT_JX only sees a single column atmosphere (or ICA) thru calling params
!
! version 6.8
!       New layout and formatting of FJX_spec.dat (required)
!              allows 1,2, or 3 T's (or P's) for all X-sects
!       VOCs mostly use Pressure interp.
!
! version 6.7  (should output same J-values as ver 6.6 EXCEPT for VOCs)   (3/12)
!       New JPL-2010 update complete through all VOCs (see notes in FJX_spec.dat)
!           Most important change is incorporation of Stern-Volmer pressure-dependent
!           and wavelength-dependent quantum yields into the Temperature interpolation.
!           Acetone now has only one pair of X sections for each pathway.
!       Redo of mapping of (external) chemical model J's onto fastJX J's
!           see examples of splitting a J with fixed branching q-yields
!           see how chem model labels are not tied by the cross section labels in fastJX
!       Changes in FX_spec.dat make it incompatible with earlier versions, although
!           all of the Xsection data have the same format.
!       Now the number of X sections read in equals the number of J's that fast-JX calculates.
!       As before, all J's except for O2, O3, O3(1D) have only pairs of data at different T.
!
! version 6.6x
!       N.B. SPECIAL FIX of tropospheric bin 11 used in FULL(W_=18) and TROP-ONLY(W_=12)
!            Attenuation of this bin is too weak since it mixes 218 nm and 288 nm
!            Thus it allows a low level of J-O2 (and related 218 nm Xsects) thru-out trop.
!            The attenuation in strat (1e-4) is OK, but must zero out in trop (>100 hPa)
!            The TROP-QUICK (W_=8) does not use bin#11 and is unaffected.
!       Redo of the 4x4 matrix inversions from subroutine to inline code (faster)
!       FJX_spec.dat:  includes new JPL2010 and solar flux, but J-VOC unchanged (a to-do).
!           The J-O2 and O3 steady-state are improved, but the NOy:N2O is too small.
!           The J-NO appears too large based on both photocomp2008 & NOy values.
!             Possibly lack of thermopsheric NO absorption, or correlation with S-R bands.
!           ACTION: J-NO X-sections multiplied by 0.6 in new FJX_spec.dat
!
! version 6.5
!      J-values are now averaged over CTM layer (and not just mid-layer point)
!        This is important when cloud optical depth is large (~1).
!
! version 6.4
!     allows for shortened, speeded up troposphere versions<<<<<<
!     STD:           W_=18
!        identical results to v-6.2 if cloud OD is consistent
!     TROP-ONLY:     W_=12
!        collapses the wavelength bins from 18 to 12 (5-6-7-8 & 11-18)
!        drops many 'stratospheric' cross-sections (denoted by 'x' in 2nd title)
!        allows use of single standard spectral data set:  FJX_spec.dat
!        results close to W_=18, largest difference is J-O2 (<1% in 13-18 km!!)
!        This is recommended as accurate for troposphere only calculations.
!     TROP-QUICK:    W_=8
!        reverts to original fast-J 7-bins (12-18) plus 1 scaled UV (5) for J-O2
!        errors in 12-18 km range for J-O2, high sun are 10%, worse above.
!     ***Photolysis of O2 in the upper tropical troposphere is an important
!        source of O3.  It needs to be included in tropospheric runs.
!        TROP-ONLY is recommended, W_=8 is a quick fix if speed essential.
!
!     Major rewrite of code to minimize calls and allow better vector-type ops.
!     loop over wavelengths internal to Mie soln.
!     Driven by profiling of CTM code, may still need optimization.
!     Wavelengths can be contracted to W_=12 (trop only) and strat-only
!        X-sections are dropped.  With parm W_=18, the std fast-JX is retrieved.
!     Many call eliminated and summations in BLKSLV and GEN_ID are explicit
!     GEN_ID replaces GEN and calculates all matrix coeff's (1:L_) at once
!     RD_XXX changed to collapse wavelengths & x-sections to Trop-only:
!           WX_ = 18 (parm_CTM.f) should match the JX_spec.dat wavelengths
!           W_ = 12 (Trop-only) or 18 (std) is set in (parm_MIE.f).
!       if W_=12 then drop strat wavels, and drop x-sects (e.g. N2O, ...)
!
! version 6.3
!     revise cloud/aerosol OD & wavelength properties for CTM link:
!         OPTICL is new sub for cloud optical properties, but it
!              now starts with cloud OD @ 600 nm and cloud NDX
!              fast-JX now uses cloud NDX to scale OD to other wavelengths
!         OPTICA & OPTICM are new subs to convert aerosol path (g/m2) to OD
!              A is std UCI scat data
!              M is U Michigan data tables for aerosols, includes Rel Hum effect
!     drop sub GAUSSP and put into Parameter statement (parm_MIE.f)
!
! version 6.2
!     corrects a long-standing problem at SZA > 89 degrees.
!     In prior versions the ray-tracing of the path (and air-mass functions)
!     back to the sun was done at the edges of the CTM layers (it was developed
!     for the grid-point J-value code at Harvard/GISS/UCI).  This left the
!     interpolation to the mid-layer (needed for J's) open.  The prior method
!     gave irregular fluctuations in the direct solar beam at mid-layer for
!     large SZA > 88.  This is now corrected with exact ray-tracing from
!     the mid-pt of each CTM layer.  For small SZA, there is no effective
!     difference, for large SZA, results could be erratic.
!   v-6.2 fix should be easy if you have migrated to v6.1, else some minor
!      caution may be needed:
!      replace sub SPHERE with SPHERE2, AMF2 report factors for mid and egdes.
!      replace sub OPMIE with new OPMIE, this uses the new AMF2 correctly.
!      replace sub PHOTOJ with new PHOTOJ, this just hands off AMF2 from
!            SPHERE2 to OPMIE.
!
! version 6.1 adds
!      6.1b simplifies calling sequences feeds solar factor, albedo, to PHOTOJ
!         and read LAT, LNG directly.  No substantive changes.
!      new read-in of scat data for clouds/aerosols to allow for UMich data
!      This has required substantial rewrite of some of the core subroutines:
!         OPMIE is now called for each wavelength and without aersol/cloud data
!              all subs below OPMIE are unchanged
!         OPTICD & OPTICM are new subs to convert path (g/m2) to OD and phase fn
!              D is std UCI scat data (re-ordered for clouds 1st)
!              M is U Michigan data tables for aerosols, includes Rel Hum effect
!         PHOTOJ now assembles the aerosol data (better for CTM implementation)
!      This version can reproduce earlier versions exactly, but the test input
!         is changed from OD and NDX to PATH (g/m2) and NDX.
!
! version 6.0 adds
!      new 200-nm scattering data so that stratospheric aerosols can be done!
! version 5.7
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<begin CTM-specific subroutines<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<fastJX codes<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<only access to external variable thru cmn_FJX.f and calls<<<<<<
!<<<<<<<<<<<<<<<<<<<<version 7.0+  (3/2013, mjp)<<<<<<<<<<<<<<<<<<<<<<<<

!
! !MODULE: FJX_INIT (j2l: adapted to fastjx_init in AM3)
!
! !DESCRIPTION: FJX70_INIT contains routines to input fast-JX data
!
!
! !INTERFACE:
! subroutines:
! 
!	INIT_FJX (TITLEJXX,NJXU,NJXX)
!         INIT_JX is called only once to read in and store all fast-JX data
!
!	SOLAR_JX(GMTAU,IDAY,YGRD,XGRD, SZA,U0,SOLF)
!         SOLAR_JX is called only once per grid-square to set U0, etc.
!
!	ACLIM_FJX (YLAT, MONTH, PPP,TTT,ZZZ,DDD,OOO, L1_)
!         ACLIM_FJX sets up climatologies for O3, T, Density and Z.
!
!	RD_XXX(JXUNIT,'FJX_spec.dat')
!         Read in fast-J X-sections (spectral data)
!
!	RD_CLD(JXUNIT,'FJX_scat-cld.dat')
!         Read in cloud scattering data
!
!	RD_MIE(JXUNIT,'FJX_scat-aer.dat')   
!         Read in aerosols scattering data
!
!	RD_UM (JXUNIT,'FJX_scat-UMa.dat')
!         Read in UMich aerosol scattering data
!
!	RD_PROF(JXUNIT,'atmos_std.dat')
!         Read in T and O3 climatology used to fill e.g. upper layers or if O3 not calc.
!
!	RD_JS_JX(9,'FJX_j2j.dat', TITLEJXX,NJXX)
!         Read in photolysis rates used in chemistry code and mapping onto FJX J's
!         CTM call:  read in J-values names and link to fast-JX names
!
! END MODULE FJX_INIT_MOD
!<<<<<<<<<<<<<<<<<<<<<begin core fast-J subroutines<<<<<<<<<<<<<<<<<<<<<
!
! MODULE FJX_SUB_MOD (j2l: adapted to fastjx_photo in AM3)
! !DESCRIPTION: JX version 7.1  (10/13)  works with FL=0 for <200 nm
!
!
! !INTERFACE:
! subroutines:
!	PHOTO_JX(U0,SZA,REFLB,SOLF, LPRTJ,PPP,ZZZ,TTT,DDD,RRR,OOO,            &
!       LWP,IWP,REFFL,REFFI, AERSP,NDXAER,L1U,ANU,  VALJXX,NJXU)
!         PHOTO_JX is the gateway to fast-JX calculations:
!          calc J's for a single column atmosphere (aka Indep Colm Atmos or ICA)
!          needs P, T, O3, clds, aersls; adds top-of-atmos layer from climatology
!          needs day-fo-year for sun distance, SZA (not lat or long)
!	
!	SPHERE2 (U0,RAD,ZZJ,ZZHT,AMF2, L1U,JXL1_)
!         calculate the optical properties (opt-depth, single-scat-alb, phase-fn(1:8))
!          at the 5 std wavelengths 200-300-400-600-999 nm for cloud+aerosols
!
!	OPTICL (OPTX,SSAX,SLEGX,  ODL,NDXL)
!         UCI CLOUD data sets, calculate scattering properties,          
!           scales OD at 600 nm to JX wavelength        
!
!	OPTICA (OPTX,SSAX,SLEGX,  PATH,RH, NAER)
!       OPTICM (OPTX,SSAX,SLEGX,  PATH,RH,-NAER)
!         UC Irvine & U Michigan aerosol data sets,calculate scattering properties
!
!	EXTRAL(OD600,L1U,L2U,N_,JTAUMX,ATAU,ATAU0, JXTRA)
!         add sub-layers (JXTRA) to thick cloud/aerosol layers 
!
!	OPMIE (DTAUX,POMEGAX,U0,RFL,AMF2,JXTRA, &
!              AVGF,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0, LU)
!         calculate mean intensity (actinic) at each CTM levels    
!          calculate fluxes and deposition (heating rates)        
!
!	JRATET(PPJ,TTJ,FFF, VALJXX, LU,NJXU)
!         mapping J-values from fast-JX onto CTM chemistry is done in main
!
!	JP_ATM(PPJ,TTJ,DDJ,OOJ,ZZJ,DTAU600,POMG600,JXTRA, LU)
!         print out atmosphere used in J-value calc.
!
!	SOLAR_JX(GMTIME,NDAY,YGRDJ,XGRDI, SZA,COSSZA,SOLFX)
!         calc SZA and Solar Flux factor for given lat/lon/UT
!
!	ACLIM_FJX (YLATD, MONTH, PPP,TTT,ZZZ,DDD,OOO, L1U)
!         Load fast-JX climatology for latitude & month given pressure grid
!
!<<<<<<<<<<<<<<<<<<<<<<<begin core scattering subroutines<<<<<<<<<<<<<<<
!                                                                       
!	MIESCT(FJ,FJT,FJB, POMEGA,FZ,ZTAU,ZFLUX,RFL,U0,ND)                  
!                                                                       
!	LEGND0 (X,PL,N)                                                     
!                                                                       
!	BLKSLV (FJ,POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,FJTOP,FJBOT,ND)          
!                                                                       
!	GEN_ID(POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,B,CC,AA,A,H,C,ND)            
!
! END MODULE FJX_SUB_MOD
!----------------------------------------------------------------------------------                                                                       

      use              fms_mod, only : file_exist,              &
                                       write_version_number,    &
                                       mpp_pe,                  &
                                       mpp_root_pE,             &
                                       close_file,              &
                                       stdlog,                  &
                                       mpp_clock_begin, mpp_clock_end, &
                                       mpp_clock_id, CLOCK_MODULE, &
                                       check_nml_error, error_mesg, &
                                       open_namelist_file, FATAL
      use           mpp_io_mod, only : mpp_open, mpp_close, MPP_RDONLY, &
                                       MPP_ASCII, MPP_SEQUENTIAL,   &
                                       MPP_MULTI, MPP_SINGLE 
      use           time_manager_mod, only : time_type, get_date
      use tropchem_types_mod, only : tropchem_opt
      implicit none
      
      private
      public :: fastjx_init, fastjx_end, fastjx_photo, JVN_

!-----------------------------------------------------------------------
!     version number and tagname.
!-----------------------------------------------------------------------
      character(len=128)            :: version     = '$Id: mo_fastjx.F90,v 19.0.4.2 2015/01/30 17:46:00 Jingyi.Li Exp $'
      character(len=128)            :: tagname     = '$Name: siena_201211 $'
!-------------------------------------------------------------------------
!     ************************
!     *  parameters overview  *
!     ************************
!     L_     =  altitude(levels) dim of CTM grid
!     L1_    =  assuming another layer above the top layer
!     L2_    =  no. levels in the Fast-JX grid that includes both layer edges and layer mid-points
!     JXL_   =  vertical(levels) dim for J-values computed within fast-JX
!     JXL1_  =  assuming another layer above the top layer
!     JXL2_  =  mx no. levels in the basic Fast-JX grid (mid-level)
!     SZAMAX =  Solar zenith angle cut-off, above which to skip calculation
!     WX_    =  dim = no. of wavelengths in input file
!     X_     =  dim = no. of X-section data sets (input data) 
!     A_     =  dim = no. of Aerosol/cloud Mie sets (input data)
!     C_     =  dim = no. of cld-data sets (input data)
!     W_     =  dim = no. of Wavelength bins:  =18 std, =12 trop only
!     JVN_   =  max no. of J-values
!     AN_    =  no of separate aerosols per layer (needs NDX for each)
!     NJX_   =  
!     N_     =  no. of levels in Mie scattering arrays
!     M_     =  no. of Gauss points used, must = 4 in fast_JX (no option)
!     M2_    =  2*M_ = 8, replaces MFIT
!     A_AM3  =  dim = no. of Aerosol/cloud Mie sets FOR AM3(input data)
!     RAD    =  Radius of Earth (cm)
!     ZZHT   =  Effective scale height above top of atmosphere (cm)
!     ATAU   =  heating rate (factor increase from one layer to the next)
!     ATAU0  =  minimum heating rate
!     JTAUMX =  maximum number of divisions (i.e., may not get to ATAUMN)
!     WBIN   =  Boundaries of wavelength bins
!     WL     =  Centres of wavelength bins - 'effective wavelength'
!     FL     =  Solar flux incident on top of atmosphere (cm-2.s-1)
!     QO2    =  O2 cross-sections
!     QO3    =  O3 cross-sections
!     Q1D    =  O3 => O(1D) quantum yield
!     QQQ    =  Supplied cross sections in each wavelength bin (cm2)
!     QRAYL  =  Rayleigh parameters (effective cross-section) (cm2)
!     TQQ    =  Temperature for supplied cross sections
!     LQQ    =  determine interpolation with T or P
!     TITLEJX=  Title for supplied cross sections, from 'FJX_spec.dat'
!     SQQ    =  Flag for supplied cross sections, from 'FJX_spec.dat'
!     NAA    =  Number of categories for scattering phase functions, is 30
!     NAA_AM3=  Number of categories for scattering phase functions for am3, is 880
!     MEE_AM3=  Aerosol mass extinction efficiency, MEE*colume mass =od
!     SAA_AM3=  Single scattering albedo
!     PAA_AM3=  Phase function: first 8 terms of expansion      
!     WAA_AM3=  Wavelengths for the NK supplied phase functions(nm)
!     TITLE0 =  Title of the first line in FJX_spec.dat
!     TITLE0_AM3 = Title of the first line in am3_scat.dat
!     TITLAA =  
!     RAA    =  Effective radius associated with aerosol type
!     DAA    =  rho
!     WAA    =  Wavelengths for the NK supplied phase functions(nm)
!     QAA    =  Aerosol scattering phase functions
!     SAA    =  Single scattering albedo
!     PAA    =  Phase function: first 8 terms of expansion      
!     QCC    =  Cloud scattering phase functions
!     WCC    =  5 Wavelengths for supplied phase functions
!     PCC    =  Phase function: first 8 terms of expansion
!     RCC    =  Effective radius associated with cloud type
!     SCC    =  Single scattering albedo
!     DCC    =  density (g/cm^3)
!     NCC    =  Number of categories for cloud scattering phase functions
!     WMM    =  U Michigan aerosol wavelengths
!     UMAER  =  U Michigan aerosol data sets
!     JFACTA =  Quantum yield (or multiplication factor) for photolysis
!     JLABEL =  Reference label identifying appropriate J-value to use
!     JIND   =  Set index arrays that map Jvalue(j) onto rates
!     NRATJ  =  number of J species
!     MASFAC = used to calc ZH
!     pfactor1 = pa => molec/cm2
!-----------------------------------------------------------------------
      integer, parameter                ::      L_      =48             !altitude(levels) dim of CTM grid
      integer, parameter                ::      L1_     =L_+1           !assuming another layer above the top layer
      integer, parameter                ::      L2_     =2*L_+2         !2*L1_ = 2*L_ + 2 = no. levels in the basic Fast-JX grid (mid-level_ 
      integer, parameter                ::      JXL_    =48             ! JXL_: vertical(levels) dim for J-values computed within fast-JX
      integer, parameter                ::      JXL1_   =JXL_+1
      integer, parameter                ::      JXL2_   =2*JXL_+2       ! JXL2_: 2*JXL_ + 2 = mx no. levels in the basic Fast-JX grid (mid-level)
      real*8,  parameter                ::      SZAMAX  =98.0d0         !Solar zenith angle cut-off, above which to skip calculation
      integer, parameter                ::      WX_     =18             !dim = no. of wavelengths in input file
!      integer, parameter                ::      X_      =72            !dim = no. of X-section data sets (input data)
      integer, parameter                ::      X_      =69             ! ++j2l modified to match AM3
      integer, parameter                ::      A_      =40             !dim = no. of Aerosol/cloud Mie sets (input data)
      integer, parameter                ::      C_      =16             ! C_   = dim = no. of cld-data sets (input data)
      integer, parameter                ::      W_      =18             !dim = no. of Wavelength bins:  =18 std, =12 trop only

!      integer, parameter                ::      JVN_    =101           ! max no. of J-values
      integer, parameter                ::      JVN_    =69             ! max no. of J-values
      integer, parameter                ::      AN_     =25             ! no of separate aerosols per layer (needs NDX for each)
!      integer, parameter                ::      NJX_    =100
      integer, parameter                ::      NJX_    =69             ! ++j2l modified to match AM3 
      integer, parameter                ::      N_      =601            !no. of levels in Mie scattering arrays
                                                                        !     = 2*NC+1 = 4*(L_+1) + 1`+ 2*sum(JADDLV)
      integer, parameter                ::      M_      =4              !no. of Gauss points used, must = 4 in fast_JX (no option)
      integer, parameter                ::      M2_     =2*M_           !2*M_ = 8, replaces MFIT
      integer, parameter                ::      A_AM3   =1000           !dim = no. of Aerosol/cloud Mie sets FOR AM3(input data)     
!-----------------------------------------------------------------------
!    4 Gauss pts = 8-stream
      real*8, dimension(M_), parameter  ::                            &
                        EMU = [.06943184420297d0, .33000947820757d0,  &
                               .66999052179243d0, .93056815579703d0]
      real*8, dimension(M_), parameter  ::                            &
                        WT  = [.17392742256873d0, .32607257743127d0,  &
                               .32607257743127d0, .17392742256873d0]
!-----INPHOT------------------------------------------------------------ 
      real*8, parameter                 ::      RAD     =6375.d5        !Radius of Earth (cm)  
      real*8, parameter                 ::      ZZHT    =5.d5           !Effective scale height above top of atmosphere (cm)
      real*8, parameter                 ::      ATAU    =1.120d0        ! ATAU: heating rate (factor increase from one layer to the next)
      real*8, parameter                 ::      ATAU0   =0.010d0        ! ATAU0: minimum heating rate
      integer, parameter                ::      JTAUMX  =(N_ - 4*L_)/2  ! JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)

!-----RD_XXX  <- FJX_spec.dat-------------------------------------------
      real*8                            ::      WBIN(WX_+1)             ! WBIN: Boundaries of wavelength bins
      real*8                            ::      WL(WX_)                 ! WL: Centres of wavelength bins - 'effective wavelength'
      real*8                            ::      FL(WX_)                 ! FL: Solar flux incident on top of atmosphere (cm-2.s-1)
      real*8                            ::      QO2(WX_,3)              ! QO2: O2 cross-sections
      real*8                            ::      QO3(WX_,3)              ! QO3: O3 cross-sections
      real*8                            ::      Q1D(WX_,3)              ! Q1D: O3 => O(1D) quantum yield
      real*8                            ::      QQQ(WX_,3,X_)           ! QQQ: Supplied cross sections in each wavelength bin (cm2)
      real*8                            ::      QRAYL(WX_+1)            ! QRAYL: Rayleigh parameters (effective cross-section) (cm2)
      real*8                            ::      TQQ(3,X_)               ! TQQ: Temperature for supplied cross sections
      integer                           ::      LQQ(X_)                 ! LQQ = 1, 2, or 3 to determine interpolation with T or P
      character*6                       ::      TITLEJX(X_)             ! TITLEJX: Title for supplied cross sections, from 'FJX_spec.dat'
      character*1                       ::      SQQ(X_)                 ! SQQ: Flag for supplied cross sections, from 'FJX_spec.dat'

!-----RD_MIE  <- FJX_scat.dat-------------------------------------------
      integer                           ::      NAA                     !Number of categories for scattering phase functions, is 30
      integer                           ::      NAA_AM3                 !Number of categories for scattering phase functions for am3, is 880
      real*8                            ::      MEE_AM3(5,A_AM3)        !Aerosol mass extinction efficiency, MEE*colume mass =od
      real*8                            ::      SAA_AM3(5,A_AM3)        !Single scattering albedo
      real*8                            ::      PAA_AM3(8,5,A_AM3)      !Phase function: first 8 terms of expansion      
      real*8                            ::      WAA_AM3(5,A_AM3)        !Wavelengths for the NK supplied phase functions(nm)
      character*78                      ::      TITLE0                  !Title of the first line in FJX_spec.dat
      character*78                      ::      TITLE0_AM3              !Title of the first line in am3_scat.dat
      character*20                      ::      TITLAA(A_)              ! A_ =40
      real*8                            ::      RAA(A_)                 !Effective radius associated with aerosol type
      real*8                            ::      DAA(A_)                 ! rho
      real*8                            ::      WAA(5,A_)               !Wavelengths for the NK supplied phase functions(nm)
      real*8                            ::      QAA(5,A_)               !Aerosol scattering phase functions
      real*8                            ::      SAA(5,A_)               !Single scattering albedo
      real*8                            ::      PAA(8,5,A_)             !Phase function: first 8 terms of expansion      
!---- Variables in file 'FJX_scat-cld.dat' (RD_CLD)

      real*8                            ::      QCC(5,C_)               ! QCC: Cloud scattering phase functions
      real*8                            ::      WCC(5,C_)               ! WCC: 5 Wavelengths for supplied phase functions
      real*8                            ::      PCC(8,5,C_)             ! PCC: Phase function: first 8 terms of expansion
      real*8                            ::      RCC(C_)                 ! RCC: Effective radius associated with cloud type
      real*8                            ::      SCC(5,C_)               ! SCC: Single scattering albedo
      real*8                            ::      DCC(C_)                 ! DCC: density (g/cm^3)
      integer                           ::      NCC                     ! NCC: Number of categories for cloud scattering phase functions

!---- Variables in file 'FJX_scat-UMa.dat' (RD_CLD)

      real*8                            ::      WMM(6)                  ! WMM: U Michigan aerosol wavelengths
      real*8                            ::      UMAER(3,6,21,33)        ! UMAER: U Michigan aerosol data sets

!-----RD_JS  <-  chem_Js.dat--------------------------------------------
      real*8                            ::      JFACTA(JVN_)            !Quantum yield (or multiplication factor) for photolysis
      character*7                       ::      JLABEL(JVN_)            !Reference label identifying appropriate J-value to use
      integer                           ::      JIND(JVN_)              !Set index arrays that map Jvalue(j) onto rates
      integer                           ::      NRATJ                   !number of J species
      logical                           ::      module_is_initialized = .false.      
      real*8,    parameter              ::      MASFAC = 100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)         !used to calc ZH
      real*8,    parameter              ::      pfactor1 = 1.e3/9.8/28.9644*6.02e23/1.e4                !pa => molec/cm2
      integer                           ::      NJX,NW1,NW2
      ! solar flux (solar constant)
      real, dimension(:,:), allocatable :: solflxtot
      real                              :: solflxtot_ann_1850, &
                                           solflxtot_ann_2300
      ! flux by band
      real, dimension(:,:,:), allocatable :: solflxband_all
      real, dimension(:),     allocatable :: solflxband_ann_1850, &
                                             solflxband_ann_2300  
      integer, parameter                ::      nbands   =  18      
      integer ::  first_yr, last_yr,   &
                  nvalues_per_year, numbands
      integer ::  years_of_data = 0

      CONTAINS
! <SUBROUTINE NAME="fastjx_init">
!   <OVERVIEW>
!     Initialize data for fastjx calculation
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine initializes the calculation of photolysis rates
!     for FAST-JX calculation
!   </DESCRIPTION>
!   <TEMPLATE>
!     call fastjx_init
!   </TEMPLATE>                                                                        
!-----------------------------------------------------------------------
      subroutine fastjx_init
!-----------------------------------------------------------------------
!  Routine to initialise photolysis rate data, called directly from the 
!  cinit routine in ASAD. Currently use it to read the JPL spectral data
!  and standard O3 and T profiles and to set the appropriate reaction in
!-----------------------------------------------------------------------
      implicit none 

      integer  NJXX
      character*6, dimension(NJX_) :: TITLEJXX

      integer  JXUNIT,J

      if (module_is_initialized) return

!----------------------------------------------------------------------
!     output version number and tagname to logfile.
!----------------------------------------------------------------------
      call write_version_number (version, tagname)

      if (W_.ne.8 .and. W_.ne.12 .and. W_.ne.18) then
         write(*,*) ' INIT_JX: invalid no. wavelengths', W_
         stop
      endif

!-->! Use channel 8 to read fastJX data files:                              
!-->      JXUNIT  = 8
 
! Read in fast-J X-sections (spectral data)
!-->      call RD_XXX(JXUNIT,'FJX_spec.dat')
      call RD_XXX('INPUT/FJX71_spec.dat')

! Read in cloud scattering data
!-->      call RD_CLD(JXUNIT,'FJX_scat-cld.dat')
      call RD_CLD('INPUT/FJX71_scat-cld.dat')

! Read in aerosols scattering data
!-->      call RD_MIE(JXUNIT,'FJX_scat-aer.dat')
!      call RD_MIE('INPUT/FJX_scat-aer.dat')

! Read in aerosol/cloud scattering data <<<<<<<<<<<<<<<<<< new fast-JX  
      call RD_MIE_AM3('INPUT/am3_scat.dat')
! Read in UMich aerosol scattering data
!-->      call RD_UM (JXUNIT,'FJX_scat-UMa.dat')
!      call RD_UM ('INPUT/FJX_scat-UMa.dat')

!-->! Read in T and O3 climatology used to fill e.g. upper layers or if O3 not calc.
!-->      call RD_PROF(JXUNIT,'atmos_std.dat')

      NJXX = NJX
      do J = 1,NJX
        TITLEJXX(J) = TITLEJX(J)
      enddo

! Read in photolysis rates used in chemistry code and mapping onto FJX J's
!---CTM call:  read in J-values names and link to fast-JX names
!-->      call RD_JS_JX(9,'FJX_j2j.dat', TITLEJXX,NJXX)
      call RD_JS_JX('INPUT/FJX71_j2j.dat', TITLEJXX,NJXX)

      call RD_SOLAR('INPUT/PHOTON_18bin.txt')
                                                                        
!-----------------------------------------------------------------------
!   mark module as initialized.
!-----------------------------------------------------------------------
         module_is_initialized = .true.                                                                     
    end subroutine fastjx_init                                         
 
                  
    subroutine fastjx_photo                                             &
          (U0, SOLF, phalf1, zhalf1, pfull1, tfull, XO3, pwt,           &
           lwc, cloudsf, REFLB, aerop, aeron,o3_column_top, ZPJ, Time,  & 
           time_varying_solarflux)
!-----------------------------------------------------------------------
!  fastJX version 7.2 (f90) - Prather notes (Jan 2013)

!---calculation of cloud optical depth in FJX-72 !!!
!---    assumes that clouds are 100% if in layer

!   IWP = ice water path (in layer, in cloud) in g/m**2
!   LWP = liquid water path (in layer, in cloud) in g/m**2
!   REFFI = effective radius of ice cloud (microns)
!   REFFL = effective radius of liquid cloud (microns)

!>>>>method for calculating cloud OD (in layer) used by FJX core or outside
!>>>>FJX core needs only the _WP and the REFF_
!>>>> note that FJX can use correct Q's and density updates if need be.
!   ODL = LWP * 0.75 * 2.10 / REFFL
!   ODI = IWP * 0.75 * 2.00 / (REFFI * 0.917)

!>>>R-effective determined by main code, not FJX
!   REFF determined by user - some recommendations below (from Neu & Prather)
!       REFFI is a simplle function of ice water content IWC (g/m3, 0.0001 to 0.1)
!          IWC = IWP / delta-Z (of layer in m, approx OK)
!          REFFI = 50. * (1. + 8.333 * IWC)
!   prefer Heymsfield++ 2003 JAM, log-log fit ext(/m) vs. IWC, Fig B1a, p.1389
!              EXT (/m) = 1.7e-3 * (IWC/0.1)**0.77
!          REFFI = 164. * IWC**0.23     (33 microns at 0.001 --- 164 at 1.0)

!          REFFL is a simple function of pressure (PCLD):
!            FACTOR = (PCLD - 610.) / 200.
!            FACTOR = min(1.0, max(0.0, FACTOR))
!          REFFL = 9.60*FACTOR + 12.68*(1.-FACTOR)

!>>>indices for cloud scattering determined by FJX core, not main code.
!   NDX = cloud scattering index based on REFFL or TCLD(for ice cloud)
!          NDXI = 6  ! ice hexag (cold)
!             if (TCLD .ge. 233.15) then
!          NDXI = 7  ! ice irreg
!             endif
!          NDXC = 1
!             do I=2,4
!             if (REFFL .gt. 0.5*(RCC(I-1)+RCC(I))) then
!          NDXC = I
!             endif
!             enddo
!>>>>>revise PHOTO_JX to drop NDXCLD, and to include
!   LWP, REFFL, IWP, REFFI only, use T, etc to get NDX.
!---------------------------------------------------------------------------

      implicit none
!--------------------key params sent to fast JX-------------------------
!---SZA = solar zenith angle, U0 = cos(SZA)
!---SOLF = solar flux factor for sun-earth distance
!---REFLB = Lambertian reflectivity at the Lower Boundary
!---ZPJ(layer, JV#) = 2D array of J-values, indexed to CTM chemistry code
!---turn internal PHOTO_JX print (unit=6) on/off
!
!---Independent Column Atmosphere data passed to PHOTO_JX:
!--- P = edge press (hPa), Z = edge alt (cm), T = layer temp (K)
!--- D = layer dens-path (# molec /cm2), O = layer O3 path (# O3 /cm2)
!---
!--- R = layer rel.hum.(fraction)
!--- CLDWP = cloud water path (g/m2), AERSP = aerosol path (g/m2)
!--- NDXCLD = cld index type (liq & size, ice & scatt-phase)
!--- NDXAER = aerosol index type
!--- NB. clouds have only a single type within an ICA layer, pick liquid or ice
!---     aerosols are dimensioned with up to AN_ different types in an ICA layer

!   PPP - pressure at CTM layer edge, PPP(L1_+1)=0 (top-of-atmos)
!---calling sequence variables
      real*8,  intent(in)                   :: U0,REFLB,SOLF
      real*8,  intent(in),dimension(L1_)    :: phalf1,&               !Pressure at boundaries (Pa) 49 layer
                                               zhalf1                 !height at layer boundaries      
      real*8,  intent(in),dimension(L_)     :: pfull1,&               !Pressure at mid-layer (Pa)  48
                                               tfull,&                !Temperature at mid-Layer (K)           
                                               XO3,&                  !O3 mixing ratio (VMR)              
                                               pwt                    !Column air density (Kg/m2)          

      real*8,  intent(in)                   ::  o3_column_top         !O3 column density, a small number
      real*8,  intent(in)                   ::  aerop(:,:)            !Aerosol mass path (g/m2) (L_,AERONUM)
      integer, intent(in)                   ::  aeron(:,:)            !(AERONUM) Aerosol type index, see FJX_scat.dat 
                                                                      !for aerosol type info 
      real*8,  intent(in)                   ::  cloudsf(:,:)          !Clouds fraction
      real*8,  intent(in)                   ::  lwc(:,:)              !(L_,ncloud=2)grid averaged liquid water content (Kg/Kg) 
      logical, intent(in)                   ::  time_varying_solarflux    ! solar cycle on fastjx?
!---reports out the JX J-values, upper level program converts to CTM chemistry J's
!      real*8, intent(out), dimension(L1_-1,NJXU)::  VALJXX
      real*8, intent(out), dimension(L1_-1,NJX_)::  ZPJ

!-----------------------------------------------------------------------
!--------key LOCAL atmospheric data needed to solve plane-parallel J----
!-----these are dimensioned JXL_, and must have JXL_ .ge. L_
      real*8, dimension(L1_)                ::  LWP,IWP,REFFL,REFFI,phalf
      real*8, dimension(L_)                 ::  pfull
      real*8, dimension(L1_+1)              ::  TTJ,DDJ,OOJ,ZZJ
      real*8, dimension(L1_+1)              ::  PPJ,RELH
      integer,dimension(L2_+1)              ::  JXTRA
      real*8, dimension(W_)                 ::  FJTOP,FJBOT,FSBOT,FLXD0,RFL
      real*8, dimension(L_, W_)             ::  AVGF, FJFLX
      real*8, dimension(L1_,W_)             ::  DTAUX, FLXD
      real*8, dimension(8,L1_,W_)           ::  POMEGAX
      real*8, dimension(L1_)                ::  DTAU600
      real*8, dimension(8,L1_)              ::  POMG600
      real*8, dimension(W_,L1_)             ::  FFX
      real*8, dimension(W_,8)               ::  FFXNET
      logical  :: LPRTJ                   ! set to false

!---flux/heating arrays (along with FJFLX,FLXD,FLXD0)
      real*8             :: FLXJ(L1_),FFX0,FXBOT,FABOT
      real*8             :: ODABS,ODRAY,ODI,ODL
      real*8             :: RFLECT,FREFS,FREFL,FREFI
      real*8             :: AMF2(2*L1_+1,2*L1_+1)
!------------key SCATTERING arrays for clouds+aerosols------------------
      real*8             :: OPTX(5),SSAX(5),SLEGX(8,5)
      real*8             :: OD(5,L1_),SSA(5,L1_),SLEG(8,5,L1_)
      real*8             :: OD600(L1_)
      real*8             :: PATH,RH,XTINCT
!------------key arrays AFTER solving for J's---------------------------
      real*8             :: FFF(W_,JXL_),VALJ(X_)
      real*8             :: FLXUP(W_),FLXDN(W_),DIRUP(W_),DIRDN(W_)

      integer            :: LU,L2U, I,J,K,L,M,KMIE,KW,NAER,NDXI,NDXL, RATIO(W_)
      real*8             :: XQO3,XQO2,DTAUC,WAVE, TTTX
!-----------------------------------------------------------------------
      real*8             :: FACTOR
      real*8             :: SZA, IWC
      real*8, dimension(L1_+1)  :: ZZZ                 ! ZZZ at boundaries
      real*8, dimension(L1_-1,NJX_)::  VALJXX
!-----------------------------------------------------------------------
      type(time_type),intent(in) :: Time
      integer :: yr, mo, day, hr, minute, sec
      real*8                 :: solar_constant
      real*8, dimension(18)  :: solflxband_now
      if(time_varying_solarflux) then
         call get_date(Time,yr,mo,day,hr,minute,sec)  !model GMT
         call get_solar_data(yr, mo, solar_constant, solflxband_now)
         FL = solflxband_now
      endif
      if (L1_ .gt. JXL1_) then
        write(*,*) ' PHOTO_JX: not enough levels in JX', L1_, JXL1_
        stop
      endif

      LU = L1_ - 1
      L2U = LU + LU + 2

      FFF(:,:) = 0.d0
      FREFI = 0.d0
      FREFL = 0.d0
      FREFS = 0.d0

!---check for dark conditions SZA > 98.0 deg => tan ht = 63 km
!                        or         99.0                 80 km
      SZA = acos(U0)*180./ 3.141592653589793d0
!      write(*,*) 'SZA ',SZA
      if (SZA .gt. 98.d0) then
          ZPJ(:,:) = 0.0
          goto 99
      endif

!-----load the amtospheric column data ( convert Pa to hPa)     
      phalf=phalf1/100.               
      pfull=pfull1/100.
!
!****NOTE!!! In fastjx, L=1 surface; In AM3 L=1 TOA
!
!-----calculate effective altitude of each CTM level edge (cm)
!        ZZZ(1) = 0.d0
!      do L = 1,L1_-1
!        SCALEH      = 1.3806d-19*MASFAC*TTT(L)
!        ZZZ(L+1) = ZZZ(L) -( LOG(PPP(L+1)/PPP(L)) * SCALEH )
!      enddo
!        ZZZ(L1_+1) = ZZZ(L1_) + ZZHT

       do L=1,L1_
          ZZZ(L)= zhalf1(L1_+1-L)*100.            ! height at boundaries (cm)
       end do
       ZZZ(L1_+1) = ZZZ(L1_) + ZZHT
!---load the amtospheric column data
!      do L = 1,L1_
!        PPJ(L) = PPP(L)
!        TTJ(L) = TTT(L)
!        DDJ(L) = DDD(L)
!        OOJ(L) = OOO(L)
!      enddo
      do L = 1,L1_
         PPJ(L) = phalf(L1_+1-L)   ! pressure at boundaries (hPa)
         if(L <L1_)then
            TTJ(L) = tfull(L1_ -L) !TJ(ILNG,JLAT,L) 
            DDJ(L) = pwt(L1_ -L)*6.023e22/28.92   !DM(ILNG,JLAT,L)  air den moleculars/cm2 
            OOJ(L) = XO3(L1_ -L)*DDJ(L)           !O3 column density (molecs/cm2) per layer
         else
            TTJ(L) = TTJ(L1_-1) + (TTJ(L1_-1)-TTJ(L1_-2))*(ZZZ(L1_+1)-ZZZ(L1_-1))/(ZZZ(L1_)-ZZZ(L1_-2))                         !K
            DDJ(L) = phalf(1)/pfactor1*100.     !hPa=>molec/cm2
            OOJ(L) = o3_column_top              !molec/cm2 
         end if
      enddo
!      write(*,*)'Assigned TTJ, DDJ, OOJ', TTJ(L1_), TTJ(L1_-1)
      PPJ(L1_+1) = 0.d0

!---calculate spherical weighting functions (AMF: Air Mass Factor)
      do L = 1,L1_+1
        ZZJ(L) = ZZZ(L)
      enddo
!      write(*,*)'Assigned ZZJ'
!>>>R-effective of clouds determined by main code, not FJX
!   REFF determined by user - some recommendations below
!       REFFI is a simple function of ice water content IWC (g/m3, 0.0001 to 0.1)
!          IWC = IWP / delta-Z (of layer in m, approx OK)
!   Heymsfield++ 2003 JAM, log-log fit ext(/m) vs. IWC, Fig B1a, p.1389
!              EXT (/m) = 1.7e-3 * (IWC/0.1)**0.77
!          REFFI = 164. * IWC**0.23     (33 microns at 0.001 --- 164 at 1.0)
!          REFFL is a simple function of pressure (PCLD):
!            FACTOR = (PCLD - 610.) / 200.
!            FACTOR = min(1.0, max(0.0, FACTOR))
!          REFFL = 9.60*FACTOR + 12.68*(1.-FACTOR)
!----set layer L1_ to zero
      LWP(:)  =  0.d0
      IWP(:)  =  0.d0
      REFFL(:) = 0.d0
      REFFI(:) = 0.d0

      do L = 1,L_
         if (TTJ(L) .gt. 253.d0) then
!           LWP(L) = CLDP(L)                    ! liquid cloud
            LWP(L) = lwc(L1_-L,1)*pwt(L1_-L)*1.e3
         else
!           IWP(L) = CLDP(L)                    ! ice cloud
            IWP(L) = lwc(L1_-L,2)*pwt(L1_-L)*1.e3 
         endif
         if (IWP(L) .gt. 1.d-5) then
            IWC = IWP(L) *100.d0 / (ZZZ(L+1)-ZZZ(L))
            IWC = max(0.001d0, IWC)
            REFFI(L) = 164.d0 * IWC**0.23d0
!       write(6,'(a,i3,3f10.4)') 'ICE:',L,IWP(L),IWC,REFFI(L)
         endif
         if (LWP(L) .gt. 1.d-5) then
!         PCLD = 0.5d0*(PPP(L)+PPP(L+1))                      ! PCLD = pressure of grid center, replaced by pfull
!         FACTOR = min(1.d0, max(0.d0, (PCLD-610.d0)/200.d0))
            FACTOR = min(1.d0, max(0.d0, (pfull(L1_-L)-610.d0)/200.d0))
            REFFL(L) = 9.60d0*FACTOR + 12.68d0*(1.-FACTOR)
         endif
      enddo


!      write(*,*)'Calculated REFFL, REFFI, IWC'
      RFLECT = REFLB

!-----------------------------------------------------------------------
      call SPHERE2 (U0,RAD,ZZJ,ZZHT,AMF2, L1_,JXL1_)
!-----------------------------------------------------------------------
!      write(*,*)'Calculated AMF2'

!---calculate the optical properties (opt-depth, single-scat-alb, phase-fn(1:8))
!---  at the 5 std wavelengths 200-300-400-600-999 nm for cloud+aerosols
      do L = 1,L1_
       do K=1,5
         OD(K,L)  = 0.d0
         SSA(K,L) = 0.d0
        do I=1,8
         SLEG(I,K,L) = 0.d0
        enddo
       enddo
      enddo

      do L = 1,L1_

!---Liquid Water Cloud
         if (LWP(L) .gt. 1.d-5 .and. REFFL(L) .gt. 0.1d0) then
!---extinction K(m2/g) = 3/4 * Q / [Reff(micron) * density(g/cm3)]
             if(L<L1_) then
               ODL = LWP(L) * 0.75d0 * 2.1d0 / REFFL(L) * cloudsf(L1_-L,1)**0.5   ! here 2.1 => QAA(4,NDCLD1) in v6.4, RHO=1
                                                                                  ! assuming approximate random overlap
             else
               ODL = 0.0d0
             endif
!             NDXL = 1
!             do I=2,4
!                if (REFFL(L) .gt. 0.5*(RCC(I-1)+RCC(I))) then
!                   NDXL = I
!                endif
!             enddo
             NDXL = 3    !   Re=12 um 
!             write(*,*)'L ODL ', L, ODL
             call OPTICL (OPTX,SSAX,SLEGX,  ODL,NDXL)
             do K=1,5
                OD(K,L)  = OD(K,L)  + OPTX(K)
                SSA(K,L) = SSA(K,L) + SSAX(K)*OPTX(K)
                do I=1,8
                   SLEG(I,K,L)=SLEG(I,K,L) + SLEGX(I,K)*SSAX(K)*OPTX(K)
!                   write(*,*)'SLEG Water: ',I,K,L, SLEG(I,K,L)
                enddo
             enddo
!>>>diagnostic print of cloud data:
!        write(*,'(a,i3,2f8.2,f8.4,f8.2,f8.4,i4)') &
!         'Liq Cld',L,PPJ(L),PPJ(L+1),LWP(L),REFFL(L),ODL,NDXL
         endif

!         write(*,*)'Called OPTICL for Liquid Water Cloud'
!---Ice Water Cloud
         if (IWP(L) .gt. 1.d-5 .and. REFFI(L) .gt. 0.1d0) then
            if(L<L1_) then
              ODI = IWP(L) * 0.75d0 * 2.0d0 / (REFFI(L) * 0.917d0) * cloudsf(L1_-L,2)**0.5  ! here 2.0 => QAA(4,NDCLD2) in v6.4, RHO = 0.917
                                                                                            ! assuming approximate random overlap
            else
              ODI = 0.0d0
            endif
            if (TTJ(L) .ge. 233.15d0) then
               NDXI = 7  ! ice irreg
            else
               NDXI = 6  ! ice hexag (cold)
            endif
!            write(*,*)'L ODI ', L, ODI, NDXI
            call OPTICL (OPTX,SSAX,SLEGX,  ODI,NDXI)
            do K=1,5
               OD(K,L)  = OD(K,L)  + OPTX(K)
               SSA(K,L) = SSA(K,L) + SSAX(K)*OPTX(K)
               do I=1,8
                  SLEG(I,K,L)=SLEG(I,K,L) + SLEGX(I,K)*SSAX(K)*OPTX(K)
!                  write(*,*)'SLEG Ice: ',I,K,L, SLEG(I,K,L)
               enddo
            enddo
!>>>diagnostic print of cloud data:
!        write(*,'(a,i3,2f8.2,f8.4,f8.2,f8.4,i4)') &
!         'Ice Cld',L,PPJ(L),PPJ(L+1),IWP(L),REFFI(L),ODI,NDXI
         endif
!         write(*,*)'Called OPTICL for Ice Water Cloud'

!---aerosols in layer: check aerosol index
!---this uses data from climatology OR from current CTM (STT of aerosols)

!---FIND useful way to sum over different aerosol types!
         do M = 1,size(aerop,2)
            if(L<L1_)then
               NAER = aeron(L1_-L,M)
               PATH = aerop(L1_-L,M) !g/m2
            else
               NAER = 0
               PATH = 0.D0
            endif
!---subroutines OPTICA & OPTICM return the same information:
!---  optical depth (OPTX), single-scat albedo (SSAX) and phase fn (SLEGX(8))
!---subs have slightly different inputs:
!---  PATH is the g/m2 in the layer, NAER in the cloud/aerosol index
!---  UMich aerosols use relative humidity (RH)
!----------------------------------------------------------------------------
! Replaced by OPTICAM3
!
            if (PATH .gt. 0.d0 .and. NAER .gt. 0) then
                call OPTICAM3 (OPTX,SSAX,SLEGX,  PATH, NAER)

               do K=1,5
                  OD(K,L)  = OD(K,L)  + OPTX(K)
                  SSA(K,L) = SSA(K,L) + SSAX(K)*OPTX(K)
                  do I=1,8
                     SLEG(I,K,L)=SLEG(I,K,L) + SLEGX(I,K)*SSAX(K)*OPTX(K)
                  enddo
               enddo
!               write(*,'(a,i3,f8.4,i4)') 'AM3 Aero',L,PATH,NAER
            endif
         enddo  ! M
!         write(*,*)'Called OPTICAM3'

         do K=1,5
            if (OD(K,L) .gt. 0.d0) then
               SSA(K,L) = SSA(K,L)/OD(K,L)
               do I=1,8
                  SLEG(I,K,L) = SLEG(I,K,L)/OD(K,L)
               enddo
            endif
         enddo

!---Include aerosol with cloud OD at 600 nm to determine added layers:
         OD600(L) = OD(4,L)

      enddo             ! end loop of L

!      write(*,*)'End of SLEG loop'
!---when combining with Rayleigh and O2-O3 abs, remember the SSA and
!---  phase fn SLEG are weighted by OD and OD*SSA, respectively.
!---Given the aerosol+cloud OD/layer in visible (600 nm) calculate how to add
!       additonal levels at top of clouds (now uses log spacing)

!-----------------------------------------------------------------------
      call EXTRAL(OD600,L1_,L2U,N_,JTAUMX,ATAU,ATAU0, JXTRA)
!-----------------------------------------------------------------------

!---set surface reflectance
        RFL(:) = max(0.d0,min(1.d0,RFLECT))

!      write(*,*)'Called EXTRAL and passed values to RFL'
!-----------------------------------------------------------------------
!---Loop over all wavelength bins to calc mean actinic flux AVGF(L)
!-----------------------------------------------------------------------

      do K = 1,W_
!         write(logunit,*) 'FL() => ', K, FL(K)
!         write(*,*) 'FL() => ', K, FL(K)
!      if (FL(K) .gt. 1.d0) then
         WAVE = WL(K)
!---Pick nearest Mie wavelength to get scattering properites------------
                                KMIE=1  ! use 200 nm prop for <255 nm
         if( WAVE .gt. 255.d0 ) KMIE=2  ! use 300 nm prop for 255-355 nm
         if( WAVE .gt. 355.d0 ) KMIE=3  ! use 400 nm prop for 355-500 nm
         if( WAVE .gt. 500.d0 ) KMIE=4
         if( WAVE .gt. 800.d0 ) KMIE=5


!---Combine: Rayleigh scatters & O2 & O3 absorbers to get optical properties
!---values at L1_=L_+1 are a pseudo/climatol layer above the top CTM layer (L_)
         do L = 1,L1_
            TTTX     = TTJ(L)
            call X_interp (TTTX,XQO2, TQQ(1,1),QO2(K,1), TQQ(2,1),QO2(K,2), &
                 TQQ(3,1),QO2(K,3), LQQ(1))
            call X_interp (TTTX,XQO3, TQQ(1,2),QO3(K,1), TQQ(2,2),QO3(K,2), &
                 TQQ(3,2),QO3(K,3), LQQ(2))
            ODABS = XQO3*OOJ(L) + XQO2*DDJ(L)*0.20948d0
            ODRAY = DDJ(L)*QRAYL(K)

            DTAUX(L,K) = OD(KMIE,L) + ODABS + ODRAY

            do I=1,8
               POMEGAX(I,L,K) = SLEG(I,KMIE,L)*OD(KMIE,L)
            enddo
               POMEGAX(1,L,K) = POMEGAX(1,L,K) + 1.0d0*ODRAY
               POMEGAX(3,L,K) = POMEGAX(3,L,K) + 0.5d0*ODRAY
            do I=1,8
               POMEGAX(I,L,K) = POMEGAX(I,L,K)/DTAUX(L,K)
            enddo
         enddo   ! L

!      endif
      enddo    ! K
!      write(*,*)'Calculated mean actinic flux AVGF'
!-----------------------------------------------------------------------

      call OPMIE (DTAUX,POMEGAX,U0,RFL,AMF2,JXTRA, &
              AVGF,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0, LU)

!       write(*,*)'Called OPMIE'
!-----------------------------------------------------------------------

!         FLXUP(:) = 0.d0
!         DIRUP(:) = 0.d0
!         FLXDN(:) = 0.d0
!         DIRDN(:) = 0.d0
!         FLXJ(:) = 0.d0
!         FFX(:,:) = 0.d0
!         FFXNET(:,:) = 0.d0

      do K = 1,W_
!      if (FL(K) .gt. 1.d0) then

!----direct(DIR) and diffuse(FLX) fluxes at top(UP) (solar = negative by convention)
!----     also at bottom (DN), does not include diffuse reflected flux.
!        FLXUP(K) =  FJTOP(K)
!        DIRUP(K) = -FLXD0(K)
!        FLXDN(K) = -FJBOT(K)
!        DIRDN(K) = -FSBOT(K)

        do L = 1,LU
          FFF(K,L) = FFF(K,L) + SOLF*FL(K)*AVGF(L,K)
!          write(*,*) K, L, FFF(K,L)
        enddo
!        FREFI = FREFI + SOLF*FL(K)*FLXD0(K)/WL(K)
!        FREFL = FREFL + SOLF*FL(K)*FJTOP(K)/WL(K)
!        FREFS = FREFS + SOLF*FL(K)/WL(K)

!---for each wavelength calculate the flux budget/heating rates:
!  FLXD(L) = direct flux deposited in layer L  [approx = MU0*(F(L+1) -F(L))]
!            but for spherical atmosphere!
!  FJFLX(L) = diffuse flux across top of layer L

!---calculate divergence of diffuse flux in each CTM layer (& t-o-a)
!---     need special fix at top and bottom:
!---FABOT = total abs at L.B. &  FXBOT = net diffusive flux at L.B.
!        FABOT = (1.d0-RFL(K))*(FJBOT(K)+FSBOT(K))
!        FXBOT = -FJBOT(K) + RFL(K)*(FJBOT(K)+FSBOT(K))
!        FLXJ(1) = FJFLX(1,K) - FXBOT
!        do L=2,LU
!          FLXJ(L) = FJFLX(L,K) - FJFLX(L-1,K)
!        enddo
!        FLXJ(LU+1) = FJTOP(K) - FJFLX(LU,K)
!---calculate net flux deposited in each CTM layer (direct & diffuse):
!        FFX0 = 0.d0
!        do L=1,L1_
!           FFX(K,L) = FLXD(L,K) - FLXJ(L)
!           FFX0 = FFX0 + FFX(K,L)
!        enddo

!  NB: the radiation level ABOVE the top CTM level is included in these budgets
!      these are the flux budget/heating terms for the column:
!  FFXNET(K,1) = FLXD0        direct(solar) flux dep into atmos (spherical)
!  FFXNET(K,2) = FSBOT        direct(solar) flux dep onto LB (surface)
!  FFXNET(K,3) = FLXD0+FSBOT  TOTAL solar into atmopshere+surface
!  FFXNET(K,4) = FJTOP        diffuse flux leaving top-of-atmos
!  FFXNET(K,5) = FFX0         diffuse flux absorbed in atmos
!  FFXNET(K,6) = FABOT        total (dir+dif) absorbed at LB (surface)
!       these are surface fluxes to compare direct vs. diffuse:
!  FFXNET(K,7) = FSBOT        direct flux dep onto LB (surface) - for srf diags
!  FFXNET(K,8) = FJBOT        diffuse flux dep onto LB (surface)

!        FFXNET(K,1) = FLXD0(K)
!        FFXNET(K,2) = FSBOT(K)
!        FFXNET(K,3) = FLXD0(K) + FSBOT(K)
!        FFXNET(K,4) = FJTOP(K)
!        FFXNET(K,5) = FFX0
!        FFXNET(K,6) = FABOT
!        FFXNET(K,7) = FSBOT(K)
!        FFXNET(K,8) = FJBOT(K)

!-----------------------------------------------------------------------
!      endif
      enddo       ! end loop over wavelength K
!-----------------------------------------------------------------------
!      FREFL = FREFL/FREFS      !calculate reflected flux (energy weighted)
!      FREFI = FREFI/FREFS

!       write(*,*)'Calculated FFF'
!---NB UVB = 280-320 = bins 12:15, UVA = 320-400 = bins 16:17, VIS = bin 18 (++)

!-----------------------------------------------------------------------
!      write(logunit,*) 'Check NJX ', NJX, ' NJX_ ', NJX_
!      write(*,*) 'Check NJX ', NJX, ' NJX_ ', NJX_
      call JRATET(pfull,tfull,FFF, VALJXX, LU,NJX_)
!-----------------------------------------------------------------------

!      write(*,*)'Called JRATET'
!---mapping J-values from fast-JX onto CTM chemistry 

      LPRTJ = .false.
      do L = 1,L_
         do J = 1,NRATJ
            if (JIND(J).gt.0) then
               ZPJ(L_+1-L,J) = VALJXX(L,JIND(J))*JFACTA(J)
               if(ZPJ(L_+1-L,J) .lt. 0.0d0) then
                 LPRTJ = .true.
               endif
            else
               ZPJ(L_+1-L,J) = 0.d0
            endif
         enddo
      enddo
!      write(*,*)'Mapped J-values'


!-----------------------------------------------------------------------
!      if (LPRTJ) then
!---diagnostics below are NOT returned to the CTM code
!         write(*,*)'fast-JX-(7.0)---PHOTO_JX internal print: Atmosphere---'
!---used last called values of DTAUX and POMEGAX, should be 600 nm
!         do L=1,L1_
!            DTAU600(L) = DTAUX(L,W_)
!            do I=1,8
!               POMG600(I,L) = POMEGAX(I,L,W_)
!            enddo
!         enddo
!
!         call JP_ATM(PPJ,TTJ,DDJ,OOJ,ZZJ,DTAU600,POMG600,JXTRA, LU)
!
!---PRINT SUMMARY of mean intensity, flux, heating rates:
!         write(*,*)
!         write(*,*)'fast-JX(7.0)---PHOTO_JX internal print: Mean Intens---'
!         write(*,'(a,5f10.4)') &
!                 ' SUMMARY fast-JX: albedo/SZA/u0/F-incd/F-refl/', &
!                 RFLECT,SZA,U0,FREFI,FREFL
!         write(*,'(a5,18i8)')   ' bin:',(K, K=NW2,NW1,-1)
!         write(*,'(a5,18f8.1)') ' wvl:',(WL(K), K=NW2,NW1,-1)
!         write(*,'(a)') ' ----  100000=Fsolar   MEAN INTENSITY per wvl bin'
!         RATIO(:) = 0.d0
!         do L = LU,1,-1
!            do K=NW1,NW2
!               if (FL(K) .gt. 1.d0) then
!                  RATIO(K) = (1.d5*FFF(K,L)/FL(K))
!               endif
!            enddo
!            write(*,'(i3,2x,18i8)') L,(RATIO(K),K=NW2,NW1,-1)
!         enddo
!
!         write(*,*)
!         write(*,*)'fast-JX(7.0)---PHOTO_JX internal print: Net Fluxes---'
!         write(*,'(a11,18i8)')   ' bin:',(K, K=NW2,NW1,-1)
!         write(*,'(a11,18f8.1)') ' wvl:',(WL(K), K=NW2,NW1,-1)
!!      write(6,'(a11,18f8.4)') ' sol in atm',(FFXNET(K,1), K=NW2,NW1,-1)
!!      write(6,'(a11,18f8.4)') ' sol at srf',(FFXNET(K,2), K=NW2,NW1,-1)
!         write(*,*) ' ---NET FLUXES--- '
!         write(*,'(a11,18f8.4)') ' sol TOTAL ',(FFXNET(K,3), K=NW2,NW1,-1)
!         write(*,'(a11,18f8.4)') ' dif outtop',(FFXNET(K,4), K=NW2,NW1,-1)
!         write(*,'(a11,18f8.4)') ' abs in atm',(FFXNET(K,5), K=NW2,NW1,-1)
!         write(*,'(a11,18f8.4)') ' abs at srf',(FFXNET(K,6), K=NW2,NW1,-1)
!         write(*,*) ' ---SRF FLUXES--- '
!         write(*,'(a11,18f8.4)') ' srf direct',(FFXNET(K,7), K=NW2,NW1,-1)
!         write(*,'(a11,18f8.4)') ' srf diffus',(FFXNET(K,8), K=NW2,NW1,-1)
!         write(*,'(2a)') '  ---NET ABS per layer:       10000=Fsolar', &
!              '  [NB: values <0 = numerical error w/clouds or SZA>90, colm OK]'
!         do L = LU,1,-1
!            do K=NW1,NW2
!               RATIO(K) = 1.d5*FFX(K,L)
!            enddo
!            write(*,'(i9,2x,18i8)') L,(RATIO(K),K=NW2,NW1,-1)
!         enddo
!         write(*,'(a)')
!         write(*,'(a)') ' fast-JX (7.0)----J-values----'
!         write(*,'(1x,a,72(a6,3x))') 'L=  ',(TITLEJX(K), K=1,NJX)
!         do L = LU,1,-1
!            write(*,'(i3,1p, 72e9.2)') L,(VALJXX(L,K),K=1,NJX)
!         enddo
!
!      endif

   99 continue
   
         
      return
      end subroutine fastjx_photo                                          
                                                                        
!<<<<<<<<<<<<<<<<<<<<<<<end CTM-fastJX linking subroutines<<<<<<<<<<<<<<
 

      subroutine fastjx_end
         implicit none
         module_is_initialized = .false.
      deallocate (solflxtot)
      deallocate (solflxband_all)
      deallocate (solflxband_ann_1850)
      deallocate (solflxband_ann_2300)
      end subroutine fastjx_end

!<<<<<<<<<<<<<<<<<<<<<<<<fast-J initialization subroutines<<<<<<<<<<<<<<

      subroutine RD_XXX(NAMFIL)
!-----------------------------------------------------------------------
!  Read in wavelength bins, solar fluxes, Rayleigh, T-dep X-sections.
!
!>>>>NEW v-6.8  now allow 1 to 3 sets of X-sects for T or P
!           LQQ = 1, 2, or 3 to determine interpolation with T or P
!           IF the temperatures TQQQ are <0, then use as pressure interp (hPa)
!           NB - the temperatures and pressures must be increasing
!>>>>NEW v-6.4  changed to collapse wavelengths and x-sections to Trop-only:
!           WX_ = 18 should match the JX_spec.dat wavelengths
!           W_ = 12 (Trop-only) or 18 (std) is set in (cmn_FXJ.f).
!       if W_=12 then drop strat wavels, and drop x-sects (e.g. N2O, ...)
!           W_ = 8, reverts to quick fix:  fast-J (12-18) plus bin (5) scaled
!
!-----------------------------------------------------------------------
!     NAMFIL   Name of spectral data file (JX_spec.dat) >> j2 for fast-J2
!     NUN      Channel number for reading data file
!
!     NJX    Number of species to calculate J-values for
!     NWWW     Number of wavelength bins, from 1:NWWW
!     WBIN     Boundaries of wavelength bins
!     WL       Centres of wavelength bins - 'effective wavelength'
!     FL       Solar flux incident on top of atmosphere (cm-2.s-1)
!     QRAYL    Rayleigh parameters (effective cross-section) (cm2)
!     QO2      O2 cross-sections
!     QO3      O3 cross-sections
!     Q1D      O3 => O(1D) quantum yield
!     TQQ      Temperature for supplied cross sections
!     QQQ      Supplied cross sections in each wavelength bin (cm2)
!-----------------------------------------------------------------------
      implicit none

!      integer, intent(in) :: NUN
      character(*), intent(in) ::  NAMFIL

      integer  I, J, JJ, K, IW, NQRD, NWWW,   LQ, NUN
      real*8  QQ2(199), TQQ2

      character*78 TITLE0
      character*6  TITLEJ2,TITLEJ3
      character*1  TSTRAT

      TQQ(:,:) = 0.d0

!----------spectral data----set for new format data------------------
!         note that X_ = max # Xsects read in
!                   NJX = # fast-JX J-values derived from this (.le. X_)

! >>>> W_ = 12 <<<< means trop-only, discard WL #1-4 and #9-10, some X-sects

!      open (NUN,FILE=NAMFIL,status='old',form='formatted')
      call mpp_open (NUN, trim(NAMFIL), MPP_RDONLY, MPP_ASCII,  &
                     MPP_SEQUENTIAL, MPP_MULTI, MPP_SINGLE)

      read (NUN,100) TITLE0

!----note that NQRD is not used any more, a read until 'endofJ' is performed
      read (NUN,101) NQRD,NWWW
         NW1 = 1
         NW2 = NWWW
!      write(*,'(1x,a)') TITLE0
!      write(*,'(i8)') NWWW
!----J-values:  1=O2, 2=O3P,3=O3D 4=readin Xsects
      read (NUN,102) (WL(IW),IW=1,NWWW)
      read (NUN,102) (FL(IW),IW=1,NWWW)
      read (NUN,102) (QRAYL(IW),IW=1,NWWW)

!---Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps)
!---NB the O3 and q-O3-O1D are at different temperatures and cannot be combined
      read (NUN,103) TITLEJX(1),TQQ(1,1), (QO2(IW,1),IW=1,NWWW)
      read (NUN,103) TITLEJ2,  TQQ(2,1), (QO2(IW,2),IW=1,NWWW)
      read (NUN,103) TITLEJ3,  TQQ(3,1), (QO2(IW,3),IW=1,NWWW)

      read (NUN,103) TITLEJX(2),TQQ(1,2), (QO3(IW,1),IW=1,NWWW)
      read (NUN,103) TITLEJ2,  TQQ(2,2), (QO3(IW,2),IW=1,NWWW)
      read (NUN,103) TITLEJ3,  TQQ(3,2), (QO3(IW,3),IW=1,NWWW)

      read (NUN,103) TITLEJX(3),TQQ(1,3), (Q1D(IW,1),IW=1,NWWW)
      read (NUN,103) TITLEJ2,  TQQ(2,3), (Q1D(IW,2),IW=1,NWWW)
      read (NUN,103) TITLEJ3,  TQQ(3,3), (Q1D(IW,3),IW=1,NWWW)

      SQQ(1) = ' '
      SQQ(2) = ' '
      SQQ(3) = ' '

      LQQ(1) = 3
      LQQ(2) = 3
      LQQ(3) = 3

!---Read remaining species:  X-sections at 1-2-3 T_s
        JJ = 3
      do I=1,9999

!--try to read in 3 X-sects per J-value (JJ)
        read (NUN,104) TITLEJ2,TSTRAT,TQQ2,(QQ2(IW),IW=1,NWWW)
        if (TITLEJ2 .eq. 'endofJ') goto 1
!---skip stratosphere only J's (denoted by 'x')if W_<18 => trop-only J's
        if (W_.eq.18 .or. TSTRAT.ne.'x') then
          if (TITLEJ2 .ne. TITLEJX(JJ)) then
               JJ = JJ+1

         if (JJ .gt. X_) then
            write(*,*)' RD_XXX:  X_ not large enough for Xsects read in', JJ, X_
            stop
         endif

               TITLEJX(JJ) = TITLEJ2
               LQQ(JJ) = 1
               SQQ(JJ) = TSTRAT
                 LQ = LQQ(JJ)
               TQQ(LQ,JJ) = TQQ2
             do IW = 1,NWWW
               QQQ(IW,LQ,JJ) = QQ2(IW)
             enddo
           else
               LQQ(JJ) = LQQ(JJ)+1
            if (LQQ(JJ) .le. 3) then
               LQ = LQQ(JJ)
               TQQ(LQ,JJ) = TQQ2
             do IW = 1,NWWW
               QQQ(IW,LQ,JJ) = QQ2(IW)
             enddo
            endif
           endif
        endif
      enddo
    1 continue
        NJX = JJ

      do J = 1,NJX
!        write(*,200) J,TITLEJX(J),SQQ(J),LQQ(J),(TQQ(I,J),I=1,LQQ(J))
!---need to check that TQQ is monotonically increasing:
        if (LQQ(J) .eq. 3) then
          if (TQQ(2,J) .ge. TQQ(3,J)) then
             write(*,*) 'TQQ out of order ', j,  TQQ(2, J), TQQ(3,J)
             stop
          endif
          if (TQQ(1,J) .ge. TQQ(2,J)) then
             write(*,*) 'TQQ out of order', j, TQQ(1,J), TQQ(2,J)
          endif
        endif
        if (LQQ(J) .eq. 2) then
          if (TQQ(1,J) .ge. TQQ(2,J)) then
            write(*,*) 'TQQ out of order',j, TQQ(1,J), TQQ(2,j)
          endif
        endif
      enddo

!---check on doingpressure interp
!---check on consolidating Qo2 and others into
!---wrte a newFJX_J2J.dat for mapping on fjx Xsects


!---truncate number of wavelengths to do troposphere-only
      if (W_ .ne. WX_) then
!---TROP-ONLY
       if (W_ .eq. 12) then
!        write(*,'(a)')  &
!         ' >>>TROP-ONLY reduce wavelengths to 12, drop strat X-sects'
        NW2 = 12
        do IW = 1,4
          WL(IW) = WL(IW+4)
          FL(IW) = FL(IW+4)
          QRAYL(IW) = QRAYL(IW+4)
         do K = 1,3
          QO2(IW,K) = QO2(IW+4,K)
          QO3(IW,K) = QO3(IW+4,K)
          Q1D(IW,K) = Q1D(IW+4,K)
         enddo
         do J = 4,NJX
!          QQQ(IW,1,J) = QQQ(IW+4,1,J)
!          QQQ(IW,2,J) = QQQ(IW+4,2,J)
! j2l follow GeosChem v10-01c fast_jx_mod.F
           do LQ=1,LQQ(J)
              QQQ(IW,LQ,J) = QQQ(IW+4,LQ,J)
           enddo
         enddo
        enddo
        do IW = 5,12
          WL(IW) = WL(IW+6)
          FL(IW) = FL(IW+6)
          QRAYL(IW) = QRAYL(IW+6)
         do K = 1,3
          QO2(IW,K) = QO2(IW+6,K)
          QO3(IW,K) = QO3(IW+6,K)
          Q1D(IW,K) = Q1D(IW+6,K)
         enddo
         do J = 4,NJX
!          QQQ(IW,1,J) = QQQ(IW+6,1,J)
!          QQQ(IW,2,J) = QQQ(IW+6,2,J)
! j2l
           do LQ=1,LQQ(J)
              QQQ(IW,LQ,J) = QQQ(IW+6,LQ,J)
           enddo
         enddo
        enddo
!---TROP-QUICK  (must scale solar flux for W=5)
       elseif (W_ .eq. 8) then
!        write(*,'(a)')   &
!         ' >>>TROP-QUICK reduce wavelengths to 8, drop strat X-sects'
        NW2 = 8
        do IW = 1,1
          WL(IW) = WL(IW+4)
          FL(IW) = FL(IW+4)  * 2.d0
          QRAYL(IW) = QRAYL(IW+4)
         do K = 1,3
          QO2(IW,K) = QO2(IW+4,K)
          QO3(IW,K) = QO3(IW+4,K)
          Q1D(IW,K) = Q1D(IW+4,K)
         enddo
         do J = 4,NJX
!          QQQ(IW,1,J) = QQQ(IW+4,1,J)
!          QQQ(IW,2,J) = QQQ(IW+4,2,J)
! j2l
           do LQ=1,LQQ(J)
              QQQ(IW,LQ,J) = QQQ(IW+4,LQ,J)
           enddo
         enddo
        enddo
        do IW = 2,8
          WL(IW) = WL(IW+10)
          FL(IW) = FL(IW+10)
          QRAYL(IW) = QRAYL(IW+10)
         do K = 1,3
          QO2(IW,K) = QO2(IW+10,K)
          QO3(IW,K) = QO3(IW+10,K)
          Q1D(IW,K) = Q1D(IW+10,K)
         enddo
         do J = 4,NJX
!          QQQ(IW,1,J) = QQQ(IW+10,1,J)
!          QQQ(IW,2,J) = QQQ(IW+10,2,J)
! j2l
           do LQ=1,LQQ(J)
              QQQ(IW,LQ,J) = QQQ(IW+10,LQ,J)
           enddo
         enddo
        enddo

       else
         write(*,*)' no. wavelengths wrong: W_ .ne. 8,12,18', W_
         stop
       endif
      endif

      call mpp_close (NUN)

  100 format(a)
  101 format(10x,5i5)
  102 format(10x,    6e10.3/(10x,6e10.3)/(10x,6e10.3))
  103 format(a6,1x,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))
  104 format(a6,a1,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))
  200 format(1x,' x-sect:',i3,a10,a4,i5,3(3x,f6.2))
  201 format(' Number of x-sections supplied to Fast-J2: ',i3,/,    &
             ' Maximum number allowed (X_) only set to: ',i3,       &
             ' - increase in cmn_FJX.f')
      return
      END SUBROUTINE RD_XXX
                                                                        
!-----------------------------------------------------------------------
      subroutine RD_CLD(NAMFIL)
!-----------------------------------------------------------------------
!-------aerosols/cloud scattering data set for fast-JX (ver 5.3+)
!  >>>>>>>>>>>>>>>>spectral data rev to J-ref ver8.5 (5/05)<<<<<<<<<<<<
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
!     NUN      Channel number for reading data file
!     NCC      Number of categories for cloud scattering phase functions
!     QCC      Cloud scattering phase functions
!     WCC      5 Wavelengths for supplied phase functions
!     PCC      Phase function: first 8 terms of expansion
!     RCC      Effective radius associated with cloud type
!     SCC      Single scattering albedo
!     DCC      density (g/cm^3)
!-----------------------------------------------------------------------
      implicit none

      character(*), intent(in) ::  NAMFIL

      integer  I, J, K, NUN
      character*78 TITLE0
      character*20 TITLAA(A_)   ! TITLAA: Title for scatering data

!      open (NUN,FILE=NAMFIL,status='old',form='formatted')
      call mpp_open (NUN, trim(NAMFIL), MPP_RDONLY, MPP_ASCII,  &
                     MPP_SEQUENTIAL, MPP_MULTI, MPP_SINGLE)


      read (NUN,'(i2,a78)') NCC,TITLE0
        if (NCC .gt. C_) then
          write(*,*)' too many cld-data sets: NCC > C_'
          stop
        endif

!      write(*,'(a,2f9.5,i5)') ' ATAU/ATAU0/JMX',ATAU,ATAU0,JTAUMX

      read (NUN,*)
      read (NUN,*)
      do J = 1,NCC
          read (NUN,'(3x,a8,1x,2f6.3)') TITLAA(J),RCC(J),DCC(J)
        do K = 1,5
          read (NUN,'(f4.0,f7.4,f7.4,7f6.3)')     &
        WCC(K,J),QCC(K,J),SCC(K,J),(PCC(I,K,J),I=2,8)
          PCC(1,K,J) = 1.d0
        enddo
      enddo

      call mpp_close(NUN)

!      write(*,'(a,9f8.1)') ' Aerosol optical: r-eff/rho/Q(@wavel):'  &
!                   ,(WCC(K,1),K=1,5)
!      write(*,*) TITLE0
!      do J=1,NCC
!      write(*,'(i3,1x,a8,7f8.3)')   &
!                    J,TITLAA(J),RCC(J),DCC(J),(QCC(K,J),K=1,5)
!      enddo
      return
      END SUBROUTINE RD_CLD

!-----------------------------------------------------------------------
      subroutine RD_MIE(NAMFIL)
!-----------------------------------------------------------------------
!-------aerosols scattering data set for fast-JX (ver 5.3+)
!  >>>>>>>>>>>>>>>>spectral data rev to J-ref ver8.5 (5/05)<<<<<<<<<<<<
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
!     NUN      Channel number for reading data file
!     NAA      Number of categories for scattering phase functions
!     QAA      Aerosol scattering phase functions
!     WAA      5 Wavelengths for the supplied phase functions
!     PAA      Phase function: first 8 terms of expansion
!     RAA      Effective radius associated with aerosol type
!     SAA      Single scattering albedo
!     DAA      density (g/cm^3)
!-----------------------------------------------------------------------
      implicit none

!      logical, intent(in)                   :: LPRTJ                   ! set to false
      character(*), intent(in) ::  NAMFIL

      integer  I, J, K, NUN
      character*78 TITLE0
      character*20 TITLAA(A_)   ! TITLAA: Title for scatering data

!      open (NUN,FILE=NAMFIL,status='old',form='formatted')
      call mpp_open (NUN, trim(NAMFIL), MPP_RDONLY, MPP_ASCII,  &
                     MPP_SEQUENTIAL, MPP_MULTI, MPP_SINGLE)

      read (NUN,'(i2,a78)') NAA,TITLE0
        if (NAA .gt. A_) then
          write(*,*)' too many aerosol-data sets: NAA > A_'
          stop
        endif

!      write(*,'(a,2f9.5,i5)') ' ATAU/ATAU0/JMX',ATAU,ATAU0,JTAUMX

      read (NUN,*)
      read (NUN,*)
      do J = 1,NAA
          read (NUN,'(3x,a8,1x,2f6.3)') TITLAA(J),RAA(J),DAA(J)
        do K = 1,5
          read (NUN,'(f4.0,f7.4,f7.4,7f6.3)')    &
        WAA(K,J),QAA(K,J),SAA(K,J),(PAA(I,K,J),I=2,8)
          PAA(1,K,J) = 1.d0
        enddo
      enddo

      call mpp_close(NUN)

!      write(*,'(a,9f8.1)') ' Aerosol optical: r-eff/rho/Q(@wavel):'  &
!                   ,(WAA(K,1),K=1,5)
!      write(*,*) TITLE0
!      do J=1,NAA
!      write(*,'(i3,1x,a8,7f8.3)')   &
!                    J,TITLAA(J),RAA(J),DAA(J),(QAA(K,J),K=1,5)
!      enddo
      return
      END SUBROUTINE RD_MIE
!-----------------------------------------------------------------------
      subroutine RD_UM(NAMFIL)
!-----------------------------------------------------------------------
!-------UMich aerosol optical data for fast-JX (ver 6.1+)
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
!     NUN      Channel number for reading data file
!-----------------------------------------------------------------------
      implicit none

      character(*), intent(in) ::  NAMFIL

      integer  I, J, K, L, NUN
      character*78 TITLE0
      character*20 TITLUM(33)   ! TITLUM: Title for U Michigan aerosol data set

!      open (NUN,FILE=NAMFIL,status='old',form='formatted')
     call mpp_open (NUN, trim(NAMFIL), MPP_RDONLY, MPP_ASCII,  &
                     MPP_SEQUENTIAL, MPP_MULTI, MPP_SINGLE)

      read (NUN,'(a78)') TITLE0
!        write(*,*) 'UMichigan Aerosols', TITLE0
      read(NUN,'(5x,10f5.0)') WMM
!        write(*,'(a,10f7.1)') ' UMIchigan aerosol wavelengths:',WMM

!---33 Different UM Aerosol Types:  SULF, SS-1,-2,-3,-4, DD-1,-2,-3,-4,
!---      FF00(0%BC), FF02, ...FF14(14%BC),  BB00, BB02, ...BB30(30%BC)
      do L=1,33
          read(NUN,'(a4)') TITLUM(L)
!---21 Rel Hum:    K=1=0%, =2=5%, ... =20=95%, =21=99%
        do K=1,21
!---6 wavelengths: J=1=200nm, 2=300nm, 3=400nm, (4'=550nm) 5=600nm, 6=1000nm
!---3 optic vars:  I=1=SSAlbedo,  =2=g,  =3=k-ext
          read(NUN,'(18f9.5)')  ((UMAER(I,J,K,L),I=1,3),J=1,6)
        enddo
      enddo

      call mpp_close(NUN)

!        write(*,'(a)') 'collapse UM wavelengths, drop 550 nm'
          WMM(4) = WMM(5)
          WMM(5) = WMM(6)
       do L=1,33
       do K=1,21
       do I=1,3
          UMAER(I,4,K,L) = UMAER(I,5,K,L)
          UMAER(I,5,K,L) = UMAER(I,6,K,L)
       enddo
       enddo
       enddo

!        write(*,'(7(i5,1x,a4))') (L,TITLUM(L), L=1,33)

      return
      END SUBROUTINE RD_UM

!-----------------------------------------------------------------------
!      subroutine RD_PROF(NJ2,NAMFIL)
!-----------------------------------------------------------------------
!  Routine to input T and O3 reference profiles 'atmos_std.dat'
!-----------------------------------------------------------------------
!      implicit none
!
!      integer, intent(in) ::  NJ2
!      character(*), intent(in) ::  NAMFIL
!
!      integer IA, I, M, L, LAT, MON, NTLATS, NTMONS, N216
!      real*8  OFAC, OFAK
!
!      character*78 TITLE0
!
!      open (NJ2,file=NAMFIL,status='old',form='formatted')
!      read (NJ2,'(A)') TITLE0
!      read (NJ2,'(2I5)') NTLATS,NTMONS
!      write(6,'(1X,A)') TITLE0
!      write(6,1000) NTLATS,NTMONS
!      N216  = min(216, NTLATS*NTMONS)
!      do IA = 1,N216
!        read (NJ2,'(1X,I3,3X,I2)') LAT, MON
!        M = min(12, max(1, MON))
!        L = min(18, max(1, (LAT+95)/10))
!        read (NJ2,'(3X,11F7.1)') (TREF(I,L,M), I=1,41)
!        read (NJ2,'(3X,11F7.4)') (OREF(I,L,M), I=1,31)
!      enddo
!      close (NJ2)
!
!  Extend climatology to 100 km
!      OFAC = exp(-2.d5/5.d5)
!      do I = 32,51
!        OFAK = OFAC**(I-31)
!        do M = 1,NTMONS
!        do L = 1,NTLATS
!          OREF(I,L,M) = OREF(31,L,M)*OFAK
!        enddo
!        enddo
!      enddo
!      do L = 1,NTLATS
!      do M = 1,NTMONS
!      do I = 42,51
!        TREF(I,L,M) = TREF(41,L,M)
!      enddo
!      enddo
!      enddo
!
! 1000 format(1x,'std atmos profiles: ',i3,' lat x ',i2,' mon')
!
!      END SUBROUTINE RD_PROF

!-----------------------------------------------------------------------
      subroutine RD_JS_JX(NAMFIL,TITLEJX,NJX)
!-----------------------------------------------------------------------
!  Read 'FJX_j2j.dat' that defines mapping of fast-JX J's (TITLEJX(1:NJX))
!    onto the CTM reactions:  react# JJ, named T_REACT, uses fast-JX's T_FJX
!    including scaling factor F_FJX.
!-----------------------------------------------------------------------
!---mapping variables stored in  block /jvchem/JFACTA,JIND,NRATJ,JLABEL,JMAP
!           real*8  JFACTA(JVN_)          integer JIND(JVN_), NRATJ
!           character*50 JLABEL(JVN_)     character*6  JMAP(JVN_)
!     JFACTA    multiplication factor for fast-JX calculated J
!     JLABEL    label(*50) of J-value used in the main chem model
!     JMAP      label(*6) of J-value used to match with fast-JX J's
!     NRATJ     number of Photolysis reactions in CTM chemistry, derived here
!                   NRATJ must be .le. JVN_
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in)                    ::  NJX
      character(*), intent(in)               ::  NAMFIL
      character*6, intent(in),dimension(NJX) :: TITLEJX
      integer   J,JJ,K,NUNIT
      character*120 CLINE
      character*50 T_REACT
      character*6  JMAP(JVN_), T_FJX
      real*8 F_FJX

! Read the FJX_j2j.dat file to map model specific J's onto fast-JX J's
! The chemistry code title describes fully the reaction (a50)
! Blank (unfilled) chemistry J's are unmapped
! The number NRATJ is the last JJ readin that is .le. JVN
!   include fractional quantum yield for the fast-JX J's

      JLABEL(:) = '------'
      JMAP(:)   = '------'
      JFACTA(:) = 0.d0

!      open (NUNIT,file=NAMFIL,status='old',form='formatted')
      call mpp_open (NUNIT, trim(NAMFIL), MPP_RDONLY, MPP_ASCII,  &
                     MPP_SEQUENTIAL, MPP_MULTI, MPP_SINGLE)


       read (NUNIT,'(a)') CLINE
!         write(*,'(a)') CLINE
      do J = 1,JVN_
       read (NUNIT,'(i4,1x,a50,4x,f5.3,2x,a6)') JJ,T_REACT,F_FJX,T_FJX
       if (JJ.gt.JVN_) goto 20
        JLABEL(JJ) = T_REACT
        JFACTA(JJ) = F_FJX
        JMAP(JJ) = T_FJX
        NRATJ = JJ
      enddo

 20   call mpp_close(NUNIT)

!---Zero / Set index arrays that map Jvalue(j) onto rates
      do K = 1,NRATJ
         JIND(K) = 0
       do J = 1,NJX
         T_FJX = TITLEJX(J)
        if (JMAP(K) .eq. TITLEJX(J)) then
         JIND(K)=J
        endif
       enddo
      enddo

!      write(*,'(a,i4,a)')' Photochemistry Scheme with',NRATJ,' J-values'
!      do K=1,NRATJ
!       if (JMAP(K) .ne. '------' ) then
!        J = JIND(K)
!        if (J.eq.0) then
!         write(*,'(i5,a50,f6.3,a,1x,a6)') K,JLABEL(K),JFACTA(K),      &
!               ' no mapping onto fast-JX',JMAP(K)
!        else
!         write(*,'(i5,a50,f6.3,a,i4,1x,a6)') K,JLABEL(K),JFACTA(K),   &
!               ' mapped to FJX:',J,TITLEJX(J)
!        endif
!       endif
!      enddo

      return
      END SUBROUTINE RD_JS_JX

      subroutine RD_SOLAR(NAMFIL)
      character(*), intent(in)               ::  NAMFIL
      integer   NUN, year, month,nyr,nv,nband

      call mpp_open (NUN, trim(NAMFIL), MPP_RDONLY, MPP_ASCII,  &
                     MPP_SEQUENTIAL, MPP_MULTI, MPP_SINGLE)
          read (NUN, FMT = '(4i8)') first_yr, last_yr,  &
                                   nvalues_per_year, numbands
          if (numbands /= nbands) then
            call error_mesg ('fastjx_mod', &
            ' number of wavelength bins is incorrect',&
                                                           FATAL)
          endif

          years_of_data  = last_yr  - first_yr  + 1
          allocate (solflxtot (years_of_data , nvalues_per_year))
          allocate (solflxband_all(years_of_data, nvalues_per_year, numbands))
          allocate (solflxband_ann_1850(numbands))
          allocate (solflxband_ann_2300(numbands))

          read (NUN, FMT = '(2i6,e12.5)') year, month, solflxtot_ann_1850
          read (NUN, FMT = '(6e12.5 )')  (solflxband_ann_1850(nband), nband =1,numbands)

          do nyr=1,years_of_data
            do nv=1,nvalues_per_year
              read (NUN, FMT = '(2i6,e12.5)') year, month, solflxtot(nyr,nv)
              read (NUN, FMT = '(6e12.5 )')  (solflxband_all  &
                                (nyr,nv,nband), nband =1,numbands)
            end do
          end do

          read (NUN, FMT = '(2i6,e12.5)') year, month, solflxtot_ann_2300
          read (NUN, FMT = '(6e12.5 )')  (solflxband_ann_2300(nband), nband =1,numbands)
      call mpp_close(NUN)

      END SUBROUTINE RD_SOLAR                                           

      subroutine get_solar_data (year, month, solar_constant, solflxband_now)
      
      integer,            intent(in)  :: year, month
      real,               intent(out) :: solar_constant
      real, dimension(:), intent(out) :: solflxband_now
      
      integer :: nband
      !--------------------------------------------------------------------
      !    returns solar constant for year and month
      !--------------------------------------------------------------------
      
      !   first some error checks
      
            if (years_of_data == 0) then
                call error_mesg ('solar_data_mod', &
                   'data may not exist', FATAL)
            endif
      
            if (size(solflxband_all,3) /= numbands) then
              call error_mesg ('solar_data_mod', &
                    'bands present in solar constant time data differs from &
                     &fastjx band number', FATAL)
            endif
      
      !--------------------------------------------------------------------
      
            if (year < first_yr) then
               solar_constant = solflxtot_ann_1850
               do nband=1,numbands
                  solflxband_now(nband) = solflxband_ann_1850(nband)
               end do
            else if (year > last_yr) then
               solar_constant = solflxtot_ann_2300
               do nband=1,numbands
                  solflxband_now(nband) = solflxband_ann_2300(nband)
               end do
            else
               solar_constant = solflxtot(year-first_yr+1, month)
               do nband=1,numbands
                  solflxband_now(nband) = solflxband_all(year-first_yr+1, month, nband)
               end do
            endif

      !--------------------------------------------------------------------
      
      end subroutine get_solar_data

!-----------------------------------------------------------------------
!      subroutine EXITC(T_EXIT)
!!-----------------------------------------------------------------------
!      character(len=*), intent(in) ::  T_EXIT
!
!      write(*,'(a)') T_EXIT
!      stop
!
!      END SUBROUTINE EXITC
           
!<<<<<<<<<<<<<<<<<<<<<<<<end fast-J initialization subroutines<<<<<<<<<<<<<<<<<<<
                                                                        
                                                                        
!<<<<<<<<<<<<<<<<<<<<<begin core fast-J special call subroutines<<<<<<<<<<<<<<<<<<<<<
                                                                        
!-----------------------------------------------------------------------
!      subroutine SOLAR_JX(GMTIME,NDAY,YGRDJ,XGRDI, SZA,COSSZA,SOLFX)
!-----------------------------------------------------------------------
!     GMTIME = UT for when J-values are wanted
!           (for implicit solver this is at the end of the time step)
!     NDAY   = integer day of the year (used for solar lat and declin)
!     YGRDJ  = laitude (radians) for grid (I,J)
!     XGDRI  = longitude (radians) for grid (I,J)
!
!     SZA = solar zenith angle in degrees
!     COSSZA = U0 = cos(SZA)
!-----------------------------------------------------------------------
!      implicit none
!      real*8, intent(in) ::   GMTIME,YGRDJ,XGRDI
!      integer, intent(in) ::  NDAY
!      real*8, intent(out) ::  SZA,COSSZA,SOLFX
!
!      real*8  PI, PI180, LOCT
!      real*8  SINDEC, SOLDEK, COSDEC, SINLAT, SOLLAT, COSLAT, COSZ
!
!      PI     = 3.141592653589793d0
!      PI180  = PI/180.d0
!      SINDEC = 0.3978d0*sin(0.9863d0*(dble(NDAY)-80.d0)*PI180)
!      SOLDEK = asin(SINDEC)
!      COSDEC = cos(SOLDEK)
!      SINLAT = sin(YGRDJ)
!      SOLLAT = asin(SINLAT)
!      COSLAT = cos(SOLLAT)
!
!      LOCT   = (((GMTIME)*15.d0)-180.d0)*PI180 + XGRDI
!      COSSZA = COSDEC*COSLAT*cos(LOCT) + SINDEC*SINLAT
!      SZA    = acos(COSSZA)/PI180
!
!      write(6,*) ' XGRDI,YGRDJ',XGRDI,YGRDJ
!      write(6,*) ' LOCT (rad)',LOCT
!      write(6,*) ' SINDEC,COSDEC', SINDEC,COSDEC
!      write(6,*) ' SINLAT,COSLAT', SINLAT,COSLAT
!      write(6,*) ' COS, SZA',COSSZA,SZA
!
!      SOLFX  = 1.d0-(0.034d0*cos(dble(NDAY-186)*2.d0*PI/365.d0))
!
!      END SUBROUTINE SOLAR_JX

!-----------------------------------------------------------------------
!      subroutine ACLIM_FJX (YLATD, MONTH, PPP,TTT,ZZZ,DDD,OOO, L1U)
!-----------------------------------------------------------------------
!  Load fast-JX climatology for latitude & month given pressure grid
!-----------------------------------------------------------------------
!      implicit none
!
!      real*8, intent(in)  :: YLATD
!      integer, intent(in) :: MONTH, L1U
!      real*8, intent(in),  dimension(L1U+1) :: PPP
!      real*8, intent(out), dimension(L1U+1) :: ZZZ
!      real*8, intent(out), dimension(L1U) :: TTT,DDD,OOO
!
!      real*8, dimension(51) :: OREF2,TREF2
!      real*8, dimension(52) :: PSTD
!      integer  K, L, M, N
!      real*8   DLOGP,F0,T0,PB,PC,XC,MASFAC,SCALEH
!
!  Select appropriate month
!      M = max(1,min(12,MONTH))
!  Select appropriate latitudinal profiles
!      N = max(1, min(18, (int(YLATD+99)/10 )))
!      do K = 1,51
!        OREF2(K) = OREF(K,N,M)
!        TREF2(K) = TREF(K,N,M)
!      enddo
!
!  Apportion O3 and T on supplied climatology z levels onto CTM levels +1
!  with mass (pressure) weighting, assuming constant mixing ratio and
!  temperature half a layer on either side of the point supplied.
!   PPP(L=1:L1_)=edge-pressure of CTM layer, PPP(L1_+1)=0 (top-of-atmos)
!  Mass factor - delta-Pressure (mbars) to delta-Column (molecules.cm-2)
!      MASFAC = 100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)

!  Set up pressure levels for O3/T climatology - assume that value
!  given for each 2 km z* level applies from 1 km below to 1 km above,
!  so select pressures at these boundaries. Surface level values at
!  1000 mb are assumed to extend down to the actual PSURF (if > 1000)
!           PSTD(1) = max(PPP(1),1000.d0)
!           PSTD(2) = 1000.d0 * 10.d0**(-1.d0/16.d0)
!           DLOGP   = 10.d0**(-2.d0/16.d0)
!      do K = 3,51
!        PSTD(K) = PSTD(K-1)*DLOGP
!      enddo
!        PSTD(52)  = 0.d0
!      do L = 1,L1U
!        F0 = 0.d0
!        T0 = 0.d0
!        do K = 1,51
!          PC   = min(PPP(L),PSTD(K))
!          PB   = max(PPP(L+1),PSTD(K+1))
!          if (PC .gt. PB) then
!            XC = (PC-PB)/(PPP(L)-PPP(L+1))
!            F0 = F0 + OREF2(K)*XC
!            T0 = T0 + TREF2(K)*XC
!          endif
!        enddo
!        TTT(L)  = T0
!        DDD(L)  = (PPP(L)-PPP(L+1))*MASFAC
!        OOO(L) = F0*1.d-6*DDD(L)
!      enddo
!
!  Calculate effective altitudes using scale height at each level
!        ZZZ(1) = 0.d0
!      do L = 1,L1U-1
!        SCALEH      = 1.3806d-19*MASFAC*TTT(L)
!        ZZZ(L+1) = ZZZ(L) -( LOG(PPP(L+1)/PPP(L)) * SCALEH )
!      enddo
!        ZZZ(L1U+1) = ZZZ(L1U) + ZZHT
!
!      END SUBROUTINE ACLIM_FJX
!
!<<<<<<<<<<<<<<<<<<<end CTM-fastJX special call subroutines<<<<<<<<<<<<<

!<<<<<<<<<<<<<<<<<<<<<<<begin core scattering subroutines<<<<<<<<<<<<<<<
!-----------------------------------------------------------------------
      subroutine OPMIE (DTAUX,POMEGAX,U0,RFL,AMF2,JXTRA, &
              FJACT,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0, LU)
!-----------------------------------------------------------------------
      implicit none

      real*8, intent(in) ::   DTAUX(JXL1_,W_),POMEGAX(8,JXL1_,W_)
      real*8, intent(in) ::   AMF2(2*JXL1_+1,2*JXL1_+1)
      real*8, intent(in) ::   U0,RFL(W_)
      integer, intent(in) ::  JXTRA(JXL2_+1), LU
      real*8, intent(out) ::FJACT(JXL_,W_),FJTOP(W_),FJBOT(W_),FSBOT(W_)
      real*8, intent(out) ::  FJFLX(JXL_,W_),FLXD(JXL1_,W_),FLXD0(W_)

      integer JNDLEV(JXL_),JNELEV(JXL1_)
      integer JADDLV(JXL2_+1),JADDTO(JXL2_+1),L2LEV(JXL2_+1)
      integer JTOTL,I,II,J,K,L,LL,IX,JK,   L2,L2L,L22,LZ,LZZ,ND
      integer L1U,L2U,   LZ0,LZ1,LZMID
      real*8   SUMT,SUMJ

      real*8  DTAU(JXL1_+1,W_),POMEGAJ(M2_,JXL2_+1,W_),TTAU(JXL2_+1,W_)
      real*8  FTAU2(JXL2_+1,W_),POMEGAB(M2_,W_)
      real*8  ATAUA,ATAUZ,XLTAU,TAUDN,TAUUP,DTAUJ,FJFLX0
      real*8, dimension(W_) :: TAUBTM,TAUTOP,FBTM,FTOP,ZFLUX
!--- variables used in mie code-----------------------------------------
      real*8, dimension(W_)         :: FJT,FJB
      real*8, dimension(N_,W_)      :: FJ,FZ,ZTAU
      real*8, dimension(M2_,N_,W_)  :: POMEGA
      real*8, dimension(2*JXL1_,W_)  :: FLXD2

!---there is a parallel correspondence:
!  dimension of JX arrays JXL_ .ge. dimension that CTM is using = L_
!  but calculation is done for L_=LU, L1_=L1U, L2_=L2U lengths of CTM
!
!  fast-J Mie code for J_s, only uses 8-term expansion, 4-Gauss pts
!
! in:
!     DTAUX(1:L1_,1:W_) = optical depth of each layer
!     POMEGAX(1:8,1:L1_,1:W_) = scattering phase fn (multiplied by s-s abledo)
!     U0  = cos (SZA)
!     RFL(1:W_) = Lambertian albedo of surface
!     AMF2(1:2*L1_+1,1:2*L1_+1) = air mass factor (I,L)=wt of layer-I to layer-L
!        AMF2 now does both edges and middle of CTM layers
!     JXTRA(1:L1_) = number 0:J = no. of additional levels to be inserted
! out:
!     FJACT(1:L_,1:W_) = mean actinic flux(diff+direct) at std CTM levels(mid-lyr)
!  (new ver 5.7 diagnostics for fluxes, deposition)  fluxes 'down' are <0
!     FJTOP(1:W_) = diffuse flux out top-of-atmosphere (TAU=0 above top model lyr)
!     FJBOT(1:W_) = diffuse flux onto surface (<0 by definition)
!     FSBOT(1:W_) = direct/solar flux onto surface  (<0 by definition)
!     FJFLX(1:L_,1:W_) = diffuse flux across top of model layer L
!        this connects with FJBOT = FJFLX(0) & FJTOP = FJFLX(L_+1) (not dim!!)
!     FLXD(1:L_+1,1:W_) = solar flux deposited in layer L (includes lyr above CTM)
!        this should take into account sphericity, and is not just = mu0
!     FLXD0(1:W_) = sum of solar flux deposited in atmos
!        does NOT include flux on lower surface, does NOT mean absorbed!
!-----------------------------------------------------------------------
!
!     DTAU     Local optical depth of each CTM level
!     TTAU     Optical depth of air vertically above each point (to top of atm)
!     FTAU2     Attenuation of solar beam
!     POMEGAJ  Scattering phase function
!
!---new ver 5.3 code adds sub-layers (# = JXTRA(L2)) using ATAU as the
!   factor increase from sub-layer to sub-layer
!
!---------------------SET UP FOR MIE CODE-------------------------------
!
!-----------------wavelength independent--------------------------------
!
!  Transpose the ascending TTAU grid to a descending ZTAU grid.
!  Double the resolution - TTAU points become the odd points on the
!  ZTAU grid, even points needed for asymm phase fn soln, contain 'h'.
!  Odd point added at top of grid for unattenuated beam   (Z='inf')
!
!  The following mapping holds for JADDLV=0
!        Surface:   TTAU(1)    ==> ZTAU(2*L2_+1)
!        Top:       TTAU(L2_)  ==> ZTAU(3)
!        Infinity:     0.0     ==> ZTAU(1)
!        index: 2*(L2_+1-L2)+1 ==> LZ
!
!  Mie scattering code only used from surface to level L2_
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
!  Insert new levels, working downwards from the top of the atmosphere
!  to the surface (down in 'LZ', up in 'L2'). This allows ztau and pomega
!  to be incremented linearly, and the flux fz to be attenuated top-down
!    (avoiding problems where lower level fluxes are zero).
!------------------------------------------------------------------------
!
!  Ascend through atmosphere transposing grid and adding extra points
!  remember L2=1 is surface of CTM, but last layer (LZ) in scattering code.
!  there are twice the number of layers in the LZ arrays (2*L2_ + 2*JADDTO + 1)
!    because we need to insert the intermediate layers (even LZ) for the
!    asymmetric scattering code.
!
!  Transfer the L2=1:L2_+1 values (TTAU,FTAU2,POMEGAJ) onto the reverse
!    order, expanded, doubled-level scatter grid.
!    Note that we need to deal with the expansion by JADD levels (L2L).
!      These JADDLV levels are skipped and need to be interpolated later.
!    Note that only odd LZ levels are filled,
!
!----------------------re-grid data---------------------------------------------
!  Calculate cumulative total and define levels we want J-values at.
!  Sum upwards for levels, and then downwards for Mie code readjustments.
!
!     JXTRA(L2)  Number of new levels to add between (L2) and (L2+1)
!           ***JXTRA(1:L2_+1) is calculated based on the aerosol+cloud OD_s
!     JADDLV(L2)  Number of new levels actually added at each wavelength
!            where JADDLV = 0 when there is effectively no FTAU2
!     JADDTO(L2)   Total number of new levels to add to and above level (L2)
!     JNDLEV(L) = L2 index that maps on CTM mid-layer L
!
!---JADDLV(L2=1:L2_) = number of levels to add between TTAU2(L2) and TTAU(L2+1)
!---    JADDLV is taken from JXTRA, which is based on visible OD.
!---    JADDTO(L2=1:L2_+1) is the cumulative number of levels to be added
!---these should be fixed for all wavelengths to lock-in the array sizes
      if (LU .gt. JXL_) then
         write(*,*) ' OPMIE:  JXL_ .lt. L_', LU, JXL_
         stop
      endif

      L1U = LU + 1
      L2U = 2*LU + 2


      do L2 = 1,L2U,1
        JADDLV(L2) = JXTRA(L2)             ! here surf=1
      enddo
        JADDTO(L2U+1) = 0
      do L2 = L2U,1,-1
        JADDTO(L2) = JADDTO(L2+1) + JADDLV(L2)  ! here surf = 1
      enddo

!---expanded grid now included CTM edge and mid layers plus expanded
!---    grid to allow for finer delta-tau at tops of clouds.
!---    DIM of new grid = L2U + JADDTO(1) + 1

!---L2LEV(L2) = L2-index for old level L2 in expanded J-grid (w/JADDLV)
!     in absence of JADDLV, L2LEV(L2) = L2
        L2LEV(1)  = 1
      do L2 = 2,L2U+1
        L2LEV(L2) = L2LEV(L2-1) + 1 + JADDLV(L2-1)    ! here surf = 1
      enddo

!---JNDLEV(L=1:L_) = L2-index in expanded grid for CTM mid-layer L
!---JNELEV(L=1:L_) = L2-index for top of layer L
      do L = 1,LU
        JNDLEV(L) = L2LEV(2*L)                        ! here surf=1
        JNELEV(L) = L2LEV(2*L+1)
      enddo
        JNELEV(LU+1) = 0  !need to set this to top-of-atmosphere

      ND = 2*L2U + 2*JADDTO(1) + 1

      if(ND .gt. N_) then
        write(*,*) ' overflow of scatter arrays: ND > N_'
        stop
      endif

!----------------begin wavelength dependent set up------------------------------

!---Reinitialize arrays
      ZTAU(:,:)     = 0.d0
      FZ(:,:)       = 0.d0
      POMEGA(:,:,:) = 0.d0

      FJACT(:,:) = 0.d0
      FJTOP(:) = 0.d0
      FJBOT(:) = 0.d0
      FSBOT(:) = 0.d0
      FJFLX(:,:) = 0.d0
      FLXD(:,:) = 0.d0
      FLXD0(:) = 0.d0
      FJT(:) = 0.d0
      FJB(:) = 0.d0
      FJ(:,:) = 0.d0
      FZ(:,:) = 0.d0
      ZTAU(:,:) = 0.d0

      do K=1,W_
!      if (FL(K) .gt. 1.d0) then

!---Set up optical depth DTAU(L)
       do L = 1,L1U
        DTAU(L,K) = DTAUX(L,K)             ! here surf = 1
       enddo
        DTAU(L1U+1,K) = 0.d0

!---Define the total scattering phase fn for each CTM layer L=1:L_+1
!---   from a DTAU-wt_d mix of aerosols, cloud & Rayleigh
!---No. of quadrature pts fixed at 4(M_), expansion of phase fn @ 8
       do L = 1,L1U
        do I = 1,M2_
          POMEGAJ(I,L,K) = POMEGAX(I,L,K)       ! here surf = 1
        enddo
       enddo

!---Calculate attenuated incident beam exp(-TTAU/U0 = DTAU * AirMassFactor)
!---      at the middle & edges of the CTM layers L=1:2*L1_+1
!---  L1_ is top-edge of CTM (ie, L=38 = 2 hPa) which has TAU > 0
!---  note that DTAU(L1_) is optical depth in the FULL CTM layer just above
        FTAU2(:,:) = 0.d0                       ! here surf = 1
        FTAU2(L2U+1,:) = 1.0d0
       do LL = 1,2*L1U+1
         L = (LL+1)/2
        if (AMF2(LL,LL) .gt. 0.0d0) then
           XLTAU = 0.0d0
         do II = 1,2*L1U+1
           I = (II+1)/2
           XLTAU = XLTAU + 0.5d0*DTAU(I,K)*AMF2(II,LL)
         enddo
         if (XLTAU .lt. 76.d0) then   ! zero out flux at 1e-33
          FTAU2(LL,K) = exp(-XLTAU)
         endif
        endif
       enddo

!---calculate direct solar flux deposited in each CTM half-layer: L=1:L2_
!---     use FSBOT for surface flux, cannot do layer above CTM (L_+1)
          FLXD2(:,:) = 0.d0                    ! here surf = 1
       do LL = 1,2*L1U
        if (AMF2(LL,LL) .gt. 0.d0) then
          FLXD2(LL,K) = (FTAU2(LL+1,K) - FTAU2(LL,K))/AMF2(LL,LL)
        endif
       enddo
        if (AMF2(1,1) .gt. 0.d0) then           ! bottom = surf = 1
          FSBOT(K) = FTAU2(1,K)/AMF2(1,1)
        else
          FSBOT(K) = 0.d0
        endif

       do LL = 2,2*L1U,2
         L=LL/2
         FLXD(L,K) = FLXD2(LL,K)+FLXD2(LL-1,K)   ! here surf = 1
       enddo

!---integrate solar flux depositied in CTM layers L=1:L_, cannot do top layer
!---  note FLXD0 .ne. (1.d0 - FTAU(L_+1))/AMF(L_+1,L_+1) with spherical atmos
        FLXD0(K) = 0.d0
       if (AMF2(2*L1U,2*L1U) .gt. 0.d0) then
        do L=1,L1U
         FLXD0(K) = FLXD0(K) + FLXD(L,K)
        enddo
       endif

!------------------------------------------------------------------------
!  Take optical properties on CTM layers and convert to a photolysis
!  level grid corresponding to layer centres and boundaries. This is
!  required so that J-values can be calculated for the centre of CTM
!  layers; the index of these layers is kept in the JNDLEV array.
!------------------------------------------------------------------------
!---Now combine the CTM layer edges (1:L_+2) with the CTM mid-layer
!---    points (1:L_) plus 1 for the mid point of added top layer.
!---combine these edge- and mid-layer points into grid of size:
!---              L2_+1 = 2*L1_+1 = 2*L_+3
!---calculate column optical depths above each level, TTAU(1:L2_+1)
!---      note that TTAU(L2_+1)=0 and TTAU(1)=total OD

        TTAU(L2U+1,K) = 0.0d0
       do L2 = L2U,1,-1
        L          = (L2+1)/2
        DTAUJ      = 0.5d0 * DTAU(L,K)
        TTAU(L2,K)   = TTAU(L2+1,K) + DTAUJ
       enddo

!----solar flux incident on lower boundary & Lambertian reflect factor:
       if (FSBOT(K) .gt. 0.d0) then
        ZFLUX(K) = FSBOT(K)*RFL(K)/(1.d0+RFL(K))
       else
        ZFLUX(K) = 0.d0
       endif

!  Calculate scattering properties, level centres then level boundaries
!>>>>>be careful of order, we are overwriting/shifting the 'POMEGAJ' upward in index
       do L2 = L2U,2,-2
        L   = L2/2
        do I = 1,M2_
          POMEGAJ(I,L2,K) = POMEGAJ(I,L,K)
        enddo
       enddo
!---lower boundary value is set (POMEGAJ(I,1)), but set upper:
       do I = 1,M2_
         POMEGAJ(I,L2U+1,K) = POMEGAJ(I,L2U,K)
       enddo
!---now have POMEGAJ filled at even points from L2=3:L2_-1
!---use inverse interpolation for correct tau-weighted values at edges
       do L2 = 3,L2U-1,2
        TAUDN = TTAU(L2-1,K)-TTAU(L2,K)
        TAUUP = TTAU(L2,K)-TTAU(L2+1,K)
        do I = 1,M2_
          POMEGAJ(I,L2,K) = (POMEGAJ(I,L2-1,K)*TAUDN + &
                 POMEGAJ(I,L2+1,K)*TAUUP) / (TAUDN+TAUUP)
        enddo
       enddo

!---at this point FTAU2(1:L2_+1) and POMEAGJ(1:8, 1:L2_+1)
!---    where FTAU2(L2_+1) = 1.0 = top-of-atmos, FTAU2(1) = surface
       do L2 = 1,L2U+1          ! L2 = index of CTM edge- and mid-layers
        L2L = L2LEV(L2)        ! L2L = index for L2 in expanded scale(JADD)
        LZ  = ND + 2 - 2*L2L  ! LZ = index for L2 in scatt arrays
          ZTAU(LZ,K) = TTAU(L2,K)
          FZ(LZ,K)   = FTAU2(L2,K)
        do I=1,M2_
          POMEGA(I,LZ,K) = POMEGAJ(I,L2,K)
        enddo
       enddo

!   Now go thru the pairs of L2 levels to see if we need JADD levels
       do L2 = 1,L2U             ! L2 = index of CTM edge- and mid-layers
         L2L = L2LEV(L2)         ! L2L = index for L2 in expanded scale(JADD)
         LZ  = ND + 2 - 2*L2L   ! LZ = index for L2 in scatt arrays
         L22 = L2LEV(L2+1) - L2LEV(L2) - 1   ! L22 = 0 if no added levels

        if (L22 .gt. 0) then
          TAUBTM(K) = TTAU(L2,K)
          TAUTOP(K) = TTAU(L2+1,K)
          FBTM(K)   = FTAU2(L2,K)
          FTOP(K)   = FTAU2(L2+1,K)
         do I = 1,M2_
          POMEGAB(I,K) = POMEGAJ(I,L2,K)
         enddo

!---to fit L22 new layers between TAUBOT > TAUTOP, calculate new 1/ATAU factor
!---  such that TAU(just above TAU-btm) = ATUAZ * TAUBTM < TAUBTM
         ATAUZ = exp(-log(TAUBTM(K)/max(TAUTOP(K),ATAU0))/float(L22+1))
         do L = 1,L22           ! add odd levels between L2LEV(L2) & L2LEV(L2+1)
          LZZ = LZ - 2*L       ! LZZ = index(odd) of added level in scatt arrays
          ZTAU(LZZ,K) = TAUBTM(K) * ATAUZ

!---fraction from TAUBTM=>TAUTOP
          ATAUA=(TAUBTM(K)-ZTAU(LZZ,K))/(TAUBTM(K)-TAUTOP(K))
!---solar flux at interp-levels: use exp(TAU/U0) if U0>0.02 (89 deg),
!---else scale by TAU
          if (U0 .gt. 0.02d0) then
            FZ(LZZ,K) = FTOP(K) * exp((TAUTOP(K)-ZTAU(LZZ,K))/U0)
          else
            if (FBTM(K) .lt. 1.d-32) then
              FZ(LZZ,K) = 0.d0
            else
              FZ(LZZ,K) = FBTM(K) * (FTOP(K)/FBTM(K))**ATAUA
            endif
          endif
          do I = 1,M2_
            POMEGA(I,LZZ,K) = POMEGAB(I,K) + &
                     ATAUA*(POMEGAJ(I,L2+1,K)-POMEGAB(I,K))
          enddo
            TAUBTM(K)    = ZTAU(LZZ,K)
            FBTM(K)      = FZ(LZZ,K)
          do I = 1,M2_
            POMEGAB(I,K) = POMEGA(I,LZZ,K)
          enddo
         enddo
        endif
       enddo

!   Now fill in the even points with simple interpolation in scatter arrays:
       do LZ = 2,ND-1,2
         ZTAU(LZ,K) = 0.5d0*(ZTAU(LZ-1,K)+ZTAU(LZ+1,K))
         FZ(LZ,K)   = sqrt(FZ(LZ-1,K)*FZ(LZ+1,K))
        do I=1,M2_
         POMEGA(I,LZ,K) = 0.5d0*(POMEGA(I,LZ-1,K)+POMEGA(I,LZ+1,K))
        enddo
       enddo

!      endif
      enddo  ! wavelength loop!

!-----------------------------------------------------------------------
       call MIESCT(FJ,FJT,FJB,POMEGA,FZ,ZTAU,ZFLUX,RFL,U0,ND)
!-----------------------------------------------------------------------

!---Move mean intensity from scatter array FJ(LZ=1:ND)
!---              to CTM mid-level array FJACT(L=1:L_)

      do K=1,W_
!      if (FL(K) .gt. 1.0d0) then

!---mean intensity at mid-layer:  4*<I> + solar
!       do L = 1,LU
!        L2L = JNDLEV(L)
!        LZ  = ND+2 - 2*L2L
!        FJACT(L,K) = 4.d0*FJ(LZ,K) + FZ(LZ,K)
!       enddo

!---mean intensity averaged throughout layer:
       do L = 1,LU
         LZ0 = ND+2 - 2*JNELEV(L)
        if (L .gt. 1) then
         LZ1 = ND+2 - 2*JNELEV(L-1)
        else
         LZ1 = ND
        endif
         SUMJ = (4.d0*FJ(LZ0,K)+FZ(LZ0,K))*(ZTAU(LZ0+2,K)-ZTAU(LZ0,K)) &
               +(4.d0*FJ(LZ1,K)+FZ(LZ1,K))*(ZTAU(LZ1,K)-ZTAU(LZ1-2,K))
         SUMT = ZTAU(LZ0+2,K)-ZTAU(LZ0,K) + ZTAU(LZ1,K)-ZTAU(LZ1-2,K)

        do LZ = LZ0+2,LZ1-2,2
         SUMJ =SUMJ+(4.d0*FJ(LZ,K)+FZ(LZ,K))*(ZTAU(LZ+2,K)-ZTAU(LZ-2,K))
         SUMT =SUMT + ZTAU(LZ+2,K)-ZTAU(LZ-2,K)
        enddo
        FJACT(L,K) = SUMJ/SUMT
       enddo

!---mean diffuse flux:  4<I*mu> (not solar) at top of layer L
!---      average (tau-wtd) the h's just above and below the L-edge
       do L = 1,LU
        L2L = JNELEV(L)
        LZ  = ND+2 - 2*L2L
        FJFLX0 = (ZTAU(LZ+1,K)-ZTAU(LZ,K))/(ZTAU(LZ+1,K)-ZTAU(LZ-1,K))
        FJFLX(L,K)=4.d0*(FJ(LZ-1,K)*FJFLX0 + FJ(LZ+1,K)*(1.d0-FJFLX0))
       enddo

!---diffuse fluxes reflected at top, incident at bottom
         FJTOP(K) = FJT(K)
         FJBOT(K) = FJB(K)
!      endif
      enddo  ! wavelength loop!
 
      return
      END SUBROUTINE OPMIE

!-----------------------------------------------------------------------
      subroutine MIESCT(FJ,FJT,FJB, POMEGA,FZ,ZTAU,ZFLUX,RFL,U0,ND)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) ::  ND
      real*8, intent(in)  ::  POMEGA(M2_,N_,W_),FZ(N_,W_),ZTAU(N_,W_) &
                             ,RFL(W_),U0,ZFLUX(W_)
      real*8, intent(out) ::  FJ(N_,W_),FJT(W_),FJB(W_)

      real*8  PM(M_,M2_),PM0(M2_)
      integer I, IM  ,K

!-----------------------------------------------------------------------
!   This is an adaption of the Prather radiative transfer code, (mjp, 10/95)
!     Prather, 1974, Astrophys. J. 192, 787-792.
!         Solution of inhomogeneous Rayleigh scattering atmosphere.
!         (original Rayleigh w/ polarization)
!     Cochran and Trafton, 1978, Ap.J., 219, 756-762.
!         Raman scattering in the atmospheres of the major planets.
!         (first use of anisotropic code)
!     Jacob, Gottlieb and Prather, 1989, J.Geophys.Res., 94, 12975-13002.
!         Chemistry of a polluted cloudy boundary layer,
!         (documentation of extension to anisotropic scattering)
!
!    takes atmospheric structure and source terms from std J-code
!    ALSO limited to 4 Gauss points, only calculates mean field! (M=1)
!-----------------------------------------------------------------------
      do I = 1,M_
       call LEGND0 (EMU(I),PM0,M2_)
       do IM = 1,M2_
         PM(I,IM) = PM0(IM)
       enddo
      enddo

       call LEGND0 (-U0,PM0,M2_)
       do IM=1,M2_
         PM0(IM) = 0.25d0*PM0(IM)
       enddo

!---BLKSLV now called with all the wavelength arrays (K=1:W_)

      call BLKSLV(FJ,POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,FJT,FJB, ND)

      return
      END SUBROUTINE MIESCT
 
!-----------------------------------------------------------------------
      subroutine LEGND0 (X,PL,N)
!-----------------------------------------------------------------------
!---Calculates ORDINARY Legendre fns of X (real)
!---   from P[0] = PL(1) = 1,  P[1] = X, .... P[N-1] = PL(N)
      implicit none
      integer, intent(in) :: N
      real*8, intent(in)  :: X
      real*8, intent(out) :: PL(N)
      integer I
      real*8  DEN
!---Always does PL(2) = P[1]
        PL(1) = 1.d0
        PL(2) = X
        do I = 3,N
         DEN = (I-1)
         PL(I) = PL(I-1)*X*(2.d0-1.0/DEN) - PL(I-2)*(1.d0-1.d0/DEN)
        enddo
      return
      END SUBROUTINE LEGND0

!-----------------------------------------------------------------------
      subroutine BLKSLV &
           (FJ,POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,FJTOP,FJBOT,ND)
!-----------------------------------------------------------------------
!  Sets up and solves the block tri-diagonal system:
!               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
!  This goes back to the old, dumb, fast version 5.3
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) ::  ND
      real*8, intent(in)  ::  POMEGA(M2_,N_,W_),FZ(N_,W_),ZTAU(N_,W_) &
                             ,PM(M_,M2_),PM0(M2_) &
                             ,RFL(W_),ZFLUX(W_)
      real*8, intent(out) ::  FJ(N_,W_),FJTOP(W_),FJBOT(W_)

      real*8, dimension(M_,N_,W_)    ::  A,C,H,   RR

      real*8, dimension(M_,M_,N_,W_) ::  B,AA,CC,  DD
      real*8, dimension(M_,M_) ::  E
      real*8  SUMB,SUMBX,SUMT
      integer I, J, K, L

      do K = 1,W_
!      if (FL(K) .gt. 1.0d0) then
!       call GEN_ID (POMEGA(1,1,K),FZ(1,K),ZTAU(1,K),ZFLUX(K),RFL(K), &
!           PM,PM0, B(1,1,1,K),CC(1,1,1,K),AA(1,1,1,K), &
!                   A(1,1,K),H(1,1,K),C(1,1,K), ND)
       call GEN_ID (POMEGA(:,:,K),FZ(:,K),ZTAU(:,K),ZFLUX(K),RFL(K), &
           PM,PM0, B(:,:,:,K),CC(:,:,:,K),AA(:,:,:,K), &
                   A(:,:,K),H(:,:,K),C(:,:,K), ND)

!      endif
      enddo

      do K = 1,W_
!      if (FL(K) .gt. 1.0d0) then
!-----------UPPER BOUNDARY L=1
       L = 1
        do J = 1,M_
         do I = 1,M_
          E(I,J) = B(I,J,1,K)
         enddo
        enddo

!---setup L & U matrices
         E(2,1) = E(2,1)/E(1,1)
         E(2,2) = E(2,2)-E(2,1)*E(1,2)
         E(2,3) = E(2,3)-E(2,1)*E(1,3)
         E(2,4) = E(2,4)-E(2,1)*E(1,4)
         E(3,1) = E(3,1)/E(1,1)
         E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
         E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
         E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
         E(4,1) = E(4,1)/E(1,1)
         E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
         E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
         E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
!---invert L
         E(4,3) = -E(4,3)
         E(4,2) = -E(4,2)-E(4,3)*E(3,2)
         E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
         E(3,2) = -E(3,2)
         E(3,1) = -E(3,1)-E(3,2)*E(2,1)
         E(2,1) = -E(2,1)
!---invert U
         E(4,4) = 1.d0/E(4,4)
         E(3,4) = -E(3,4)*E(4,4)/E(3,3)
         E(3,3) = 1.d0/E(3,3)
         E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
         E(2,3) = -E(2,3)*E(3,3)/E(2,2)
         E(2,2) = 1.d0/E(2,2)
         E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
         E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
         E(1,2) = -E(1,2)*E(2,2)/E(1,1)
         E(1,1) = 1.d0/E(1,1)
!---multiply U-invers * L-inverse
         E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
         E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
         E(1,3) = E(1,3)+E(1,4)*E(4,3)
         E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
         E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
         E(2,3) = E(2,3)+E(2,4)*E(4,3)
         E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
         E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
         E(3,3) = E(3,3)+E(3,4)*E(4,3)
         E(4,1) = E(4,4)*E(4,1)
         E(4,2) = E(4,4)*E(4,2)
         E(4,3) = E(4,4)*E(4,3)

        do J = 1,M_
         do I = 1,M_
          DD(I,J,1,K) = -E(I,1)*CC(1,J,1,K)-E(I,2)*CC(2,J,1,K) &
                        -E(I,3)*CC(3,J,1,K)-E(I,4)*CC(4,J,1,K)
         enddo
          RR(J,1,K) = E(J,1)*H(1,1,K)+E(J,2)*H(2,1,K) &
                    +E(J,3)*H(3,1,K)+E(J,4)*H(4,1,K)
        enddo

!----------CONTINUE THROUGH ALL DEPTH POINTS ID=2 TO ID=ND-1
       do L = 2,ND-1

        do J = 1,M_
         do I = 1,M_
          B(I,J,L,K) = B(I,J,L,K) + A(I,L,K)*DD(I,J,L-1,K)
         enddo
          H(J,L,K) = H(J,L,K) - A(J,L,K)*RR(J,L-1,K)
        enddo

        do J = 1,M_
         do I = 1,M_
          E(I,J) = B(I,J,L,K)
         enddo
        enddo

!---setup L & U matrices
         E(2,1) = E(2,1)/E(1,1)
         E(2,2) = E(2,2)-E(2,1)*E(1,2)
         E(2,3) = E(2,3)-E(2,1)*E(1,3)
         E(2,4) = E(2,4)-E(2,1)*E(1,4)
         E(3,1) = E(3,1)/E(1,1)
         E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
         E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
         E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
         E(4,1) = E(4,1)/E(1,1)
         E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
         E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
         E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
!---invert L
         E(4,3) = -E(4,3)
         E(4,2) = -E(4,2)-E(4,3)*E(3,2)
         E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
         E(3,2) = -E(3,2)
         E(3,1) = -E(3,1)-E(3,2)*E(2,1)
         E(2,1) = -E(2,1)
!---invert U
         E(4,4) = 1.d0/E(4,4)
         E(3,4) = -E(3,4)*E(4,4)/E(3,3)
         E(3,3) = 1.d0/E(3,3)
         E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
         E(2,3) = -E(2,3)*E(3,3)/E(2,2)
         E(2,2) = 1.d0/E(2,2)
         E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
         E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
         E(1,2) = -E(1,2)*E(2,2)/E(1,1)
         E(1,1) = 1.d0/E(1,1)
!---multiply U-invers * L-inverse
         E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
         E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
         E(1,3) = E(1,3)+E(1,4)*E(4,3)
         E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
         E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
         E(2,3) = E(2,3)+E(2,4)*E(4,3)
         E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
         E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
         E(3,3) = E(3,3)+E(3,4)*E(4,3)
         E(4,1) = E(4,4)*E(4,1)
         E(4,2) = E(4,4)*E(4,2)
         E(4,3) = E(4,4)*E(4,3)

        do J = 1,M_
         do I = 1,M_
          DD(I,J,L,K) = - E(I,J)*C(J,L,K)
         enddo
          RR(J,L,K) = E(J,1)*H(1,L,K)+E(J,2)*H(2,L,K) &
                    + E(J,3)*H(3,L,K)+E(J,4)*H(4,L,K)
        enddo

       enddo

!---------FINAL DEPTH POINT: L=ND
       L = ND
        do J = 1,M_
         do I = 1,M_
          B(I,J,L,K) = B(I,J,L,K) &
           + AA(I,1,L,K)*DD(1,J,L-1,K) + AA(I,2,L,K)*DD(2,J,L-1,K) &
           + AA(I,3,L,K)*DD(3,J,L-1,K) + AA(I,4,L,K)*DD(4,J,L-1,K)
         enddo
          H(J,L,K) = H(J,L,K) &
           - AA(J,1,L,K)*RR(1,L-1,K) - AA(J,2,L,K)*RR(2,L-1,K) &
           - AA(J,3,L,K)*RR(3,L-1,K) - AA(J,4,L,K)*RR(4,L-1,K)
        enddo

        do J = 1,M_
         do I = 1,M_
          E(I,J) = B(I,J,L,K)
         enddo
        enddo

!---setup L & U matrices
         E(2,1) = E(2,1)/E(1,1)
         E(2,2) = E(2,2)-E(2,1)*E(1,2)
         E(2,3) = E(2,3)-E(2,1)*E(1,3)
         E(2,4) = E(2,4)-E(2,1)*E(1,4)
         E(3,1) = E(3,1)/E(1,1)
         E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
         E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
         E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
         E(4,1) = E(4,1)/E(1,1)
         E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
         E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
         E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
!---invert L
         E(4,3) = -E(4,3)
         E(4,2) = -E(4,2)-E(4,3)*E(3,2)
         E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
         E(3,2) = -E(3,2)
         E(3,1) = -E(3,1)-E(3,2)*E(2,1)
         E(2,1) = -E(2,1)
!---invert U
         E(4,4) = 1.d0/E(4,4)
         E(3,4) = -E(3,4)*E(4,4)/E(3,3)
         E(3,3) = 1.d0/E(3,3)
         E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
         E(2,3) = -E(2,3)*E(3,3)/E(2,2)
         E(2,2) = 1.d0/E(2,2)
         E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
         E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
         E(1,2) = -E(1,2)*E(2,2)/E(1,1)
         E(1,1) = 1.d0/E(1,1)
!---multiply U-invers * L-inverse
         E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
         E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
         E(1,3) = E(1,3)+E(1,4)*E(4,3)
         E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
         E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
         E(2,3) = E(2,3)+E(2,4)*E(4,3)
         E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
         E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
         E(3,3) = E(3,3)+E(3,4)*E(4,3)
         E(4,1) = E(4,4)*E(4,1)
         E(4,2) = E(4,4)*E(4,2)
         E(4,3) = E(4,4)*E(4,3)

        do J = 1,M_
         RR(J,L,K) = E(J,1)*H(1,L,K)+E(J,2)*H(2,L,K) &
                    +E(J,3)*H(3,L,K)+E(J,4)*H(4,L,K)
        enddo

!-----------BACK SOLUTION
       do L = ND-1,1,-1
        do J = 1,M_
         RR(J,L,K) = RR(J,L,K) &
          + DD(J,1,L,K)*RR(1,L+1,K) + DD(J,2,L,K)*RR(2,L+1,K) &
          + DD(J,3,L,K)*RR(3,L+1,K) + DD(J,4,L,K)*RR(4,L+1,K)
        enddo
       enddo

!----------mean J & H
       do L = 1,ND,2
        FJ(L,K) = RR(1,L,K)*WT(1) + RR(2,L,K)*WT(2) &
                + RR(3,L,K)*WT(3) + RR(4,L,K)*WT(4)
       enddo
       do L = 2,ND,2
        FJ(L,K) = RR(1,L,K)*WT(1)*EMU(1) + RR(2,L,K)*WT(2)*EMU(2) &
                + RR(3,L,K)*WT(3)*EMU(3) + RR(4,L,K)*WT(4)*EMU(4)
       enddo

!---FJTOP = scaled diffuse flux out top-of-atmosphere (limit = mu0)
!---FJBOT = scaled diffuse flux onto surface:
!---ZFLUX = reflect/(1 + reflect) * mu0 * Fsolar(lower boundary)
!---SUMBX = flux from Lambert reflected I+
       SUMT = RR(1, 1,K)*WT(1)*EMU(1) + RR(2, 1,K)*WT(2)*EMU(2) &
            + RR(3, 1,K)*WT(3)*EMU(3) + RR(4, 1,K)*WT(4)*EMU(4)
       SUMB = RR(1,ND,K)*WT(1)*EMU(1) + RR(2,ND,K)*WT(2)*EMU(2) &
            + RR(3,ND,K)*WT(3)*EMU(3) + RR(4,ND,K)*WT(4)*EMU(4)
       SUMBX = 4.d0*SUMB*RFL(K)/(1.0d0 + RFL(K)) + ZFLUX(K)

       FJTOP(K) = 4.d0*SUMT
       FJBOT(K) = 4.d0*SUMB - SUMBX

!      endif
      enddo
      return
      END SUBROUTINE BLKSLV

!-----------------------------------------------------------------------
      subroutine GEN_ID(POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0 &
                    ,B,CC,AA,A,H,C,  ND)
!-----------------------------------------------------------------------
!  Generates coefficient matrices for the block tri-diagonal system:
!               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) ::  ND
      real*8, intent(in)  ::  POMEGA(M2_,N_),PM(M_,M2_),PM0(M2_)
      real*8, intent(in)  ::  ZFLUX,RFL
      real*8, intent(in),dimension(N_) :: FZ,ZTAU

      real*8, intent(out),dimension(M_,M_,N_) ::  B,AA,CC
      real*8, intent(out),dimension(M_,N_) ::  A,C,H

      integer I, J, K, L1,L2,LL
      real*8  SUM0, SUM1, SUM2, SUM3
      real*8  DELTAU, D1, D2, SURFAC
!
      real*8, dimension(M_,M_) :: S,T,U,V,W
!---------------------------------------------

!---------upper boundary:  2nd-order terms
       L1 = 1
       L2 = 2
!       print*, 'TEST GEN_ID: ', POMEGA(1,L1), POMEGA(5,L1)
       do I = 1,M_
        SUM0 = &
         POMEGA(1,L1)*PM(I,1)*PM0(1) + POMEGA(3,L1)*PM(I,3)*PM0(3) &
       + POMEGA(5,L1)*PM(I,5)*PM0(5) + POMEGA(7,L1)*PM(I,7)*PM0(7)
        SUM2 = &
         POMEGA(1,L2)*PM(I,1)*PM0(1) + POMEGA(3,L2)*PM(I,3)*PM0(3) &
       + POMEGA(5,L2)*PM(I,5)*PM0(5) + POMEGA(7,L2)*PM(I,7)*PM0(7)
        SUM1 = &
         POMEGA(2,L1)*PM(I,2)*PM0(2) + POMEGA(4,L1)*PM(I,4)*PM0(4) &
       + POMEGA(6,L1)*PM(I,6)*PM0(6) + POMEGA(8,L1)*PM(I,8)*PM0(8)
        SUM3 = &
         POMEGA(2,L2)*PM(I,2)*PM0(2) + POMEGA(4,L2)*PM(I,4)*PM0(4) &
       + POMEGA(6,L2)*PM(I,6)*PM0(6) + POMEGA(8,L2)*PM(I,8)*PM0(8)
         H(I,L1) = 0.5d0*(SUM0*FZ(L1) + SUM2*FZ(L2))
         A(I,L1) = 0.5d0*(SUM1*FZ(L1) + SUM3*FZ(L2))
       enddo

       do I = 1,M_
        do J = 1,I
         SUM0 = &
         POMEGA(1,L1)*PM(I,1)*PM(J,1) + POMEGA(3,L1)*PM(I,3)*PM(J,3) &
       + POMEGA(5,L1)*PM(I,5)*PM(J,5) + POMEGA(7,L1)*PM(I,7)*PM(J,7)
         SUM2 = &
         POMEGA(1,L2)*PM(I,1)*PM(J,1) + POMEGA(3,L2)*PM(I,3)*PM(J,3) &
       + POMEGA(5,L2)*PM(I,5)*PM(J,5) + POMEGA(7,L2)*PM(I,7)*PM(J,7)
         SUM1 = &
         POMEGA(2,L1)*PM(I,2)*PM(J,2) + POMEGA(4,L1)*PM(I,4)*PM(J,4) &
       + POMEGA(6,L1)*PM(I,6)*PM(J,6) + POMEGA(8,L1)*PM(I,8)*PM(J,8)
         SUM3 = &
         POMEGA(2,L2)*PM(I,2)*PM(J,2) + POMEGA(4,L2)*PM(I,4)*PM(J,4) &
       + POMEGA(6,L2)*PM(I,6)*PM(J,6) + POMEGA(8,L2)*PM(I,8)*PM(J,8)
         S(I,J) = - SUM2*WT(J)
         S(J,I) = - SUM2*WT(I)
         T(I,J) = - SUM1*WT(J)
         T(J,I) = - SUM1*WT(I)
         V(I,J) = - SUM3*WT(J)
         V(J,I) = - SUM3*WT(I)
         B(I,J,L1) = - 0.5d0*(SUM0 + SUM2)*WT(J)
         B(J,I,L1) = - 0.5d0*(SUM0 + SUM2)*WT(I)
        enddo
       enddo

       do I = 1,M_
         S(I,I)   = S(I,I)   + 1.0d0
         T(I,I)   = T(I,I)   + 1.0d0
         V(I,I)   = V(I,I)   + 1.0d0
         B(I,I,L1)= B(I,I,L1) + 1.0d0

         C(I,L1)= S(I,1)*A(1,L1)/EMU(1) + S(I,2)*A(2,L1)/EMU(2) &
                + S(I,3)*A(3,L1)/EMU(3) + S(I,4)*A(4,L1)/EMU(4)
       enddo

       do I = 1,M_
        do J = 1,M_
         W(J,I) = S(J,1)*T(1,I)/EMU(1) + S(J,2)*T(2,I)/EMU(2) &
                + S(J,3)*T(3,I)/EMU(3) + S(J,4)*T(4,I)/EMU(4)
         U(J,I) = S(J,1)*V(1,I)/EMU(1) + S(J,2)*V(2,I)/EMU(2) &
                + S(J,3)*V(3,I)/EMU(3) + S(J,4)*V(4,I)/EMU(4)
        enddo
       enddo
!-------------upper boundary, 2nd-order, C-matrix is full (CC)
         DELTAU = ZTAU(L2) - ZTAU(L1)
         D2 = 0.25d0*DELTAU
       do I = 1,M_
        do J = 1,M_
         B(I,J,L1) = B(I,J,L1) + D2*W(I,J)
         CC(I,J,L1) = D2*U(I,J)
        enddo
         H(I,L1) = H(I,L1) + 2.0d0*D2*C(I,L1)
         A(I,L1) = 0.0d0
       enddo
       do I = 1,M_
        D1 = EMU(I)/DELTAU
        B(I,I,L1)  = B(I,I,L1) + D1
        CC(I,I,L1) = CC(I,I,L1) - D1
       enddo

!------------intermediate points:  can be even or odd, A & C diagonal
!---mid-layer h-points, Legendre terms 2,4,6,8
       do LL=2,ND-1,2
        DELTAU = ZTAU(LL+1) - ZTAU(LL-1)
        do I = 1,M_
          A(I,LL) = EMU(I)/DELTAU
          C(I,LL) = -A(I,LL)
          H(I,LL) = FZ(LL)*( &
           POMEGA(2,LL)*PM(I,2)*PM0(2) + POMEGA(4,LL)*PM(I,4)*PM0(4) &
         + POMEGA(6,LL)*PM(I,6)*PM0(6) + POMEGA(8,LL)*PM(I,8)*PM0(8))
        enddo
        do I = 1,M_
         do J=1,I
          SUM0 = &
           POMEGA(2,LL)*PM(I,2)*PM(J,2) + POMEGA(4,LL)*PM(I,4)*PM(J,4) &
          +POMEGA(6,LL)*PM(I,6)*PM(J,6) + POMEGA(8,LL)*PM(I,8)*PM(J,8)
          B(I,J,LL) =  - SUM0*WT(J)
          B(J,I,LL) =  - SUM0*WT(I)
         enddo
        enddo
        do I = 1,M_
          B(I,I,LL) = B(I,I,LL) + 1.0d0
        enddo
       enddo

!---odd-layer j-points, Legendre terms 1,3,5,7
       do LL=3,ND-2,2
        DELTAU = ZTAU(LL+1) - ZTAU(LL-1)
        do I = 1,M_
          A(I,LL) = EMU(I)/DELTAU
          C(I,LL) = -A(I,LL)
          H(I,LL) = FZ(LL)*( &
           POMEGA(1,LL)*PM(I,1)*PM0(1) + POMEGA(3,LL)*PM(I,3)*PM0(3) &
         + POMEGA(5,LL)*PM(I,5)*PM0(5) + POMEGA(7,LL)*PM(I,7)*PM0(7))
        enddo
        do I = 1,M_
         do J=1,I
          SUM0 = &
           POMEGA(1,LL)*PM(I,1)*PM(J,1) + POMEGA(3,LL)*PM(I,3)*PM(J,3) &
          +POMEGA(5,LL)*PM(I,5)*PM(J,5) + POMEGA(7,LL)*PM(I,7)*PM(J,7)
          B(I,J,LL) =  - SUM0*WT(J)
          B(J,I,LL) =  - SUM0*WT(I)
         enddo
        enddo
        do I = 1,M_
          B(I,I,LL) = B(I,I,LL) + 1.0d0
        enddo
       enddo

!---------lower boundary:  2nd-order terms
       L1 = ND
       L2 = ND-1
       do I = 1,M_
        SUM0 = &
         POMEGA(1,L1)*PM(I,1)*PM0(1) + POMEGA(3,L1)*PM(I,3)*PM0(3) &
       + POMEGA(5,L1)*PM(I,5)*PM0(5) + POMEGA(7,L1)*PM(I,7)*PM0(7)
        SUM2 = &
         POMEGA(1,L2)*PM(I,1)*PM0(1) + POMEGA(3,L2)*PM(I,3)*PM0(3) &
       + POMEGA(5,L2)*PM(I,5)*PM0(5) + POMEGA(7,L2)*PM(I,7)*PM0(7)
        SUM1 = &
         POMEGA(2,L1)*PM(I,2)*PM0(2) + POMEGA(4,L1)*PM(I,4)*PM0(4) &
       + POMEGA(6,L1)*PM(I,6)*PM0(6) + POMEGA(8,L1)*PM(I,8)*PM0(8)
        SUM3 = &
         POMEGA(2,L2)*PM(I,2)*PM0(2) + POMEGA(4,L2)*PM(I,4)*PM0(4) &
       + POMEGA(6,L2)*PM(I,6)*PM0(6) + POMEGA(8,L2)*PM(I,8)*PM0(8)
         H(I,L1) = 0.5d0*(SUM0*FZ(L1) + SUM2*FZ(L2))
         A(I,L1) = 0.5d0*(SUM1*FZ(L1) + SUM3*FZ(L2))
       enddo

       do I = 1,M_
        do J = 1,I
         SUM0 = &
          POMEGA(1,L1)*PM(I,1)*PM(J,1) + POMEGA(3,L1)*PM(I,3)*PM(J,3) &
        + POMEGA(5,L1)*PM(I,5)*PM(J,5) + POMEGA(7,L1)*PM(I,7)*PM(J,7)
         SUM2 = &
          POMEGA(1,L2)*PM(I,1)*PM(J,1) + POMEGA(3,L2)*PM(I,3)*PM(J,3) &
        + POMEGA(5,L2)*PM(I,5)*PM(J,5) + POMEGA(7,L2)*PM(I,7)*PM(J,7)
         SUM1 = &
          POMEGA(2,L1)*PM(I,2)*PM(J,2) + POMEGA(4,L1)*PM(I,4)*PM(J,4) &
        + POMEGA(6,L1)*PM(I,6)*PM(J,6) + POMEGA(8,L1)*PM(I,8)*PM(J,8)
         SUM3 = &
          POMEGA(2,L2)*PM(I,2)*PM(J,2) + POMEGA(4,L2)*PM(I,4)*PM(J,4) &
        + POMEGA(6,L2)*PM(I,6)*PM(J,6) + POMEGA(8,L2)*PM(I,8)*PM(J,8)
         S(I,J) = - SUM2*WT(J)
         S(J,I) = - SUM2*WT(I)
         T(I,J) = - SUM1*WT(J)
         T(J,I) = - SUM1*WT(I)
         V(I,J) = - SUM3*WT(J)
         V(J,I) = - SUM3*WT(I)
         B(I,J,L1) = - 0.5d0*(SUM0 + SUM2)*WT(J)
         B(J,I,L1) = - 0.5d0*(SUM0 + SUM2)*WT(I)
        enddo
       enddo

       do I = 1,M_
         S(I,I)   = S(I,I)   + 1.0d0
         T(I,I)   = T(I,I)   + 1.0d0
         V(I,I)   = V(I,I)   + 1.0d0
         B(I,I,L1)= B(I,I,L1) + 1.0d0

         C(I,L1)= S(I,1)*A(1,L1)/EMU(1) + S(I,2)*A(2,L1)/EMU(2) &
                + S(I,3)*A(3,L1)/EMU(3) + S(I,4)*A(4,L1)/EMU(4)
       enddo

       do I = 1,M_
        do J = 1,M_
         W(J,I) = S(J,1)*T(1,I)/EMU(1) + S(J,2)*T(2,I)/EMU(2) &
                + S(J,3)*T(3,I)/EMU(3) + S(J,4)*T(4,I)/EMU(4)
         U(J,I) = S(J,1)*V(1,I)/EMU(1) + S(J,2)*V(2,I)/EMU(2) &
                + S(J,3)*V(3,I)/EMU(3) + S(J,4)*V(4,I)/EMU(4)
        enddo
       enddo

!------------lower boundary, 2nd-order, A-matrix is full (AA)
         DELTAU = ZTAU(L1) - ZTAU(L2)
         D2 = 0.25d0*DELTAU
         SURFAC = 4.0d0*RFL/(1.0d0 + RFL)
       do I = 1,M_
          D1 = EMU(I)/DELTAU
          SUM0 = D1 + D2*(W(I,1)+W(I,2)+W(I,3)+W(I,4))
          SUM1 = SURFAC*SUM0
        do J = 1,M_
         AA(I,J,L1) = - D2*U(I,J)
         B(I,J,L1) = B(I,J,L1) + D2*W(I,J) - SUM1*EMU(J)*WT(J)
        enddo
         H(I,L1) = H(I,L1) - 2.0d0*D2*C(I,L1) + SUM0*ZFLUX
       enddo

       do I = 1,M_
          D1 = EMU(I)/DELTAU
        AA(I,I,L1) = AA(I,I,L1) + D1
        B(I,I,L1)  = B(I,I,L1) + D1
        C(I,L1) = 0.0d0
       enddo
      return
      END SUBROUTINE GEN_ID

!<<<<<<<<<<<<<<<<<<<<<<<<<<end fastJX scattering code<<<<<<<<<<<<<<<<<<<<<<<<<

!<<<<<begin fastJX subroutines called from PHOTO_JX or OPMIE<<<<<<<<<<<<

!------------------------------------------------------------------------------
      subroutine OPTICL (OPTD,SSALB,SLEG, ODCLD,NDCLD)
!------------------------------------------------------------------------------
!---set CLOUD fast-JX properties at the std 5 wavelengths:200-300-400-600-999nm
!      ----FJX ver 7.0+ clouds separated
! 01 W_C02   (C1/Deir)GAMMA:r-m=2.0/alf=6 n=1.335   reff=3.000___G=19.55_rho=1.000
! 02 W_C04   (C1/Deir)GAMMA:r-m=4.0/alf=6 n=1.335   reff=6.000___G=78.19_rho=1.000
! 03 W_C08   (C1/Deir)GAMMA:r-m=8.0/alf=2 n=1.335   reff=12.00___G=301.1_rho=1.000
! 04 W_C13   (C1/Deir)GAMMA:r-m=13./alf=2 n=1.335   reff=20.00___G=472.9_rho=1.000
! 05 W_L06   (W/Lacis)GAMMA:r-m=5.5/alf=11/3        reff=10.00___G=183.9_rho=1.000
! 06 Ice-Hexagonal (Mishchencko)                    reff=50.00___G=999.9_rho=0.917
! 07 Ice-Irregular (Mishchencko)                    reff=50.00___G=999.9_rho=0.917
! 08 S-Bkg   LOGN:r=.090 s=.600 n=1.514/.../1.435   reff=0.221___G=.0523_rho=1.630
! 09 S-Vol   LOGN:r=.080 s=.800 n=1.514/.../1.435   reff=0.386___G=.0721_rho=1.630
!
      implicit none

      real*8, intent(out)::    OPTD(5)    ! optical depth of layer
      real*8, intent(out)::    SSALB(5)   ! single-scattering albedo
      real*8, intent(out)::    SLEG(8,5)  ! scatt phase fn (Leg coeffs)
      real*8, intent(in)::     ODCLD      ! optical depth of cloud layer @600 nm
      integer,intent(inout)::  NDCLD      ! index of cloud layer:  4:13

      integer I,J
      real*8  XTINCT, REFF,RHO

!---later versions should allow for interpolation, averaging of R-eff

!---default cloud type C1, Reff = 12 microns
      if (NDCLD .lt. 1 .or. NDCLD .gt. 9) then
         NDCLD = 3
      endif

!--rescale OD by Qext at 600 nm (J=4)
      do J=1,5
         OPTD(J) = ODCLD * QCC(J,NDCLD)/QCC(4,NDCLD)
         SSALB(J) = SCC(J,NDCLD)
        do I=1,8
         SLEG(I,J) =  PCC(I,J,NDCLD)
        enddo
      enddo
      return
      END SUBROUTINE OPTICL



!------------------------------------------------------------------------------
      subroutine OPTICA (OPTD,SSALB,SLEG, PATH,RELH,K)
!------------------------------------------------------------------------------
!---for the UCI aerosol data sets, calculates optical properties at fast-JX's
!              std 5 wavelengths:200-300-400-600-999nm
!
!---UCI aersols optical data  v-7.0+
! 01 S-Bkg   LOGN:r=.090 s=.600 n=1.514/.../1.435   reff=0.221___G=.0523_rho=1.630
! 02 S-Vol   LOGN:r=.080 s=.800 n=1.514/.../1.435   reff=0.386___G=.0721_rho=1.630
! 03 UT-sulfate LOGN:r=0.05 s=.693 n=1.44           reff=0.166___G=.0205_rho=1.769
! 04 UT-sulfate LOGN:r=0.05 s=.693 n=1.46           reff=0.166___G=.0205_rho=1.769
! 05 UT-sulfatM LOGN:r=.050 s=.642 n=1.53           reff=0.140___G=.0179_rho=1.769
! 06 UM-BC1     LOGN:r=.050 s=.642 n=1.80+0.50i     reff=0.140___G=.0179_rho=1.500
! 07 UM-BC2     LOGN:r=.080 s=.501 n=1.80+0.50i     reff=0.150___G=.0332_rho=1.500
! 08 UM-BB08 (%BC)LOGN:r=.080 s=.500 n=1.552+0.04i  reff=0.149___G=.0331_rho=1.230
! 09 UM-FF04 (%BC) LOGN:r=.050 s=.642 n=1.541+0.02i reff=0.140___G=.0179_rho=1.212
! 10 UM-FF10 (%BC)LOGN:r=.050 s=.642 n=1.557+0.05i  reff=0.140___G=.0179_rho=1.230
! 11 MDust.15  (R.V. Martin generated phase fns)    reff=0.150___G=1.000_rho=2.600
! 12 MDust.25  (R.V. Martin generated phase fns)    reff=0.250___G=1.000_rho=2.600
! 13 MDust.40  (R.V. Martin generated phase fns)    reff=0.400___G=1.000_rho=2.600
! 14 MDust.80  (R.V. Martin generated phase fns)    reff=0.800___G=1.000_rho=2.6001
! 15 MDust1.5  (R.V. Martin generated phase fns)    reff=1.500___G=1.000_rho=2.600
! 16 MDust2.5  (R.V. Martin generated phase fns)    reff=2.500___G=1.000_rho=2.600
! 17 MDust4.0  (R.V. Martin generated phase fns)    reff=4.000___G=1.000_rho=2.600
!
      implicit none

      real*8, intent(out)::    OPTD(5)    ! optical depth of layer
      real*8, intent(out)::    SSALB(5)   ! single-scattering albedo
      real*8, intent(out)::    SLEG(8,5)  ! scatt phase fn (Leg coeffs)
      real*8, intent(in)::     PATH       ! path (g/m2) of aerosol/cloud
      real*8, intent(in)::     RELH       ! relative humidity (0.00->1.00+)
      integer,intent(inout)::     K       ! index of cloud/aerosols

      integer I,J
      real*8  XTINCT, REFF,RHO

      if (K .gt. NAA .or. K .lt. 1) then
         write(*,*) ' aerosol index out-of-range: K/NAA',K,NAA
         K = 18
      endif

         REFF = RAA(K)
         RHO = DAA(K)
      do J=1,5
!---extinction K(m2/g) = Q(wvl) / [4/3 * Reff(micron) * aerosol-density(g/cm3)]
         XTINCT = 0.75d0*QAA(J,K)/(REFF*RHO)
         OPTD(J) = PATH*XTINCT
         SSALB(J) = SAA(J,K)
       do I=1,8
         SLEG(I,J) =  PAA(I,J,K)
       enddo
      enddo
      return
      END SUBROUTINE OPTICA

!------------------------------------------------------------------------------
      subroutine OPTICM (OPTD,SSALB,SLEG, PATH,RELH,LL)
!------------------------------------------------------------------------------
!---for the U Michigan aerosol data sets, this generate fast-JX data formats.
!---NB Approximates the Legendre expansion(L) of the scattering phase fn
!---           as (2*L+1)*g**L
!---UMAER(I,J,K,L):
!   I=1:3 = [SSAbldeo, g, k-ext(m2/g)]
!   J=1:5 = [200, 300, 400, (550,) 600 , 1000 nm]
!   K=1:21= [0, 5, 10, 15, ..., 90, 95, 99 %RelHum]
!   L=1:33= UM aerosol types [SULF, SS-1,-2,-3,-4, DD-1,-2,-3,-4, FF00(0%BC),
!                      FF02, ...FF14(14%BC), BB00, BB02, ...BB30(30%BC)]
      implicit none

      real*8, intent(out)::    OPTD(5)    ! optical depth of layer
      real*8, intent(out)::    SSALB(5)   ! single-scattering albedo
      real*8, intent(out)::    SLEG(8,5)  ! scatt phase fn (Leg coeffs)
      real*8, intent(in)::     PATH       ! path (g/m2) of aerosol/cloud
      real*8, intent(in)::     RELH       ! relative humidity (0.00->1.00)
      integer,intent(in)::     LL         ! index of cloud/aerosols

      integer KR,J,L
      real*8  R,FRH, GCOS, XTINCT

!---calculate fast-JX properties at the std 5 wavelengths:200-300-400-600-999nm
!---extrapolate phase fn from first term (g)
      L = LL
      if (L .lt. 1 .or. L .gt. 33) then
!ccc         write(6,*) ' UM aer index too large: L',L
         L = 1
      endif

!---pick nearest Relative Humidity
      KR =  20.d0*RELH  + 1.5d0
      KR = max(1, min(21, KR))

      do J=1,5
       SSALB(J) = UMAER(1,J,KR,L)
         XTINCT = UMAER(3,J,KR,L)
       OPTD(J) = PATH*XTINCT
         GCOS   = UMAER(2,J,KR,L)
       SLEG(1,J) =  1.d0
       SLEG(2,J) =  3.d0*GCOS
       SLEG(3,J) =  5.d0*GCOS**2
       SLEG(4,J) =  7.d0*GCOS**3
       SLEG(5,J) =  9.d0*GCOS**4
       SLEG(6,J) = 11.d0*GCOS**5
       SLEG(7,J) = 13.d0*GCOS**6
       SLEG(8,J) = 15.d0*GCOS**7
      enddo
      return
      END SUBROUTINE OPTICM

!-----------------------------------------------------------------------
      subroutine JRATET(PPJ,TTJ,FFF, VALJL,LU,NJXU)
!-----------------------------------------------------------------------
! in:
!        PPJ(L_+1) = pressure profile at edges
!        TTJ(L_+1) = = temperatures at mid-level
!        FFF(K=1:NW, L=1:L_) = mean actinic flux
! out:
!        VALJL(L_,JX_)  JX_ = no of dimensioned J-values in CTM code
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in)  :: LU,NJXU
      real*8, intent(in)  ::  PPJ(LU),TTJ(LU)             ! use pfull and tfull
      real*8, intent(inout)  ::  FFF(W_,LU)
      real*8, intent(out), dimension(LU,NJXU) ::  VALJL

      real*8  VALJ(NJX_)
      real*8  QO2TOT, QO3TOT, QO31DY, QO31D, QQQT, TFACT
      real*8  TT,PP,DD,TT200,TFACA,TFAC0,TFAC1,TFAC2,QQQA,QQ2,QQ1A,QQ1B
      integer J,K,L, IV
      if (NJXU .lt. NJX) then
        write(*,*) ' JRATET:  CTM has not enough J-values dimensioned'
        stop
      endif

      do L = 1,LU
!---need temperature, pressure, and density at mid-layer (for some quantum yields):
          TT   = TTJ(L1_-L)
          PP   = PPJ(L1_-L)
!         if (L .eq. 1) then
!          PP = PPJ(1)
!         else
!          PP  = (PPJ(L)+PPJ(L+1))*0.5d0
!         endif
          DD = 7.24e18*PP/TT

!---if W_=18/12, must zero bin-11/5 below 100 hPa, since O2 e-fold is too weak
!        and does not represent the decay of 215.5-221.5 nm sunlight.
!         write(logunit,*) 'LI===: ', TT, PP, DD
!         write(*,*) 'LI===: ', TT, PP, DD
        if (PP .gt. 100.d0) then
          if (W_ .eq. 18) then
            FFF(11,L) = 0.d0
          elseif (W_ .eq. 12) then
            FFF(5,L) = 0.d0
          endif
        endif

!        write(logunit,*) 'JRATET NJX: ', NJX
!        write(*,*) 'JRATET NJX: ', NJX
        do J = 1,NJX
          VALJ(J) = 0.d0
        enddo

         do K = 1,W_
          call X_interp (TT,QO2TOT, TQQ(1,1),QO2(K,1), &
                  TQQ(2,1),QO2(K,2), TQQ(3,1),QO2(K,3), LQQ(1))
          call X_interp (TT,QO3TOT, TQQ(1,2),QO3(K,1), &
                  TQQ(2,2),QO3(K,2), TQQ(3,2),QO3(K,3), LQQ(2))
          call X_interp (TT,QO31DY, TQQ(1,3),Q1D(K,1), &
                  TQQ(2,3),Q1D(K,2), TQQ(3,3),Q1D(K,3), LQQ(3))
            QO31D  = QO31DY*QO3TOT
           VALJ(1) = VALJ(1) + QO2TOT*FFF(K,L)
           VALJ(2) = VALJ(2) + QO3TOT*FFF(K,L)
           VALJ(3) = VALJ(3) + QO31D*FFF(K,L)
!           write(logunit,*) 'LI==> ', K, 'VALJ(1-3) ', VALJ(1), VALJ(2), VALJ(3) 
         enddo

        do J = 4,NJX
         do K = 1,W_
!---also need to allow for Pressure interpolation if SQQ(J) = 'p'
          if (SQQ(J) .eq.'p') then
           call X_interp (PP,QQQT, TQQ(1,J),QQQ(K,1,J), &
              TQQ(2,J),QQQ(K,2,J), TQQ(3,J),QQQ(K,3,J), LQQ(J))
          else
           call X_interp (TT,QQQT, TQQ(1,J),QQQ(K,1,J), &
              TQQ(2,J),QQQ(K,2,J), TQQ(3,J),QQQ(K,3,J), LQQ(J))
          endif
            VALJ(J) = VALJ(J) + QQQT*FFF(K,L)
!            write(logunit,*) 'LI==> ', J, K, 'VALJ(J) ', VALJ(J)
            if(VALJ(J) .gt. 1.0d0) then
               write(*,*) 'Problems ', J, K, 'VALJ(J) ', VALJ(J), QQQT, FFF(K,L)
            endif
         enddo
        enddo

        do J=1,NJX
          VALJL(L,J) = VALJ(J)
        enddo

      enddo
      return
      END SUBROUTINE JRATET

!-----------------------------------------------------------------------
      subroutine X_interp (TINT,XINT, T1,X1, T2,X2, T3,X3, L123)
!-----------------------------------------------------------------------
!  up-to-three-point linear interpolation function for X-sections
!-----------------------------------------------------------------------
      implicit none

      real*8, intent(in)::  TINT,T1,T2,T3, X1,X2,X3
      integer,intent(in)::  L123
      real*8,intent(out)::  XINT

      real*8  TFACT

      if (L123 .le. 1) then
           XINT = X1
      elseif (L123 .eq. 2) then
             TFACT = max(0.d0,min(1.d0,(TINT-T1)/(T2-T1) ))
           XINT = X1 + TFACT*(X2 - X1)
      else
        if (TINT.le. T2) then
             TFACT = max(0.d0,min(1.d0,(TINT-T1)/(T2-T1) ))
           XINT = X1 + TFACT*(X2 - X1)
        else
             TFACT = max(0.d0,min(1.d0,(TINT-T2)/(T3-T2) ))
           XINT = X2 + TFACT*(X3 - X2)
        endif
      endif
      return
      END SUBROUTINE X_interp


!-----------------------------------------------------------------------
      subroutine JP_ATM(PPJ,TTJ,DDJ,OOJ,ZZJ,DTAU6,POMEG6,JXTRA,LU)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!---the CTM has L_ = LU layers and fast-JX adds layer LU+1
!---the pressure and altitude(Z) are on layer edge (LU+2)

      integer,intent(in)                  :: LU
      real*8, intent(in), dimension(LU+2) :: PPJ,ZZJ
      real*8, intent(in), dimension(LU+1) :: TTJ,DDJ,OOJ,DTAU6
      real*8, intent(in), dimension(8,LU+1) :: POMEG6
      integer,intent(in), dimension(LU+LU+3) :: JXTRA
!-----------------------------------------------------------------------
      integer  I,J,K,L
      real*8   COLO2,COLO3,ZKM,DELZ,ZTOP

!      write(*,'(4a)') '   L z(km)     p      T   ', &
!       '    d(air)   d(O3)','  col(O2)  col(O3)     d-TAU   SS-alb', &
!       '  g(cos) CTM lyr=>'

      L = LU+2
!      write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)') &
!            L,ZZJ(L)*1.d-5,PPJ(L)

          COLO2 = 0.d0
          COLO3 = 0.d0
          ZTOP = ZZJ(LU+2)

        do L = LU+1,1,-1
          COLO2 = COLO2 + DDJ(L)*0.20948d0
          COLO3 = COLO3 + OOJ(L)
          DELZ = ZTOP-ZZJ(L)
          ZTOP = ZZJ(L)
          ZKM = ZZJ(L)*1.d-5

!      write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)') &
!            L,ZKM,PPJ(L),TTJ(L),DDJ(L)/DELZ,OOJ(L)/DELZ, &
!            COLO2,COLO3,DTAU6(L),POMEG6(1,L),POMEG6(2,L)/3.d0, &
!            JXTRA(L+L),JXTRA(L+L-1)

        enddo
      return
      END SUBROUTINE JP_ATM

!-----------------------------------------------------------------------
      subroutine JP_ATM0(PPJ,TTJ,DDJ,OOJ,ZZJ, LU)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!---the CTM has L_ = LU layers and fast-JX adds layer LU+1
!---the pressure and altitude(Z) are on layer edge (LU+2)

      integer,intent(in)                  :: LU
      real*8, intent(in), dimension(LU+2) :: PPJ,ZZJ
      real*8, intent(in), dimension(LU+1) :: TTJ,DDJ,OOJ
!-----------------------------------------------------------------------
      integer  I,J,K,L
      real*8   COLO2,COLO3,ZKM,DELZ,ZTOP
!      write(*,'(4a)') '   L z(km)     p      T   ', &
!       '    d(air)   d(O3)','  col(O2)  col(O3)     d-TAU   SS-alb', &
!       '  g(cos) CTM lyr=>'
      L = LU+2
!      write(*,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)') &
!            L,ZZJ(L)*1.d-5,PPJ(L)
          COLO2 = 0.d0
          COLO3 = 0.d0
          ZTOP = ZZJ(LU+2)
        do L = LU+1,1,-1
          COLO2 = COLO2 + DDJ(L)*0.20948d0
          COLO3 = COLO3 + OOJ(L)
          DELZ = ZTOP-ZZJ(L)
          ZTOP = ZZJ(L)
          ZKM = ZZJ(L)*1.d-5
!      write(*,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)') &
!            L,ZKM,PPJ(L),TTJ(L),DDJ(L)/DELZ,OOJ(L)/DELZ, &
!            COLO2,COLO3
        enddo
      return
      END SUBROUTINE JP_ATM0

!-----------------------------------------------------------------------
      subroutine SPHERE2(U0,RAD,ZHL,ZZHT,AMF2, L1U,LJX1U)
!-----------------------------------------------------------------------
!----new v6.2: does AirMassFactors for mid-layer, needed for SZA ~ 90
!  This new AMF2 does each of the half-layers of the CTM separately,
!     whereas the original, based on the pratmo code did the whole layers
!     and thus calculated the ray-path to the CTM layre edges, NOT the middle.
!  Since fast-JX is meant to calculate the intensity at the mid-layer, the
!     solar beam at low sun (interpolated between layer edges) was incorrect.
!  This new model does make some approximations of the geometry of the layers:
!     the CTM layer is split evenly in mass (good) and in height (approx).
!
!  Calculation of spherical geometry; derive tangent heights, slant path
!  lengths and air mass factor for each layer. Not called when
!  SZA > 98 degrees.  Beyond 90 degrees, include treatment of emergent
!  beam (where tangent height is below altitude J-value desired at).
!-----------------------------------------------------------------------
! in:
!     U0      cos(solar zenith angle)
!     RAD     radius of Earth mean sea level (cm)
!     ZHL(L)  height (cm) of the bottome edge of CTM level L
!     ZZHT    scale height (cm) used above top of CTM (ZHL(L_+1))
!     L1U     dimension of CTM = levels +1 (L+1 = above-CTM level)
! out:
!     AMF2(I,J) = air mass factor for CTM level I for sunlight reaching J
!          these are calcualted for both layer middle and layer edge
!-----------------------------------------------------------------------
      implicit none
      integer, intent(in) ::   L1U, LJX1U
      real*8, intent(in)  ::   U0,RAD,ZHL(L1U+1),ZZHT
      real*8, intent(out) ::   AMF2(2*LJX1U+1,2*LJX1U+1)

      integer, parameter  ::  LSPH_ = 100

!     RZ      Distance from centre of Earth to each point (cm)
!     RQ      Square of radius ratios
!     SHADHT  Shadow height for the current SZA
!     XL      Slant path between points

      integer  I, J, K, II, L2
      real*8   XMU1,XMU2,XL,DIFF,SHADHT,RZ(LSPH_+1)
      real*8   RZ2(2*LSPH_+1),RQ2(2*LSPH_+1)

!--- must have top-of-atmos (NOT top-of-CTM) defined
!      ZHL(L1U+1) = ZHL(L1U) + ZZHT

      if (L1U .gt. LSPH_) then
        write(*,*)' SPHERE2: temp arrays not large enough', L1U, LSPH_
        stop
      endif

        RZ(1) = RAD + ZHL(1)
      do II = 2,L1U+1
        RZ(II)   = RAD + ZHL(II)
      enddo

!---calculate heights for edges of split CTM-layers
     L2 = 2*L1U
      do II = 2,L2,2
        I = II/2
        RZ2(II-1) = RZ(I)
        RZ2(II) = 0.5d0*(RZ(I)+RZ(I+1))
      enddo
        RZ2(L2+1) = RZ(L1U+1)
      do II = 1,L2
        RQ2(II) = (RZ2(II)/RZ2(II+1))**2
      enddo

!---shadow height for SZA > 90
      if (U0 .lt. 0.0d0)  then
        SHADHT = RZ2(1)/dsqrt(1.0d0 - U0**2)
      else
        SHADHT = 0.d0
      endif

!---up from the surface calculating the slant paths between each level
!---  and the level above, and deriving the appropriate Air Mass Factor
         AMF2(:,:) = 0.d0

      do 16 J = 1,2*L1U+1

!  Air Mass Factors all zero if below the tangent height
        if (RZ2(J) .lt. SHADHT) goto 16
!  Ascend from layer J calculating AMF2s
        XMU1 = abs(U0)
        do I = J,2*L1U
          XMU2     = dsqrt(1.0d0 - RQ2(I)*(1.0d0-XMU1**2))
          XL       = RZ2(I+1)*XMU2 - RZ2(I)*XMU1
          AMF2(I,J) = XL / (RZ2(I+1)-RZ2(I))
          XMU1     = XMU2
        enddo
!--fix above top-of-atmos (L=L1U+1), must set DTAU(L1U+1)=0
          AMF2(2*L1U+1,J) = 1.d0
!
!  Twilight case - Emergent Beam, calc air mass factors below layer
        if (U0 .ge. 0.0d0) goto 16

!  Descend from layer J
          XMU1       = abs(U0)
         do II = J-1,1,-1
          DIFF        = RZ2(II+1)*sqrt(1.0d0-XMU1**2)-RZ2(II)
          if (II.eq.1)  DIFF = max(DIFF,0.d0)   ! filter
!  Tangent height below current level - beam passes through twice
          if (DIFF .lt. 0.0d0)  then
            XMU2      = sqrt(1.0d0 - (1.0d0-XMU1**2)/RQ2(II))
            XL        = abs(RZ2(II+1)*XMU1-RZ2(II)*XMU2)
            AMF2(II,J) = 2.d0*XL/(RZ2(II+1)-RZ2(II))
            XMU1      = XMU2
!  Lowest level intersected by emergent beam
          else
            XL        = RZ2(II+1)*XMU1*2.0d0
            AMF2(II,J) = XL/(RZ2(II+1)-RZ2(II))
            goto 16
          endif
         enddo

   16 continue
       
     return
     END SUBROUTINE SPHERE2


!-----------------------------------------------------------------------
      subroutine EXTRAL(DTAUX,L1X,L2X,NX,JTAUMX,ATAU,ATAU0, JXTRA)
!-----------------------------------------------------------------------
!
!    new version 6.1, add sub-layers (JXTRA) to thick cloud/aerosol layers
!    this version sets up log-spaced sub-layers of increasing thickness ATAU
!
!     DTAUX(L=1:L1X) = Optical Depth in layer L (generally 600 nm OD)
!        This can be just cloud or cloud+aerosol, it is used only to set
!        the number in levels to insert in each layer L
!        Set for log-spacing of tau levels, increasing top-down.
!
!     N.B. the TTAU, etc calculated here are NOT used elsewhere

!---The log-spacing parameters have been tested for convergence and chosen
!---  to be within 0.5% for ranges OD=1-500, rflect=0-100%, mu0=0.1-1.0
!---  use of ATAU = 1.18 and min = 0.01, gives at most +135 pts for OD=100
!---  ATAU = 1.12 now recommended for more -accurate heating rates (not J's)
!-----------------------------------------------------------------------
!
      implicit none
      integer, intent(in) ::  JTAUMX,L1X,L2X  !index of cloud/aerosol
      integer, intent(in) ::  NX              !Mie scattering array size
      real*8,  intent(in) ::  DTAUX(L1X)      !cloud+3aerosol OD in each layer
      real*8,  intent(in) ::  ATAU,ATAU0
      integer, intent(out)::  JXTRA(L2X+1)    !number of sub-layers to be added
!
      integer JTOTL,I,L,L2
      real*8  TTAU(L2X+1),DTAUJ, ATAU1,ATAULN,ATAUM,ATAUN1
!
!---Reinitialize arrays
      TTAU(:)  = 0.d0
      JXTRA(:) = 0
!
!---combine these edge- and mid-layer points into grid of size:
!---              L2X+1 = 2*L1X+1 = 2*L_+3
!---calculate column optical depths above each level, TTAU(1:L2X+1)
!---      note that TTAU(L2X+1)=0 and TTAU(1)=total OD
!
!---Divide thick layers to achieve better accuracy in the scattering code
!---In the original fast-J, equal sub-layers were chosen, this is wasteful
!---and this new code (ver 5.3) uses log-scale:
!---        Each succesive layer (down) increase thickness by ATAU > 1
!---        e.g., if ATAU = 2, a layer with OD = 15 could be divided into
!---        4 sub-layers with ODs = 1 - 2 - 4 - 8
!---The key parameters are:
!---        ATAU = factor increase from one layer to the next
!---        ATAUMN = the smallest OD layer desired
!---        JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)
!---These are hardwired below, can be changed, but have been tested/optimized

      ATAU1  = ATAU - 1.d0
      ATAULN = log(ATAU)
        TTAU(L2X+1)  = 0.0d0
      do L2 = L2X,1,-1
        L         = (L2+1)/2
        DTAUJ     = 0.5d0 * DTAUX(L)
        TTAU(L2)  = TTAU(L2+1) + DTAUJ
!---Now compute the number of log-spaced sub-layers to be added in
!---   the interval TTAU(L2) > TTAU(L2+1)
!---The objective is to have successive TAU-layers increasing by factor ATAU >1
!---the number of sub-layers + 1
        if (TTAU(L2) .lt. ATAU0) then
          JXTRA(L2) = 0
        else
          ATAUM    = max(ATAU0, TTAU(L2+1))
          ATAUN1 = log(TTAU(L2)/ATAUM) / ATAULN
          JXTRA(L2) = min(JTAUMX, max(0, int(ATAUN1 - 0.5d0)))
        endif
      enddo

!---check on overflow of arrays, cut off JXTRA at lower L if too many levels
      JTOTL    = L2X + 2
      do L2 = L2X,1,-1
        JTOTL  = JTOTL + JXTRA(L2)
        if (JTOTL .gt. NX/2)  then
!          write(*,'(A,2I5,F9.2)') 'N_/L2_/L2-cutoff JXTRA:',NX,L2X,L2
          do L = L2,1,-1
            JXTRA(L) = 0
          enddo
          go to 10
        endif
      enddo
  10  continue

      return
      END SUBROUTINE EXTRAL

!<<<<<<<end fastJX subroutines called from PHOTO_JX or OPMIE<<<<<<<<<<<<


                                                                       
!-----------------------------------------------------------------------
      subroutine OPTICAM3 (OPTD,SSALB,SLEG, PATH,L) 
!-----------------------------------------------------------------------
!---for the AM3 aerosol data sets, calculates optical properties at fast-JX
!              std 5 wavelengths:200-300-400-600-999nm                  :                                    
!                                                                         
!                    
!                                                                        
      implicit none 
!      include 'parm_CTM.f' 
!      include 'parm_MIE.f' 
!      include 'cmn_JVdat.f' 
                                                                        
                                          ! optical depth of layer      
      real*8, intent(out)::    OPTD(5) 
                                          ! single-scattering albedo    
      real*8, intent(out)::    SSALB(5) 
                                          ! scatt phase fn (Leg coeffs) 
      real*8, intent(out)::    SLEG(8,5) 
                                          ! path (g/m2) of aerosol/cloud
      real*8, intent(in)::     PATH 
                                          ! index of cloud/aerosols     
      integer,intent(inout)::     L 
                                                                        
      integer I,J 
!      real*8  XTINCT, REFF,RHO 
                                                                        
      if (L .gt. NAA_AM3 ) then 
         write(*,*) 'ATMOS:fastjx_photo:OPTICAM3:FATAL: aerosol index out-of-range: L/NAA',L,NAA_AM3
         call error_mesg ('ATMOS: fastjx_init:OPTICAM3','aerosol index out-of-range.', FATAL)
      endif 
                                                                        
!          MEE1 = MEE(L)                                                          
!         REFF = RAA(L) 
!         RHO = DAA(L) 
      do J=1,5 
!---extinction K(m2/g) = Q(wvl) / [4/3 * Reff(micron) * aerosol-density(
!         XTINCT = 0.75d0*QAA(J,L)/(REFF*RHO) 
         OPTD(J) = PATH*MEE_AM3(J,L) 
         SSALB(J) = SAA_AM3(J,L) 
       do I=1,8 
         SLEG(I,J) =  PAA_AM3(I,J,L) 
       enddo 
      enddo 
                                                                        
      return 
      END  subroutine OPTICAM3                                         
 
!-----------------------------------------------------------------------
                                                                        
                                                                        
!-----------------------------------------------------------------------
      subroutine RD_MIE_AM3(NAMFIL) 
!-----------------------------------------------------------------------
!-------aerosols/cloud scattering data set for fast-JX (ver 5.3+)       
!  >>>>>>>>>>>>>>>>spectral data rev to J-ref ver8.5 (5/05)<<<<<<<<<<<< 
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)        
!     NJ1      Channel number for reading data file                     
!     NAA      Number of categories for scattering phase functions      
!     QAA      Aerosol scattering phase functions                       
!     NK       Number of wavelengths at which functions supplied (set as
!     WAA      Wavelengths for the NK supplied phase functions          
!     PAA      Phase function: first 8 terms of expansion               
!     RAA      Effective radius associated with aerosol type            
!     SAA      Single scattering albedo                                 
!-----------------------------------------------------------------------
      implicit none 
!      include 'parm_CTM.f' 
!      include 'parm_MIE.f' 
!      include 'cmn_JVdat.f' 
                                                                        
      character(*), intent(in) ::  NAMFIL 
      integer  I, J, K, NJ1 
      character*90 :: str1
!jul++
!   NAA,TITLE0
!   JTAUMX,ATAU,ATAU0 
!   TITLAA(J),RAA(J),DAA(J) 
!   WAA(K,J),QAA(K,J),SAA(K,J),PAA(I,K,J)
!
!jul--                                                                        
      call mpp_open (NJ1, trim(NAMFIL), MPP_RDONLY, MPP_ASCII,  &
                     MPP_SEQUENTIAL, MPP_MULTI, MPP_SINGLE)
                                                                        
      read (NJ1,'(i4,a78)') NAA_AM3, TITLE0_AM3 
      if (NAA_AM3 .gt. A_AM3) then 
          write(*,*) 'ATMOS:fastjx_init: too many scat-data sets for AM3:', NAA_AM3, A_AM3 
          stop 
      endif 
!      read (NJ1,'(5x,i5,2f10.5)') JTAUMX,ATAU,ATAU0 
!      if (mpp_pe() == mpp_root_pe()) then 
!               write(*,*)'ATMOS:fastjx_init: RD_MIE_AM3:'
!               write(*,'(a,2f9.5,i5)') ' ATAU/ATAU0/JMX',ATAU,ATAU0,JTAUMX 
!      end if
        
!      read (NJ1,*) 
      do J = 1,NAA_AM3
          read (NJ1,'(a90)') str1
!          if (mpp_pe() == mpp_root_pe()) then 
!             write (*,'(a12,a90)')              &
!                   'fastjx_init:',str1               
!          endif       
!          read (NJ1,'(3x,a20,32x,f5.3,15x,f5.3)')                       &
!     &           TITLAA(J),RAA(J),DAA(J)                                
                       ! ver 6.0 extend to 5 ref wavelengths for mie-sca
        do K = 1,5 
          read (NJ1,'(f4.0,2e14.6,7f6.3,1x,f7.3,f8.4)')              &
                   WAA_AM3(K,J),MEE_AM3(K,J),SAA_AM3(K,J),(PAA_AM3(I,K,J),I=2,8)                   
          PAA_AM3(1,K,J) = 1.d0 
  
          if (mpp_pe() == mpp_root_pe()) then 
!!!              write (*,'(a12,i4,1x,f4.0,2e14.6,8f6.3,1x,f7.3,f8.4)')              &
!!!                  'fastjx_init:',J,WAA_AM3(K,J),MEE_AM3(K,J),SAA_AM3(K,J),(PAA_AM3(I,K,J),I=1,8)                   
          endif
        enddo 
      enddo 
                                                                        
       call mpp_close (NJ1)
!      if (mpp_pe() == mpp_root_pe()) then 
!        write(*,'(a,9f8.1)') ' ATMOS:fastjx_init: RD_MIE: Aerosol optical: r-eff/rho/Q(@wavel):'     &
!                  ,(WAA_AM3(K,1),K=1,5)                                    
!!        write(*,*) TITLE0_AM3 
!      end if                                                                        
!        if (mpp_pe() == mpp_root_pe()) then 

!      do J=1,NAA_AM3 
!          write(*,*) ' ATMOS:fastjx_init: RD_MIE_AM3:'
!           write(*,'(i3,1x,a8,7f8.3)')                                 &
!                   J,TITLAA_AM3(J),(QAA(K,J),K=1,5)   
!        end if          
!      enddo 
                                                                        
      return 
      END  subroutine RD_MIE_AM3                                         
 
                                                                        
!-----------------------------------------------------------------------
      real*8 FUNCTION FLINT (TINT,T1,T2,T3,F1,F2,F3) 
!-----------------------------------------------------------------------
!  Three-point linear interpolation function                            
!-----------------------------------------------------------------------
      real*8  TINT,T1,T2,T3,F1,F2,F3 
      if (TINT .le. T2)  then 
        if (TINT .le. T1)  then 
          FLINT = F1 
        else 
          FLINT = F1 + (F2 - F1)*(TINT -T1)/(T2 -T1) 
        endif 
      else 
        if (TINT .ge. T3)  then 
          FLINT = F3 
        else 
          FLINT = F2 + (F3 - F2)*(TINT -T2)/(T3 -T2) 
        endif 
      endif 
      return 
      END  FUNCTION FLINT                                         
end module MO_FASTJX_MOD
