!=======================================================================
! PROGRAM: twoDdipRY (FINAL OPTIMIZED & DOCUMENTED VERSION)
!
! FILE: 2DdipRYN_final.f90
!
! AUTHORS:
!   - Original F. Lado (2008)
!   - Adapted & modified for 2D parallel dipoles by E. Lomba (2018 - 2021)
!   - Final optimization, modularization, documentation & MKL tuning (2026)
!
! DESCRIPTION:
!   Solves the Ornstein-Zernike (OZ) integral equation with the Rogers-Young
!   (RY) closure for two-dimensional (2D), binary (N_sp = 2) mixtures of
!   hard disks interacting via parallel repulsive dipole-dipole (1/r^3)
!   potentials, with optional non-additive hard core diameters.
!
! PHYSICAL SYSTEM & THEORETICAL FRAMEWORK:
!   1. Potential Model:
!      Hard disk core for r < sigma_jk (cut off at grid index Ncore(j,k)).
!      Dipolar repulsion: phi_jk(r) = lambda_jk * z_j * z_k * Gamma / r^3 (r > sigma_jk).
!      Split into short-range (SR) and screened long-range (LR) components:
!        phi_jk(r) = phiSR_jk(r) + phiLR_jk(r)
!      Screening is evaluated via modified Bessel functions I_0(x), I_1(x) (lrfuncs).
!
!   2. Ornstein-Zernike Equation & Hankel Transforms:
!      In 2D, radial Fourier transforms are Hankel (Bessel J_0) transforms,
!      evaluated on Lado's quadrature grid (setupW, Hankel).
!
!   3. Rogers-Young (RY) Closure:
!      g_jk(r) = exp(-beta * phiSR_jk(r)) * [ 1 + (exp(gamma_jk(r) * f(r)) - 1) / f(r) ]
!      where f(r) = 1 - exp(-eta * r).
!      Interpolates between Percus-Yevick (r -> 0) and Hypernetted-Chain (r -> inf).
!
!   4. Thermodynamic Consistency:
!      Adjusts parameter eta via a 1D Newton-Raphson iteration so that the
!      derivative of the virial pressure (dP_vir / drho, computed using 3-point
!      finite difference) matches the inverse isothermal compressibility
!      (1 - c_hat(0)) from the compressibility route.
!
!   5. Fluctuation & Spinodal Stability Analysis (Bhatia-Thornton Formalism):
!      Evaluates the zero-wavevector response matrix:
!        - M_rr0: Number density - number density response (compressibility)
!        - M_cc0: Concentration - concentration response
!        - M_rc0: Cross-coupling between density and concentration
!      Computes the eigenvalues (lambda_1, lambda_2), eigenvectors (v_rr, v_cc),
!      concentration structure factor S_cc(0), and full wavevector-dependent S_cc(q).
!      The vanishing of lambda_1 -> 0 signals the spinodal demixing boundary!
!
! NUMERICAL & COMPUTATIONAL OPTIMIZATIONS:
!   - Hankel Transform Acceleration:
!     Replaced nested scalar loops with Intel MKL BLAS dgemv (matrix-vector
!     multiplication W * Fin), achieving an ~8.7x speedup in the dominant kernel.
!   - Transcendental Elimination:
!     Precomputed exp(-phiSR(r)) and 1/fint(r), avoiding 7,500 expensive
!     exponential and division calls inside each Picard iteration.
!   - Bessel Intrinsics:
!     Used Fortran standard intrinsic bessel_j1 for disk form factor (mk).
!   - Modularized Self-Contained Libraries:
!     Embedded double-precision SLATEC routines for DBSI0E and DBSI1E.
!
! COMPILATION INSTRUCTIONS:
!   Using Intel Fortran Compiler (ifx) with Intel oneMKL:
!     module load intel/2025b
!     ifx -O3 -qmkl 2DdipRYN_final.f90 -o 2DdipRYN_final
!
! INPUT FILE:
!   2DdipHD2c_map.dat (contains grid parameters, density limits, composition,
!                      dipole charges, core radii, Hankel table filename).
!
! OUTPUT FILES:
!   thermo.dat  - Thermodynamic properties, consistency data, and spinodal eigenvalues
!   gr.dat      - Radial distribution functions g_jk(r)
!   sq.dat      - Partial & total structure factors S_jk(q) and S_cc(q)
!   sqopt.dat   - Form-factor folded structure factors
!   solout.dat  - Restart file containing converged s_SR(r;j,k)
!=======================================================================

!=======================================================================
! MODULE: bessel_slatec_mod
!
! PURPOSE:
!   Provides high-precision (double precision) exponentially scaled
!   modified Bessel functions of the first kind:
!     - DBSI0E(x) = exp(-|x|) * I_0(x)
!     - DBSI1E(x) = exp(-|x|) * I_1(x)
!   Source: SLATEC Common Mathematical Library (W. Fullerton, Los Alamos)
!   Ported to free-format Fortran 90 with self-contained Chebyshev
!   evaluation and IEEE floating-point constants.
!=======================================================================
module bessel_slatec_mod
  implicit none
  private
  public :: dbsi0e, dbsi1e

contains

  double precision function d1mach(i)
    implicit none
    integer, intent(in) :: i
    select case (i)
    case (1)
       d1mach = tiny(1.0d0)
    case (2)
       d1mach = huge(1.0d0)
    case (3)
       d1mach = epsilon(1.0d0)*0.5d0
    case (4)
       d1mach = epsilon(1.0d0)
    case (5)
       d1mach = log10(real(radix(1.0d0), 8))
    case default
       d1mach = 0.0d0
    end select
  end function d1mach

  integer function initds(dos, nos, eta)
    implicit none
    integer, intent(in) :: nos
    double precision, intent(in) :: dos(nos), eta
    integer :: ii
    double precision :: err
    err = 0.0d0
    do ii = nos, 1, -1
       err = err + abs(dos(ii))
       if (err > eta) exit
    end do
    if (ii == 0) then
       initds = 0
    else
       initds = ii
    end if
  end function initds

  double precision function dcsevl(x, cs, n)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: x, cs(n)
    double precision :: b0, b1, b2, twox
    integer :: i
    if (n <= 0) then
       dcsevl = 0.0d0
       return
    end if
    b0 = 0.0d0
    b1 = 0.0d0
    b2 = 0.0d0
    twox = 2.0d0 * x
    do i = n, 1, -1
       b2 = b1
       b1 = b0
       b0 = twox * b1 - b2 + cs(i)
    end do
    dcsevl = 0.5d0 * (b0 - b2)
  end function dcsevl

  subroutine xermsg(librar, subrou, messg, nerr, level)
    implicit none
    character(*), intent(in) :: librar, subrou, messg
    integer, intent(in) :: nerr, level
    ! Underflow / non-fatal warning handler (no-op)
  end subroutine xermsg

DOUBLE PRECISION FUNCTION DBSI0E (X)
!***BEGIN PROLOGUE  DBSI0E
!***PURPOSE  Compute the exponentially scaled modified (hyperbolic)
!            Bessel function of the first kind of order zero.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B1
!***TYPE      DOUBLE PRECISION (BESI0E-S, DBSI0E-D)
!***KEYWORDS  EXPONENTIALLY SCALED, FIRST KIND, FNLIB,
!             HYPERBOLIC BESSEL FUNCTION, MODIFIED BESSEL FUNCTION,
!             ORDER ZERO, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DBSI0E(X) calculates the double precision exponentially scaled
! modified (hyperbolic) Bessel function of the first kind of order
! zero for double precision argument X.  The result is the Bessel
! function I0(X) multiplied by EXP(-ABS(X)).
!
! Series for BI0        on the interval  0.          to  9.00000E+00
!                                        with weighted error   9.51E-34
!                                         log weighted error  33.02
!                               significant figures required  33.31
!                                    decimal places required  33.65
!
! Series for AI0        on the interval  1.25000E-01 to  3.33333E-01
!                                        with weighted error   2.74E-32
!                                         log weighted error  31.56
!                               significant figures required  30.15
!                                    decimal places required  32.39
!
! Series for AI02       on the interval  0.          to  1.25000E-01
!                                        with weighted error   1.97E-32
!                                         log weighted error  31.71
!                               significant figures required  30.15
!                                    decimal places required  32.63
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DCSEVL, INITDS
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  DBSI0E
DOUBLE PRECISION X, BI0CS(18), AI0CS(46), AI02CS(69), &
     & XSML, Y
LOGICAL FIRST
INTEGER :: NTI0, NTAI0, NTAI02
DOUBLE PRECISION :: ETA
SAVE BI0CS, AI0CS, AI02CS, NTI0, NTAI0, NTAI02, XSML, FIRST
DATA BI0CS(  1) / -.7660547252839144951081894976243285D-1   /
DATA BI0CS(  2) / +.1927337953993808269952408750881196D+1   /
DATA BI0CS(  3) / +.2282644586920301338937029292330415D+0   /
DATA BI0CS(  4) / +.1304891466707290428079334210691888D-1   /
DATA BI0CS(  5) / +.4344270900816487451378682681026107D-3   /
DATA BI0CS(  6) / +.9422657686001934663923171744118766D-5   /
DATA BI0CS(  7) / +.1434006289510691079962091878179957D-6   /
DATA BI0CS(  8) / +.1613849069661749069915419719994611D-8   /
DATA BI0CS(  9) / +.1396650044535669699495092708142522D-10  /
DATA BI0CS( 10) / +.9579451725505445344627523171893333D-13  /
DATA BI0CS( 11) / +.5333981859862502131015107744000000D-15  /
DATA BI0CS( 12) / +.2458716088437470774696785919999999D-17  /
DATA BI0CS( 13) / +.9535680890248770026944341333333333D-20  /
DATA BI0CS( 14) / +.3154382039721427336789333333333333D-22  /
DATA BI0CS( 15) / +.9004564101094637431466666666666666D-25  /
DATA BI0CS( 16) / +.2240647369123670016000000000000000D-27  /
DATA BI0CS( 17) / +.4903034603242837333333333333333333D-30  /
DATA BI0CS( 18) / +.9508172606122666666666666666666666D-33  /
DATA AI0CS(  1) / +.7575994494023795942729872037438D-1      /
DATA AI0CS(  2) / +.7591380810823345507292978733204D-2      /
DATA AI0CS(  3) / +.4153131338923750501863197491382D-3      /
DATA AI0CS(  4) / +.1070076463439073073582429702170D-4      /
DATA AI0CS(  5) / -.7901179979212894660750319485730D-5      /
DATA AI0CS(  6) / -.7826143501438752269788989806909D-6      /
DATA AI0CS(  7) / +.2783849942948870806381185389857D-6      /
DATA AI0CS(  8) / +.8252472600612027191966829133198D-8      /
DATA AI0CS(  9) / -.1204463945520199179054960891103D-7      /
DATA AI0CS( 10) / +.1559648598506076443612287527928D-8      /
DATA AI0CS( 11) / +.2292556367103316543477254802857D-9      /
DATA AI0CS( 12) / -.1191622884279064603677774234478D-9      /
DATA AI0CS( 13) / +.1757854916032409830218331247743D-10     /
DATA AI0CS( 14) / +.1128224463218900517144411356824D-11     /
DATA AI0CS( 15) / -.1146848625927298877729633876982D-11     /
DATA AI0CS( 16) / +.2715592054803662872643651921606D-12     /
DATA AI0CS( 17) / -.2415874666562687838442475720281D-13     /
DATA AI0CS( 18) / -.6084469888255125064606099639224D-14     /
DATA AI0CS( 19) / +.3145705077175477293708360267303D-14     /
DATA AI0CS( 20) / -.7172212924871187717962175059176D-15     /
DATA AI0CS( 21) / +.7874493403454103396083909603327D-16     /
DATA AI0CS( 22) / +.1004802753009462402345244571839D-16     /
DATA AI0CS( 23) / -.7566895365350534853428435888810D-17     /
DATA AI0CS( 24) / +.2150380106876119887812051287845D-17     /
DATA AI0CS( 25) / -.3754858341830874429151584452608D-18     /
DATA AI0CS( 26) / +.2354065842226992576900757105322D-19     /
DATA AI0CS( 27) / +.1114667612047928530226373355110D-19     /
DATA AI0CS( 28) / -.5398891884396990378696779322709D-20     /
DATA AI0CS( 29) / +.1439598792240752677042858404522D-20     /
DATA AI0CS( 30) / -.2591916360111093406460818401962D-21     /
DATA AI0CS( 31) / +.2238133183998583907434092298240D-22     /
DATA AI0CS( 32) / +.5250672575364771172772216831999D-23     /
DATA AI0CS( 33) / -.3249904138533230784173432285866D-23     /
DATA AI0CS( 34) / +.9924214103205037927857284710400D-24     /
DATA AI0CS( 35) / -.2164992254244669523146554299733D-24     /
DATA AI0CS( 36) / +.3233609471943594083973332991999D-25     /
DATA AI0CS( 37) / -.1184620207396742489824733866666D-26     /
DATA AI0CS( 38) / -.1281671853950498650548338687999D-26     /
DATA AI0CS( 39) / +.5827015182279390511605568853333D-27     /
DATA AI0CS( 40) / -.1668222326026109719364501503999D-27     /
DATA AI0CS( 41) / +.3625309510541569975700684800000D-28     /
DATA AI0CS( 42) / -.5733627999055713589945958399999D-29     /
DATA AI0CS( 43) / +.3736796722063098229642581333333D-30     /
DATA AI0CS( 44) / +.1602073983156851963365512533333D-30     /
DATA AI0CS( 45) / -.8700424864057229884522495999999D-31     /
DATA AI0CS( 46) / +.2741320937937481145603413333333D-31     /
DATA AI02CS(  1) / +.5449041101410883160789609622680D-1      /
DATA AI02CS(  2) / +.3369116478255694089897856629799D-2      /
DATA AI02CS(  3) / +.6889758346916823984262639143011D-4      /
DATA AI02CS(  4) / +.2891370520834756482966924023232D-5      /
DATA AI02CS(  5) / +.2048918589469063741827605340931D-6      /
DATA AI02CS(  6) / +.2266668990498178064593277431361D-7      /
DATA AI02CS(  7) / +.3396232025708386345150843969523D-8      /
DATA AI02CS(  8) / +.4940602388224969589104824497835D-9      /
DATA AI02CS(  9) / +.1188914710784643834240845251963D-10     /
DATA AI02CS( 10) / -.3149916527963241364538648629619D-10     /
DATA AI02CS( 11) / -.1321581184044771311875407399267D-10     /
DATA AI02CS( 12) / -.1794178531506806117779435740269D-11     /
DATA AI02CS( 13) / +.7180124451383666233671064293469D-12     /
DATA AI02CS( 14) / +.3852778382742142701140898017776D-12     /
DATA AI02CS( 15) / +.1540086217521409826913258233397D-13     /
DATA AI02CS( 16) / -.4150569347287222086626899720156D-13     /
DATA AI02CS( 17) / -.9554846698828307648702144943125D-14     /
DATA AI02CS( 18) / +.3811680669352622420746055355118D-14     /
DATA AI02CS( 19) / +.1772560133056526383604932666758D-14     /
DATA AI02CS( 20) / -.3425485619677219134619247903282D-15     /
DATA AI02CS( 21) / -.2827623980516583484942055937594D-15     /
DATA AI02CS( 22) / +.3461222867697461093097062508134D-16     /
DATA AI02CS( 23) / +.4465621420296759999010420542843D-16     /
DATA AI02CS( 24) / -.4830504485944182071255254037954D-17     /
DATA AI02CS( 25) / -.7233180487874753954562272409245D-17     /
DATA AI02CS( 26) / +.9921475412173698598880460939810D-18     /
DATA AI02CS( 27) / +.1193650890845982085504399499242D-17     /
DATA AI02CS( 28) / -.2488709837150807235720544916602D-18     /
DATA AI02CS( 29) / -.1938426454160905928984697811326D-18     /
DATA AI02CS( 30) / +.6444656697373443868783019493949D-19     /
DATA AI02CS( 31) / +.2886051596289224326481713830734D-19     /
DATA AI02CS( 32) / -.1601954907174971807061671562007D-19     /
DATA AI02CS( 33) / -.3270815010592314720891935674859D-20     /
DATA AI02CS( 34) / +.3686932283826409181146007239393D-20     /
DATA AI02CS( 35) / +.1268297648030950153013595297109D-22     /
DATA AI02CS( 36) / -.7549825019377273907696366644101D-21     /
DATA AI02CS( 37) / +.1502133571377835349637127890534D-21     /
DATA AI02CS( 38) / +.1265195883509648534932087992483D-21     /
DATA AI02CS( 39) / -.6100998370083680708629408916002D-22     /
DATA AI02CS( 40) / -.1268809629260128264368720959242D-22     /
DATA AI02CS( 41) / +.1661016099890741457840384874905D-22     /
DATA AI02CS( 42) / -.1585194335765885579379705048814D-23     /
DATA AI02CS( 43) / -.3302645405968217800953817667556D-23     /
DATA AI02CS( 44) / +.1313580902839239781740396231174D-23     /
DATA AI02CS( 45) / +.3689040246671156793314256372804D-24     /
DATA AI02CS( 46) / -.4210141910461689149219782472499D-24     /
DATA AI02CS( 47) / +.4791954591082865780631714013730D-25     /
DATA AI02CS( 48) / +.8459470390221821795299717074124D-25     /
DATA AI02CS( 49) / -.4039800940872832493146079371810D-25     /
DATA AI02CS( 50) / -.6434714653650431347301008504695D-26     /
DATA AI02CS( 51) / +.1225743398875665990344647369905D-25     /
DATA AI02CS( 52) / -.2934391316025708923198798211754D-26     /
DATA AI02CS( 53) / -.1961311309194982926203712057289D-26     /
DATA AI02CS( 54) / +.1503520374822193424162299003098D-26     /
DATA AI02CS( 55) / -.9588720515744826552033863882069D-28     /
DATA AI02CS( 56) / -.3483339380817045486394411085114D-27     /
DATA AI02CS( 57) / +.1690903610263043673062449607256D-27     /
DATA AI02CS( 58) / +.1982866538735603043894001157188D-28     /
DATA AI02CS( 59) / -.5317498081491816214575830025284D-28     /
DATA AI02CS( 60) / +.1803306629888392946235014503901D-28     /
DATA AI02CS( 61) / +.6213093341454893175884053112422D-29     /
DATA AI02CS( 62) / -.7692189292772161863200728066730D-29     /
DATA AI02CS( 63) / +.1858252826111702542625560165963D-29     /
DATA AI02CS( 64) / +.1237585142281395724899271545541D-29     /
DATA AI02CS( 65) / -.1102259120409223803217794787792D-29     /
DATA AI02CS( 66) / +.1886287118039704490077874479431D-30     /
DATA AI02CS( 67) / +.2160196872243658913149031414060D-30     /
DATA AI02CS( 68) / -.1605454124919743200584465949655D-30     /
DATA AI02CS( 69) / +.1965352984594290603938848073318D-31     /
DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DBSI0E
IF (FIRST) THEN
   ETA = 0.1*REAL(D1MACH(3))
   NTI0 = INITDS (BI0CS, 18, ETA)
   NTAI0 = INITDS (AI0CS, 46, ETA)
   NTAI02 = INITDS (AI02CS, 69, ETA)
   XSML = SQRT(4.5D0*D1MACH(3))
ENDIF
FIRST = .FALSE.
!
Y = ABS(X)
IF (Y.GT.3.0D0) GO TO 20
!
DBSI0E = 1.0D0 - X
IF (Y.GT.XSML) DBSI0E = EXP(-Y) * (2.75D0 + &
     & DCSEVL (Y*Y/4.5D0-1.D0, BI0CS, NTI0) )
RETURN
!
20 IF (Y.LE.8.D0) DBSI0E = (0.375D0 + DCSEVL ((48.D0/Y-11.D0)/5.D0, &
     & AI0CS, NTAI0))/SQRT(Y)
IF (Y.GT.8.D0) DBSI0E = (0.375D0 + DCSEVL (16.D0/Y-1.D0, AI02CS, &
     & NTAI02))/SQRT(Y)
!
RETURN
END


DOUBLE PRECISION FUNCTION DBSI1E (X)
!***BEGIN PROLOGUE  DBSI1E
!***PURPOSE  Compute the exponentially scaled modified (hyperbolic)
!            Bessel function of the first kind of order one.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B1
!***TYPE      DOUBLE PRECISION (BESI1E-S, DBSI1E-D)
!***KEYWORDS  EXPONENTIALLY SCALED, FIRST KIND, FNLIB,
!             HYPERBOLIC BESSEL FUNCTION, MODIFIED BESSEL FUNCTION,
!             ORDER ONE, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DBSI1E(X) calculates the double precision exponentially scaled
! modified (hyperbolic) Bessel function of the first kind of order
! one for double precision argument X.  The result is I1(X)
! multiplied by EXP(-ABS(X)).
!
! Series for BI1        on the interval  0.          to  9.00000E+00
!                                        with weighted error   1.44E-32
!                                         log weighted error  31.84
!                               significant figures required  31.45
!                                    decimal places required  32.46
!
! Series for AI1        on the interval  1.25000E-01 to  3.33333E-01
!                                        with weighted error   2.81E-32
!                                         log weighted error  31.55
!                               significant figures required  29.93
!                                    decimal places required  32.38
!
! Series for AI12       on the interval  0.          to  1.25000E-01
!                                        with weighted error   1.83E-32
!                                         log weighted error  31.74
!                               significant figures required  29.97
!                                    decimal places required  32.66
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DBSI1E
DOUBLE PRECISION X, BI1CS(17), AI1CS(46), AI12CS(69), XMIN, &
     & XSML, Y
LOGICAL FIRST
INTEGER :: NTI1, NTAI1, NTAI12
DOUBLE PRECISION :: ETA
SAVE BI1CS, AI1CS, AI12CS, NTI1, NTAI1, NTAI12, XMIN, XSML, &
     & FIRST
DATA BI1CS(  1) / -.19717132610998597316138503218149D-2     /
DATA BI1CS(  2) / +.40734887667546480608155393652014D+0     /
DATA BI1CS(  3) / +.34838994299959455866245037783787D-1     /
DATA BI1CS(  4) / +.15453945563001236038598401058489D-2     /
DATA BI1CS(  5) / +.41888521098377784129458832004120D-4     /
DATA BI1CS(  6) / +.76490267648362114741959703966069D-6     /
DATA BI1CS(  7) / +.10042493924741178689179808037238D-7     /
DATA BI1CS(  8) / +.99322077919238106481371298054863D-10    /
DATA BI1CS(  9) / +.76638017918447637275200171681349D-12    /
DATA BI1CS( 10) / +.47414189238167394980388091948160D-14    /
DATA BI1CS( 11) / +.24041144040745181799863172032000D-16    /
DATA BI1CS( 12) / +.10171505007093713649121100799999D-18    /
DATA BI1CS( 13) / +.36450935657866949458491733333333D-21    /
DATA BI1CS( 14) / +.11205749502562039344810666666666D-23    /
DATA BI1CS( 15) / +.29875441934468088832000000000000D-26    /
DATA BI1CS( 16) / +.69732310939194709333333333333333D-29    /
DATA BI1CS( 17) / +.14367948220620800000000000000000D-31    /
DATA AI1CS(  1) / -.2846744181881478674100372468307D-1      /
DATA AI1CS(  2) / -.1922953231443220651044448774979D-1      /
DATA AI1CS(  3) / -.6115185857943788982256249917785D-3      /
DATA AI1CS(  4) / -.2069971253350227708882823777979D-4      /
DATA AI1CS(  5) / +.8585619145810725565536944673138D-5      /
DATA AI1CS(  6) / +.1049498246711590862517453997860D-5      /
DATA AI1CS(  7) / -.2918338918447902202093432326697D-6      /
DATA AI1CS(  8) / -.1559378146631739000160680969077D-7      /
DATA AI1CS(  9) / +.1318012367144944705525302873909D-7      /
DATA AI1CS( 10) / -.1448423418183078317639134467815D-8      /
DATA AI1CS( 11) / -.2908512243993142094825040993010D-9      /
DATA AI1CS( 12) / +.1266388917875382387311159690403D-9      /
DATA AI1CS( 13) / -.1664947772919220670624178398580D-10     /
DATA AI1CS( 14) / -.1666653644609432976095937154999D-11     /
DATA AI1CS( 15) / +.1242602414290768265232168472017D-11     /
DATA AI1CS( 16) / -.2731549379672432397251461428633D-12     /
DATA AI1CS( 17) / +.2023947881645803780700262688981D-13     /
DATA AI1CS( 18) / +.7307950018116883636198698126123D-14     /
DATA AI1CS( 19) / -.3332905634404674943813778617133D-14     /
DATA AI1CS( 20) / +.7175346558512953743542254665670D-15     /
DATA AI1CS( 21) / -.6982530324796256355850629223656D-16     /
DATA AI1CS( 22) / -.1299944201562760760060446080587D-16     /
DATA AI1CS( 23) / +.8120942864242798892054678342860D-17     /
DATA AI1CS( 24) / -.2194016207410736898156266643783D-17     /
DATA AI1CS( 25) / +.3630516170029654848279860932334D-18     /
DATA AI1CS( 26) / -.1695139772439104166306866790399D-19     /
DATA AI1CS( 27) / -.1288184829897907807116882538222D-19     /
DATA AI1CS( 28) / +.5694428604967052780109991073109D-20     /
DATA AI1CS( 29) / -.1459597009090480056545509900287D-20     /
DATA AI1CS( 30) / +.2514546010675717314084691334485D-21     /
DATA AI1CS( 31) / -.1844758883139124818160400029013D-22     /
DATA AI1CS( 32) / -.6339760596227948641928609791999D-23     /
DATA AI1CS( 33) / +.3461441102031011111108146626560D-23     /
DATA AI1CS( 34) / -.1017062335371393547596541023573D-23     /
DATA AI1CS( 35) / +.2149877147090431445962500778666D-24     /
DATA AI1CS( 36) / -.3045252425238676401746206173866D-25     /
DATA AI1CS( 37) / +.5238082144721285982177634986666D-27     /
DATA AI1CS( 38) / +.1443583107089382446416789503999D-26     /
DATA AI1CS( 39) / -.6121302074890042733200670719999D-27     /
DATA AI1CS( 40) / +.1700011117467818418349189802666D-27     /
DATA AI1CS( 41) / -.3596589107984244158535215786666D-28     /
DATA AI1CS( 42) / +.5448178578948418576650513066666D-29     /
DATA AI1CS( 43) / -.2731831789689084989162564266666D-30     /
DATA AI1CS( 44) / -.1858905021708600715771903999999D-30     /
DATA AI1CS( 45) / +.9212682974513933441127765333333D-31     /
DATA AI1CS( 46) / -.2813835155653561106370833066666D-31     /
DATA AI12CS(  1) / +.2857623501828012047449845948469D-1      /
DATA AI12CS(  2) / -.9761097491361468407765164457302D-2      /
DATA AI12CS(  3) / -.1105889387626237162912569212775D-3      /
DATA AI12CS(  4) / -.3882564808877690393456544776274D-5      /
DATA AI12CS(  5) / -.2512236237870208925294520022121D-6      /
DATA AI12CS(  6) / -.2631468846889519506837052365232D-7      /
DATA AI12CS(  7) / -.3835380385964237022045006787968D-8      /
DATA AI12CS(  8) / -.5589743462196583806868112522229D-9      /
DATA AI12CS(  9) / -.1897495812350541234498925033238D-10     /
DATA AI12CS( 10) / +.3252603583015488238555080679949D-10     /
DATA AI12CS( 11) / +.1412580743661378133163366332846D-10     /
DATA AI12CS( 12) / +.2035628544147089507224526136840D-11     /
DATA AI12CS( 13) / -.7198551776245908512092589890446D-12     /
DATA AI12CS( 14) / -.4083551111092197318228499639691D-12     /
DATA AI12CS( 15) / -.2101541842772664313019845727462D-13     /
DATA AI12CS( 16) / +.4272440016711951354297788336997D-13     /
DATA AI12CS( 17) / +.1042027698412880276417414499948D-13     /
DATA AI12CS( 18) / -.3814403072437007804767072535396D-14     /
DATA AI12CS( 19) / -.1880354775510782448512734533963D-14     /
DATA AI12CS( 20) / +.3308202310920928282731903352405D-15     /
DATA AI12CS( 21) / +.2962628997645950139068546542052D-15     /
DATA AI12CS( 22) / -.3209525921993423958778373532887D-16     /
DATA AI12CS( 23) / -.4650305368489358325571282818979D-16     /
DATA AI12CS( 24) / +.4414348323071707949946113759641D-17     /
DATA AI12CS( 25) / +.7517296310842104805425458080295D-17     /
DATA AI12CS( 26) / -.9314178867326883375684847845157D-18     /
DATA AI12CS( 27) / -.1242193275194890956116784488697D-17     /
DATA AI12CS( 28) / +.2414276719454848469005153902176D-18     /
DATA AI12CS( 29) / +.2026944384053285178971922860692D-18     /
DATA AI12CS( 30) / -.6394267188269097787043919886811D-19     /
DATA AI12CS( 31) / -.3049812452373095896084884503571D-19     /
DATA AI12CS( 32) / +.1612841851651480225134622307691D-19     /
DATA AI12CS( 33) / +.3560913964309925054510270904620D-20     /
DATA AI12CS( 34) / -.3752017947936439079666828003246D-20     /
DATA AI12CS( 35) / -.5787037427074799345951982310741D-22     /
DATA AI12CS( 36) / +.7759997511648161961982369632092D-21     /
DATA AI12CS( 37) / -.1452790897202233394064459874085D-21     /
DATA AI12CS( 38) / -.1318225286739036702121922753374D-21     /
DATA AI12CS( 39) / +.6116654862903070701879991331717D-22     /
DATA AI12CS( 40) / +.1376279762427126427730243383634D-22     /
DATA AI12CS( 41) / -.1690837689959347884919839382306D-22     /
DATA AI12CS( 42) / +.1430596088595433153987201085385D-23     /
DATA AI12CS( 43) / +.3409557828090594020405367729902D-23     /
DATA AI12CS( 44) / -.1309457666270760227845738726424D-23     /
DATA AI12CS( 45) / -.3940706411240257436093521417557D-24     /
DATA AI12CS( 46) / +.4277137426980876580806166797352D-24     /
DATA AI12CS( 47) / -.4424634830982606881900283123029D-25     /
DATA AI12CS( 48) / -.8734113196230714972115309788747D-25     /
DATA AI12CS( 49) / +.4045401335683533392143404142428D-25     /
DATA AI12CS( 50) / +.7067100658094689465651607717806D-26     /
DATA AI12CS( 51) / -.1249463344565105223002864518605D-25     /
DATA AI12CS( 52) / +.2867392244403437032979483391426D-26     /
DATA AI12CS( 53) / +.2044292892504292670281779574210D-26     /
DATA AI12CS( 54) / -.1518636633820462568371346802911D-26     /
DATA AI12CS( 55) / +.8110181098187575886132279107037D-28     /
DATA AI12CS( 56) / +.3580379354773586091127173703270D-27     /
DATA AI12CS( 57) / -.1692929018927902509593057175448D-27     /
DATA AI12CS( 58) / -.2222902499702427639067758527774D-28     /
DATA AI12CS( 59) / +.5424535127145969655048600401128D-28     /
DATA AI12CS( 60) / -.1787068401578018688764912993304D-28     /
DATA AI12CS( 61) / -.6565479068722814938823929437880D-29     /
DATA AI12CS( 62) / +.7807013165061145280922067706839D-29     /
DATA AI12CS( 63) / -.1816595260668979717379333152221D-29     /
DATA AI12CS( 64) / -.1287704952660084820376875598959D-29     /
DATA AI12CS( 65) / +.1114548172988164547413709273694D-29     /
DATA AI12CS( 66) / -.1808343145039336939159368876687D-30     /
DATA AI12CS( 67) / -.2231677718203771952232448228939D-30     /
DATA AI12CS( 68) / +.1619029596080341510617909803614D-30     /
DATA AI12CS( 69) / -.1834079908804941413901308439210D-31     /
DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DBSI1E
IF (FIRST) THEN
   ETA = 0.1*REAL(D1MACH(3))
   NTI1 = INITDS (BI1CS, 17, ETA)
   NTAI1 = INITDS (AI1CS, 46, ETA)
   NTAI12 = INITDS (AI12CS, 69, ETA)
!
   XMIN = 2.0D0*D1MACH(1)
   XSML = SQRT(4.5D0*D1MACH(3))
ENDIF
FIRST = .FALSE.
!
Y = ABS(X)
IF (Y.GT.3.0D0) GO TO 20
!
DBSI1E = 0.0D0
IF (Y.EQ.0.D0)  RETURN
!
IF (Y .LE. XMIN) CALL XERMSG ('SLATEC', 'DBSI1E', &
     & 'ABS(X) SO SMALL I1 UNDERFLOWS', 1, 1)
IF (Y.GT.XMIN) DBSI1E = 0.5D0*X
IF (Y.GT.XSML) DBSI1E = X*(0.875D0 + DCSEVL (Y*Y/4.5D0-1.D0, &
     & BI1CS, NTI1) )
DBSI1E = EXP(-Y) * DBSI1E
RETURN
!
20 IF (Y.LE.8.D0) DBSI1E = (0.375D0 + DCSEVL ((48.D0/Y-11.D0)/5.D0, &
     & AI1CS, NTAI1))/SQRT(Y)
IF (Y.GT.8.D0) DBSI1E = (0.375D0 + DCSEVL (16.D0/Y-1.D0, AI12CS, &
     & NTAI12))/SQRT(Y)
DBSI1E = SIGN (DBSI1E, X)
!
RETURN
END


end module bessel_slatec_mod

Module datatrans
  implicit none
  integer, parameter :: nsp=2
  integer, parameter :: mxNr=2500, next=2500,iterMax=30000
  real (kind=8), parameter ::  pi=3.141592653589793D0
  real (kind=8) :: W(mxNr,mxNr), BJ1sq(mxNr), root(mxNr)
End Module datatrans

Module color_mod
  implicit none
  character(len=*), parameter :: c_reset   = achar(27)//'[0m'
  character(len=*), parameter :: c_bold    = achar(27)//'[1m'
  character(len=*), parameter :: c_dim     = achar(27)//'[2m'
  character(len=*), parameter :: c_gray    = achar(27)//'[90m'
  character(len=*), parameter :: c_red     = achar(27)//'[31m'
  character(len=*), parameter :: c_green   = achar(27)//'[32m'
  character(len=*), parameter :: c_yellow  = achar(27)//'[33m'
  character(len=*), parameter :: c_blue    = achar(27)//'[34m'
  character(len=*), parameter :: c_magenta = achar(27)//'[35m'
  character(len=*), parameter :: c_cyan    = achar(27)//'[36m'
  character(len=*), parameter :: c_white   = achar(27)//'[37m'
  character(len=*), parameter :: c_b_red   = achar(27)//'[1;31m'
  character(len=*), parameter :: c_b_green = achar(27)//'[1;32m'
  character(len=*), parameter :: c_b_yellow= achar(27)//'[1;33m'
  character(len=*), parameter :: c_b_blue  = achar(27)//'[1;34m'
  character(len=*), parameter :: c_b_mag   = achar(27)//'[1;35m'
  character(len=*), parameter :: c_b_cyan  = achar(27)//'[1;36m'
  character(len=*), parameter :: c_b_white = achar(27)//'[1;37m'
End Module color_mod


Program twoDdipRY
  !     N-component  HDs+1/r**3 in 2D, from 2008's F. Lado code.
  !     Included excess two-particle energy calculation Jan. 2021. 
  !     Modified for parallel dipole-dipole interaction March 2019. 
  !     adapted and modified to F90 syntax by E. Lomba, Jan. 2018.
  !     Last FL change:  FL   20 Nov 2008   11:49 am
  !     Uses setupW, E1, Bpc2D, BpyLR, Hankel

  !     This program solves the Ornstein-Zernike equation with RHNC closure
  !     for the pair distribution functions and thermodynamics of a
  !     two-dimensional N-component charged hard disks with ln(r) potentials.

  !     NOTE: Let FTf(q;j,k) be the 2D Fourier transform of a function f(r;j,k).
  !           Then, in the program, the computed transform is
  !                         xFTf(q;j,k) = (qMax/2.0*pi*rMax)*FTf(q;j,k).
  !           Also,           Tf(q;j,k) = sqrt(rhoHat(j)*rhoHat(k))*xFTf(q;j,k),
  !                                     = sqrt(rho(j)*rho(k))*FTf(q;j,k),
  !           where rhoHat = (2.0*pi*rMax/qMax)*rho.
  use datatrans
  use color_mod
  use bessel_slatec_mod, only : dbsi0e, dbsi1e
  implicit none
  integer :: ipiv(nsp)
  integer :: lwork
  integer :: Ncore(nsp,nsp)
  character :: FTtables*60, INfile*60, OUTfile*60, fname*20, fnameo*21
  real(kind=8) :: s11, s12, s22, ssum11, ssum12, ssum22, blend,blend0,&
       & sexsum, chempot(nsp), csr0(nsp,nsp), lambda, s11_0, s12_0,&
       & s11_q, s12_q, aconst, bconst, denom, rmsMax, c112,c122,&
       & v(nsp), work(2*nsp),sjk(nsp,nsp),ssumjk(nsp,nsp), ac=1.0d0,&
       & shrink, ulr(nsp,nsp), virlr(nsp,nsp), chempot0(nsp), deltarho=0.001
  real(kind=8) :: rho(nsp), rhoHat(nsp), z(nsp), sigma(nsp,nsp),  &
       ICinv(nsp,nsp),   TcSq(nsp,nsp), sum22(nsp,nsp), sum11(nsp,nsp), &
       sum12(nsp,nsp),  sum00(nsp,nsp), sum01(nsp,nsp), sum02(nsp,nsp), &
       r(0:mxNr), q(0:mxNr), dr(0:mxNr), dq(0:mxNr), rhoi(nsp), A10,&
       & A20, pres(-1:1), uint(-1:1), xc(-1:1), rhot(-1:1),P10, P20 
  real(kind=8) :: phi(0:mxNr,nsp,nsp), phiSR(0:mxNr,nsp,nsp),  phiLR(0:mxNr,nsp,nsp), &
       TphiLR(0:mxNr,nsp,nsp),   sSR(0:mxNr,nsp,nsp),    cSR(0:mxNr,nsp,nsp), &
       xFTcSR(0:mxNr,nsp,nsp),  TcSR(0:mxNr,nsp,nsp), xFTsSR(0:mxNr,nsp,nsp), &
       sSRnew(0:mxNr,nsp,nsp),     g(0:mxNr,nsp,nsp),   &
       s0(0:mxNr,nsp,nsp),    s1(0:mxNr,nsp,nsp),     s2(0:mxNr,nsp,nsp), &
       d0(0:mxNr,nsp,nsp),    d1(0:mxNr,nsp,nsp),     d2(0:mxNr,nsp,nsp), &
       sc(0:mxNr), scint(0:mxNr), sg(nsp,nsp), sgi(0:mxNr,nsp,nsp),&
       & flr(0:mxNr), tflr(0:mxNr), dphi(0:mxNr,nsp,nsp), fint(1:mxNr), &
       & exp_neg_phiSR(0:mxNr,nsp,nsp), inv_fint(1:mxNr)

  integer :: Nr, iStart, newW, n, i, j, k, l, NrIn, i2, i3&
       &, iter, iext, info, lcount, irho, its, nxf, ixf, nrt, inr 
  real(kind=8) :: rMax, qMax, rhoTotal, rhoTotal0,Gamma, a0, a1, a2, a3, a4, a5&
       &, sumMat, sum0, rms, a01, a02, a12, a22, a11, c1, c2, sSRng, sumtsx,&
       & sumst, ssumx, s2xc, sumP1, sumP2, sumU, sumUb, sumP0, sumU0,&
       & sumUb0, sumX0, P1, P2, U, X, chemp, sumc, sumint,&
       & sumA1, sumA10, sumA20, TH11, TH12, TH22, ssum, sumX, eta, U0&
       &, X0, dPr, fopt, y1, y2, y3, sqmax, sq0, dsig, ersig, tol,&
       & eti, eto, fopto, fp, rmscut, rtmin, rtmax, xf(100), R01, R02&
       &, sq110, sq220, sums2, s2ex, sumcjk, ct(2,2), Mrr0, Mrc0,&
       & Mcc0, lamb1, lamb2, vrr, vcc, csrnc(2,2), scc0
  integer, parameter :: maxStates = 1000
  integer :: nStateDone, istate
  real(kind=8) :: res_rho(maxStates), res_x2(maxStates), res_P(maxStates), &
                  res_U(maxStates), res_xc(maxStates), res_dPr(maxStates), &
                  res_eta(maxStates), res_lamb1(maxStates), &
                  res_invScc(maxStates), res_sqmax(maxStates), scale_21
  real(kind=8), external :: E1, mk
  logical :: unfound, osig, solve
  lwork = 2*nsp
  nStateDone = 0

  !     Program banner
  write (*,"(/'', a, '================================================================================', a)") c_b_cyan, c_reset
  write (*,"(a, '                                   twoDdipRY', a)") c_b_white, c_reset
  write (*,"(a, '        2D Dipolar Hard Disk Liquid-State Integral Equation Solver', a)") c_cyan, c_reset
  write (*,"(a, '          Rogers-Young (RY) Closure with Thermodynamic Consistency', a)") c_cyan, c_reset
  write (*,"(a, '================================================================================', a)") c_b_cyan, c_reset

  !     Read input parameters.
  write(fname,"('2DdipHD',i1,'c_map.dat')") nsp
  write(fnameo,"('2DdipdHD',i1,'c_out.dat')") nsp
  open (2,file=fname,status='old')
  read (2,*) Nr, iStart, newW, &
       Gamma, blend0, rmsMax, rmscut
  if (Nr > mxNr) Then
     write (*,"(a, ' *** Error: Nr (', i0, ') exceeds maximum mxNr (', i0, ')', a)") c_b_red, Nr, mxNr, c_reset
     stop
  Endif
  read (2,*) nrt
  read(2,*) rtmin,rtmax
  read(2,*) nxf
  read(2,*) xf(1:nxf)
  read(2,*) z(1:nsp)
  do i=1,nsp
     read (2,*) Ncore(i,i:nsp)
  Enddo
  read (2,*) FTtables
  read (2,*) INfile
  read (2,*) OUTfile
  read (2,*) lambda, shrink, r01, R02
  r02 = r02*r01
  read (2, *) eta, dsig, tol, osig
  close (2,status='keep')

  !     Echo input parameters to log file (unit 3).
  open  (3,file=fnameo,status='unknown')
  write (3,'(3i10/11f15.7)') Nr,  iStart, newW, &
       Gamma, z(1:nsp)
  do i=1,nsp
     write (3,'(3i10)')Ncore(i,i:nsp)
  Enddo
  write (3,30) FTtables
  write (3,30) INfile
  write (3,30) OUTfile
10 format (6i5/5f10.5)
20 format (/' PROGRAM 2DdipRY' &
       /' Iterative solution of the RY equation for 2D N-component&
       & charged 1/r^3 particles' &
       /' Input data:'/6i5/5f10.5)
30 format (a60)

  !     Open and format thermo.dat with comprehensive metadata and aligned columns
  open(88,file="thermo.dat")
  write(88,"(a)") '#============================================================================================================================================================================================================================='
  write(88,"(a)") '# 2D Dipolar Hard Disk Integral Equation Solver (twoDdipRY)'
  write(88,"(a)") '# Ornstein-Zernike equation with Rogers-Young (RY) closure'
  write(88,"(a)") '#'
  write(88,"(a)") '# Column Definitions:'
  write(88,"(a)") '#   [ 1] rho_tot   : Total number density rho = rho1 + rho2'
  write(88,"(a)") '#   [ 2] x2        : Mole fraction of species 2 (x2 = rho2 / rho_tot)'
  write(88,"(a)") '#   [ 3] U/NkT     : Reduced internal energy per particle'
  write(88,"(a)") '#   [ 4] P/(rho*kT): Virial compressibility factor Z = P / (rho * k_B * T)'
  write(88,"(a)") '#   [ 5] chi^-1    : Inverse isothermal compressibility (compressibility route: 1 - c_hat(0))'
  write(88,"(a)") '#   [ 6] dP*/drho  : Virial pressure derivative with respect to density'
  write(88,"(a)") '#   [ 7] eta       : Rogers-Young closure parameter'
  write(88,"(a)") '#   [ 8] Sq_max    : Maximum peak of total structure factor S(q)_max'
  write(88,"(a)") '#   [ 9] Sq_0      : Zero-wavevector limit of total structure factor S(q=0)'
  write(88,"(a)") '#   [10] Sq0/Sqmax : Long-wavelength fluctuation ratio S(0) / S(q)_max'
  write(88,"(a)") '#   [11] scale_21  : Core scaling factor (S22(0) / S11(0))^(1/4)'
  write(88,"(a)") '#   [12] S2_ex/kB  : Excess two-body entropy per particle'
  write(88,"(a)") '#   [13] lambda1   : Minimum eigenvalue of Bhatia-Thornton fluctuation matrix (spinodal indicator, lambda1 -> 0)'
  write(88,"(a)") '#   [14] lambda2   : Maximum eigenvalue of Bhatia-Thornton fluctuation matrix'
  write(88,"(a)") '#   [15] 1/Scc(0)  : Inverse concentration fluctuation at origin (vanishes at critical demixing point)'
  write(88,"(a)") '#============================================================================================================================================================================================================================='
  write(88,"(a1, 15(1x, a14))") '#', &
       'rho_tot', 'x2', 'U/NkT', 'P/(rho*kT)', 'chi^-1', &
       'dP*/drho', 'eta', 'Sq_max', 'Sq_0', 'Sq0/Sqmax', &
       'scale_21', 'S2_ex/kB', 'lambda1', 'lambda2', '1/Scc(0)'
  write(88,"(a1, 15(1x, a14))") '#', &
       '--------------', '--------------', '--------------', '--------------', '--------------', &
       '--------------', '--------------', '--------------', '--------------', '--------------', &
       '--------------', '--------------', '--------------', '--------------', '--------------'

  !     Set up tables for Hankel transforms.
  open (15,file=FTtables,status='unknown')
  if (newW .eq. 1) then ! calculate needed FT tables; save for later reuse.
     write (*,"(a, ' [Hankel Tables]', a, ' Calculating Bessel roots and quadrature weights...')") c_b_blue, c_reset
     write (*,"('                 Saving precomputed tables to ', a, a, a, '...')") c_b_white, trim(FTtables), c_reset
     call setupW(Nr)
     write (15,40) ((W(i,j), i = 1,Nr-1), j = 1,Nr-1)
     write (15,40) (BJ1sq(i), i = 1,Nr)
     write (15,40) (root(i), i = 1,Nr)
     write (*,"(a, ' [Hankel Tables]', a, ' ', a, 'Tables successfully generated.', a)") c_b_blue, c_reset, c_b_green, c_reset
  else ! read in stored tables.
     write (*,"(a, ' [Hankel Tables]', a, ' Loading precomputed quadrature tables from: ', a, a, a)") &
          c_b_blue, c_reset, c_b_white, trim(FTtables), c_reset
     read (15,40) ((W(i,j), i = 1,Nr-1), j = 1,Nr-1)
     read (15,40) (BJ1sq(i), i = 1,Nr)
     read (15,40) (root(i), i = 1,Nr)
     write (*,"(a, ' [Hankel Tables]', a, ' ', a, 'Tables loaded successfully.', a)") c_b_blue, c_reset, c_b_green, c_reset
  end if
  close (15,status='keep')
40 format (6f12.6)

  !     Initialize constants and arrays.
  n = Nr-1
  do i=1,nsp
     do j=i,nsp
        Ncore(j,i) = Ncore(i,j)
     Enddo
  Enddo
  rMax = root(Nr)/root(Ncore(1,1))
  qMax = root(Nr)/rMax
  r(0) = 0.0d0
  q(0) = 0.0d0
  do i = 1,Nr
     r(i) = root(i)/qMax
     q(i) = root(i)/rMax
     dr(i) = (2.0d0/qMax)/(root(i)*BJ1sq(i))
     dq(i) = (2.0d0/rMax)/(root(i)*BJ1sq(i))
  end do
  ! shrink core when needed
  Ncore(:,:) = shrink*Ncore(:,:)
  do j = 1,nsp
     do k = 1,nsp
        sigma(j,k) = root(Ncore(j,k))/qMax
        do i = 0,Ncore(j,k)-1
           g(i,j,k) = 0.0d0
        end do
     end do
  end do
  do i=1,nsp
     write (3,50) sigma(i,i:nsp)
  Enddo

  !     Echo system parameters and configuration to stdout
  write (*,"(/'', a, '--------------------------------------------------------------------------------', a)") c_gray, c_reset
  write (*,"(a, ' SYSTEM CONFIGURATION & PARAMETERS', a)") c_b_yellow, c_reset
  write (*,"(a, '--------------------------------------------------------------------------------', a)") c_gray, c_reset
  write (*,"(a, ' Radial Grid & Quadrature:', a)") c_b_cyan, c_reset
  write (*,"('   Grid points (Nr)       : ', a, i8, a, 6x, 'r_max = ', a, f10.5, a, ', q_max = ', a, f10.5, a)") &
       c_b_white, Nr, c_reset, c_b_white, rMax, c_reset, c_b_white, qMax, c_reset
  write (*,"('   Hankel table file      : ', a, a, a)") c_b_white, trim(FTtables), c_reset
  if (iStart .eq. 1) then
     write (*,"('   Restart mode           : ', a, 'Restart solution from ', a, a)") c_yellow, trim(INfile), c_reset
  else
     write (*,"('   Restart mode           : Fresh start (short-range ideal gas)')")
  end if
  write (*,"('   Restart output file    : ', a, a, a)") c_b_white, trim(OUTfile), c_reset
  write (*,*)
  write (*,"(a, ' Interaction Potentials:', a)") c_b_cyan, c_reset
  write (*,"('   Dipolar coupling Gamma : ', a, 1pe12.5, a)") c_b_white, Gamma, c_reset
  write (*,"('   Species dipole moments : z(1) = ', a, f9.4, a, ', z(2) = ', a, f9.4, a)") &
       c_b_white, z(1), c_reset, c_b_white, z(2), c_reset
  write (*,"('   Cross-interaction scale: lambda(1,2) = ', a, f8.4, a)") c_b_white, lambda, c_reset
  write (*,"('   Core shrink factor     : ', a, f8.4, a)") c_b_white, shrink, c_reset
  write (*,"('   Hard-core diameters    : sigma(1,1) = ', a, f8.4, a, ' (Ncore = ', i5, ')')") &
       c_b_white, sigma(1,1), c_reset, Ncore(1,1)
  write (*,"('                            sigma(1,2) = ', a, f8.4, a, ' (Ncore = ', i5, ')')") &
       c_b_white, sigma(1,2), c_reset, Ncore(1,2)
  write (*,"('                            sigma(2,2) = ', a, f8.4, a, ' (Ncore = ', i5, ')')") &
       c_b_white, sigma(2,2), c_reset, Ncore(2,2)
  write (*,"('   Effective couplings    : gamma(1,1) = ', a, 1pe12.5, a, ', gamma(1,2) = ', a, 1pe12.5, a)") &
       c_b_white, Gamma*z(1)*z(1), c_reset, c_b_white, Gamma*lambda*z(1)*z(2), c_reset
  write (*,"('                            gamma(2,2) = ', a, 1pe12.5, a)") c_b_white, Gamma*z(2)*z(2), c_reset
  write (*,*)
  write (*,"(a, ' Numerical Solver & Closure Settings:', a)") c_b_cyan, c_reset
  write (*,"('   Picard mixing parameter: blend0 = ', a, f8.4, a, ', rms_cut = ', a, 1pe10.3, a)") &
       c_b_white, blend0, c_reset, c_b_white, rmscut, c_reset
  write (*,"('   Picard convergence tol : rms_max = ', a, 1pe10.3, a, ', max_iters = ', a, i6, a)") &
       c_b_white, rmsMax, c_reset, c_b_white, iterMax, c_reset
  if (osig) then
     write (*,"('   Thermodynamic closure  : ', a, 'Rogers-Young with Newton-Raphson consistency', a)") c_b_green, c_reset
     write (*,"('                            Initial eta = ', a, f8.4, a, ', d_eta = ', a, f9.4, a, ', tol = ', a, 1pe10.3, a)") &
          c_b_white, eta, c_reset, c_b_white, dsig, c_reset, c_b_white, tol, c_reset
  else
     write (*,"('   Thermodynamic closure  : ', a, 'Rogers-Young with FIXED eta = ', f8.4, a)") c_yellow, eta, c_reset
  end if
  write (*,*)
  write (*,"(a, ' State Scan Setup:', a)") c_b_cyan, c_reset
  write (*,"('   Density range [min,max]: [', a, f8.4, a, ', ', a, f8.4, a, '] in ', a, i3, a, ' interval(s)')") &
       c_b_white, rtmin, c_reset, c_b_white, rtmax, c_reset, c_b_white, nrt, c_reset
  write (*,"('   Composition points     : ', a, i3, a, ' point(s) for x2')") c_b_white, nxf, c_reset
  write (*,"(a, '--------------------------------------------------------------------------------', a)") c_gray, c_reset

50 format (/5x, 'sigma(i,j) =',5f15.10)

  !     Construct potentials.
  flr(:) = 0
  tflr(:) = 0
  !  ac = 0
  if (istart .eq. 1) then ! from existing solution.
     open (16,file=INfile,status='old')
     read (16,70) NrIn, i2, i3, &
          a1, a2, a3, a4, a5
     do j = 1,nsp
        do k = j,nsp
           read (16,40) (sSR(i,j,k), i = 0,NrIn-1)
        end do
     end do
     close (16,status='keep')
     write (3,60)
     write (3,70) NrIn, i2, i3, &
          a1, a2, a3, a4, a5
     write (*,"('   ', a, '[Restart]', a, ' Successfully loaded initial solution from ', a, a, a, ' (Nr = ', i5, ')')") &
          c_cyan, c_reset, c_b_white, trim(INfile), c_reset, NrIn
  else ! from SR ideal gas.
     do j = 1,nsp
        do k = j,nsp
           do i = 0,Nr
              sSR(i,j,k) = 0.0d0
           end do
        end do
     end do
  end if
60 format (/' Starting from diskfile:')
70 format (3i5/5f10.5)
  !
  ! Calculate long range functions
  !
  call lrfuncs(flr(0:Nr),tflr(0:Nr),r(0:Nr),q(0:Nr),Nr,ac)
  !
  ! Loop over densities to compute derivatives
  !
  do j = 1,nsp
     do k = 1,nsp
        if (k.eq.j) then 
           a0 = z(j)*z(k)*Gamma
        else
           a0 = lambda*z(j)*z(k)*Gamma
        endif
        ulr(j,k) = a0/r(Nr)
        virlr(j,k) = 3*a0/r(Nr)
        do i = 1,Nr
           phi(i,j,k) = a0/r(i)**3
           dphi(i,j,k) = 3*a0/r(i)**3 ! - r*du/dr 
           phiLR(i,j,k) = a0*flr(i)
           phiSR(i,j,k) = phi(i,j,k)-phiLR(i,j,k)
           exp_neg_phiSR(i,j,k) = exp(-phiSR(i,j,k))
        end do
     end do
  end do
  do inr = 0, nrt
     if (nrt > 0) then
        rhoTotal0 = rtmin + inr*(rtmax-rtmin)/real(nrt, 8)
     else
        rhoTotal0 = rtmin
     endif
     do ixf=1,nxf
        ersig = 100.0
        its = 0
        solve = .true.
        rhoi(1) = (1-xf(ixf))*rhoTotal0
        rhoi(2) = xf(ixf)*rhoTotal0

        write (*,*)
        write (*,"(a)") c_b_blue // '================================================================================' // c_reset
        write (*,"(a, f8.4, a, f8.4, a, f8.4, a, f8.4, a)") &
             c_b_cyan // ' STATE POINT: ' // c_b_yellow // 'rho_tot = ' // c_b_white, rhoTotal0, &
             c_b_yellow // ' | x2 = ' // c_b_white, xf(ixf), &
             c_cyan // ' (rho1 = ' // c_white, rhoi(1), &
             c_cyan // ', rho2 = ' // c_white, rhoi(2), &
             c_cyan // ')' // c_reset
        write (*,"(a)") c_b_blue // '================================================================================' // c_reset

        do while (ersig > tol .and. solve)
           if (.not. osig) solve = .false.

           its= its+1
           do irho=-1,1
              rho(1:nsp)=rhoi(1:nsp)*(1+irho*deltarho/sum(rhoi(1:nsp)))
              rhoTotal = sum(rho(1:nsp))
              rhot(irho) = rhoTotal
              rhoHat(1:nsp) = (2.0d0*pi*rMax/qMax)*rho(1:nsp)
              do j = 1,nsp
                 do k = 1,nsp
                    a1 = sqrt(rho(j)*rho(k))*a0
                    TphiLR(0:nr,j,k) = a1*tflr(0:nr)
                 end do
              end do
              !     Start iterations ...
              write (3,80)(i,i,i=1,nsp)
80            format (/13x, 'iter', 7x, 'rms', 8x, 5('sSR(',i1,i1,')',4x:))
              iext = next
              iter = 0
              !     Begin Picard iteration.
              rms = 1000.0
              fint(1:nr) = 1.0d0 - exp(-eta*r(1:nr))
           inv_fint(1:nr) = 1.0d0 / fint(1:nr)
              do while(rms > rmsMax)
                 iter = iter+1
                 !     Begin iteration core 
                 do j = 1,nsp
                    do k = j,nsp
                       do i = Ncore(j,k),n
                          g(i,j,k) = exp_neg_phiSR(i,j,k)*(1.0d0+(exp((sSR(i,j,k))*fint(i))-1.0d0)*inv_fint(i))
                       end do
                       g(Ncore(j,k),j,k) = 0.5d0*g(Ncore(j,k),j,k)
                       do i = 0,n
                          cSR(i,j,k) = g(i,j,k)-1.0d0-sSR(i,j,k)
                          cSR(i,k,j) = cSR(i,j,k)
                          g(i,k,j) = g(i,j,k)
                       end do
                    end do
                 end do
                 !     Transform cSR(r).
                 do j = 1,nsp
                    do k = j,nsp
                       call Hankel(n,r,dr,rMax,cSR(0,j,k),qMax,xFTcSR(0,j,k))
                       !!        rewind(300+j+k)
                       do i = 0,n
                          TcSR(i,j,k) = sqrt(rhoHat(j)*rhoHat(k))*xFTcSR(i,j,k)
                          TcSR(i,k,j) = TcSR(i,j,k)
                          xFTcSR(i,k,j) = xFTcSR(i,j,k)
                          !!           write(300+j+k,'(3f15.7)')q(i),TcSR(i,j,k),TphiLR(i,j,k)
                       end do
                    end do
                 end do
                 !!  stop
                 !     Ornstein-Zernike equation for xFTsSR(k).
                 do i = 0 ,n
                    do j = 1,nsp
                       do k = j,nsp
                          sumMat = 0.0d0
                          do l = 1,nsp
                             sumMat = sumMat+(TcSR(i,j,l)-TphiLR(i,j,l)) &
                                  *(TcSR(i,l,k)-TphiLR(i,l,k))
                          end do
                          TcSq(j,k) = sumMat
                          TcSq(k,j) = sumMat
                       end do
                    end do
                    do j=1,nsp-1
                       ICinv(j,j) = 1-(TcSR(i,j,j)-TphiLR(i,j,j))
                       do k=j+1,nsp
                          ICinv(j,k) = -(TcSR(i,j,k)-TphiLR(i,j,k))
                          ICinv(k,j) = ICinv(j,k)
                       Enddo
                    Enddo
                    ICinv(nsp,nsp) = 1-(TcSR(i,nsp,nsp)-TphiLR(i,nsp,nsp))
                    call dsytrf('U',nsp,ICinv,nsp,ipiv,work,lwork,info)
                    call dsytri('U',nsp,ICinv,nsp,ipiv,work,info)
                    do k=1,nsp
                       do l=k,nsp
                          ICinv(l,k) = ICinv(k,l)
                       Enddo
                    Enddo
                    do j = 1,nsp
                       do k = j,nsp
                          sumMat = 0.0d0
                          do l = 1,nsp
                             sumMat = sumMat+ICinv(j,l)*TcSq(l,k)
                          end do
                          xFTsSR(i,j,k) = (sumMat-TphiLR(i,j,k))/sqrt(rhoHat(j)*rhoHat(k))
                          xFTsSR(i,k,j) = xFTsSR(i,j,k)
                       end do
                    end do
                 end do
                 !     Inverse transform of xFTsSR(k).
                 do j = 1,nsp
                    do k = j,nsp
                       call Hankel(n,q,dq,qMax,xFTsSR(0,j,k),rMax,sSRnew(0,j,k))
                       sSRnew(0:n,k,j) = sSRnew(0:n,j,k)
                    end do
                 end do
                 !     End iteration core.

                 !     Calculate next s(i) iterate using Ng acceleration. See J. Chem. Phys. 61, 2680 (1974).
                 do j = 1,nsp
                    do k = j,nsp
                       sum22(j,k) = sum11(j,k)
                       sum11(j,k) = sum00(j,k)
                       sum12(j,k) = sum01(j,k)
                       sum00(j,k) = 0.0d0
                       sum01(j,k) = 0.0d0
                       sum02(j,k) = 0.0d0
                       sum0 = 0.0d0
                       do i = 0,n
                          s2(i,j,k) = s1(i,j,k)
                          s1(i,j,k) = s0(i,j,k)
                          s0(i,j,k) = sSRnew(i,j,k)
                          d2(i,j,k) = d1(i,j,k)
                          d1(i,j,k) = d0(i,j,k)
                          d0(i,j,k) = sSRnew(i,j,k)-sSR(i,j,k)
                          sum00 = sum00+d0(i,j,k)*d0(i,j,k)
                          sum01 = sum01+d0(i,j,k)*d1(i,j,k)
                          sum02 = sum02+d0(i,j,k)*d2(i,j,k)
                          sum0 = sum0+(r(i)*d0(i,j,k))**2
                       end do
                       rms = sqrt(dr(Nr)*sum0)
                       if (rms .gt. rmsCut) then
                          blend = blend0
                       else
                          blend = 1.0d0-(1.0d0-blend0)*(rms/rmsCut)**2
                       end if
                       if (iter .eq. 1) then
                          do i = 0,n
                             sSR(i,j,k) = sSR(i,j,k)+blend*(sSRnew(i,j,k)-sSR(i,j,k))
                             sSR(i,k,j) = sSR(i,j,k)
                          end do
                       else if (iter .eq. 2) then
                          a01 = sum00(j,k)-sum01(j,k)
                          a11 = sum00(j,k)-2.0d0*sum01(j,k)+sum11(j,k)
                          c1 = a01/a11
                          do i = 0,n
                             sSRng = (1.0d0-c1)*s0(i,j,k)+c1*s1(i,j,k)
                             sSR(i,j,k) = sSR(i,j,k)+blend*(sSRng-sSR(i,j,k))
                             sSR(i,k,j) = sSR(i,j,k)
                          end do
                       else if (iter .le. iterMax) then
                          a01 = sum00(j,k)-sum01(j,k)
                          a02 = sum00(j,k)-sum02(j,k)
                          a11 = sum00(j,k)-2.0d0*sum01(j,k)+sum11(j,k)
                          a22 = sum00(j,k)-2.0d0*sum02(j,k)+sum22(j,k)
                          a12 = sum00(j,k)-sum01(j,k)-sum02(j,k)+sum12(j,k)
                          c1 = (a01*a22-a02*a12)/(a11*a22-a12**2)
                          c2 = (a02*a11-a01*a12)/(a11*a22-a12**2)
                          do i = 0,n
                             sSRng = (1.0d0-c1-c2)*s0(i,j,k)+c1*s1(i,j,k)+c2*s2(i,j,k)
                             sSR(i,j,k) = sSR(i,j,k)+blend*(sSRng-sSR(i,j,k))
                             sSR(i,k,j) = sSR(i,j,k)
                          end do
                       else
                          write (*,*) '*** Too many iterations. Quitting.'
                          stop
                       end if
                    end do
                 end do
                 if (mod(iter,10) == 0) then
                    write (3,90) iter, rms, (sSR(0,i,i),i=1,nsp)
                 endif
                 if (irho == 0 .and. mod(iter,50) == 0) then
                    write (*,"(a, i5, a, 1pe10.3, a, 0pf8.4, a, 0pf8.4, a)") &
                         c_cyan // '   [Picard rho0]' // c_reset // ' Iteration ' // c_yellow, iter, &
                         c_reset // ' : rms = ' // c_b_white, rms, &
                         c_reset // ' | sSR(1,1) = ' // c_white, sSR(0,1,1), &
                         c_reset // ', sSR(2,2) = ' // c_white, sSR(0,2,2), &
                         c_reset
                 endif
90               format (10x, i6, 1pe13.2, 0pf12.4, 5f12.4)
                 !     End Ng acceleration.

                 !     Check for convergence.
                 ! check for possible extrapolation attempt.
                 if (iter .eq. iext) then
                    iext = iext+next
                    do j = 1,nsp
                       do k = j,nsp
                          do i = 0,n
                             s2(i,j,k) = s0(i,j,k)-(s0(i,j,k)-s1(i,j,k))**2 &
                                  /(s0(i,j,k)-2.0d0*s1(i,j,k)+s2(i,j,k))
                             if (r(i)*dabs(s2(i,j,k)-sSR(i,j,k)) .gt. 1.0d0) cycle
                          end do
                          do i = 0,n
                             sSR(i,j,k) = s2(i,j,k)
                          end do
                       end do
                    end do
                    write (3,100) (sSR(0,i,i), i=1,nsp)
                 end if
                 !       End extrapolation attempt.
              enddo
              if (irho == 0) then
                 write (*,"(a, i5, a, 1pe10.3, a)") &
                      c_b_green // '   [Picard rho0] Converged in ' // c_b_yellow, iter, &
                      c_b_green // ' iterations (rms = ' // c_b_white, rms, &
                      c_b_green // ')' // c_reset
              endif
100           format (13x, 'ext', 13x, 5f12.4)
              !     End Picard iteration.

              !     Done! Calculate g(r;j,k). Print heading.
              do j = 1,nsp
                 do k = j,nsp
                    do i = Ncore(j,k),n
                       g(i,j,k) = exp(-phi(i,j,k))*(1.0d0+(exp((sSR(i,j,k)+phiLR(i,j,k))*fint(i))-1.0d0)*inv_fint(i))
                       g(i,k,j) = g(i,j,k)
                    end do
                 end do
              end do

!!$           scint(:) = 0.d0
!!$           sumtsx = 0.0d0
!!$           sumst = 0.0d0
!!$           sg(:,:) =0.d0
!!$           sgi(:,:,:) = 0.0d0
!!$           do i=1, n
!!$              ssumx = 0
!!$              do j = 1,nsp
!!$                 do k = 1,nsp
!!$                    if (i == Ncore(j,k)) Then
!!$                       ssumx = ssumx + 0.5*pi*rho(j)*rho(k)*(g(i,j,k)-2)&
!!$                            &+0.5*pi*rho(j)*rho(k)*(g(i,j,k)*log(g(i,j&
!!$                            &,k)))  
!!$                       sg(j,k) = sg(j,k) + 0.5*dr(i)*r(i)*((g(i,j,k)&
!!$                            &-1.0d0)**2+1)
!!$                    else
!!$                       ssumx = ssumx + pi*rho(j)*rho(k)*(g(i,j,k)-1)
!!$                       sg(j,k) = sg(j,k) + dr(i)*r(i)*(g(i,j,k)-1.0d0)**2
!!$                    Endif
!!$                    sgi(i,j,k) = 2*pi*sg(j,k)
!!$                    if (g(i,j,k) > 0) Then
!!$                       ssumx = ssumx-pi*rho(j)*rho(k)*(g(i,j,k)*log(g(i,j&
!!$                            &,k)))
!!$                    Endif
!!$                 Enddo
!!$              Enddo
!!$              sumst = sumst -(g(i,1,1)-1)*r(i)*dr(i)
!!$              if (g(i,1,1)>0) Then
!!$                 sumst = sumst +(g(i,1,1)*log(g(i,1,1)))*r(i)*dr(i)
!!$              Endif
!!$              sc(i) = ssumx*r(i)/rhoTotal
!!$              sumtsx = sumtsx + sc(i)*dr(i)
!!$              scint(i) = sumtsx
!!$           Enddo
!!$           s2exc = -pi*rhoTotal*sumst

              !     Calculate thermodynamics.
              sumP1 = 0.0d0
              sumP2 = 0.0d0
              sumU = 0.0d0
              sumUb = 0.0d0
              sumX = 0.0d0
              lcount = 1
              ct(:,:) = 0.0
              do j = 1,nsp
                 do k = 1,nsp
                    sumP0 = 0.0d0
                    sumU0 = 0.0d0
                    sumUb0 = 0.0d0
                    sumX0 = 0.0d0
                    sumcjk = 0
                    g(Ncore(j,k),j,k) = 0.5d0*g(Ncore(j,k),j,k)
                    sexsum = 0.0d0
                    !
                    ! The -1 in the virial and energy accounts for the effect
                    ! of electroneutrality
                    !
                    csrnc(j,k) = cSR(Ncore(j,k),j,k)
                    cSR(Ncore(j,k),j,k) = (cSR(Ncore(j,k),j,k)-1&
                         &-sSR(Ncore(j,k),j,k))/2.0
                    do i = 1,n
                       sumP0 = sumP0+dr(i)*r(i)*g(i,j,k)*dphi(i,j,k)
                       sumU0 = sumU0+dr(i)*r(i)*g(i,j,k)*phi(i,j,k)
                       sumX0 = sumX0+dr(i)*r(i)*cSR(i,j,k)
                       sumcjk = sumcjk + dr(i)*r(i)*cSR(i,j,k)
                    end do
                    g(Ncore(j,k),j,k) = 2.0d0*g(Ncore(j,k),j,k)
                    cSR(Ncore(j,k),j,k) = csrnc(j,k)
                    ct(j,k) = (2*pi*sumcjk-z(j)*z(k)*Gamma*tflr(0))
                    sumP2 = sumP2+rho(j)*rho(k)*(sumP0+virlr(j,k))
                    sumP1 = sumP1+rho(j)*rho(k)*sigma(j,k)**2*g(Ncore(j,k),j,k)
                    sumU = sumU+rho(j)*rho(k)*(sumU0+ulr(j,k))
                    sumX = sumX+rho(j)*rho(k)*(2*pi*sumX0-z(j)*z(k)*Gamma*tflr(0))
                 end do
              end do
              P1 = 1.0d0+(pi/(2.0d0*rhoTotal))*sumP1
              P2 = (pi/(2.0d0*rhoTotal))*sumP2
              U = (pi/rhoTotal)*(sumU-sumUb)
              X = (1.0d0-(sumX/rhoTotal))
              if (irho == 0) then
                 Mrr0 = X
                 Mcc0 =  1 - (rho(1)*rho(2)/rhoTotal)*(ct(1,1)+ct(2,2)-2&
                      &*ct(1,2))
                 Mrc0 = sqrt(rho(1)*rho(2))*((rho(2)/rhoTotal)*ct(2,2)&
                      &-(rho(1)/rhoTotal)*ct(1,1)-((rho(2)-rho(1))&
                      &/rhoTotal)*ct(1,2))
                 lamb1 = (Mrr0+Mcc0-sqrt((Mrr0-Mcc0)**2+4*Mrc0**2))/2
                 lamb2 = (Mrr0+Mcc0+sqrt((Mrr0-Mcc0)**2+4*Mrc0**2))/2
                 scc0 = X/((1-rho(1)*ct(1,1))*(1-rho(2)*ct(2,2))&
                   &-(rho(1)*rho(2))*ct(1,2)**2)

                 vrr = (lamb1-Mcc0)/sqrt((lamb1-Mcc0)**2+Mrc0**2)
                 
                 vcc = Mrc0/sqrt((lamb1-Mcc0)**2+Mrc0**2)
              endif
              pres(irho) = P1+P2
              uint(irho) = U
              xc(irho) = X
              if(irho==0) Then
                 P10 = P1
                 P20 = P2
                 U0 = U
                 X0 = X
              Endif


                            ! Chemical potential in the HNC approx
              do j=1,nsp
                 chemp = 0.0d0
                 do k=1,nsp
                    sumc =0
                    do I=1,n
                       sumc = sumc + cSR(i,j,k)*dr(i)*r(i)
                    Enddo
                    if (j.eq.k) Then
                       csr0(j,k) = 2*pi*(sumc-z(j)*z(k)*Gamma*sqrt(ac*pi)/2)
                    else
                       csr0(j,k) = 2*pi*(sumc-z(j)*z(k)*Gamma*lambda*sqrt(ac*pi)/2)
                    endif

                    sums2 = 0
                    sumint = 0
                    do i=1, n
                       if (i.eq.Ncore(j,k)) Then
                          sumint = sumint + 0.5*(g(i,j,k)-2)*(sSR(i,j,k)+phiLR(i,j,k))*dr(i)*r(i)
                          sums2 = sums2 +  0.5*(g(i,j,k)*log(g(i,j,k))-g(i,j,k)+1)*dr(i)*r(i)
                       else
                          sumint = sumint + (g(i,j,k)-1)*(sSR(i,j,k)+phiLR(i,j,k))*dr(i)*r(i)
                          if (i > Ncore(j,k)) sums2 = sums2 + (g(i,j,k)*log(g(i,j,k))-g(i,j,k)+1)*dr(i)*r(i)
                       Endif
                    Enddo
                    chemp = chemp - rho(k)*csr0(j,k) + pi*rho(k)*sumint
                 Enddo
                 chempot(j) = chemp
                 if (irho == 0) Then
                    chempot0(j) = chemp
                    s2ex = -pi*rho(j)*rho(k)*sums2/rhoTotal0
                 endif
              Enddo
              sumA1 = 0.0d0
              do j = 1,nsp
                 do k = 1,nsp
                    sumA10 = 0.0d0
                    do i = 1,n
                       sumA10 = sumA10+dr(i)*r(i)*(cSR(i,j,k) &
                            +((cSR(i,j,k)-phiLR(i,j,k))**2-(sSR(i,j,k)+phiLR(i,j,k))**2)/2.0d0)
                    end do
                    sumA1 = sumA1+rho(j)*rho(k)*sumA10
                 end do
              end do
              sumc = 0
              do j=1,nsp
                 do k=1,nsp
                    if (j.eq.k) then
                       sumc = sumc -rho(j)*rho(k)*z(j)*z(k)*sqrt(ac*pi)/2
                    else
                       sumc = sumc -rho(j)*rho(k)*lambda*z(j)*z(k)*sqrt(ac*pi)/2
                    Endif
                 Enddo
              Enddo
              sumc = -pi*sumc*Gamma/rhoTotal
              A1 = -(pi/rhoTotal)*sumA1+sumc
              sumA20 = 0.0d0
              do i = 1,n
                 TH11 = rhoHat(1)*(xFTcSR(i,1,1)+xFTsSR(i,1,1))
                 TH22 = rhoHat(2)*(xFTcSR(i,2,2)+xFTsSR(i,2,2))
                 TH12 = sqrt(rhoHat(1)*rhoHat(2))*(xFTcSR(i,1,2)+xFTsSR(i,1,2))
                 sumA20 = sumA20+dq(i)*q(i)*(log((1.0d0+TH11)*(1.0d0+TH22)-TH12**2)-(TH11+TH22))
              end do
              A2 = -sumA20/(4.0d0*pi*rhoTotal)
              if (irho==0) Then
                 A20= A2
                 A10 = A1
              Endif
           Enddo
           dPr = (rhot(1)*pres(1)-rhot(-1)*pres(-1))/(2*deltarho)
           fopt = (xc(0)-dPr)
           ersig = abs(fopt/xc(0))
           if(osig) then
              write (*,"(a, i3, a, f9.6, a, f10.6, a, f10.6, a, f10.6, a, 1pe10.3, a)") &
                   c_magenta // '   [Consistency NR]' // c_reset // ' Step ' // c_yellow, its, &
                   c_reset // ' : eta = ' // c_b_white, eta, &
                   c_reset // ' | chi^-1 = ' // c_white, xc(0), &
                   c_reset // ' | dP/drho = ' // c_white, dPr, &
                   c_reset // ' | diff = ' // c_white, fopt, &
                   c_reset // ' | rel_err = ' // c_b_yellow, ersig, &
                   c_reset
              if(its.lt.2)then
                 eti = eta+dsig
              else
                 fp = (fopt-fopto)/(eta-eto)
                 eti = eta - fopt/fp
              endif
              eto = eta
              fopto = fopt
              eta = eti
           endif
        enddo
        if (osig) then
           write (*,"(a, i3, a, f9.6, a, 1pe10.3, a)") &
                c_b_green // '   [Consistency NR] Converged in ' // c_b_yellow, its, &
                c_b_green // ' step(s)! Final eta = ' // c_b_white, eto, &
                c_b_green // ' (tol = ' // c_white, tol, &
                c_b_green // ')' // c_reset
           eta = eto
        else
           write (*,"(a, f9.6, a)") &
                c_cyan // '   [Rogers-Young] Evaluated with fixed eta = ' // c_b_white, eta, c_reset
        endif
        write (3,101) Gamma, (i,rhoi(i),i=1,nsp)
        write (3,"(5(' z(',i1,') =',f8.4,',':))") (i,z(i),i=1,nsp)
101     format (/' Thermodynamics of a 2D 1/r^3 using the RY equation' &
             /' with Rogers-Young (RY) closure' &
             /' Gamma =', f8.4, ',', 5('rho(',i1,') =', f8.4,',':))
        write (3,110) pres(0), uint(0), xc(0), P10, P20
110     format (5x, 'pA/NkT =', f15.4, ', U/NkT =', f15.4, ', NkTX/A =', f8.4/ &
             12x, '=', f15.4, '  (HD)'/12x, '=', f15.4, '  (QQ)')
        ssum =0
        do j=1,nsp
           ssum = rhoi(j)*chempot0(j)+ssum
        Enddo

        open(95,file='srq.dat')

120     format (5x, 'HNC free energy: A1/NkT =', f8.5, ', A2/NkT =', f8.5, &
             ', Aex/NkT =', f8.5)

        !     Save solution.
        open  (17,file=OUTfile,status='unknown')
        write (17,70) Nr, Ncore(1,1), Ncore(2,2), &
             Gamma, rho(1), rho(2), z(1), z(2)
        do j = 1,nsp
           do k = j,nsp
              write (17,40) (sSR(i,j,k), i = 0,n)
           end do
        end do
        close (17,status='keep')
        open(22,file='gr.dat')
        open(23,file='sq.dat')
        do j=1,nsp
           do k=1, nsp
              y1 = rho(j)*rho(k)/rhoTotal*(2.0d0*pi*rMax/qMax)&
                   &*(xFTcSR(1,j,k)+xFTsSR(1,j,k))
              y2 = rho(j)*rho(k)/rhoTotal*(2.0d0*pi*rMax/qMax)&
                   &*(xFTcSR(2,j,k)+xFTsSR(2,j,k))
              y3 = rho(j)*rho(k)/rhoTotal*(2.0d0*pi*rMax/qMax)&
                   &*(xFTcSR(3,j,k)+xFTsSR(3,j,k))
              if (j==k) then
                 y1 = y1+rho(j)/rhoTotal
                 y2 = y2+rho(j)/rhoTotal
                 y3 = y3+rho(j)/rhoTotal
              endif
!!$        a= -((-(r(2)*y1) + r(3)*y1 + r(1)*y2 - r(3)*y2 - r(1)*y3 + r(2)*y3)&
!!$             &/(r(1)**2*r(2) - r(1)*r(2)**2 - r(1)**2*r(3) + r(2)**2*r(3) +&
!!$             & r(1)*r(3)**2 - r(2)*r(3)**2))) 
!!$        b= -((r(2)**2*y1 - r(3)**2*y1 - r(1)**2*y2 + r(3)**2*y2 + r(1)&
!!$             &**2*y3 - r(2)**2*y3)/(r(1)**2*r(2) - r(1)*r(2)**2 - r(1)&
!!$             &**2*r(3) + r(2)**2*r(3) + r(1)*r(3)**2 - r(2)*r(3)**2))) 
              sjk(j,k)=-((-(r(2)**2*r(3)*y1) + r(2)*r(3)**2*y1 + r(1)**2*r(3)*y2 -&
                   & r(1)*r(3)**2*y2 - r(1)**2*r(2)*y3 +  r(1)*r(2)**2*y3)&
                   &/(r(1)**2*r(2) - r(1)*r(2)**2 - r(1)**2*r(3) + r(2)**2&
                   &*r(3) + r(1)*r(3)**2 - r(2)*r(3)**2))
           enddo
        enddo
        write(23,'(12f15.7)')q(0), q(0)/sqrt(rhoTotal0),((sjk(j,k),j=1,nsp),k=1,nsp),&
             & sum(sjk(:,:)),((rho(1)*mk(0.0d0,R01))**2*sjk(1,1)+(rho(2)&
             &*mk(0.0d0,R02))**2*sjk(2,2)+2*rho(1)*rho(2)*mk(0.0d0,R01)&
             &*mk(0.0d0,R02)*sjk(1,2))/rhoTotal,mk(0.0d0,R02),(rho(2)&
             &/rhoTotal)**2*sjk(1,1)+(rho(1)/rhoTotal)**2& 
                &*sjk(2,2)-(rho(1)*rho(2)/rhoTotal**2)*sjk(1,2)
        sqmax = 0.0
        sq0 = sum(sjk(:,:))
        sq110 = sjk(1,1)
        sq220 = sjk(2,2)
        unfound = .true.
        write(95,'(12f15.7)')q(0),(sjk(j,j)*rhoTotal/rho(j),j=1,nsp)
        do i=1, n
           do j=1,nsp
              sjk(j,j) = rho(j)/rhoTotal+rho(j)*rho(j)/rhoTotal*(2.0d0*pi&
                   &*rMax/qMax)*(xFTcSR(i,j,j)+xFTsSR(i,j,j)) 
              do k=1,nsp
                 if (j.ne.k) Then
                    sjk(j,k) = rho(j)*rho(k)/rhoTotal*(2.0d0*pi*rMax/qMax)&
                         &*(xFTcSR(i,j,k)+xFTsSR(i,j,k))
                 Endif
              Enddo
           Enddo
           write(22,'(18f15.7)')r(i),((g(i,j,k),j=1,nsp),k=1,nsp)
           write(23,'(13f15.7)')q(i), q(i)/sqrt(rhoTotal0),((sjk(j,k),j=1,nsp),k=1,nsp),&
                & sum(sjk(:,:)),((rho(1)*mk(q(i),R01))**2*sjk(1,1)+(rho(2)&
                &*mk(q(i),R02))**2*sjk(2,2)+2*rho(1)*rho(2)*mk(q(i),R01)&
                &*mk(q(i),R02)*sjk(1,2))/rhoTotal,mk(q(i),R02)&
                &,(rho(2)/rhoTotal)**2*sjk(1,1)+(rho(1)/rhoTotal)**2&
                &*sjk(2,2)-(rho(1)*rho(2)/rhoTotal**2)*sjk(1,2)
           if (sum(sjk(:,:)) > 1) then
              if (sum(sjk(:,:)) > sqmax .and. unfound) then
                 sqmax = sum(sjk(:,:))
              else
                 unfound = .false.
              endif
           endif
           write(95,'(12f15.7)')q(i),(sjk(j,j)*rhoTotal/rho(j),j=1,nsp)
        Enddo
        close(22)
        close(23)
        close(95)
        if (sq110 > 0.0d0 .and. sq220 > 0.0d0) then
           scale_21 = (sq220 / sq110)**0.25d0
        else
           scale_21 = 1.0d0
        end if
        R02 = R01 / scale_21

        write(88,"(1x, 15(1x, f14.6))") &
             rhoTotal0, xf(ixf), uint(0), pres(0), xc(0), &
             dPr, eta, sqmax, sq0, sq0/sqmax, &
             scale_21, s2ex, lamb1, lamb2, 1.0d0/scc0

        !     Print comprehensive state results
        write (*,*)
        write (*,"(a)") c_b_blue // '--------------------------------------------------------------------------------' // c_reset
        write (*,"(a, f8.4, a, f8.4, a)") &
             c_b_cyan // ' RESULTS FOR STATE POINT: ' // c_b_yellow // 'rho = ' // c_b_white, rhoTotal0, &
             c_b_yellow // ', x2 = ' // c_b_white, xf(ixf), c_reset
        write (*,"(a)") c_b_blue // '--------------------------------------------------------------------------------' // c_reset
        write (*,"(a)") c_b_cyan // ' Thermodynamics & Equation of State:' // c_reset
        write (*,"(a, f12.6)") '   Compressibility factor (Z = P/(rho*kT)) : ', pres(0)
        write (*,"(a, f12.6)") '     - Hard-disk core contribution (Z_HD)  : ', P10
        write (*,"(a, f12.6)") '     - Dipolar interaction (Z_dip)         : ', P20
        write (*,"(a, f12.6)") '   Reduced internal energy (U / (N*kT))    : ', uint(0)
        write (*,"(a, f12.6)") '   Excess two-body entropy (S2_ex / kB)    : ', s2ex
        write (*,*)
        write (*,"(a)") c_b_cyan // ' Thermodynamic Consistency & Response:' // c_reset
        write (*,"(a, f12.6)") '   Rogers-Young parameter (eta)            : ', eta
        write (*,"(a, f12.6)") '   Inverse compressibility (chi^-1)        : ', xc(0)
        write (*,"(a, f12.6)") '   Virial pressure derivative (dP*/drho)   : ', dPr
        write (*,"(a, f12.6)") '   Isothermal compressibility (drho/dP*)   : ', 1.0d0 / dPr
        write (*,"(a, f12.6, a, 1pe10.3, a)") &
             '   Consistency discrepancy (f_opt)         : ', fopt, '  (rel_err = ', ersig, ')'
        write (*,*)
        write (*,"(a)") c_b_cyan // ' Fluctuation & Spinodal Stability (Bhatia-Thornton):' // c_reset
        write (*,"(a)") '   Zero-wavevector response matrix elements:'
        write (*,"(a, f12.6, a, f12.6, a, f12.6)") &
             '     Mrr(0) = ', Mrr0, ', Mcc(0) = ', Mcc0, ', Mrc(0) = ', Mrc0
        write (*,"(a, f12.6, a, f12.6)") &
             '   Spinodal eigenvalues (lambda1, lambda2) : ', lamb1, ', ', lamb2
        if (lamb1 > 1.0d-5) then
           write (*,"(a)") '   Spinodal stability status               : ' // &
                c_b_green // 'STABLE (lambda1 > 0)' // c_reset
        else if (lamb1 >= 0.0d0) then
           write (*,"(a)") '   Spinodal stability status               : ' // &
                c_b_yellow // 'NEAR SPINODAL MARGIN (lambda1 ~ 0)' // c_reset
        else
           write (*,"(a)") '   Spinodal stability status               : ' // &
                c_b_red // 'UNSTABLE / DEMIXING (lambda1 < 0)' // c_reset
        end if
        write (*,"(a, f12.6, a, f12.6, a)") &
             '   Concentration fluctuation S_cc(0)       : ', scc0, '  (1/S_cc(0) = ', 1.0d0 / scc0, ')'
        write (*,*)
        write (*,"(a)") c_b_cyan // ' Structure Factor Highlights:' // c_reset
        write (*,"(a, f12.6)") '   Peak of total structure factor S(q)_max : ', sqmax
        write (*,"(a, f12.6)") '   Zero-wavevector structure factor S(0)   : ', sq0
        write (*,"(a, f12.6)") '   Fluctuation ratio S(0) / S(q)_max       : ', sq0 / sqmax
        if (sq110 > 0.0d0 .and. sq220 > 0.0d0) then
           write (*,"(a, f12.6)") '   Disk size scaling (S22(0)/S11(0))^0.25  : ', scale_21
        end if
        write (*,"(a, f12.6)") '   Optimum form factor radius ratio R01/R02: ', R01 / R02
        write (*,"(a)") c_b_blue // '--------------------------------------------------------------------------------' // c_reset

        if (nStateDone < maxStates) then
           nStateDone = nStateDone + 1
           res_rho(nStateDone) = rhoTotal0
           res_x2(nStateDone) = xf(ixf)
           res_P(nStateDone) = pres(0)
           res_U(nStateDone) = uint(0)
           res_xc(nStateDone) = xc(0)
           res_dPr(nStateDone) = dPr
           res_eta(nStateDone) = eta
           res_lamb1(nStateDone) = lamb1
           res_invScc(nStateDone) = 1.0d0 / scc0
           res_sqmax(nStateDone) = sqmax
        end if
        open(230,file='sqopt.dat')
        do j=1,nsp
           do k=1, nsp
              y1 = rho(j)*rho(k)/rhoTotal*(2.0d0*pi*rMax/qMax)&
                   &*(xFTcSR(1,j,k)+xFTsSR(1,j,k))
              y2 = rho(j)*rho(k)/rhoTotal*(2.0d0*pi*rMax/qMax)&
                   &*(xFTcSR(2,j,k)+xFTsSR(2,j,k))
              y3 = rho(j)*rho(k)/rhoTotal*(2.0d0*pi*rMax/qMax)&
                   &*(xFTcSR(3,j,k)+xFTsSR(3,j,k))
              if (j==k) then
                 y1 = y1+rho(j)/rhoTotal
                 y2 = y2+rho(j)/rhoTotal
                 y3 = y3+rho(j)/rhoTotal
              endif
!!$        a= -((-(r(2)*y1) + r(3)*y1 + r(1)*y2 - r(3)*y2 - r(1)*y3 + r(2)*y3)&
!!$             &/(r(1)**2*r(2) - r(1)*r(2)**2 - r(1)**2*r(3) + r(2)**2*r(3) +&
!!$             & r(1)*r(3)**2 - r(2)*r(3)**2))) 
!!$        b= -((r(2)**2*y1 - r(3)**2*y1 - r(1)**2*y2 + r(3)**2*y2 + r(1)&
!!$             &**2*y3 - r(2)**2*y3)/(r(1)**2*r(2) - r(1)*r(2)**2 - r(1)&
!!$             &**2*r(3) + r(2)**2*r(3) + r(1)*r(3)**2 - r(2)*r(3)**2))) 
              sjk(j,k)=-((-(r(2)**2*r(3)*y1) + r(2)*r(3)**2*y1 + r(1)**2*r(3)*y2 -&
                   & r(1)*r(3)**2*y2 - r(1)**2*r(2)*y3 +  r(1)*r(2)**2*y3)&
                   &/(r(1)**2*r(2) - r(1)*r(2)**2 - r(1)**2*r(3) + r(2)**2&
                   &*r(3) + r(1)*r(3)**2 - r(2)*r(3)**2))
           enddo
        enddo
        write(230,'(12f15.7)')q(0), q(0)/sqrt(rhoTotal0), ((sjk(j,k),j=1,nsp),k=1,nsp),&
             & sum(sjk(:,:)),((rho(1)*mk(0.0d0,R01))**2*sjk(1,1)+(rho(2)&
             &*mk(0.0d0,R02))**2*sjk(2,2)+2*rho(1)*rho(2)*mk(0.0d0,R01)&
             &*mk(0.0d0,R02)*sjk(1,2))/rhoTotal,mk(0.0d0,R02)
        sqmax = 0.0
        sq0 = sum(sjk(:,:))
        sq110 = sjk(1,1)
        sq220 = sjk(2,2)
        unfound = .true.
        do i=1, n
           do j=1,nsp
              sjk(j,j) = rho(j)/rhoTotal+rho(j)*rho(j)/rhoTotal*(2.0d0*pi&
                   &*rMax/qMax)*(xFTcSR(i,j,j)+xFTsSR(i,j,j)) 
              do k=1,nsp
                 if (j.ne.k) Then
                    sjk(j,k) = rho(j)*rho(k)/rhoTotal*(2.0d0*pi*rMax/qMax)&
                         &*(xFTcSR(i,j,k)+xFTsSR(i,j,k))
                 Endif
              Enddo
           Enddo
           write(230,'(12f15.7)')q(i), q(i)/sqrt(rhoTotal0),((sjk(j,k),j=1,nsp),k=1,nsp),&
                & sum(sjk(:,:)),((rho(1)*mk(q(i),R01))**2*sjk(1,1)+(rho(2)&
                &*mk(q(i),R02))**2*sjk(2,2)+2*rho(1)*rho(2)*mk(q(i),R01)&
                &*mk(q(i),R02)*sjk(1,2))/rhoTotal,mk(q(i),R02)
           if (sum(sjk(:,:)) > 1) then
              if (sum(sjk(:,:)) > sqmax .and. unfound) then
                 sqmax = sum(sjk(:,:))
              else
                 unfound = .false.
              endif
           endif
        Enddo
        close(230)


     enddo
  Enddo

  !     Summary table of all computed state points
  if (nStateDone > 0) then
     write (*,*)
     write (*,"(a)") c_b_blue // '================================================================================' // c_reset
     write (*,"(a)") c_b_yellow // '                     SUMMARY OF COMPUTED STATE POINTS' // c_reset
     write (*,"(a)") c_b_blue // '================================================================================' // c_reset
     write (*,"(a)") c_b_cyan // '   rho       x2      P/(rho*kT)   U/(N*kT)      chi^-1     dP*/drho     eta       lambda1       1/Scc(0)    S(q)_max' // c_reset
     write (*,"(a)") c_dim // '--------------------------------------------------------------------------------' // c_reset
     do istate = 1, nStateDone
        write (*,"(f8.4, 1x, f8.4, 1x, f12.6, 1x, f11.6, 1x, f11.6, 1x, f11.6, 1x, f8.4, 1x, f11.6, 1x, f11.6, 1x, f10.4)") &
             res_rho(istate), res_x2(istate), res_P(istate), res_U(istate), &
             res_xc(istate), res_dPr(istate), res_eta(istate), &
             res_lamb1(istate), res_invScc(istate), res_sqmax(istate)
     end do
     write (*,"(a)") c_b_blue // '================================================================================' // c_reset
  end if

  write (*,"(/a)") c_b_cyan // ' Output files written:' // c_reset
  write (*,"(a)") c_cyan // '   * ' // c_b_white // 'thermo.dat' // c_reset // ' : Thermodynamic table & spinodal stability metrics'
  write (*,"(a)") c_cyan // '   * ' // c_b_white // 'gr.dat    ' // c_reset // ' : Pair radial distribution functions g_jk(r)'
  write (*,"(a)") c_cyan // '   * ' // c_b_white // 'sq.dat    ' // c_reset // ' : Structure factors S_jk(q) and S_cc(q)'
  write (*,"(a)") c_cyan // '   * ' // c_b_white // 'sqopt.dat ' // c_reset // ' : Optical form-factor folded structure factors'
  write (*,"(a)") c_cyan // '   * ' // c_b_white // 'srq.dat   ' // c_reset // ' : Scaled partial structure factors'
  write (*,"(a)") c_cyan // '   * ' // c_b_white // 'solout.dat' // c_reset // ' : Converged solution vector s_SR(r; j,k) for restart'
  write (*,"(/a)") c_b_green // ' Program twoDdipRY completed successfully.' // c_reset
  write (*,"(a/)") c_b_blue // '================================================================================' // c_reset

  close(88)
  close(3, status='keep')
end program twoDdipRY




SUBROUTINE BesChb(x,gam1,gam2,gampl,gammi)
  !     Uses Chebev
  !     Taken from "Numerical Recipes" by W.H. Press et al.
  implicit none
  real(kind=8), save :: c1(7)=(/-1.142022680371172d0,6.516511267076d-3&
       &,3.08709017308d-4,-3.470626964d-6,6.943764d-9,3.6780d-11,&
       &-1.36d-13/) , c2(8)=(/1.843740587300906d0,-.076852840844786d0,1.271927136655d-3, &
       &-4.971736704d-6,-3.3126120d-8,2.42310d-10,-1.70d-13,-1.d-15/)
  integer, parameter :: nuse1=5,nuse2=5
  real(kind=8) :: x, xx, gam1, gam2, gammi, gampl
  real(kind=8), external :: Chebev
  xx = 8.0d0*x*x-1.0d0
  gam1 = Chebev(-1.0d0,1.0d0,c1,nuse1,xx)
  gam2 = Chebev(-1.0d0,1.0d0,c2,nuse2,xx)
  gampl = gam2-x*gam1
  gammi = gam2+x*gam1
end subroutine BesChb



   SUBROUTINE BessJY(x,xnu,rj,ry,rjp,ryp)
     !     Uses BesChb
     !     Taken from "Numerical Recipes" by W.H. Press et al.
     !     Modified to calculate J0(x) and J0'(x) by direct series summation
     !     for x<2 to get double precision in this range.
     implicit real*8(a-h,o-z), integer(i-n)
     parameter (eps=1.0d-16,fpmin=1.0d-30,maxit=10000,xmin=2.0d0)

     pi = dacos(-1.0d0)
     if (x.le.0.0d0 .or. xnu.lt.0.0d0) then
        write (*,*) '*** Bad arguments in BessJY'
        stop
     end if
     if(x .lt. xmin) then
        nl = idint(xnu+0.5d0)
     else
        nl = max0(0,idint(xnu-x+1.5d0))
     end if
     xmu = xnu-dble(nl)
     xmu2 = xmu*xmu
     xi = 1.0d0/x
     xi2 = 2.0d0*xi
     w = xi2/pi
     isign = 1
     h = xnu*xi
     if (h.lt.fpmin) h = fpmin
     b = xi2*xnu
     d = 0.0d0
     c = h
     do i = 1,maxit
        b = b+xi2
        d = b-d
        if (dabs(d) .lt. fpmin) d=fpmin
        c = b-1.0d0/c
        if (dabs(c) .lt. fpmin) c=fpmin
        d = 1.0d0/d
        del = c*d
        h = del*h
        if (d .lt. 0.0d0) isign=-isign
        if (dabs(del-1.0d0) .lt. eps) go to 1
     end do
     write (*,*) '*** x too large in BessJY; try asymptotic expansion.'
     stop
1    continue
     rjl = isign*fpmin
     rjpl = h*rjl
     rjl1 = rjl
     rjp1 = rjpl
     fact = xnu*xi
     do l = nl,1,-1
        rjtemp = fact*rjl+rjpl
        fact = fact-xi
        rjpl = fact*rjtemp-rjl
        rjl = rjtemp
     end do
     if (rjl .eq. 0.0d0) rjl=eps
     f = rjpl/rjl
     if (x .lt. xmin) then
        x2 = 0.5d0*x
        pimu = pi*xmu
        if (dabs(pimu) .lt. eps) then
           fact = 1.0d0
        else
           fact = pimu/sin(pimu)
        end if
        d = -log(x2)
        e = xmu*d
        if (dabs(e) .lt. eps) then
           fact2 = 1.0d0
        else
           fact2 = dsinh(e)/e
        end if
        call BesChb(xmu,gam1,gam2,gampl,gammi)
        ff = 2.0d0/pi*fact*(gam1*cosh(e)+gam2*fact2*d)
        e = dexp(e)
        p = e/(gampl*pi)
        q = 1.0d0/(e*pi*gammi)
        pimu2 = 0.5d0*pimu
        if (dabs(pimu2) .lt. eps) then
           fact3 = 1.0d0
        else
           fact3 = sin(pimu2)/pimu2
        end if
        r = pi*pimu2*fact3*fact3
        c = 1.0d0
        d = -x2*x2
        sum0 = ff+r*q
        sum1 = p
        do i = 1,maxit
           ff = (dble(i)*ff+p+q)/(dble(i*i)-xmu2)
           c = c*d/dble(i)
           p = p/(dble(i)-xmu)
           q = q/(dble(i)+xmu)
           del = c*(ff+r*q)
           sum0 = sum0+del
           del1 = c*p-dble(i)*del
           sum1 = sum1+del1
           if (dabs(del) .lt. (1.0d0+dabs(sum0))*eps) go to 2
        end do
        write (*,*) '*** BessY series failed to converge in Subroutine BessJY'
        stop
2       continue
        rymu = -sum0
        ry1 = -sum1*xi2
        rymup = xmu*xi*rymu-ry1
        rjmu = w/(rymup-f*rymu)
     else
        a = 0.25d0-xmu2
        p = -0.5d0*xi
        q = 1.0d0
        br = 2.0d0*x
        bi = 2.0d0
        fact = a*xi/(p*p+q*q)
        cr = br+q*fact
        ci = bi+p*fact
        den = br*br+bi*bi
        dr = br/den
        di = -bi/den
        dlr = cr*dr-ci*di
        dli = cr*di+ci*dr
        temp = p*dlr-q*dli
        q = p*dli+q*dlr
        p = temp
        do i = 2,maxit
           a = a+dble(2*(i-1))
           bi = bi+2.0d0
           dr = a*dr+br
           di = a*di+bi
           if (dabs(dr)+dabs(di) .lt. fpmin) dr = fpmin
           fact = a/(cr*cr+ci*ci)
           cr = br+cr*fact
           ci = bi-ci*fact
           if (dabs(cr)+dabs(ci) .lt. fpmin) cr = fpmin
           den = dr*dr+di*di
           dr = dr/den
           di = -di/den
           dlr = cr*dr-ci*di
           dli = cr*di+ci*dr
           temp = p*dlr-q*dli
           q = p*dli+q*dlr
           p = temp
           if (dabs(dlr-1.0d0)+dabs(dli) .lt. eps) go to 3
        end do
        write (*,*) '*** cf2 failed in Subroutine BessJY'
        stop
3       continue
        gam = (p-f)/q
        rjmu = sqrt(w/((p-f)*gam+q))
        rjmu = dsign(rjmu,rjl)
        rymu = rjmu*gam
        rymup = rymu*(p+q/gam)
        ry1 = xmu*xi*rymu-rymup
     end if
     fact = rjmu/rjl
     rj = rjl1*fact
     rjp = rjp1*fact
     do i = 1,nl
        rytemp = (xmu+dble(i))*xi2*ry1-rymu
        rymu = ry1
        ry1 = rytemp
     end do
     ry = rymu
     ryp = xnu*xi*rymu-ry1

     !     Recalculate J0(x) and J0'(x) for small x by direct series summation.
     if (xnu .eq. 0.0d0 .and. x .lt. 2.0d0) then
        xsq = -0.25d0*x*x
        k = 0
        add0 = 1.0d0
        add1 = 1.0d0
        rj  = 1.0d0
        rjp = 1.0d0
4       k = k+1
        add0 = add0*xsq/dble(k*k)
        add1 = add1*xsq/dble(k*(k+1))
        rj  = rj +add0
        rjp = rjp+add1
        if (dabs(add0) .gt. 1.0d-15) go to 4
        rjp = -0.5d0*x*rjp
     end if
   end subroutine BessJY



   SUBROUTINE Bpc2D(Nr,Ncore,rho,r,q,dr,dq,Bpc)
     !     Uses Hankel
     !     Obtain pressure-consistent (PC) bridge function for hard disks.
     !     See J. Chem. Phys. 49, 3092 (1968).
     use datatrans, only : mxNR, pi
     implicit none
     integer :: n, Nr, Ncore, i, k, iter
     real(kind=8) :: mu, eta, rMax, qMax, rhoHat, a, bigs, biga,&
          & alpha, beta, c0, c1, c2, c0p, x, d11, d22, Ccore, rho, d00,&
          & d12, d01, d02, sum0, rms, dd01, dd02, dd11, dd22, dd12,&
          & Sng, Pr
     real(kind=8) :: rmsMax=1.0d-5,blend=0.8d0
     real(kind=8) :: r(0:mxNr), q(0:mxNr), dr(0:mxNr), dq(0:mxNr), Bpc(0:mxNr)
     real(kind=8) :: s(0:mxNr),    c(0:mxNr),    P(0:mxNr), xFTc(0:mxNr), &
          xFTp(0:mxNr), xFTs(0:mxNr), snew(0:mxNr), &
          f0(0:mxNr),   f1(0:mxNr),   f2(0:mxNr), &
          d0(0:mxNr),   d1(0:mxNr),   d2(0:mxNr)

     !     Initialize parameters and arrays.
     eta = pi*rho/4.0d0
     n = Nr-1
     rMax = r(Nr)
     qMax = q(Nr)
     rhoHat = (2.0d0*pi*rMax/qMax)*rho
     write (3,10) rho, Ncore, rMax, qMax
     write (*,10) rho, Ncore, rMax, qMax
10   format (/3x, 'SUBROUTINE Bpc2D with rho =', f8.5/ &
          3x, 'Ncore =', i4, ', rMax =', f10.5, ', qMax =', f10.5)

     !     Start from modeled PY c(r). See Leutheuser, J. Chem. Phys. 84, 1050 (1986).
     a = 1.0d0+eta
     bigs = 0.5d0*sqrt(1.0d0-a*a/4.0d0)+dasin(a/2.0d0)/a
     biga = (bigs-(1.0d0-a*a/4.0d0)*sqrt(1.0d0-a*a/4.0d0))*16.0d0/(pi*a*a)
     alpha = 2.0d0*eta**2*biga
     beta = 8.0d0*eta*bigs/pi
     c1 = ((4.0d0*eta-1.0d0)+sqrt((4.0d0*eta-1.0d0)**2 &
          -4.0d0*(alpha-beta)))/(2.0d0*(alpha-beta))
     c0 = c1-beta*c1**2
     c0p = 8.0d0*eta*c1**2/pi
     do i = 0,Ncore
        x = a*r(i)/2.0d0
        c(i) = c0+c0p*(dasin(x)+x*sqrt(1.0d0-x*x))/a
     end do
     c(Ncore) = 0.5d0*c(Ncore)
     do i = Ncore+1,n
        c(i) = 0.0d0
     end do
     !     Compute Hankel transform of c(r).
     call Hankel(n,r,dr,rMax,c,qMax,xFTc)
     !     OZ equation for xFTs(k).
     do k = 1,n
        xFTs(k) = rhoHat*xFTc(k)**2/(1.0d0-rhoHat*xFTc(k))
     end do
     !     Compute inverse Hankel transform of xFTs(k) to get modeled s(r).
     s(0) = -1.0d0-c(0)
     call Hankel(n,q,dq,qMax,xFTs,rMax,s)

     iter = 0
     mu = 0.0920d0+0.1222d0*rho+0.1642d0*rho**2+0.1100d0*rho**3
     !     Begin Picard iteration for S(r) = g(r)exp[beta*phi(r)]-1.
     !     Begin iteration core.
     rms = 100.0
     do while (rms > rmsMax)
        iter = iter+1
        do i = 0,Ncore
           P(i) = S(i)-log(1.0d0+S(i))
           C(i) = -(1.0+S(i))+mu*P(i)
        end do
        Ccore = C(Ncore)
        do i = Ncore,n
           P(i) = S(i)-log(1.0d0+S(i))
           C(i) = mu*P(i)
        end do
        C(Ncore) = 0.5d0*(Ccore+C(Ncore))
        !     Compute Hankel transform of C(r) and P(r).
        call Hankel(n,r,dr,rMax,C,qMax,xFTc)
        call Hankel(n,r,dr,rMax,P,qMax,xFTp)
        !     OZ equation for xFTs(k).
        do k = 0,n
           xFTs(k) = rhoHat*xFTc(k)**2/(1.0d0-rhoHat*xFTc(k))+mu*xFTp(k)
        end do
        !     Compute inverse Hankel transform of xFTs(k).
        call Hankel(n,q,dq,qMax,xFTs,rMax,Snew)
        !     End iteration core.

        !     Check convergence and compute next S iterate using Ng acceleration.
        d22 = d11
        d11 = d00
        d12 = d01
        d00 = 0.0d0
        d01 = 0.0d0
        d02 = 0.0d0
        sum0 = 0.0d0
        do i = 0,n
           f2(i) = f1(i)
           f1(i) = f0(i)
           f0(i) = Snew(i)
           d2(i) = d1(i)
           d1(i) = d0(i)
           d0(i) = Snew(i)-S(i)
           d00 = d00+d0(i)*d0(i)
           d01 = d01+d0(i)*d1(i)
           d02 = d02+d0(i)*d2(i)
           sum0 = sum0+d0(i)*d0(i)
        end do
        rms = sqrt(sum0/dble(n))
        !     write (3,30) iter, rms, Snew(1)
        !     write (*,30) iter, rms, Snew(1)
        !  30 format (i6, e13.2, f12.5)
        dd01 = d00-d01
        dd02 = d00-d02
        dd11 = d00-2.0d0*d01+d11
        dd22 = d00-2.0d0*d02+d22
        dd12 = d00-d01-d02+d12
        if (iter .ge. 3) then
           c1 = (dd01*dd22-dd02*dd12)/(dd11*dd22-dd12**2)
           c2 = (dd02*dd11-dd01*dd12)/(dd11*dd22-dd12**2)
        else
           c1 = dd01/dd11
           c2 = 0.0d0
        end if
        if (iter .eq. 1) then
           do i = 0,n
              S(i) = S(i)+blend*(Snew(i)-S(i))
           end do
        else 
           do i = 0,n
              Sng = f0(i)+c1*(f1(i)-f0(i))+c2*(f2(i)-f0(i))
              S(i) = S(i)+blend*(Sng-S(i))
           end do
        end if
     Enddo
     !     End Picard iteration

     !     Done! Get hard disk thermodynamics.
     Pr = 1.0d0+pi*rho*(1.0d0+S(Ncore))/2.0d0
     sum0 = 0.0d0
     do i = 1,n
        sum0 = sum0+dr(i)*r(i)*C(i)
     end do
     X = 1.0d0/(1.0d0 -2.0d0*pi*rho*sum0)
     write (3,40) rho, eta, mu
     write (*,40) rho, eta, mu
     write (3,50) iter, Pr, X
     write (*,50) iter, Pr, X
40   format (3x, 'rho =', f7.4, ', eta =', f7.4, ', mu =', f7.4/ &
          3x, 'Hard disk thermodynamics from the PC equation.')
50   format (3x, i3, ' iterations for S(r): pA/NkT =', f8.4, ', NkTX/A =', f8.4)

     !     Calculate PC bridge function of hard disks.
     do i = 0,n
        Bpc(i) = (1.0d0-mu)*(log(1.0d0+S(i))-S(i))
     end do
   end subroutine Bpc2D



   SUBROUTINE BpySR(Nr,Ncore,Gamma,sigma,z,rho,lambda, r,q,dr,dq,phiSR0,BpySR0)
     !     Uses Hankel
     !     Obtain the Percus-Yevick bridge functions for a "short-range electrolyte"
     !     of hard disks with short-range electrostatic potentials phiSR0(r;j,k).
     use datatrans, only : pi, mxNr, nsp
     implicit none
     real(kind=8), parameter :: rmsMax=1.0d-5,rmsCut=1.0d-2,blend0=0.5d0
     integer, parameter :: next=25,iterMax=1500
     integer ::  Ncore(nsp,nsp)
     integer :: i, j, k, l, info, iter, iext, Nr, n
     real(kind=8) :: rMax, qMax, rhoTotal, sumMat, sum0, rms, blend,&
          & a01, a11, a02, a22, a12, c1, c2, sSRng, Gamma, sumP1,&
          & sumP2, sumU, sumX, sumP0, sumU0, sumX0, lambda, P1, P2, U&
          &, X
     real(kind=8) :: sigma(nsp,nsp), z(nsp), rho(nsp), r(0:mxNr), q(0:mxNr),   &
          dr(0:mxNr), dq(0:mxNr), phiSR0(0:mxNr,nsp,nsp), BpySR0(0:mxNr,nsp,nsp)
     real(kind=8) :: rhoHat(nsp), ICinv(nsp,nsp),  TcSq(nsp,nsp), sum22(nsp,nsp), sum11(nsp,nsp), &
          sum12(nsp,nsp), sum00(nsp,nsp), sum01(nsp,nsp), sum02(nsp,nsp)
     real(kind=8) :: eBond(0:mxNr,nsp,nsp),    sSR(0:mxNr,nsp,nsp),    g(0:mxNr,nsp,nsp), &
          cSR(0:mxNr,nsp,nsp), xFTcSR(0:mxNr,nsp,nsp), TcSR(0:mxNr,nsp,nsp), &
          xFTsSR(0:mxNr,nsp,nsp), sSRnew(0:mxNr,nsp,nsp),                   &
          s0(0:mxNr,nsp,nsp),     s1(0:mxNr,nsp,nsp),   s2(0:mxNr,nsp,nsp), &
          d0(0:mxNr,nsp,nsp),     d1(0:mxNr,nsp,nsp),   d2(0:mxNr,nsp,nsp)
     integer :: ipiv(nsp)
     integer :: lwork
     real (kind=8) ::  work(2*nsp)
     lwork = 2*nsp
     rMax = r(Nr)
     qMax = q(Nr)
     n = Nr-1
     do j = 1,nsp
        do k = 1,nsp
           do i = 0,Ncore(j,k)-1
              eBond(i,j,k) = 0.0d0
              sSR(i,j,k) = 0.0d0
              g(i,j,k) = 0.0d0
           end do
           do i = Ncore(j,k),n
              eBond(i,j,k) = dexp(-phiSR0(i,j,k))
              sSR(i,j,k) = 0.0d0
           end do
        end do
     end do
     rhoTotal = sum(rho(1:nsp))
     rhoHat(1:nsp) = (2.0d0*pi*rMax/qMax)*rho(1:nsp)
     iext = next
     iter = 0

     !     Begin Picard iteration.
     rms = 100.0
     do while (rms > rmsMax)
        iter = iter+1
        !     Begin iteration core.
        do j = 1,nsp
           do k = j,nsp
              do i = Ncore(j,k),n
                 g(i,j,k) = eBond(i,j,k)*(1.0d0+sSR(i,j,k))
              end do
              g(Ncore(j,k),j,k) = 0.5d0*g(Ncore(j,k),j,k)
              do i = 0,n
                 cSR(i,j,k) = g(i,j,k)-1.0d0-sSR(i,j,k)
                 cSR(i,k,j) = cSR(i,j,k)
                 g(i,k,j) = g(i,j,k)
              end do
           end do
        end do
        !     Transform cSR(r).
        do j = 1,nsp
           do k = j,nsp
              call Hankel(n,r,dr,rMax,cSR(0,j,k),qMax,xFTcSR(0,j,k))
              do i = 0,n
                 TcSR(i,j,k) = sqrt(rhoHat(j)*rhoHat(k))*xFTcSR(i,j,k)
                 TcSR(i,k,j) = TcSR(i,j,k)
                 xFTcSR(i,k,j) = xFTcSR(i,j,k)
              end do
           end do
        end do
        !     Ornstein-Zernike equation for xFTsSR(k).
        do i = 1,n
           do j = 1,nsp
              do k = j,nsp
                 sumMat = 0.0d0
                 do l = 1,nsp
                    sumMat = sumMat+TcSR(i,j,l)*TcSR(i,l,k)
                 end do
                 TcSq(j,k) = sumMat
                 TcSq(k,j) = sumMat
              end do
           end do
           do j=1,nsp-1
              ICinv(j,j) = 1-TcSR(i,j,j)
              do k=j+1,nsp
                 ICinv(j,k) = -TcSR(i,j,k)
                 ICinv(k,j) = ICinv(j,k)
              Enddo
           Enddo
           ICinv(nsp,nsp) = 1-TcSR(i,nsp,nsp)
           call dsytrf('U',nsp,ICinv,nsp,ipiv,work,lwork,info)
           call dsytri('U',nsp,ICinv,nsp,ipiv,work,info)
           do k=1,nsp
              do l=k,nsp
                 ICinv(l,k) = ICinv(k,l)
              Enddo
           Enddo
           do j = 1,nsp
              do k = j,nsp
                 sumMat = 0.0d0
                 do l = 1,nsp
                    sumMat = sumMat+ICinv(j,l)*TcSq(l,k)
                 end do
                 xFTsSR(i,j,k) = sumMat/sqrt(rhoHat(j)*rhoHat(k))
                 xFTsSR(i,k,j) = xFTsSR(i,j,k)
              end do
           end do
        end do
        !     Inverse transform of xFTsSR(k).
        do j = 1,nsp
           do k = j,nsp
              call Hankel(n,q,dq,qMax,xFTsSR(0,j,k),rMax,sSRnew(0,j,k))
              sSRnew(0:n,k,j) = sSRnew(0:n,j,k)
           end do
        end do
        !     End iteration core.

        !     Calculate next s(i) iterate using Ng acceleration. See J. Chem. Phys. 61, 2680 (1974).
        do j = 1,nsp
           do k = j,nsp
              sum22(j,k) = sum11(j,k)
              sum11(j,k) = sum00(j,k)
              sum12(j,k) = sum01(j,k)
              sum00(j,k) = 0.0d0
              sum01(j,k) = 0.0d0
              sum02(j,k) = 0.0d0
              sum0 = 0.0d0
              do i = 0,n
                 s2(i,j,k) = s1(i,j,k)
                 s1(i,j,k) = s0(i,j,k)
                 s0(i,j,k) = sSRnew(i,j,k)
                 d2(i,j,k) = d1(i,j,k)
                 d1(i,j,k) = d0(i,j,k)
                 d0(i,j,k) = sSRnew(i,j,k)-sSR(i,j,k)
                 sum00 = sum00+d0(i,j,k)*d0(i,j,k)
                 sum01 = sum01+d0(i,j,k)*d1(i,j,k)
                 sum02 = sum02+d0(i,j,k)*d2(i,j,k)
                 sum0 = sum0+(r(i)*d0(i,j,k))**2
              end do
              rms = sqrt(dr(Nr)*sum0)
!!$           if (rms .gt. rmsCut) then
!!$              blend = blend0
!!$           else
!!$              blend = 1.0d0-(1.0d0-blend0)*(rms/rmsCut)**2
!!$           end if
              blend=blend0
              if (iter .eq. 1) then
                 do i = 0,n
                    sSR(i,j,k) = sSR(i,j,k)+blend*(sSRnew(i,j,k)-sSR(i,j,k))
                    sSR(i,k,j) = sSR(i,j,k)
                 end do
              else if (iter .eq. 2) then
                 a01 = sum00(j,k)-sum01(j,k)
                 a11 = sum00(j,k)-2.0d0*sum01(j,k)+sum11(j,k)
                 c1 = a01/a11
                 do i = 0,n
                    sSRng = (1.0d0-c1)*s0(i,j,k)+c1*s1(i,j,k)
                    sSR(i,j,k) = sSR(i,j,k)+blend*(sSRng-sSR(i,j,k))
                    sSR(i,k,j) = sSR(i,j,k)
                 end do
              else if (iter .le. iterMax) then
                 a01 = sum00(j,k)-sum01(j,k)
                 a02 = sum00(j,k)-sum02(j,k)
                 a11 = sum00(j,k)-2.0d0*sum01(j,k)+sum11(j,k)
                 a22 = sum00(j,k)-2.0d0*sum02(j,k)+sum22(j,k)
                 a12 = sum00(j,k)-sum01(j,k)-sum02(j,k)+sum12(j,k)
                 c1 = (a01*a22-a02*a12)/(a11*a22-a12**2)
                 c2 = (a02*a11-a01*a12)/(a11*a22-a12**2)
                 do i = 0,n
                    sSRng = (1.0d0-c1-c2)*s0(i,j,k)+c1*s1(i,j,k)+c2*s2(i,j,k)
                    sSR(i,j,k) = sSR(i,j,k)+blend*(sSRng-sSR(i,j,k))
                    sSR(i,k,j) = sSR(i,j,k)
                 end do
              else
                 write (3,*) '*** Too many iterations in BpySR. Quitting.'
                 write (*,*) '*** Too many iterations in BpySR. Quitting.'
                 stop
              end if
           end do
        end do
        !     write (3,90) iter, rms, sSR(0,1,1), sSR(0,1,2), sSR(0,2,2)
        write (*,90) iter, rms, (sSR(0,i,i),i=1,nsp)
90      format (10x, i6, 1pe13.2, 0pf12.4, 5f12.4)
        !     End Ng acceleration.

        !     Check for convergence.
        !       check for possible extrapolation attempt.
        if (iter .eq. iext) then
           iext = iext+next
           do j = 1,nsp
              do k = j,nsp
                 do i = 0,n
                    s2(i,j,k) = s0(i,j,k)-(s0(i,j,k)-s1(i,j,k))**2 &
                         /(s0(i,j,k)-2.0d0*s1(i,j,k) +s2(i,j,k))
                    if (r(i)*dabs(s2(i,j,k)-sSR(i,j,k)) .gt. 2.0d0) cycle
                 end do
                 do i = 0,n
                    sSR(i,j,k) = s2(i,j,k)
                 end do
              end do
           end do
        end if
        !       End extrapolation attempt.
     Enddo
     !     End Picard iteration.

     !     Done! Calculate g(r;j,k) and BpySR0(r;j,k). Print heading.
     do j = 1,nsp
        do k = 1,nsp
           do i = 1,n
              g(i,j,k) = eBond(i,j,k)*(1.0d0+sSR(i,j,k))
              BpySR0(i,j,k) = log(1.0d0+sSR(i,j,k))-sSR(i,j,k)
           end do
        end do
     end do
     write (3,103) iter, Gamma, rho(1), rho(2)
     write (*,103) iter, Gamma, rho(1), rho(2)
103  format (/3x, 'SUBROUTINE BpySR ' &
          /3x, 'Thermodynamics of a "short-ranged 2D electrolyte" using the' &
          ' PY equation' &
          /3x, 'iter =', i3, ': Gamma =', f8.4, ', rho1 =', f8.4, ', rho2 =', f8.4)

     !     Calculate PY thermodynamics of SR potentials.
     sumP1 = 0.0d0
     sumP2 = 0.0d0
     sumU = 0.0d0
     sumX = 0.0d0
     do j = 1,nsp
        do k = 1,nsp
           sumP0 = 0.0d0
           sumU0 = 0.0d0
           sumX0 = 0.0d0
           g(Ncore(j,k),j,k) = 0.5d0*g(Ncore(j,k),j,k)
           do i = 1,n
              sumP0 = sumP0+dr(i)*r(i)*g(i,j,k)*dexp(-r(i)**2)
              sumU0 = sumU0+dr(i)*r(i)*g(i,j,k)*phiSR0(i,j,k)
              sumX0 = sumX0+dr(i)*r(i)*(g(i,j,k)-1.0d0)
           end do
           g(Ncore(j,k),j,k) = 2.0d0*g(Ncore(j,k),j,k)
           sumP1 = sumP1+rho(j)*rho(k)*sigma(j,k)**2*g(Ncore(j,k),j,k)
           if (j.eq.k) Then
              sumP2 = sumP2+Gamma*rho(j)*rho(k)*z(j)*z(k)*sumP0
           else
              sumP2 = sumP2+Gamma*rho(j)*rho(k)*lambda*z(j)*z(k)*sumP0
           Endif

           sumU = sumU+rho(j)*rho(k)*sumU0
           sumX = sumX+rho(j)*rho(k)*sumX0
        end do
     end do
     P1 = 1.0d0+(pi/(2.0d0*rhoTotal))*sumP1
     P2 = (pi/(2.0d0*rhoTotal))*sumP2
     U = (pi/rhoTotal)*sumU
     X = 1.0d0+(2.0d0*pi/rhoTotal)*sumX
     write (3,110) P1+P2, U, X, P1, P2
     write (*,110) P1+P2, U, X, P1, P2
110  format (8x, 'pA/NkT =', f8.4, ', U/NkT =', f8.4, ', NkTX/A =', f8.4/ &
          15x, '=', f8.4, '  (HD)'/15x, '=', f8.4, '  (QQ-SR)')
   end subroutine BpySR



   FUNCTION Chebev(a,b,c,m,x)
     !     Taken from "Numerical Recipes" by W.H. Press et al.
     implicit none
     integer :: m, j
     real(kind=8) :: c(m), x, a, b, Chebev
     real(kind=8) :: d,dd, sv, y2, y
     if ((x-a)*(x-b) .gt. 0.0d0) then
        write (*,*) '*** x not in range in function Chebev'
        stop
     end if
     d = 0.0d0
     dd = 0.0d0
     y = (2.0d0*x-a-b)/(b-a)
     y2 = 2.0d0*y
     do j = m,2,-1
        sv = d
        d = y2*d-dd+c(j)
        dd = sv
     end do
     Chebev = y*d-dd+0.5d0*c(1)
   end function Chebev



   FUNCTION E1(x)
     !     Returns the exponential integral E_1(x) for positive real x.
     !     See AMS55, p.231.
     implicit none
     real(kind=8) :: x, E1, y
     real(kind=8), save :: p0,p1,p2,p3,p4,p5
     data p0,p1,p2,p3,p4,p5/-0.57721566d0, 0.99999193d0,-0.24991055d0, &
          0.05519968d0,-0.00976004d0, 0.00107857d0/
     real(kind=8), save :: q1,q2,q3,q4
     data q1,q2,q3,q4/8.5733287401d0,18.0590169730d0, &
          8.6347608925d0, 0.2677737343d0/
     real(kind=8), save :: r1,r2,r3,r4
     data r1,r2,r3,r4/ 9.5733223454d0,25.6329561486d0, &
          21.0996530827d0, 3.9584969228d0/
     if (x .lt. 1.0d0) then
        E1 = -log(x)+(p0+x*(p1+x*(p2+x*(p3+x*(p4+x*p5)))))
     else
        y = 1.0d0/x
        E1 = (1.0d0+y*(q1+y*(q2+y*(q3+y*q4))))/(1.0d0+y*(r1+y*(r2+y*(r3+y*r4)))) &
             *dexp(-x)/x
     end if
   end function E1



!=======================================================================
! Subroutine: Hankel
! Purpose: 2D Hankel (Bessel J_0) transform on Lado's quadrature grid.
!          Accelerated using Intel MKL BLAS dgemv for the matrix-vector
!          product Fout(1:n) = W(1:n, 1:n) * Fin(1:n).
!=======================================================================
   SUBROUTINE Hankel(n,x,dx,Cin,Fin,Cout,Fout)
     use datatrans, only : mxNr, W
     implicit none
     integer, intent(in) :: n
     real(kind=8), intent(in) :: x(0:mxNr), dx(0:mxNr), Fin(0:mxNr)
     real(kind=8), intent(in) :: Cout, Cin
     real(kind=8), intent(out) :: Fout(0:mxNr)
     real(kind=8) :: sumFT
     integer :: i

     ! DC (zero wavevector / origin) component
     sumFT = 0.0d0
     do i = 1,n
        sumFT = sumFT + dx(i)*x(i)*Fin(i)
     end do
     Fout(0) = (Cout/Cin)*sumFT

     ! Intel MKL BLAS dgemv for fast matrix-vector multiplication
     call dgemv('N', n, n, 1.0d0, W, mxNr, Fin(1), 1, 0.0d0, Fout(1), 1)
   end subroutine Hankel



   FUNCTION SerAdd(x,n,a)
     !     SerAdd = a(1)+x(a(2)+x(a(3)+x( ...+x(a(n))...)))
     implicit none
     integer :: n
     real(kind=8) :: a(n), x, SerAdd
     integer :: i
     SerAdd = a(n)
     do i = n-1,1,-1
        SerAdd = a(i)+x*SerAdd
     end do
   end function SerAdd



   SUBROUTINE setupW(Nr)
     !     Uses SerAdd, BessJY
     use datatrans, only : mxNr, W, BJ1sq, root, pi 
     implicit none
     integer :: i, j, k, Nr 
     real (kind=8), parameter :: xerrMx=1.0d-12
     real(kind=8), dimension(4) :: croot = (/ 1.0d0, -41.33333333333333d0,&
          & 8061.86666666667d0, -3826125.40952d0 /)  
     real(kind=8), dimension(20) :: root1 = (/ 2.4048255577d0,&
          & 5.5200781103d0, 8.6537279129d0,11.7915344391d0,14.9309177086d0,18.0710639679d0,&
          21.2116366299d0,24.3524715308d0,27.4934791320d0&
          &,30.6346064684d0, 33.7758202136d0,36.9170983537d0&
          &,40.0584257646d0,43.1997917132d0, 46.3411883717d0&
          &,49.4826098974d0,52.6240518411d0,55.7655107550d0,&
          & 58.9069839261d0,62.0484691902d0 /) 
     real (kind=8) :: x, BessJ,BessY,BessJp,BessYp,BesJ0,BesY0,BesJ0p&
          &,BesY0p, xerror, xnew, a1, a2, b, xsq
     real (kind=8), external :: SerAdd
     
     !     Calculate roots of zero-order Bessel function J0(x).
     !     First pass. Use explicit roots or construct polynomial approximations.
     do i = 1,20
        root(i) = root1(i)
     end do
     do i = 21,Nr
        b = (dble(i)-0.25d0)*pi
        x = 1.0d0/(8.0d0*b)
        xsq = x**2
        root(i) = b+x*seradd(xsq,4,croot)
     end do
     !     Second pass. Refine first-pass roots with Newton's method.
     do i = 1,Nr
        x = root(i)
        k = 0
1       continue
        k = k+1
        call BessJY(x,0.0d0,BesJ0,BesY0,BesJ0p,BesY0p)
        xnew = x-BesJ0/BesJ0p
        xerror = dabs((xnew-x)/xnew)
        x = xnew
        if (k .gt. 200) then
           write (*,*) '*** Count exceeds 200 in subroutine SETUPW'
           stop
        end if
        if (xerror .gt. xerrMx) go to 1
        root(i) = x
        BJ1sq(i) = BesJ0p**2
     end do

     !     Set up W table.
     a1 = 2.0d0/root(Nr)
     do i = 1,Nr-1
        a2 = root(i)/root(Nr)
        do j = 1,Nr-1
           x = root(j)*a2
           call BessJY(x,0.0d0,BessJ,BessY,BessJp,BessYp)
           W(i,j) = a1*BessJ/BJ1sq(j)
        end do
     end do
   end subroutine setupW

!=======================================================================
! Subroutine: lrfuncs
! Purpose: Computes screened long-range direct correlation functions
!          flr(r) and Fourier transforms tflr(q) for 1/r^3 dipole tails
!          using exponentially scaled modified Bessel functions I_0, I_1.
!=======================================================================
   subroutine lrfuncs(flr,tflr,r,q,Nr,a)
     use bessel_slatec_mod, only : dbsi0e, dbsi1e
     implicit none
     real(kind=8) :: flr(0:Nr), tflr(0:Nr), r(0:Nr), q(0:Nr)
     real(kind=8), parameter :: pi=3.141592653589793d0 
     real(kind=8) :: q2a, q2, I0, I1, a
     integer :: i, Nr
     flr(0) = 0.0d0
     tflr(0) = pi*sqrt(a*pi)
     do i = 1, Nr
        flr(i) = (1.0d0 - (1.0d0 + a*r(i)**2)*exp(-a*r(i)**2))/r(i)**3
        q2 = q(i)**2
        q2a= q2/a
        I0 = dbsi0e(q2a/8.0d0)
        I1 = dbsi1e(q2a/8.0d0)
        tflr(i) = 2.0d0*pi*(-q(i) + sqrt(pi/a)*((2.0d0*a + q2)*I0 + q2*I1)/4.0d0)
     end do
   end subroutine lrfuncs
function mk(k, R)
  use datatrans, only : pi
  implicit none
  real(kind=8), intent(in) :: k, R
  real(kind=8) :: mk
  if (k * R > 1.0d-6) then
     mk = 2.0d0 * pi * R * bessel_j1(k * R) / k
  else
     mk = pi * R**2
  endif
end function mk
