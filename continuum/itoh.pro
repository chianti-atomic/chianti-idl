
PRO itoh, temp, wvl, int, no_setup=no_setup, sumt=sumt, dem_int=dem_int, $
          photons=photons, min_abund=min_abund, kev=kev, em_int=em_int

;+
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory (USA), the University of Florence (Italy), the
;       University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
; NAME
;
;    ITOH
;
; EXPLANATION
;
;    Calculates the relativistic free-free continuum using the fitting 
;    formula of Itoh et al. (ApJS 128, 125, 2000).
;
;    To calculate the continuum spectrum the user should specify
;    whether the plasma structure is specified through a differential
;    emission measure (DEM) curve or through emission measure (EM)
;    values. For the former, the DEM values are specified through the
;    DEM_INT keyword and for the latter the EM values are specified
;    through the EM_INT keyword. If neither are specified EM_INT
;    values of unity are assumed.
;
;    If DEM_INT values are specifed they must be specified on a
;    temperature scale that has constant intervals in log10(T). EM_INT
;    values do not have this restriction.
;
;    Note that EM_INT must be used for isothermal plasmas.
;
; INPUTS
;
;    TEMP    Temperature (in K). If the keyword DEM_INT= is set, then the
;            temperatures must be specified with a fixed spacing in
;            logT. E.g., TEMP=10.^[4.0,4.1,4.2].
;
;    WVL     Wavelengths in angstroms. Can be a scalar or vector. If /keV is
;            set, then WVL is assumed to contain energies in keV units.
;
; OUTPUTS
;
;    INT     Free-free continuum emissivity in units of 
;            10^-40 erg cm^3 / s / sr / Angstrom per unit emission 
;            measure [ integral(N_e N_H dh) in cm^-5 ]. If T is given as 
;            a 1-D array, then RAD will be output as a 2-D array, 
;            with one element for each temperature and wavelength 
;            (but also see SUMT).
;
;            If the keyword /keV is set, then the units of INT will be 
;            10^-40 erg cm^3/s/sr/keV
;
; OPTIONAL INPUTS
;
;    DEM_INT   This should be the same size as TEMP and contain
;              differential emission measure values. Specifying
;              DEM_INT places a restriction on TEMP (see above). Note
;              that DEM_INT values are multiplied by the factor dT
;              within the routine, whereas EM_INT values are not.
;
;    EM_INT    This should be the same size as TEMP and contain
;              emission measure values. Note
;              that DEM_INT values are multiplied by the factor dT
;              within the routine, whereas EM_INT values are not. If
;              neither EM_INT or DEM_INT are specified, then a EM_INT
;              vector of same size as TEMP with values of unity is
;              assumed.
;
;    MIN_ABUND This keyword allows the specification of a minimum abundance, 
;              such that any elements with an abundance (relative to 
;              hydrogen) less than MIN_ABUND will not be included in the 
;              calculation. E.g., MIN_ABUND=1e-5.
;
; KEYWORDS
;
;    NO_SETUP By default the routine asks the user which ion balance 
;             and abundance files to use via pop-up widgets. If 
;             /no_setup is used then this data is taken from the common 
;             block.
;
;    SUMT     The default is to output the intensity array as an array 
;             of size (nwvl x nT). Setting this keyword performs a sum 
;             over the temperatures to yield a vector of same size as 
;             the input wavelengths, thus producing the complete 
;             free-free spectrum.
;
;    PHOTONS  Gives output emissivity in photon units rather than ergs.
;
;    KEV      If set, then WVL is assumed to contain energies in keV units
;             rather than wavelengths in Angstroms.
;
; CALLS
;
;    SETUP_ELEMENTS
;
; PROGRAMMING NOTES
;
;    In the Itoh et al. paper they state that their fitting formula is 
;    valid for ion charges (Z_j) of 1-28. The data-file they provide 
;    actually goes up to Z_j=30, and so you will see that the loop over 
;    z (see below) goes up to 30.
;
;    There is no restriction on the elements for which the fitting 
;    formula is valid, and so all elements are allowed to be included 
;    (subject to their abundances being non-zero).
;
; HISTORY
;
;    Ver.1, 3-Dec-2001, Peter Young
;
;    Ver.2, 22-May-2002, Peter Young
;            Added MIN_ABUND optional input.
;
;    Ver.3, 28-May-2002, Peter Young
;            Corrected way in which DEM is handled.
;
;    Ver.4, 9-Mar-2006, Peter Young
;            Added /keV keyword.
;
;    Ver.5, 10-Sep-2009, Peter Young
;            Routine now checks that temperatures are evenly spaced
;            when DEM_INT is specified. Also 0.1 dex spacing is no
;            longer assumed.
;
;    Ver.6, 6-Oct-2009, Peter Young
;            Introduced EM_INT keyword for specifying emission measure
;            values for isothermal plasmas.
;-


COMMON elements,abund,abund_ref,ioneq,ioneq_t,ioneq_ref

IF n_params() LT 3 THEN BEGIN
  print,'Use:  IDL> itoh, temp, wvl, int [, /no_setup, /sumt, /kev, '
  print,'                        /photons, min_abund= , dem_int= , '
  print,'                        em_int= ]'
  return
ENDIF


temp=double(temp)
wvl=double(wvl)

;
; If /keV specified, then the input WVL needs to be converted to Angstroms
; in order to calculate the free-free continuum. At the end of the routine
; WVL gets converted back to keV using WVL_SAVE.
;
IF keyword_set(kev) THEN BEGIN
  wvl_save=wvl
  efact=12.398/wvl^2   ; hc/E^2 factor to multiply continuum by (remember
                       ; wvl is energies in keV at this stage)
  efact=reverse(efact)
  wvl=12.398/wvl   ; convert energies to angstroms
  wvl=reverse(wvl)
ENDIF ELSE BEGIN
  efact=1d0
ENDELSE

nt=n_elements(temp)
nw=n_elements(wvl)

cc=2.997825d18    ; ang/s
hkt=4.1356692d-15 / 8.61735d-5 / temp


;
; DEM_INT and EM_INT implementation
; ---------------------------------
;  If DEM_INT is specified then it's necessary to multiply the
;  DEM values by the factor dT (=ln(10)*T*d(log10(T)).
;
;  I implement this by introducing a keyword dT. If dT=1 then a DEM is
;  assumed, while dt=0 implies EM values. In the latter case DEM_INT
;  is set to be EM_INT.
;
;  EM_VALS is set to be DEM_INT or EM_INT and is used in the rest of
;  the code.
;
CASE 1 OF 
  n_elements(dem_int) NE 0 AND n_elements(em_int) NE 0: BEGIN
    print,'%ITOH: Please specify only one of DEM_INT and EM_INT. Returning...'
    int=0.
    return
  END
 ;
  n_elements(dem_int) EQ 0 AND n_elements(em_int) NE 0: BEGIN
    em_vals=em_int
    dt=0
  END 
 ;
  n_elements(dem_int) NE 0 AND n_elements(em_int) EQ 0: BEGIN
    em_vals=dem_int
    dt=1
  END 
 ;
  n_elements(dem_int) EQ 0 AND n_elements(em_int) EQ 0: BEGIN
    em_vals=dblarr(nt)+1.
    dt=0
  END 
ENDCASE


;
; Make sure EM_VALS has same size as TEMP
;
em_vals=double(em_vals)
IF n_elements(em_vals) NE nt THEN BEGIN
  print,'%ITOH: Warning, number of elements of DEM_INT (EM_INT) must match the'
  print,'  number of temperatures. Returning...'
  int=0.
  return
ENDIF 


;
; If /DT set then need to multiply by dT for each element of EM_VALS.
;
IF keyword_set(dt) THEN BEGIN
  IF nt LT 2 THEN BEGIN
    print,'%ITOH: Two or more temperatures must be specified if DEM_INT values are specified.'
    print,'       Returning...'
    int=0.
    return
  ENDIF 
 ;
  lt1=alog10(temp[0:nt-2])
  lt2=alog10(temp[1:*])
  mean_lt=mean(lt2-lt1)
  max_LT=max(lt2-lt1)
  min_LT=min(lt2-lt1)
  IF max_LT-mean_LT GE 0.05*mean_LT OR mean_LT-min_LT GE 0.05*mean_LT THEN BEGIN
    print,'%ITOH: Error, LOG10 TEMP needs to be evenly spaced if DEM_INT is specified.'
    print,'       E.g., LOG10 TEMP should be spaced at 0.1 dex intervals.'
    int=0.
    return
  ENDIF 
  demfactor=temp*mean_lt/alog10(exp(1))
END ELSE BEGIN
  demfactor=dblarr(nt)+1.
ENDELSE 
;
; Create dem_arr
;
dem_arr= (em_vals*demfactor) # (dblarr(nw) + 1. )



;
; chk is used to flag locations where the temperature and wavelength 
; are outside of the range of validity of the Itoh fitting formula
;
chk=bytarr(n_elements(temp),nw)

t=(alog10(temp)-7.25)/1.25d0 # (dblarr(nw) +1.)
ind=where( (t LT -1.0) OR (t GT 1.0) )
IF ind[0] NE -1 THEN chk[ind]=1

little_u=hkt # (cc/wvl)
u=(alog10(little_u)+1.5)/2.5d0    ; this is the big-U of Itoh et al.
ind=where( (u LT -1.0) OR (u GT 1.0) )
IF ind[0] NE -1 THEN chk[ind]=1

;
; Note that wvl has been converted to Angstroms stage, so no correction for
; the /kev keyword is required.
; 
IF keyword_set(photons) THEN BEGIN
  erg2phot=( (dblarr(nt) + 1.) # wvl ) / 1.9864d-8
ENDIF ELSE BEGIN
  erg2phot=(dblarr(nt) + 1.) # (dblarr(nw) + 1. )
ENDELSE

hkt=hkt # (dblarr(nw) +1.)

IF NOT keyword_set(no_setup) THEN setup_elements

openr,lun,!xuvtop+'/continuum/itoh.dat',/get_lun
data=dblarr(121,30)    ; 30 is the number of elements (H to Zn)
readf,lun,data
free_lun,lun

int=dblarr(nt,nw)

zmax=max(where(abund gt 0.))+1
;
IF n_elements(min_abund) EQ 0. THEN min_abund=0.
elt_i=where(abund[0:zmax-1] GE min_abund)
nei=n_elements(elt_i)

FOR z=1,30 DO BEGIN
  a=data[*,z-1]
  a=[0d0,a]
  gf= A(1) +A(2)*U +A(3)*U^2 +A(4)*U^3 +A(5)*U^4  $
       +A(6)*U^5+A(7)*U^6 +A(8)*U^7 +A(9)*U^8   $
       +A(10)*U^9 +A(11)*U^10   $
       +( A(12)+A(13)*U+A(14)*U^2+A(15)*U^3+A(16)*U^4  $
          +A(17)*U^5+A(18)*U^6+A(19)*U^7 +A(20)*U^8   $
          +A(21)*U^9 +A(22)*U^10 )*T  $
       +( A(23)+A(24)*U+A(25)*U^2+A(26)*U^3+A(27)*U^4   $
          +A(28)*U^5+A(29)*U^6+A(30)*U^7 +A(31)*U^8   $
          +A(32)*U^9 +A(33)*U^10 )*T^2  $
       +( A(34)+A(35)*U+A(36)*U^2+A(37)*U^3+A(38)*U^4  $
          +A(39)*U^5+A(40)*U^6+A(41)*U^7 +A(42)*U^8   $
          +A(43)*U^9 +A(44)*U^10 )*T^3  $
       +( A(45)+A(46)*U+A(47)*U^2+A(48)*U^3+A(49)*U^4  $
          +A(50)*U^5+A(51)*U^6+A(52)*U^7 +A(53)*U^8   $
          +A(54)*U^9 +A(55)*U^10 )*T^4  $
       +( A(56)+A(57)*U+A(58)*U^2+A(59)*U^3+A(60)*U^4  $
          +A(61)*U^5+A(62)*U^6+A(63)*U^7 +A(64)*U^8   $
          +A(65)*U^9 +A(66)*U^10 )*T^5  $
       +( A(67)+A(68)*U+A(69)*U^2+A(70)*U^3+A(71)*U^4  $
          +A(72)*U^5+A(73)*U^6+A(74)*U^7 +A(75)*U^8   $
          +A(76)*U^9 +A(77)*U^10 )*T^6  $
       +( A(78)+A(79)*U+A(80)*U^2+A(81)*U^3+A(82)*U^4  $
          +A(83)*U^5+A(84)*U^6+A(85)*U^7 +A(86)*U^8   $
          +A(87)*U^9 +A(88)*U^10 )*T^7  $
       +( A(89)+A(90)*U+A(91)*U^2+A(92)*U^3+A(93)*U^4  $
          +A(94)*U^5+A(95)*U^6+A(96)*U^7 +A(97)*U^8   $
          +A(98)*U^9 +A(99)*U^10 )*T^8  $
       +( A(100)+A(101)*U+A(102)*U^2+A(103)*U^3+A(104)*U^4  $
          +A(105)*U^5+A(106)*U^6+A(107)*U^7+A(108)*U^8  $
          +A(109)*U^9+A(110)*U^10)*T^9  $
       +( A(111)+A(112)*U+A(113)*U^2+A(114)*U^3+A(115)*U^4  $
          +A(116)*U^5+A(117)*U^6+A(118)*U^7+A(119)*U^8  $
          +A(120)*U^9+A(121)*U^10)*T^10

  FOR j=0,nei-1 DO BEGIN

    ii=ioneq[*,elt_i[j],z]

    ti=where(ii NE 0.)

    IF ti[0] NE -1 THEN BEGIN

      ind=where( (alog10(temp) LE max(ioneq_t[ti])) AND $ 
                 (alog10(temp) GE min(ioneq_t[ti])) )

      IF ind[0] NE -1 THEN BEGIN
        yy=ii[ti]
        xx=ioneq_t[ti]
        xi=alog10(temp[ind])
        y2=spl_init(xx,yy)
        yi=spl_interp(xx,yy,y2,xi)
        
        contrib=gf[ind,*]*z^2*1.426d-27*abund[elt_i[j]]* $
             (yi # (dblarr(nw) + 1.)) * $
             (sqrt(temp[ind]) # (dblarr(nw) + 1.)) * $
             exp(-little_u[ind,*])*hkt[ind,*]*cc/ $
             ((dblarr(n_elements(yi)) + 1.) # wvl)^2 * $
             erg2phot[ind,*] * dem_arr[ind,*]

        int[ind,*]=int[ind,*]+contrib

      ENDIF

    ENDIF

  ENDFOR
ENDFOR

ind=where(chk NE 0)
IF ind[0] NE -1 THEN int[ind]=0d0

int=transpose(int)*1d40/4d0/!pi

;
; multiply int by the hc/E^2 factor if /kev set
;
IF keyword_set(kev) THEN BEGIN
  e_arr=efact#(dblarr(nt)+1d0)
  int=int*e_arr
 ;
 ; reverse the wavelength dimension to match energy ordering rather than
 ; wavelength ordering
  int=int[reverse(indgen(nw)),*]
 ;
 ; and set wavelengths back to energy units
  wvl=wvl_save
ENDIF

siz=size(int)

IF keyword_set(sumt) AND (siz[0] GT 1) THEN int=total(int,2)

;IF dem_tst EQ 1 THEN junk=temporary(dem_tst)

END
