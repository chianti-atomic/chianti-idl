PRO freebound, temp, wvl, int, sumt=sumt, photons=photons, $
               noverner=noverner, iz=iz, ion=ion, no_setup=no_setup, $
               min_abund=min_abund, dem_int=dem_int, verbose=verbose, $
               kev=kev, em_int=em_int

;+
; PROJECT:  CHIANTI
;
;      CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;      Astrophysical Plasmas. It is a collaborative project involving the Naval
;      Research Laboratory (USA), the University of Florence (Italy), the
;      University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
; NAME
;
;      FREEBOUND
;
; PURPOSE:
;
;      Calculates the free-bound (radiative recombination) continuum.
;
;      A report giving detailed information on how the CHIANTI
;      free-bound is calculated is available at:
;
;      http://www.pyoung.org/chianti/freebound.pdf
;
;      To calculate the continuum spectrum the user should specify
;      whether the plasma structure is specified through a differential
;      emission measure (DEM) curve or through emission measure (EM)
;      values. For the former, the DEM values are specified through the
;      DEM_INT keyword and for the latter the EM values are specified
;      through the EM_INT keyword. If neither are specified EM_INT
;      values of unity are assumed.
;
;      If DEM_INT values are specifed they must be specified on a
;      temperature scale that has constant intervals in log10(T). EM_INT
;      values do not have this restriction.
;
;      Note that EM_INT must be used for isothermal plasmas.
;
; INPUTS
;
;      TEMP    Temperature in K (can be an array). If the input
;              DEM_INT is used then the temperatures must be spaced at
;              equal intervals of 
;              log10(T). E.g., log T=6.0, 6.1, 6.2, 6.3, ...
;
;      WVL     Wavelength in angstroms (can be an array). If /keV is
;              set, then WVL is assumed to contain energies in keV units.
;
; OUTPUTS
;
;      INT     Free-bound continuum intensity in units of 
;              10^-40 erg cm^3/s/sr/Angstrom 
;              ( integral(N_H N_e dh) in cm^-5) if a DEM is not defined. 
;
;              If DEM values are defined, it is assumed that they are given
;              as N_H N_e dh/dT.  The units are 10^-40 erg/cm^2/s/srAngstrom 
;
;              If T is given as a 1-D array, then the output will be a 
;              2-D array, with one element for each temperature and 
;              wavelength (but also see SUMT).
;
;              If the keyword /keV is set, then the units of INT will be 
;              10^-40 erg cm^3/s/sr/keV.
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
;      IZ     Only calculate continuum for the element with atomic 
;             number IZ
;
;      ION    (To be used in conjunction with IZ.) Calculated continuum 
;             for a single ion (IZ, ION).
;
; KEYWORDS
;
;      NO_SETUP If the procedure setup_elements has already been called 
;               then the keyword /nosetup should be set to avoid 
;               repeating this step
;
;      MIN_ABUND If set, calculates the continuum only from those 
;                elements which have an abundance greater than 
;                min_abund.  Can speed up the calculations.  For 
;                example:
;                   abundance (H)  = 1.
;                   abundance (He) = 0.085
;                   abundance (C)  = 3.3e-4
;                   abundance (Si) = 3.3e-5
;                   abundance (Fe) = 3.9e-5
;
;      PHOTONS  The output spectrum is given in photon units rather 
;               than ergs.
;
;      SUMT     When a set of temperatures is given to FREEBOUND, the 
;               default is to output INTENSITY as an array of size 
;               (nwvl x nT). With this keyword set, a summation over 
;               the temperatures is performed.
;
;      VERBOSE  Output information from FREEBOUND.
;
;      KEV      If set, then WVL is assumed to contain energies in keV units
;               rather than wavelengths in Angstroms.
;
; COMMON BLOCKS
;
;      ELEMENTS
;
; CALLS
;
;      FREEBOUND_ION, SETUP_ELEMENTS, READ_KLGFB, GET_IEQ
;
; HISTORY
;
;      Ver.1, 24-Jul-2002, Peter Young
;
;      Ver.2, 26-Jul-2002, Peter Young
;           revised call to freebound_ion; corrected ion fraction problem
;
;       V 3, 25-May-2005, GDZ 
;                  corrected routine header.
;
;      Ver.4, 10-Mar-2006, Peter Young
;           added keV keyword
;
;      Ver.5, 8-Sep-2009, Peter Young
;           routine now checks to make sure temperatures are evenly
;           spaced when the DEM_INT keyword is used; if an error is
;           found due to incorrect temperature specificiations, then
;           the continuum intensity is returned as zero.
;
;      Ver.6, 6-Oct-2009, Peter Young
;            Introduced EM_INT keyword for specifying emission measure
;            values for isothermal plasmas.
;
;      Ver.7, 5-Sep-2012, Peter Young
;            Updated header (no change to code).
;
; VERSION     :  7, 5-Spe-2012
;-

COMMON elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref

IF n_params() LT 3 THEN BEGIN
  print,'Use:  IDL> freebound, temp, wvl, int [/sumt, /photons, /noverner,'
  print,'                 /iz, /ion, /no_setup,  min_abund= , /dem_int, /verbose, '
  print,'                  /kev, /em_int]'
  return
ENDIF


IF n_elements(ion) NE 0 AND n_elements(iz) EQ 0 THEN BEGIN
  print,'%FREEBOUND:  Error, please specify IZ when using ION'
  return
ENDIF

IF NOT keyword_set(no_setup) THEN BEGIN
  setup_elements
  no_setup=1
ENDIF

t1=systime(1)

IF n_elements(min_abund) EQ 0 THEN min_abund=0.

wvl=double(wvl)
temp=double(temp)
nt=n_elements(temp)
nwvl=n_elements(wvl)
int=dblarr(nwvl,nt)

ident_wvl=make_array(nwvl,val=1.)
ident_t=make_array(nt,val=1.)

read_ip,!xuvtop+'/ip/chianti.ip',ionpot,ipref

read_klgfb,pe,klgfb
ksize=size(klgfb)
max_ngf=ksize(2)

IF NOT keyword_set(noverner) THEN BEGIN
  vdata=dblarr(10,465)
  dir=concat_dir(!xuvtop,'continuum')
  fname=concat_dir(dir,'verner_short.txt')
  openr,lun,fname,/get_lun
  readf,lun,vdata
  free_lun,lun
ENDIF


IF n_elements(iz) NE 0 THEN BEGIN
 ;
  ab=abund[iz-1]
  IF ab LT min_abund THEN GOTO,lbl2
  IF n_elements(ion) NE 0 THEN BEGIN
    ieq=get_ieq(temp,iz,ion+1,ioneq_logt=ioneq_logt,ioneq_frac=ioneq)
    ip=ionpot[iz-1,ion-1]
    IF (total(ieq) NE 0.) THEN BEGIN
      freebound_ion,temp,wvl,int,iz,ion,ip=ip, $
           vdata=vdata,pe=pe,klgfb=klgfb,noverner=noverner, kev=kev
      int=int*ab*(ident_wvl#ieq)
    ENDIF
 ;
  ENDIF ELSE BEGIN
    FOR ion=1,iz DO BEGIN
      ieq=get_ieq(temp,iz,ion+1,ioneq_logt=ioneq_logt,ioneq_frac=ioneq)
      ip=ionpot[iz-1,ion-1]
      IF total(ieq) NE 0. THEN BEGIN
        freebound_ion,temp,wvl,rad,iz,ion,ip=ip, $
             vdata=vdata,pe=pe,klgfb=klgfb,noverner=noverner, kev=kev
        rad=rad*ab*(ident_wvl#ieq)
        int=int+rad
      ENDIF
    ENDFOR
  ENDELSE
 ;
ENDIF ELSE BEGIN
 ;
  FOR iz=1,30 DO BEGIN
    ab=abund[iz-1]
    IF ab LT min_abund THEN GOTO,lbl1
    FOR ion=1,iz DO BEGIN
      ieq=get_ieq(temp,iz,ion+1,ioneq_logt=ioneq_logt,ioneq_frac=ioneq)
      ip=ionpot[iz-1,ion-1]
      IF total(ieq) NE 0. THEN BEGIN
        freebound_ion,temp,wvl,rad,iz,ion,ip=ip, $
             vdata=vdata,pe=pe,klgfb=klgfb,noverner=noverner, kev=kev
        rad=rad*ab*(ident_wvl#ieq)
        int=int+rad
      ENDIF
    ENDFOR
    lbl1:
  ENDFOR
 ;
ENDELSE

lbl2:

;
; DEM_INT and EM_INT implementation
; ---------------------------------
;  If DEM_INT is specified then it's necessary to multiply by the
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
    print,'%FREEBOUND: Please specify only one of DEM_INT and EM_INT. Returning...'
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
  print,'%FREEBOUND: Warning, number of elements of DEM_INT (EM_INT) must match the'
  print,'  number of temperatures. Returning...'
  int=0.
  return
ENDIF 


;
; If /DT set then need to multiply by dT for each element of EM_VALS.
;
IF keyword_set(dt) THEN BEGIN
  IF nt LT 2 THEN BEGIN
    print,'%FREEBOUND: Two or more temperatures must be specified if DEM_INT is specified.'
    print,'            Returning...'
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
    print,'%FREEBOUND: Error, LOG10 TEMP needs to be evenly spaced if DEM_INT is specified.'
    print,'            E.g., LOG10 TEMP should be spaced at 0.1 dex intervals.'
    int=0.
    return
  ENDIF 
  demfactor=temp*mean_lt/alog10(exp(1))
END ELSE BEGIN
  demfactor=dblarr(nt)+1.
ENDELSE 
;
int=int* (ident_wvl # (em_vals*demfactor))


;
; If /keV specified, then the conversion to photons (from ergs) is different
;
IF keyword_set(photons) THEN BEGIN
  IF keyword_set(kev) THEN BEGIN
    int=int*((12.398/wvl/1.9865d-8)#ident_t)
  ENDIF ELSE BEGIN
    int=int*((wvl/1.9865d-8)#ident_t)
  ENDELSE
ENDIF


siz=size(int)
IF keyword_set(sumt) AND siz[0] EQ 2 THEN int=total(int,2)

t2=systime(1)

IF keyword_set(verbose) THEN $
  print,format='("% FREEBOUND:  continuum calculated in ",f8.1," seconds")',t2-t1


END
