PRO freebound, temp, wvl, int, sumt=sumt, photons=photons, $
               noverner=noverner, iz=iz, ion=ion, no_setup=no_setup, $
               min_abund=min_abund, dem_int=dem_int, verbose=verbose, $
               kev=kev, em_int=em_int, abund_file=abund_file, $
               ioneq_file=ioneq_file, element=element, sngl_ion=sngl_ion

;+
; NAME:
;     FREEBOUND
;
; PURPOSE:
;     Calculates the free-bound (radiative recombination)
;     continuum. Please see CHIANTI Technical Report No. 12 for more
;     details about this routine.
;
; CATEGORY:
;     CHIANTI; continuum; free-bound.
;
; CALLING SEQUENCE:
;     FREEBOUND, Temp, Wvl, Int
;
; INPUTS:
;     Temp:   Temperature in K (can be an array). If the input
;             DEM_INT is used then the temperatures must be spaced at
;             equal intervals of 
;             log10(T). E.g., log T=6.0, 6.1, 6.2, 6.3, ...
;     Wvl:    Wavelength in angstroms (can be an array). If /keV is
;             set, then WVL is assumed to contain energies in keV
;             units. In both cases WVL should be monotonically
;             increasing. 
;
; OPTIONAL INPUTS:
;     Dem_Int:  This should be the same size as TEMP and contain
;               differential emission measure values. Specifying
;               DEM_INT places a restriction on TEMP (see above). 
;     Em_Int:   This should be the same size as TEMP and contain
;               emission measure values. The emission measure is
;               N_e*N_H*h, where N_e is electron number density, N_H
;               is the hydrogen number density and h the plasma column
;               depth. Units are cm^-5.
;     Sngl_Ion: A string or string array specifying individual ions
;               for which you want the continuum. Ions should be
;               specified in the CHIANTI format, for example, 'fe_13'
;               for Fe XIII. The name of the *recombined* ion is
;               specified. For example, if you want the hydrogen
;               continuum, use 'h_1' (not 'h_2').
;     Element:  An array specifying which elements are used for the
;               calculation. An element can be specified either
;               with an integer (e.g., 26 for iron) or a string
;               (e.g., 'fe' for iron).
;     Abund_File: The name of a CHIANTI format element abundance
;               file. If not specified then a widget will appear
;               asking you to choose a file. The CHIANTI default file
;               can be specified by giving !abund_file.
;     Ioneq_File: The name of a CHIANTI format ionization equilibrium
;               file. If not specified then the CHIANTI default file
;               (!ioneq_file) will be used.
;     Iz:    Only calculate continuum for the element with atomic 
;            number IZ. **It is recommended that you use the ELEMENT=
;            input instead.** 
;     Ion    (To be used in conjunction with IZ.) Calculate continuum 
;            for a single ion (IZ, ION). **It is recommended that you
;            use the SNGL_ION= input instead.**
;	
; KEYWORD PARAMETERS:
;     MIN_ABUND: If set, calculates the continuum only for those 
;                elements that have an abundance greater than 
;                min_abund. Typically the most abundant elements
;                have abundances > 10^-5.
;     PHOTONS: The output spectrum is given in photon units rather 
;              than ergs.
;     SUMT:    When a set of temperatures is given to FREEBOUND, the 
;              default is to output INTENSITY as an array of size 
;              (nwvl x nT). With this keyword set, a summation over 
;              the temperatures is performed.
;     VERBOSE: Output information from FREEBOUND.
;     KEV:     If set, then WVL is assumed to contain energies in keV units
;              rather than wavelengths in Angstroms.
;     NOVERNER  If set, then the Verner & Yakovlev cross-sections will 
;               not be used. Instead the Karzas & Latter (1961) gaunt
;               factors will be used.
;     NO_SETUP: Does not do anything. Retained only for backwards
;               compatibility. 
;
; OUTPUTS:
;     Int:    Free-bound continuum intensity in units of 
;             10^-40 erg cm^-2/s/sr/Angstrom. That is, it is necessary
;             to multiply by 10^40 to get the specific intensity.
;
;             Note that if neither
;             EM_INT or DEM_INT are specified, then it is assumed that
;             EM_INT=1 for each of the specified temperatures.
;
;             If the keyword /keV is set, then the units of INT will be 
;             10^-40 erg cm^-2/s/sr/keV.
;
; CALLS:
;      FREEBOUND_ION, GET_IEQ, READ_ABUND, READ_IOENQ
;
; EXAMPLE:
;     IDL> wvl=findgen(100)+1.
;     IDL> freebound,1e7,wvl,int
;
; MODIFICATION HISTORY:
;     Ver.1, 24-Jul-2002, Peter Young
;     Ver.2, 26-Jul-2002, Peter Young
;           revised call to freebound_ion; corrected ion fraction problem
;      V 3, 25-May-2005, GDZ 
;                  corrected routine header.
;     Ver.4, 10-Mar-2006, Peter Young
;           added keV keyword
;     Ver.5, 8-Sep-2009, Peter Young
;           routine now checks to make sure temperatures are evenly
;           spaced when the DEM_INT keyword is used; if an error is
;           found due to incorrect temperature specificiations, then
;           the continuum intensity is returned as zero.
;     Ver.6, 6-Oct-2009, Peter Young
;            Introduced EM_INT keyword for specifying emission measure
;            values for isothermal plasmas.
;     Ver.7, 5-Sep-2012, Peter Young
;            Updated header (no change to code).
;     Ver.8, 11-Mar-2020, Peter Young
;        Added ELEMENT=, SNGL_ION=, ABUND_FILE=, IONEQ_FILE= optional
;        inputs; removed common block and disabled /no_setup keyword;
;        tidied up header.
;     Ver.9, 18-Jun-2020, Peter Young
;        Fixed bug in SNGL_ION implementation - needs to be converted
;        to the recombining ion (RECOMB_ION).
;     Ver.10, 23-Jun-2020, Peter Young
;        Updated header to refer to CTR No. 12; no change to code.
;     Ver.11, 12-Feb-2021, Peter Young
;        Removed call to read_klgfb as no longer needed.
;-



IF n_params() LT 3 THEN BEGIN
  print,'Use:  IDL> freebound, temp, wvl, int [/sumt, /photons, /noverner,'
  print,'                      min_abund= , dem_int=, /verbose, abund_file='
  print,'                      /kev, em_int=, element=, sngl_ion=, ioneq_file= ]'
  return
ENDIF



IF n_elements(abund_file) NE 0 THEN BEGIN
  chck=file_search(abund_file,count=count)
  IF count EQ 0 THEN BEGIN
    print,'% FREEBOUND: the specified abund_file was not found. Returning...'
    return
  ENDIF 
ENDIF ELSE BEGIN
  abund_file=''
ENDELSE 
read_abund,abund_file,abund,ab_ref, element=element

;
; Always use the default ioneq file, unless specified by user.
;
IF n_elements(ioneq_file) EQ 0 THEN BEGIN
  ioneq_file=!ioneq_file
ENDIF ELSE BEGIN
  chck=file_search(ioneq_file,count=count)
  IF count EQ 0 THEN BEGIN
    print,'% FREEBOUND: the specified ioneq_file was not found. Returning...'
    return
  ENDIF 
ENDELSE

;
; This converts the input SNGL_ION (specifying the recombined ions) to
; the recombining ions.
;
IF n_elements(sngl_ion) NE 0 THEN BEGIN
  nion=n_elements(sngl_ion)
  recomb_ion=strarr(nion)
  FOR i=0,nion-1 DO BEGIN
    convertname,sngl_ion[i],a,b
    zion2name,a,b+1,rec_ion
    recomb_ion[i]=trim(rec_ion)
  ENDFOR 
ENDIF

read_ioneq,ioneq_file,ioneq_logt,ioneq, ioneq_ref, element=element, sngl_ion=recomb_ion



IF n_elements(ion) NE 0 AND n_elements(iz) EQ 0 THEN BEGIN
  print,'%FREEBOUND:  Error, please specify IZ when using ION'
  return
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

;; read_klgfb,pe,klgfb
;; ksize=size(klgfb)
;; max_ngf=ksize(2)

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
           vdata=vdata,noverner=noverner, kev=kev
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
