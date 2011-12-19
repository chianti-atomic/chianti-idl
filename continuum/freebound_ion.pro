
PRO freebound_ion, temp, wvl, int, iz, ion, ip=ip, vdata=vdata, pe=pe, $
                   klgfb=klgfb, noverner=noverner, kev=kev

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
; EXPLANATION
;
;      Calculates the free-bound (radiative recombination) continuum 
;      from a single ion. Note that the output does not contain the ion 
;      fraction, element abundance or differential emission measure term.
;
; INPUTS
;
;      TEMP    Temperature in K (can be an array).
;
;      WVL     Wavelength in angstroms (can be an array). If the keyword
;              /keV is set, then wvl should instead be given in units of
;              keV.
;
;      IZ      Atomic number of ion (e.g., 26 = Fe)
;
;      ION     Spectroscopic number of ion (e.g., 13 = XIII)
;
; OUTPUTS
;
;      INT     Free-bound continuum intensity. Needs to be multiplied by 
;              element abundance, ion fraction and DEM to obtain the final 
;              continuum intensity.
;
; OPTIONAL INPUTS
;
;      IP      The ionization potential of the ion.
;
;      VDATA   An array containing the Verner & Yakovlev data array.
;
;      PE      An array containing the PE data from READ_FBLVL
;
;      KLGFB   An array containing the KLGFB data from READ_FBLVL
;
;      [Note: the above 3 inputs are used when calling freebound_ion from 
;       freebound]
;
; KEYWORDS
;
;      NOVERNER  If set, then the Verner & Yakovlev cross-sections will 
;                not be used.
;
;      KEV      If set, then WVL is assumed to contain energies in keV units
;               rather than wavelengths in Angstroms.
;
; COMMON BLOCKS
;
;      None.
;
; CALLS
;
;      READ_FBLVL, ZION2FILENAME, VERNER_XS, KARZAS_XS, CONCAT_DIR,
;      READ_IP, READ_KLGFB, FILE_EXIST
;
; PROGRAMMING NOTES
;
;      The way I treat the exponential function in the expression for the 
;      emissivity may seem strange, but it saves a bit of time in the 
;      calculation. Basically calculating exp(E-IP+E_i/T) for each level 
;      was time consuming so I split it into exp(E-IP/T)*exp(E_i/T). The 
;      first term comes out of the for loop, while the second term is the 
;      exponential of a vector rather than an array. This made a ~30% 
;      time-saving for freebound.
;
; HISTORY
;
;      Ver.1, 24-Jul-2002, Peter Young
;
;      Ver.2, 26-Jul-2002, Peter Young
;          Added /noverner keyword.
;
;      Ver.3, 30-Jul-2002, Peter Young
;          Speeded up routine by modifying treatment of exponential.
;
;      Ver.4, 10-Mar-2006, Peter Young
;          Added /keV keyword
;
;       Ver. 5 8-Jun-2010, Ken Dere
;           added check for n_params()
;
;-
if n_params() lt 5 then begin
    print,' > freebound_ion, temp, wvl, int, iz, ion, ip=ip, vdata=vdata, pe=pe, '
    print,'        klgfb=klgfb, noverner=noverner, kev=kev'
    print,'  '
    return
endif

const=3.0992d-46

nt=n_elements(temp)
nwvl=n_elements(wvl)
rad=dblarr(nwvl,nt)
int=rad
ident_wvl=make_array(nwvl,val=1,/double)

;
; * keV units *
;
IF keyword_set(kev) THEN BEGIN
  wvl_save=wvl
  efact=12.398/wvl^2   ; hc/E^2 factor to multiply continuum by
  efact=reverse(efact)
  wvl=12.398/wvl   ; convert energies to angstroms
  wvl=reverse(wvl)
ENDIF ELSE BEGIN
  efact=1d0
ENDELSE


IF (n_elements(vdata) EQ 0) AND NOT keyword_set(noverner) THEN BEGIN
  vdata=dblarr(10,465)
  dir=concat_dir(!xuvtop,'continuum')
  fname=concat_dir(dir,'verner_short.txt')
  openr,lun,fname,/get_lun
  readf,lun,vdata
  free_lun,lun
ENDIF

IF n_elements(pe) EQ 0 AND n_elements(klgfb) EQ 0 THEN read_klgfb,pe,klgfb

IF n_elements(ip) EQ 0 THEN BEGIN
  dir=concat_dir(!xuvtop,'ip')
  fname=concat_dir(dir,'chianti.ip')
  read_ip,fname,ionpot,ref
  ip=ionpot[iz-1,ion-1]
ENDIF

ewvl=1d8/wvl

;
; read .fblvl file for recombined ion
;
zion2filename,iz,ion,filename
ename=filename+'.fblvl'
IF file_exist(ename) THEN BEGIN
  read_fblvl,ename,l1,conf,pqn,ll,spd,mult,ecm,ecmth,eref
  mult=float(mult)
  ind=where(ecm EQ 0. AND l1 NE 1)
  IF ind[0] NE -1 THEN ecm[ind]=ecmth[ind]
  wecm=1.e+8/(ip-ecm)
  necm=n_elements(ecm)
ENDIF ELSE BEGIN
  IF keyword_set(kev) THEN wvl=wvl_save
  return
ENDELSE


;
; if the maximum energy in WVL is not enough to ionize the ion, then don't
; bother calculating emissivity
;
IF max(ewvl) LT min(ip-ecm) THEN BEGIN
  IF keyword_set(kev) THEN wvl=wvl_save
  return
ENDIF


;
; read .fblvl file for recombining ion
;
zion2filename,iz,ion+1,filenamer
enamer=filenamer+'.fblvl'
IF file_exist(enamer) THEN BEGIN
  read_fblvl,enamer,l1r,confr,pqnr,llr,spdr,multr,ecmr,erefr
  multr=float(multr)
ENDIF ELSE BEGIN
  multr=1.
ENDELSE

ind=where( (wecm GE 0.) AND (max(ewvl) GE ecm) )
ni=n_elements(ind)
IF ind[0] EQ -1 THEN BEGIN
  IF keyword_set(kev) THEN wvl=wvl_save
  return
ENDIF

FOR j=0,ni-1 DO BEGIN
  i=ind[j]

  IF (i EQ 0) AND NOT keyword_set(noverner) THEN BEGIN
    xs=verner_xs(iz,ion,wvl,data=vdata)
  ENDIF ELSE BEGIN
    xs=karzas_xs(wvl,pqn[i],ll[i],ip-ecm[i],pe=pe,klgfb=klgfb)
  ENDELSE

  expf=ident_wvl # exp(-1.43873 * ecm[i] / temp)
  rad_i=const * mult[i]/multr[0] * ((ewvl^5*xs) # (1d0/temp^(1.5))) * expf
  rad=rad+rad_i

ENDFOR

ind=where(rad NE 0.)
IF ind[0] NE -1 THEN BEGIN
  lrad=alog(rad[ind])
  f=-1.43873 * (ewvl-ip) # (1d0/temp)
  int[ind]=exp(lrad+f[ind])*1d40/4d0/!pi
ENDIF

;
; multiply int by the hc/E^2 factor if /kev set
;
IF keyword_set(kev) THEN BEGIN
  e_arr=efact#(dblarr(nt)+1d0)
  int=int*e_arr
 ;
 ; reverse the wavelength dimension to match energy ordering rather than
 ; wavelength ordering
  int=int[reverse(indgen(nwvl)),*]
 ;
 ; and set wavelengths back to energy units
  wvl=wvl_save
ENDIF

END
