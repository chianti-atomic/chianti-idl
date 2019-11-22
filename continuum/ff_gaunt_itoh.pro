
FUNCTION ff_gaunt_itoh, big_t, big_z, nonrel=nonrel

;+
; NAME:
;     FF_GAUNT_ITOH
;
; PURPOSE:
;     Returns the frequency-integrated free-free Gaunt
;     factor, as computed by Itoh et al. (2002, A&A, 382, 722). The
;     routine checks the input parameters and returns the relativistic
;     or non-relativistic Gaunt factor accordingly
;
; CATEGORY:
;     CHIANTI; freefree; radiative losses.
;
; CALLING SEQUENCE:
;     Result = FF_GAUNT_ITOH( Big_T, Big_Z )
;
; INPUTS:
;     Big_T:   Temperature in Kelvin (scalar).
;     Big_Z:   Electric charge of ion (scalar).
;
; KEYWORD PARAMETERS:
;     NONREL:  If set, then the routine computes the non-relativistic
;              Gaunt factor using Itoh's fit parameters. 
;
; OUTPUTS:
;     The relativistic, frequency-integrated Gaunt factor. If
;     there's a problem, then a value of -1 is returned.
;
; RESTRICTIONS:
;     The relativistic Gaunt factor formula is only valid for
;     temperatures between logT=6.0 and 8.5, and elements up to
;     Ni. However, for CHIANTI we'll assume it works up to Zn
;     too. Below logT=6.0 we use Itoh's non-relativistic
;     formula which is valid for log gamma^2 from -3 to +2. gamma^2 is
;     Z^2 Ry/kT. If the T,Z pair results in gamma outside of this
;     range, then -1 is returned.
;
; EXAMPLE:
;     IDL> g=ff_gaunt_itoh( 1d7, 20)
;     IDL> g=ff_gaunt_itoh( 1d7, 20, /nonrel)
;
; MODIFICATION HISTORY:
;     Ver.1, 26-Apr-2019, Peter Young
;        Note: checked results against Figure 1 of Itoh et al. (2002).
;     Ver.2, 30-Apr-2019, Peter Young
;        Added non-relativistic formula and confirmed that results are
;        close to relativistic results.
;-


IF n_params() LT 2 THEN BEGIN
  print,'Use:  IDL> g = ff_gaunt_itoh( big_t, big_z [, /nonrel ])'
  return,-1.
ENDIF 

datadir=concat_dir(!xuvtop,'continuum')

;
; The relativistic formula is only valid for temperatures greater
; than 1 MK
;
IF (alog10(big_t) LT 6.0 OR alog10(big_t) GT 8.5) OR keyword_set(nonrel) THEN BEGIN
  datafile=concat_dir(datadir,'itoh_integrated_gaunt_nonrel.txt')
  data=dblarr(11)
  openr,lin,datafile,/get_lun
  readf,lin,data
  free_lun,lin
 ;
  gamma2=(big_z^2*1.579e5)/big_t
  big_gamma=(alog10(gamma2) + 0.5)/2.5
  IF alog10(gamma2) LT -3 OR alog10(gamma2) GT 2 THEN BEGIN
    gaunt=-1.
  ENDIF ELSE BEGIN 
    gaunt=0.0
    FOR i=0,10 DO gaunt=gaunt + data[i]*big_gamma^i
  ENDELSE 
  
ENDIF ELSE BEGIN  ;
 ; Note that data[i,k] corresponds to a_ik in the paper.
 ;
  datafile=concat_dir(datadir,'itoh_integrated_gaunt.txt')
;
  data=dblarr(11,11)
  openr,lin,datafile,/get_lun
  readf,lin,data
  free_lun,lin
 ;
  small_z=(big_z-14.5)/13.5
  small_t=(alog10(big_t)-7.25)/1.25
 ;
  z_arr=dblarr(11,11)
  FOR i=0,10 DO z_arr[*,i]=small_z^i
  t_arr=dblarr(11,11)
  FOR i=0,10 DO t_arr[i,*]=small_t^i
 ;
  gaunt=total(data*z_arr*t_arr)
ENDELSE

return,gaunt

END
