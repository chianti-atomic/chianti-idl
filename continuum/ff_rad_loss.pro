
pro ff_rad_loss,t,rad_loss,no_setup=no_setup,min_abund=min_abund, $
                abund_file=abund_file, ioneq_file=ioneq_file, quiet=quiet, $
                element=element, sngl_ion=sngl_ion

;+
; NAME:
;     FF_RAD_LOSS
;
; PURPOSE:
;     Calculate the free-free radiative energy losses.
;
; CATEGORY:
;     CHIANTI; continuum; radiative losses.
;
; CALLING SEQUENCE:
;     FF_RAD_LOSS, T, Rad_Loss
;
; INPUTS:
;     None.
;
; OPTIONAL INPUTS:
;     Abund_File: The name of a CHIANTI format element abundance
;                 file. If not specified, then !abund_file is used.
;     Min_Abund:  A number specifying the minimum abundance to
;                 consider. For example 1e-5 implies only elements
;                 with abundances greater than this will be included
;                 in the calculation.
;     Ioneq_File: The name of a CHIANTI format ionization equilibrium
;                 file. If not set, then !ioneq_file is used.
;     Element:    Either the atomic number of an element, or the
;                 symbol name (e.g., 'Fe' for iron). The loss rate
;                 will only be computed for this element.
;     Sngl_ion:   To calculate the loss spectrum for a single
;                 ion. Ions specified in CHIANTI format (e.g.,
;                 'o_6'). Can be an array of multiple ions. 
;	
; KEYWORD PARAMETERS:
;     QUIET:      If set, then no information is printed to the screen.
;     NO_SETUP:   This keyword is obsolete but retained for backwards
;                 compatibility. 
;
; OUTPUTS:
;     T:          Temperatures in degrees Kelvin. Taken from the
;                 CHIANTI ioneq file.
;     Rad_Loss:   The radiative loss function for free-free
;                 emission. Units: erg s^-1 cm^3.  
;
; EXAMPLE:
;     IDL> ff_rad_loss, t, rad_loss
;     IDL> ff_rad_loss, t, rad_loss, abund_file=''proto_solar_2009_lodders.abund'
;     IDL> ff_rad_loss, t, rad_loss, element='fe'
;     IDL> ff_rad_loss, t, rad_loss, element=26
;     IDL> ff_rad_loss, t, rad_loss, sngl_ion=['fe_8','fe_17','fe_26']
;
; MODIFICATION HISTORY:
;     Ver.1, Apr-2000, Ken Dere
;     Ver.2, 8-Aug-2017, Peter Young
;           added abund_file and ioneq_file optional inputs; added
;           /quiet and check on no_setup
;     Ver.3, 24-Apr-2019, Peter Young
;           removed common block; added element and sngl_ion keywords;
;           updated header.
;     Ver.4, 30-Apr-2019, Peter Young
;           modified to call ff_gaunt_sutherland to get the Gaunt
;           factor rather than compute it internally.
;-


if n_params() lt 2 then begin
   print,' '
   print,' type> ff_rad_loss,temperature,loss_rate [,/no_setup, min_abund= '
   print,'                    abund_file=, ioneq_file=, sngl_ion=, element= ]'
   print,' '
   return
endif
;
kb=1.38062d-16   ;  erg deg-1
ryd=2.17992d-11  ; erg
factor=1.42554d-27

;
; Read element abundances
;
IF n_elements(abund_file) EQ 0 THEN BEGIN
  abund_file=ch_choose_abund()
ENDIF
read_abund,abund_file,abund,abund_ref

;
; Apply min_abund (if set).
;
IF n_elements(min_abund) NE 0 THEN BEGIN
  k=where(abund LE min_abund,nk)
  IF nk NE 0 THEN abund[k]=0.
ENDIF ELSE BEGIN
  min_abund=0.
ENDELSE 

;
; Handle the input 'element'.
;
elt_iz=-1
IF n_elements(element) NE 0 THEN BEGIN
  IF datatype(element) EQ 'STR' THEN BEGIN
    z2element,indgen(30)+1,elt,/symbol,/lower_CASE
    k=where(strlowcase(element) EQ elt,nk)
    IF nk NE 0 THEN elt_iz=k[0]+1
  ENDIF ELSE BEGIN 
    elt_iz=element[0]
  ENDELSE
  ab_save=abund
  abund=abund*0.
  abund[elt_iz-1]=ab_save[elt_iz-1]
ENDIF 


IF n_elements(ioneq_file) EQ 0 THEN ioneq_file=!ioneq_file
read_ioneq,ioneq_file,ioneq_t,ioneq,ioneq_ref
t=10.^ioneq_t
n_ioneq_t=n_elements(ioneq_t)


IF n_elements(sngl_ion) NE 0 THEN BEGIN
  n=n_elements(sngl_ion)
  ioneq_save=ioneq
  ioneq=ioneq*0.
  FOR i=0,n-1 DO BEGIN 
    convertname,sngl_ion[i],iz,ion
    ioneq[*,iz-1,ion-1]=ioneq_save[*,iz-1,ion-1]
  ENDFOR 
  ioneq_save=0.
ENDIF

;
; The following could be speeded up by allowing ff_gaunt_sutherland to
; take array inputs, but ff_rad_loss runs quickly so there's no need.
;
nz=30
z=findgen(nz)+1.
gaunt=fltarr(n_ioneq_t,nz)
FOR i=0,n_ioneq_t-1 DO BEGIN
  FOR j=0,nz-1 DO BEGIN
    gaunt[i,j]=ff_gaunt_sutherland(t[i],z[j])
  ENDFOR
ENDFOR

;
; Now compute the rad_loss contribution from each ion.
; 
zmax=max(where(abund gt 0.))+1
rad_loss=fltarr(n_ioneq_t)
FOR iz=1,zmax DO BEGIN
  this_abund=abund[iz-1]
  IF this_abund GT min_abund THEN BEGIN
    FOR ion=2,iz+1 DO BEGIN
      z=float(ion-1)
      this_ioneq=reform(ioneq[*,iz-1,ion-1])
      k=where(gaunt[*,ion-2] NE -1)
      rad_loss[k]=rad_loss[k]+this_abund*this_ioneq[k]*factor*sqrt(t[k])*z^2*gaunt[k]
    ENDFOR 
  ENDIF
ENDFOR 


END 
