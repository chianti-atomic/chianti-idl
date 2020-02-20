
FUNCTION ff_gaunt_sutherland, big_t, big_z

;+
; NAME:
;     FF_GAUNT_SUTHERLAND
;	ROUTINE_NAME
;
; PURPOSE:
;     Computes the non-relativistic, frequency-averaged free-free
;     Gaunt factor using the approximation of Sutherland (1998, MNRAS,
;     300, 321). Note that there seems to be an error in
;     Sutherland's formulation (see Programming Notes below). 
;
; CATEGORY:
;     CHIANTI; continuum; free-free.
;
; CALLING SEQUENCE:
;     Result = FF_Gaunt_Sutherland( Big_T, Big_Z )
;
; INPUTS:
;     Big_T:   Temperature in K (scalar).
;     Big_Z:   Electric charge of ion (scalar).
;
; OUTPUTS:
;     The relativistic, frequency-integrated Gaunt factor. If
;     there's a problem, then a value of -1 is returned.
;
; PROGRAMMING NOTES:
;     Eq. 18 of the paper refers to gamma_eff^2, which is given as
;     I/kT, where I is the ionization potential of the ion. However,
;     implementing this gives a Gaunt factor that isn't
;     consistent with the non-relativistic one of Itoh et
;     al. (2002). If the expression Z^2 Ry/kT is used instead, then we
;     do get good agreement so I've implemented this instead.
;
; EXAMPLE:
;     IDL> g=ff_gaunt_sutherland( 1e7, 20 )
;
; MODIFICATION HISTORY:
;     Ver.1, 30-Apr-2019, Peter Young
;       Adapted from the code that was previously in the routine
;       ff_rad_loss.pro. 
;-

IF n_params() LT 2 THEN BEGIN
  print,'Use:  IDL> g = ff_gaunt_sutherland( big_t, big_z )'
  return,-1.
ENDIF 

kb=1.38062d-16   ;  erg deg-1
ryd=2.17992d-11  ; erg


read_gffint,g2,gff,s1,s2,s3

;
; This defines log gamma_eff^2, but I'm using Z^2 Ryd/kT rather
; than the expression used by Sutherland.
;
lng2eff=alog10(big_z^2*ryd/(kb*big_t))

IF lng2eff LT -4 OR lng2eff GT 4 THEN BEGIN
  gtot=-1
ENDIF ELSE BEGIN 
  idx=fix((lng2eff+4.)/.2)
  idx=((idx>0)<40)
  delta=lng2eff-g2(idx)
 ;
  gtot=gff(idx)+delta*(s1(idx)+delta*(s2(idx)+delta*s3(idx)))
ENDELSE

return,gtot

END 
