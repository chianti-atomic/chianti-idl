
FUNCTION r2w, r

;+
; NAME
;
;     R2W()
;
; EXPLANATION
;
;     Converts a distance above a star into a dilution factor.
;
; INPUTS
;
;     R     Distance above a star's surface, measured from the star's 
;           center, in stellar radii units. E.g., R=1 is the surface.
;
; OUTPUT
;
;     Radiation dilution factor - a number between 0 and 0.5
;
; HISTORY
;
;     Ver.1, 8-Aug-2001, Peter Young
;     Ver.2, 12-Nov-2001, Peter Young
;              catches error if r < 1.
;
; CONTACT
;
;     Peter Young, CfA, pyoung@cfa.harvard.edu
;-

r=double(r)

IF r LT 1d0 THEN BEGIN
  print,'%R2W: Warning -- RPHOT must be >= 1. Dilution factor set to zero.'
  return,0d0
ENDIF

return, 0.5d0 * (1d0 - SQRT(1d0 - 1d0/r^2))

END
