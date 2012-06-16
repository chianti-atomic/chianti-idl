
function W2R, W

;+
; NAME
;
;     W2R()
;
; EXPLANATION
;
;     Converts a dilution factor into a distance above a star.
;
; INPUTS
;
;     W     Radiation dilution factor.
;
; OUTPUT
;
;     Distance from star -- a number greater than 1.
;
; HISTORY
;
;     Ver.1, 7-Dec-2001, Peter Young
;
; CONTACT
;
;     Peter Young, CfA, pyoung@cfa.harvard.edu
;-

w=double(w)

IF (w GT 0.5) THEN BEGIN
  print,'%R2W: Warning -- W must be =< 0.5. Dilution factor set to 0.5.'
  w=0.5
ENDIF

IF (w LT 0.0) THEN BEGIN
  print,'%R2W: Warning -- W must be >= 0.0. Dilution factor set to 0.'
  return,0.
ENDIF

IF (w EQ 0.0) THEN return,1d50

return,SQRT( 1. / (1-(1-2.*w)^2) )

END
