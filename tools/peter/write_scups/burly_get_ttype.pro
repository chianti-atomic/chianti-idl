
FUNCTION burly_get_ttype, gf, de

;+
; NAME:
;    BURLY_GET_TTYPE
;
; PURPOSE:
;    Given the energy and oscillator strength for a transition this
;    routine automatically works out the transition type.
;
; CATEGORY:
;    CHIANTI; data assessment.
;
; CALLING SEQUENCE:
;    Result = BURLY_GET_TTYPE( Gf, De )
;
; INPUTS:
;    Gf:   Weighted oscillator strength.
;    De:   Transition energy (Rydbergs).
;
; OUTPUTS:
;    The transition type, i.e., an integer with one of three values:
;      1  Allowed transition
;      2  Forbidden transition
;      4  Allowed transition with small gf value
;    If an error occurs, then a value of -1 is returned.
;
; EXAMPLE:
;    IDL> ttype=burly_get_ttype(gf,de)
;
; MODIFICATION HISTORY:
;    Ver.1, 2-Jul-2010, Peter Young
;    Ver.2, 22-Jul-2010, Peter Young
;      added parameter check
;    Ver.3, 3-Feb-2017, Peter Young
;      tidied up header.
;-



IF n_params() LT 2 THEN BEGIN
  print,'Use:  IDL> ttype=burly_get_ttype(gf,de)'
  return,-1
ENDIF 

IF gf NE 0 THEN BEGIN
  IF gf GE 1e-3 THEN t_type=1 ELSE t_type=4
ENDIF ELSE BEGIN
  t_type=2
ENDELSE 

return,t_type

END
