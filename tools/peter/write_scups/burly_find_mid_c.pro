FUNCTION burly_find_mid_c, temp, de, t_type

;+
; NAME
;
;    BURLY_FIND_MID_C
;
; PROJECT
;
;    CHIANTI
;
; EXPLANATION
;
;    Works out the scaling parameter that places the middle
;    temperature of the array TEMP at a scaled temperature of 0.5.
;
; INPUTS
;
;    TEMP    Temperature array.
;
;    DE      Energy of transition.
;
;    T_TYPE  CHIANTI transition type.
;
; OUTPUTS
;
;    A scalar giving the scaling parameter that places the middle
;    temperature of TEMP at a scaled temperature value of 0.5.
;
; HISTORY
;
;    Ver.1, 2-Jul-2010, Peter Young
;    Ver.2, 26-May-2013, Peter Young
;       added type 3
;-

IF t_type GT 4 THEN BEGIN
  print,'%BURLY_FIND_MID_C: This routine only works for t_type=1-4.'
  print,'                   Returning value of -1.'
  return,-1
ENDIF 

t=temp
n=n_elements(t)/2
t=t[n]
kt_e=t*8.61735d-5/(de*13.606)

st=0.5

CASE t_type OF
  1: return,sqrt(0.25+kt_e)+0.5
  2: return, kt_e*(1./st -1.)
  3: return, kt_e
  4: return,sqrt(0.25+kt_e)+0.5
ENDCASE

END
