
FUNCTION ion_frac_interp, temp, ioneq_logt, ioneq_frac

;+
; NAME
;
;   ION_FRAC_INTERP
;
; PROJECT
;
;   CHIANTI
;
; EXPLANATION
;
;   Performs interpolation of an ion fraction array (such as contained in
;   the CHIANTI .ioneq files) to yield ion fractions at temperatures TEMP.
;
;   Note that the ion must already have been selected before calling this
;   routine, since IONEQ_FRAC is a 1D array. The routine GET_IEQ is an
;   alternative if the ion has not been selected.
;
; INPUTS
;
;   TEMP        Temperature(s) for which ion fractions are required. 
;               Units: K.
;
;   IONEQ_LOGT  Log (base 10) of temperature (units: K). A 1D array of
;               temperatures for which the ion fraction is defined.
;
;   IONEQ_FRAC  A 1D array of ion fractions at the temperatures defined by
;               ioneq_logt
;
; OUTPUT
;
;   The ion fraction defined at the temperatures TEMP. If values of TEMP lie
;   outside the range of the ion fraction or the range of IONEQ_LOGT, then
;   values of 0 are returned at these temperatures.
;
; HISTORY
;
;   Ver.1, 10-Jun-2005, Peter Young
;-

nt=n_elements(temp)
ltemp=alog10(temp)

answer=dblarr(nt)

i=where(ioneq_frac NE 0.)
IF i[0] EQ -1 THEN return,dblarr(nt)

x=ioneq_logt[i]
y=alog10(ioneq_frac[i])

ind=where(ltemp GE min(x) AND ltemp LE max(x))
IF ind[0] EQ -1 THEN return,dblarr(nt)

xi=ltemp[ind]

y2=spl_init(x,y)
yi=spl_interp(x,y,y2,xi)
answer[ind]=10.^yi

IF n_elements(answer) EQ 1 THEN answer=answer[0]

return,answer

END
