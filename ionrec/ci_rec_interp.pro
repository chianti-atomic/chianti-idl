
FUNCTION ci_rec_interp, TEMP, RATE_TEMP, RATE, EXTRAP_ABOVE=EXTRAP_ABOVE, $
                        EXTRAP_BELOW=EXTRAP_BELOW

;+
; NAME
;
;   CI_REC_INTERP()
;
; PROJECT
;
;   CHIANTI
;
; EXPLANATION
;
;   For including ionization and recombination into the level balance, it's
;   necessary to interpolate the data stored in the CHIANTI files and, in
;   addition, perform extrapolation to lower or higher temperatures.
;
;   This routine performs the interpolation and extrapolation for general
;   input arrays.
;
; INPUTS
;
;   TEMP   Temperature at which rate required. Units: K. Must be a scalar
;          quantity.
;
;   RATE_TEMP  Temperatures at which RATE is tabulated. Input as Log (base 10)
;              values.
;
;   RATE   The rate coefficient, tabulated at temperatures RATE_TEMP.
;
; KEYWORDS
;
;   EXTRAP_ABOVE Extrapolation to higher temperatures will only take place
;                if this keyword is set.
;
;   EXTRAP_BELOW Extrapolation to lower temperatures will only take place
;                if this keyword is set.
;
; OUTPUT
;
;   The value of the rate coefficient at the temperature TEMP. If no data is
;   available at the specified temperature then a value of zero is returned.
;
; HISTORY
;
;   Ver.1, 29-Jun-2005, Peter Young
;       adapted from original code by Enrico Landi
;-

ltemp=alog10(temp)

;
; Perform extrapolation to higher temperatures
; I've changed to just calculating the linear parameters myself rather than
; using poly & poly_fit as they were generating divide by zero errors for
; some reason.
;
IF ltemp GT max(rate_temp) AND keyword_set(extrap_above) THEN BEGIN
  nrec=n_elements(rate)
  t1=rate_temp[nrec-2:nrec-1] & rat1=alog10(rate[nrec-2:nrec-1])
 ;
  gg=(rat1[1]-rat1[0])/(t1[1]-t1[0])
  cc=rat1[0] - gg*t1[0]
 ;
  out_rate=10.^( gg*ltemp+cc )
  return,out_rate
ENDIF

;
; Perform extrapolation to lower temperatures
;
IF ltemp LT min(rate_temp) AND keyword_set(extrap_below) THEN BEGIN
  t1=rate_temp[0:1] & rat1=alog10(rate[0:1])
 ;
  gg=(rat1[1]-rat1[0])/(t1[1]-t1[0])
  cc=rat1[0] - gg*t1[0]
 ;
  out_rate=10.^( gg*ltemp+cc )
  return,out_rate
ENDIF

;
; Perform interpolation. Note that I use spl_init and spl_interp instead of
; spline as my experience is that they're more reliable.
;
IF ltemp GE min(rate_temp) AND ltemp LE max(rate_temp) THEN BEGIN
  y2=spl_init(rate_temp,alog10(rate))
  out_rate=spl_interp(rate_temp,alog10(rate),y2,alog10(temp))
  return,10.^out_rate
ENDIF

; if we get to this point, then just return zero
;
return,0d0

END
