
FUNCTION ci_rec_interp, TEMP, RATE_TEMP, RATE, EXTRAP_ABOVE=EXTRAP_ABOVE, $
                        EXTRAP_BELOW=EXTRAP_BELOW

;+
; NAME:
;   CI_REC_INTERP()
;
; PURPOSE:
;   For including ionization and recombination into the level balance, it's
;   necessary to interpolate the data stored in the CHIANTI files and, in
;   addition, perform extrapolation to lower or higher temperatures.
;
;   This routine performs the interpolation and extrapolation for general
;   input arrays.
;
; CATEGORY:
;   CHIANTI; ionization; recombination; interpolation.
;
; CALLING SEQUENCE:
;   Result = CI_REC_INTERP( Temp, Rate_Temp, Rate )
;
; INPUTS:
;   Temp:  Temperature at which rate required. Units: K. Must be a scalar
;          quantity.
;   Rate_Temp:  Temperatures at which RATE is tabulated. Input as Log (base 10)
;               values.
;   Rate:  The rate coefficient, tabulated at temperatures RATE_TEMP.
;	
; KEYWORD PARAMETERS:
;   EXTRAP_ABOVE Extrapolation to higher temperatures will only take place
;                if this keyword is set. Note that extrapolation is
;                performed in logT-logRate space. If the extrapolated
;                points have higher rates than the last tabulated
;                point in the data file, then they will be replaced by
;                the last tabulated point. As of 2018, the
;                /extrap_above option is only used for the
;                recombination rates.
;
;   EXTRAP_BELOW Extrapolation to lower temperatures will only take place
;                if this keyword is set. As of 2018, the
;                /extrap_below option is only used for the
;                ionization rates.
;
; OUTPUTS:
;   The value of the rate coefficient at the temperature TEMP. If no data is
;   available at the specified temperature then a value of zero is returned.
;
; EXAMPLE:
;   Perform interpolation for a Fe XVIII transition, with
;   extrapolation used for high temperature points.
;     IDL> zion2filename,26,18,fname
;     IDL> read_ionrec,fname,rec_rate,ci_rate,temp_ionrec,luprec,lupci,status,rec_ref,ci_ref
;     IDL> temp=10.^(findgen(21)/10.+6.0)
;     IDL> rate=ci_rec_interp(temp,temp_ionrec,reform(rec_rate[50,*]),/extrap_above)
;
; MODIFICATION HISTORY:
;   Ver.1, 29-Jun-2005, Peter Young
;       Adapted from original code by Enrico Landi
;   Ver.2, 5-Nov-2018, Peter Young
;       Modified behavior of /EXTRAP_ABOVE so that the extrapolated
;       value can not be higher than the last point in the rate
;       tabulation.  This is a safety net in case there are anomalously
;       high rates in the data tabulation. Also updated header, and
;       added check on inputs.
;-


IF n_params() LT 3 THEN BEGIN
  print,'Use:  IDL> output=ci_rec_interp(temp,rate_temp,rate [,/extrap_above,/extrap_below])'
  return,0.
ENDIF


;
; Make sure the input temperature is a scalar.
;
nt=n_elements(temp)
IF nt GT 1 THEN BEGIN
  print,'% CI_REC_INTERP: the input TEMP must be a scalar. Returning...'
  return,0.
ENDIF 

ltemp=alog10(temp)

;
; Perform extrapolation to higher temperatures
; I've changed to just calculating the linear parameters myself rather than
; using poly & poly_fit as they were generating divide by zero errors for
; some reason.
;
; PRY, 5-Nov-2018: Now force the output rate to be no larger than the
; last point in the data tabulation.
;
IF ltemp GT max(rate_temp) AND keyword_set(extrap_above) THEN BEGIN
  nrec=n_elements(rate)
  t1=rate_temp[nrec-2:nrec-1] & rat1=alog10(rate[nrec-2:nrec-1])
 ;
  gg=(rat1[1]-rat1[0])/(t1[1]-t1[0])
  cc=rat1[0] - gg*t1[0]
 ;
  out_rate=10.^( gg*ltemp+cc )
  out_rate=min([out_rate,rate[nrec-1]])  ; 5-Nov-2018
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
