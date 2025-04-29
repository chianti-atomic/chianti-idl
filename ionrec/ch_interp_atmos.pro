
FUNCTION ch_interp_atmos, param, temp, output_temp, plot=plot, quiet=quiet, $
                          index_good=index_good

;+
; NAME:
;     CH_INTERP_ATMOS
;
; PURPOSE:
;     Interpolate one of the parameters in the CHIANTI model atmosphere
;     file onto a new temperature scale. The interpolation is done in
;     log(temp)-log(param) space.
;
; CATEGORY:
;     CHIANTI; model atmosphere; interpolation.
;
; CALLING SEQUENCE:
;     Result = CH_INTERP_ATMOS( Param, Temp, Output_Temp )
;
; INPUTS:
;     Param:  An array containing the parameter to be interpolated.
;     Temp:   The temperature array on which PARAM is defined (K).
;     Output_Temp: The temperature array for which the parameter values
;                  are needed.
;
; KEYWORD PARAMETERS:
;     QUIET:  If set, then informational messages will not be printed to the
;             IDL window.
;     PLOT:   If set, then a plot showing the interpolation result will be
;             shown.
;
; OUTPUTS:
;     An array of same size as OUTPUT_TEMP containing the interpolated
;     parameter values. If some of the OUTPUT_TEMP values lie outside the
;     range of TEMP, then the parameter values are set to zero at these
;     locations.
;
;     If a problem is found, then a value of -1 is returned.
;
; OPTIONAL OUTPUTS:
;     Index_Good: An integer array containing the indices of OUTPUT_TEMP
;                 that correspond to temperatures within the range of TEMP.
;                 This is used when deriving the ion fraction of he_3 (for
;                 example) from the ion fractions of he_1 and he_2 (see
;                 example below). 
;
; EXAMPLE:
;     IDL> atmos_params=ch_read_atmos()   ; choose a file
;     IDL> ltemp=findgen(41)/10.+4.
;     IDL> h_elec=ch_interp_atmos(atmos_params.h_elec,atmos_params.temp,10.^ltemp)
;
;     IDL> he1_frac=ch_interp_atmos(atmos_params.he1_frac,atmos_params.temp,10.^ltemp)
;     IDL> he2_frac=ch_interp_atmos(atmos_params.he1_frac,atmos_params.temp,10.^ltemp,index_good=k)
;     IDL> he3_frac=he2_frac-he2_frac
;     IDL> he3_frac[k]=1d0-he1_frac[k]-he2_frac[k]
;
; MODIFICATION HISTORY:
;     Ver.1, 23-Apr-2025, Peter Young
;       Adapted from code in ch_adv_model_setup (GDZ,RPD), however I've switched
;       to using spline interpolation (instead of linear interpolation).
;-


IF n_params() LT 3 THEN BEGIN
  print,'Use:  IDL> newparam = ch_interp_atmos( param, temp, output_temp [, /plot, /quiet ] )'
  return,-1
ENDIF 

nt=n_elements(temp)
np=n_elements(param)

IF nt NE np THEN BEGIN
  IF NOT keyword_set(quiet) THEN BEGIN
    message,/info,/cont,'TEMP and PARAM do not have the same numbers of elements. Returning...'
  ENDIF 
  return,-1
ENDIF

n=n_elements(output_temp)
outparam=dblarr(n)

;
; The input temperature needs to be monotonically increasing for
; interpolation to work, so restrict array to lie between min and max
; values.
;
getmin=min(temp,imin)
getmax=max(temp,imax)
;
temp_use=temp[imin:imax]
param_use=param[imin:imax]

;
; Only perform interpolations for temperature region where PARAM is defined.
;
k=where(output_temp LE max(temp) AND output_temp GE min(temp),nk)
out_temp=output_temp[k]
index_good=k

;
; The spline option for interpol requires at least 4 data points. For less
; than this, the default (linear interpolation) is used.
;
nt=n_elements(temp_use)
IF nt GE 4 THEN spline=1 ELSE spline=0
outparam[k]=10.0d0^interpol(alog10(param_use),alog10(temp_use),alog10(out_temp), $
                           spline=spline)

IF keyword_set(plot) THEN BEGIN
  plot,alog10(temp),alog10(param),psym=5,symsize=2, $
       charsize=2, $
       xtitle='Log!d10!n ( Temperature (K) )', $
       ytitle='Log!d10!n ( Parameter )', $
       /ynozero
  
  oplot,alog10(output_temp[k]),alog10(outparam[k]),th=2
ENDIF 

return,outparam

END
