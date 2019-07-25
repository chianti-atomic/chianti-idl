
FUNCTION ch_tot_recomb,gname,temperature, $
                       filename=filename, quiet=quiet

;+
; NAME:
;      CH_TOT_RECOMB
;
; PURPOSE:
;      Read the CHIANTI format total recombination files (i.e., the
;      rates include both radiative and dielectronic recombination)
;      and return the rates at the specified temperatures.
;
; CATEGORY:
;      CHIANTI; atomic data; recombination.
;
; CALLING SEQUENCE:
;	Result = CH_TOT_RECOMB( Ion_Name, Temperature )
;
; INPUTS:
;       Gname:  The name of the ion in CHIANTI format, e.g., 'fe_13'
;               for Fe XIII.
;       Temperature:  The temperature(s) in K at which the rates are
;                     required. 
;
; OPTIONAL OUTPUTS:
;       Filename:  This directly specifies the total recombination file
;                  to be read. If specified, then Gname is ignored. 
;
; KEYWORD PARAMETERS:
;       QUIET:   If set, then no information will be printed to the
;                screen. 
;
; OUTPUTS:
;       The recombination rate in cm^3 s^-1 at the specified
;       temperatures. If a problem occurs, then -1 is returned. 
;
; EXAMPLE:
;       IDL> logt=findgen(41)/10.+4.0
;       IDL> t=10.^t
;       IDL> rate=ch_tot_recomb('o_6',t)
;
; MODIFICATION HISTORY:
;       Ver.1, 25-Oct-2016, Peter Young
;          Code extracted from recomb_rate.pro. 
;-



IF n_params() LT 2 THEN BEGIN
  print,'Use:  IDL> rate = ch_tot_recomb( Ion_Name, Temperature )'
  return,-1.
ENDIF 

IF n_elements(filename) EQ 0 THEN BEGIN
  convertname,gname,iz,ion
  zion2filename,iz,ion,fname
  trfile=fname+'.trparams'
ENDIF ELSE BEGIN
  trfile=filename
ENDELSE

chck=file_search(trfile,count=count)
IF count EQ 0 THEN BEGIN
  IF NOT keyword_set(quiet) THEN print,'%CH_TOT_RECOMB:  The total recombination file does not exist. Returning...'
  return,-1.
ENDIF 

t=temperature

IF NOT keyword_set(quiet) THEN print,'%CH_TOT_RECOMB: found trparams file = ',trfile
openr,lur,trfile,/get_lun
readf,lur,nt,format='(i5)'
totalt=fltarr(nt)
totalr=fltarr(nt)
FOR it=0,nt-1 DO BEGIN
  readf,lur,tt,tr,format='(2e12.4)'
  totalt(it)=tt
  totalr(it)=tr
ENDFOR
free_lun,lur
recomb=exp(spline(alog(totalt),alog(totalr),alog(t)))

return,recomb

END
