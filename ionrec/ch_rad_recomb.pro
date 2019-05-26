
FUNCTION ch_rad_recomb,gname,temperature, $
                       filename=filename, quiet=quiet

;+
; NAME:
;      CH_RAD_RECOMB
;
; PURPOSE:
;      Read the CHIANTI format radiative recombination files and
;      return the rates at the specified temperatures.
;
; CATEGORY:
;      CHIANTI; atomic data; radiative recombination.
;
; CALLING SEQUENCE:
;	Result = CH_RAD_RECOMB( Ion_Name, Temperature )
;
; INPUTS:
;       Gname:  The name of the ion in CHIANTI format, e.g., 'fe_13'
;               for Fe XIII.
;       Temperature:  The temperature(s) in K at which the rates are
;                     required. 
;
; OPTIONAL OUTPUTS:
;       Filename:  This directly specifies the rad. recombination file
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
;       IDL> rate=ch_rad_recomb('o_6',t)
;
; MODIFICATION HISTORY:
;       Ver.1, 25-Oct-2016, Peter Young
;          Code extracted from recomb_rate.pro. 
;-



IF n_params() LT 2 THEN BEGIN
  print,'Use:  IDL> rate = ch_rad_recomb( Ion_Name, Temperature )'
  return,-1.
ENDIF 

IF n_elements(filename) EQ 0 THEN BEGIN
  convertname,gname,iz,ion
  zion2filename,iz,ion,fname
  rrfile=fname+'.rrparams'
ENDIF ELSE BEGIN
  rrfile=filename
ENDELSE

chck=file_search(rrfile,count=count)
IF count EQ 0 THEN BEGIN
  IF NOT keyword_set(quiet) THEN print,'%CH_RAD_RECOMB:  The radiative recombination file does not exist. Returning...'
  return,-1.
ENDIF 

t=temperature


openr,lur,rrfile,/get_lun
if NOT keyword_set(quiet) then print, ' using rrparams file = ',rrfile
rrtype=1
readf,lur,rrtype
;
if rrtype eq 1 then begin
; a Badnell type
;
  z=1 
  n=1
  m=1
  w=1
  a=1.
  b=1.
  t0=1.
  t1=1.
  fmt1a='(3i5,e12.4,f10.5,2e12.4)'
  readf,lur,z,ion,m,a,b,t0,t1,format=fmt1a
  free_lun,lur
  recomb=a/(sqrt(t/t0)*(1.+sqrt(t/t0))^(1.-b)*(1.+sqrt(t/t1))^(1.+b))
;
endif else if rrtype eq 2 then begin
; a badnell type
;
  z=1 
  n=1
  m=1
  w=1
  a=1.
  b=1.
  t0=1.
  t1=1.
  fmt2a='(3i5,e12.4,f10.5,2e11.4,f10.5,e12.4)'
  readf,lur,z,ion,m,a,b,t0,t1,c,t2,format=fmt2a
  b=b+c*exp(-t2/t)
  recomb=a/(sqrt(t/t0)*(1.+sqrt(t/t0))^(1.-b)*(1.+sqrt(t/t1))^(1.+b))
  free_lun,lur
;
endif else if rrtype eq 3 then begin
;  shull type
  z=1
  ion=1
  arad=1.
  xrad=1.
  readf,lur,z,ion,arad,xrad,format='(2i5,2e12.4)'
  recomb=arad/(temperature/1.e+4)^xrad
  free_lun,lur
endif else begin
  print,' rrtype not defined'
  recomb=1.
ENDELSE

return,recomb

END
