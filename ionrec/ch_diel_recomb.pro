
FUNCTION ch_diel_recomb,gname,temperature, $
                        filename=filename, quiet=quiet

;+
; NAME:
;      CH_DIEL_RECOMB
;
; PURPOSE:
;      Read the CHIANTI format dielectronic recombination files and
;      return the rates at the specified temperatures.
;
; CATEGORY:
;      CHIANTI; atomic data; dielectronic recombination.
;
; CALLING SEQUENCE:
;	Result = CH_DIEL_RECOMB( Ion_Name, Temperature )
;
; INPUTS:
;       Gname:  The name of the ion in CHIANTI format, e.g., 'fe_13'
;               for Fe XIII.
;       Temperature:  The temperature(s) in K at which the rates are
;                     required. 
;
; OPTIONAL OUTPUTS:
;       Filename:  This directly specifies the diel. recombination file
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
;       IDL> rate=ch_diel_recomb('o_6',t)
;
; MODIFICATION HISTORY:
;       Ver.1, 25-Oct-2016, Peter Young
;          Code extracted from recomb_rate.pro. 
;       Ver.2 13-Nov-2017, Ken Dere
;          For Badnell data, ions have 9 parameters rather than the 8 previously assumed.
;          code now checks for number of parameters and reads them accordingly
;          pointed out by Will Barnes (Rice U.)
;-



IF n_params() LT 2 THEN BEGIN
  print,'Use:  IDL> rate = ch_diel_recomb( Ion_Name, Temperature [, /quiet, filename= ] )'
  return,-1.
ENDIF 

IF n_elements(filename) EQ 0 THEN BEGIN
  convertname,gname,iz,ion
  zion2filename,iz,ion,fname
  drfile=fname+'.drparams'
ENDIF ELSE BEGIN
  drfile=filename
ENDELSE

chck=file_search(drfile,count=count)
IF count EQ 0 THEN BEGIN
  IF NOT keyword_set(quiet) THEN print,'%CH_DIEL_RECOMB:  The dielectronic recombination file does not exist. Returning...'
  return,-1.
ENDIF 

str=''
t=temperature


z=1
ion=1
;
dr1=fltarr(n_elements(temperature))
openr,lur,drfile,/get_lun
if NOT keyword_set(quiet) then print, ' using drparams file = ',drfile
;
drtype=1
dstr = ''
;
readf,lur,drtype
;
if drtype eq 1 then begin
                        ; a Badnell type
  if NOT keyword_set(quiet) then print, ' drtype = ',drtype, ' a Badnell type 1 file'
  readf,lur,dstr
  dsplit = strsplit(dstr, ' ')
  if n_elements(dsplit) eq 10 then begin
     ncoef=8
     ener=fltarr(ncoef)
     coef=fltarr(ncoef)
     format10 = '(2i5,8e12.4)'
     reads,dstr,z,ion,ener,format = format10
     readf,lur,z,ion,coef,format = format10
  endif  else if n_elements(dsplit) eq 11 then begin
     ncoef=9
     ener=fltarr(ncoef)
     coef=fltarr(ncoef)
     format11 = '(2i5,9e12.4)'
     reads,dstr,z,ion,ener,format=format11
     readf,lur,z,ion,coef,format=format11
  endif
  for icoef=0,ncoef-1 do begin
    dr1(0)=dr1+coef[icoef]*exp(-ener[icoef]/temperature)
  endfor
  dr1=dr1*temperature^(-1.5)
  if NOT keyword_set(quiet)  and keyword_set(dr) then begin
    for icoef=0,ncoef-1 do begin
      print,z,ion,coef[icoef],ener[icoef],format='(2i5,8e12.4)'
    endfor
    print ,' max of dr = ',max(dr1)
  endif
  readf,lur,str
  dr_ref=''
  while not eof(lur) do begin
    readf,lur,str
    dr_ref=dr_ref+str
  endwhile
  free_lun,lur
endif else if drtype eq 2 then begin
  if NOT keyword_set(quiet) then print, ' drtype = ',drtype, ' a Shull type 1 file'
   ;   a Shull type file
  readf,lur,z,ion,adi,bdi,t0,t1,format='(2i5,4e12.4)'
  dr1=adi*exp(-t0/temperature)*(1.+bdi*exp(-t1/temperature))/temperature^1.5
  readf,lur,str
  dr_ref=''
  while not eof(lur) do begin
    readf,lur,str
    dr_ref=dr_ref+str
  endwhile
  free_lun,lur
endif else begin
  if NOT keyword_set(quiet) then print,' drtype not understood'
  dr1=-1.
endelse

return,dr1

END
