
FUNCTION ch_diel_recomb,gname,temperature, $
            filename=filename,quiet=quiet,level_resolved=level_resolved

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
; OPTIONAL INPUTS:
;       Filename:  This directly specifies the diel. recombination file
;                  to be read. If specified, then Gname is ignored. 
;
; KEYWORD PARAMETERS:
;       QUIET:   If set, then no information will be printed to the
;                screen. 
;       LEVEL_RESOLVED:  Calculates rate coefficients for ground and 
;                        metastable levels.
;
; OUTPUTS:
;       The recombination rate in cm^3 s^-1 at the specified
;       temperatures. If a problem occurs, then -1 is returned. 
;
; EXAMPLE:
;       IDL> logt=findgen(41)/10.+4.0
;       IDL> t=10.^t
;       IDL> rate=ch_diel_recomb('o_6',t)
;       IDL> rate=ch_diel_recomb('s_2',t,/level,/quiet)
;
; MODIFICATION HISTORY:
;       Ver.1, 25-Oct-2016, Peter Young
;          Code extracted from recomb_rate.pro. 
;       Ver.2 13-Nov-2017, Ken Dere
;          For Badnell data, ions have 9 parameters rather than the 8 previously assumed.
;          code now checks for number of parameters and reads them accordingly
;          pointed out by Will Barnes (Rice U.)
;       Ver.3, 29 Aug 2024, Roger Dufresne
;          Added option to calculate recombination rate coefficients for metastable levels.
;          In this case it uses .drcoeffs files, which contain the fitting coefficients.
;          The format of the .drcoeffs file is different than the .drparams files.
;-



IF n_params() LT 2 THEN BEGIN
  print,'Use:  IDL> rate = ch_diel_recomb( Ion_Name, Temperature [, /quiet, filename= ] )'
  return,-1.
ENDIF 

IF n_elements(filename) EQ 0 THEN BEGIN
  convertname,gname,iz,ion
  zion2filename,iz,ion,fname
  if keyword_set(level_resolved) then drfile=fname+'.drcoeffs' $
    else drfile=fname+'.drparams'
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
openr,lur,drfile,/get_lun
if NOT keyword_set(quiet) then print, ' using drcoeffs file = ',drfile
;
drtype=1
dstr = ''
;
if keyword_set(level_resolved) then begin
  lvl=1
  wght=1
  readf,lur,dstr
  drspl=strsplit(dstr,' ',/extract)

  z=fix(drspl[0])
  ion=fix(drspl[1])
  ; read in number of levels
  nmeta=fix(drspl[2])
  dr1=dblarr(n_elements(temperature),nmeta)

  ; all level resolved files are currently a Badnell type 1 file
  if NOT keyword_set(quiet) then print, ' drtype = ',drtype, ' a Badnell type 1 file'
  
  for im=0,nmeta-1 do begin
    ncoef=9
    ener=fltarr(ncoef)
    coef=fltarr(ncoef)
    format12='(3i5,9e12.4)'
    readf,lur,lvl,wght,drtype,ener,format=format12
    readf,lur,lvl,wght,drtype,coef,format=format12

    for icoef=0,ncoef-1 do begin
      dr1[*,im]=dr1[*,im]+coef[icoef]*exp(-ener[icoef]/temperature)
    endfor
    dr1[*,im]=dr1[*,im]*temperature^(-1.5)

    if NOT keyword_set(quiet) and keyword_set(dr) then begin
      for icoef=0,ncoef-1 do begin
      	print,z,ion,im+1,coef[icoef],ener[icoef],format='(3i5,8e12.4)'
      endfor
      print ,' max of dr for level '+(im+1)+' = ',max(dr1[*,im])
    endif
  endfor
  
  readf,lur,str
  dr_ref=''
  while not eof(lur) do begin
    readf,lur,str
    dr_ref=dr_ref+str
  endwhile
  free_lun,lur

endif else begin
  dr1=fltarr(n_elements(temperature))
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

endelse

return,dr1

END
