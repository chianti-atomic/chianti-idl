
FUNCTION ch_rad_recomb,gname,temperature, $
                       filename=filename,quiet=quiet,level_resolved=level_resolved

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
; OPTIONAL INPUTS:
;       Filename:  This directly specifies the rad. recombination file
;                  to be read. If specified, then Gname is ignored. If
;                  the file extension is "rrcoeffs", then the file is
;                  read as a level-resolved file.
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
;       IDL> rate=ch_rad_recomb('o_6',t)
;       IDL> rate=ch_rad_recomb('s_2',t,/level,/quiet)
;
; MODIFICATION HISTORY:
;       Ver.1, 25-Oct-2016, Peter Young
;          Code extracted from recomb_rate.pro. 
;       Ver.2, 29 Aug 2024, Roger Dufresne
;          Added option to calculate recombination rate coefficients for metastable levels.
;          In this case it uses .rrcoeffs files, which contain the fitting coefficients.
;          The format of the .rrcoeffs files is different than the .rrparams files.
;       Ver.3, 20-Sep-2024, Peter Young
;          If filename is specified then the routine now checks if the extension is
;          "rrcoeffs" and treats it as a level-resolved file; removed an unnecessary print
;          statement.
;-



IF n_params() LT 2 THEN BEGIN
  print,'Use:  IDL> rate = ch_rad_recomb( Ion_Name, Temperature [, /level_resolved, /quiet, '
  print,'                                 filename= ] )'
  return,-1.
ENDIF

ext=['rrparams','rrcoeffs']

IF n_elements(filename) EQ 0 THEN BEGIN
  convertname,gname,iz,ion
  zion2filename,iz,ion,fname
  if n_elements(level_resolved) gt 0 then rrfile=fname+'.'+ext[1] $
    else rrfile=fname+'.'+ext[0]
ENDIF ELSE BEGIN
 ;
 ; If the specified file has the rrcoeffs extension, then it is
 ; treated as a level resolved file. Otherwise file is treated as
 ; rrparams file, even if /level_resolved has been set.
 ;
  filebasename=file_basename(filename)
  fileparts=filebasename.split('\.')
  IF fileparts[-1] EQ ext[1] THEN BEGIN
    level_resolved=1
    IF NOT keyword_set(quiet) THEN message,/info,/cont,'Specified file will be read as rrcoeffs file.'
  ENDIF ELSE BEGIN
    level_resolved=0
    IF NOT keyword_set(quiet) THEN message,/info,/cont,'Specified file will be read as rrparams file.'
  ENDELSE 
  rrfile=filename
ENDELSE

chck=file_search(rrfile,count=count)
IF count EQ 0 THEN BEGIN
  IF NOT keyword_set(quiet) THEN message,/info,/cont,'The file '+rrfile+' does not exist. Returning...'
  return,-1.
ENDIF 

t=temperature


openr,lur,rrfile,/get_lun

if keyword_set(level_resolved) then begin

  z=1
  ion=1
  nmeta=1
  lvl=1
  wght=1
  rrtype=1
  rstr=''

  readf,lur,rstr
  rrspl=strsplit(rstr,' ',/extract)
  z=fix(rrspl[0])
  ion=fix(rrspl[1])
  nmeta=fix(rrspl[2])

  recomb=dblarr(n_elements(temperature),nmeta)

  for im=0,nmeta-1 do begin
    a=1.
    b=1.
    t0=1.
    t1=1.
    lstr=''
    
    readf,lur,lstr
    rrtype=fix(trim(strmid(lstr,10,5)))

    if rrtype eq 1 then begin
    ; a Badnell type 1
      fmt1a='(3i5,e12.4,f10.5,2e12.4)'

      reads,lstr,lvl,wght,rrtype,a,b,t0,t1,format=fmt1a
      br=b
    endif else if rrtype eq 2 then begin
    ; a Badnell type 2
      c=1.
      t2=1.
      fmt2b='(3i5,e12.4,f10.5,2e11.4,f10.5,e12.4)'

      reads,lstr,lvl,wght,rrtype,a,b,t0,t1,c,t2,format=fmt2b
      br=b+c*exp(-t2/t)
    endif else begin
      free_lun,lur
      print,rrtype
      message,' rrtype not defined'
    ENDELSE

    recomb[*,im]=a/(sqrt(t/t0)*(1.+sqrt(t/t0))^(1.-br)*(1.+sqrt(t/t1))^(1.+br))

  endfor
  free_lun,lur

endif else begin

  rrtype=1
  readf,lur,rrtype

  if rrtype eq 1 then begin
  ; a Badnell type
  ;
    z=1 
    ion=1
    m=1
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
    ion=1
    m=1
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
  endif else BEGIN
    print,rrtype
    print,' rrtype not defined'
    recomb=1.
  ENDELSE

endelse

return,recomb

END
