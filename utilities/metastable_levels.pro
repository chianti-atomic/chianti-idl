
PRO metastable_levels, ionname, meta, cutoff=cutoff, path=path, quiet=quiet

;+
; NAME:
;    METASTABLE_LEVELS
;
; PURPOSE:
;    Given a CHIANTI .wgfa file this routine works out which levels
;    will be metastable and stores the result in META. A metastable
;    level is defined by whether it has an A-value greater than
;    CUTOFF. If yes, then the level is not metastable; if no then it
;    is metastable. The autoionization rate is included in the check
;    if necessary.
;
; CATEGORY:
;    CHIANTI; utility.
;
; CALLING SEQUENCE:
;    METASTABLE_LEVELS, Ionname, Meta
;
; INPUTS:
;    Ionname:   Ion name in CHIANTI format. E.g., 'fe_13' for Fe XIII.
;
; OPTIONAL INPUTS:
;    Cutoff:   Metastable levels are identified according to whether they
;              have at least one A-value that is > CUTOFF. The default for
;              CUTOFF is 10^5 s^-1.
;    Path:     By default the routine looks in the user's
;              CHIANTI distribution for the .wgfa file. Setting PATH
;              allows you to directly specify a directory containing
;              the .wgfa file.
;	
; KEYWORD PARAMETERS:
;    QUIET:    If set, then no information is printed to the screen.
;
; OUTPUTS:
;    Meta:     A byte array with N elements, where N is the number of
;              levels in the ion model (as judged from the .wgfa file).
;              Metastable levels are identified in this array with the
;              value 1.
;
; CALLS:
;    READ_WGFA_STR, READ_ELVLC, CONVERTNAME, ZION2FILENAME
;
; EXAMPLE:
;    IDL> metastable_levels, 'fe_13', meta
;
; MODIFICATION HISTORY:
;    Ver.1, 4-Feb-2009, Peter Young
;    Ver.2, 13-Feb-2009, Peter Young
;      changed input to ionname; added print out of levels; added
;      /quiet keyword.
;    Ver.3, 26-Jul-2019, Peter Young
;      fixed error when a dielectronic ion is specified; updated
;      header format.
;    Ver.4, 2-Aug-2019, Peter Young
;      moved read_elvlc inside if statement; now check whether the sum
;      of A-value + autoionization rate is greater than cutoff.
;-


IF n_params() LT 1 THEN BEGIN
  print,'Use: IDL> metastable_levels, ionname, meta [, cutoff=, path=, /quiet]'
  return
ENDIF

convertname,ionname,iz,ion,diel=diel
IF n_elements(path) EQ 0 THEN BEGIN 
  zion2filename,iz,ion,filename,diel=diel
  wgfaname=filename+'.wgfa'
  elvlcname=filename+'.elvlc'
ENDIF ELSE BEGIN
  wgfaname=concat_dir(path,ionname+'.wgfa')
  elvlcname=concat_dir(path,ionname+'.elvlc')
ENDELSE 

IF n_elements(cutoff) EQ 0 THEN cutoff=1e5
chck=file_search(wgfaname)
IF chck[0] EQ '' THEN BEGIN
  meta=-1
  print,'%METASTABLE_LEVELS: .wgfa file can not be found. Returning...'
  return
ENDIF
read_wgfa_str,wgfaname,str,ref

n=max(str.lvl2)

meta=bytarr(n)

;
; PRY, 2-Aug-2019: I now sum aval and auto to fix a problem for fe_11d
; whereby levels 11-20 had no A-values, only auto rates.
;
FOR i=1,n DO BEGIN
  k=where(str.lvl2 EQ i)
  IF k[0] EQ -1 THEN BEGIN
    meta[i-1]=1b
  ENDIF ELSE BEGIN
    amax=max(str[k].aval+str[k].auto)
    IF amax LT cutoff THEN meta[i-1]=1b
  ENDELSE
ENDFOR


IF NOT keyword_set(quiet) THEN BEGIN 
  read_elvlc,elvlcname,l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,ref
  k=where(meta EQ 1,nk)
  print,'Metastable levels are: '
  FOR i=0,nk-1 DO BEGIN
    print,format='(i3,a25)',l1[k[i]],term[k[i]]
  ENDFOR
ENDIF 

END
