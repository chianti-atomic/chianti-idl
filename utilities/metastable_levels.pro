
PRO metastable_levels, ionname, meta, cutoff=cutoff, path=path, quiet=quiet

;+
; NAME
;
;    METASTABLE_LEVELS
;
; PROJECT
;
;    CHIANTI
;
; EXPLANATION
;
;    Given a CHIANTI .wgfa file this routine works out which levels
;    will be metastable and stores the result in META. A metastable
;    level is defined by whether it has an A-value greater than
;    CUTOFF. If yes, then the level is not metastable; if no then it
;    is metastable.
;
; INPUTS
;
;    IONNAME   Ion name in CHIANTI format. E.g., 'fe_13' for Fe XIII.
;
; OPTIONAL INPUTS
;
;    CUTOFF    Metastable levels are identified according to whether they
;              have at least one A-value that is > CUTOFF. The default for
;              CUTOFF is 10^5 s^-1.
;
;    PATH      By default the routine looks in the user's
;              CHIANTI distribution for the .wgfa file. Setting PATH
;              allows you to directly specify a directory containing
;              the .wgfa file.
;
; KEYWORDS
;
;    QUIET     If set, then no information is printed to the screen.
;
; OUTPUTS
;
;    META      A byte array with N elements, where N is the number of
;              levels in the ion model (as judged from the .wgfa file).
;              Metastable levels are identified in this array with the
;              value 1.
;
; CALLS
;
;    READ_WGFA_STR, READ_ELVLC
;
; HISTORY
;
;    Ver.1, 4-Feb-2009, Peter Young
;    Ver.2, 13-Feb-2009, Peter Young
;      changed input to ionname; added print out of levels; added
;      /quiet keyword.
;-

IF n_params() LT 1 THEN BEGIN
  print,'Use: IDL> metastable_levels, ionname, meta [, cutoff=, path=, /quiet]'
  return
ENDIF

convertname,ionname,iz,ion
IF n_elements(path) EQ 0 THEN BEGIN 
  zion2filename,iz,ion,filename
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

FOR i=1,n DO BEGIN
  k=where(str.lvl2 EQ i)
  IF k[0] EQ -1 THEN BEGIN
    meta[i-1]=1b
  ENDIF ELSE BEGIN
    amax=max(str[k].aval)
    IF amax LT cutoff THEN meta[i-1]=1b
  ENDELSE
ENDFOR

read_elvlc,elvlcname,l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,ref

IF NOT keyword_set(quiet) THEN BEGIN 
  k=where(meta EQ 1,nk)
  print,'Metastable levels are: '
  FOR i=0,nk-1 DO BEGIN
    print,format='(i3,a25)',l1[k[i]],term[k[i]]
  ENDFOR
ENDIF 

END
