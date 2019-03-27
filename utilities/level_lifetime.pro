
PRO level_lifetime, ionname, level, lifetime, path=path, quiet=quiet

;+
; NAME
;
;      LEVEL_LIFETIME
;
; PROJECT
;
;      CHIANTI
;
; EXPLANATION
;
;      Calculates the lifetime of an ion's level using atomic
;      data in the CHIANTI atomic database.
;
; INPUTS
;
;      IONNAME  Name of the ion in CHIANTI format. E.g., 'fe_13' for
;               Fe XIII.
;
;      LEVEL    CHIANTI level index. Check the CHIANTI .elvlc file to
;               find this.
;
; OPTIONAL INPUTS
;
;      PATH     By default the routine uses the data files in
;               the user's CHIANTI distribution. By setting
;               PATH the user can directly the specify the directory
;               containing the atomic data files.
;
; KEYWORDS
;
;      QUIET    If set, then nothing is printed to the screen.
;
; OUTPUTS
;
;      LIFETIME The level's lifetime in seconds.
;
; EXAMPLE
;
;      IDL> level_lifetime,'fe_13',3
;
;      Level:     3s2.3p2 3P2
;      Lifetime (seconds):  1.027e-01
;
; CALLS
;
;      CONVERTNAME, ZION2FILENAME, READ_WGFA_STR, READ_ELVLC
;
; HISTORY
;
;      Ver.1, 13-Feb-2009, Peter Young
;-

IF n_params() LT 2 THEN BEGIN
  print,'Use  IDL> level_lifetime, ionname, level, lifetime [, path=, /quiet]'
  return
ENDIF 

convertname,ionname,iz,ion
IF n_elements(path) EQ 0 THEN BEGIN
  zion2filename,iz,ion,fname
  wgfaname=fname+'.wgfa'
  elvlcname=fname+'.elvlc'
ENDIF ELSE BEGIN
  wgfaname=concat_dir(path,ionname+'.wgfa')
  elvlcname=concat_dir(path,ionname+'.elvlc')
ENDELSE 

read_wgfa_str,wgfaname,struc

read_elvlc,elvlcname,l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,ref

i=where(struc.lvl2 EQ level,ni)

IF ni EQ 0 THEN BEGIN
  print,'%LEVEL_LIFETIME:  There are no radiative decays from this level!  Returning...'
  return
ENDIF 

aval=total(struc[i].aval)

lifetime=1./aval

IF NOT keyword_set(quiet) THEN BEGIN 
  print,''
  print,'Level:     ',term[level-1]
  print,format='("Lifetime (seconds): ",e10.3)',lifetime
  print,''
ENDIF 

END
