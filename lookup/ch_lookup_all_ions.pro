
PRO ch_lookup_all_ions,ldens_start=ldens_start,ldens_end=ldens_end, execute=execute

;+
; NAME:
;     CH_LOOKUP_ALL_IONS
;
; PURPOSE:
;     Compute population lookup tables for all ions in CHIANTI.
;
; CATEGORY:
;     CHIANTI; level populations.
;
; CALLING SEQUENCE:
;     CH_LOOKUP_ALL_IONS
;
; INPUTS:
;     None.
;
; OPTIONAL INPUTS:
;     Ldens_Start: The start electron density number (cm^-3),
;                  specified as log10. Default is 7.0
;     Ldens_End:   The end electron density number (cm^-3),
;                  specified as log10. Default is 13.0
;	
; KEYWORD PARAMETERS:
;     EXECUTE:  This is required to actually write the lookup
;               tables. If not set, then only an information message
;               is given.
;
; OUTPUTS:
;     If the /EXECUTE keyword is given, then the lookup tables are
;     written to the directory $CHIANTI_LOOKUP. If not given, then
;     only an information message is printed to the screen. Note that
;     the $CHIANTI_LOOKUP variable must be set.
;
; RESTRICTIONS:
;     Must define the environment variable $CHIANTI_LOOKUP.
;
; EXAMPLE:
;     IDL> ch_lookup_all_ions
;     IDL> ch_lookup_all_ions, /execute
;     IDL> ch_lookup_all_ions, ldens_start=1.0, ldens_end=10.0
;
; MODIFICATION HISTORY:
;     Ver.1, 26-Jul-2019, Peter Young
;-

outdir=getenv('CHIANTI_LOOKUP')
IF outdir EQ '' THEN BEGIN
  print,'% CH_LOOKUP_ALL_IONS: please define the environment variable $CHIANTI_LOOKUP to point to where the lookup files should be written.'
  print,'                      Returning...'
  return
ENDIF 

mlistdir=concat_dir(!xuvtop,'masterlist')
mlistfile=concat_dir(mlistdir,'masterlist.ions')
read_masterlist,mlistfile,mlist

IF n_elements(ldens_start) EQ 0 THEN ldens_start=7.0
IF n_elements(ldens_end) EQ 0 THEN ldens_end=13.0

n=n_elements(mlist)

IF NOT keyword_set(execute) THEN BEGIN
  print,''
  print,'* This routine will create lookup tables for all ions in CHIANTI and write them to the directory'
  print,'         ',expand_path(outdir)
  print,''
  print,'* The density range for the calculations will be:'
  print,format='(9x,"ldens_start=",f8.2,"  ldens_end=",f8.2)',ldens_start,ldens_end
  print,''
  print,'* The number of ion tables to be generated is '+trim(n)+'.'
  print,''
  print,'* To perform the calculation, please run this routine again but with the /EXECUTE keyword set.'
  print,''
  return
ENDIF 

FOR i=0,n-1 DO BEGIN
  print,'Writing lookup table for '+trim(mlist[i])+' (ion '+trim(i+1)+' of '+trim(n)+')...'
  ch_write_pop_lookup_table,mlist[i],dir_lookup=outdir, $
                            ldens_start=ldens_start,ldens_end=ldens_end
ENDFOR

END
