
PRO wgfa_tidy_masterlist, overwrite=overwrite

;+
; NAME:
;      WGFA_TIDY_MASTERLIST
;
; PURPOSE:
;      This runs wgfa_tidy on every ion listed in the CHIANTI
;      masterlist file, and optionally puts the new files back into
;      the database.
;
; CATEGORY:
;      CHIANTI; file formatting.
;
; CALLING SEQUENCE:
;      WGFA_TIDY_MASTERLIST
;
; INPUTS:
;      None.
;	
; KEYWORD PARAMETERS:
;      OVERWRITE: If set, then the new .wgfa files overwrite the files
;                 currently in the database. This keyword should be
;                 used with great care!
;
; OUTPUTS:
;      Writes new versions of the .wgfa files in the current working
;      directory. A text file called 'wgfa_tidy_masterlist_log.txt' is
;      created in the current working directory and summarizes the
;      results of wgfa_tidy for each ion.
;
; MODIFICATION HISTORY:
;      Ver.1, 14-Jan-2019, Peter Young
;-


mlistfile=concat_dir(!xuvtop,'masterlist/masterlist.ions')
read_masterlist,mlistfile,mlist

;
; Remove the "d" ions
;
;; chck=strpos(mlist,'d')
;; k=where(chck LT 0)
;; mlist=mlist[k]

n=n_elements(mlist)

logfile='wgfa_tidy_masterlist_log.txt'
chck=file_info(logfile)
IF chck.exists THEN file_delete,logfile
openw,lout,logfile,/get_lun

FOR i=0,n-1 DO BEGIN
  printf,lout,'Processing '+mlist[i]
 ;
  convertname,mlist[i],iz,ion
  zion2filename,iz,ion,fname
  oldwgfaname='original.wgfa'
  newwgfaname=file_basename(fname+'.wgfa')
  elvlcname=fname+'.elvlc'
  file_copy,fname+'.wgfa',oldwgfaname,/overwrite
  file_copy,elvlcname,'.',/overwrite
 ;
  wgfa_tidy,oldwgfaname,newwgfaname,enfile=elvlcname, /quiet, text_output=text_output
  file_delete,oldwgfaname
  file_delete,file_basename(elvlcname)
 ;
  IF keyword_set(overwrite) THEN BEGIN
    file_copy,newwgfaname,file_dirname(fname),/overwrite
  ENDIF
 ;
  nt=n_elements(text_output)
  FOR j=0,nt-1 DO printf,lout,text_output[j]
ENDFOR

free_lun,lout

END
