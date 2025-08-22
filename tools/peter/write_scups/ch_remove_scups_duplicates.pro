
PRO ch_remove_scups_duplicates, file

;+
; NAME:
;     CH_REMOVE_SCUPS_DUPLICATES
;
; PURPOSE:
;     Removes duplicate transitions from a CHIANTI scups file. It is
;     assumed that the last occurrence of the transition in the file is
;     the one to be kept.
;
; CATEGORY:
;     CHIANTI; scups; validation.
;
; CALLING SEQUENCE:
;     CH_REMOVE_SCUPS_DUPLICATES, File
;
; INPUTS:
;     File:  The name of a CHIANTI scups file.
;
; OUTPUTS:
;     Writes the file [file]_remove_duplicates to the working directory.
;     This is identical to the input file except duplicate transitions
;     are removed. The last occurrence of the duplicated transition is
;     retained.
;
; EXAMPLE:
;     IDL> ch_remove_scups_duplicates, 'ar_3.scups'
;
; MODIFICATION HISTORY:
;     Ver.1, 22-Aug-2025, Peter Young
;-

IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> ch_remove_scups_duplicates, scups_file'
  return
ENDIF 

chck=file_info(file)
IF chck.exists EQ 0 THEN BEGIN
  message,/cont,/info,'The specified file does not exist! Returning...'
  return
ENDIF 

read_scups,file,str
d=str.data
n=n_elements(d)

s=ch_duplicate_transitions(file,status=status)

;
; The code below for writing out the transition data is copied from
; write_scups.pro.
;
IF status EQ 1 THEN BEGIN
  outfile=file+'_remove_duplicates'
  openw,lout,outfile,/get_lun
  ;
  missing=str.info.missing
  ;
  FOR i=0,n-1 DO BEGIN
    k=where(d.lvl1 EQ d[i].lvl1 AND d.lvl2 EQ d[i].lvl2,nk)
    IF nk EQ 1 OR (nk GT 1 AND i EQ max(k)) THEN BEGIN
      ;
      chck=where(d[i].spl NE missing,nt)
      IF nt EQ 1 THEN stop
      ;
      ;  Write out first line of data
      ;
      c_ups=d[i].c_ups
      IF c_ups GE 1e4 OR c_ups LT 0.1 THEN cform='e12.4' ELSE cform='f12.5'
      IF d[i].lim GE 0 THEN limform='(e12.3)' ELSE limform='(i12)'
      format_str='(2i7,2e12.3,'+limform+',2i5,'+cform+')'
      
      printf,lout,format=format_str, $
             d[i].lvl1,d[i].lvl2,d[i].de,d[i].gf,d[i].lim, $
             nt,d[i].t_type,c_ups
      ;
      ;  Write out second line of data (scaled temperatures)
      ;
      format_str='('+trim(nt)+'e12.3)'
      printf,lout,format=format_str,d[i].stemp[0:nt-1]
      ;
      ;  Write out third line of data (scaled upsilons)
      ;
      printf,lout,format=format_str,d[i].spl[0:nt-1]
    ENDIF
  ENDFOR
  printf,lout,' -1'
  nc=n_elements(str.info.comments)
  FOR i=0,nc-1 DO printf,lout,str.info.comments[i]
  printf,lout,'File processed with ch_remove_scups_duplicates, '+systime()+'.'
  printf,lout,' -1'
  free_lun,lout
  message,/info,/cont,'The updated scups file is: '+outfile
ENDIF ELSE BEGIN
  message,/info,/cont,'No duplicates found in this file.'
ENDELSE 

END
