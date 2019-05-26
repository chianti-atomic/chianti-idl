
PRO write_splstr, splstr, filename=filename, overwrite=overwrite

;+
; NAME
;
;      WRITE_SPLSTR
;
; PROJECT
;
;      CHIANTI
;
; EXPLANATION
;
;      Given a CHIANTI spline structure, this routine writes out the
;      SCUPS format file.
;
; INPUTS
;
;      SPLSTR   A structure containing data from the CHIANTI SCUPS
;               format files. The format of the structure should be
;               the same as that produced by READ_SCUPS.
;
; OPTIONAL INPUTS
;
;      FILENAME  The name of the output filename. If not specified,
;                then the file ION_NAME.scups is written, where
;                ION_NAME is the CHIANTI format for the ion name
;                (e.g., 'o_6' for O VI).
;
; KEYWORDS
;
;      OVERWRITE  If the SCUPS file already exists, then it will not
;                 be over-written unless the /OVERWRITE keyword is
;                 set.
;
; OUTPUTS
;
;      A file is written to the current working directory.
;
; HISTORY
;
;      Ver.1, 28-May-2013, Peter Young
;-

IF n_elements(filename) EQ 0 THEN BEGIN
  ion_name=splstr.info.ion_name
  filename=ion_name+'.scups'
ENDIF 

chck=file_search(filename)
IF chck[0] NE '' AND NOT keyword_set(overwrite) THEN BEGIN
  print,'%WRITE_SPLSTR: the file '+trim(filename)+' already exists. Please use the keyword /OVERWRITE to over-write this file.'
  return
ENDIF 

openw,lout,filename,/get_lun

nt=splstr.info.ntrans
FOR i=0,nt-1 DO BEGIN 
;
;  Write out first line of data
;
  IF splstr.data[i].lim GE 0 THEN format_str='(2i7,3e12.3,2i5,e12.3)' ELSE format_str='(2i7,2e12.3,i12,2i5,e12.3)'
  printf,lout,format=format_str, $
         splstr.data[i].lvl1, $
         splstr.data[i].lvl2, $
         splstr.data[i].de, $
         splstr.data[i].gf, $
         splstr.data[i].lim, $
         splstr.data[i].nspl, $
         splstr.data[i].t_type, $
         splstr.data[i].c_ups 
;
  nspl=splstr.data[i].nspl
  stemp=splstr.data[i].stemp
  spl=splstr.data[i].spl
;
;  Write out second line of data (scaled temperatures)
;
  format_str='('+trim(nspl)+'e12.3)'
  printf,lout,format=format_str,stemp[0:nspl-1]
;
;  Write out third line of data (scaled upsilons)
;
  printf,lout,format=format_str,spl[0:nspl-1]
ENDFOR

;
; Write out comments.
;
printf,lout,' -1'
comments=splstr.info.comments
n=n_elements(comments)
FOR i=0,n-1 DO printf,lout,comments[i]

free_lun,lout

END
