
PRO write_scups_struc, scupstr, filename=filename, overwrite=overwrite

;+
; NAME:
;      WRITE_SCUPS_STRUC
;
; CATEGORY:
;      CHIANTI; file input/output.
;
; PURPOSE:
;      For a structure containing scaled upsilon data in the CHIANTI
;      format (see CHIANTI Technical Report No. 13), this routine the
;      data into the CHIANTI SCUPS file format.
;
; INPUTS:
;      Scupstr: A structure containing data from the CHIANTI SCUPS
;               format files. The format of the structure should be
;               the same as that produced by READ_SCUPS.
;
; OPTIONAL INPUTS:
;      Filename: The name of the output filename. If not specified,
;                then the file ION_NAME.scups is written, where
;                ION_NAME is the CHIANTI format for the ion name
;                (e.g., 'o_6' for O VI).
;
; KEYWORDS PARAMETERS:
;      OVERWRITE: If the SCUPS file already exists, then it will not
;                 be over-written unless the /OVERWRITE keyword is
;                 set.
;
; OUTPUTS:
;      A file is written to the current working directory.
;
; MODIFICATION HISTORY:
;      Ver.1, 28-May-2013, Peter Young
;      Ver.2, 27-Jan-2017, Peter Young
;         Renamed to WRITE_SCUPS_STR; now writes an additional -1 at
;         the end of the file.
;      Ver.3, 3-Feb-2017, Peter Young
;         Adjusted output format.
;      Ver.4, 24-Jun-2020, Peter Young
;         Now checks that a comment string isn't empty before printing.
;-

IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> write_scups_str, scupstr [, filename=, /overwrite]'
  return
ENDIF 

IF n_elements(filename) EQ 0 THEN BEGIN
  ion_name=scupstr.info.ion_name
  filename=ion_name+'.scups'
ENDIF 

chck=file_search(filename)
IF chck[0] NE '' AND NOT keyword_set(overwrite) THEN BEGIN
  print,'%WRITE_SCUPS_STRUC: the file '+trim(filename)+' already exists. Please use the keyword /OVERWRITE to over-write this file.'
  return
ENDIF 

openw,lout,filename,/get_lun

nt=scupstr.info.ntrans
FOR i=0,nt-1 DO BEGIN 
;
;  Write out first line of data
;
  IF scupstr.data[i].c_ups GE 1e4 OR scupstr.data[i].c_ups LT 0.1 THEN cform='e12.4' ELSE cform='f12.5'
    IF scupstr.data[i].lim GE 0 THEN limform='(e12.3)' ELSE limform='(i12)'
  format_str='(2i7,2e12.3,'+limform+',2i5,'+cform+')'

  printf,lout,format=format_str, $
         scupstr.data[i].lvl1, $
         scupstr.data[i].lvl2, $
         scupstr.data[i].de, $
         scupstr.data[i].gf, $
         scupstr.data[i].lim, $
         scupstr.data[i].nspl, $
         scupstr.data[i].t_type, $
         scupstr.data[i].c_ups 
;
  nspl=scupstr.data[i].nspl
  stemp=scupstr.data[i].stemp
  spl=scupstr.data[i].spl
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
; Write out comments. Don't print out any empty lines, though.
;
printf,lout,'-1'
comments=scupstr.info.comments
n=n_elements(comments)
FOR i=0,n-1 DO BEGIN
  IF trim(comments[i]) NE '' THEN printf,lout,comments[i]
ENDFOR 
printf,lout,'-1'

free_lun,lout

END
