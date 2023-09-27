
FUNCTION ch_masterlist, filename=filename, count=count, comments=comments, $
                        check=check

;+
; NAME:
;     CH_MASTERLIST
;
; PURPOSE:
;     Reads the CHIANTI masterlist of ions.
;
; CATEGORY:
;     CHIANTI; read.
;
; CALLING SEQUENCE:
;     Result = CH_MASTERLIST( )
;
; INPUTS:
;     None.
;
; OPTIONAL INPUTS:
;     Filename: The name of the masterlist file. If not specified, then
;               the routine looks for !xuvtop/masterlist/masterlist.ions.
;
; KEYWORD PARAMETERS:
;     CHECK:  If set, then the routine checks to see if the SCUPS file
;             is present for all of the ions in the database. If not,
;             then the list of ions for which the file is missing is
;             listed to the screen.
;
; OUTPUTS:
;     A string array containing a list of ions in the CHIANTI format. For
;     example, 'fe_13' for Fe XIII. The "dielectronic" ions have a "d"
;     appended to their name.
;
; OPTIONAL OUTPUTS:
;     Count:  The number of ions in the masterlist.
;     Comments: String array containing comments at the bottom of the file
;               (if available).
;
; EXAMPLE:
;     IDL> mlist=ch_masterlist()
;     IDL> mlist=ch_masterlist(file='my_masterlist.ions')
;     IDL> mlist=ch_masterlist(comments=comments,count=count)
;     IDL> mlist=ch_masterlist(/check)
;
; MODIFICATION HISTORY:
;     Ver.1, 27-Sep-2023, Peter Young
;-


count=0

IF n_elements(filename) EQ 0 THEN BEGIN
  dir=concat_dir(!xuvtop,'masterlist')
  filename=concat_dir(dir,'masterlist.ions')
ENDIF

chck=file_search(filename,count=ct)
IF ct EQ 0 THEN BEGIN
  message,/cont,/info,'The masterlist file was not found. Returning...'
  return,''
ENDIF

openr,lin,filename,/get_lun

str1=''
mlist=''
comments=''
swtch=0b
WHILE eof(lin) EQ 0 DO BEGIN
  readf,lin,str1
  IF trim(str1) EQ '-1' THEN BEGIN
    swtch=1b
    CONTINUE
  ENDIF
 ;
  ionname=trim(strmid(str1,0,7))
 ;
  IF swtch THEN BEGIN
    comments=[comments,trim(str1)]
  ENDIF ELSE BEGIN
    mlist=[mlist,ionname]
  ENDELSE 
ENDWHILE

mlist=mlist[1:*]
IF n_elements(comments) GT 1 THEN comments=comments[1:*]

free_lun,lin

count=n_elements(mlist)

;
; The section below checks for the existence of the SCUPS files.
;
IF keyword_set(check) THEN BEGIN
  chck_all=bytarr(count)
  FOR i=0,count-1 DO BEGIN
    convertname,mlist[i],iz,ion
    zion2filename,iz,ion,fname
    fname=fname+'.scups'
    chck=file_info(fname)
    chck_all[i]=chck.exists
  ENDFOR
 ;
  k=where(chck_all EQ 0b,nk)
  IF nk EQ 0 THEN BEGIN
    message,/info,/cont,'The SCUPS file is present for all the ions in the masterlist.'
  ENDIF ELSE BEGIN
    message,/info,/cont,'The SCUPS file is missing for '+trim(nk)+' ions in the database:'
    FOR j=0,nk-1 DO print,'    '+mlist[k[j]]
  ENDELSE
 ;
ENDIF 

return,mlist

END
