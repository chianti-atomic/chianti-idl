
FUNCTION ch_get_version, ssw=ssw, online=online, idl=idl

;+
; NAME:
;      CH_GET_VERSION
;
; PURPOSE:
;      Returns a string containing the version of the CHIANTI database
;      being used. If /idl is set, then the version of the IDL
;      software is returned.
;
; CATEGORY:
;      CHIANTI; utility.
;
; CALLING SEQUENCE:
;      Result = CH_GET_VERSION()
;      
; INPUTS:
;      None.
;
; KEYWORD PARAMETERS:
;      SSW:  If set, then the routine checks the version of the
;            Solarsoft (SSW) version of the database. Requires an
;            internet connection.
;      IDL:  If set, then the routine reads the VERSION file for the
;            IDL routines.
;      ONLINE:  This is keyword is obsolete, but retained for
;            backwards compatibility.
;
; OUTPUTS:
;      A string containing the CHIANTI version (either database or
;      idl) 
;
; MODIFICATION HISTORY:
;      Ver.1, 30-Nov-2015, Peter Young
;      Ver.2, 4-Feb-2019, Peter Young
;        Added /ssw and /online keywords.
;      Ver.3, 19-Jul-2019, Peter Young
;        Changed address of SSW file (ftp->https).
;      Ver.4, 2-Jul-2020, Peter Young
;        Added /idl keyword and disabled /online keyword.
;-

fname='VERSION'

IF keyword_set(idl) THEN dir='idl' ELSE dir='dbase'

CASE 1 OF
  keyword_set(ssw): BEGIN
    url='https://sohoftp.nascom.nasa.gov/solarsoft/packages/chianti/'+dir+'/'+fname
    out_dir=getenv('IDL_TEMPDIR')
    sock_get,url,out_dir=out_dir,local_file=verfile,status=status
    IF status EQ 0 THEN BEGIN
      print,'% CH_GET_VERSION: file not found online.'
      return,''
    ENDIF 
  END
  keyword_set(online): BEGIN
    print,'% CH_GET_VERSION: this keyword is obsolete. Returning...'
    return,''
  END
  ELSE: BEGIN
    IF keyword_set(idl) THEN BEGIN
     ;
     ; A user may keep the CHIANTI IDL routines separate from the
     ; database, so the code backs out the directory location using
     ; the location of ch_get_version.
     ;
      which,'ch_get_version',outfile=outfile,/quiet
      dir1=file_dirname(outfile)
      dir2=file_dirname(dir1)
      verfile=concat_dir(dir2,fname)
    ENDIF ELSE BEGIN
      verfile=concat_dir(!xuvtop,fname)
    ENDELSE 
  END 
ENDCASE 

str1=''
openr,lin,verfile,/get_lun
readf,lin,str1
free_lun,lin

return,trim(str1)

END
