
FUNCTION ch_get_version, ssw=ssw, online=online

;+
; NAME:
;      CH_GET_VERSION
;
; PURPOSE:
;      Returns a string containing the version of CHIANTI being used.
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
;      ONLINE:  If set, then the routine checks the version of the
;            database is maintained by the CHIANTI team at
;            http://chiantidatabase.org. Requires an internet
;            connection. 
;
; OUTPUTS:
;      A string containing the CHIANTI version.
;
; MODIFICATION HISTORY:
;      Ver.1, 30-Nov-2015, Peter Young
;      Ver.2, 4-Feb-2019, Peter Young
;        Added /ssw and /online keywords.
;      Ver.3, 19-Jul-2019, Peter Young
;        Changed address of SSW file (ftp->https).
;-

fname='VERSION'

CASE 1 OF
  keyword_set(ssw): BEGIN
    url='https://sohoftp.nascom.nasa.gov/solarsoft/packages/chianti/dbase/'+fname
    out_dir=getenv('IDL_TEMPDIR')
    sock_get,url,out_dir=out_dir,local_file=verfile,status=status
    IF status EQ 0 THEN BEGIN
      print,'% CH_GET_VERSION: file not found online.'
      return,''
    ENDIF 
  END
  keyword_set(online): BEGIN
    url='http://www.chiantidatabase.org/sswidl/dbase/'+fname
    out_dir=getenv('IDL_TEMPDIR')
    sock_get,url,out_dir=out_dir,local_file=verfile, status=status
    IF status EQ 0 THEN BEGIN
      print,'% CH_GET_VERSION: file not found online.'
      return,''
    ENDIF 
  END
  ELSE: verfile=concat_dir(!xuvtop,fname)
  
ENDCASE 

str1=''
openr,lin,verfile,/get_lun
readf,lin,str1
free_lun,lin

return,trim(str1)

END
