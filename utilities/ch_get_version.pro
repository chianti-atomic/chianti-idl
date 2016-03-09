
FUNCTION ch_get_version

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
; OUTPUTS:
;      A string containing the CHIANTI version.
;
; MODIFICATION HISTORY:
;      Ver.1, 30-Nov-2015, Peter Young
;-

verfile=concat_dir(!xuvtop,'VERSION')

str1=''
openr,lin,verfile,/get_lun
readf,lin,str1
free_lun,lin

return,trim(str1)

END
