
FUNCTION ch_conf_rem_parent, input, nparent=nparent, quiet=quiet

;+
; NAME:
;     CH_CONF_REM_PARENT
;
; PURPOSE:
;     Removes any parent term designations within the CHIANTI
;     configuration descriptor. Parent terms are indicated by pairs of
;     brackets. For example '(3P)'.
;
; CATEGORY:
;     CHIANTI; string processing.
;
; CALLING SEQUENCE:
;     Result = CH_CONF_REM_PARENT( Input )
;
; INPUTS:
;     Input:  A string containing the CHIANTI configuration
;             description. 
;
; KEYWORD PARAMETERS:
;     QUIET:  If set, then information messages are not printed. 
;
; OUTPUTS:
;     A string containing the configuration description, but with any
;     parent terms removed. The routine checks to make sure there are
;     equal numbers of open and closed brackets, and that the open
;     brackets occur before teh closed brackets. If problems are
;     found, the input string is returned.
;
; OPTIONAL OUTPUTS:
;     Nparent:  The number of parent terms that were found.
;
; EXAMPLE:
;     IDL> config='3d3.(2D)4s.(3D)4p'
;     IDL> print,ch_conf_rem_parent(config)
;     3d3.4s.4p
;
;     IDL> zion2filename,26,2,fname
;     IDL> read_elvlc,fname+'.elvlc',elvlc=elvlc
;     IDL> print,ch_conf_rem_parent(elvlc.data[0].conf)
;     3d6.4s
;
; MODIFICATION HISTORY:
;     Ver.1, 23-Jul-2020, Peter Young
;-

IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> result = ch_conf_rem_parent( input )'
  return,input
ENDIF 

;
; I add a space to input in case the close bracket is at the end of the string.
;
open=strsplit(input+' ','(',/extract,count=n1)
close=strsplit(input+' ',')',/extract,count=n2)

IF n1 EQ 1 AND n2 EQ 1 THEN BEGIN
  IF NOT keyword_set(quiet) THEN print,'% CH_CONF_REM_PARENT: No parent terms found. Returning...'
  nparent=0
  return,input
ENDIF

IF n1 NE n2 THEN BEGIN
  print,'% CH_CONF_REM_PARENT: There is a mismatch in the numbers of open and closed brackets: '+input
  print,'                      Returning...'
  return,input
ENDIF

n=n1-1

output=input
FOR i=0,n-1 DO BEGIN
  len=strlen(output)
  a=strpos(output,'(')
  b=strpos(output,')')
  IF b GT a THEN BEGIN
    output=strmid(output,0,a)+strmid(output,b+1)
  ENDIF ELSE BEGIN
    print,'% CH_CONF_REM_PARENT: The closed bracket occurs before the open bracket! '+input
    print,'                      Returning...'
    return,input
  ENDELSE 
ENDFOR 

nparent=n

return,output

END
