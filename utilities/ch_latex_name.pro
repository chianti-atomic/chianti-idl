

FUNCTION ch_latex_name, ionname, apj=apj, diel=diel, quiet=quiet

;+
; NAME:
;     CH_LATEX_NAME
;
; PURPOSE:
;     Converts a CHIANTI format ion name to a latex string.
;
; CATEGORY:
;     CHIANTI; text; format.
;
; CALLING SEQUENCE:
;	Result = CH_LATEX_NAME( IonName )
;
; INPUTS:
;     Ionname:  A string containing the name of an ion in CHIANTI format.
;               For example, 'o_6' for O VI.
;
; KEYWORD PARAMETERS:
;     APJ:  Returns the latex string in the format used by The Astrophysical
;           Journal.
;     QUIET:  If set, then no informational messages are printed.
;
; OUTPUTS:
;     A string containing the latex format of the ion name. For example,
;     "fe_13" will be converted to "\ion{Fe}{xiii}". If the keyword /apj is
;     given, then the ion name is given as "\ion{Fe}{13}". For a dielectronic
;     ion (e.g., "fe_17d"), the d is stripped out and the output diel is set
;     to 1.
;
;     If a problem is found, then an empty string is returned.
;
; OPTIONAL OUTPUTS:
;     Diel:  If IONNAME ends with a "d" (indicating it is a dielectronic ion),
;            then diel is returned as 1. Otherwise it is returned as 0.
;
; EXAMPLE:
;     IDL> print,ch_latex_name('fe_17')
;     \ion{Fe}{xvii}
;     IDL> print,ch_latex_name('fe_17',/apj)
;     \ion{Fe}{17}
;     IDL> print,ch_latex_name('fe_17d',diel=diel)
;     \ion{Fe}{17}
;     IDL> print,diel
;            1
;
; MODIFICATION HISTORY:
;     Ver.1, 20-Nov-2023, Peter Young
;-


IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> latex_name = ch_latex_name( ionname [,/apj, diel=, /quiet] )'
  return,''
ENDIF 

diel=0
IF ionname.EndsWith('d') THEN BEGIN
  IF NOT keyword_set(quiet) THEN message,/info,/cont,'This is a dielectronic ion. The "d" is stripped out and the keyword diel is set to 1.'
  i=ionname.LastIndexOf('d')
  iname=ionname.Substring(0,i-1)
  bits=iname.split('_')
  diel=1
ENDIF ELSE BEGIN 
  bits=ionname.split('_')
ENDELSE 



IF n_elements(bits) NE 2 THEN BEGIN
  IF NOT keyword_set(quiet) THEN message,/info,/cont,'There is an error in the format of the ion name. Returning...'
  return,''
ENDIF 

eltname=bits[0].CapWords()

ion_num=fix(bits[1])

IF keyword_set(apj) THEN BEGIN
  latex_string='\ion{'+eltname+'}{'+trim(ion_num)+'}'
  return,latex_string
ENDIF

;
; The method below works as long as ion_num is less than 50.
;
nx=ion_num/10
roman_num=''
IF nx GE 1 THEN BEGIN
  FOR i=0,nx-1 DO roman_num=roman_num+'x'
ENDIF

ix=ion_num-nx*10
rom_all=['','i','ii','iii','iv','v','vi','vii','viii','ix']
roman_num=roman_num+rom_all[ix]

latex_string='\ion{'+eltname+'}{'+roman_num+'}'

return,latex_string

END
