
FUNCTION ch_element_ions, element, exist=exist

;+
; NAME:
;      CH_ELEMENT_IONS()
;
; CATEGORY:
;      CHIANTI; formatting.
;
; PURPOSE:
;      For the specified element, this routine returns the list of all
;      ions of this element.
;
; CALLING SEQUENCE:
;	Result = CH_ELEMENT_IONS( Element )
;
; INPUTS:
;      Element:  Can be a string, e.g., 'fe' or 'Fe', or can be an
;                integer giving atomic number of element.
;
; OUTPUTS:
;      A string array giving the list of ions belonging to the element
;      in CHIANTI format. E.g., 'fe_1', 'fe_2', etc.
;
;      If the input is incorrect, then an empty string is returned.
;
; KEYWORDS:
;      EXIST:  If set, then only those ions that actually exist in
;              CHIANTI will be returned.
;
; CALLS:
;      CH_ION_EXIST, CONVERTNAME, Z2ELEMENT
;
; EXAMPLES:
;      IDL> ions=ch_element_ions(8)     ; for oxygen
;      IDL> ions=ch_element_ions('zn',/exist)
;
; MODIFICATION HISTORY:
;      Ver.1, 13-Aug-2016, Peter Young
;      Ver.2, 15-Aug-2016, Peter Young
;          Added /exist keyword.
;-

IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> ions = ch_element_ions( element )'
  print,''
  print,"   e.g., ions=ch_element_ions('fe')"
  print,"   e.g., ions=ch_element_ions(26)"
  print,''
ENDIF 

IF datatype(element) EQ 'STR' THEN BEGIN
  convertname,element,iz
  elt=trim(element)
  IF iz EQ 0 THEN return,''
ENDIF ELSE BEGIN
  iz=element
  z2element,element,elt,/symbol
  elt=trim(elt)
  IF elt EQ '' THEN return,''
ENDELSE 

output=''
FOR i=0,iz-1 DO BEGIN
  ionname=strlowcase(elt)+'_'+trim(i+1)
  chck=ch_ion_exist(ionname)
  IF chck EQ 1 THEN output=[output,ionname]
ENDFOR
output=output[1:*]

return,output

END
