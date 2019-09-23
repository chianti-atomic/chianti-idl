
PRO z2element,iz,element, symbol=symbol, english=english, lower_CASE=lower_case


;+
; NAME:
;      Z2ELEMENT
;
; PURPOSE:
;      Convert the atomic number for an element to either the element
;      or the chemical symbol for the element.
;
; CATEGORY:
;      CHIANTI; notation conversion.
;
; CALLING SEQUENCE:
;      Z2ELEMENT, Iz, Name
;
; INPUTS:
;      Iz:    Atomic number for element. Can be an array.
;
; OUTPUTS:
;      Name:  A string array of same size as IZ containing the
;             (American) names of the elements corresponding to IZ.  
;
; KEYWORDS:
;      SYMBOL:  If set, then the chemical symbol for the element is
;               returned instead.
;      ENGLISH: If set, then the English names of the elements will be
;               returned.
;      LOWER_CASE: If set, then names or symbols will be returned in
;               lower case format.
;
; EXAMPLES:
;      IDL> z2element,2,name
;      IDL> print,name
;      helium
;
;      IDL> z2element,[8,26],name,/symbol
;      IDL> print,name
;      O Fe
;
;      IDL> z2element,[13,16],name,/english,/lower_case
;      IDL> print,name
;      aluminium sulphur
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	March 1996:     Version 2.0
;
;       Ver. 3, 26-Jun-2012, Peter Young
;         added /symbol keyword.
;
;       Ver.4, 25-Jan-2018, Peter Young
;         The input IZ can now be an array; added /english and
;         /lower_case keywords; the element names now start with a
;         capital letter; removed extra space from names of some
;         symbols. 
;-

if n_params() lt 2 then begin
   print,'Use:  IDL> z2element, z, element [, /symbol, /english, /lower_case ]'
   print,'  e.g., z2element,5,element'
   print,'         returns element = ''boron'' '
   print,'        z2element,5,element,/symbol'
   print,'         returns element = ''B'' '
   return
endif

niz=n_elements(iz)
element=strarr(niz)
k=where(iz GE 1 AND iz LE 30,nk)
IF nk EQ 0 OR nk NE niz THEN BEGIN 
  print,'% Z2ELEMENT: The input IZ must take values between 1 and 30 (inclusive).'
  print,'             Other values will be returned as empty strings.'
ENDIF
IF nk EQ 0 THEN BEGIN
  element=''
  return
ENDIF 

elt_short=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si',$
           'P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni',$
           'Cu','Zn']

;
ele=['Hydrogen','Helium','Lithium','Beryllium','Boron','Carbon','Nitrogen','Oxygen','Fluorine']
ele=[ele,'Neon','Sodium','Magnesium','Aluminum','Silicon','Phosphorus','Sulfur','Chlorine','Argon']
ele=[ele,'Potassium','Calcium','Scandium','Titanium','Vanadium','Chromium','Manganese','Iron','Cobalt']
ele=[ele,'Nickel','Copper','Zinc']

IF keyword_set(english) THEN BEGIN
  ele[12]='Aluminium'
  ele[15]='Sulphur'
ENDIF 

IF keyword_set(symbol) THEN BEGIN
  element[k]=elt_short[iz[k]-1]
ENDIF ELSE BEGIN 
  element[k]=ele[iz[k]-1]
ENDELSE

;
; Make sure ELT is not returned as an array if there's only element.
;
IF n_elements(element) EQ 1 THEN element=element[0]

IF keyword_set(lower_CASE) THEN element=strlowcase(element)

END 
