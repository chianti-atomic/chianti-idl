;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an atomic database package for the calculation of
;       continuum and emission line spectra from astrophysical plasmas. It is a 
;       collaborative project involving the Naval Research Laboratory
;       (Washington DC, USA), the Arcetri Observatory (Firenze, Italy), and the
;       Cambridge University (United Kingdom).
;
;
; NAME:
;	Z2ELEMENT
;
; PURPOSE:
;
;	Convert the atomic number for an element to either the element
;       or the chemical symbol for the element.
;
; CATEGORY:
;
;	database.
;
; CALLING SEQUENCE:
;
;       Z2ELEMENT, Iz, Name , [ /symbol ]
;
; INPUTS:
;
;	Iz:  nuclear charge of ion of interest, i.e. 26 for Fe
;
; OUTPUTS:
;
;	Name:  A string giving the element name. If /short is set,
;              then the chemical symbol for the element is returned.
;
; KEYWORDS:
;
;       SYMBOL If set then the chemical symbol for the element is
;              returned. 
;
; EXAMPLES:
;
;             > z2element,2,name
;             > print,name
;             > helium  
;
;             > z2element,2,name,/symbol
;             > print,name
;             > He
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	March 1996:     Version 2.0
;
;       Ver. 3, 26-Jun-2012, Peter Young
;         added /symbol keyword.
;
;-
pro z2element,iz,element, symbol=symbol
;
if n_params() lt 2 then begin
   print,' '
   print,'  IDL> z2element, z, element [, /symbol ]'
   print,'  e.g., z2element,5,element'
   print,'         returns element = ''boron'' '
   print,'        z2element,5,element,/symbol'
   print,'         returns element = ''B'' '
   print,' '
   return
endif
;

elt_short=['H','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si',$
           'P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni',$
           'Cu','Zn']

;
ele=['hydrogen','helium','lithium','beryllium','boron','carbon','nitrogen','oxygen','fluorine']
ele=[ele,'neon','sodium','magnesium','aluminum','silicon','phosphorus','sulfur','chlorine','argon']
ele=[ele,'potassium','calcium','scandium','titanium','vanadium','chromium','manganese','iron','cobalt']
ele=[ele,'nickel','copper','zinc']
;
IF iz LE 30 THEN BEGIN
  IF keyword_set(symbol) THEN BEGIN
    element=elt_short[iz-1]
  ENDIF ELSE BEGIN 
    element=ele[iz-1]
  ENDELSE
  return
ENDIF ELSE  BEGIN
  element=' '
  print, ' ascii characters not tabulated for iz=',iz
  return 
ENDELSE 
;
END 
