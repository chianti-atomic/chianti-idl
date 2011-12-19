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
;
;	provide identification strings
;
; CATEGORY:
;
;	database.
;
; CALLING SEQUENCE:
;
;       Z2ELEMENT, Iz, Name
;
;
; INPUTS:
;
;	Iz:  nuclear charge of ion of interest, i.e. 26 for Fe
;
;
; OUTPUTS:
;
;	Name:  a string identifying the element
;
;
;
; EXAMPLE:
;
;             > z2element,2,name
;             > print,name
;             > He  
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	March 1996:     Version 2.0
;
;-
pro z2element,iz,element
;
if n_params() lt 2 then begin
   print,' '
   print,'  IDL> z2element,z,element'
   print,' i.e.> z2element,5,element'
   print,'         returns element = ''boron'' '
   print,' '
   return
endif
;
;
ele=['hydrogen','helium','lithium','beryllium','boron','carbon','nitrogen','oxygen','fluorine']
ele=[ele,'neon','sodium','magnesium','aluminum','silicon','phosphorus','sulfur','chlorine','argon']
ele=[ele,'potassium','calcium','scandium','titanium','vanadium','chromium','manganese','iron','cobalt']
ele=[ele,'nickel','copper','zinc']
;
if iz le 30 then begin
   element=ele(iz-1)
   return
endif else begin
   element=' '
   print, ' ascii characters not tabulated for iz=',iz
   return
endelse
;
end
