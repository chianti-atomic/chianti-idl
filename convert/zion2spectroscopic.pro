;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory (USA), the University of Florence (Italy), the
;       University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
;
; NAME:
;	ZION2SPECTROSCOPIC
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
;       ZION2SPECTROSCOPIC, Iz, Ion, Name
;
;
; INPUTS:
;
;	Iz:  nuclear charge of ion of interest, i.e. 26 for Fe
;       Ion:   charge state of ion of interest, i.e. 2 for Fe II	
;
;
; OUTPUTS:
;
;	Name:  the spectroscopic notation for the ion, i.e. 'Fe II'
;
;
;
; EXAMPLE:
;
;             > zion2spectroscopic,26,2,name
;             > print,name
;             > Fe II   
;
; WRITTEN     :  Ken Dere
;
; MODIFICATION HISTORY:
;
;	March 1996:     Version 2.0
;
;       V.3, 25-May-2002, Giulio Del Zanna (GDZ)
;            added the DIELECTRONIC keyword.
;
; VERSION     : 3, 25-May-2002
;
;-
pro zion2spectroscopic,iz,ion,snote, dielectronic=dielectronic
;
;  convert z, ion to spectroscopic notation
;
if n_params() lt 2 then begin
   print,' '
   print,'   IDL> zion2spectroscopic,z,ion,spectroscopic, diel=diel'
   print,'  i.e.> zion2spectroscopic,7,3,spectroscopic'
   print,'            yields spectroscopic = ''N III''  '
   print,' '
   return
end
;
element=['H','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si',$
         'P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni',$
         'Cu','Zn']
;
ionstage=['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII',$
         'XIV','XV','XVI','XVII','XVIII','XIX','XX','XXI',' XXII','XXIII','XXIV',$
         'XXV','XXVI','XXVII','XXVIII','XXIX','XXX','XXXI']
;
snote=element(iz-1)+' '+ionstage(ion-1)

IF keyword_set(dielectronic) THEN snote = snote+' d'
;
end
