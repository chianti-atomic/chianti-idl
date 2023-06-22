
PRO zion2spectroscopic,iz,ion,snote, dielectronic=dielectronic


;+
; NAME:
;	ZION2SPECTROSCOPIC
;
; PURPOSE:
;       Convert an atomic number and spectroscopic number pair to an ion name
;       in spectroscopic format.
;
; CATEGORY:
;       CHIANTI; convert.
;
; CALLING SEQUENCE:
;	ZION2SPECTROSCOPIC, IZ, ION
;
; INPUTS:
;	Iz:  nuclear charge of ion of interest, i.e. 26 for Fe
;       Ion:   charge state of ion of interest, i.e. 2 for Fe II	
;
; OPTIONAL INPUTS:
;	Parm2:	Describe optional inputs here. If you don't have any, just
;		delete this section.
;	
; KEYWORD PARAMETERS:
;       DIELECTRONIC:  If set, then a "d" is appended to the ion name to
;                      denote that it is a dielectronic ion.
;
; OUTPUTS:
;       The name of the ion in spectroscopic format. For example,
;       'Fe XXII' or 'O VI'.
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;
; EXAMPLE:
;       IDL> zion2spectroscopic,26,2,name
;       IDL> zion2spectroscopic,26,2,name,/diel
;
; MODIFICATION HISTORY:
;	V.2, Mar-1996, Ken Dere
;
;       V.3, 25-May-2002, Giulio Del Zanna (GDZ)
;            added the DIELECTRONIC keyword.
;
;       V.4, 12-Jun-2023, Peter Young
;            Used the trim routine to get rid of redundant spaces.
;-


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
element=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si',$
         'P','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni',$
         'Cu','Zn']
;
ionstage=['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII',$
         'XIV','XV','XVI','XVII','XVIII','XIX','XX','XXI',' XXII','XXIII','XXIV',$
         'XXV','XXVI','XXVII','XXVIII','XXIX','XXX','XXXI']
;
snote=trim(element(iz-1))+' '+trim(ionstage(ion-1))

IF keyword_set(dielectronic) THEN snote = snote+' d'
;
end
