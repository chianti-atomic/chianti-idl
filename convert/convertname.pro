;+
;
; PROJECT:  CHIANTI
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an atomic database package for the calculation of
;       continuum and emission line spectra from astrophysical plasmas. It is a 
;       collaborative project involving the Naval Research Laboratory
;       (Washington DC, USA), the Arcetri Observatory (Firenze, Italy), and the
;       Cambridge University (United Kingdom).
;
; NAME:
;	convertname
;
; PURPOSE:
;	Ion names as character strings are converted into
;	numerical values (note c_2 is C II or C^+1
;	in spectroscopic or atomic notation)
;
; CATEGORY:
;	
;	naming utility
;
; CALLING SEQUENCE:
;
;       CONVERTNAME,Name,Iz,Ion
;
;
; INPUTS:
;	Name:   such as 'c_2'
;
;
; OUTPUTS:
;
;	Iz:  nuclear charge Z  (6 for 'c_2', the equivalent of C II)
;       Ion:  ionization stage:  (2 for 'c_2')
;
; OPTIONAL OUTPUTS
;
;       DIELECTRONIC   Set to 1 if NAME has a 'd' appended to it 
;                      (indicating dielectronic recombination data) else 
;                      set to 0
;
; EXAMPLE:
;
;                     > convertname,'c_2',iz,ion
;                     > print,iz,ion
;                     > 6,2
;
;                     > convertname,'o_6d',iz,ion
;                     > print,iz,ion
;                     > 8,6
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	March 1996:     Version 2.0
;       October 1999:   Version 3.  by kpd
;
;       Ver.4, 11-Dec-01, Peter Young
;           Revised routine, removing ch_repstr call.
;           Added DIELECTRONIC optional output.
;
;-
pro convertname,name,iz,ion,dielectronic=dielectronic
;
;to convert names such as c_2 to iz=6 and ion=2
;
if n_params() lt 3 then begin
    print,' '
    print,'   IDL> convertname,name,iz,ion'
    print,'  i.e.> convertname,''c_5'',iz,ion'
    print,'          giving iz = 6 and ion = 5, for C V'
    print,' '
    return
endif
;
;
tname=strtrim(name,2)
;
;
zlabl=['H','He','Li','Be','B','C','N','O','F','Ne','Na',$
       'Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr', $
       'Mn','Fe','Co','Ni','Cu','Zn']
;
zlabl=strupcase(zlabl)


bits=str_sep(tname,'_')
IF n_elements(bits) NE 2 THEN BEGIN
  iz=0
  ion=0
ENDIF ELSE BEGIN
  first=strupcase(bits[0])
  ind=where(first EQ zlabl)
  IF ind[0] EQ -1 THEN BEGIN
    print,'%CONVERTNAME: Error, ion not in database  ',name
    iz=0 &  ion=0
    return
  ENDIF
  iz=ind[0]+1
 ;
  bits2=str_sep(bits[1],'d')
  ion=fix(bits2[0])
  IF n_elements(bits2) EQ 2 THEN dielectronic=1 ELSE dielectronic=0
ENDELSE

END
