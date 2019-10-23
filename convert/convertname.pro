;+
; PROJECT:  CHIANTI
;
;       CHIANTI is an atomic database package for the calculation of
;       continuum and emission line spectra from astrophysical plasmas. It is a 
;       collaborative project involving George Mason University USA),
;        the University of Michigan (USA), and Cambridge University (UK).

; NAME:
;	convertname
;
; PURPOSE:
;	Ion names as character strings are converted into
;	numerical values (note c_2 is C II or C^+1
;	in spectroscopic or atomic notation)
;
; CATEGORY:
;	naming utility
;
; CALLING SEQUENCE:
;       CONVERTNAME,Name,Iz,Ion
;
; INPUTS:
;	Name:   Name of ion in CHIANTI format, e.g., 'c_2' for C
;               II. Can also be just an element.
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
;       Ver.5, 13-Aug-2016, Peter Young
;           Modified so that if only the element is given, then the
;           routine returns the correct IZ value.
;-
pro convertname,name,iz,ion,dielectronic=dielectronic
;
;
if n_params() lt 2 then begin
    print,' '
    print,'   IDL> convertname,name,iz,ion,diel=diel'
    print,'  e.g.> convertname,''c_5'',iz,ion,diel'
    print,'          giving iz = 6 and ion = 5, diel=0, for C V'
    print,'  e.g.> convertname,''c_4d'',iz,ion,diel'
    print,'          giving iz = 6 and ion = 4, diel=1, for the C IV dielectronic transitions'
    print,"  e.g.> convertname,'fe',iz"
    print,'          giving iz = 26.'
    print,''
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


;
; Extract element from string.
;
first=strupcase(bits[0])
ind=where(first EQ zlabl)
IF ind[0] EQ -1 THEN BEGIN
  print,'%CONVERTNAME: Error, element not in database  ',name
  iz=0 &  ion=0
  return
ENDIF
iz=ind[0]+1

;
; Extract ion number (and optionally dielectronic id).
;
IF n_elements(bits) EQ 2 THEN BEGIN 
  bits2=str_sep(bits[1],'d')
  ion=fix(bits2[0])
  IF n_elements(bits2) EQ 2 THEN dielectronic=1 ELSE dielectronic=0
ENDIF ELSE BEGIN
  ion=0
ENDELSE 
  

END
