;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an atomic database package for the calculation of
;       astrophysical emission line spectra.  It is a collaborative project
;       involving Ken Dere (Naval Research Laboratory, Washington DC), 
;       Brunella Monsignori-Fossi and Enrico Landi (Arcetri Observatory, 
;       Florence), and Helen Mason and Peter Young (DAMTP, Cambridge Univ.).
;
;
; NAME:
;	ZION2NAME
;
; PURPOSE:
;
;       Help locate CHIANTI database files.
;
; CATEGORY:
;
;	Database.
;
; CALLING SEQUENCE:
;
;       ZION2NAME, Iz, Ion, Name [, /diel]
;
;
; INPUTS:
;
;	Iz:  nuclear charge of ion of interest, i.e. 26 for Fe
;       Ion:   charge state of ion of interest, i.e. 2 for Fe II
;
; KEYWORDS:
;
;	diel:   set if ion produced by dielectronic recombination	
;
;
; OUTPUTS:
;
;	Name:  the generic filename for a CHIANTI database file.
;
; EXAMPLE:
;
;             > zion2name,26,2,name
;             > print,name
;             > fe_2   
;
; MODIFICATION HISTORY:
;
;     Ver.1, Mar-96, Ken Dere
;     Ver.2, Sep-99, Ken Dere
;              added /diel keyword
;
;-

pro zion2name,z,ion,name,diel=diel

IF n_params() LT 3 THEN BEGIN
  print,'Use: IDL> zion2name, z, ion, name [, /diel]'
  return
ENDIF


z_lbl=['h','he','li','be','b','c','n','o','f','ne','na','mg','al','si', $
       'p','s','cl','ar','k','ca','sc','ti','v','cr','mn','fe','co','ni', $
       'cu','zn']

IF z GT 0 THEN BEGIN
  name=strtrim(z_lbl(z-1),2)+'_'+strtrim(string(ion,'(i2)'),2)
  IF keyword_set(diel) THEN name=name+'d'
ENDIF ELSE name=''

END
