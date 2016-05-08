;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving
;       George Mason University USA), the University of Michigan (USA),
;        and Cambridge University (UK).
;
;
; NAME:
;	ZION2FILENAME
;
; PURPOSE:
;
;	help locate CHIANTI database files
;
; CATEGORY:
;
;	database.
;
; CALLING SEQUENCE:
;
;       ZION2FILENAME, Iz, Ion, Filename
;
;
; INPUTS:
;
;	Iz:  nuclear charge of ion of interest, i.e. 26 for Fe
;       Ion:   charge state of ion of interest, i.e. 2 for Fe II	
;
;  KEYWORDS:
;       
;       diel:  set if excitation of this ion is by dielectronic
;		recombination
;
; OUTPUTS:
;
;	Filename:  the complete filename and path specification for generic
;                  CHIANTI database file, i.e. '.elvlc' suffix is not included
;
;
; RESTRICTIONS:
;
;	!xuvtop must be set
;
;
; EXAMPLE:
;
;             > zion2filename,26,2,filename
;             > print,filename
;             > /data1/xuv/fe/fe_2/fe_2   assuming !xuvtop = '/data1/xuv'
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	March 1996:     Version 2.0
;       Sept  1996:     Modified for use with VMS
;       December 1998:  Modified to diel keyword
;
;       V.5, 29-May-2002, Giulio Del Zanna (GDZ) 
;                   generalized directory concatenation to work for
;                   Unix, Windows  and VMS. Added keyword name to output just
;                   the name of the file and changed the dielectronic  keyword. 
;
; VERSION     :   5, 29-May-2002 
;
;-
pro zion2filename,z,ion,filename,name=name, dielectronic=dielectronic
;
;  convert Z, ionization state to the filename
;  i.e.     z=26  ion=24 > !xuvtop/fe/fe_24
;
if n_params() lt 3 then begin
   print,' '
   print,'   IDL> zion2filename,z,ion,filename, [/diel]
   print,'  i.e.> zion2filename,26,11,filename
   print,'         yields filename = ''!xuvtop/fe/fe_11/fe_11''   '
   print,'  or '
   print,'  i.e.> zion2filename,26,11,filename,/diel
   print,'         yields filename = ''!xuvtop/fe/fe_11d/fe_11d''   '
   print,' '
   return
endif
;
z_lbl=['h','he','li','be','b','c','n','o','f','ne','na','mg','al','si',$
       'p','s','cl','ar','k','ca','sc','ti','v','cr','mn','fe','co','ni',$
      'cu','zn']

IF z GT 0 THEN BEGIN 
dir=concat_dir(!xuvtop,strtrim(z_lbl(z-1),2))

name=strtrim(z_lbl(z-1),2)+'_'+strtrim(string(ion,'(i2)'),2)

IF  keyword_set(dielectronic) THEN name=name+'d' 

filename=concat_dir(concat_dir(dir, name), name)

ENDIF ELSE BEGIN 
filename=''
name = ''
END 

END 


