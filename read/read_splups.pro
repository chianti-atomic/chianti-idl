;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an atomic database package for the calculation of
;       continuum and emission line spectra from astrophysical plasmas. It 
;       is a 
;       collaborative project involving the Naval Research Laboratory
;       (Washington DC, USA), the Arcetri Observatory (Firenze, Italy), and the
;       Cambridge University (United Kingdom).
;
;
; NAME:
;	READ_SPLUPS
;
; PURPOSE:
;
;	to read file containing spline fits to the Burgess-Tully scaled
;       collision strengths
;
; CATEGORY:
;
;	science.
;
; CALLING SEQUENCE:
;
;       READ_SPLUPS, File, Splstr, Splref
;
;
; INPUTS:
;
;	File:	the name of the input file, i.e. !xuvtop/si/si_4/si_4.splups
;
;
; OUTPUTS:
;
;       SPLSTR  Structure containing the data from the file. The tags are 
;               as follows:
;
;               .lvl1   lower level index
;               .lvl2   upper level index
;               .t_type transition type
;               .gf     gf value
;               .de     Delta-E for transition (rydbergs)
;               .c_ups  the scaling parameter
;               .nspl   
;               .spl    Vector of length 9, containing spline points
;
;
; OPTIONAL OUTPUTS
;
;       SPLREF  String array containing references.
;
; KEYWORDS
;
;       PROT    Allows reading of .psplups files for proton rates.
;
;
; PROCEDURE:
;
;	see Burgess and Tully, 1992, Astronomy and Astrophysics, 254, 436.
;
; EXAMPLE:
;
;       > read_splups, !xuvtop+'/si/si_4/si_4.splups',splstr,splref
;
; PROGRAMMING NOTES
;
;       This routine is marginally quicker (20-25%) reading the .splups 
;       files than Ken's original routine. The improvement in speed is 
;       through minimising the lines of code in the WHILE loop.
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	March 1996:     Version 2.0
;
;       Ver.3, 23-Jan-01, Peter Young
;                completely revised. Now reads into a structure and 
;                handles 9 point spline fits.
;
;       Ver.4, 26-Jan-01, Peter Young
;                speeded up routine
;
;       Ver.5, 22-Mar-01, Peter Young
;                now checks if file exists
;
;-
pro read_splups,name,splstr,splref,prot=prot

IF n_params(0) LT 2 THEN BEGIN
   print,' '
   print,'    IDL> read_splups,upsname,splstr,splref [, /prot]'
   print,' '
   return
endif

result=findfile(expand_path(name))
IF result[0] EQ '' THEN BEGIN
  splstr=-1
  splref='File does not exist!'
  return
ENDIF


IF keyword_set(prot) THEN BEGIN
  m=0 
  format='(3i3,12e10.0)'
ENDIF ELSE BEGIN
  m=2
  format='(5i3,12e10.0)'
ENDELSE


data1=dblarr(11+m)
data2=dblarr(15+m)

dd1=dblarr(16+m)

openr,lur,name,/get_lun

string1=''
tst1=0


WHILE tst1 EQ 0 DO BEGIN
  readf,lur,string1
  IF strcompress(string1,/rem) NE '-1' THEN BEGIN
    IF strlen(string1) LT 100 THEN BEGIN
      reads,format=format,string1,data1
      dd1=[[dd1],[data1,-1,-1,-1,-1,5]]
    ENDIF ELSE BEGIN
      reads,format=format,string1,data2
      dd1=[[dd1],[data2,9]]
    ENDELSE
  ENDIF ELSE BEGIN
    tst1=1
  ENDELSE
ENDWHILE

splref=''
WHILE tst1 EQ 1 DO BEGIN
  readf,lur,string1
  IF (strcompress(string1,/rem) NE '-1') THEN splref=[splref,string1] $
  ELSE tst1=2
ENDWHILE

free_lun,lur

dd1=dd1[*,1:*]
siz=size(dd1)
IF siz[0] EQ 1 THEN rep=1 ELSE rep=siz[2]

str={ lvl1: 0, lvl2: 0, t_type: 0, gf: 0., de: 0., c_ups: 0., $
      nspl:0, spl: fltarr(9)}
splstr=replicate(str,rep)
;
splstr.lvl1=reform(round(dd1[0+m,*]))
splstr.lvl2=reform(round(dd1[1+m,*]))
splstr.t_type=reform(round(dd1[2+m,*]))
splstr.gf=reform(dd1[3+m,*])
splstr.de=reform(dd1[4+m,*])
splstr.c_ups=reform(dd1[5+m,*])
splstr.nspl=reform(round(dd1[15+m,*]))
splstr.spl=dd1[6+m:14+m,*]


END
