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
;       Ver.6, 17-Jul-2012, Ken Dere
;               can handle files with 5 through 9 spline points
;
;       Ver.7, 10-Mar-2014 Ken Dere
;               corrected for reading lines with trailing white space
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
  status = 0
  return
ENDIF Else begin
  status = 1
endelse


lvl1 = -1
lvl2 = -1
t_type = -1
gf = -1.
de = -1.
c_ups = -1.
nspl = -1
spl = fltarr(9)



openr,lur,name,/get_lun

string1=''
tst1=0


WHILE tst1 EQ 0 DO BEGIN
  readf,lur,string1
  nchar = strlen(strtrim(string1))
;   so far, the number of splups values is limited to 9
  IF strcompress(string1,/rem) NE '-1' THEN BEGIN
    if keyword_set(prot) then begin
		nsplups = (nchar - 3*3 - 3*10)/10 < 9
		nvar = (nchar - 3*3)/10
        splformat = '(3i3,'+string(nvar)+'e10.0)'
        datav = dblarr(nsplups+3+3)
        reads,format=splformat,string1,datav
        lvl1 = [lvl1,fix(datav[0])]
        lvl2 = [lvl2,fix(datav[1])]
        t_type = [t_type,fix(datav[2])]
        gf = [gf, datav[3]]
        de = [de, datav[4]]
        c_ups = [c_ups,datav[5]]
        nspl = [nspl, nsplups]
        spl1 = datav[6:nsplups+6-1]
		; nsplups should equal 9 for a psplups file
        nodata = 9 - nsplups
        for inot=0,nodata-1 do begin
            spl1 = [spl1,-1.]
        endfor
        spl = [spl, spl1]
    endif else begin
  		nsplups = (nchar - 5*3 - 3*10)/10 < 9
  		nvar = (nchar - 5*3)/10
        splformat = '(5i3,'+string(nvar)+'e10.0)'
        datav = dblarr(nsplups+3+5)
        reads,format=splformat,string1,datav
        lvl1 = [lvl1,fix(datav[2])]
        lvl2 = [lvl2,fix(datav[3])]
        t_type = [t_type,fix(datav[4])]
        gf = [gf, datav[5]]
        de = [de, datav[6]]
        c_ups = [c_ups,datav[7]]
        nspl = [nspl, nsplups]
        spl1 = datav[8:nsplups+8-1]
        nodata = 9 - nsplups
        for inot=0,nodata-1 do begin
            spl1 = [spl1,-1.]
        endfor
        spl = [spl, spl1]
    endelse
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

str={ lvl1: 0, lvl2: 0, t_type: 0, gf: 0., de: 0., c_ups: 0., $
nspl:0, spl: fltarr(9)}
splstr=replicate(str,n_elements(lvl1)-1)
;
splstr.lvl1 = lvl1[1:*]
splstr.lvl2 = lvl2[1:*]
splstr.t_type = t_type[1:*]
splstr.gf = gf[1:*]
splstr.de = de[1:*]
splstr.c_ups = c_ups[1:*]
splstr.nspl = nspl[1:*]
nlines = n_elements(spl)/9
spl = reform(spl,9,nlines)
spl = spl[*,1:*]
splstr.spl = spl

END
