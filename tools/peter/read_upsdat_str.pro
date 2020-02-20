;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an atomic database package for the calculation of
;       astrophysical emission line spectra.  It is a collaborative project
;       involving the Naval Research Laboratory, Washington DC, Arcetri 
;       Observatory, Florence, Cambridge University, Rutherford
;       Appleton Laboratory, and the University of Central Lancashire.
;
;
; NAME:
;	READ_UPSDAT_STR
;
; PURPOSE:
;
;	This is a modified version of Ken's READ_UPSDAT, which reads the 
;	upsdat data into a structure.
;
;	Note that the routine does not output the file references.
;
; CATEGORY:
;
;	science.
;
; CALLING SEQUENCE:  
;
;        READ_UPSDAT_STR, FILE, OMDAT_STRUC
;
; INPUTS:
;
;	File:	the name of the input file, i.e. !xuvtop/si/si_4/si_4.omdat
;
;		 
; OUTPUTS:
;
;	UPSDAT_STRUC	A structure with the following tags:
;
;		.lvl1	The lower level of the transition
;		.lvl2	The upper level
;		.de	The Delta-E for the transition (Ryd)
;		.gf	The gf value (for allowed transitions)
;		.temp	The temperatures at which the upsilons are tabulated
;		.ups	The upsilons
;
; AUTHOR
;
;	Peter Young, Rutherford Appleton Laboratory
;
; HISTORY
;
;	Ver.1	PRY, 26-Jan-00
;       Ver.2   PRY, 29-Apr-05
;            Now outputs references to REF.
;       Ver.3   PRY, 13-Jun-06
;            Changed name to read_upsdat_str
;       Ver.4,  PRY, 2-Feb-2009
;            Increase size of upsilon structure to 30000.
;-


pro read_upsdat_str, name, upsdat_struc, ref


; open .upsdat file...
;
openr,lur,name,/get_lun


; get no. of temperatures
;
n_t=1
readf,lur,n_t


; Create upsdat_struc
;
str={ lvl1: 0, $
      lvl2: 0, $
      de:   0., $
      gf:   0., $
      temp:   fltarr(n_t), $
      ups:   fltarr(n_t) }
upsdat_struc=REPLICATE(str,30000)

l1=1 & l2=1
;
de1=0.  &  e1=fltarr(n_t)  &  gf1=0.  &  om1=fltarr(n_t)
string1=' '
i=0
swtch=0
WHILE swtch EQ 0 DO BEGIN
  readf,lur,string1
  IF trim(string1) EQ '-1' THEN BEGIN
    swtch=swtch+1
  ENDIF ELSE BEGIN
    READS,string1,l1,l2,de1,e1,format='(2i5,70f10.6)'
    READF,lur,string1
    READS,string1,l1,l2,gf1,om1,format='(2i5,70f10.6)'
   ;
    upsdat_struc(i).lvl1=l1
    upsdat_struc(i).lvl2=l2
    upsdat_struc(i).de=de1
    upsdat_struc(i).gf=gf1
    upsdat_struc(i).temp=e1
    upsdat_struc(i).ups=om1
   ;
    i=i+1    
  ENDELSE
ENDWHILE

upsdat_struc=upsdat_struc(0:i-1)

ref=''
WHILE swtch EQ 1 DO BEGIN
  readf,lur,string1
  IF (strcompress(string1,/rem) NE '-1') THEN ref=[ref,string1] $
  ELSE swtch=2
ENDWHILE
;
ref=ref[1:*]

FREE_LUN, lur

END
