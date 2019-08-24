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
; NAME:
;	READ_GFFINT
;
; PURPOSE:
;
;	Read gffint file containing integrated free-free gaunt factors of 
;             R. S. Sutherland, 1998, MNRAS, 300, 321
;
;
; CALLING SEQUENCE:
;
;       READ_GFFINT,g2,gff,s1,s2,s3
;
;
; INPUTS:
;
;	None	
;
;	
;
; OUTPUTS:
;
;	g2,gff,s1,s2,s3 defined in the paper by Sutherland
;
;

;
; COMMON BLOCKS:
;
;	None
;
;
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	April 2000:     Version 3.0
;
;-
pro read_gffint,g2,gff,s1,s2,s3
;
;   read gffint.dat file of R. S. Sutherland, 1998, MNRAS, 300, 321.
;      integrated free-free gaunt factors 
;
openr,lur,!xuvtop+'/continuum/gffint.dat',/get_lun
;
str=''
f5=fltarr(5)
ngamma=41
g2=fltarr(ngamma)
gff=g2
s1=g2
s2=g2
s3=g2
;
for i=0,3 do readf,lur,str
;
for i=0,40 do begin
   readf,lur,f5
   g2(i)=f5(0)
   gff(i)=f5(1)
   s1(i)=f5(2)
   s2(i)=f5(3)
   s3(i)=f5(4)
endfor
;
free_lun,lur
;
end
