;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an atomic database package for the calculation of
;       continuum and emission line spectra from astrophysical plasmas. It is 
;       a collaborative project involving the Naval Research Laboratory
;       (Washington DC, USA), the Arcetri Observatory (Firenze, Italy), and the
;       Cambridge University (United Kingdom).
;
; NAME:
;       READ_GFFGU
;
; PURPOSE:
;
;       Read gffgu.dat file containing free-free gaunt factors of 
;             R. S. Sutherland, 1998, MNRAS, 300, 321
;             note:  the file available from the web site are mislabelled
;
;
; CALLING SEQUENCE:
;
;       READ_GFFGU,g2,u,gff
;
;
; INPUTS:
;
;       None    
;       
;
; OUTPUTS:
;
;       g2,u,gff defined in the paper by Sutherland
;
;
;
; COMMON BLOCKS:
;
;       None
;
;
;
; MODIFICATION HISTORY:
;       Written by:     Ken Dere
;       April 2000:     Version 3.0
;       October 2001:   Version 4 - Ken Dere -
;                       corrected for the fact that the original
;                       Sutherland file was mislabelled
;       November 2001:  Corrected the address of the gffgu.dat file - Enrico Landi
;
;-
pro read_gffgu,g2,u,gff
;
;   read gffgu.dat file of R. S. Sutherland, 1998, MNRAS, 300, 321.
;      free-free gaunt factors as a function of gamma^2 and u
;
openr,lur,!xuvtop+'/continuum/gffgu.dat',/get_lun
;
str=''
f5=fltarr(5)
ngamma=41
nu=81
g2=fltarr(ngamma)
u=fltarr(nu)
gff=fltarr(ngamma,nu)
f=fltarr(3)
;
for i=0,4 do readf,lur,str
;
for ig2=0,ngamma-1 do begin
;
      for iu=0,nu-1 do begin
;
      readf,lur,f
      g2(ig2)=f(0)
      u(iu)=f(1)
      gff(ig2,iu)=f(2)
   endfor
endfor
;
free_lun,lur
;
end
