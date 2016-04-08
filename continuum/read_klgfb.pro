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
;
; NAME:
;	READ_KLGFB
;
; PURPOSE:
;
;	to read CHIANTI files file containing the free-bound gaunt factors for n=1-6
;       from Karzas and Latter, 1961, ApJSS, 6, 167
;
; CATEGORY:
;
;	science.
;
; CALLING SEQUENCE:
;
;       READ_KLGFB,PE,GFB
;
;
; INPUTS:
;
;	File:  none
;
;
;	
;
; OUTPUTS:
;
;	PE:   an array of natural log of 41 values of Photon Energy values relative to 
;             the recombination edge i.e., at the recombination edge, PE=0.
;       GFB:  an array containing the natural log of free-bound gaunt factors indexed 
;             by energy, n and l
;
;
;
;  PROCEDURE:  reads the file:  '!xuvtop/continuum/klgfb.dat'

;   
;
; EXAMPLE:
;
;             > read_klgfb,pe,klgfb
;                  the free-bound gaunt factor vs. energy for recombination to a hydrogen 2p state
;                  (n=2 and l=1) is given by klgfb(*,n-1,l-1)
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	July 2002:     Version 1.0
;             Version 2, 8-Aug-02, GDZ
;             corrected a typo.
;
; VERSION     :   V.2, 8-Aug-02
;
;-
pro read_klgfb,pe,gfb
;
; read karzas and latter files containing log(photon energy),log (fb gaunt factor)
;
openr,lur,concat_dir(concat_dir(!xuvtop,'continuum'),'klgfb.dat'),/get_lun
;
ngfb=1
nume=1
;
readf,lur,ngfb,nume    ; number of n values, number of photon energy values
;
pe=fltarr(nume)
readf,lur,pe
;
g=fltarr(nume)
gfb=dblarr(nume,ngfb,ngfb)
;
i1=1 & i2=1
;
for n=0,ngfb-1 do begin
   for l=0,n do begin
      readf,lur,i1,i2,g
      gfb(0,n,l)=g
   endfor
endfor
;
free_lun,lur
;
end


