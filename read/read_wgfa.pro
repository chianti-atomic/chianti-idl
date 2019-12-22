

pro read_wgfa,filename,lvl1,lvl2,wvl,gf,a_value,ref

;+
; PROJECT:  CHIANTI
;
;       CHIANTI is an atomic database package for the calculation of
;       continuum and emission line spectra from astrophysical plasmas.
;
; NAME:
;	READ_WGFA
;
; PURPOSE:
;
;	To read CHIANTI .wgfa files, which contain wavelength,
;	oscillator strengths and radiative decay rates.
;
; CATEGORY:
;
;	Science
;
; CALLING SEQUENCE:
;
;       READ_WGFA, File, Lvl1, Lvl2, Wvl, Gf, A_value, Ref
;
;
; INPUTS:
;
;	File:  name of the file containing the radiative data
;                i.e. !xuvtop/c/c_4/c_4.wgfa
;
;
; OUTPUTS:
;
;	Lvl1:  1D array of indices of the lower level (starting at 1)
;       Lvl2:  1D array of indices of the upper level (starting at 1)
;       Wvl:   2D array of transition wavelengths in Angstroms
;       Gf:    2D array of weighted oscillator strength gf
;       A_value:  2D array of radiative transition probability (s^-1)
;       Ref:   1D string array of references to the data in the scientific literature
;
;       Note that the 2D arrays have size (N,N) where N is the maximum
;       value of [L1,L2].  If the user is interested in the radiative
;       decay 20->1, then the value will be stored in A_value[0,19].
;
; EXAMPLE:
;
;       > read_wgfa,!xuvtop+'/c/c_4/c_4.wgfa',lvl1,lvl2,wvl,gf,a,ref
;             
;
; CALLS
;
;       READ_WGFA2
;
; MODIFICATION HISTORY:
;
; 	Written by:	Ken Dere
;	March 1996:     Version 2.0
;
;       v.3, 6-Mar-02, Peter Young
;           Changed the way the reference string is read in order to 
;           prevent '-1' problems.
;
;       v.4, 12-Mar-02, Peter Young
;           Corrected bug following above change.
;       
;       v.5, 11-Feb-2012, Peter Young
;           This routine has been completely revised so that it simply
;           calls read_wgfa2 and then re-formats the output to match
;           the old read_wgfa routine.
;-

if n_params() lt 6 then begin
   print,' '
   print,'   IDL> read_wgfa,filename,lvl1,lvl2,wvl,gf,a_value,ref'
   print,'  i.e.> read_wgfa,!xuvtop+''/c/c_4/c_4.wgfa'',lvl1,lvl2,wvl,gf,a,ref'
   print,'         to read the gf and A values for C IV'
   print,' '
   return
endif

print,'%READ_WGFA: This routine is now obsolete. Please consider using READ_WGFA2.'

read_wgfa2,filename,lvl1,lvl2,ww,gg,aa,ref

n=max([lvl1,lvl2])

wvl=fltarr(n,n)
gf=fltarr(n,n)
a_value=fltarr(n,n)

nl=n_elements(lvl1)

FOR i=0,nl-1 DO BEGIN
  lev1=lvl1[i]
  lev2=lvl2[i]
  wvl[lev1-1,lev2-1]=ww[i]
  gf[lev1-1,lev2-1]=gg[i]
  a_value[lev1-1,lev2-1]=aa[i]
ENDFOR

END
