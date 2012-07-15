
FUNCTION karzas_xs, wvl, n, l, ioniz_en, pe=pe, klgfb=klgfb

;+
; PROJECT:  CHIANTI
;
;      CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;      Astrophysical Plasmas. It is a collaborative project involving the Naval
;      Research Laboratory (USA), the University of Florence (Italy), the
;      University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
; NAME
;
;    KARZAS_XS()
;
; EXPLANATION
;
;    Outputs the photoionization cross-section at the wavelengths WVL for 
;    ionization of an N,L electron with ionization energy IONIZ_EN, 
;    calculated using the Karzas & Latter (ApJS 6, 167, 1961) formulation. 
;    The bound-free gaunt factor is derived from the tables of Karzas & 
;    Latter.
;
; INPUTS
;
;    WVL    Wavelengths (in angstroms) for which cross-sections are 
;           required (1-D array).
;
;    N      Principal quantum number of the electron being removed in the 
;           ionization.
;
;    L      The orbital angular momentum of the electron being removed in the 
;           ionization.
;
;    IONIZ_EN  The ionization energy (in cm^-1) of the electron being 
;              removed.
;
; OPTIONAL INPUTS
;
;    PE      Allows the PE array from READ_KLGFB to be directly input to 
;            KARZAS_XS thus avoiding the need to read the K&L data 
;            repeatedly for many ions. Requires KLGFB to also be input.
;
;    KLGFB   Allows the KLGFB array from READ_KLGFB to be directly input to 
;            KARZAS_XS thus avoiding the need to read the K&L data 
;            repeatedly for many ions. Requires PE to also be input.
;
; OUTPUT
;
;    The photoionization cross-section for the ionization of the outer 
;    electron in units of mega-barns (Mb = 10^-18 cm^2) at the input 
;    wavelengths. E.g., for Fe XIII (ground configuration 
;    1s2.2s2.2p6.3s2.3p2) it is the cross-section for the removal of the 
;    3p electron.
;
; CALLS
;
;    READ_KLGFB
;
; HISTORY
;
;    Ver.1, 24-Jul-2002, Peter Young
;-

IF n_elements(pe) EQ 0 AND n_elements(klgfb) EQ 0 THEN read_klgfb,pe,klgfb

ksize=size(klgfb)
max_ngf=ksize(2)

ewvl=1d8/wvl

;
; fbgf is the free-bound gaunt factor. Note that if the N,L are not covered
; by the K&L table, then the gaunt factor is set to 1.
;
IF n LE max_ngf THEN BEGIN
  wrat=alog(ewvl/ioniz_en)>0.
  klgfb_i=klgfb[*,n-1,l]
  y2=spl_init(pe,klgfb_i)
  fbgf=spl_interp(pe,klgfb_i,y2,wrat)
  fbgf=exp(fbgf)
ENDIF ELSE BEGIN
  fbgf=replicate(1.d,n_elements(wvl))
ENDELSE

xs=1.077294d-1*8065.54d3*(ioniz_en)^2*fbgf/n/ewvl^3
ind=where(ewvl LT ioniz_en)
IF ind[0] NE -1 THEN xs[ind]=0d0

return,xs

END
