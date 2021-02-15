
FUNCTION karzas_xs, wvl, n, l, ioniz_en, z, pe=pe, klgfb=klgfb

;+
; NAME:
;     KARZAS_XS()
;
; PURPOSE:
;     Outputs the photoionization cross-section at the wavelengths WVL for 
;     ionization of an N,L electron with ionization energy IONIZ_EN, 
;     calculated using the Karzas & Latter (ApJS 6, 167, 1961) formulation. 
;     The bound-free gaunt factor is derived from the tables of Karzas & 
;     Latter (1961).
;
; CATEGORY:
;     CHIANTI; continuum; freebound.
;
; CALLING SEQUENCE:
;     Result = KARZAS_XS( Wvl, N, L, Ioniz_En )
;
; INPUTS:
;     Wvl:   A 1D array of wavelengths for which cross-section is
;            required. 
;     N:     Principal quantum number of the electron being removed in the 
;            ionization. Integer between 1 and 6.
;     L:     The orbital angular momentum of the electron being removed in the 
;            ionization. Integer between 0 and N-1.
;     Ioniz_En: The ionization energy (in cm^-1) of the electron being 
;               removed.
;
; OPTIONAL INPUTS:
;     Pe:    This input is now obsolete and does not do anything. Kept
;            for backwards compatibility.
;     Klgfb: This input is now obsolete and does not do anything. Kept
;            for backwards compatibility.
;
; OUTPUTS:
;     The photoionization cross-section for the ionization of the outer 
;     electron in units of mega-barns (Mb = 10^-18 cm^2) at the input 
;     wavelengths. E.g., for Fe XIII (ground configuration 
;     1s2.2s2.2p6.3s2.3p2) it is the cross-section for the removal of the 
;     3p electron.
;
; CALLS:
;     KARZAS_GFB
;
; EXAMPLE:
;     Compute cross-section for the first excited configuration of Fe
;     XIII.
;
;     IDL> zion2filename,26,13
;     IDL> read_fblvl,fname+'.fblvl',l1,conf,pqn,ll,spd,mult,ecm,ecmth,ref
;     IDL> wvl=findgen(100)+1.
;     IDL> ip=ch_ip('fe_13',/cm)
;     IDL> xs=karzas_xs(wvl,pqn[1],ll[1],ip-ecm[1])
;
; MODIFICATION HISTORY:
;     Ver.1, 24-Jul-2002, Peter Young
;     Ver.2, 12-Feb-2021, Peter Young
;-


IF n_params() LT 4 THEN BEGIN
   print,'Use:  IDL> gfb=karzas_xs( wvl, pqn, level, ioniz_en )'
   return,-1
ENDIF 

;
; Convert wavelength to energy in cm^-1
;
ewvl=1d8/wvl

gfb=karzas_gfb(wvl,n,l,ioniz_en)

xs=1.077294d-1*8065.54d3*(ioniz_en)^2*gfb/n/ewvl^3

return,xs

END
