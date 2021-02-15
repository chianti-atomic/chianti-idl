
FUNCTION karzas_gfb, wvl, pqn, level, ioniz_en

;+
; NAME:
;     KARZAS_GFB
;
; PURPOSE:
;     Gives the free-bound Gaunt factor, interpolated from the Tables
;     of Karzas & Latter (1961).
;
; CATEGORY:
;     CHIANTI; continuum; freebound.
;
; CALLING SEQUENCE:
;     Result = KARZAS_GFB( Wvl, Pqn, Level, Ioniz_En )
;
; INPUTS:
;     Wvl:  1D array of wavelengths in angstroms.
;     Pqn:  Principal quantum number (n). Integer from 1 to 6.
;     Level: Orbital angular momentum for orbital. Integer from 0
;            to PQN-1. For example, level=0 corresponds to 's', 1 to
;            'p', 2 to 'd'.
;     Ioniz_En:  Energy in cm^-1 to ionize the ion from the specified
;                shell. 
;
; OUTPUTS:
;     A 1D array of same size as WVL containing the free-bound Gaunt
;     factor. If a problem is found, then a value of -1 is returned.
;
; CALLS:
;     READ_KLGFB
;
; EXAMPLE:
;     Get Gaunt factor for the Ne V 2s2.2p.3s configuration
;
;     IDL> zion2filename,10,5,fname
;     IDL> read_fblvl,fname+'.fblvl',l1,conf,pqn,ll,spd,mult,ecm,ecmth,ref
;     IDL> wvl=findgen(100)+1.
;     IDL> ip=ch_ip('ne_5',/cm)
;     IDL> gfb=karzas_gfb(wvl,pqn[1],ll[1],ip-ecm[1])
;
; MODIFICATION HISTORY:
;     Ver.1, 12-Feb-2021, Peter Young
;-


IF n_params() LT 4 THEN BEGIN
   print,'Use:  IDL> gfb=karzas_gfb( wvl, pqn, level, ioniz_en )'
   return,-1
ENDIF 


IF pqn GT 6 THEN BEGIN
   print,'% KARZAS_GFB: The Karzas data are only tabulated for PQN < 7. Returning...'
   return,-1.
ENDIF 

read_klgfb,pe,klgfb,pqn
IF n_elements(pe) EQ 0 THEN return,-1.

nwvl=n_elements(wvl)

IF level GE pqn THEN BEGIN
   print,'% KARZAS_GFB: The input LEVEL must be less than or equal to PQN. Returning...'
   return,-1.
ENDIF 

;
; The original data is in reverse energy ordering, so I need to
; reverse it here.
;
; The Karzas tabulated energies are for a purely Rydberg set of
; levels, and the minimum energy corresponds to the ionization energy
; from the N level. For CHIANTI we will use actual ionization
; energies. For this reason EN is normalized to the ionization energy
; so that the minimum energy is 1
;
gfb=reverse(reform(klgfb[*,level]))
en=reverse(pe)/min(pe)


;
; I need to convert WVL to cm^-1 and then normalize by ionization energy.
;
wvl_norm=reverse(1d8/wvl)/ioniz_en

output=dblarr(nwvl)


;
; Only do interpolation over tabulated range of data, and for energies
; above the ionization edge.
;
k=where(wvl_norm GE min(en) AND wvl_norm LE max(en) AND wvl_norm GE 1.0,nk)
IF nk NE 0 THEN BEGIN 
   x=alog10(en)
   y=alog10(gfb)
   xi=alog10(wvl_norm[k])
  ;
   y2=spl_init(x,y)
   yi=spl_interp(x,y,y2,xi)
  ;
   output[k]=10.^yi
ENDIF


;
; For energies greater than tabulated range, I perform a linear fit to
; the five highest energy points of the data. I then use the fit to
; obtain the gfb values.
;
k=where(wvl_norm GT max(en),nk)
IF nk GT 0 THEN BEGIN
   n=n_elements(en)
   x=alog10(en[n-5:n-1])
   y=alog10(gfb[n-5:n-1])
   c=linfit(x,y,/double)
   yi=c[0]+c[1]*alog10(wvl_norm[k])
   output[k]=10.^yi
ENDIF 

output=reverse(output)

return,output

END
