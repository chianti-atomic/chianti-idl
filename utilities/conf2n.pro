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
;	CONF2N
;
; PURPOSE:
;
;	Extract the highest principal quantum number from the configuration
;
;
; CALLING SEQUENCE:
;
;       CONF2N,conf,n
;
;
; INPUTS:
;
;	Conf:  the configuration returned from read_elvlc_direct
;
;	
; KEYWORD PARAMETERS:
;
;	None
;
; OUTPUTS:
;
;	N:  the principal quantum number
;
;
;
;
; EXAMPLE:
;
;             > conf2n,'2s2.3p 2P1.0',n
;             > print,n
;             > 3
;
; CALLS:
;
;       None
;
; RESTRICTIONS
;
;       If the principal quantum number is 10 or greater, it will not be 
;       picked up.
;
;       If the configurations are written in upper case (e.g., 3S2.3P2), 
;       then the principal quantum numbers will not be picked up (a value 
;       of 0 will be returned).
;
; MODIFICATION HISTORY:
;
;       Ver.1, April-2000, Ken Dere
;
;       Ver.2, 17-Oct-2000, Peter Young
;                removed call to str_index
;
;-

PRO conf2n,conf,n

pqn=0
spd=['s','p','d','f','g','h','i','j','k']
nspd = n_elements(spd)
for ispd=0,nspd-1 do BEGIN
   idx = -1
   ii = strpos(conf,spd[ispd])
   IF ii NE -1 THEN BEGIN
      idx = [idx,ii]
      ctest = 0
      WHILE ctest EQ 0 DO BEGIN
         ii = strpos(conf,spd[ispd],ii+1)
         IF ii EQ -1 THEN ctest = 1 ELSE idx = [idx,ii]
      ENDWHILE
      idx = idx[1:*]
   ENDIF ELSE BEGIN
      idx = -1
   ENDELSE

;   print,spd[ispd],idx
;   idx=str_index(conf,spd(ispd))
;   print,spd[ispd],idx

   IF max(idx) GE 0 THEN BEGIN
      nidx=n_elements(idx)
      FOR i=0,nidx-1 DO BEGIN
         pqn=[pqn,fix(strmid(conf,idx[i]-1,1))]
      ENDFOR
   ENDIF
ENDFOR
n=max(pqn)

END
