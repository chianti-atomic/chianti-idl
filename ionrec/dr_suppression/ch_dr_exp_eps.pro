
FUNCTION ch_dr_exp_eps, n, q, x, t

;+
; NAME:
;      CH_DR_EPS
;
; PURPOSE:
;      Computes the exponential term in  Eq. 14 of
;      Nikolic et al. (2018, ApJS, 237, 41). The routine is intended 
;      to be called from CH_DR_SUPPRESS.
;
; CATEGORY:
;      CHIANTI; ionization balance; dielectronic recombination
;
; CALLING SEQUENCE:
;	Result = CH_DR_EXP_EPS(N,Q)
;
; INPUTS:
;      N:    The index of isoelectronic sequence (1=hydrogen,
;            2=helium, etc).
;      q:    The charge on the ion.
;      X:    Logarithm of electron number density (units: cm^-3). Can
;            be a scalar or a 1D array. For the latter, it must be the
;            same size as T.
;      T:    Temperature (units: K). Can be a 1D array.
;
; KEYWORD PARAMETERS:
;      None.
;
; OUTPUTS:
;      The quantity exp (- epsilon_N(q)/10kT) from Eq. 14 of Nikolic
;      et al. (2013, ApJ, 768, 82). Note that if this quantity is 0,
;      then there is no suppression, whereas if it is 1 then
;      suppression is unmodified from the high temperature case. The
;      output is a 1D array of same size as the input T.
;
;      If a problem is found then the scalar '-1' is returned.
;
; MODIFICATION HISTORY:
;      Ver.1, 6-Nov-2017, Peter Young
;      Ver.2, 15-Jun-2020, Peter Young
;        Updated header; no change to code.
;      Ver.3, 05-Sep-2025, Peter Young
;        Added special cases for N=1,2,10 and 14.
;-

IF n_params() LT 4 THEN BEGIN
  print,'Use:  IDL> output=ch_dr_exp_eps(n,q,x,t)'
  return,-1.
ENDIF 

nx=n_elements(x)
nt=n_elements(t)

IF nx NE 1 AND nt NE nx THEN BEGIN
  message,/info,/cont,'The input X must be either a scalar or an array with the same size as T. Returning value of 1.'
  return,1.0
ENDIF 

kt=8.61735e-5*t   ; in eV

output=fltarr(nt)

;
; Load the coefficients p_j for select isoelectronic sequences. For
; the other sequences the p_j are zero.
;
p=fltarr(6)
CASE n OF
  3: p=[1.963e0,2.030e1,-9.710e-1,8.545e-1,1.355e-1,2.401e-2]
  4: p=[5.789e0,3.408e1,1.517e0,-1.212e0,7.756e-1,-4.100e-3]
  7: p=[1.137e1,3.622e1,7.084e0,-5.168e0,2.451e0,-1.696e-1]
  11: p=[2.248e0,2.228e1,-1.123e0,9.027e-1,-3.860e-2,1.468e-2]
  12: p=[2.745e0,1.919e1,-5.432e-1,7.868e-1,-4.249e-2,1.357e-2]
  14: p=[17.6874,0.,0.,0.,0.,0.]
  15: p=[1.428e0,3.908e0,7.312e-1,-1.914e0,1.051e0,-8.992e-2]
  ELSE: 
ENDCASE 

;
; For Si-like sequence (N=14), only S2+ is non-zero.
;
IF n EQ 14 AND q NE 2 THEN p[0]=0.

;
; Special case for H, He and Li-like ions - see Table 5 of Nikolic
; et al. (2018).
; Note that xa0 is fixed for all sequences (Sect. 2 of Nikolic et
; al. 2018).
;
IF n EQ 1 OR n EQ 2 OR n EQ 10 THEN BEGIN
  xa0=10.1821
  p[0]=20.*erfc(2.*(x-xa0))
ENDIF 

eps=0.0
FOR i=0,5 DO eps=eps+p[i]*(q/10.0)^i

output=exp( - eps / (10.*kt) )
return,output

END
