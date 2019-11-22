FUNCTION ch_dr_a_n, n, q, t


;+
; NAME:
;      CH_DR_A_N
;
; PURPOSE:
;      Computes the function A(N) from the paper of Nikolic et
;      al. (2018, ApJS, 237,41). 
;
; CATEGORY:
;      CHIANTI; atomic data; recombination.
;
; CALLING SEQUENCE:
;      Result = CH_DR_A_N( N )
;
; INPUTS:
;      N:   The atomic number of the species for which the function is
;           required.
;
; OPTIONAL INPUTS:
;      q:   The charge on the ion. Only needed for N <= 5.
;      T:   Temperature (units: K). Can be a 1D array.
;           Only needed for N <= 5.
;
; OUTPUTS:
;      The quantity A(N) is returned, as computed using Eq. 9 from
;      Nikolic et al. (2018). A(N) is a scalar unless q and T are
;      input, in which case A will be a 1D array of same size as T.
;
; CALLS:
;      CH_PTABLE_POS, CH_DR_PSI, CH_DR_PSI_SEC.
;
; MODIFICATION HISTORY:
;      Ver.1, 28-Jul-2017, Peter Young
;      Ver.2, 7-Nov-2017, Peter Young
;         Removed the theta input, and added q and t.
;      Ver.3, 17-Nov-2017, Peter Young
;         Updated to call ch_dr_psi_sec.
;      Ver.4, 15-Jun-2020, Peter Young
;         Updated header; no change to code.
;-

IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> a_n=ch_dr_a_n(n [, q, T])'
  return,-1.
ENDIF

IF n LE 5 AND n_params() LT 3 THEN BEGIN
  print,'% CH_DR_A_N: for N <= 5, you must specify the extra inputs q and T. Please check routine header. Returning...'
  return,-1.
ENDIF 


;
; Get N1, N2 values (Eq. 13 of Nikolic et al. 2013).
;
ch_ptable_pos,n,row=row
CASE row OF
  2: BEGIN
    n1=3
    n2=10
  END
  3: BEGIN
    n1=11
    n2=18
  END
  4: BEGIN
    n1=19
    n2=36
  END
  ELSE: BEGIN  ; set arbitrary values
    n1=10
    n2=1
  END 
ENDCASE
;
; Compute A(N) - Eq. 12 of Niklic et al. (2013).
;
a=12.+10.*n1+(10.*n1-2.*n2)/(n1-n2)*(n-n1)

;
; Separately deal with N <= 5 cases.
;
IF n LE 5 THEN BEGIN
  CASE n OF
    1: a=16.
    2: a=18.
    3: a=66.
    4: a=66.
    5: a=52.
  ENDCASE 
ENDIF

IF n_elements(t) NE 0 THEN BEGIN 
  nt=n_elements(t)
  output=fltarr(nt)+a
ENDIF ELSE BEGIN
  output=a
ENDELSE 

;
; 17-Nov-2017: I updated the following to call ch_dr_psi_sec following
; latest info from D. Nikolic.
;
CASE 1 OF
  n LE 4: output=output*ch_dr_psi(n,q,t)
  n GE 5 AND n LE 6: output=output*ch_dr_psi_sec(n,q,t)
  n GE 13 AND n LE 14: output=output*ch_dr_psi_sec(n,q,t)
  ELSE: 
ENDCASE

return,output
  

END
