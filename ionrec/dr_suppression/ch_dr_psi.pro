
FUNCTION ch_dr_psi, n, q, t, simple=simple

;+
; NAME:
;      CH_DR_PSI
;
; PURPOSE:
;      Computes the Psi function used in the Nikolic et
;      al. (2018, ApJS, 237, 41) paper.
;
; CATEGORY:
;      CHIANTI; ionization balance; dielectronic recombination.
;
; CALLING SEQUENCE:
;      Result = CH_DR_PSI(N,Q,T)
;
; INPUTS:
;      N:    The index of isoelectronic sequence (1=hydrogen,
;            2=helium, etc). Must take be between 1 and 5.
;      q:    The charge on the ion.
;      T:    Temperature (units: K). Can be a 1D array.
;
; KEYWORD PARAMETERS:
;      SIMPLE: Implements the "simplified" approximation of Table 2
;              of Nikolic et al. (2018).
;
; OUTPUTS:
;      The quantity Psi^N(q,T), as defined by Nikolic et al. It has
;      the same size as the input T.
;
;      If a problem is found then -1 is returned.
;
; MODIFICATION HISTORY:
;      Ver.1, 7-Nov-2017, Peter Young
;      Ver.2, 15-Jun-2020, Peter Young
;        Updated header and comments; no change to code.
;      Ver.3, 08-Sep-2025, Peter Young
;        Added /simple keyword.
;-


IF n_params() LT 3 THEN BEGIN
  print,'Use:  IDL> psi = ch_dr_psi(N,q,T [,/simple] )'
  print,'   N - isoelectronic sequence index (1-hydrogen, etc)'
  print,'   q - charge on recombining ion'
  print,'   T - temperature (K)'
  return,-1.
ENDIF

IF n GT 5 THEN BEGIN
  message,/info,/cont,'The input N must be between 1 and 5. Returning...'
  return,-1.
ENDIF 

;
; These numbers are typed in from Table 2 of Nikolic et al. (2018).
; I'm using 1e10 instead of the infinities in the table.
;
data1=[ 4.7902, 0.32456, 0.97838, 24.78084, $
        -0.0327, 0.13265, 0.29226, 1e10, $
        -0.66855, 0.28711, 0.29083, 6.65275, $
        6.23776, 0.11389, 1.24036, 25.79559, $
        0.33302, 0.00654, 5.67945, 0.92602, $
        -0.75788, 1.75669, -0.63105, 184.82361 ]
d1=reform(data1,4,6)

data2=[ 4.82857, 0.3, 1.04558, 19.6508, $
        -0.50889, 0.6, 0.17187, 47.19496, $
        -1.03044, 0.35, 0.3586, 39.4083, $
        6.14046, 0.15, 1.46561, 10.17565, $
        0.08316, 0.08, 1.37478, 8.54111, $
        -0.19804, 0.4, 0.74012, 2.54024]
d2=reform(data2,4,6)

data3=[ 4.55441, 0.08, 1.11864, 1e10, $
        0.3, 2.0, -2.0, 67.36368, $
        -0.4, 0.38, 1.62248, 2.78841, $
        4.00192, 0.58, 0.93519, 21.28094, $
        0.00198, 0.32, 0.84436, 9.73494, $
        0.55031, -0.32251, 0.75493, 19.89169]
d3=reform(data3,4,6)

data4=[ 2.79861, 1.0, 0.82983, 18.05422, $
        -0.01897, 0.05, 1.34569, 10.82096, $
        -0.56934, 0.68, 0.78839, 2.77582, $
        4.07101, 1.0, 0.7175, 25.89966, $
        0.44352, 0.05, 3.54877, 0.94416, $
        -0.57838, 0.68, 0.08484, 6.70076]
d4=reform(data4,4,6)

data5=[ 6.75706, -3.77435, 0.0, 4.59785, $
        0.0, 0.08, 1.34923, 7.36394, $
        -0.63, 0.06, 2.65736, 2.11946, $
        7.74115, -4.82142, 0.0, 4.04344, $
        0.26595, 0.09, 1.29301, 6.81342, $
        -0.39209, 0.07, 2.27233, 1.9958]
d5=reform(data5,4,6)

CASE n OF
  1: data=d1
  2: data=d2
  3: data=d3
  4: data=d4
  5: data=d5
  ELSE: return,-1.
ENDCASE 

c=sqrt(2.5e4*q^2/t)
IF keyword_set(simple) THEN BEGIN
  psi=2.0/(1.0 + exp(-c) )
ENDIF ELSE BEGIN 
  pi=data[0,*]+data[1,*]*q^data[2,*]*exp(-q/data[3,*])

  a=(alog10(t)-pi[0])/sqrt(2.0)/pi[1]
  b=(alog10(t)-pi[3])/sqrt(2.0)/pi[4]
  
  psi=2.0 * ( 1 + pi[2]*exp(-a^2) +pi[5]*exp(-b^2) ) / $
      ( 1.0 + exp(-c) )
ENDELSE 
  
return,psi

END
