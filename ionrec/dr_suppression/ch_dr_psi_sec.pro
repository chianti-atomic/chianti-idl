
FUNCTION ch_dr_psi_sec, n, q, t

;+
; NAME:
;      CH_DR_PSI
;
; PURPOSE:
;      Computes the Psi_sec function used in the Nikolic et
;      al. (2018, ApJS, 237, 41) paper.
;
; CATEGORY:
;      CHIANTI; ionization balance; dielectronic recombination.
;
; CALLING SEQUENCE:
;      Result = CH_DR_PSI_SEC(N,Q,T)
;
; INPUTS:
;      N:    The index of the isoelectronic sequence. Must take a
;            value of 5, 6, 13 or 14.
;      q:    The charge on the ion.
;      T:    Temperature (units: K). Can be a 1D array.
;
; OUTPUTS:
;      The quantity Psi_sec^N(q,T), as defined by Nikolic et al. It has
;      the same size as the input T.
;
;      If a problem is found then -1 is returned.
;
; MODIFICATION HISTORY:
;      Ver.1, 17-Nov-2017, Peter Young
;      Ver.2, 15-Jun-2020, Peter Young
;        Updated header and comments; no change to code.
;-


IF n_params() LT 3 THEN BEGIN
  print,'Use:  IDL> psi_sec = ch_dr_psi_sec(N,q,T)'
  print,'   N - isoelectronic sequence index (1-hydrogen, etc)'
  print,'   q - charge on recombining ion'
  print,'   T - temperature (K)'
  return,-1.
ENDIF

IF (n LT 5) OR (n GT 6 AND n LT 13) OR (n GT 14) THEN BEGIN
  print,'% CH_DR_PSI_SEC: the input N must be 5, 6, 13 or 14. Returning...'
  return,-1.
ENDIF 

;
; These numbers are typed in from Table 2 of Nikolic et al. (2018)
;
data1=[ 6.91078, -1.6385, 2.18197, 1.45091, $
        0.4959, -0.08348, 1.24745, 8.55397, $
        -0.27525, 0.132, 1.15443, 3.79949, $
        7.45975, -2.6722, 1.7423, 1.19649, $
        0.51285, -0.60987, 5.15431, 0.49095, $
        -0.24818, 0.125, 0.59971, 8.34052 ]
d1=reform(data1,4,6)

data2=[ 5.90184, -1.2997, 1.32018, 2.10442, $
        0.12606, 0.009, 8.33887, 0.44742, $
        -0.28222, 0.018, 2.50307, 3.83303, $
        6.96615, -0.41775, 2.75045, 1.32394, $
        0.55843, 0.45, 0.0, 2.06664, $
        -0.17208, -0.17353, 0.0, 2.57406]

d2=reform(data2,4,6)

data3=[ 6.59628, -3.03115, 0.0, 10.519821, $
        1.20824, -0.85509, 0.21258, 25.56, $
        -0.34292, -0.06013, 4.09344, 0.90604, $
        7.92025, -3.38912, 0.0, 10.02741, $
        0.06976, 0.6453, 0.24827, 20.94907, $
        -0.34108, -0.17353, 0.0, 6.0384]
d3=reform(data3,4,6)

data4=[ 5.54172, -1.54639, 0.01056, 3.24604, $
        0.39649, 0.8, 3.19571, 0.642068, $
        -0.35475, -0.08912, 3.55401, 0.73491, $
        6.88765, -1.93088, 0.23469, 3.23495, $
        0.58577, -0.31007, 3.30137, 0.83096, $
        -0.14762, -0.16941, 0.0, 18.53007]
d4=reform(data4,4,6)


CASE n OF
  5: data=d1
  6: data=d2
  13: data=d3
  14: data=d4
  ELSE: return,-1.
ENDCASE 



gamma=data[0,*]+data[1,*]*q^data[2,*]*exp(-q/data[3,*])

a=(alog10(t)-gamma[0])/sqrt(2.0)/gamma[1]
b=(alog10(t)-gamma[3])/sqrt(2.0)/gamma[4]

psi_sec= 1 + gamma[2]*exp(-a^2) +gamma[5]*exp(-b^2) 

return,psi_sec

END
