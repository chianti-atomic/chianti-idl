
;+
; NAME:
;      CH_DR_SUPPRESS
;
; PURPOSE:
;      Computes the total dielectronic recombination rate for the
;      specified ion, including the density-dependent suppression
;      factor as specified by Nikolic et al. (2018, ApJS, 237, 41). 
;
; CATEGORY:
;      CHIANTI; atomic data; recombination.
;
; CALLING SEQUENCE:
;	Result = CH_DR_SUPPRESS( ion_name, temperature )
;
; INPUTS:
;       Ion_Name:  The name of the recombining ion in CHIANTI
;                  format. For example, 'o_5' for O V.
;       Temperature: An array of temperatures (units: K).
;
; OPTIONAL INPUTS:
;       Density:  An electron number density (units: cm^-3) for which
;                 the suppression factor is calculated.
;
; KEYWORD PARAMETERS:
;       NO_LOWT:  If set, then the low-T correction factor is not
;                 applied. See Sect. 2.3 of Nikolic et al.
;       QUIET:    If set, then information is not printed to the IDL
;                 window. 
; 
; OUTPUTS:
;       The dielectronic recombination rate (units: cm^3 s^-1) for the
;       specified recombining ion tabulated for specified temperatures.
;
; OPTIONAL OUTPUTS:
;       S:    The suppression factor (Eq. 2 of Nikolic et al. 2018).
;       Q0:   Q0 value (Eq. 5 of Nikolic et al. 2018).
;       XA:   Activation density (Eq. 3 of Nikolic et al. 2018).
;
; CALLS:
;       CH_DIEL_RECOMB, CH_DR_A_N, CONVERTNAME, CH_DR_EXP_EPS
;
; EXAMPLE:
;       IDL> log_temp=findgen(101)/20.+4.0
;       IDL> rate=ch_dr_suppress('si_4',10.^log_temp)
;
; MODIFICATION HISTORY:
;      Ver.1, 25-Apr-2017, Peter Young
;      Ver.2, 27-Jul-2017, Peter Young
;        If neither pressure or density are provided, then just return
;        the zero-density rates.
;      Ver.3, 7-Nov-2017, Peter Young
;        Updated call to ch_dr_a_n, and added call to ch_dr_exp_eps
;        (for low temperature modification to DR rates).
;      Ver.4, 12-Jun-2020, Peter Young
;        Updated header, changing references to Nikolic et al. (2018)
;        instead of Nikolic et al. (2013); code unchanged.
;      Ver.5, 08-Sep-2025, Peter Young
;        Updated value of w to 5.64586 from 5.64548. The latter was
;        from the 2013 paper, the former is from the 2018 paper.
;-



FUNCTION ch_dr_suppress, gname, temperature, density=density, $
                         q0=q0, xa=xa, s=s, quiet=quiet,filename=filename, $
                         pressure=pressure, no_lowt=no_lowt


IF n_params() LT 2 THEN BEGIN
  print,'Use:  IDL> rate=ch_dr_suppress(gname, temperature [, density=, pressure=, q0=, xa=, s='
  print,'                               /quiet, filename=, /no_lowt])'
  return,-1
ENDIF 

;
; Get zero density rates
;
rate0=ch_diel_recomb(gname,temperature,quiet=quiet,filename=filename)
IF rate0[0] EQ -1 THEN return,-1

;
; If neither density or pressure have been specified, then just return
; the zero-density rates from ch_diel_recomb.
;
IF n_elements(density) EQ 0 AND n_elements(pressure) EQ 0 THEN return,rate0


nt=n_elements(temperature)

;
; Note that pressure over-rides density if both are input. 
;
IF n_elements(pressure) NE 0 THEN BEGIN
  density=pressure/temperature
  nd=nt
ENDIF ELSE BEGIN
  nd=1
ENDELSE 

;
; Compute activation density
;
convertname,gname,iz,ion
q=ion-1     ; charge
N=iz-ion+1  ; iso-electronic sequence number (hydrogen=1)
;theta=temperature/q^2

;
; q0 - Eq. 5 of Nikolic et al. (2013).
;
; Note that q0 is an array of same size as temperature
;
q0=( 1.-sqrt(2./(3*q)) )*ch_dr_a_n(n,q,temperature)/sqrt(q)

;
; T0 - Eq. 4 of Nikolic et al. (2013).
;
t0=5e4*q0^2

;
; Activation density - Eq. 3 of Nikolic et al. (2013).
;
xa0=10.1821
w=5.64586
xa=xa0+alog10( (q/q0)^7 * sqrt(temperature/t0) )

;
; Define x
;
x=alog10(density)

;
; Compute suppression factor, S, using Eq. 2 of Nikolic et
; al. (2013). 
;
IF nd EQ 1 THEN BEGIN 
  s=dblarr(nt)+1.0
  k=where(x GT xa,nk)
  IF nk GT 0 THEN s[k]=exp( - (x-xa[k])^2/w^2*alog(2.) )
ENDIF ELSE BEGIN
  s=dblarr(nt)+1.0
  k=where(x GT xa,nk)
  IF nk GT 0 THEN s[k]=exp( - (x[k]-xa[k])^2/w^2*alog(2.) )
ENDELSE 

;
; Now compute the modification for low temperatures (Eq. 14 of
; Nikolic).
;
IF NOT keyword_set(no_lowt) THEN BEGIN 
  exp_eps=ch_dr_exp_eps(n,q,x,temperature)
  s=1.0 - (1.0-s)*exp_eps
ENDIF 

;
; Return the suppressed rate.
;
return,rate0*s

END
