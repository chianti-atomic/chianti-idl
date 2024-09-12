;+
; PROJECT:  CHIANTI
;
;        This program was developed as part of CHIANTI-VIP. 
;        CHIANTI-VIP (Virtual IDL and Python) is  a member of the CHIANTI family
;        mantained by Giulio Del Zanna, to develop additional features and
;        provide them to the astrophysics community.
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the
;       University of Cambridge, Goddard Space Flight Center, and University of Michigan.
;
; NAME:
;
;      CH_NIKOLIC_DR_SUPPRESSION
;
; PURPOSE:
;
;      Computes the  density-dependent suppression
;      factor as specified by Nikolic et al. (2018, ApJS, 237, 41). 
;      NOTE: these factors are obtained with several approximations and 
;            are not always accurate.
;
; 
; CATEGORY:
;
;      CHIANTI; atomic data; recombination.
;
;
; CALLING SEQUENCE:
;
;      IDL> temperature=10.^(findgen(11)*0.1+5.5)
;   	 IDL> result=ch_nikolic_dr_suppression('fe_9',temperature,dens=1.e9)
;
;
; INPUTS:
;
;       Ion_Name:  The name of the recombining ion in CHIANTI
;                  format. For example, 'o_5' for O V.
;
;       Temperature: An array of temperatures (units: K).
;
;
; KEYWORD PARAMETERS:
;
;       NO_LOWT:  If set, then the low-T correction factor is not
;                 applied. See Sect. 2.3 of Nikolic et al.
;
;       QUIET:    If set, then information is not printed to the IDL
;                 window. 
; 
;
; OPTIONAL INPUTS:
;
;       DENSITY:  An electron number density (units: cm^-3) for which
;                 the suppression factor is calculated. Usually a scalar, but
;                 an array can be used if the values correspond to the temperature
;                 array.
;
;       PRESSURE: Specifies the pressure (N_e*T; units cm^-3 K). This is a scalar
;                 which will be use to create the corresponding density array.
;
;
; OUTPUTS:
;
;       S:    The suppression factor (Eq. 2 of Nikolic et al. 2018).
;
; OPTIONAL OUTPUTS:
;
;       Q0:   Q0 value (Eq. 5 of Nikolic et al. 2018).
;
;       XA:   Activation density (Eq. 3 of Nikolic et al. 2018).
;
;
; CALLS:
;
;       CH_DR_A_N, CONVERTNAME, CH_DR_EXP_EPS
;
;
; PREVIOUS HISTORY:
;
;       Modification of CH_DR_SUPPRESS written by PRY to output just the
;       suppression factor
;
;
; WRITTEN:
;         
;       Giulio Del Zanna (GDZ) and Roger Dufresne (RPD)
;       DAMTP, University of Cambridge, 14 Sept 2023 
;
;
; MODIFICATION HISTORY:
;      
;       v.2, 16 Sept 2023  RPD, changed density input to allow for an array
;
;       v.3, 27-AUG-2024  RPD, corrected error messages for density and pressure inputs
;
; VERSION     : 3
;-

FUNCTION ch_nikolic_dr_suppression, gname, temperature, density=density, $
                         q0=q0, xa=xa, quiet=quiet,filename=filename, $
                         pressure=pressure, no_lowt=no_lowt


IF n_params() LT 2 THEN BEGIN
  print,'Use:  IDL> rate=ch_dr_suppress(gname, temperature [, density=, pressure=, q0=, xa=, s='
  print,'                               /quiet, filename=, /no_lowt])'
  return,-1
ENDIF 


;
; If neither density or pressure have been specified, then just return
; the zero-density rates from ch_diel_recomb.
;
IF n_elements(density) EQ 0 AND n_elements(pressure) EQ 0 THEN $
  message,'%Error: please specify density or pressure'

IF n_elements(density) GT 0 AND n_elements(pressure) GT 0 THEN $
  message,'%Error: please specify either density or pressure'


nt=n_elements(temperature)

;
; Note that pressure over-rides density if both are input. 
;
IF n_elements(pressure) EQ 1 THEN BEGIN
  density=pressure/temperature
  nd=nt
ENDIF ELSE IF n_elements(pressure) GT 1 THEN BEGIN
  message,'%Error: only a single, constant pressure can be specified'
ENDIF ELSE IF n_elements(density) GT 1 THEN BEGIN
  IF n_elements(density) EQ nt THEN nd=nt ELSE $
    message,'%Error: density array is not the same length as temperature array'
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
w=5.64548
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
; Return the suppression factor.
;
return, s

END
