;+
;
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
;
; NAME:
;
;      CH_BURGCHID_RATE
;
;
; PURPOSE:
;
;      Calculate ionisation rate coefficients from 
;      Burgess and Chidichimo (1983) formula and return the rates
;      at the specified temperatures.
;
;
; CATEGORY:
;
;      CHIANTI; atomic data; electron impact ionisation.
;
;
; CALLING SEQUENCE:
;
;      For ionisation rate coefficient from n=2 shell of neutral nitrogen:
;      IDL> logt=findgen(41)/10.+3.0
;      IDL> t=10.^logt
;      IDL> rate=burgchid_rate(0,3,14.534,t)+burgchid_rate(0,2,20.340,t)
;
;
; INPUTS:
;
;      Charge:    electric charge of the ion, e.g., O III has electric
;                 charge of +2.
;
;      Eff_elec:  Effective number of electrons in shell being ionised,
;                 e.g., ionisation of 2p shell in 1s2 2s2 2p3 has 3
;                 effective electrons.
;
;      Shell_ip:  Ionisation potential in eV of shell being ionised.
;
;      Temp:      The temperature(s) in K at which the rates are
;                 required.
;
;
; KEYWORD PARAMETERS:
;
;      NONE
; 
;
; OPTIONAL INPUTS:
;
;      Burg_c:    Value of the constant used in the Burgess and Chidichimo
;                 formula. This can be used to fit the cross section to a known
;                 value. Otherwise, this will be set to the default value 2.3
;                 given in the paper. 
;
;
; OUTPUTS:
;
;      The ionisation rate in cm^3 s^-1 at the specified temperatures.
;
;
; OPTIONAL OUTPUTS:
;
;       NONE
;
;
; CALLS:
;
;       NONE
;
;
; PREVIOUS HISTORY:
;
;       NONE
;
;
; WRITTEN:
;         
;       Ver.1, 07-Aug-2023, Roger Dufresne
;
;
; MODIFICATION HISTORY:
;
;       NONE
;
;
; VERSION:  1
;
;- 

function ch_burgchid_rate, chg, eff_elec, shell_ip, temp, burg_c=burg_c


k_ev=8.61734d-5
ryd=13.605698d0

if n_elements(burg_c) eq 0 then burg_c=2.3d0

wannier_beta=0.25d0*(sqrt((100.0d0*chg+91.0d0)/(4.0d0*chg+3.d0))-5.0d0)

reduced_t=k_ev*temp/shell_ip
wannier_rt=(alog(1.0d0+reduced_t))^(wannier_beta/(1.0d0 + reduced_t))

rate=2.1715d-8*burg_c*eff_elec*sqrt((ryd/shell_ip)^3.0d0)* $
  expint(1,1.0d0/reduced_t)*wannier_rt/sqrt(reduced_t)


return,rate


end
