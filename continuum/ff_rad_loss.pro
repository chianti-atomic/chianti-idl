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
;           FF_RAD_LOSS
;	
;
; PURPOSE:
;
;	Calculate the free-free radiative energy losses losses.
;       Uses the free-free integrated gaunt factor calculations of 
;       Sutherland, 1998, MNRAS, 300, 321
;
;
; CALLING SEQUENCE:
;
;      FF_RAD_LOSS,Temperature,LossRate
;
;
; INPUTS:
;
;       None
;
; OPTIONAL INPUTS:
;	Abund_File:  Specifies the element abundance file to be
;                    used.
;       Ioneq_File:  Specifies the ionization equilibrium file.
;	
; KEYWORD PARAMETERS:
;
;	NO_SETUP:   If the procedure setup_elements has already been called then
;                   the keyword /no_setup should be set to avoid
;                   repeating this step 
;
;       QUIET:      If set, then no information messages are printed.
;
;       MIN_ABUND:  If set, calculates the continuum only from those elements which 
;                   have an abundance greater than min_abund.  Can speed up the 
;                   calculations.  For example:
;                   abundance (H)  = 1.
;                   abundance (He) = 0.085
;                   abundance (C)  = 3.3e-4
;                   abundance (Si) = 3.3e-5
;                   abundance (Fe) = 3.9e-5
;
;
; OUTPUTS:
;
;	Temperature:	Temperatures in degrees Kelvin. Taken from the
;                       ionization equilibrium file selected by user.
;       LossRate:       Radiative energy loss rate in erg s^-1 cm^3
;                       (radiative loss rate per emission measure (N_e
;                       N_H V) calculated for the set of temperatures.
;
;
;
; COMMON BLOCKS:
;
;	common elements,abund,abund_ref,ioneq,ioneq_t,ioneq_ref
;
;
;
; EXAMPLE:
;             > ff_rad_loss,t,rad
;             > ff_rad_loss,t,rad,min_abund=3.e-5
;             > ff_rad_loss,t,rad,/no_setup,min_abund=1.e-6
;             > ff_rad_loss,t,rad,abund_file=!abund_file,ioneq_file=!ioneq_file.
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	April 2000:     Version 3.0
;
;       Ver.2, 8-Aug-2017, Peter Young
;           added abund_file and ioneq_file optional inputs; added
;           /quiet and check on no_setup
;
;-
pro ff_rad_loss,t,rad_loss,no_setup=no_setup,min_abund=min_abund, $
                abund_file=abund_file, ioneq_file=ioneq_file, quiet=quiet
;
;  to calculate the total free-free radiative losses vs. temperature
;  temperature array determined by ionization equilibrium calculations
;
common elements,abund,abund_ref,ioneq,ioneq_t,ioneq_ref
;
;
;
if n_params() lt 2 then begin
   print,' '
   print,' type> ff_rad_loss,temperature,loss_rate [,/no_setup, min_abund= '
   print,'                    abund_file=, ioneq_file= ]'
   print,' '
   return
endif
;
kb=1.38062d-16   ;  erg deg-1
ryd=2.17992d-11  ; erg
factor=1.42554d-27
;
;
; Read elemental abundances, ionization equilibrium.
;  PRY, 8-Aug-2017: if /no_setup is set but abund is empty then I use
;  setup_elements to prevent a crash.
;
chck=n_elements(abund)
if not keyword_set(no_setup) OR chck EQ 0 then begin
  IF keyword_set(no_setup) AND chck EQ 0 THEN BEGIN
    IF NOT keyword_set(quiet) THEN print,'%FF_RAD_LOSS: /no_setup was used, but common block was empty.'
  ENDIF 
  setup_elements, abund_file=abund_file, ioneq_file=ioneq_file
endif
;
   read_ip,!xuvtop+'/ip/chianti.ip',ionpot,ipref
;
zmax=max(where(abund gt 0.))+1
;   print,zmax,abund(zmax-1)
;
n_ioneq_t=n_elements(ioneq_t)
;
read_gffint,g2,gff,s1,s2,s3
;
rad_loss=fltarr(n_ioneq_t)
;
;
;
for iz=1,zmax do begin
   this_abund=abund(iz-1)
;   print,iz,'  this_abund = ',this_abund
;
    for ion=2,iz+1 do begin
        ;  iz=1 & ion=2
        z=float(ion-1)
        dielectronic=0
        this_ioneq=ioneq(*,iz-1,ion-1+dielectronic)
        ip=ionpot(iz-1,ion-2)
        ;
        ;print,' iz, ion = ',iz,ion,' max ioneq = ',max(this_ioneq)
        ;print,' ionization potential = ',ip
        ;
        test1=max(this_ioneq gt 0.)
        if keyword_set(min_abund) then begin
                test2 = this_abund gt min_abund
        endif else test2 = this_abund gt 0.
        test3=ip gt 0.
        ;
        ;print,' test1,2,3 = ',test1, test2, test3
        ;
        if test1 and test2 and test3 then begin
        ;
            ;ipz=z^2*ryd
            ;ipnist=ip*ryd/(109737.3d)
            ;print,' ipz , ipnist = ',ipz,ipnist
            ;
            for it=0,n_ioneq_t-1 do begin
                t=10.^ioneq_t
                g2x=z^2*ryd/(kb*t(it))
                g2eff=ip*ryd/(109737.3d*kb*t(it))
                ;  print,' g2 = ',g2x,'  g2eff = ',g2eff
                lng2x=alog10(z^2*ryd/(kb*t(it)))
                lng2eff=alog10(g2eff)
                idx=fix((lng2eff+4.)/.2)
                idx=((idx>0)<40)
                delta=lng2eff-g2(idx)
                gtot=gff(idx)+delta*(s1(idx)+delta*(s2(idx)+delta*s3(idx)))
                ;  print,idx,s1(idx),s2(idx),s3(idx)
                ;  print,t(it),lng2x,idx,fix(idx),g2(idx),delta,gtot
                rad_loss(it)=rad_loss(it)+this_abund*this_ioneq(it)*factor*sqrt(t(it))*z^2*gtot
            endfor ;  nt
        ;
        endif  ;  various tests
        ;
    endfor;  ion
    ;
endfor;  iz
;
;for it=0,n_ioneq_t-1 do begin
;   print,t(it),rad_loss(it)
;endfor
;
end
