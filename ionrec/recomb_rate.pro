
;+
; PROJECT     :  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving George Mason University, 
;       the University of Michigan and the University of Cambridge.
;       
;
; NAME:
;     RECOMB_RATE
;
; PURPOSE:
;     This routine computes the total recombination rate
;     coefficient of an ion in units cm^3 s^-1. The rates are taken
;     from a variety of sources and are summarized in Dere et
;     al. (2009, A&A, 498, 915).
;
; INPUTS:
;    GNAME     A string, specifying the name of the recombining ion in
;              CHIANTI style.  For example 'o_6' specifies O VI or O+5.
;
;    TEMPERATURE    Temperature (units: K).
;
;
; OUTPUTS
;
;    The total (radiative + dielectronic) recombination rate
;    coefficient (units: cm^3 s^-1). If there is a problem, then an
;    array of the same size as TEMPERATURE will be returned, but
;    containing all zeros.
;
; OPTIONAL INPUTS:
;    Z, ION  Z specified the nuclear charge and ion specifies the ionization state in
;            Chianti style.  Both must be set.  For example, Z=8 and ION=6 specifies the
;            ion 'o_6'
;    Pressure:  The electron pressure (N_e*T; units cm^-3 K) at which
;               the rates are required.
;    Density:   The electron number density (units: cm^-3) at which
;               the rates are required.
;    Path:   Allows the directory containing the recombination files
;            to be directly specified. If not set, then the files are
;            obtained from the user's CHIANTI directory.
;
; KEYWORDS:
;    RADIATIVE  If set, then only the radiative recombination rate is
;               returned. 
;
;    DIELECTRONIC  If set, then only the dielectronic recombination
;                  rate is returned.
;
;    TOTAL       If set, then the routine will look for a total
;                recombination rate file (extension .trparams) and
;                read this (if it exists).
;    
; PROGRAMMING NOTES:
;    The recombination rates are stored in files with the extensions:
;       .rrparams   Radiative recombination rates  (RR)
;       .drparams   Dielectronic recombination rates  (DR)
;       .trparams   Total recombination rates  (TR)
;    These files are stored in the ion directores, along with the
;    .elvlc, .wgfa and .scups files. Generally atomic physicists
;    calculate the RR and DR rates separately. The group of Nahar &
;    Pradhan, however, compute total recombination rates
;
; MODIFICATION HISTORY:
;    Ver.1, 17-Nov-2006, Ken Dere
;    Ver.2, 1-Aug-2013,  Ken Dere
;           corrected radiative recombination calculation
;    Ver.3, 18-May-2015, Peter Young
;           added /radiative and /dielectronic keywords and updated
;           header
;    Ver.4, 25-Oct-2016, Peter Young
;           introduced /TOTAL keyword so that total rate files
;           (.trparams) are only read if this keyword is
;           set. Previously the .trparams files had the highest
;           priority, and so were read instead of .rrparams and
;           .drparams if they existed (although they only existed for
;           a few ions). Added PATH= input.  Also rates are now
;           obtained with the new routines ch_rad_recomb,
;           ch_diel_recomb and ch_tot_recomb.
;    Ver.5, 20-Feb-2017, Peter Young
;           gone back to the old method of prioritizing the .trparams
;           files over the separate .rrparams and .drparams.
;    Ver.6, 21-Aug-2018, Peter Young
;           Fixed typo in header (O VII -> VI); no change to code.
;    Ver.7, 15-Jun-2020, Peter Young
;           Added DENSITY and PRESSURE optional inputs, and added call
;           to ch_dr_suppress.pro. These allow the density suppression
;           of DR to be included.
;
; VERSION     :  7, 15-Jun-2020
;
;-
FUNCTION recomb_rate,gname,temperature,z=z,ion=ion,verbose=verbose, $
                     radiative=radiative, dielectron=dielectronic, $
                     total=total, path=path, density=density, $
                     pressure=pressure
; 
;   
;
if n_params() lt 2 then begin
   print,' '
   print,' > rate = recomb_rate(gname,temperature,[z=z,ion=ion]) '
   print,'    calculate the recombination rate coefficient in cm^3 s^-1 '
   print,'    as a function of temperature (K) '
   print,' '
   return,-1
endif
;
if not keyword_set(iz) and not keyword_set(ion) then convertname,gname,z,ion
;
str=''
;
t=temperature

;
; The following creates the path to the recombination files, with the
; ion name appended. For example /chianti/dbase/fe/fe_13/fe_13.
;
IF n_elements(path) NE 0 THEN BEGIN
  fname=concat_dir(path,gname)
ENDIF ELSE BEGIN 
  zion2filename,z,ion,fname
ENDELSE

;
; Create the names of the three files
;
rrfile=fname+'.rrparams'
drfile=fname+'.drparams'
trfile=fname+'.trparams'


quiet=1-keyword_set(verbose)

rate=fltarr(n_elements(temperature))

;
; The routines ch_rad_recomb, ch_diel_recomb and ch_tot_recomb return
; a -1 if there is a problem (e.g., file missing).
;
rr=ch_rad_recomb(gname,temperature,quiet=quiet,filename=rrfile)
IF rr[0] EQ -1 THEN rr_rate=rate ELSE rr_rate=rr

dr=ch_dr_suppress(gname,temperature,density=density,quiet=quiet,filename=drfile,pressure=pressure)
IF dr[0] EQ -1 THEN dr_rate=rate ELSE dr_rate=dr

tr=ch_tot_recomb(gname,temperature,quiet=quiet,filename=trfile)
IF tr[0] EQ -1 THEN tr_rate=rate ELSE tr_rate=tr

;
; If one of these keywords is set, then only that rate will be returned.
;
IF keyword_set(radiative) THEN return,rr_rate
IF keyword_set(dielectronic) THEN return,dr_rate
IF keyword_set(total) THEN return,tr_rate

;
; PRY, 20-Feb-2017
; The code below prioritizes the total recombination rate data over
; the separate RR+DR data. The commented out code is the situation
; where the separate RR+DR data are prioritized (v.4 of routine).
;
IF tr[0] NE -1 THEN BEGIN
  rate=tr_rate
ENDIF ELSE BEGIN
  rate=rr_rate+dr_rate
ENDELSE 

;; IF rr[0] EQ -1 AND dr[0] EQ -1 THEN BEGIN
;;   rate=tr_rate
;; ENDIF ELSE BEGIN
;;   rate=rr_rate+dr_rate
;; ENDELSE 

return,rate



END

