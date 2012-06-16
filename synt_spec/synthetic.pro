;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory (USA), the University of Florence (Italy), the
;       University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
;
; NAME:
;	SYNTHETIC
;
; PURPOSE:
;
;       calculates a synthetic spectrum
;
;
; PROCEDURE:
;
;                 Calculations are done assuming either constant density or
;                 constant pressure. See CH_SYNTHETIC for details.
;
;
; CALLING SEQUENCE:
;
;       SYNTHETIC,Wmin, Wmax, Fwhm, Pressure= , Lambda, Spectrum ,List_wvl, List_ident
;                 ,[/all, density=, /cont, min_abund=]
;
;
; INPUTS:
;
;
;	Wmin:   lower limit of the wavelength/energy range of interest (Angstroms)
;
;	Wmax:   upper limit of the wavelength/energy range of interest (Angstroms)
;
;       Pressure:  pressure in emitting region (cm^-3 K), or 
;       Density:   density in emitting region (cm^-3).
;
;       Fwhm:  gaussian full width at half maximum of the resolution of the output 
;                  spectrum, for example, to correspond to an observed spectrum
;
;
; OPTIONAL INPUTS:
;
;	SNGL_ION:  specifies  a single ion (e.g. SNGL_ION='Fe_10' to include
;                 only Fe X lines) or an array (e.g. SNGL_ION=['Fe_10','Fe_11']
;                 to include only Fe X and Fe XI lines) of ions to be used
;                 instead of the complete set of ions specified in
;                 !xuvtop/masterlist/masterlist.ions 
;
;       MASTERLIST: string of a specific masterlist file (full path). 
;                   If defined as a keyword (i.e. MASTERLIST=1 or /MASTERLIST)
;                   then a widget allows the user to select a  user-defined
;                   masterlist file. Shortcut for SNGL_ION.   
;
;
;       MIN_ABUND:  If set, calculates the continuum only from those elements which 
;                   have an abundance greater than min_abund. 
;
;	DEM_NAME:  Name of the DEM file to used.  If not passed, then the user
;		   is prompted for it.
;
;	ABUND_NAME:  Name of the abundance file to use.  If not passed, then
;		     the user is prompted for it.
;
;	IONEQ_NAME:  Name of the ionization equilization name to use.  If not
;		     passed, then the user is prompted for it.
;
;       RPHOT   Distance from the centre of the star in stellar radius units.
;               I.e., RPHOT=1 corresponds to the star's surface. (Default is
;               infinity, i.e., no photoexcitation.)
;
;       RADTEMP The blackbody radiation field temperature (default 6000 K).
; 
;
; OUTPUTS:
;
;       Lambda:  wavelength array of calculated synthetic spectrum
;       Spectrum:  intensity array (erg cm^-2 s^-1 str^-1 Ang^-1),
;                  unless keyword photons is set then output is is
;                  photons cm^-2 s^-1 str^-1 Ang^-1
;       List_wvl:  a list of wavelengths for use with synthetic_plot.pro
;       List_ident:  a list of line identifications for use with 
;                        synthetic_plot.pro
;
; OPTIONAL OUTPUTS:
;
;	
; KEYWORD PARAMETERS:
;
;     
;	ALL:  if set, then all lines are included.  This means that lines for which
;             only an approximate wavelength is known (only theoretical energy
;             levels are known) are included.
;
;	SNGL_ION:  specifies  a single ion (e.g. SNGL_ION='Fe_10' to include
;                 only Fe X lines) or an array (e.g. SNGL_ION=['Fe_10','Fe_11']
;                 to include only Fe X and Fe XI lines) of ions to be used
;                 instead of the complete set of ions specified in
;                 !xuvtop/masterlist/masterlist.ions 
;
;       MASTERLIST: string of a specific masterlist file (full path). 
;                   If defined as a keyword (i.e. MASTERLIST=1 or /MASTERLIST)
;                   then a widget allows the user to select a  user-defined
;                   masterlist file. Shortcut for SNGL_ION.   
;
;       CONTINUUM:   if set, then the continuum (free-free, free-bound and
;                  two-photon) are  included 
;
;       MIN_ABUND:  If set, calculates the continuum only from those elements which 
;                   have an abundance greater than min_abund.  Can speed up the 
;                   calculations.  For example, from Allen (1973):
;                   abundance (H)  = 1.
;                   abundance (He) = 0.085
;                   abundance (C)  = 3.3e-4
;                   abundance (O)  = 6.6e-4
;                   abundance (Si) = 3.3e-5
;                   abundance (Fe) = 3.9e-5
;
;
;
;	PHOTONS:  if set, intensities are in photons cm^-2 s^-1 sr^-1 Ang^-1
;
;	DEM_NAME:  Name of the DEM file to used.  If not passed, then the user
;		   is prompted for it.
;
;	ABUND_NAME:  Name of the abundance file to use.  If not passed, then
;		     the user is prompted for it.
;
;	IONEQ_NAME:  Name of the ionization equilization name to use.  If not
;		     passed, then the user is prompted for it.
;
;       NOPROT   If set, then proton rates are not included.
;
;	RADTEMP	Specify background radiation temperature (default: 6000 K)
;
;	RPHOT   Distance from the centre of the star in stellar radius units.
;               I.e., RPHOT=1 corresponds to the star's surface. (Default is
;               infinity, i.e., no photoexcitation.)
;
; CALLS:
;       
;       CH_SYNTHETIC, MAKE_CHIANTI_SPEC, READ_ABUND, STRPAD
;
;
; COMMON BLOCKS: None
;
;
; RESTRICTIONS:
;
; SIDE EFFECTS:
;
;
; EXAMPLE:
;
;       > synthetic,100.,200.,.1, pressure=1.e+16,lambda,spectrum,list_wvl,list_ident
;
;
; CATEGORY:
;
;	spectral synthesis.
;
; WRITTEN     : 
;
;       Version 1, 8-Nov-01, Giulio Del Zanna (GDZ). 
;
;       Rewritten as a wrapper routine using the new procedures.
;
;       Compared to the previous SYNTHETIC, these are the main changes:
;
;       1-Now the PRESSURE value is a keyword as the DENSITY value
;       2-The keyword CONT is now renamed CONTINUUM
;       3-Added keywords PHOTONS, DEM_NAME, ABUND_NAME, IONEQ_NAME
;       4-MASTERLIST can now be used both as an input string or as a keyword.
;       5-The description of the line details now has the spectroscopic 
;         designation at the end.
;
;
; MODIFICATION HISTORY:
;
;       Version 2, 18-Nov-01, Peter Young
;           Added /noprot, rphot and radtemp keywords.
;
;       Version 3, 11-Dec-01, Peter Young
;           Changed call to ch_strpad to strpad.
;
;       Version 4, 28-Apr-02, GDZ, changed the call to make_chianti_spec and the
;       continuum keyword.
;
;       V. 5, 22-May-2002 GDZ.  Removed const_net definitions.
;
;       V.6, 14-Feb-2003 GDZ.
;             Fixed a bug (keyword PHOTONS was not active). 
;             
; VERSION     : 6, 14-Feb-2003
;
;
;-

PRO  synthetic, wmin, wmax, fwhm, pressure=pressure, density=density,$
                lambda, spectrum, list_wvl, list_ident, $
                all=all, sngl_ion=sngl_ion, min_abund=min_abund, $
                continuum=continuum, $
                elapsed=elapsed, masterlist=masterlist,$
                dem_name=dem_name,abund_name=abund_name, $
                ioneq_name=ioneq_name, photons=photons, $
                noprot=noprot, rphot=rphot, radtemp=radtemp


if n_params() lt 7 then begin
   print,' type> synthetic,wmin,wmax,fwhm,pressure= ,lambda,spectrum,list_wvl,list_ident, '
   print,'         [density= , /all, sngl_ion= , /cont, min_abund=3.e-5, masterlist= , '
   print,'          radtemp= , rphot= , /noprot ]'
   print,'   '
   print,'     continuum calculations (/cont set) are fairly slow, setting higher min_abund values'
   print,'        can help but these also affect the spectral line calculations'
   print,'     Allen abundances for some elements:'
   print,'     abund(H)  = 1.'
   print,'     abund(He) = 8.5e-2'
   print,'     abund(C)  = 3.3e-4'
   print,'     abund(N)  = 9.1e-5'
   print,'     abund(Mg) = 2.6e-5'
   print,'     abund(Na) = 1.8e-6'
   print,'     abund(Si) = 3.3e-5'
   print,'     abund(S)  = 1.5e-5'
   print,'     abund(Ca) = 2.0e-6'
   print,'     abund(Fe) = 4.0e-5'
   return
endif

IF  keyword_set(pressure)  AND  keyword_set(density)  THEN  begin
   print,' You  have  to decide if you want the intensities calculated '+$
     'for constant pressure or not ! '
   print,'(e.g.  pressure=1.e16) or '+$
     'for constant electron density (e.g.  density=1.e9) !!!'
   err_msg ='Error '
   return
endif

IF n_elements(density) NE 0 THEN BEGIN 
   print,' using constant density = ',density
ENDIF ELSE IF  n_elements(pressure) NE 0 THEN BEGIN 
   print,' using constant Pressure = ', pressure 
ENDIF ELSE message, 'You must specify either the Electron Density (DENSITY) or'+$
  ' the Pressure (PRESSURE = cm^-3 K)'


;
; check on how long the procedure takes
;
IF  keyword_set(elapsed) THEN  time0=systime(1)  
;


ch_synthetic, wmin, wmax, output=TRANSITIONS, err_msg=err_msg, $
  pressure=pressure, density=density,all=all,sngl_ion=sngl_ion, $
  photons=photons,  masterlist=masterlist,  $
  verbose=verbose,$
      logt_isothermal=logt_isothermal,  logem_isothermal=logem_isothermal,$
  ioneq_name=ioneq_name, dem_name=dem_name, $
  noprot=noprot, rphot=rphot, radtemp=radtemp


IF err_msg NE '' THEN BEGIN 
   print, err_msg
   return 
ENDIF 

pixel = 0
KEV = 0
INSTR_FWHM=fwhm
WRANGE = [wmin, wmax]

;set binsize = 1/10 the FHWM
BIN_SIZE = FWHM/10.


MAKE_CHIANTI_SPEC, TRANSITIONS, X, Y, BIN_SIZE=BIN_SIZE, $
  INSTR_FWHM=INSTR_FWHM,  $
  WRANGE=WRANGE, ALL=ALL, continuum=continuum, photons=photons, $
  ABUND_NAME=ABUND_NAME, MIN_ABUND=MIN_ABUND,  BINSIZE=BINSIZE

lambda = X
spectrum = Y.spectrum

;line wavelengths (A):

list_wvl = TRANSITIONS.lines[*].wvl

;read abundance file

read_abund,abund_name,abund,abund_ref

intensity = TRANSITIONS.lines[*].int * abund[TRANSITIONS.lines[*].iz-1]


;line ID in ascii form:

; IF keyword_set(ALL) THEN 

list_ident = strpad(TRANSITIONS.lines[*].snote, 14,/after)+$
  '  Int='+string(intensity[*], '(e10.2)')+$
  '  Tmax='+string(TRANSITIONS.lines[*].tmax,'(f4.1)')+$
  ' '+ strtrim(string(TRANSITIONS.lines[*].ident), 2)


;sort them 

isort = sort(list_wvl) 

list_wvl = list_wvl(isort)
list_ident =list_ident(isort)


IF  keyword_set(elapsed) THEN  BEGIN 
   time1=systime(1)
   dt=time1-time0
   print,' elapsed seconds in synthetic = ',dt
ENDIF 


END 


