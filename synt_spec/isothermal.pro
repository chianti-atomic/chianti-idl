
;+
; NAME
;
;       ISOTHERMAL
;
; PURPOSE:
;       Computes spectra from isothermal plasmas. A number of isothermal 
;       plasmas can be included.
;       Note that this routine has a number of unique features that 
;       distinguish it from the other CHIANTI synthetic spectra routines. 
;       See the Programming Notes section.
;
; INPUTS:
;
;       WMIN      Minimum of desired wavelength range in Angstroms.
;
;       WMAX      Maximum of desired wavelength range in Angstroms.
;
;       WAVESTEP  Bin size of spectrum (in Angstroms)
;
;       TEMP      Electron temperature (or array).
;
; OPTIONAL INPUTS
;
;       PRESSURE  Electron pressure in units of cm^-3 K.
;
;       EDENSITY  Electron density in units of cm^-3.
;
;       EM        Emission measure. The units of EM govern the intensity 
;                 units of the emission lines (i.e., column or volume 
;                 emission measure units). If EM is not specified, then the 
;                 emission measure is set to (N_e * N_H) where N_e is 
;                 derived from the user-specified PRESSURE or EDENSITY, 
;                 and N_H is derived from the routine PROTON_DENS.PRO.
;
;       SNGL_ION  Rather than include the entire list of CHIANTI ions in 
;                 the calculation, this input can be used to select a 
;                 single ion, or a number of different ions. E.g., 
;                 SNGL_ION='s_2' or SNGL_ION=['s_2','s_3','s_4'].
;
;       RADTEMP   The blackbody radiation field temperature (default 6000 K).
;
;       RPHOT    Distance from the centre of the star in stellar radius units.
;                I.e., RPHOT=1 corresponds to the star's surface. (Default is
;                infinity, i.e., no photoexcitation.)
;
;       MASTERLIST  The list of ions that will be considered for the 
;                   spectrum is contained in the masterlist file in the 
;                   CHIANTI directories. The user can specify his own file 
;                   through this keyword. E.g., 
;                   masterlist='/user/home/masterlist.ions'
;
;
;	ABUND_NAME:  Name of the abundance file to use.  If not passed, then
;		     the user is prompted for it.
;
;	IONEQ_NAME:  Name of the ionization equilization name to use.  If not
;		     passed, then the user is prompted for it.
;
;
;  KEYWORDS
;
;       NOPROT     Switch off the inclusion of proton rates in the level 
;                  balance.
;
;       ERGS       The units of the output spectrum are by default in photons. 
;                  Setting /ERGS switches to erg units.
;
;       CONT       Adds continuum (free-free, free-bound, two-photon) to 
;                  spectrum.
;  
;       ALL        Include all lines, i.e. also those for which wavelengths 
;                  are only theoretical and not observed. 
;
;  OUTPUTS:
;
;        LAMBDA   Wavelength array of calculated synthetic spectrum.
;
;        SPECTRUM Intensity array. The units depend on the user inputs to 
;                 ISOTHERMAL -- see note below. 
;
;        LIST_WVL A list of wavelengths for use with synthetic_plot.pro
;
;        LIST_IDENT A list of line identifications for use with 
;                   synthetic_plot.pro
;
; PROGRAMMING NOTES
;
;        Intensity Units
;        ---------------
;        The units are of the form photons cm^3 s^-1 sr^-1 * (units of EM), 
;        changing to ergs if the /ergs keyword is set.
;
;        The volume emission measure (N_e*N_H*V) has units cm^-3.
;
;        The column emission measure (N_e*N_H*h) has units cm^-5.
;
;
;        Unique features
;        ---------------
;        The emission lines in the final spectrum have no width and thus 
;        each occupies a single pixel of the spectrum. The size of the 
;        pixels are set by WAVESTEP.
;
;        As stated above, the units of the output spectrum are 
;        photons cm^3 s^-1 sr^-1, i.e., there is no "per angstrom" term. 
;        This means that (i) the height of the emission lines in the 
;        spectrum does not change with varying WAVESTEP, and (ii) the height
;        of continuum does change with WAVESTEP.
;
; COMMON BLOCKS
;
;        ELEMENTS
;
; CALLS
;
;        CH_SYNTHETIC, READ_ABUND, CH_GET_FILE, CONCAT_DIR, FREEFREE, 
;        FREEBOUND, TWO_PHOTON
;
; HISTORY
;        Ver.1, 8-Apr-02, Peter Young  Rutherford Appleton Laboratory,
;        p.r.young@rl.ac.uk 
;        Tries to replicate the behaviour of the original ISOTHERMAL which 
;        was found in earlier versions of CHIANTI (v.3 and earlier). 
;
; MODIFICATION HISTORY
;
;       Ver. 2, Giulio Del Zanna (GDZ), 28-Apr-02 
;               Added abund_name,ioneq_name keywords.
;               Also, added photons keyword in call to MAKE_CHIANTI_SPEC.
;
;       Ver. 3, Peter Young, 24-May-02
;                 Modified to produce arrays of spectra when an array of 
;                 temperatures is given
;
;       V.4, GDZ, 28-May-02 
;              Added a couple of checks on file existence and modified the call
;              to ch_synthetic and make_chianti_spec  due to change of keyword
;              names.  
;
;       V.5, Peter Young, 16-Jul-02
;              Restructured routine to avoid crashes when a large number of 
;              temperatures is input.
;
;       V.6, 8-Aug-02 GDZ
;              Added one error checking
;
;       V.7, 18-Aug-03, Peter Young
;              Added EM= keyword.
; 
;       V.8, 14-Sept-2005 GDZ 
;              Added ALL keyword and modified header, error message.
;
;       V.9, 3-Oct-2005, GDZ
;              Now the FOR loop accepts more than 32000 lines.
;
;       V.10, 7-Oct-2009, Peter Young
;              The keyword EM_INT is now used in the call to the
;              continuum routines.
;
;       V.11, 17-Aug-2015, Peter Young
;              Changed 'lambda' to 'lmbda' as IDL have introduced a
;              function called lambda.
;
;       V.12, 5-Dec-2018,  GDZ
;              replaced noverbose with verbose.
;
;       V.13, 20-Feb-2019, Peter Young
;              Added check on TEMP to make sure it's in
;              ascending temperature order; now uses !ioneq_file for
;              the ion balance file if ioneq_name is not set.
;
; VERSION     : 
;       Version 13, 20-Feb-2019
;        
;-


PRO isothermal, wmin,wmax,wavestep,temp,lmbda,spectrum,list_wvl,list_ident,$
      pressure=pressure,edensity=edensity,photons=photons, ergs=ergs, $
      sngl_ion=sngl_ion, abund_name= abund_name , ioneq_name=ioneq_name, $
      verbose=verbose, min_abund=min_abund, cont=cont, $
      masterlist=masterlist, noprot=noprot, radtemp=radtemp, $
      rphot=rphot, em=em,all=all

COMMON elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref

IF n_params() LT 6 THEN BEGIN
  print,'Use: IDL> isothermal,wmin,wmax,wavestep,temp,lambda,spectrum,list_wvl,'
  print,'                 list_ident, pressure= , edensity= , /photons,'
  print,'                 /verbose, sngl_ion= , min_abund= , /cont,'
  print,'                 masterlist= ,abund_name= , ioneq_name= ,/noprot]'
  print,'                 radtemp= ,rphot= ,em= , /all '
  print,''
  print,'E.g., IDL> isothermal,170,180,0.1,1e6,lambda,spectrum,list_wvl,'+$
       'list_ident,'
  print,'                  edensity=1e10'
  return
ENDIF


;
; PRY, 24-Jul-2018
;  Check to make sure temp is in ascending order (if it's an array)
;
chck_temp=temp[sort(temp)]
k=where(chck_temp NE temp,nk)
IF nk GT 0 THEN BEGIN
  print,'% ISOTHERMAL: the TEMP array must be in ascending temperature order. Returning...'
  return
ENDIF 


IF keyword_set(ergs) THEN photons=0 ELSE photons=1

IF n_elements(spectrum) NE 0 THEN junk=temporary(spectrum)

IF n_elements(min_abund) EQ 0 THEN min_abund=0.

IF (n_elements(pressure) NE 0 AND n_elements(edensity) NE 0) OR $
   (n_elements(pressure) EQ 0 AND n_elements(edensity) EQ 0) THEN BEGIN
  print,'%ISOTHERMAL: Please specify either EDENSITY or PRESSURE'
  return
ENDIF


n=n_elements(edensity)
nt=n_elements(temp)

IF n_elements(em) NE 0 THEN BEGIN
  IF n_elements(em) EQ nt THEN BEGIN
    logem_isothermal=alog10(em)
  ENDIF ELSE BEGIN
    logem_isothermal=alog10(em[0])
  ENDELSE
ENDIF ELSE BEGIN
  logem_isothermal=dblarr(nt)
ENDELSE


;
; PRY, 24-Jul-2018
; Get rid of widget selection for ioneq file
;
IF n_elements(ioneq_name) EQ 0 THEN ioneq_name=!ioneq_file
;
;; IF n_elements(ioneq_name) EQ 0 THEN BEGIN 
;;   dir=concat_dir(!xuvtop,'ioneq')
;;   ioneq_name=ch_get_file(path=dir,filter='*.ioneq', $
;;                          title='Select Ionization Equilibrium File')
;; IF NOT file_exist(ioneq_name) THEN message, 'Error, no file selected -- EXIT'
;; END


IF n_elements(abund_name) EQ 0 THEN BEGIN 
  dir=concat_dir(!xuvtop,'abundance')
  abund_name=ch_get_file(path=dir,filter='*.abund', $
                         title='Select Abundance File')
IF NOT file_exist(abund_name) THEN message, 'Error, no file selected -- EXIT'
END

IF nt NE 1 THEN no_sum_int=1

spectrum=-1

ch_synthetic,wmin,wmax,out=out, press=pressure, err_msg=err_msg, $
     sngl_ion=sngl_ion,ioneq_name=ioneq_name,dem_name=dem_name, $
     photons=photons, masterlist=masterlist, noprot=noprot, $
     radtemp=radtemp, rphot=rphot, verbose=verbose, progress=progress, $
     density=edensity, no_sum_int=no_sum_int, logt_isothermal=alog10(temp), $
     logem_isothermal=logem_isothermal,all=all

IF err_msg [0]  NE  '' THEN BEGIN 
  print, '% ISOTHERMAL: Error - EXIT'
  return
END 

nl=n_elements(out.lines)

format='('+strtrim(string(nt),2)+'e10.2)'

read_abund,abund_name,abund,ref

nw=fix((wmax-wmin)/wavestep)+1
lmbda=findgen(nw)*wavestep+wmin
spectrum=dblarr(nw,nt)

str1=out.lines.snote+' '+out.lines.ident
FOR i=0L, long(nl-1) DO BEGIN
  int=out.lines[i].int*abund[out.lines[i].iz-1]
  str1[i]=strpad(str1[i],50,/after)+'  Int='+string(format=format,int)
  getmin=min(abs(lmbda-out.lines[i].wvl),ind)
  spectrum[ind,*]=spectrum[ind,*]+int
ENDFOR

list_ident=str1
list_wvl=out.lines.wvl

i=sort(list_wvl)
list_wvl=list_wvl[i]
list_ident=list_ident[i]

IF keyword_set(cont) THEN BEGIN
 ;
  em_int=double(10.)^logem_isothermal
 ;
  IF n_elements(edensity) EQ 0 THEN edensity=pressure/temp
  freebound, temp,lmbda,fb,/no_setup,min_abund=min_abund, $
       photons=photons, em_int=em_int
  freefree,temp,lmbda,ff,/no_setup,min_abund=min_abund, $
       photons=photons, em_int=em_int
  two_photon, temp,  lmbda, two_phot,/no_setup,min_abund=min_abund, $
       edensity=edensity, photons=photons, em_int=em_int
  totcont=(fb+ff+two_phot)/1d40*wavestep
  spectrum=spectrum+totcont
ENDIF

print, '% ISOTHERMAL: calculation finished'

END
