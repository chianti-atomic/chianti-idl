
FUNCTION ch_continuum, temp, wvl, photons=photons, density=density, $
                       pressure=pressure, $
                       abund_file=abund_file, em_int=em_int, dem_int=dem_int, $
                       sumt=sumt, ioneq_file=ioneq_file, $
                       sngl_ion=sngl_ion, element=element, $
                       advanced_model=advanced_model, min_abund=min_abund, $
                       kev=kev

;+
; NAME:
;     CH_CONTINUUM
;
; PURPOSE:
;     Derives the total continuum spectrum from CHIANTI. The routine is
;     a wrapper to the freefree, freebound and two_photon routines.
;
; CATEGORY:
;     CHIANTI; continuum.
;
; CALLING SEQUENCE:
;     Result = CH_CONTINUUM( Temp, Wvl )
;
; INPUTS:
;     Temp:   Temperature in K (can be an array). If the input
;             DEM_INT is used then the temperatures must be spaced at
;             equal intervals of 
;             log10(T). E.g., log T=6.0, 6.1, 6.2, 6.3, ...
;     Wvl:    Wavelength in angstroms (can be an array). If /keV is
;             set, then WVL is assumed to contain energies in keV
;             units. In both cases WVL should be monotonically
;             increasing. 
;
; OPTIONAL INPUTS:
;     Dem_Int:  This should be the same size as TEMP and contain
;               differential emission measure values. Specifying
;               DEM_INT places a restriction on TEMP (see above). 
;     Em_Int:   This should be the same size as TEMP and contain
;               emission measure values. The emission measure is
;               N_e*N_H*h, where N_e is electron number density, N_H
;               is the hydrogen number density and h the plasma column
;               depth. Units are cm^-5.
;     Sngl_Ion: A string or string array specifying individual ions
;               for which you want the continuum. Ions should be
;               specified in the CHIANTI format, for example, 'fe_13'
;               for Fe XIII. The name of the *recombined* ion is
;               specified. For example, if you want the hydrogen
;               continuum, use 'h_1' (not 'h_2').
;     Element:  An array specifying which elements are used for the
;               calculation. An element can be specified either
;               with an integer (e.g., 26 for iron) or a string
;               (e.g., 'fe' for iron).
;     Abund_File: The name of a CHIANTI format element abundance
;               file. If not specified then a widget will appear
;               asking you to choose a file. The CHIANTI default file
;               can be specified by giving !abund_file.
;     Ioneq_File: The name of a CHIANTI format ionization equilibrium
;               file. If not specified then the CHIANTI default file
;               (!ioneq_file) will be used.
;     Density:  Scalar giving density (cm^-3) for the call to
;               two_photon. The default is 10^10.
;     Pressure: Scalar giving pressure (K cm^-3). The value is used
;               to calculate the density using the input temperature
;               array, and the density is input to two_photon.
;     Min_Abund: Sets a minimum abundance value such that only
;                elements with abundances higher than this value will
;                be included in the calculation. For example
;                min_abund=1e-6 will exclude the minor elements.
;	
; KEYWORD PARAMETERS:
;     PHOTONS: The output spectrum is given in photon units rather 
;              than ergs.
;     SUMT:    When a set of temperatures is given to FREEBOUND, the 
;              default is to output INTENSITY as an array of size 
;              (nwvl x nT). With this keyword set, a summation over 
;              the temperatures is performed.
;     ADVANCED_MODEL: By default, this is set to 1 and thus a new
;              ionization balance is calculated using the CHIANTI
;              advanced models. To switch it off, set
;              advanced_model=0 and either the default ion balance
;              (!ioneq_file) will be used, or the one specified by
;              ioneq_file=.
;     KEV:     If set, then the units of WVL are assumed to be keV.
;
; OUTPUTS:
;     An IDL structure with the following tags:
;      .wvl   Input wavelength (WVL).
;      .int   Continuum intensity.
;      .temp  Input temperature (TEMP).
;      .dem_int  Input DEM values (set to -1 if not set).
;      .em_int  Input EM values (set to -1 if not set).
;      .sumt  Indicates if /sumt was set.
;      .pressure Value of input pressure (set to -1 if not set).
;      .density Value of input density (set to 10^10 if not set and
;               pressure is not set).
;      .kev   Set to 1 if /kev was set.
;      .time_stamp  String giving time at which spectrum was created.
;      .int_fb  The freebound component of the intensity.
;      .int_ff  The freefree component of the intensity.
;      .int_pp  The two-photon component of the intensity.
;      .chianti_version  Version of CHIANTI used to generate spectrum.
;
;     The intensity is given in units erg/cm2/s/sr/Ang. If /photons
;     was set, then the units are photon/cm2/s/sr/Ang. If /keV was set,
;     then the units will be /keV rather than /Ang.
;
;     The output intensity will have dimensions (NWVL, NTEMP). If
;     /sumt is set, then the intensity array is summed over
;     temperature to give an array of size NWVL.
;
; EXAMPLE:
;     IDL> wvl=findgen(100)+1.
;     IDL> temp=10.^(findgen(11)/10.+6.5)
;     IDL> c=ch_continuum(wvl,temp)
;
;     IDL> c=ch_continuum(wvl,temp,advanced=0,density=1e8)
;
; MODIFICATION HISTORY:
;     Ver.1, 25-Nov-2025, Peter Young
;-

IF n_params() LT 2 THEN BEGIN
  print,'Use:  IDL> output=ch_continuum( temp, wvl [, /photons, pressure=, abund_file=, em_int=, '
  print,'                                dem_int=, /sumt, ioneq_file=, sngl_ion=, element=, '
  print,'                                advanced_model= , density=] )'
  return,-1
ENDIF 

IF n_elements(pressure) NE 0 AND n_elements(density) NE 0 THEN BEGIN
  message,/info,/cont,'Both PRESSURE and DENSITY have been specified. Please use only one. Returning...'
  return,-1.
ENDIF

IF n_elements(pressure) GT 1 or n_elements(density) GT 1 THEN BEGIN
  message,/info,/cont,'PRESSURE or DENSITY must be scalars. Returning...'
  return,-1.
ENDIF

;
; Need to define density for call to two_photon.
;
density_scalar=-1.
pressure_scalar=-1.
;
IF n_elements(pressure) NE 0 THEN BEGIN
  density=pressure/temp
  pressure_scalar=pressure
  pressure_ioneq=pressure_scalar
ENDIF ELSE BEGIN
  ;
  IF n_elements(density) EQ 0 THEN density=1e10
  density_scalar=density
  density_ioneq=density_scalar
  ;
ENDELSE 

IF n_elements(abund_file) EQ 0 THEN abund_file=!abund_file

IF n_elements(advanced_model) EQ 0 THEN advanced_model=1b
;
IF NOT keyword_set(advanced_model) THEN BEGIN
  IF n_elements(ioneq_file) EQ 0 THEN ioneq_file=!ioneq_file
  swtch=0b
ENDIF ELSE BEGIN
  ioneq_file=concat_dir(getenv('IDL_TMPDIR'),'chianti.ioneq')
  log_temp=findgen(101)/20.+4.0
  ioneq_data=ch_calc_ioneq(10.^log_temp, dens=density_ioneq, press=pressure_ioneq, /adv, $
                           outname=ioneq_file)
  junk=temporary(ioneq_data)
  swtch=1b
ENDELSE 

freebound,temp,wvl,int_fb,abund_file=abund_file,photons=photons, $
          em_int=em_int,dem_int=dem_int,ioneq_file=ioneq_file, $
          element=element, sngl_ion=sngl_ion,min_abund=min_abund, $
          kev=kev
freefree,temp,wvl,int_ff,abund_file=abund_file,photons=photons, $
          em_int=em_int,dem_int=dem_int,ioneq_file=ioneq_file, $
          element=element, sngl_ion=sngl_ion,min_abund=min_abund, $
          kev=kev
two_photon,temp,wvl,int_pp,edensity=density,abund_file=abund_file,photons=photons, $
          em_int=em_int,dem_int=dem_int,/lookup,ioneq_file=ioneq_file, $
          element=element, sngl_ion=sngl_ion,min_abund=min_abund, $
          kev=kev

;
; Get current time and convert to UTC for use as time_stamp.
;
jd=systime(/julian,/utc)
mjd=jd-2400000.5d
mjd_str={ mjd: floor(mjd), time: (mjd-floor(mjd))*8.64d7 }
time_stamp=anytim2utc(/ccsds,mjd_str)

int=(int_fb+int_ff+int_pp)*1d-40
IF keyword_set(sumt) THEN int=total(int,2)

IF n_elements(dem_int) EQ 0 THEN dem_int=-1.
IF n_elements(em_int) EQ 0 THEN em_int=-1.

;
; This forces the output temperature to be a scalar if there is only
; one temperature.
;
IF n_elements(temp) EQ 1 THEN temp_output=temp[0] ELSE temp_output=temp

output={wvl: wvl, $
        int: int, $
        temp: temp_output, $
        dem_int: dem_int, $
        em_int: em_int, $
        sumt: keyword_set(sumt), $
        density: density_scalar, $
        pressure: pressure_scalar, $
        kev: keyword_set(kev), $
        time_stamp: time_stamp, $
        abund_file: abund_file, $
        ioneq_file: ioneq_file, $
        int_fb: int_fb*1d-40, $
        int_ff: int_ff*1d-40, $
        int_pp: int_pp*1d-40, $
        chianti_version: ch_get_version() }

;
; Delete the temporary ioneq file.
;
IF swtch THEN file_delete,ioneq_file

return,output

END
