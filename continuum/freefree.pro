
PRO freefree, temp, wvl, int, sumt=sumt, dem_int=dem_int, $
              photons=photons, min_abund=min_abund, verbose=verbose, $
              kev=kev, em_int=em_int, sngl_ion=sngl_ion, element=element, $
              abund_file=abund_file, ioneq_file=ioneq_file


;+
; NAME:
;     FREEFREE
;
; PURPOSE:
;     Computes the free-free (bremsstrahlung) continuum spectrum for
;     an ionized plasma.
;
; CATEGORY:
;     CHIANTI; continuum.
;
; CALLING SEQUENCE:
;     FREEFREE, Temp, Wvl, Int
;
; INPUTS:
;     Temp:   Temperature (in K). Can be an array. If the input
;             DEM_INT= is used then the temperatures must be spaced at
;             equal intervals of 
;             log10(T). E.g., log T=6.0, 6.1, 6.2, 6.3, ...
;
;     Wvl:    Wavelengths in angstroms. Can be a scalar or vector. If /keV is
;             set, then WVL is assumed to contain energies in keV units.
;
; OPTIONAL INPUTS:
;     Dem_Int:  This should be the same size as TEMP and contain
;               differential emission measure values. Specifying
;               DEM_INT places a restriction on TEMP (see above). Note
;               that DEM_INT values are multiplied by the factor dT
;               within the routine, whereas EM_INT values are not.
;
;     Em_Int:   This should be the same size as TEMP and contain
;               emission measure values. Note
;               that DEM_INT values are multiplied by the factor dT
;               within the routine, whereas EM_INT values are not. If
;               neither EM_INT or DEM_INT are specified, then a EM_INT
;               vector of same size as TEMP with values of unity is
;               assumed.
;
;     Min_Abund:This keyword allows the specification of a minimum abundance, 
;               such that any elements with an abundance (relative to 
;               hydrogen) less than MIN_ABUND will not be included in the 
;               calculation. E.g., MIN_ABUND=1e-5.
;
;     Abund_File: The name of an element abundance file to use. If not
;                 specified, then a widget will appear allowing you to
;                 choose one of CHIANTI abundance files.
;
;     Ioneq_File: The name of an ionization equilibrium file to use. If
;                 not specified, then the default CHIANTI ionization
;                 file (!ioneq_file) is used.
;
;     Sngl_Ion: A string containing the name of an ion in CHIANTI
;               format (e.g., 'fe_21' for Fe XXI). Can be an array. The
;               continuum will only be calculated for the specified
;               ions. 
;
;     Element:  Either the atomic number of an element, or the
;               symbol name (e.g., 'Fe' for iron). If set, then the
;               continuum will only be computed for the specified
;               element. This can not be an array.
;	
; KEYWORD PARAMETERS:
;     SUMT:    The default is to output the intensity array as an array 
;              of size (nwvl x nT). Setting this keyword performs a sum 
;              over the temperatures to yield a vector of same size as 
;              the input wavelengths, thus producing the complete 
;              free-free spectrum.
;
;     PHOTONS: Gives output emissivity in photon units rather than ergs.
;
;     KEV:     If set, then WVL is assumed to contain energies in keV units
;              rather than wavelengths in Angstroms.
;
;     VERBOSE: If set, then information will be printed to the screen.
;
; OUTPUTS:
;     Int:    Free-free continuum intensity in units of 
;             10^-40 erg cm^3/s/sr/Angstrom  
;             [ integral(N_H N_e dh) in cm^-5 ] if a DEM is not defined. 
;
;             If DEM values are defined, it is assumed that they are given
;             as N_H N_e dh/dT.  The units are 10^-40 erg/cm^2/s/sr/Angstrom. 
;
;             If T is given as a 1-D array, then the output will be a 2-D array,
;             with one element for each temperature and wavelength 
;             (but also see SUMT).
;
;             If the keyword /keV is set, then the units of INT will be 
;             10^-40 erg cm^3/s/sr/keV
;
; CALLS:
;     SUTHERLAND, ITOH
;
; PROGRAMMING NOTES:
;     This routine calls out to sutherland.pro and itoh.pro to
;     actually calculate the continuum. 
;
;     The Itoh fitting formula is only valid for (6.0 LE logT LE 8.5). 
;     For temperatures below this, we thus switch to the Sutherland 
;     fitting formula. There is very little (<1%) between the two at 
;     logT=6.
;
;     Itoh also has a constraint on the quantity u=hc/kTl (l=wavelength), 
;     such that (-4 LE log u LE 1.0). The upper limit corresponds to the 
;     continuum being cut-off prematurely at low wavelengths. E.g., for 
;     T=10^6 the cutoff is at 14.39 angstroms. For these low wavelengths 
;     we also use the Sutherland data to complete the continuum. Note that 
;     the continuum at these wavelengths is very weak
;
;     To calculate the continuum spectrum the user should specify
;     whether the plasma structure is specified through a differential
;     emission measure (DEM) curve or through emission measure (EM)
;     values. For the former, the DEM values are specified through the
;     DEM_INT keyword and for the latter the EM values are specified
;     through the EM_INT keyword. If neither are specified EM_INT
;     values of unity are assumed.
;
;     If DEM_INT values are specifed they must be specified on a
;     temperature scale that has constant intervals in log10(T). EM_INT
;     values do not have this restriction.
;
;     Note that EM_INT must be used for isothermal plasmas.
;
; EXAMPLE:
;     IDL> abfile=ch_choose_abund()
;     IDL> wvl=findgen(100)+1.
;     IDL> freefree,1e7,wvl,int,abund_file=abfile
;     IDL> freefree,1e7,wvl,int,abund_file=abfile,element='fe'
;
; MODIFICATION HISTORY:
;    Ver.1, 5-Dec-2001, Peter Young
;         Completely revised to call the separate itoh.pro and 
;         sutherland.pro routines.
;
;    V. 2, 21-May-2002,  Giulio Del Zanna (GDZ),
;          Corrected the description of the  units.
;          Added verbose keyword and a printout.
;
;    V. 3, 22-May-2002,  Peter Young (PRY)
;          Added MIN_ABUND optional input.
;          Changed ioneq_t to ioneq_logt (GDZ).
;
;    V 4, 25-May-2005, GDZ 
;                  corrected routine header.
;
;    V. 5, 9-Mar-2006, Peter Young
;          Added /kev keyword
;
;    V. 6, 9-May-2008, Peter Young & Kevin Tritz
;          Modifed code so that temperatures can be specified in a unordered
;          list.
;
;    V. 7, 7-Oct-2009, Peter Young
;          Added EM_INT= keyword; simplified how the routine decides
;          which calculation to use (Itoh or Sutherland)
;
;    v. 8, 8-Jun-2010, Ken Dere
;           added print statement following check for params
;
;    v. 9, 5-Sep-2012, Peter Young
;           updated header (no change to code).
;
;    v.10, 26-Apr-2019, Peter Young
;          Removed common block; added ABUND_FILE, IONEQ_FILE, ELEMENT
;          and SNGL_ION inputs.
;
;    v.11, 19-Aug-2019, Peter Young
;          Updated header; deleted some old, commented-out code;
;          removed references to no_setup.
;-


IF n_params() lt 3 then begin
  print,'Use: IDL> freefree, temp, wvl, int [, /sumt, dem_int=, em_int= '
  print,'                    /photons, min_abund=, /verbose, /kev, '
  print,'                    sngl_ion=, element=, ioneq_file='
  print,'                    abund_file= ] '
  print,'  '
  return
ENDIF 


t1=systime(1)

temp=double(temp)
wvl=double(wvl)

;
; * keV units *
; For the calculation of the continuum, it's simply necessary to pass the
; /kev keyword on to itoh.pro and sutherland.pro. However, there's also a
; check in freefree.pro based on wavelength that needs to be corrected. I
; thus define 'chckwvl' for this, with a different definition in the case
; of /kev.
; 
chckwvl=wvl
IF keyword_set(kev) THEN chckwvl=12.398/wvl


IF n_elements(dem_int) NE 0 THEN BEGIN
   IF n_elements(dem_int) NE n_elements(temp) THEN BEGIN
      print,'% FREEFREE: Warning, number of elements of DEM_INT must match'
      print,'  the number of temperatures. Freefree continuum not calculated.'
      return
   ENDIF
ENDIF


;
; I've modified the code below to simply calculate the Itoh
; continuum for the complete range of wavelengths and temperatures,
; and then replace any zero values with the continuum from Sutherland.
; 
; Previously I had split the temperature range into above and below
; 10^6 K, using Sutherland below and Itoh above. I then also had to
; check for wavelengths outside of Itoh's range and replace
; these with Sutherland's data. A potential problem was if
; DEM_INT was specified and there was only one temperature below 10^6
; (e.g., logT=5.9). The routine sutherland.pro would then reject the
; single temperature value as two or more are required for a DEM
; calculation. The new method thus fixes this.
;
; The new method will be slower than the old but since computers are
; now a lot faster this isn't a major problem.
;
itoh,temp,wvl,int,photons=photons, $
     dem_int=dem_int,min_abund=min_abund,kev=kev,em_int=em_int, $
     abund_file=abund_file, ioneq_file=ioneq_file, sngl_ion=sngl_ion, $
     element=element
;
k=where(int EQ 0.,nk)
IF nk NE 0 THEN BEGIN
  sutherland,temp,wvl,int1,photons=photons, $
             dem_int=dem_int,min_abund=min_abund,kev=kev,em_int=em_int, $
             abund_file=abund_file, ioneq_file=ioneq_file, sngl_ion=sngl_ion, $
             element=element
  int[k]=int1[k]
ENDIF 

siz=size(int)

IF keyword_set(sumt) AND (siz[0] GT 1) THEN int=total(int,2)

IF keyword_set(verbose) THEN BEGIN 
  t2=systime(1)
  print,format='("% FREEFREE: continuum calculated in ",f8.1," seconds")',t2-t1
END  

END
