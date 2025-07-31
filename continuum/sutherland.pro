;+
; NAME:
;	SUTHERLAND
;
; PURPOSE:
;	Calculate the free-free continuum from an hot, low density plasma
;
;       Uses the free-free gaunt factor calculations of Sutherland, 1998, 
;       MNRAS, 300, 321.
;
;       Note that Sutherland's Eq.(15) has units of erg/cm^3/s. Comparing 
;       with Rybicki & Lightman's Eq.5.14(a) (in their book 'Radiative 
;       Processes in Astrophysics'), suggests that Sutherland's units 
;       should be erg/cm^3/s/sr/Hz. We are assuming the latter to be 
;       correct in this routine.
;
;       When using the DEM_INT optional input, FREEFREE expects the 
;       differential emission measure to have been derived from a product 
;       of N_e*N_H rather than N_e*N_e. This can be important when dealing 
;       with a regime (typically T < 10^4.5) where H and He are not fully 
;       ionised.
;
;       To calculate the continuum spectrum the user should specify
;       whether the plasma structure is specified through a differential
;       emission measure (DEM) curve or through emission measure (EM)
;       values. For the former, the DEM values are specified through the
;       DEM_INT keyword and for the latter the EM values are specified
;       through the EM_INT keyword. If neither are specified EM_INT
;       values of unity are assumed.
;
;       If DEM_INT values are specifed they must be specified on a
;       temperature scale that has constant intervals in log10(T). EM_INT
;       values do not have this restriction.
;
;       Note that EM_INT must be used for isothermal plasmas.
;
; CATEGORY:
;       CHIANTI; continuum.
;
; CALLING SEQUENCE:
;       SUTHERLAND,temperature, wavelength, intensity
;
; INPUTS:
;	T	   Temperature in degrees Kelvin, can be a 1 D
;                  array. If the keyword DEM_INT= is given, then the
;                  temperatures must be specified with a fixed spacing
;                  in logT. E.g., TEMP=10.^[4.0,4.1,4.2]. 
;
;       WVL        Wavelength in Angstroms. If /keV is set, then WVL is
;                  assumed to contain energies in keV units.
;
;
; OPTIONAL INPUTS:
;     Sngl_Ion:  A string or string array containing a CHIANTI ion
;                name (or names). If specified, then only the
;                specified ions will be included in the calculation. 
;     Element:   Either the atomic number of an element, or the
;                symbol name (e.g., 'Fe' for iron). If specified, then
;                only the specified element will be included in the
;                calculation. 
;     Abund_File:  The name of CHIANTI abundance file. If not set,
;                then the routine asks you to choose a file with a
;                widget.
;
;     Ioneq_File: The name of an ionization equilibrium file. If not
;                set, then the CHIANTI default, !ioneq_file, is used.
;
;     MIN_ABUND: If set, calculates the continuum only from those 
;                elements which have an abundance greater than min_abund.  
;                  This can speed up the calculations. ;
;
;     DEM_INT:  This should be the same size as TEMP and contain
;               differential emission measure values. Specifying
;               DEM_INT places a restriction on TEMP (see above). Note
;               that DEM_INT values are multiplied by the factor dT
;               within the routine, whereas EM_INT values are not.
;
;     EM_INT:   This should be the same size as TEMP and contain
;               emission measure values. Note
;               that DEM_INT values are multiplied by the factor dT
;               within the routine, whereas EM_INT values are not. If
;               neither EM_INT or DEM_INT are specified, then a EM_INT
;               vector of same size as TEMP with values of unity is
;               assumed.
;
; KEYWORD PARAMETERS:
;	NO_SETUP:  If the procedure setup_elements has already been called 
;                  then the keyword /no_setup should be set to avoid 
;                  repeating this step.
;
;       PHOTONS    The output spectrum is given in photon units rather 
;                  than ergs.
;
;       SUMT       The default is to output the intensity array as an array 
;                  of size (nwvl x nT). Setting this keyword performs a sum 
;                  over the temperatures to yield a vector of same size as 
;                  the input wavelengths, thus producing the complete 
;                  free-free spectrum.
;	
;       KEV        If set, then WVL is assumed to contain energies in keV
;                  units rather than wavelengths in Angstroms.
;
; OUTPUTS:
;	RAD        Free-free continuum intensity in units
;                  10^-40 erg cm^3 s^-1 str^-1 angstrom^-1 per unit 
;                  emission measure [ integral (N_e N_H dh) in cm^-5 ]. 
;                  If T is given as a 1-D array, then RAD will be output 
;                  as a 2-D array, with one element for each temperature 
;                  and wavelength (but also see SUMT).
;
;                  If the keyword /keV is set, then the units of INT will be 
;                  10^-40 erg cm^3/s/sr/keV
;
; PROGRAMMING NOTES:
;     The gaunt factors from Sutherland (MNRAS 300, 321, 1998) are a 
;     function of uu and gg (see his Eq. 14). uu is a function of both 
;     wavelength and T, while gg is a function of T only.
;
;     The gaunt factor (gff) is tabulated for values of uu and gg at 
;     fixed intervals in log(uu) and log(gg). The log(uu) values go from 
;     -4 to 4 in 0.2 steps; the log(gg) values go from -4 to 4 in 0.1 steps.
;
;     A particular (input) temperature and wavelength give rise to values 
;     uu0 and gg0, the logs of which lie between -4 and 4. To derive the 
;     corresponding gff0 value, I use the IDL routine BILINEAR. 
;
;     BILINEAR requires, not uu and gg values as input, but indices. 
;     E.g., the indices corresponding to the tabulated values of uu are 
;     0 (=-4.0), 1 (=-3.8), 2 (=-3.6), etc. Thus, if log(uu0)=-3.76, then 
;     i_uu0=1.20 is the index of uu0.
;
;     In order to make significant time-savings, I give BILINEAR all of 
;     the wavelengths and temperatures in the same call for a particular 
;     ion. To do this, I make my i_uu and i_gg values 2-D arrays of size 
;     (nwvl x nT), and BILINEAR then produces a (nwvl x nT) array 
;     containing the gff values.
;
;     A problem occurred if nT=1, as BILINEAR will turn the input uu and 
;     gg vectors into 2-D arrays of size (nwvl x nwvl). If there are a 
;     large number of wavelengths, this uses a lot of memory. To solve 
;     this I make a 2 element temperature vector whose values are 
;     identical, and then change this back to a 1 element vector after 
;     BILINEAR has been called. See the parts of the code where I use 
;     'tst1'.
;
; CALLS:
;       READ_IP, READ_GFFGU, READ_IONEQ, READ_ABUND, CH_CHOOSE_ABUND
;
; EXAMPLES:
;       IDL> freefree,1.e+6,wvl,int
;       IDL> freefree,1.e+6,wvl,int,min_abund=3.e-5
;       IDL> freefree,1.e+6,wvl,int,/no_setup,min_abund=1.e-6
;
;       IDL> wvl=findgen(5001)/10. + 50.
;       IDL> temp=10.^(findgen(41)/10. +4.)
;       IDL> freefree,temp,wvl,int
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	March 1999:     Version 2.0
;       September 1999:  Version 3.0
;
;       Ver.3.1, 11-Aug-00, Peter Young
;           Improved call to bilinear, allowing routine to solve for 
;           all temperatures in one go. This makes the routine quicker, 
;           and also lowers the memory usage of the routine when dealing 
;           with many wavelengths.
;
;       Ver.3.2, 16-Aug-00, Peter Young
;           Corrected expression for 'gg', replacing ip with z^2.
;
;       Ver.3.3, 16-Oct-00, Peter Young
;           Now deals with dem_int correctly
;
;       Ver.3.4, 10-Oct-01, Ken Dere
;           Corrected for labelling errors in Sutherland's gffgu.dat file
;           No longer reads ionization potential file
;
;       Ver.3.5, 5-Dec-01, Peter Young
;           Corrected expression for gamma^2
;           Renamed routine sutherland.pro
;           Restructured code to make it run quicker.
;
;       Ver.3.6, 22-May-01, Peter Young
;           Re-instated the MIN_ABUND optional input. 
;           Changed ioneq_t to ioneq_logt (GDZ).
;
;       Ver.3.7, 18-Aug-03, Peter Young
;           Activated /PHOTONS keyword
;
;       Ver.3.8, 5-Nov-03, Peter Young
;           Corrected bug found by Jim McTiernan when multiple temperatures
;           were input. The quantity 'newfactor' was not being calculated
;           correctly due to indexing problems.
;
;       Ver. 4, 9-Mar-2006, Peter Young
;           Added /keV keyword.
;
;       Ver.5, 10-Sep-2009, Peter Young
;            Routine now checks that temperatures are evenly spaced
;            when DEM_INT is specified. Also 0.1 dex spacing is no
;            longer assumed.
;
;       Ver.6, 6-Oct-2009, Peter Young
;            Introduced EM_INT keyword for specifying emission measure
;            values for isothermal plasmas.
;
;       Ver.7, 5-Sep-2012, Peter Young
;            Updated header (no change to code).
;
;       Ver.8, 26-Apr-2019, Peter Young
;            Added SNGL_ION, ELEMENT, IONEQ_FILE and ABUND_FILE
;            optional inputs; removed common block.
;
;       Ver.9, 19-Aug-2019, Peter Young
;            Removed references to /no_setup; tidied up header. 
;
;       Ver.10, 07-Dec-2023, Peter Young
;            For conversion to keV, I changed indgen to lindgen due to
;            errors if there are too many wavelength bins.
;
;       Ver.11, 31-Jul-2025, Peter Young
;            Now calls get_ieq for interpolating the ion fractions.
;-


pro sutherland,t,wvl,rad,min_abund=min_abund, $
               dem_int=dem_int, photons=photons, sumT=sumT, kev=kev, em_int=em_int, $
               sngl_ion=sngl_ion, element=element, $
               abund_file=abund_file, ioneq_file=ioneq_file

if n_params(0) lt 3 then begin
  print,'Use: IDL> freefree,temperature,wavelength,intensity  [,/no_setup, min_abund='
  print,'                     dem_int= , em_int= , /photons, /sumT, /kev, abund_file='
  print,'                     ioneq_file=, element=, sngl_ion= ]'
  return
endif


kb=1.38062d-16   ; erg/K
h=6.6262d-27     ; erg s
c=2.997925d+10   ; cm/s
ryd=2.17992d-11  ; erg
factor=5.44436d-39  ; this is f_k of Eq.15 of Sutherland
rescale=1.d+40
;

;
; Use a widget to select the abundance file if not specified by abund_file.
;
IF n_elements(abund_file) EQ 0 THEN BEGIN
  abund_file=ch_choose_abund()
ENDIF
read_abund,abund_file,abund,ab_ref, element=element

;
; Always use the default ioneq file, unless specified by user.
;
IF n_elements(ioneq_file) EQ 0 THEN BEGIN
  ioneq_file=!ioneq_file
ENDIF
read_ioneq,ioneq_file,ioneq_logt,ioneq, ioneq_ref, element=element, sngl_ion=sngl_ion



zmax=max(where(abund gt 0.))+1

IF n_elements(min_abund) EQ 0 THEN min_abund=0.

n_ioneq_logt=n_elements(ioneq_logt)
;dlnt=alog(10.^(ioneq_logt(1)-ioneq_logt(0)))
;
;  read free-free gaunt factors
;
read_gffgu,g2,u,gff
;

t = double(t)
wvl = double(wvl)

;
; * keV units *
;
IF keyword_set(kev) THEN BEGIN
  wvl_save=wvl
  efact=12.398/wvl^2   ; hc/E^2 factor to multiply continuum by
  efact=reverse(efact)
  wvl=12.398/wvl   ; convert energies to angstroms
  wvl=reverse(wvl)
ENDIF ELSE BEGIN
  efact=1d0
ENDELSE

;
; If t and wvl only have 1 element, then need to convert them to arrays
;
IF n_elements(wvl) EQ 1 THEN wvl = [wvl]
IF n_elements(t) EQ 1 THEN t = [t]

nwvl=n_elements(wvl)
nt=n_elements(t)

;
; Note that wvl has been converted to Angstroms stage, so no correction for
; the /kev keyword is required.
; 
IF keyword_set(photons) THEN BEGIN
  erg2phot=( wvl # (dblarr(nt) + 1.) ) / 1.9864d-8
ENDIF ELSE BEGIN
  erg2phot= (dblarr(nwvl) + 1. ) # (dblarr(nt) + 1.)
ENDELSE

ident_t = make_array(max([nt,2]),val=1.,/double)
ident_wvl = make_array(nwvl,val=1.,/double)

;
; DEM_INT and EM_INT implementation
; ---------------------------------
;  If DEM_INT is specified then it's necessary to multiply by the
;  DEM values by the factor dT (=ln(10)*T*d(log10(T)).
;
;  I implement this by introducing a keyword dT. If dT=1 then a DEM is
;  assumed, while dt=0 implies EM values. In the latter case DEM_INT
;  is set to be EM_INT.
;
;  EM_VALS is set to be DEM_INT or EM_INT and is used in the rest of
;  the code.
;
CASE 1 OF 
  n_elements(dem_int) NE 0 AND n_elements(em_int) NE 0: BEGIN
    print,'%SUTHERLAND: Please specify only one of DEM_INT and EM_INT. Returning...'
    rad=0.
    return
  END
 ;
  n_elements(dem_int) EQ 0 AND n_elements(em_int) NE 0: BEGIN
    em_vals=em_int
    dt=0
  END 
 ;
  n_elements(dem_int) NE 0 AND n_elements(em_int) EQ 0: BEGIN
    em_vals=dem_int
    dt=1
  END 
 ;
  n_elements(dem_int) EQ 0 AND n_elements(em_int) EQ 0: BEGIN
    em_vals=dblarr(nt)+1.
    dt=0
  END 
ENDCASE

;
; Make sure EM_VALS has same size as TEMP
;
em_vals=double(em_vals)
IF n_elements(em_vals) NE nt THEN BEGIN
  print,'%SUTHERLAND: Warning, number of elements of DEM_INT (EM_INT) must match the'
  print,'  number of temperatures. Returning...'
  int=0.
  return
ENDIF 



;
; Evaluate the quantity 'demfactor', this is simply the factor dT and
; so if /dt is not set then demfactor is set to 1.
;
IF keyword_set(dt) THEN BEGIN
  IF nt LT 2 THEN BEGIN
    print,'%SUTHERLAND: Two or more temperatures must be specified if DEM_INT is specified.'
    print,'             Returning...'
    rad=0.
    return 
  ENDIF 
 ;
  lt1=alog10(t[0:nt-2])
  lt2=alog10(t[1:*])
  mean_lt=mean(lt2-lt1)
  max_LT=max(lt2-lt1)
  min_LT=min(lt2-lt1)
  IF max_LT-mean_LT GE 0.05*mean_LT OR mean_LT-min_LT GE 0.05*mean_LT THEN BEGIN
    print,'%SUTHERLAND: Error, LOG10 TEMP needs to be evenly spaced if DEM_INT is specified.'
    print,'             E.g., LOG10 TEMP should be spaced at 0.1 dex intervals.'
    rad=0.
    return
  ENDIF 
  demfactor=t*mean_lt/alog10(exp(1))
END ELSE BEGIN
  demfactor = make_array(max([nt,2]),val=1.,/double) 
ENDELSE 



IF nt EQ 1 THEN rad = dblarr(nwvl,2) ELSE rad=dblarr(nwvl,nt)

IF n_elements(t) EQ 1 THEN temp=[t,t] ELSE temp=t

uu = h*c*1.d8/kb*(1/wvl)#(1/temp)
i_uu = (ALOG10(uu)+4.d)*10.d

elt_i=where(abund[0:zmax-1] GE min_abund)
nei=n_elements(elt_i)

FOR z=1,30 DO BEGIN
  gg = z^2*ryd/kb* ident_wvl #(1/temp)
  i_gg = (ALOG10(gg)+4.d)*5.d
  big = bilinear(gff,i_gg,i_uu)
  newbig=big*exp(-uu)*z^2*( (1/wvl^2)#ident_t )

  FOR iz=0,nei-1 DO BEGIN
    this_abund=abund[elt_i[iz]]

    yi=get_ieq(temp,elt_i[iz]+1,z+1,ioneq_logt=ioneq_logt, $
               ioneq_frac=ioneq)

    newfactor=this_abund*yi*factor*c*1d8*rescale/sqrt(temp)* $
              em_vals * demfactor

    rad=rad+(ident_wvl#newfactor)*newbig*erg2phot

  ENDFOR

ENDFOR

;
; multiply int by the hc/E^2 factor if /kev set
;
IF keyword_set(kev) THEN BEGIN
  e_arr=efact#(dblarr(nt)+1d0)
  rad=rad*e_arr
 ;
 ; reverse the wavelength dimension to match energy ordering rather than
 ; wavelength ordering
  rad=rad[reverse(lindgen(nwvl)),*]
 ;
 ; and set wavelengths back to energy units
  wvl=wvl_save
ENDIF

IF n_elements(t) EQ 1 THEN BEGIN
  rad=rad[*,0]
ENDIF ELSE BEGIN
  IF keyword_set(sumt) THEN rad=total(rad,2)
ENDELSE


END
