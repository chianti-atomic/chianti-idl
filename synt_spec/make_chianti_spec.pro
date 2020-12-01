PRO make_chianti_spec, TRANSITIONS, LAMBDA, OUTPUT, BIN_SIZE=BIN_SIZE,  $
        INSTR_FWHM=INSTR_FWHM,  BINSIZE=BINSIZE, $ 
        WRANGE=WRANGE, ALL=ALL, continuum=continuum, $
        ABUND_NAME=ABUND_NAME, MIN_ABUND=MIN_ABUND, $
        photons=photons, file_effarea=file_effarea, $
        err_msg=err_msg,  verbose=verbose, kev=kev, $
        no_thermal_width=no_thermal_width, lookup=lookup

on_error,0
;+
; PROJECT     :  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory (USA), the University of Florence (Italy), the
;       University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
;                   
; NAME        : MAKE_CHIANTI_SPEC
;     		          
; PURPOSE     : 
;              To create a CHIANTI synthetic spectrum 
;               
; CALLING SEQUENCE:
;
;       IDL> make_chianti_spec, TRANSITIONS,  LAMBDA, SPECTRUM,$ 
;                    [BIN_SIZE= ,  ,INSTR_FWHM= , PIXEL=PIXEL, BINSIZE = BINSIZE, $
;                    WRANGE= , ALL=ALL, continuum=continuum, $
;                    ABUND_NAME= , MIN_ABUND=, photons=photons,
;                    effarea=effarea, $
;                    /no_thermal_width
;
;
; PROCEDURE : 
; 		
;     From information contained in the structure TRANSITIONS, constructs 
;     a synthetic spectrum
;
;     By default, routine assumes thermal widths for lines. This can
;     be switched off with /no_thermal_width
;
;   PROGRAMMING NOTES
;
;     The line profile is constructed using the IDL gaussint routine. 
;     For a wavelength pixel centred at l and with width dl, gaussint 
;     is used to integrate the Gaussian up to l-dl/2 and up to l+dl/2. 
;     The difference between the two is the intensity in this pixel.
;
;
;    
; INPUTS      : 
;		
;               TRANSITIONS, the structure created by ch_synthetic.pro.
;               
; OPT. INPUTS : 
;
;     LAMBDA   Array of wavelengths (X-values). If not defined as input, it is
;              calculated on the basis of BIN_SIZE, and returned as an output. 
;              If defined as input, the routine checks that there are at least
;              10 points in the wavelength range defined by WRANGE. If there
;              are, the corresponding subset of LAMBDA is returned, otherwise
;              the routine exits with an error. If /keV is set, then
;              LAMBDA should be given in keV units. IMPORTANT: LAMBDA
;              is re-defined by make_chianti_spec, so (to avoid
;              unexpected results) a "fresh" LAMBDA
;              should always be created prior to calling MAKE_CHIANTI_SPEC.
;
;     BIN_SIZE      Bin size  in Angstroms of the spectrum to be created. A linear
;              spectrum is assumed. In case an effective area file is used, the
;              wavelenghts in the file (that might not be linear) are used to
;              create the spectrum, and this bin size looses any meaning.
;
;     WRANGE   Allows a subset of the wavelength range to be turned into 
;              a spectrum. When using syn_plot, this speeds up the routine 
;              a lot if you've elected to zoom in on a region.
;
;     INSTR_FWHM Instrumental FWHM (Angstroms). 
;                In case an effective area file is used, The line intensities
;                contributing to each bin are summed, and subsequently convolved
;                with a gaussian of full-width-half-maximum FWHM if FWHM is not
;                set = 0 . Please make sure that the FWHM value (if not set to
;                zero) is larger than  the bin size. 
;
;     ABUND_NAME  A CHIANTI abundance file name can be set. 
;                It allows the abundance file given in transitions.abund_name
;                (if present)   to be over-ridden. Note that it also used for
;                the continuum calculation.
;
;     MIN_ABUND: If set, calculates line intensities only from those elements
;                  which  have an abundance greater than min_abund. 
;                  Can speed up the calculations. For example, from Allen
;                  (1973):
;
;                   abundance (H)  = 1.
;                   abundance (He) = 0.085
;                   abundance (C)  = 3.3e-4
;                   abundance (O)  = 6.6e-4
;                   abundance (Si) = 3.3e-5
;                   abundance (Fe) = 3.9e-5
;
;     FILE_EFFAREA
;                   The Effective Area File (TWO COLUMNS - wavelengths in
;                   Angstroms and cm^2 values) If defined, then the spectrum is
;                   multiplied with these values.  Any input  LAMBDA value is
;                   over-written by the wavelenghts in the file (that might not
;                   be linear) and  are used to create the spectrum.
;		    Note that this option only works well if a sufficient number
;		    of bins is given. The line intensities contributing to each
;		    bin are summed, and  subsequently convolved with a gaussian
;		    of full-width-half-maximum FWHM, if FWHM is not set = 0.
;                   Please note that the convolution might not work if a small
;                   number of  bins is defined. 
;
;               
; OUTPUTS     : 
;		
;		LAMBDA  Array of wavelengths (X-values).
;              If not defined as input, it is
;              calculated on the basis of BIN_SIZE, and returned as an output. 
;              If defined as input, the routine checks that there are at least
;              10 points in the wavelength range defined by WRANGE. If there
;              are, the corresponding subset of LAMBDA is returned, otherwise
;              the routine exits with an error.
;
;
;               SPECTRUM  A structure containing all the information:
;
;                     LAMBDA      The array of X-values
;                     SPECTRUM    The array of Y-values
;                     UNITS       The units of LAMBDA, SPECTRUM
;                     INSTR_FWHM  The Instrumental FWHM
;                     BIN_SIZE    Width of the Bins  (fixed) in angstroms
;                     ABUND_NAME  The CHIANTI abundance file name             
;                     ABUND       The abundance values
;                     MIN_ABUND   The minimum abundance value used                 
;                     ABUND_REF   The references
;                     CONTINUUM   The values of the continuum (if calculated)
;                     
;                     FILE_EFFAREA The Effective Area File used (optional)
;                     EFFAREA       The array of effective area values
;                                 (optional - same size of LAMBDA)
;
;                    .LINES      An array of structures, for all the lines used               
;                                to calculate the SPECTRUM. 
;                                The tags are the same as those created by 
;                                CH_SYNTHETIC, plus
;                       .PEAK    The peak intensity of the line in the spectrum
;                                (approx. value) 
;	
; OPT. OUTPUTS:
;		
;     BINSIZE  If BIN_SIZE  is not  specified, then the spectrum 
;              bin-sizes are computed automatically, and the size of the 
;              bin returned in BINSIZE.
;
;
; KEYWORDS    : 
;
;     PIXEL    The spectrum is given in /pixel units rather /ang
;        (DISABLED)
;      
;     ALL      Add  lines that originally had negative wavelengths  
;               
;     PHOTONS  If set=1, the output intensities will be in photons instead of 
;                  ergs.
;
;     CONTINUUM
;              If set, then the  continuum is added to the 
;              spectrum.
;
;     VERBOSE  If set, then print out information while the routine is
;              running.
;
;     KEV      If set, then the output spectrum is in units 
;              erg cm^-2 s^-1 keV^-1.
;
;     NO_THERMAL_WIDTH  If set, then the lines will not be broadened
;                       by the thermal width. If both
;                       /no_thermal_width and instr_fwhm=0, then an
;                       emission line will be concentrated in a
;                       single pixel within the spectrum.
;
;     LOOKUP:   This keyword is passed on to two_photon in order to
;               speed the calculation
;
; CALLS       : 
;		
;		PRY:     	GET_ATOMIC_WEIGHTS
;		Chianti: FREEBOUND, FREEFREE
;		
; COMMON       (with freefree freebound and two_photon):
; 		
;		elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref
;
; RESTRICTIONS: The input structure has to be of the type created by ch_synthetic.
;               The LAMBDA, EFFAREA values must be ordered in wavelength and the
;               LAMBDA values must be in Angstroms.
;               
; SIDE EFFECTS: None known yet.
;               
;
; EXAMPLES    : 
;		
;		make_chianti_spec, output_struct,  LAMBDA, SPECTRUM,$
;		 bin_size=0.01, instr=0.1 
;
; CATEGORY    : 
;               spectral synthesis.
;
; PREV. HIST. :
;
;      
; WRITTEN     : 
;           Peter Young , CfA, pyoung@cfa.harvard.edu  1-Sept-2000
;
; MODIFICATION HISTORY:  
; 
;               Version 1, PRY 1-Sept-2000
;
;               Version 2, Giulio Del Zanna (GDZ)  10-Oct-2000
;
;               put ALL keyword, removed the FWHM obsolete and
;               confusing call. Reorganised various minor things.
;
;               Version 3, PRY 19-Oct-2000
;                 Corrected the way continuum is treated for an isothermal 
;                 spectrum.
;
;               V. 4, 2-Nov-2001 GDZ. Now MIN_ABUND is effective not only in the
;               continuum calculation but also in the line spectrum.
;               Modified for the use of logt_isothermal
;
;               V.5, GDZ, added EFFAREA keyword: an ascii file with lambdas and
;                                              effective areas can be read
;                                              in. The line intensities are
;                                              calculated in a different way.
;                        Also, changed the output.
;
;               V.6, GDZ, 28-Apr-02 redefined completely the OUTPUT structure. 
;                    Major revision (added two_photon verbose).
;
;               V.7, GDZ, 3-May-2002
;                    fixed  a bug, when negative angpix values occur.
;                     
;               V.8, GDZ, 22 May 2002,  changed some tags of the output, and
;                    added min_abund in the continuum call.
;
;               V.9, GDZ, 30-May-02 replaced fix() with round() 
;
;               V. 10, 15-July-2002 , GDZ 
;                    changed the location of Effective area files.
;
;               V.11 14-Aug-02, GDZ 
;                    speeded up the routine, by changing the way the PEAK tag is
;                    added to the structure. The drawback is that only the
;                    'standard v.4 tags' are allowed, and future additions will
;                    have to be dealt properly.
;  
;               v.12 2-Dec-2002, GDZ. 
;                   Fixed a bug:  Removed the plotting of the window with the effective areas.
;
;               v.13 26-Apr-2005, Enrico Landi (EL)
;                   Fixed a minor bug: if the lines were more than 32768 (2^15), the main
;                                loop crashed.
;
;               v.14 22-Jul-2005 GDZ 
;                 -fixed a bug. When the routine was run once without
;                 defining the lambdas, and then with the lambdas
;                 defined (the units were switched to photons)
;                 -fixed a bug. When the effective areas were used,
;                 all lines were used to create the spectrum.
;                 -added hard-wired switch to photons when using
;                 effective area files.
;
;                 -added the keV option
;
;                 -now can output a spectrum only with the continuum
;                 (i.e. even if no emission lines are present). 
;
;               v.15, 2-Aug-2005, GDZ 
;                 Added a check on the input structure. If it was
;                 calculated with ch_synthetic and the keyword
;                 /no_sum_int, it cannot be used here.
;
;               v.16, 31-Aug-2005, GDZ
;                 Fixed a typo (keyowrd_set).
;
;               v.17, 10-Mar-2006, Peter Young
;                 The /keV keyword is now passed to the continuum routines,
;                 rather than make_chianti_spec performing the conversion
;                 itself.
;
;               v.18, 4-May-2006, Peter Young
;                 Changed LAMBDA to be a double array in order to correct a
;                 problem with line profiles for very small bin sizes.
;                 Problem pointed out by Matthew West (Imperial).
;
;               v.19, 16-Nov-2007, Vincenzo Andretta, INAF/OACN (andretta@oacn.inaf.it)
;                 Fixed a problem with units when LOGT_ISOTHERMAL is set.
;
;               v.20, 7-Oct-2009, Peter Young
;                 Modified call to continuum routines to use new
;                 EM_INT keyword (for isothermal plasmas).
;                 
;               v.21, 6-Oct-2010, Peter Young
;                 Corrected bug related to the v.20 modification.
;
;               v.22, 25-Jul-2011, Peter Young
;                 Changed units for the isothermal case back to
;                 "cm-2 sr-1 s-1" instead of "cm+3 sr-1 s-1". This
;                 reverses the change made in version 19.
;
;               v.23, 26-Jun-2012, Peter Young
;                 Added /no_thermal_width option.
;
;               v.24, 27-Jul-2012, Peter Young
;                 Changed 'bigpickfile' call to dialog_pickfile.
;
;               v.25, 23-Aug-2012, Peter Young
;                 Fixed bug: dialogpickfile -> dialog_pickfile
;
;               v.26, 17-Aug-2015, Peter Young
;                 Changed 'lambda()' to 'lambda[]' to prevent crashes
;                 due to new IDL lambda.pro routine (IDL 8.4).
;
;               v.27, 26-Nov-2018, Peter Young
;                 Added information to header (about LAMBDA keyword);
;                 no change to code.
;
;               v.28, 26-Apr-2019, Peter Young
;                 Modified call to freefree (no longer uses common
;                 block).
;
;               v.29 9-Aug-2019, Peter Young
;                 Modified call to two_photon (no longer uses common
;                 block).
;
;               v.30, 11-Mar-2020, Peter Young
;                 Added ioneq_file= and abund_file= in call to
;                 freebound.
;
;               v.31, 17-Nov-2020, Peter Young
;                 Added keyword /lookup.
;-


;common with the continuum routines:
;common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref

IF  n_params() lt 3 then begin
   print,' type> make_chianti_spec, output_struct,  LAMBDA, SPECTRUM,$ '
   print,'         [BIN_SIZE= , KEV= ,INSTR_FWHM= , PIXEL=PIXEL, $'
   print,'       WRANGE= , ALL=ALL, continuum=continuum, $'
   print,'     ABUND_NAME= , MIN_ABUND=, BINSIZE = BINSIZE effarea=, /LOOKUP] '
END


IF n_elements(lambda) NE 0 THEN lambda=double(lambda)

t0=systime(1)

defsysv,'!xuvtop', EXISTS = EXISTS 
IF NOT EXISTS THEN $
  message, 'system variable !xuvtop must be set  '

xuvtop = !xuvtop

;check that we have a correct input structure:

result = ch_check_str(TRANSITIONS, /int)
if not result then begin 
      err_msg='% MAKE_CHIANTI_SPEC: incorrect input structure! - EXIT'
      return 
  endif

;check that we did not use the /no_sum_int keyword

if n_elements(TRANSITIONS.lines[0].int) ne 1 then begin 
      err_msg='% MAKE_CHIANTI_SPEC: incorrect input structure! - EXIT'
      return 
end 

ang2kev=12.39854

units = strarr(2)

IF N_ELEMENTS(bin_size) EQ 0 THEN bin_size = 0.
IF N_ELEMENTS(kev) EQ 0 THEN kev = 0.

IF N_ELEMENTS(INSTR_FWHM) EQ 0 THEN INSTR_FWHM = 0.

IF N_ELEMENTS(min_abund) EQ  0  THEN  min_abund = 0.

IF n_elements(photons) EQ 0 THEN photons = 0

nlines = N_ELEMENTS(transitions.lines)

;GDZ- kev option
IF N_ELEMENTS(wrange) EQ 0 THEN begin 
   wrange = transitions.wvl_limits

; convert  angstroms to keV

   if keyword_set(kev) then begin 
      wrange=reverse(ang2kev/wrange)
   end 
endif else begin 

;minimal check :

   if not valid_num(wrange[0]) or  not valid_num(wrange[1]) then begin 
      err_msg='% MAKE_CHIANTI_SPEC: X-ranges are not valid numbers! - EXIT'
      return 
   endif else begin   
      if wrange[0] ge wrange[1] then begin 
         err_msg='% MAKE_CHIANTI_SPEC: X-ranges are not valid! - EXIT'
         return 
      end 
   end 
end  

IF  n_elements(file_effarea) GT 0 THEN BEGIN 


   IF NOT file_exist(file_effarea) THEN BEGIN 

      dir = concat_dir(concat_dir(!xuvtop,'ancillary_data'), 'instrument_responses')

      IF DIR_EXIST(dir) THEN path = dir  ELSE BEGIN 
         cd, current=dir
         path = dir
      END 

      
      ptitle='Select Effective Area File (TWO COLUMNS - wavelengths in Angstroms and cm^2 values )'
      file_effarea =  dialog_pickfile(title=ptitle, path=path,filter='*area')
      ;; file_effarea =  bigpickfile(title=$
      ;;                             'Select Effective Area File (TWO COLUMNS - wavelengths in Angstroms and cm^2 values )', $
      ;;                             path= path, filter='*area' )

   END 

   IF file_exist(file_effarea) THEN BEGIN 

      data = read_ascii (file_effarea)

      lambda_effarea = reform( data.field1(0, *))
      lambda = lambda_effarea
      effarea = reform( data.field1(1, *))

; lambda must be sorted in increasing order.
      lambda=lambda[sort(lambda)]
      effarea =effarea(sort(lambda))

; convert to energy.
      if keyword_set(kev) then begin 
         lambda=ang2kev/lambda
         lambda=lambda[sort(lambda)]
         effarea =effarea(sort(lambda))
      end 

      bin_size = 0.

   ENDIF  ELSE BEGIN  
      err_msg='% MAKE_CHIANTI_SPEC: Wrong  effarea file -- EXIT'
      return 
   END  

   if not keyword_set(photons) then begin 
      print,'% MAKE_CHIANTI_SPEC: switching to photons'
      photons=1
   end 

END         



IF keyword_set(kev)  THEN BEGIN

   units(0) = 'keV'

; I assume for now that the units in the input structure are 'Angstroms'
; IF transitions.wvl_units EQ 'keV' THEN stop

;   wlimits = REVERSE(1.d8/wlimits*8065.54d)
   
ENDIF ELSE units(0) = 'Angstroms'

;
; angpix controls whether the spectrum is given in /ang or /pix units.
;

IF n_elements(LAMBDA) LE 1 THEN BEGIN 

;define the X-array LAMBDA
;--------------------------

   IF bin_size EQ 0. THEN BEGIN
;
; the following turns a binsize of, e.g., 0.0347862 into 0.03
;
      bs = (wrange[1]-wrange[0])/200.
      bs = STRING(format='(e6.0)',bs)
      bs = FLOAT(bs)

   ENDIF ELSE bs = bin_size

   lambda = dindgen ( ROUND((wrange[1]-wrange[0])/bs +1 ) )*bs $
     + wrange[0]

   in = where(lambda GE wrange[0] AND lambda LE  wrange[1], ng)
   lambda = lambda[in]

   angpix = dblarr( n_elements(LAMBDA))
   binsize = dblarr( n_elements(LAMBDA))

   binsize(*) = bs

   IF KEYWORD_SET(pixel) THEN angpix(*) = 1. ELSE angpix(*) = binsize

ENDIF ELSE BEGIN 

;LAMBDA is already defined.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;check the ranges.

   in = where(lambda GE wrange[0] AND lambda LE  wrange[1], ng)

;allow at least say 2 points:

   IF ng GT   2 THEN BEGIN 
      lambda = lambda[in]

      IF  n_elements(file_effarea) GT 0 THEN  effarea =  effarea(in)

   ENDIF  ELSE BEGIN 
      delvarx, lambda,  effarea
      err_msg='%MAKE_CHIANTI_SPEC: Conflict between X-axis and X-range, less than 2 points in the spectrum -- EXIT'
      return 
   END  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;check if bin etc were defined. OVERRIDES 

   angpix = fltarr( n_elements(LAMBDA))
   binsize = fltarr( n_elements(LAMBDA))

;all the bin-beginning values and the final bin-ending values:

   grid =mid2bound(lambda)

;watch for negative values:
   angpix = 2.*(grid(1:*)-lambda)


;define a finer, linear  grid:

   bad = where(angpix LE 0, nb)
   good = where(angpix GT  0)

   IF nb GT 0 THEN BEGIN 
      dlambda = average(angpix[good])/5.
;angpix(bad) = 0.001
   ENDIF ELSE    dlambda = min(angpix)/5. 

   dlambda = dlambda < 0.005

   npoints = round((max(lambda)-min(lambda))/dlambda )+1

   nlambda = min(lambda)+findgen(npoints)*dlambda
   ngrid = mid2bound(nlambda)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;GDZ
   binsize[0:n_elements(LAMBDA)-2]= lambda[1+indgen(n_elements(LAMBDA)-1)]-lambda[indgen(n_elements(LAMBDA)-1)]
;   FOR i=0L, n_elements(LAMBDA)-2 DO $ 
;     binsize[i] =  lambda[i+1]-lambda[i]

   binsize(n_elements(LAMBDA)-1) = binsize(n_elements(LAMBDA)-2)

END  

spectrum = lambda-lambda


;
; Abundances. 'line_abunds' contains the element abundance for each line 
; in the spectrum.
;

IF n_elements(abund_name) EQ 0 THEN BEGIN 

   IF tag_exist(transitions, 'abund_name') THEN BEGIN 

      abund_name = transitions.abund_name 
      abund = transitions.abund
      abund_ref = transitions.abund_ref

   ENDIF  ELSE BEGIN 

; Choose abundance file
;
      dir=concat_dir(!xuvtop,'abundance')
      abund_name=ch_get_file(path=dir,filter='*.abund',title='Select Abundance File')

      ff = findfile(abund_name)
      IF   ff(0) NE '' THEN $
        read_abund,abund_name,abund,abund_ref ELSE BEGIN 
         print, 'Error,  no abundance  file found !'
         return
      END 

   END  
ENDIF ELSE BEGIN 
   ff = findfile(abund_name)
   IF   ff(0) NE '' THEN $
     read_abund,abund_name,abund,abund_ref ELSE BEGIN 
      print, 'Error,  no abundance  file found !'
      return
   END 
END 


line_abunds = abund[transitions.lines.iz-1]
wvl = transitions.lines[*].wvl
intensities = transitions.lines[*].int * line_abunds[*] 
flag=transitions.lines[*].flag
tmax=transitions.lines[*].tmax
iz=transitions.lines[*].iz

IF keyword_set(kev)  THEN    wvl=ang2kev/wvl

;    wvl=wvl[sort(wvl)]
;    line_abunds=line_abunds[sort(wvl)] 
;    intensities =intensities[sort(wvl)] 
;    flag=flag[sort(wvl)]
;    tmax=tmax[sort(wvl)]
;    iz=iz[sort(wvl)]
;end 

IF keyword_set(ALL) THEN $
  index=where((wvl LE wrange[1]) AND $
              (wvl GE wrange[0]) AND  (line_abunds[*] GT min_abund), nl) ELSE $
  index=where((flag EQ 0) AND (wvl LE wrange[1]) AND $
              (wvl GE wrange[0]) AND  (line_abunds[*] GT min_abund) , nl)



;NOW we switch intensity units between ERGS and PHOTONS

IF photons THEN BEGIN 
 ; PRY, 25-Jul-2011: I've switched back to the same units for both DEM and
 ; isothermal as it's not clear why this change was made. (The
 ; spectrum units should be the same for both the DEM and isothermal
 ; case.) 
   int_units = 'photons cm-2 sr-1 s-1'
   ;;(16/Nov/2007 [VA] - START)
;   IF tag_exist(transitions,'DEM')  THEN BEGIN 
;     int_units = 'photons cm-2 sr-1 s-1'
;   ENDIF ELSE BEGIN
;     int_units = 'photons cm+3 sr-1 s-1'
;   ENDELSE
   If STRPOS(STRLOWCASE(transitions.int_units),'erg') GE 0 THEN BEGIN
;convert the intensities from erg to photons
     IF keyword_set(kev)  THEN intensities = intensities / (1.986e-8/(ang2kev/wvl)) else $
       intensities = intensities / (1.986e-8/wvl)
   ENDIF
   ;;CASE transitions.int_units OF
   ;;   'erg cm-2 sr-1 s-1': BEGIN 
;convert the intensities from erg to photons
   ;;      IF keyword_set(kev)  THEN intensities = intensities / (1.986e-8/(ang2kev/wvl)) else $
   ;;        intensities = intensities / (1.986e-8/wvl)
   ;;
   ;;   END 
   ;;   'photons cm-2 sr-1 s-1': BEGIN 
   ;;   END 
   ;;ENDCASE 
ENDIF ELSE BEGIN 
 ; PRY, 25-Jul-2011: I've switched back to the same units for both DEM and
 ; isothermal as it's not clear why this change was made. (The
 ; spectrum units should be the same for both the DEM and isothermal
 ; case.) 
   int_units = 'erg cm-2 sr-1 s-1'
;;    IF tag_exist(transitions,'DEM')  THEN BEGIN 
;;      int_units = 'erg cm-2 sr-1 s-1'
;;    ENDIF ELSE BEGIN
;;      int_units = 'erg cm+3 sr-1 s-1'
;;    ENDELSE
   IF STRPOS(STRLOWCASE(transitions.int_units),'erg') LT 0 THEN BEGIN
;convert the intensities
     IF keyword_set(kev)  THEN intensities = intensities*(1.986e-8/(ang2kev/wvl)) else $ 
       intensities = intensities*(1.986e-8/wvl)
   ENDIF
   ;;CASE transitions.int_units OF
   ;;   'erg cm-2 sr-1 s-1': BEGIN 
   ;;   END 
   ;;   'photons cm-2 sr-1 s-1': BEGIN 
;convert the intensities
   ;;      IF keyword_set(kev)  THEN intensities = intensities*(1.986e-8/(ang2kev/wvl)) else $ 
   ;;        intensities = intensities*(1.986e-8/wvl)
   ;;
   ;;   END 
   ;;ENDCASE 
   ;;(16/Nov/2007 [VA] - END)
END 

;remove the 4!pi factor:
;IF  n_elements(file_effarea) GT 0 THEN intensities = intensities*4.0*!pi


;if there are lines:

IF nl NE  0  THEN BEGIN 

   peak = fltarr(nl)

;	rebin to form theoretical line spectrum
; we assume that the lambda is the same as the eff


; lambda must be sorted in increasing order.

;all the bin-beginning values and the final bin-ending values:
;grid =mid2bound(lambda)



   IF n_elements(effarea) GT 0 THEN BEGIN 

;convolve with the PSF:

;spectrum in units/bin :
;GDZ: fixed bug- only lines included in [index]
;
      gspectrum =  hastogram(wvl[index], ngrid, wts=intensities[index])

      IF INSTR_FWHM GT 0 THEN BEGIN 

         sigma = INSTR_FWHM/SQRT(8.*ALOG(2))

;get the  sigma in pixels:

         sigma = sigma/(ngrid[1]-ngrid[0]) ;   average(angpix)

         nx = 2*round(7*sigma/2.)+1
         x0 = round(7*sigma/2.)


         IF x0+3*sigma GT n_elements(ngrid) THEN BEGIN 
            err_msg = 'error, FWHM value too large'
            return
         END 


         broad = gauss_put(nx,0., 1.0, x0, sigma)

         raw = gspectrum
         nel = n_elements(gspectrum)

         temp = [fltarr(x0),raw,fltarr(x0)]
         temp = convol(temp, broad,  total(broad))

;now we have to rebin spectrum over lambda:

         spectrum = rebinw(temp(x0+1:x0+1+nel-1), nlambda, grid, /perbin)

;GDZ: fixed bug- only lines included in [index]:
      ENDIF  ELSE spectrum = hastogram(wvl[index], grid, wts=intensities[index])
;spectrum = rebinw(gspectrum, nlambda, lambda, /perbin)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;spectrum in units/Angstrom:
;   spectrum = spectrum/angpix

      spectrum = spectrum*effarea
      IF nb GT 0 THEN spectrum(bad) = 0.0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      FOR jj=0L,long(nl)-1 DO BEGIN

         i = index(jj)

;get the nearest effective area value to find the 
; approximate new intensity values.

         dummy = min(abs(lambda-wvl[i]), ix)
         ix = ix[0]
         IF ix EQ -1 THEN BEGIN 
            intensities[i]=0.
            peak[jj]= 0.
         ENDIF ELSE BEGIN 

            intensities[i]=intensities[i] * effarea[ix]

            IF INSTR_FWHM GT 0 THEN BEGIN 

               tm = 10.^(tmax[i])
               pp = get_atomic_weight(iz[i])

               IF keyword_set(kev)  THEN begin 
                  wavelength=ang2kev/wvl[i] 
                  x_fwhm=ang2kev/INSTR_FWHM
               endif   else begin 
                  wavelength=wvl[i]
                  x_fwhm=INSTR_FWHM
               end 

; fwhm in Angstroms:
               IF keyword_set(no_thermal_width) THEN BEGIN 
                 fwhm=x_fwhm
               ENDIF ELSE BEGIN 
                 fwhm = SQRT( 5.093d-13 * tm * wavelength^2 / $
                              pp + X_FWHM^2)
               ENDELSE 

               IF keyword_set(kev)  THEN fwhm=ang2kev/fwhm


               peak[jj]= intensities[i] / fwhm / 1.0645 * angpix[ix]

            ENDIF ELSE      peak[jj]= intensities[i] 

         END  
      ENDFOR  


; end of effective area case.
   ENDIF  ELSE BEGIN 

      cc = SQRT(2.*!pi)/2./SQRT(2.*ALOG(2))

      FOR jj=0L,long(nl)-1 DO BEGIN

;include here ONLY the lines in index:

         i = index(jj)

         IF keyword_set(kev)  THEN begin 
            wavelength=ang2kev/wvl[i] 
            x_fwhm=ang2kev/INSTR_FWHM
         endif   else begin 
            wavelength=wvl[i]
            x_fwhm=INSTR_FWHM
         end 
         tm = 10.^(tmax[i])
         pp = get_atomic_weight(iz[i])

; fwhm in Angstroms:
         IF keyword_set(no_thermal_width) THEN BEGIN 
           fwhm=x_fwhm
         ENDIF ELSE BEGIN 
           fwhm = SQRT( 5.093d-13 * tm * wavelength^2 / $
                        pp + X_FWHM^2)
         ENDELSE 

         IF keyword_set(kev)  THEN fwhm=ang2kev/fwhm

         ind = WHERE( ABS(lambda-wvl[i]) LT 2.5*fwhm )

;; The if statements below check if the line extends over more than one 
;; pixel. If it doesn't, then the way the line is added to spectrum is 
;; different.


         IF (ind[0] EQ -1) OR (N_ELEMENTS(ind) EQ 1) THEN BEGIN

            getmin = MIN(ABS(lambda-wvl[i]),ind)
            spectrum[ind] = spectrum[ind] + $
              intensities[i]/angpix[ind]

            peak[jj]= intensities[i] /angpix[ind]

         ENDIF ELSE BEGIN

            peak[jj]= intensities[i] / fwhm / 1.0645

            int_lambda = [lambda[ind[0]] - (binsize[ind[0]])/2.,lambda[ind]+binsize[ind]/2.]

            int_lambda = (int_lambda-wvl[i])* SQRT(8.*ALOG(2)) / fwhm
            nind = N_ELEMENTS(ind)+1
            ints = DBLARR(nind)
            amplitude = intensities[i]/fwhm/cc
            FOR j=0,nind-1 DO BEGIN
               ints[j] = amplitude*cc*fwhm*GAUSSINT(int_lambda[j])
            ENDFOR 
            line_profile = ints[1:nind-1]-ints[0:nind-2]
            spectrum[ind] = spectrum[ind] + line_profile/angpix[ind]
         ENDELSE

      ENDFOR 

   END 
ENDIF   else     print, '% MAKE_CHIANTI_SPEC: no lines ??! '

;
; Add the continuum. Note that the temperature intervals are assumed to 
; be 0.1 dex. This requires access to the DEM and ABUND specifications in the
; structure.

;
IF KEYWORD_SET(continuum) THEN BEGIN

   ioneq_name=transitions.ioneq_name
   IF file_exist(ioneq_name) THEN $
     read_ioneq, ioneq_name , ioneq_logt,ioneq,ioneq_ref ELSE BEGIN 
      err_msg = 'error, could not find Ion. Fraction file: '+ioneq_name
      return
   END 

;   print, 'Calculating continuum .... '

; *** PRY, 6-Oct-2010: I've commented out the lines below as
; they caused the routine to crash. They're now placed inside
; the following IF statement
;
;;    CASE transitions.model_name OF
;;       'Constant pressure': BEGIN 
;;          edensity = transitions.model_pe/temp
;;       END 
;;       'Constant density': edensity = transitions.model_ne
;;       'Function': BEGIN 

;;          edensity =10.^SPLINE(alog10(transitions.model_te),$
;;                               alog10(transitions.model_ne), alog10(TEMP) )

;;       END 

;;    ENDCASE

   IF tag_exist(transitions,'DEM')  THEN BEGIN 

      dem_name = transitions.dem_name
      dem_logt = transitions.dem_logt
      dem = transitions.dem

; convert the DEM distribution to match the temperatures in ioneq_logt
      IF transitions.model_name EQ 'Function' THEN BEGIN 

         gdt=WHERE((ioneq_logt GE MIN(dem_logt)) AND $
                   (ioneq_logt LE MAX(dem_logt)) AND $
                   (ioneq_logt GE min(alog10(transitions.model_te))) AND $
                   (ioneq_logt LE  max(alog10(transitions.model_te))) , nn)

      ENDIF ELSE BEGIN 

         gdt=WHERE((ioneq_logt GE MIN(dem_logt)) AND $
                   (ioneq_logt LE MAX(dem_logt)) , nn) 
      END

      dem_int1=10.^(SPLINE(dem_logt,dem,ioneq_logt[gdt] ))

      dem_int=FLTARR(n_elements(ioneq_logt))
      ngt=n_elements(gdt)
      FOR igt=0,ngt-1 DO  dem_int[gdt[igt]]=dem_int1[igt] >  0.

      dem_int =  dem_int[gdt]
      temp = 10.^ioneq_logt[gdt]

     ;
     ; Compute electron density (edensity) for 2 photon calculation
     ;
      CASE transitions.model_name OF
        'Constant pressure': BEGIN 
          edensity = transitions.model_pe/temp
        END 
        'Constant density': edensity = transitions.model_ne
        'Function': BEGIN 
          edensity =10.^SPLINE(alog10(transitions.model_te),$
                              alog10(transitions.model_ne), alog10(TEMP) )
        END 
      ENDCASE


; ergs cm^3 s-1 str^-1 Ang^-1 per unit emission measure

    ; Note that dem_int= is used to specify the DEM values for the
    ; continuum routines (PRY, 7-Oct-2009)
    ;
      freebound, temp,lambda,fb,/no_setup,min_abund=min_abund, $
                 photons=photons, dem_int=dem_int,/sumt, VERBOSE=VERBOSE, kev=kev, $
               ioneq_file=ioneq_name, abund_file=abund_name

      freefree,temp,lambda,ff,min_abund=min_abund, $
               photons=photons, dem_int=dem_int,/sumt, VERBOSE=VERBOSE, kev=kev, $
               ioneq_file=ioneq_name, abund_file=abund_name

      two_photon, temp,  lambda, two_phot,min_abund=min_abund, $
                  edensity=edensity, photons=photons, dem_int=dem_int,/sumt, $
                  VERBOSE=VERBOSE, kev=kev, abund_file=abund_name, $
                  ioneq_file=ioneq_name, lookup=transitions.lookup

   ENDIF   ELSE BEGIN 

      em_int=10.^transitions.logem_isothermal
      temp = 10.^transitions.logt_isothermal      

     ;
     ; Compute electron density (edensity) for 2 photon calculation
     ;
      CASE transitions.model_name OF
        'Constant pressure': BEGIN 
          edensity = transitions.model_pe/temp
        END 
        'Constant density': edensity = transitions.model_ne
        'Function': BEGIN 
          edensity =10.^SPLINE(alog10(transitions.model_te),$
                              alog10(transitions.model_ne), alog10(TEMP) )
        END 
      ENDCASE

    ; Note that em_int= is used to specify the EM values for the
    ; continuum routines(PRY, 7-Oct-2009)
    ;
      freebound, temp,lambda,fb,/no_setup,min_abund=min_abund, $
                 photons=photons, em_int=em_int,/sumt, VERBOSE=VERBOSE, kev=kev, $
               ioneq_file=ioneq_name, abund_file=abund_name

      freefree,temp,lambda,ff,min_abund=min_abund, $
               photons=photons, em_int=em_int,/sumt, VERBOSE=VERBOSE, kev=kev, $
               ioneq_file=ioneq_name, abund_file=abund_name

      two_photon, temp,  lambda, two_phot,min_abund=min_abund, $
                  edensity=edensity, photons=photons, em_int=em_int,/sumt, $
                  VERBOSE=VERBOSE, kev=kev, ioneq_file=ioneq_name, $
                  abund_file=abund_name, lookup=transitions.lookup

   ENDELSE 


;freefree  freebound and two_photon give as output units/Angstroms.
   continuum=1.d-40*(fb+ff+two_phot)


   IF n_elements(effarea) GT 0 THEN BEGIN 

; from photons cm^-2 s-1 str^-1 Ang^-1 to  photons cm^-2 s-1 bin-1 :
      continuum = continuum*angpix ;*4.0*!pi 

; to counts s-1 bin-1:
      continuum = continuum*effarea
   ENDIF 

   spectrum = spectrum+continuum

ENDIF   



IF n_elements(effarea) GT 0 THEN  BEGIN 
   int_units = 'counts s-1' 
   units(1) = 'counts s-1 bin-1' 
ENDIF ELSE  units(1) = int_units+' '+units(0)+'-1'


OUTPUT = TRANSITIONS 
lines=OUTPUT.lines

IF nl NE  0  THEN BEGIN 

;include only the lines effectively used to build the spectrum:
;**************************************************************

;over-write the structure:

;    IF keyword_set(kev)  THEN  lines=lines[sort(wvl)]
   lines=lines[index]

;add a tag 'peak' to each structure: 

   str = lines[0] 
   str = add_tag(str, 0.0,'peak')

   new_lines = REPLICATE(str, n_elements(lines))

   IF  n_elements(tag_names(new_lines)) NE 12 THEN BEGIN 
      err_msg = '% MAKE_CHIANTI_SPEC: error, line intensities structure not compatible '
      print, err_msg
      return
   END 

   new_lines[*].peak = peak

; these are wavelengths (A) or energy units (keV):

   new_lines[*].wvl = wvl[index]
   new_lines[*].int = intensities[index]

;??
   OUTPUT.wvl_limits = [min(new_lines.wvl), max(new_lines.wvl)]

   new_lines[*].iz = lines[*].iz
   new_lines[*].ion = lines[*].ion
   new_lines[*].ident = lines[*].ident
   new_lines[*].ident_latex = lines[*].ident_latex
   new_lines[*].snote = lines[*].snote
   new_lines[*].flag = lines[*].flag
   new_lines[*].lvl1 = lines[*].lvl1
   new_lines[*].lvl2 = lines[*].lvl2
   new_lines[*].tmax = lines[*].tmax


;lines = add_tag(lines,0.0,'peak')

endif else begin 
; no lines 

   IF NOT  KEYWORD_SET(continuum) THEN  BEGIN 
      err_msg = '% MAKE_CHIANTI_SPEC: error, no lines, no continuum in the requested range! '
      print, err_msg
      return
   END 

; - just add one line 

   new_lines=lines[0]
   new_lines = add_tag(new_lines, 0.0,'peak')
   new_lines.int=0.
   new_lines.wvl=0.
   new_lines.iz = 0
   new_lines.ion = 0
   new_lines.ident = ''
   new_lines.ident_latex = ''
   new_lines.snote = ''
   new_lines.flag = 0
   new_lines.lvl1 = 0
   new_lines.lvl2 = 0
   new_lines.tmax = 0

end 

OUTPUT = rem_tag(OUTPUT, 'lines')
OUTPUT = add_tag(OUTPUT, new_lines, 'lines')


;intensity units: 'counts s-1' if effarea is used, otherwise retain the original.
OUTPUT.int_units = int_units    


OUTPUT = CREATE_STRUCT(OUTPUT,'LAMBDA', lambda )
OUTPUT = CREATE_STRUCT(OUTPUT,'SPECTRUM',spectrum )

OUTPUT = CREATE_STRUCT(OUTPUT,'UNITS',units ) ;spectrum units  
OUTPUT.WVL_UNITS = units[0]

OUTPUT = CREATE_STRUCT(OUTPUT,'INSTR_FWHM',INSTR_FWHM)
OUTPUT = CREATE_STRUCT(OUTPUT,'BIN_SIZE',bin_size )

OUTPUT = CREATE_STRUCT(OUTPUT,'ABUND_NAME',abund_name )
OUTPUT = CREATE_STRUCT(OUTPUT,'ABUND',abund )
OUTPUT = CREATE_STRUCT(OUTPUT,'MIN_ABUND',min_abund )
OUTPUT = CREATE_STRUCT(OUTPUT,'ABUND_REF',abund_ref )

IF KEYWORD_SET(continuum) THEN OUTPUT = CREATE_STRUCT(OUTPUT, 'CONTINUUM', continuum) 
IF n_elements(effarea) GT 0 THEN BEGIN 
   OUTPUT = CREATE_STRUCT(OUTPUT,'FILE_EFFAREA', file_effarea)
   OUTPUT = CREATE_STRUCT(OUTPUT,'EFFAREA', effarea)
END 

IF keyword_set(verbose) THEN $
  print,format='("% MAKE_CHIANTI_SPEC: spectrum calculated in ",f8.1," seconds")',$
  systime(1)-t0


END
