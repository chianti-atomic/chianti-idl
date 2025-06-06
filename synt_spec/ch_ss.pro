;+
; PROJECT
;
;        Part of this program was developed within CHIANTI-VIP. 
;        CHIANTI-VIP (Virtual IDL and Python) is  a member of the CHIANTI family
;        mantained by Giulio Del Zanna, to develop additional features and
;        provide them to the astrophysics community.
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving  the
;       University of Cambridge,  Goddard Space Flight Center, and University of Michigan. 
;
; NAME
;
;     CH_SS
;
; PURPOSE:
;       Widget-based multi-purpose routine 
;       to calculate CHIANTI line intensities and continua, to create a 
;       synthetic spectrum, to make tables of lines, etc. 
;
; CALLING SEQUENCE:
;
;       IDL> ch_ss
;
;
; PROCEDURE:
;
;     This routine calculates a synthetic spectrum by merging line 
;     intensities and continua. 
;
;
;     The widget is organised into four Sections: 
;
;     SECTION 1:
;
;      -The Calculation of the  CHIANTI line intensities. 
;      
;       This can be  done in two ways:
;
;       1-Restore a save file with the CHIANTI line intensities already
;       calculated.   
;
;       2-Calculate CHIANTI line intensities  with a call to CH_SYNTHETIC. 
;
;       In this case, A series of parameters must be set:
;
;       - Minimum and maximum wavelengths in Angstroms
;
;       - The model used for the calculation. Three are the options:
;          1) a constant  density (cm^-3) 
;          2) a constant pressure (cm^-3 K) 
;          3) a general (Te,Ne) model. In this case, a file will be read.
;             This file should have two columns, one with the Te (K)
;             values, and one with the Ne (cm^-3) values.
;
;       - The  ionization fraction file to be used.  "*.ioneq"  files
;          can be selected from  either  the CHIANTI database, the 
;          working directory or selected via a widget.
;
;       - All ions ?  If set to yes (default), then all the ions present in the
;                     database will be included.
; 
;                     If set to no, then it is possible to select a list of ions
;                     with a widget 
;
;       - All lines ? If set to no (default), only the lines for which there are
;                      observed energy levels are included 
;       
;                     If set to yes, also the lines that do not have
;                     corresponding observed energy levels are included. In this
;                     case, the wavelengths are calculated from the theoretical
;                     energy levels, and might not be very accurate.
;       
;       - Isothermal ?  If set to no (default), a DEM file must be selected. 
;                         "*.dem"  files (i.e. files with a .dem extension)
;                       can be selected from either  the CHIANTI database, the 
;                       working directory or selected via a widget.
;       
;                       If set to yes, then the user is requested to enter one
;                       or more temperatures (as logarithmic values - Log T )
;                       and correspondent column emission measures EM
;                       logarithmic values.
;                       NOTE: if more than one value is entered, then the
;                             sequence must be  separated by commas (e.g.: 6.0,
;                             6.5, 7.), and both Log T and Log EM must have the
;                             same number of values 
;
;       - Photoexcitation ?
;                       If set to yes, you have to define:
;                       Trad: The blackbody radiation field temperature
;                       R/Ro: Distance from the centre of the star in stellar
;                             radius units 
;
;       Units:  Photons or Ergs'
;       
;       Protons: If set to Yes, the proton data are used to calculate the level population    
;
;
;     Once all the parameters have been defined, the user should click on the
;     "Calculate intensities" button to start the calculation (which calls 
;     CH_SYNTHETIC). 
;
;     Once the calculation is finished, an IDL  structure is loaded into
;     memory. It is then possible to save it for later use by clicking  on the
;     "SAVE" button. 
;      The RESTORE button is to restore previously saved files into an IDL
;      structure in memory.
;
;     Once the IDL structure with the line intensities is in the memory, it is
;     then possible to calculate and plot a spectrum (SECTION 2).
;
;     SECTION 2:
;
;              This section controls the parameters that are needed to fold the
;              line intensities and the continua into a synthetic
;              spectrum. These parameters are used by MAKE_CHIANTI_SPEC.
;
;              Before this is done, a set of line intensities MUST be in the
;              program memory. This is done either by calculating the
;              intensities  or by restoring a save file with
;              previously calculated values (SECTION 1). 
;
;              Setting the parameters:
;
;              -Minimum and maximum wavelengths in Angstroms
;
;              -spectrum bin size in Angstroms. Disallowed if an Effective area
;                 file is used.
;              
;              -instrumental FWHM: Setting this to a non-zero value broadens
;                                  each of the spectral lines with a Gaussian of
;                                  the specified FWHM (in Angstroms) so
;                                  mimicking the effects of instrumental
;                                  broadening. 
;              
;              -continuum: Add continua to the binned spectrum:
;                          free-free, free-bound and two-photon.
;                          Please note that the continuum calculation takes some
;                          time and you may want to define a minimum abundance
;                          value to speed the calculations.
;
;              - All lines ? If set to no (default), only the lines for which there are
;                           observed energy levels are included.
;                           If set to yes, the "unobserved lines" will be added, but
;                           only if they are present in the structure.
;
;
;              -elemental abundances
;                          "*.abund"  files (i.e. files with a .abund
;                         extension) can be selected either from the CHIANTI database,
;                          the  working directory, or via a widget. 
;
;              -select a minimum abundance value
;                         If set not null, only the lines of those elements
;                         which have an abundance greater than the value set are
;                         selected. Also, the continuum is calculated only for
;                         those elements which  have an abundance greater than
;                         the value set. This can significantly speed up the
;                         calculations. By default, the minimum value in the
;                         selected abundance file is used. To have an idea of
;                         what minimum abundance should be set, the abundances
;                         of Allen (1973) give: 
;
;                         abundance (H)  = 1.
;                         abundance (He) = 0.085
;                         abundance (C)  = 3.3e-4
;                         abundance (O)  = 6.6e-4
;                         abundance (Si) = 3.3e-5
;                         abundance (Fe) = 3.9e-5
;
;              
;              Eff. Area: Yes/No 
;
;              If you want to fold the spectrum with an effective area.  
;              If set to Yes, you are requested to choose an input ascii file
;              with two columns, the wavelength and the effective area values
;              (cm^2). 
;              The wavelenghts in the file (that might not be linear)  are used
;              to create the spectrum, that is multiplied with the effective
;              area values.
;		    Note that this option only works well if a sufficient number
;		    of bins is given. The line intensities contributing to each
;		    bin are summed, and  subsequently convolved with a gaussian
;		    of full-width-half-maximum FWHM, if FWHM is not set = 0.
;                   Please note that the convolution might not work if a small
;                   number of  bins is defined. 
;
;                Also note that to have the correct output units  (counts s-1
;                bin-1) the appropiately scaled DEM (or EM) values must be provided.
;
;
;              After this, by clicking on the "Calculate and plot" button the
;               program calculates and plots the synthetic spectrum. 
;              
;              Once the spectrum is displayed, it is then possible to
;              view the details of the lines by clicking with the mouse in the
;              plot window, and to  perform various operations by clicking on
;              the buttons in SECTION 3
;              
;     SECTION 3:
;
;              This Section allows the user to select a few parameters for the
;              plotting, and to create different types of OUTPUT.
;              
;              Labels ? : Setting this to yes plots a vertical line for each
;                         spectral line in the spectrum, and also writes a label
;                         above the strongest lines indicating the ion from
;                         which the line arises. 
;              
;              Min.:      Only lines which have an intensity greater than  the
;                         value set here will be listed and, if requested,
;                         labelled and selected for inclusion in the various
;                         outputs.  Setting the value=0.  will result in all
;                         lines being listed and written in the outputs.
;              
;              X,Y, XOOM, UNZOOM: It si possible to select a region of the
;                                 spectrum, by zooming with the use of the mouse
;                                 or by setting the X,Y ranges.
;
;                                NOTE that only the line details and portion of
;                                the spectrum shown will be output.
;
;              LINEAR/LOG  To plot the spectrum in linear or log scale
; 
;              Create PS file: A postscript file is created. 
;
;              Hardcopy: the postscript file "idl.ps" is created and sent to the
;                        default printer. 
;              
;              Save Line details (latex): The  details of the lines shown in the
;                                         plot will be  saved in a latex file.
;              
;              Save Line details (ascii): The  details of the lines shown in the
;                                         plot will be  saved in an ascii file.
;              
;              
;
;              Save Spectrum (ascii): The  X,Y values of the plot are saved in
;                                     an ascii file.
;
;              Save Spectrum (IDL/FITS): The details of all the lines and the arrays
;                                   of the X,Y values of the plot are saved into
;                                   an IDL or FITS file. The  IDL structure 
;                                   has the following tags: 
;
;            .LAMBDA:   The array of wavelength X values
;
;            .SPECTRUM: The array of spectrum Y values
;
;            .UNITS       The units of LAMBDA, SPECTRUM 
;            .INSTR_FWHM  The Instrumental FWHM
;            .BIN_SIZE         Width of the Bins  (fixed) in angstroms
;            .ABUND_NAME  The CHIANTI abundance file name
;            .ABUND       The abundance values
;            .MIN_ABUND   The minimum abundance value used
;            .ABUND_REF   The references
;            .CONTINUUM   The values of the continuum (if
;                                                calculated)
;            .EFFAREA       The array of effective area
;                                      values (optional)
;            .FILE_EFFAREA  The name of the effective area file used (optional).
;
;
;            .IONEQ_NAME     The ion balance file used (full path).
;            .IONEQ_LOGT        The Log10 T values associated.
;            .IONEQ_REF      The references.
;
;            .DEM_NAME       The differential emission measure file eventually  used
;                            (full path).
;            .DEM            The Log10 DEM values 
;            .DEM_LOGT          The Log10 T values associated.
;            .DEM_REF        The references.
;
;            .MODEL_NAME    A string indicating the model used 
;                     (e.g. constant density or constant pressure).
;
;            .MODEL_NE    the Ne value.
;            .MODEL_PE    the Pe value.
;
;
;            .WVL_UNITS  The wavelength units.
;
;            .WVL_LIMITS    The wavelength limits specified by the user.
;
;            .INT_UNITS  The intensity units
;
;            .LOGT_ISOTHERMAL
;                       The Log10(T) values used. 
;
;            .LOGEM_ISOTHERMAL
;                       The Log10(EM) values used. 
;
;            .TIME      The date and time when the structure was created.
;
;            .VERSION   The version number of the CHIANTI database used.
;
;            .ADD_PROTONS 
;                       A flag (0/1) to indicate whether proton data were used (1)
;                       or not (0) to calculate the level population.
;
;            .PHOTOEXCITATION
;                       A flag (0/1) to indicate if photoexcitation was included (1)
;                       or not (0).
;
;            .RADTEMP 
;                      The blackbody radiation field temperature used (if
;                      photoexcitation was included).
;
;            .RPHOT
;                   Distance from the centre of the star in stellar radius units  
;                   (if photoexcitation was included).
;
;                                       THEN, FOR EACH LINE USED TO CALCULATE THE
;                                       SPECTRUM:
;
;            .LINES     A structure containing information about the lines. 
;                       Its size is the number of lines in the spectrum. The 
;                       tags are:
;
;                  .peak   The peak intensity value 
;
;                  .iz     The atomic number of the elements (e.g., 26=Fe)
;
;                  .ion    The ionisation stage (e.g., 13=XIII)
;
;                  .snote  The identification of the ion (e.g., 'Fe XXIV d')
;
;                  .ident  The identification of the transition, configuration
;                           and terms in text form.
;
;                  .ident_latex
;                          The identification of the transition, configuration
;                           and terms in latex form.
;
;                  .lvl1   The lower level of the transition (see .elvlc 
;                          file for ion)
;
;                  .lvl2   The upper level for transition.
;
;                  .tmax   The temperature of maximum emission of the line 
;                          (i.e., the temperature at which the product of 
;                          the emissivity and the ion fraction has its 
;                          maximum). Rounded to nearest 0.1, and zero in case
;                          the isothermal approximation is used. 
;
;                  .fwhm 
;
;                  .wvl    Wavelength of the transition, in Angstroms.
;
;                  .flag   A flag, =-1 if the line has only theoretical energy
;                          levels. Otherwise flag=0.
;
;                  .int    Intensity of line  (with the  abundance factor multiplied)
;
;
;
;              Save Spectrum (FITS): The entire information contained in the
;                                    IDL structure is stored in a FITS file. 
;                                    
;     SECTION 4:
;              Here, text information messages are printed. 
;
; INPUTS
;
;     None.
;
; OPTIONAL INPUTS:
;     The font 
;
; OUTPUTS:
;     Many.
; 
; KEYWORD PARAMETERS:
;     
;     FONT  the font to be used. Can be useful to customize the appearance of
;           the widget.
;
; CALLS:
;
;     CH_SYNTHETIC, CH_LINE_LIST, CH_DRAWBOX, MAKE_CHIANTI_SPEC, CH_XMENU_SEL,
;     plus many other CHIANTI and SolarSoft routines.
;
;
; PROGRAMMING NOTES
;
;     Within CH_SS, there are several other routines which are:
;
;     OPLOT_LINES   This overplots lines and a label on the displayed 
;                   spectrum.
;     SYN_CURSOR    When the mouse is clicked when on the spectrum window, 
;                   this routine prints out the list of nearby lines and 
;                   their IDs in the text window.
;
;     CALC_SYN_SPECTRUM
;                   Calculates line intensities with a call to CH_SYNTHETIC
;
;     PLOT_SYN_SPECTRUM  This calls make_chianti_spec to produce the 
;                        intensity vs. wavelength plot.
;     SYN_MAIN_EVENT  This handles the widget operations
;
;     SYN_WID       This creates the widgets.     
;
; COMMON BLOCKS:
;        many
;
; RESTRICTIONS:
;
; SIDE EFFECTS:
;
; EXAMPLE:
;
;     IDL> ch_ss
;
; CATEGORY:
;	
;	spectral synthesis.
;
;
; WRITTEN     : 
;
;       Ver.1, 7-Nov-01, Giulio Del Zanna (GDZ) and Peter Young (PRY)  
;
; MODIFICATION HISTORY:
;
;       V.2, 7-Nov-01, GDZ . Fixed a small bug (now the spectrum plot is always
;       plotted within the widget), and modified the option to add continua.
;       Changed the suggested  names of the outputs. 
;       Corrected a bug when creating an IDL save file with the spectrum, when
;       no line details are present.
;
;       V.3 28-Jan-02 GDZ
;           fixed a bug in the density text widget, added a few buttons 
;           and options, including the effective area.
;            Added noprot, rphot, radtemp keywords to the call to ch_synthetic
;
;       V 4, 18-Apr-2002, GDZ 
;           Added  photoexcitation, changed IDL save files to FITS files,
;
;       V.5, 21-May-2002, GDZ 
;        fixed a few small bugs: checking min_abund before calculating the
;                                spectrum; checking the ioneq file when
;                                restoring the structure; changed the status of
;                                all lines;  chnaged the font system.
;                               generalized directory concatenation to work for 
;                               Unix, Windows  and VMS.
;
;       V.6, 15-July-2002, GDZ - New major revision.
;
;           Changed the chianti top directory (for Effective areas).
;           Changed  Linear/Log button. 
;
;            Rearranged the sizes of the buttons and added a special cursor to
;            highlight the area where  details of the lines will be given. Works
;            only in linear scale.
;
;           Added quite a lot of new checks to avoid crashes and
;           fixed the problem with the zoom/unzoom/change units.
;
;       V.7, 2-Aug-02, GDZ
;
;           Modified the output labels on the plot, inside and on the axis.
;           Also modified a few minor things like the appearance of the Log T,EM
;           values. 
;           Fixed a bug when creating the latex output.
;           Now it restores at the end  previous colors and settings. 
; 
;       V.8, 8-Aug-02, GDZ
;           Changed the CHIANTI system variables. Fixed.
;           Also fixed a problem with the element ab. file.
;
;      V.9, 13-Aug-02, GDZ 
;
;           Restored the correct use of ch_line_int, now only the lines in the
;           plot window are listed, and the ALL keyword is in use. 
;           Now the correct xrange is loaded into COMMON when line int. are
;           restored. Now it checks if all ions were in the structure, when
;           restoring the line intensities, and flags the widget button accordingly.
;           Added a device,decomposed=0. to remove problems with colors.
;           Corrected the use of the DEM, IONEQ and ABUND pulldown menus,
;           avoiding conflicts between files in the working and CHIANTI
;           directory having the same name.
;           Added printing of references for ancillary files, and a check on the
;           element abundances vs. the elements present in the structure.
;
;       V.10, 7-Nov-03  GDZ
;
;          Modified format e8.2 to e9.2 for Windows compatibility.
;          Replaced f9.4 with f11.4 format for the wavelengths.
;          Some minor modifications to the widget.
;
;          Added extended details in the ascii output spectrum.
;
;          Added more explanations in the HELP buttons.
;
;       V.11,  22-Jul-2005  GDZ 
;
;          -Added keV option and a few more extra checks.
;
;       V.12,  2-Aug-2005 GDZ
;           put RETAIN=2 in the main plotting window.
;
;       V.13, 3-Oct-2005 GDZ
;          Replaced FOR i=0, calls with FOR i=0L, calls, so 
;          the routine does not crash with a large number of lines.
;
;       V.14, 26-Jun-2012, Peter Young
;          Now allow the case of FWHM=0 (this was blocked
;          previously). Note that this means the lines will broadened
;          by their thermal width (determined from Tmax value).
;
;       V.15, 9-Jul-2012, Peter Young
;          Converted all list() calls to list[] for compatibility with
;          IDL 8; modified help message for FWHM box.
;
;       V.16, 26-Jul-2012, Peter Young
;          I've swapped all the 'bigpickfile' calls for
;          'dialog_pickfile' in an attempt to fix a bug reported by
;          Ryan Milligan.
;
;       V.17, 17-Aug-2015, Peter Young
;          Replaced 'lambda()' with 'lambda[]' to prevent crashes due
;          to new IDL lambda.pro routine (IDL 8.4).
;
;       V.18, 11-Jun-2020, Peter Young
;          Added widget for using population lookup tables; now
;          automatically choose the chianti.ioneq file; removed some
;          of the zeros from the wavelength range text boxes.
;       V.19, 25-Sep-2020, Peter Young
;          Fixed bug when $CHIANTI_LOOKUP not defined.
;       V.20, 19-Nov-2020, Peter Young
;          Modified behavior when switching to log plots to better
;          display the spectrum.
;       V.21, 11-Dec-2020, Peter Young
;          Changed dem_int to double-precision.
;       V.22, 23-Jun-2023, Peter Young
;          The Tmax given in the text widget is now given to two decimal places.
;       v.23, 18-Jan-2024, Peter Young
;          Fixed bug when reading abundance files in the new archive directory.
;
;       v.24, 10 May 2024, Giulio Del Zanna
;          added the level indices of the transitions.
;       v.25 13 June 2024  Giulio Del Zanna
;          Major rewrite for version 11.
;       v.26 19 July 2024 GDZ, fixed the bug when calculating the isothermal option.
;       v.27, 23 Sept 2024 GDZ, fixed a bug: the CT option was not passed on.
;       v.28, 16-Apr-2025, Peter Young
;          If the lookup tables are enabled, then the photoexcitation options are
;          now disabled (since the two are not compatible).
;
;
; TO DO LIST:
;           Control the range of Angstroms when clicking
;           kev
;           Allow plots in intensities instead of intensities A-1
;
; VERSION     :  V.28
;
;-
PRO restore_spectrum

COMMON line_data, tran, state
COMMON plot_spec,kev_yn, ang, inst, o_lines, xrange, yrange, sty, theor_lines, $
  cont_ind, lambda, effarea, file_effarea, spectrum, fold_yn, log
COMMON abundance, abfile, abdir, abund, abund_ref,min_abund 
COMMON calc_int, int_xrange, const_names, const_nt, const_value,$
  ioneqfile, ioneqdir, demfile, demdir, isothermal,iso_logem, isothermal_flag, $
  all_ions_yn, list_ions, units, noprot,photoexcitation, rphot, radtemp

ang = spectrum.BIN_SIZE
WIDGET_CONTROL,state.ang_read,set_value=trim(ang)

inst = spectrum.INSTR_FWHM
WIDGET_CONTROL,state.inst_read,set_value=trim(inst)

o_lines = 1
WIDGET_CONTROL,state.oplot_lines,set_value=0

log = 0
; spectrum has been 'zoomed' (1) or not (0)
sty = 0

;reset X,Y ranges:
xrange = [min(spectrum.lambda), max(spectrum.lambda)]
yrange = [0.,max(spectrum.spectrum)*1.2]

WIDGET_CONTROL,state.xmin_base,  set_value=strtrim(string(format='(f11.4)',xrange(0)),2)
WIDGET_CONTROL,state.xmax_base,  set_value=strtrim(string(format='(f11.4)',xrange(1)),2)
WIDGET_CONTROL,state.ymin_base,  set_value=strtrim(string(format='(e9.2)',yrange(0)),2)
WIDGET_CONTROL,state.ymax_base,  set_value=strtrim(string(format='(e9.2)',yrange(1)),2)

index1 = where(spectrum.lines[*].flag  EQ -1, nlines1)
IF nlines1 EQ 0 THEN theor_lines = 0 ELSE theor_lines =1
WIDGET_CONTROL,state.rem_theor2,  set_value=theor_lines

; 
cont_ind=  tag_exist(spectrum, 'CONTINUUM') 
WIDGET_CONTROL,state.continua,set_value=cont_ind

lambda = spectrum.lambda


IF tag_exist(spectrum, 'FILE_EFFAREA' ) THEN BEGIN 

   widget_control,state.fold_info, set_val='Eff. Area: YES'
   fold_yn = 1
   WIDGET_CONTROL,state.ang_read,set_value=' -- '

   IF file_exist(file_effarea) THEN BEGIN 
      break_file,file_effarea, disk, dir, f, ext
      f = f+ext
      widget_control,state.fold_show, set_val=f+ext
   ENDIF ELSE BEGIN 
      widget_control,state.fold_show, set_val=' -- '
      result = DIALOG_MESSAGE(['Error, could not find  Eff. Area file: ', $
                               file_effarea] ,/info)
   END 

ENDIF ELSE BEGIN 
   widget_control,state.fold_info, set_val='Eff. Area: NO'
   fold_yn = 0
   widget_control,state.fold_show, set_val=' -- '
END   


abund_name = spectrum.abund_name

IF  file_exist(abund_name) THEN BEGIN 
   break_file, abund_name, disk, abdir, abfile, ext
   abdir = disk+abdir

   abfile = abfile+ext
   WIDGET_CONTROL,state.ab_show,set_value=abfile
ENDIF ELSE BEGIN 
   WIDGET_CONTROL,state.ab_show,set_value=' -- '
   result = DIALOG_MESSAGE(['Error, could not find abundance   file: ', $
                            abund_name] ,/info)
END 

abund = spectrum.abund 
abund_ref = spectrum.abund_ref
min_abund = spectrum.min_abund

WIDGET_CONTROL, state.min_abund_ev,$
  set_value=strtrim(string(MIN_ABUND, format='(e9.2)'),2)


widget_control, state.show_lines , $
  set_val='spectrum restored and parameters reset !'


END 

PRO restore_line_int

COMMON line_data, tran, state
COMMON calc_int, int_xrange, const_names, const_nt, const_value,$
  ioneqfile, ioneqdir, demfile, demdir, isothermal,iso_logem, isothermal_flag, $
  all_ions_yn, list_ions, units, noprot,photoexcitation, rphot, radtemp
COMMON plot_spec,kev_yn, ang, inst, o_lines, xrange, yrange, sty, theor_lines, $
  cont_ind, lambda, effarea, file_effarea, spectrum, fold_yn, log


;now  a lot of resets are needed:
;reset X,Y ranges:

int_xrange = tran.wvl_limits
xrange = int_xrange

units(0) = tran.WVL_UNITS

IF tran.WVL_UNITS  EQ 'keV' THEN BEGIN 
   widget_control,state.kev_info, set_val='keV [ / '+string(197B)+' ]'
   kev_yn = 1
   int_xrange = reverse(12.39854/int_xrange)
ENDIF ELSE IF tran.WVL_UNITS EQ 'Angstroms' THEN  BEGIN 
   widget_control,state.kev_info, set_val=string(197B)+' [ / keV ]'
   kev_yn = 0
END 


yrange = [0.,0.]

WIDGET_CONTROL,state.wminw1,  set_value=strtrim(string(format='(f11.4)',int_xrange(0)),2)
WIDGET_CONTROL,state.wmaxw1,  set_value=strtrim(string(format='(f11.4)',int_xrange(1)),2)


WIDGET_CONTROL,state.wminw2,  set_value=strtrim(string(format='(f13.6)',xrange(0)),2)
WIDGET_CONTROL,state.wmaxw2,  set_value=strtrim(string(format='(f13.6)',xrange(1)),2)

WIDGET_CONTROL,state.xmin_base,  set_value=strtrim(string(format='(f11.4)',xrange(0)),2)
WIDGET_CONTROL,state.xmax_base,  set_value=strtrim(string(format='(f11.4)',xrange(1)),2)
WIDGET_CONTROL,state.ymin_base,  set_value=strtrim(string(format='(e9.2)',yrange(0)),2)
WIDGET_CONTROL,state.ymax_base,  set_value=strtrim(string(format='(e9.2)',yrange(1)),2)


WIDGET_CONTROL, state.plot_rat, GET_VALUE=plot_rat_id,  sensitive=1 
wset,plot_rat_id
tverase


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
CASE  tran.model_name OF 
   'Constant density': BEGIN 
      const_nt = 0
      WIDGET_CONTROL,state.const_widg,set_value=trim(const_nt)
      const_value=tran.model_NE
      WIDGET_CONTROL,state.const_read, sensitive=1, $
        set_value=strtrim(string(format='(e9.2)',const_value),2)
   END

   'Constant pressure':BEGIN 
      const_nt =1
      WIDGET_CONTROL,state.const_widg,set_value=trim(const_nt)
      const_value=tran.model_pe
      WIDGET_CONTROL,state.const_read, sensitive=1, $
        set_value=strtrim(string(format='(e9.2)',const_value),2)
   END
   'Function':BEGIN 
      const_nt =2

      delvarx, const_value

      WIDGET_CONTROL,state.const_widg,set_value=trim(const_nt)

      IF tran.model_file NE ' ' THEN BEGIN 

         IF file_exist(tran.model_file) THEN BEGIN 
            break_file, tran.model_file, disk, dir, file,ext
            const_value = file+ext
            WIDGET_CONTROL,state.const_read, sensitive=0, set_value=const_value

         ENDIF ELSE BEGIN 
            WIDGET_CONTROL,state.const_read, sensitive=1, set_value=' - '
            result = DIALOG_MESSAGE(['Error, could not find (Te,Ne) file: ', $
                                     tran.model_file] ,/info)
         END 
      ENDIF ELSE BEGIN 
         WIDGET_CONTROL,state.const_read, sensitive=1,  set_value=' - '
         result = DIALOG_MESSAGE(['Error, could not find (Te,Ne) file: ', $
                                  tran.model_file ],/info)
      END 

   END

;else do nothing.
ENDCASE 


;all ions ?

list_ions = strarr(n_elements(tran.lines))

FOR i=0L, long(n_elements(tran.lines)-1) DO BEGIN 
   spectroscopic2ion, tran.lines[i].snote, ion
   list_ions[i]=ion
END 

;remove multiple occurrencies

list_ions = list_ions(rem_dup(list_ions))

;check current masterlist:

mname=concat_dir(concat_dir(!xuvtop,'masterlist'), 'masterlist.ions')
read_masterlist,mname,list

list_spectroscopic = strarr(n_elements(list))
include = bytarr(n_elements(list))

FOR ilist=0,n_elements(list)-1 DO BEGIN
   ion2spectroscopic,list[ilist],snote
   list_spectroscopic(ilist) = snote

   IF n_elements(list_ions) GT 0 THEN BEGIN 
      index = where(list_ions EQ list[ilist], nn)
      IF nn GE  1 THEN include(ilist) = 1b
   ENDIF 
ENDFOR 

IF min(include) EQ 0 THEN all_ions_yn = 0 ELSE all_ions_yn = 1

WIDGET_CONTROL,state.all_ions_ev, set_value=all_ions_yn


; all lines ? 
; check if there are lines with flags = -1 

index = where(tran.lines[*].flag EQ -1, nn)
IF nn GT 0 THEN theor_lines = 1 ELSE theor_lines = 0

WIDGET_CONTROL,state.rem_theor1,  set_value=theor_lines

;keep as default not to have all lines.
WIDGET_CONTROL,state.rem_theor2,  set_value=0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; ioneq

ioneq_name = tran.ioneq_name

IF file_exist(ioneq_name) THEN BEGIN 
   read_ioneq, ioneq_name , ioneq_logt,ioneq,ioneq_ref

;            ioneq_ref= tran.ioneq_ref

   break_file, ioneq_name,  disk, ioneqdir, ioneqfile,ext
   ioneqfile = ioneqfile+ext
   WIDGET_CONTROL,state.ioneq_show,set_value=ioneqfile
ENDIF ELSE BEGIN 
   WIDGET_CONTROL,state.ioneq_show,set_value=' -- '
   result = DIALOG_MESSAGE(['Error, could not find Ion Fractions file: ',$
                            ioneq_name] ,/info)
END 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; isothermal_flag  

IF  required_tags(tran, 'DEM') THEN BEGIN 

   isothermal_flag = 0

   dem_name = tran.dem_name
;               dem_ref= tran.dem_ref
   IF file_exist(dem_name) THEN BEGIN 

      read_dem, dem_name ,dem_logt,dem,dem_ref
      break_file, dem_name,  disk, demdir, demfile,ext
      demfile = demfile+ext

      WIDGET_CONTROL,state.temp_base, set_value=1

      WIDGET_CONTROL,state.isothermal_base1, map=0
      WIDGET_CONTROL,state.dem_base, map=1
      WIDGET_CONTROL,state.dem_show,set_value=demfile

   ENDIF ELSE BEGIN 
      WIDGET_CONTROL,state.dem_show,set_value=' -- '
      result = DIALOG_MESSAGE(['Error, could not find DEM file: ', $
                               dem_name] ,/info)
   END 

ENDIF   ELSE BEGIN 
   
   IF  required_tags(tran, 'logt_isothermal,logem_isothermal') THEN BEGIN 

      isothermal_flag = 1 

      iso_logt = tran.logt_isothermal
      iso_logem = tran.logem_isothermal

      WIDGET_CONTROL,state.temp_base, set_value=0
      WIDGET_CONTROL,state.isothermal_base1, map=1
      WIDGET_CONTROL,state.dem_base, map=0

      dummy1 = iso_logt 
      dummy2 = iso_logem

      dummy1 = trim(float(dummy1))
      dummy2 = trim(float(dummy2))

      WIDGET_CONTROL,state.iso_logt_ev,$
        set_value= arr2str(dummy1, ',',/trim)

      WIDGET_CONTROL,state.iso_logem_ev,$
        set_value=arr2str(dummy2, ',',/trim)

   ENDIF ELSE $
     result = DIALOG_MESSAGE('Error in the structure ! ',/info)

END  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; UNITS: 

CASE  tran.int_units OF
   'erg cm-2 sr-1 s-1':BEGIN 
; ytitle = 'Intensity [ erg/cm!u2!n/sr/s/'+angstrom+' ]'
      units(1) = 'ergs'
      widget_control, state.unit_info, set_val='Units: ERGS'
      widget_control, state.unit_info2, set_val='Units: ERGS'
   END 
   'photons cm-2 sr-1 s-1':BEGIN 
; ytitle = 'Intensity [ photons/cm!u2!n/sr/s/'+angstrom+' ]'
      units(1) = 'photons'
      widget_control, state.unit_info, set_val='Units: PHOTONS'
      widget_control, state.unit_info2, set_val='Units: PHOTONS'
   END
   ELSE: BEGIN 

   END 
ENDCASE 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  WE  ADD HERE THE PROTON RATES switch

CASE  tran.add_protons OF 
   1:BEGIN 
      widget_control,state.proton_info, set_val='Protons: YES'
      noprot = 0
   END 
   0:BEGIN 
      widget_control,state.proton_info, set_val='Protons: NO'
      noprot = 1
   END
   ELSE: BEGIN 
   END 
ENDCASE  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

CASE  tran.photoexcitation OF 
   0:BEGIN 
      widget_control,state.photoexcitation_info, set_val='PE: NO'
      photoexcitation = 0
      WIDGET_CONTROL,state.photoexcitation_base, map=0
   END 
   1:BEGIN 
      widget_control,state.photoexcitation_info, set_val='PE: YES'
      photoexcitation = 1
      WIDGET_CONTROL,state.photoexcitation_base, map=1

      rphot = tran.rphot
      radtemp = tran.radtemp

      WIDGET_CONTROL,state.photoexcitation_rphot_ev, $
        set_value=strtrim(string(format='(f6.3)',rphot), 2)

      WIDGET_CONTROL,state.photoexcitation_radtemp_ev, $
        set_value=strtrim(string(format='(f6.1)',radtemp), 2)


   END 
ENDCASE 

;add the photoexcitation


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


widget_control, state.show_lines , $
  set_val=trim(n_elements(tran.lines))+$
  ' line intensities restored and parameters reset. Now calculate spectrum !'



END 

PRO OPLOT_LINES

;
; need to add a fwhm tag to tran.lines in order for this to work (need to 
; divide intensity of line by FWHM) to get peak of line
;
; NEED TO ADD A CHECK ON THE X-RANGE ?
;

COMMON line_data, tran, state
COMMON plot_spec,kev_yn, ang, inst, o_lines, xrange, yrange, sty, theor_lines, $
  cont_ind, lambda, effarea, file_effarea, spectrum, fold_yn, log
COMMON wind_data, wxsiz, wysiz, o_strength
COMMON abundance, abfile, abdir, abund, abund_ref,min_abund 


;
IF o_lines EQ 1 THEN BEGIN

   ind = WHERE(spectrum.lines[*].peak GT  o_strength AND $
               (spectrum.lines[*].wvl LT  xrange[1]) AND  $
               (spectrum.lines[*].wvl GT  xrange[0]), nn)


   IF nn GT 0 THEN BEGIN 

      continuum = spectrum.lambda-spectrum.lambda
      IF tag_exist(spectrum,  'continuum') THEN continuum = spectrum.continuum

      FOR i=0L,nn-1 DO BEGIN

         peak = spectrum.lines[ind[i]].peak

;find the nearest  position 

         dummy = min(abs(spectrum.lines[ind[i]].wvl-spectrum.lambda), i_cont) 
         i_cont = i_cont[0] &  IF i_cont EQ -1 THEN i_cont=0


         IF log THEN oplot,[1.,1.]*spectrum.lines[ind[i]].wvl,  $
           [peak+ continuum(i_cont), continuum(i_cont) > 1e-10],th=2 ELSE $
           oplot,[1.,1.]*spectrum.lines[ind[i]].wvl, $
           [peak+ continuum(i_cont), continuum(i_cont) >0.],th=2

         snote = spectrum.lines[ind[i]].snote
         IF spectrum.lines[ind[i]].flag EQ -1 THEN snote = snote+' *'

         xyouts,spectrum.lines[ind[i]].wvl,peak+ continuum(i_cont),snote

      ENDFOR
   ENDIF 
END  
END  


PRO SYN_CURSOR, EVENT

;
; When the mouse button is clicked, all lines that are within +/- 1/20 
; of the X-range of the plot will be displayed
;

COMMON line_data, tran, state
COMMON plot_spec,kev_yn, ang, inst, o_lines, xrange, yrange, sty, theor_lines, $
  cont_ind, lambda, effarea, file_effarea, spectrum, fold_yn, log
COMMON wind_data, wxsiz, wysiz, o_strength
COMMON abundance, abfile, abdir, abund, abund_ref,min_abund 

angstrom = string(197B)

Ywindow=[!y.crange[0], !y.crange[1]]

;  Xwindow=[0.,1.]
;  Ywindow=[0.,1.]
;coord=convert_coord(Xwindow,Ywindow, /NORMAL, /TO_DEVICE)
;Xwindow(*)=coord(0,0:1)
;Ywindow(*)=coord(1,0:1)

dataxy = CONVERT_COORD(event.x,event.y,/device,/to_data)
x = dataxy[0] & y=dataxy[1]

range = (!x.crange[1]-!x.crange[0])/20.

Xc_old = x & Yc_old=y

device,GET_GRAPHICS_FUNCTION=orig_mode

;tvcrs,Xc_old,Yc_old

test = (n_elements(spectrum) GT 0) AND  (x GT   !x.crange[0]) AND  (x LT   !x.crange[1]) $
  AND   (y GT   !y.crange[0]) AND  (y LT   !y.crange[1]) AND NOT log


IF  test THEN BEGIN 

   device, SET_GRAPHICS_FUNCTION=6

   catch,error
   if error ne 0 then begin
;back to default error handling
      catch,/CANCEL
      device,SET_GRAPHICS_FUNCTION=orig_mode
      device,/CURSOR_ORIGINAL
      return
   ENDIF

   device,CURSOR_IMAGE=indgen(16)

   repeat begin


      oplot,[Xc_old-range, Xc_old-range], Ywindow, NOCLIP=0
      oplot,[Xc_old+range, Xc_old+range], Ywindow, NOCLIP=0

      cursor,Xc,Yc,/CHANGE, /data

      button=!mouse.button

      IF  button ne 0 then BEGIN 

;      dataxy = CONVERT_COORD(event.x,event.y,/device,/to_data)
;      x = dataxy[0]
         x = Xc

         ind = WHERE((ABS(spectrum.lines[*].wvl - x) LE range) AND spectrum.lines[*].peak GT  o_strength)

         IF ind[0] NE -1 THEN BEGIN

            nlines = N_ELEMENTS(ind)
            line_str = STRARR(nlines)

            FOR i=0L,nlines-1 DO BEGIN

               snote = trim(spectrum.lines[ind[i]].snote)+'  '

               IF spectrum.lines[ind[i]].flag EQ -1 THEN snote = snote+' *'

               line_str[i] = strpad(STRING(format='(f11.4)',$
                                           spectrum.lines[ind[i]].wvl ),10,/after)+'  ( '+$
 spectrum.units[0] +' )  '+$
                 strpad(snote,14,/after)+$
                 '    Int='+ strtrim(STRING(format='(e12.2)', $
                                            spectrum.lines[ind[i]].int), 2)+'  ( '+$
 spectrum.int_units +' )  '

               IF spectrum.lines[ind[i]].tmax NE 0.  THEN line_str[i] = line_str[i]+  $
                 '  log Tm [K]='+trim(STRING(format='(f5.2)', spectrum.lines[ind[i]].tmax))

; GDZ: add the line indices -                
               line_str[i] = line_str[i] + '    '+trim(spectrum.lines[ind[i]].lvl1)+$
                             '-'+trim(spectrum.lines[ind[i]].lvl2) +$
                 '   '+strtrim(spectrum.lines[ind[i]].ident, 2)


               IF  spectrum.lines[ind[i]].flag EQ -1 THEN $
                 line_str[i] = line_str[i]+ '  ("unobserved") '

            ENDFOR

            newind = SORT(spectrum.lines[ind].wvl)
            WIDGET_CONTROL,state.show_lines,set_value=line_str[newind]
         ENDIF ELSE BEGIN
            WIDGET_CONTROL,state.show_lines,set_value=' '
         ENDELSE
      END   

;plots, 
      oplot,[Xc_old-range, Xc_old-range], Ywindow, NOCLIP=0
      oplot,[Xc_old+range, Xc_old+range], Ywindow, NOCLIP=0


      Xc=float(Xc)
      Yc=float(Yc)
      Xc_old=Xc
      Yc_old=Yc

   endrep UNTIL (xc LE  !x.crange[0]) OR (xc GE  !x.crange[1]) $
     OR  (yc LE  !y.crange[0]) OR (yc GE  !y.crange[1]) OR button ne 0 

;  plots,[!x.crange[0], !x.crange[0]],[!y.crange[0], !y.crange[0]]
device,SET_GRAPHICS_FUNCTION=orig_mode
   device,/CURSOR_ORIGINAL


ENDIF 

END

PRO calc_syn_spectrum, err_msg=err_msg

COMMON line_data, tran, state
COMMON plot_spec,kev_yn, ang, inst, o_lines, xrange, yrange, sty, theor_lines, $
  cont_ind, lambda, effarea, file_effarea, spectrum, fold_yn, log
COMMON abundance, abfile, abdir, abund, abund_ref, min_abund
COMMON calc_int, int_xrange, const_names, const_nt, const_value, $
  ioneqfile, ioneqdir, demfile, demdir, iso_logt, iso_logem, isothermal_flag, $
  all_ions_yn, list_ions, units, noprot,photoexcitation, rphot, radtemp

;
; PRY, 18-Jan-2024: commented out line below as it prevents files in the
; abundance/archive folder being read.
;
;IF  file_exist(concat_dir(concat_dir(!xuvtop,'abundance'), abfile)) THEN $
;  abdir = concat_dir(!xuvtop, 'abundance') ELSE cd, current=abdir
abund_name=concat_dir(abdir,abfile)
read_abund,abund_name,abund,abund_ref


IF theor_lines THEN all = 1 ELSE all = 0


delvarx, err_msg

IF units(1) EQ 'photons' THEN photons = 1 ELSE IF $
  units(1) EQ 'ergs' THEN photons =0 ELSE photons =0

;IF units(0) EQ 'Angstroms' THEN kev = 0 ELSE IF units(0) EQ 'keV' THEN kev = 1 ELSE kev = 0

;clear the lambda array:

delvarx, lambda

make_chianti_spec,tran,lambda, spectrum, BIN_SIZE=ang,$ 
  photons=photons, $
  INSTR_FWHM=inst, wrange=xrange, $
  ALL=ALL, continuum=cont_ind , $
  ABUND_NAME=ABUND_NAME, binsize=binsize, $
  MIN_ABUND=MIN_ABUND, file_effarea=file_effarea,$
  err_msg=err_msg, /VERBOSE,kev=kev_yn

;IF ang EQ 0. THEN ang = binsize


END


PRO plot_syn_spectrum

COMMON line_data, tran, state
COMMON plot_spec,kev_yn, ang, inst, o_lines, xrange, yrange, sty, theor_lines, $
  cont_ind, lambda, effarea, file_effarea, spectrum, fold_yn, log
COMMON abundance, abfile, abdir, abund, abund_ref, min_abund

;
; sty
;      This parameter is either 0 or 1 and indicates whether the displayed 
; spectrum has been 'zoomed' (1) or not (0). It affects how the X and Y 
; range of the plot are displayed. It is set to 0 by the 'Unzoom' and 
; 'Make line spectrum' buttons, and set to 1 by the 'Zoom' button.
;

angstrom = string(197B)


;IF sty EQ 0 THEN yrange = [0,max(spectrum)*1.2]
;IF sty EQ 1 THEN BEGIN
;  IF ABS(yrange[0]-0.) LE (!y.crange[1]-!y.crange[0])/20. THEN yrange[0] = 0.
;  IF xrange[0] LT MIN(lambda) THEN xrange[0] = MIN(lambda)
;  IF xrange[1] GT MAX(lambda) THEN xrange[1] = MAX(lambda)
;ENDIF

IF n_elements(spectrum) EQ 0 THEN BEGIN 
   tverase
   return
END 

IF spectrum.units(0) EQ 'Angstroms' THEN  xtitle = 'Wavelength ('+angstrom+')' ELSE $
  xtitle = spectrum.units(0)

;replace some things for better look:
ytitle = spectrum.units(1)
ytitle = str_replace(ytitle, 'Angstroms', angstrom)
ytitle = str_replace(ytitle,'-1', '!U-1!N')
ytitle = str_replace(ytitle,'-2', '!U-2!N')
ytitle = str_replace(ytitle,'-3', '!U-3!N')


IF log THEN BEGIN
  ;
  ; PRY, 19-Nov-2020: added the following to make log plots display
  ; better.
  ;
   IF yrange[0] EQ 0. THEN BEGIN
      yrange[0]=min(spectrum.spectrum)
      widget_control,state.ymin_base,set_value=string(format='(e9.2)',yrange[0])
   ENDIF 
   plot_io, spectrum.lambda, spectrum.spectrum,psym=10, $
            xmargin=[15,4], $
            xsty=1,xrange=xrange,yrange=[1e-10 > yrange(0), yrange(1)],ysty=1, $
            xticklen=-0.025, $
            ytitle=ytitle,xtitle=xtitle, title='CHIANTI - Version '+spectrum.version
ENDIF ELSE BEGIN 
   plot, spectrum.lambda, spectrum.spectrum,psym=10, $
         xmargin=[15,4], $
         xsty=1,xrange=xrange,yrange=yrange,ysty=1, $
         xticklen=-0.025,yticklen=-0.015, $
         ytitle=ytitle,xtitle=xtitle, title='CHIANTI - Version '+spectrum.version
ENDELSE 

;xyouts some precious information

CASE  spectrum.model_name OF 
   'Constant density': dummy = 'Calculated with '+spectrum.model_name+ '= '+$
     string(spectrum.model_ne,'(e9.2)')+ ' (cm-3) '

   'Constant pressure':dummy = 'Calculated with '+spectrum.model_name+ '= '+$
     string(spectrum.model_pe,'(e9.2)')+ ' (cm-3 K) '

   'Function': BEGIN 
      break_file, spectrum.model_file, disk,dir,f,ext
      dummy = 'Calculated with (Te, Ne) file: '+f+ext
   END

   ELSE: dummy=''
ENDCASE 

IF dummy NE '' THEN $
  xyouts,0.300,0.9,/normal, dummy ,   chars=1.

break_file, spectrum.ioneq_name, disk,dir,f,ext
f = f+ext
xyouts,0.300,0.87,/normal, 'Ioniz. Frac. file : '+f,chars=1. 

IF tag_exist(spectrum, 'DEM') THEN BEGIN 
   isothermal_flag = 0
   break_file, spectrum.dem_name, disk,dir,f,ext
   dummy = 'DEM file: '+f+ext
ENDIF ELSE BEGIN 
   dummy1 = spectrum.logt_isothermal
   dummy2 = spectrum.logem_isothermal

         dummy1 = trim(float(dummy1))
         dummy2 = trim(float(dummy2))

   dummy = 'Isothermal approx.: log T='+arr2str(dummy1, ',',/trim)+$
     ' log EM='+ arr2str(dummy2, ',',/trim)

END 


xyouts,0.300,0.84,/normal,  dummy, chars=1. 

break_file, abfile, disk,dir,f,ext
f = f+ext
xyouts,0.300,0.81,/normal, 'Abundance file : '+f,chars=1. 

IF  min_abund GT 0. THEN BEGIN 
   xyouts,0.300,0.78,/normal, 'Minimum abundance: '+strtrim(string(MIN_ABUND, format='(e9.2)'),2),chars=1. 
ENDIF 

IF fold_yn THEN BEGIN 
   break_file,spectrum.file_effarea, disk,dir,f,ext
   f = f+ext
   xyouts,0.300,0.75,/normal, 'Effective area: '+f,chars=1. 
END

IF o_lines EQ 1 THEN oplot_lines

END



;------------------------------------------------------------------------------
PRO syn_MAIN_Event, Event
;------------------------------------------------------------------------------

  COMMON base_com, syn_main_base , isothermal_base1
  COMMON line_data, tran, state
  COMMON plot_spec,kev_yn, ang, inst, o_lines, xrange, yrange, sty, theor_lines, $
     cont_ind, lambda, effarea, file_effarea, spectrum, fold_yn, log
  COMMON wind_data, wxsiz, wysiz, o_strength
  COMMON abundance, abfile, abdir, abund, abund_ref, min_abund
  COMMON calc_int, int_xrange, const_names, const_nt, const_value, $
     ioneqfile, ioneqdir, demfile, demdir, iso_logt, iso_logem, isothermal_flag, $
     all_ions_yn, list_ions, units, noprot,photoexcitation, rphot, radtemp


  angstrom = string(197B) 
;
;  list of elements
;

  elements=['H','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg',$
            'Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr',$
            'Mn','Fe','Co','Ni', 'Cu','Zn']


  WIDGET_CONTROL,Event.top, get_uvalue=state ;,/no_copy


  CASE 1 OF 

     event.id EQ state.calc_int_butt: BEGIN 


        str=['This is a short HELP. For more details on single buttons please click the buttons where the " - HELP" sign is. ', $
             ' ', $
             'PLEASE READ THE CHIANTI USER GUIDES AND THE DOCUMENTATION IN THE HEADERS OF THE PROCEDURES to understand how to use this software. ', $
             '  ', $
             ' Send comments to Giulio Del Zanna if the details in the documentation are not clear or if you have problems in running the software.', $
             '  ', $
             'This section of the program calculates line intensities (with a call to CH_SYNTHETIC).', $
             'A series of parameters must be set: ', $
             ' - Minimum and maximum wavelengths in Angstroms ', $
             ' - The model used for the calculation. Three are the options:', $
             '   1) a constant  density (cm^-3)', $
             '   2) a constant pressure (cm^-3 K) ', $
             '   3) a general (Te,Ne) model. In this case, a file will be read. ', $
             '      This file should have two columns, one with the Te (K) values, and one with the Ne (cm^-3) values. ', $
             ' ', $
             ' - The ionization fraction file to be used or the temperature range and step for a calculation on the fly with the version 11 advanced models. ', $
             ' ',$
             ' - All ions ?  If set to yes (default), then all the ions present in the database will be included.', $
             '               If set to no, then it is possible to select a list of ions with a widget', $
             ' ', $
             ' - All lines ? If set to no, only the lines for which there are observed energy levels are included. ', $
             '               If set to yes (default), also the lines that do not have corresponding observed energy levels are included. ', $
             '               In this case, the wavelengths are calculated from the theoretical energy levels, and might not be very accurate.', $
             ' ', $
             ' - Isothermal ?   If set to no (default), a DEM file must be selected. ', $
             '                  Files either in the CHIANTI database or in the working directory can be chosen. Otherwise selected with a widget.', $
             '                 If set to yes, then Log T and Log EM values must be defined (see HELP button).', $
             ' - Photoexcitation ? If set to yes, you have to define:', $
             '   Trad: The blackbody radiation field temperature', $
             '   R/Ro: Distance from the centre of the star in stellar radius units ', $
             ' - Units:  Photons or Ergs', $
             ' -Protons: If set to Yes, the proton data are used to calculate the level population', $
             ' ', $
             ' Once all the parameters have been defined, the user should click on the "Calculate intensities" button ', $
             ' to start the calculation (with a call to CH_SYNTHETIC). ', $
             ' ', $
             ' Once the calculation is finished, an IDL CHIANTI LINE INTENSITY STRUCTURE is loaded into memory. ', $
             ' It is then possible to save it for later use by clicking  on the "SAVE" button, with two options: ', $
             '    1- an IDL file (this is done by calling SAVEGEN ).', $
             '       NOTE: to restore the above file into an IDL stucture STRUCT, use (at the command line):  IDL> restgen,file=filename,struct=STRUCT ', $
             '    2- a FITS binary table (this is done by calling CH_WRITE_FITS).',$
             '       NOTE: to restore the above file into an IDL stucture STRUCT, use (at the command line):  IDL> ch_read_fits, filename, STRUCT ', $

             '    The RESTORE option is to restore previously-saved files into  memory. ', $
             ' ', $
             ' Once the IDL structure with the line intensities is in the memory, it is then possible to ', $
             ' calculate and plot a spectrum. Please read the HELP below. ', $
             ' ', $
             'Giulio Del Zanna ']
        xpopup,str, xsize=130
     END

;-----------------------------------------------------------------------

     event.id EQ state.wminw1 :BEGIN 
        val = 0.
        widget_control,state.wminw1, get_v=val 
        IF valid_num(val(0)) AND val(0) GT 0. THEN BEGIN 
           int_xrange(0) = float(val(0)) 
        ENDIF ELSE $
           WIDGET_CONTROL,state.wminw1,  set_value=strtrim(string(0.0, format='(f11.4)'),2)

     END

     event.id EQ state.wmaxw1 :BEGIN 
        val = 0.
        widget_control,state.wmaxw1, get_v=val 
        IF valid_num(val(0)) AND val(0) GT 0. THEN BEGIN 
           int_xrange(1) = float(val(0)) 
        ENDIF ELSE $
           WIDGET_CONTROL,state.wmaxw1,  set_value=strtrim(string(0.0, format='(f11.4)'),2)
     END


;-----------------------------------------------------------------------

     event.id EQ state.const_widg: BEGIN

        const_nt = fix(event.value-1)

        CASE const_nt OF 
           0:BEGIN 
              WIDGET_CONTROL,state.const_read,sensitive=1, set_value=strtrim(string(0.0, format='(e9.2)'),2)
           END 
           1:BEGIN 
              WIDGET_CONTROL,state.const_read, sensitive=1, set_value=strtrim(string(0.0, format='(e9.2)'),2)
           END 
           2:BEGIN 

              widget_control, state.show_lines , $
                              set_val=' Will be using a Functional (Ne,Te) form '

              ptitle='Select (Te, Ne) File (TWO COLUMNS - Temperatures and Densities)'
              const_value=dialog_pickfile(title=ptitle)

              ;; const_value =  bigpickfile(title= 'Select (Te, Ne) File (TWO COLUMNS - Temperatures and Densities  )', $
              ;;                            filter='*' )

              IF file_exist(const_value) THEN BEGIN 

                 break_file, const_value,  disk, dir, file,ext
                 file = file+ext

                 widget_control, state.const_read, sensitive=0, set_value=file

              ENDIF ELSE BEGIN 

                 WIDGET_CONTROL,state.const_read, sensitive=1, set_value=' - '
                 result = DIALOG_MESSAGE('Error, no  (Ne,Te) file defined ! ',/info)

              END  
           END 
        ENDCASE    
     END 

;-----------------------------------------------------------------------

     event.id EQ state.const_read: BEGIN

        WIDGET_CONTROL,state.const_widg,get_uvalue=bob
        const_nt = fix(bob[0])

        CASE const_nt OF 
           0:BEGIN 
              WIDGET_CONTROL,state.const_read,get_value=bob
              const_value = float(bob[0]) 
              widget_control, state.show_lines , $
                              set_val='Will be using constant '+const_names(const_nt)+'='+string(const_value)
           END
           1:BEGIN 
              WIDGET_CONTROL,state.const_read,get_value=bob
              const_value = float(bob[0]) 
              widget_control, state.show_lines , $
                              set_val='Will be using constant '+const_names(const_nt)+'='+string(const_value)
           END
           2:BEGIN 
;      WIDGET_CONTROL,state.const_read,get_value=bob
           END 
        END 
     END  

;-----------------------------------------------------------------------

     event.id EQ state.ioneq_help_button: BEGIN 
        str=['With the advanced ionization equilibrium option (the default), developed by Giulio Del Zanna and Roger Dufresne', $
             'and implemented in version 11, the default is to calculate on the fly the ion charge states ', $
             'in equilibrium, for a set of elements/ions that are available. For the other ions the previous ',$
             'zero-density approximation is used. To speed up the calculations, the user can choose the required',$
             'temperature range and step in log T.  ',$
             'The ion fractions are stored in a file for later use within the programs called by ch_ss.',$
             'This file can be used later on as an input. The file is written in the working directory, with the name adv.ioneq,',$
             'but if such file exists, a different file with a time stamp added (in Julian days and fractions) is chosen by default. ',$
             ' ',$
             'The addition of charge transfer (CT) is recommended especially for silicon and oxygen ions. In this ',$
             'case, a selection of atmosphere files can be chosen, or a user-defined supplied. ',$
             ' ',$
             'The alternative is to choose as in previous versions a precompiled CHIANTI .ioneq file with the charge states.',$
             'This file could have been pre-calculated with the advanced model options using ch_ss or the program CH_CALC_IONEQ,',$
             'but the user should make sure that for consistency the same parameters are used for both the ion abundances',$
             'and the calculation of the line emissivities;  as before, three options are available, constant Ne, ',$
             'constant Pe, or a tabulated list (in a file) of Te, Ne ',$
             ' ']
        result = DIALOG_MESSAGE(str,/info)
     END

;-----------------------------------------------------------------------   

     
     event.id EQ state.ioneq_pdmenu: BEGIN

        WIDGET_CONTROL,state.ioneq_pdmenu,  get_uvalue=bob

        ioneqfile = ''
        file = ''
        
        CASE  bob OF 

           '0': BEGIN 
              WIDGET_CONTROL,state.ioneq_show,set_value=''
           END 

           '1':BEGIN 
              ptitle='Select appropriate Ion Fraction file to READ'
              file=dialog_pickfile(filter='*ioneq', tit=ptitle)
;            file = BIGPICKFILE(filter='*ioneq', tit='Select appropriate Ion Fraction file to READ')
              IF file NE '' THEN BEGIN 
                 break_file,file,  disk, ioneqdir, ioneqfile,ext
                 ioneqfile = ioneqfile+ext
                 ioneqdir = concat_dir(disk, ioneqdir)
                 WIDGET_CONTROL,state.ioneq_show,set_value=ioneqfile
              END
           END 

           ELSE:BEGIN 
              file = bob
              break_file,file,  disk, ioneqdir, ioneqfile,ext
              ioneqfile = ioneqfile+ext
              ioneqdir = concat_dir(disk, ioneqdir)
              WIDGET_CONTROL,state.ioneq_show,set_value=ioneqfile
           END   
        ENDCASE  

        IF file NE '' THEN BEGIN 
           read_ioneq, file, ioneq_logt,ioneq,ioneq_ref
           widget_control, state.show_lines,set_val= ''

           FOR  i=0,n_elements(ioneq_ref)-1 DO $
              widget_control, state.show_lines,/append , $
                              set_val=ioneq_ref[i]
        END 

     END  

;-----------------------------------------------------------------------
; event for the calculation of ion charge states
     
     event.id EQ state.ibase_sel:BEGIN 

        WIDGET_CONTROL,state.ibase_sel,get_uvalue=bob,get_value=fred
        calc_ioneq_flag = bob[fred]
        
        if calc_ioneq_flag eq 1 then begin 
           WIDGET_CONTROL,state.input_ibase1, map=1
           WIDGET_CONTROL,state.input_ct_base, map=1
           WIDGET_CONTROL,state.output_ioneq_show, map=1
           
; by default do not run CT:          
           WIDGET_CONTROL,state.ct_base, map=0
           WIDGET_CONTROL,state.ct_base_sel,set_value=1
           
           WIDGET_CONTROL,state.ioneq_base,map=0
        endif else begin
           WIDGET_CONTROL,state.output_ioneq_show, map=0
           WIDGET_CONTROL,state.input_ibase1, map=0
           WIDGET_CONTROL,state.input_ct_base, map=0
           WIDGET_CONTROL,state.ct_base, map=0
           WIDGET_CONTROL,state.ioneq_base,map=1
        end       
     end 
;-----------------------------------------------------------------------

     event.id EQ state.ct_base_sel:BEGIN 

        WIDGET_CONTROL,state.ct_base_sel,get_uvalue=bob,get_value=fred
        ct_flag = bob[fred]

        
        if ct_flag eq 1 then begin

           WIDGET_CONTROL,state.ct_base, map=1
           
           atmosphere_file=ch_get_file(path=concat_dir(concat_dir($
                           concat_dir(!xuvtop,'ancillary_data'),$
                           'advanced_models'),'model_atmospheres'),$
                                       filter='*.dat',tit='Select an atmosphere file', $
                                       /modal,  group=SYN_MAIN_base)
           
           if file_exist(atmosphere_file) then begin
              break_file, atmosphere_file,  disk, dir, file,ext
              
              WIDGET_CONTROL,state.ct_show,set_value=file
              
           endif  else begin
              WIDGET_CONTROL,state.ct_show,set_value=''
              result = DIALOG_MESSAGE(['Error, could not find atmosphere file '] ,/info)
           end 
           
        endif  else begin
           atmosphere_file=''
           WIDGET_CONTROL,state.ct_show,set_value=''
           WIDGET_CONTROL,state.ct_base, map=0
        end 
        
     end 
;-----------------------------------------------------------------------
; check input min and max log T
     event.id EQ state.min_logt_ev:BEGIN
        dummy=''
        WIDGET_CONTROL,state.min_logt_ev, get_v=dummy
        IF valid_num(dummy(0))  THEN BEGIN 
           state.min_logt= float(dummy[0])
           if state.min_logt ge state.max_logt then begin
              widget_control, state.show_lines , $
                              set_val='min log T  greater than max log T ??? '
              WIDGET_CONTROL,state.min_logt_ev,$ 
                             set_value= strtrim(string(format='(f3.1)', 4.0 ))
              state.min_logt=4.0
           endif 
        ENDIF ELSE BEGIN 
           widget_control, state.show_lines , $
                           set_val='min log T  not a  valid number !'
           WIDGET_CONTROL,state.min_logt_ev,$ 
                          set_value= strtrim(string(format='(f3.1)', 4.0 ))
           state.min_logt=4.0
        END 
     END 

     event.id EQ state.max_logt_ev:BEGIN
        dummy=''
        WIDGET_CONTROL,state.max_logt_ev, get_v=dummy
        IF valid_num(dummy(0))  THEN BEGIN 
           state.max_logt= float(dummy[0])
           if state.max_logt le state.min_logt then begin
              widget_control, state.show_lines , $
                              set_val='max log T smaller  than min log T ??? '
              WIDGET_CONTROL,state.max_logt_ev,$ 
                             set_value= strtrim(string(format='(f3.1)', 8.0 ))
              state.max_logt=8.0
           endif 
        ENDIF ELSE BEGIN 
           widget_control, state.show_lines , $
                           set_val='max log T  not a  valid number !'
           WIDGET_CONTROL,state.max_logt_ev,$ 
                          set_value= strtrim(string(format='(f3.1)', 8.0 ))
           state.max_logt=8.0
        END 
     END 

     event.id EQ state.dlogt_ev:BEGIN
        dummy=''
        WIDGET_CONTROL,state.dlogt_ev, get_v=dummy
        IF valid_num(dummy(0))  THEN BEGIN 
           state.dlogt= float(dummy[0])
           if state.dlogt le 0 or  state.dlogt ge  state.min_logt then begin
              widget_control, state.show_lines , $
                              set_val='wrong d log T  ??? '
              WIDGET_CONTROL,state.dlogt_ev,$ 
                             set_value= strtrim(string(format='(f3.1)', 0.1 ))
              state.dlogt=0.1
           endif 
        ENDIF ELSE BEGIN 
           widget_control, state.show_lines , $
                           set_val='d log T  not a  valid number !'
           WIDGET_CONTROL,state.dlogt_ev,$ 
                          set_value= strtrim(string(format='(f3.1)', 0.1 ))
           state.dlogt=0.1
        END 
     END 
     
;-----------------------------------------------------------------------
; check the name of the output .ioneq file

     event.id EQ state.output_ioneq_show:BEGIN
        dummy=''
        WIDGET_CONTROL,state.output_ioneq_show, get_v=dummy
        dummy=trim(dummy)
        if file_exist(dummy) then begin

           state.output_ioneq_file='adv_'+$
                                   trim(string(SYSTIME(/JULIAN, /UTC), format='(f14.4)'))+'.ioneq'
           
           WIDGET_CONTROL,state.output_ioneq_show, set_v=state.output_ioneq_file
           
           widget_control, state.show_lines , $
                           set_val='Output ioneq file '+dummy+' exist ! - setting a new file name: '+state.output_ioneq_file

        endif else state.output_ioneq_file=dummy
        
     END
     
;-----------------------------------------------------------------------

     event.id EQ state.all_ions_butt: BEGIN
        str=['Calculate all the ions in the database.', $
             'If this is set, the routine uses the file masterlist.ions that is located ', $
             'in the masterlist/ directory' ]
        result = DIALOG_MESSAGE(str,/info)
     END 
;-----------------------------------------------------------------------         
     
     event.id EQ state.all_ions_ev: BEGIN

        WIDGET_CONTROL,state.all_ions_ev,get_uvalue=bob,get_value=fred
        all_ions_yn= bob[fred]

        IF all_ions_yn EQ 0 THEN BEGIN 

           mname=concat_dir(concat_dir(!xuvtop,'masterlist'), 'masterlist.ions')
           read_masterlist,mname,list

           list_spectroscopic = strarr(n_elements(list))
           include = bytarr(n_elements(list))
           list_iz = intarr(n_elements(list))
           list_ion = intarr(n_elements(list))

           FOR ilist=0,n_elements(list)-1 DO BEGIN
              ion2spectroscopic,list[ilist],snote
              list_spectroscopic(ilist) = snote

              convertname,list[ilist],iz,ion
              list_iz(ilist) = iz
              list_ion(ilist) = ion

              IF n_elements(list_ions) GT 0 THEN BEGIN 
                 index = where(list_ions EQ list[ilist], nn)
                 IF nn EQ 1 THEN include(ilist) = 1b
              ENDIF 
           ENDFOR 

;sort the elements

           fsort = sort(list_iz)

;
;  get a list of elements 
;

           dummy = list_iz[fsort]
           index = rem_dup(dummy)
           list_iz_only = dummy(index)
           isort = -1


           FOR i=0, n_elements(list_iz_only)-1 DO BEGIN 

              ii = where(list_iz EQ list_iz_only(i))

              isort = [isort, ii(sort(list_ion[ii]))]

           ENDFOR 

           isort = isort[1:*]

           list_spectroscopic = list_spectroscopic[isort]
           include=include[isort]
           list = list[isort]

;         index = ch_xmenu_sel(list_spectroscopic, include=include, nl=30,
;         group=SYN_MAIN_BASE)
           index = ch_xmenu_sel(list_spectroscopic, include=include,$
                                tit=' Select ions', group=SYN_MAIN_BASE)

           IF index[0]  NE -1 THEN BEGIN 
              all_ions_yn=0
              list_ions=list[index] 

           ENDIF ELSE BEGIN 

              widget_control, state.show_lines , $
                              set_val='Error, resetting the use of all ions'

              all_ions_yn= 1
              WIDGET_CONTROL,state.all_ions_ev, set_value=all_ions_yn

           ENDELSE 
        END 
     END  

;-----------------------------------------------------------------------

     event.id EQ state.rem_theor1_butt: BEGIN
        str=['Add lines that have only theoretical energy levels ? ', $
             ' ', $
             'Please note that  the wavelengths of these  "unobserved lines"', $
             'might not be very accurate']
        result = DIALOG_MESSAGE(str,/info)
     END

;-----------------------------------------------------------------------

     event.id EQ state.rem_theor1: BEGIN
        WIDGET_CONTROL,state.rem_theor1,get_uvalue=bob,get_value=fred
        theor_lines = bob[fred]

        WIDGET_CONTROL,state.rem_theor2,  set_value=theor_lines
     END 

;-----------------------------------------------------------------------

     event.id EQ state.lookup_butt: BEGIN
        str=['Use lookup tables to compute level populations? ', $
             ' ', $
             'This results in a much faster calculation, but with a slight loss of accuracy', $
             '(typically 1% or less).']
        result = DIALOG_MESSAGE(str,/info)
     END

     event.id EQ state.lookup_ev: BEGIN
       widget_control,state.lookup_ev,set_value=event.value
       IF event.value EQ 1 THEN BEGIN
         IF photoexcitation THEN BEGIN 
           widget_control,state.photoexcitation_info, set_val='PE: NO'
           photoexcitation = 0
           WIDGET_CONTROL,state.photoexcitation_base, map=0
         ENDIF 
         widget_control,state.pexc_base,sens=0
       ENDIF ELSE BEGIN
         widget_control,state.pexc_base,sens=1
       ENDELSE 
     END
;-----------------------------------------------------------------------
     
     event.id EQ state.temp_base:BEGIN 

        WIDGET_CONTROL,state.temp_base,get_uvalue=bob,get_value=fred
        isothermal_flag = bob[fred]

        IF isothermal_flag EQ 1 THEN BEGIN 
           WIDGET_CONTROL,state.isothermal_base1, map=1
           WIDGET_CONTROL,state.dem_base, map=0
        ENDIF ELSE BEGIN 
           WIDGET_CONTROL,state.isothermal_base1, map=0
           WIDGET_CONTROL,state.dem_base, map=1
        END 
     END 

;-----------------------------------------------------------------------

     event.id EQ state.temp_butt: BEGIN 
        str=['The program can calculate line intensities either with an ', $
             '1) isothermal approximation', $
             '2) or by using a DEM(T) file.', $
             '', $
             'In the first case, the user is requested to enter one or more temperatures', $
             '(as logarithmic values - Log T )', $
             'and correspondent column emission measures EM logarithmic values', $
             '', $
             'In the second case, the user can select any *.dem file that is either ', $
             'in the standard CHIANTI database or in the working directory', $
             '', $
             'NOTE: if more than one value is entered, then the sequence must be ', $
             'separated by commas (e.g.: 6.0, 6.5, 7.), and both Log T and Log EM must have the same number of values']
        result = DIALOG_MESSAGE(str,/info)
     END

;-----------------------------------------------------------------------

     event.id EQ state.iso_logt_ev:BEGIN 
        dummy = ''
        WIDGET_CONTROL,state.iso_logt_ev, get_v=dummy
        iso_logt = float(str_sep(dummy, ',', /remove))

;      widget_control, state.show_lines , /append, $
;        set_val=' Will be using an isothermal log T='+string(iso_logt)

     END  

;-----------------------------------------------------------------------

     event.id EQ state.iso_logem_ev:BEGIN 
        dummy = ''
        WIDGET_CONTROL,state.iso_logem_ev, get_v=dummy
        iso_logem = float(str_sep(dummy, ',', /remove))


;      widget_control, state.show_lines , $
;        set_val=' Will be using  isothermal log EM='

     END 

;-----------------------------------------------------------------------

     event.id EQ state.dem_pdmenu: BEGIN

;      fix(event.value-1)

        WIDGET_CONTROL,state.dem_pdmenu,  get_uvalue=bob

        demfile =''
        file = ''

        CASE  bob OF 

           '0': BEGIN 
              WIDGET_CONTROL,state.dem_show,set_value=''
           END 

           '1':BEGIN 
              ptitle='Select appropiate DEM file to READ'
              file = dialog_pickfile(filter='*dem', tit=ptitle)
;            file = BIGPICKFILE(filter='*dem', tit='Select appropiate DEM file to READ')
              IF file NE '' THEN BEGIN 
                 break_file,file,  disk, demdir, demfile,ext
                 demfile = demfile+ext
                 demdir = concat_dir(disk,demdir)
                 WIDGET_CONTROL,state.dem_show,set_value=demfile
              END 
           END 

           ELSE:BEGIN 
              file = bob
              break_file,file,  disk, demdir, demfile,ext
              demfile = demfile+ext
              demdir = concat_dir(disk,demdir)
              WIDGET_CONTROL,state.dem_show,set_value=demfile
           END   
        ENDCASE  

        IF file NE '' THEN BEGIN 

           read_dem, file ,dem_logt,dem,dem_ref
           
           widget_control, state.show_lines,set_val= ''

           FOR  i=0,n_elements(dem_ref)-1 DO $
              widget_control, state.show_lines,/append , $
                              set_val=dem_ref[i]


           IF n_elements(ioneq_name) EQ 0 THEN ioneq_name = !ioneq_file

           read_ioneq, ioneq_name,ioneq_logt,ioneq,ioneq_ref

           n_ioneq_logt=n_elements(ioneq_logt)

           gdt=WHERE((ioneq_logt GE MIN(dem_logt)) AND  (ioneq_logt LE MAX(dem_logt))) 
           dem_int1=10.^(SPLINE(dem_logt,dem,ioneq_logt[gdt], 10))
           dem_int=dblarr(n_ioneq_logt)
           ngt=n_elements(gdt)
           FOR igt=0,ngt-1 DO  dem_int[gdt[igt]]=dem_int1[igt] >  0.


           window,0,xs=600,ys=600
           wset,0

           nb = strpos(dem_ref(0),':')

           plot, ioneq_logt, alog10(dem_int),psym=-5,syms=1.2,$
                 tit=strmid(dem_ref(0),nb+1,100),xtit='Log T (K)',$
                 ytit='Log DEM (cm!S!E-5!N K!S!E-1!N)',chars=1.4, /yno

;      plot_oo,10.^dem_logt,10.^dem,psym=-5,syms=1.4,$

           wait, 2
           wshow,!d.window,0

        END

     END

;-----------------------------------------------------------------------

     event.id EQ state.unit_info: BEGIN

;
;  Change UNITS
;
;first remove any spectrum structure still around:
        delvarx, spectrum

        if units(1) EQ 'photons'  then begin
           widget_control, state.unit_info, set_val='Units: ERGS'
           widget_control, state.unit_info2, set_val='Units: ERGS'
           units(1) = 'ergs'
        endif else begin
           widget_control, state.unit_info, set_val='Units: PHOTONS'
           widget_control, state.unit_info2, set_val='Units: PHOTONS'
           units(1) ='photons'
        endelse
     END 

;-----------------------------------------------------------------------

     event.id EQ state.proton_info: BEGIN

        IF noprot THEN BEGIN 
           widget_control,state.proton_info, set_val='Protons: YES'
           noprot = 0
        ENDIF ELSE BEGIN 
           widget_control,state.proton_info, set_val='Protons: NO'
           noprot = 1
        ENDELSE 
     END  

;-----------------------------------------------------------------------

     event.id EQ state.pe_help_button: BEGIN 
        str=['Photo-excitation (PE) by a black-body can be added to the calculations.', $
             'The two parameters are the distance from the Sun centre in solar radii ', $
             'and the radiation temperature of the photosphere, which is about 5800 K  ', $
             'for the visible and 6100 K for the near infrared. ',$
             'Note that a user-defined radiation spectrum could be input in the command-line version.',$
             'Important note: the PE is added as a correction to the A-value. It only works if the ',$
             'spectrum is a continuum, otherwise see Del Zanna (2024). ']
        result = DIALOG_MESSAGE(str,/info)
     END

     
;-----------------------------------------------------------------------
;  Switch
;
     event.id EQ state.photoexcitation_info: BEGIN
        IF photoexcitation THEN BEGIN 
           widget_control,state.photoexcitation_info, set_val='PE: NO'
           photoexcitation = 0
           WIDGET_CONTROL,state.photoexcitation_base, map=0
        ENDIF ELSE BEGIN 
           widget_control,state.photoexcitation_info, set_val='PE: YES'
           photoexcitation = 1
           WIDGET_CONTROL,state.photoexcitation_base, map=1
        END 
     END 

;-----------------------------------------------------------------------

     event.id EQ state.photoexcitation_rphot_ev:BEGIN 
        dummy = ''
        WIDGET_CONTROL,state.photoexcitation_rphot_ev, get_v=dummy
        IF valid_num(dummy(0)) AND dummy(0) GE  1. THEN BEGIN 
           rphot= float(dummy) & rphot=rphot(0)
        ENDIF ELSE BEGIN 
           widget_control, state.show_lines , $
                           set_val='R/Ro not a  valid number !'
           WIDGET_CONTROL,state.photoexcitation_rphot_ev, $
                          set_value= strtrim(string(format='(f6.3)', 1.0 ))
        END 
     END 

;-----------------------------------------------------------------------

     event.id EQ state.photoexcitation_radtemp_ev:BEGIN 
        dummy = ''
        WIDGET_CONTROL,state.photoexcitation_radtemp_ev, get_v=dummy
        IF valid_num(dummy(0)) AND dummy(0) GT 0. THEN BEGIN 
           radtemp= float(dummy) & radtemp=radtemp(0)
        ENDIF ELSE BEGIN 
           widget_control, state.show_lines , $
                           set_val='Trad not a  valid number !'
           WIDGET_CONTROL,state.photoexcitation_radtemp_ev,$ 
                          set_value= strtrim(string(format='(f6.1)', 6000.0 ))
        END 
     END 

;-----------------------------------------------------------------------
;;;;;;;;;;;;;;DO LINE CALCULATION ;;;;;;;;;;;;;;;;

     event.id EQ state.calc_lines : BEGIN

        err_mess = ''

; DO SOME CHECKS:
;----------------

; check min and max wavelengths

        minw1 = xrange(0)
        widget_control,state.wminw1, get_v=minw1 

        IF valid_num(minw1(0)) AND minw1(0) GT 0. THEN BEGIN 
           int_xrange(0) = float(minw1(0) > 0.) 
        ENDIF ELSE BEGIN 
           err_mess = [err_mess, 'No min. wavelength defined']
           WIDGET_CONTROL,state.wminw1,  set_value=strtrim(string(0.0, format='(f11.4)'),2)
        END 

        widget_control,state.wmaxw1, get_v=maxw1 
        IF valid_num(maxw1(0)) AND maxw1(0) GT 0. THEN BEGIN 
           int_xrange(1) = float(maxw1(0) ) 
        ENDIF ELSE BEGIN 
           err_mess = [err_mess, 'No max. wavelength defined']
           WIDGET_CONTROL,state.wmaxw1,  set_value=strtrim(string(0.0, format='(f11.4)'),2)
        END 

        IF int_xrange(1) LE int_xrange(0) THEN BEGIN 
           err_mess = [err_mess, 'Wrong  wavelength range defined']
           WIDGET_CONTROL,state.wminw1,  set_value=strtrim(string(0.0, format='(f11.4)'),2)
           WIDGET_CONTROL,state.wmaxw1,  set_value=strtrim(string(0.0, format='(f11.4)'),2)
        END 

;define pressure  density  temperature  :
        delvarx, pressure, density, temperature , model_file 


        WIDGET_CONTROL,state.const_widg,get_uvalue=bob 
        const_nt = fix(bob[0])

        CASE const_nt OF 
           0: BEGIN 
              WIDGET_CONTROL,state.const_read,get_value=bob
              const_value = float(bob[0]) 

              IF NOT valid_num(const_value) OR const_value LE  0 THEN BEGIN 
                 err_mess = [err_mess, 'Wrong density value']
                 WIDGET_CONTROL,state.const_read, sens=1, set_value=strtrim(string(0.0, format='(e9.2)'),2)
              ENDIF ELSE density=const_value

           END 
           1:BEGIN 
              WIDGET_CONTROL,state.const_read,get_value=bob
              const_value = float(bob[0]) 
              IF NOT valid_num(const_value) OR const_value LE  0 THEN BEGIN 
                 err_mess = [err_mess, 'Wrong pressure value']
                 WIDGET_CONTROL,state.const_read, sens=1,  set_value=strtrim(string(0.0, format='(e9.2)'),2)
              ENDIF ELSE pressure = const_value 

           END 
           2:BEGIN 

              model_file = const_value

           END 
        ENDCASE 

;-------------------------------------------------------
        

        WIDGET_CONTROL,state.all_ions_ev,get_uvalue=bob,get_value=fred
        all_ions_yn= bob[fred]

        WIDGET_CONTROL,state.rem_theor1,get_uvalue=bob,get_value=fred
        theor_lines = bob[fred]

        
        WIDGET_CONTROL,state.temp_base,get_uvalue=bob,get_value=fred
        isothermal_flag = bob[fred]

        IF  isothermal_flag EQ 1 THEN BEGIN 

           dummy = ''
           WIDGET_CONTROL,state.iso_logt_ev, get_v=dummy
           dummy1 = str_sep(dummy, ',', /remove)
           dummy = ''
           WIDGET_CONTROL,state.iso_logem_ev, get_v=dummy
           dummy2 = str_sep(dummy, ',', /remove)

           IF n_elements(dummy1) NE n_elements(dummy2) THEN $
              err_mess = [err_mess, 'Different number of Log T and Log Em values!'] ELSE BEGIN 

              iso_logt =fltarr(n_elements(dummy1))
              iso_logem =fltarr(n_elements(dummy1))

              FOR i=0, n_elements(dummy1)-1 DO BEGIN 
                 IF valid_num(dummy1[i]) THEN iso_logt[i] = float(dummy1[i]) ELSE $
                    err_mess = [err_mess, 'Not a  valid number !']
                 IF valid_num(dummy2[i]) THEN  iso_logem[i]  = float(dummy2[i]) ELSE $
                    err_mess = [err_mess, 'Not a  valid number !']
              ENDFOR 
           END 

           dummy1 = trim(float(dummy1))
           dummy2 = trim(float(dummy2))


           WIDGET_CONTROL,state.iso_logt_ev,$
                          set_value= arr2str(dummy1, ',',/trim)

           WIDGET_CONTROL,state.iso_logem_ev,$
                          set_value=arr2str(dummy2, ',',/trim)

;         IF min(iso_logt)  LE  0. THEN err_mess = [err_mess, 'Wrong Log T - Negative ? '] 

           
        ENDIF  ELSE BEGIN 

           IF n_elements(demfile) EQ 1 THEN BEGIN 
              dem_name = concat_dir(demdir, demfile)
              IF NOT  file_exist(dem_name) THEN err_mess = [err_mess, 'Check DEM file defined']
           ENDIF ELSE err_mess = [err_mess, 'No DEM file defined']

        ENDELSE 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        
        IF photoexcitation THEN BEGIN 

           dummy = ''
           WIDGET_CONTROL,state.photoexcitation_rphot_ev, get_v=dummy

           IF valid_num(dummy(0)) AND dummy(0) GE   1. THEN BEGIN 
              rphot= float(dummy) & rphot=rphot(0)
           ENDIF ELSE BEGIN 
              err_mess = [err_mess, 'Not a  valid number !']
              WIDGET_CONTROL,state.photoexcitation_rphot_ev, $
                             set_value= strtrim(string(format='(f6.3)', 1.0 ),2)
           END 

           dummy = ''
           WIDGET_CONTROL,state.photoexcitation_radtemp_ev, get_v=dummy

           IF valid_num(dummy(0)) AND dummy(0) GT 0. THEN BEGIN 
              radtemp= float(dummy) & radtemp=radtemp(0)
           ENDIF ELSE BEGIN 
              err_mess = [err_mess, 'Not a  valid number !']
              WIDGET_CONTROL,state.photoexcitation_radtemp_ev,$ 
                             set_value= strtrim(string(format='(f6.1)', 6000.0 ))
           END 

        END

;-------------------------------------------------------
; read the advanced model inputs      
;-------------------------------------------------------

        
; Advanced model options:
        WIDGET_CONTROL,state.ibase_sel,get_uvalue=bob,get_value=fred
        state.calc_ioneq_flag = bob[fred]

        if state.calc_ioneq_flag then begin
           
           advanced_model=1

; check file name first.
           
           dummy=''
           WIDGET_CONTROL,state.output_ioneq_show, get_v=dummy
           dummy=trim(dummy)
           if file_exist(dummy) then begin
              state.output_ioneq_file=''
              WIDGET_CONTROL,state.output_ioneq_show, set_v=state.output_ioneq_file
              err_mess = [err_mess, 'Error - output ioneq file exists - change name! ']
           endif else begin
              state.output_ioneq_file=dummy
              ioneq_name=state.output_ioneq_file
           end

           
           WIDGET_CONTROL,state.min_logt_ev, get_v=dummy
           IF valid_num(dummy(0))  THEN BEGIN 
              state.min_logt= float(dummy[0])
           ENDIF ELSE BEGIN 
              widget_control, state.show_lines , $
                              set_val='min log T  not a  valid number !'
              WIDGET_CONTROL,state.min_logt_ev,$ 
                             set_value= strtrim(string(format='(f3.1)', 4.0 ))
              state.min_logt=4.0
              err_mess = [err_mess, 'Error']
           END

           dummy=''
           WIDGET_CONTROL,state.max_logt_ev, get_v=dummy
           IF valid_num(dummy(0))  THEN BEGIN 
              state.max_logt= float(dummy[0])
           ENDIF ELSE BEGIN 
              widget_control, state.show_lines , $
                              set_val='max log T  not a  valid number !'
              WIDGET_CONTROL,state.max_logt_ev,$ 
                             set_value= strtrim(string(format='(f3.1)', 8.0 ))
              state.max_logt=8.0
              err_mess = [err_mess, 'Error']
           END

           if state.max_logt le state.min_logt then begin
              widget_control, state.show_lines , $
                              set_val='max log T smaller  than min log T ??? '
              err_mess = [err_mess, 'Error']
           end
           if state.min_logt ge state.max_logt then begin
              widget_control, state.show_lines , $
                              set_val='min log T  greater than max log T ??? '
              err_mess = [err_mess, 'Error']
           end

           dummy=''
           WIDGET_CONTROL,state.dlogt_ev, get_v=dummy
           IF valid_num(dummy(0))  THEN BEGIN 
              state.dlogt= float(dummy[0])
              if state.dlogt le 0 or  state.dlogt ge  state.min_logt then begin
                 widget_control, state.show_lines , $
                                 set_val='wrong d log T  ??? '          
                 WIDGET_CONTROL,state.dlogt_ev,$ 
                                set_value= strtrim(string(format='(f3.1)', 0.1 ))
                 state.dlogt=0.1
                 err_mess = [err_mess, 'Error']
              endif 
           ENDIF ELSE BEGIN 
              widget_control, state.show_lines , $
                              set_val='d log T  not a  valid number !'
              WIDGET_CONTROL,state.dlogt_ev,$ 
                             set_value= strtrim(string(format='(f3.1)', 0.1 ))
              state.dlogt=0.1
              err_mess = [err_mess, 'Error']
           END
           
           ioneq_logt= state.min_logt+ state.dlogt*( indgen((state.max_logt-state.min_logt)/state.dlogt)+1)

           ioneq_name= state.output_ioneq_file


           WIDGET_CONTROL,state.ct_base_sel,get_uvalue=bob,get_value=fred
           ct_flag = bob[fred]

           ct=ct_flag

           state.ct_flag=ct_flag

           if ct then begin
                                ;check for atmos[phere file
              WIDGET_CONTROL,state.ct_show, get_value=fred

              fred=fred[0]
              atmosphere_file=concat_dir(concat_dir(concat_dir(concat_dir(!xuvtop,'ancillary_data'),$
                                                               'advanced_models'),'model_atmospheres'), fred)+'.dat'

              if not file_exist(atmosphere_file) then begin 
                 err_mess = [err_mess, 'Atmosphere file for charge transfer not found - EXIT']

                 ct=0
              endif 

           end  


           
           IF  isothermal_flag EQ 1 THEN BEGIN

              if min(iso_logt) lt min(ioneq_logt) or max(iso_logt) gt max(ioneq_logt) then $
                 err_mess = [err_mess, 'Isothermal temperatures chosen are inconsistent with the ion fraction temperatures - EXIT ']

           ENDIF 

           
           
        endif    else begin
           
           advanced_model=0

           IF n_elements(ioneqfile) EQ 1 THEN BEGIN 

              ioneq_name = concat_dir(ioneqdir, ioneqfile)
              IF NOT  file_exist(ioneq_name) THEN err_mess = [err_mess, 'No Ionization Fraction file defined'] 

           ENDIF ELSE   err_mess = [err_mess, 'No Ioniz. Frac. file defined'] 

           IF  isothermal_flag EQ 1 THEN BEGIN
; need to check if the requested log T is within the ranges of the 
; ion fraction temperatures :

              read_ioneq, ioneq_name , logt_i,ioneq_i, ref
              if min(iso_logt) lt min(logt_i) or max(iso_logt) gt max(logt_i) then $
                 err_mess = [err_mess, 'Isothermal temperatures chosen are inconsistent with the ion fraction temperatures in the chosen file - EXIT ']

           ENDIF
           
        end 
        

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        IF n_elements(err_mess) GT 1 THEN BEGIN 

           err_mess = ['Definition for the calculation of the line intensities is incomplete ! ', err_mess]

           dummy = WIDGET_MESSAGE(err_mess, /info)


        ENDIF ELSE BEGIN 



           IF all_ions_yn THEN delvarx, list_ions 

           IF  isothermal_flag EQ 1 THEN delvarx, dem_name,ioneq_logt  ELSE $
              delvarx, iso_logt, iso_logem

           WIDGET_CONTROL, state.plot_rat, GET_VALUE=plot_rat_id ,  sensitive=1 
           wset,plot_rat_id
           wshow,plot_rat_id
           tverase

           widget_control, state.show_lines , $
                           set_val='Please wait, calculating line intensities....'


           WIDGET_CONTROL,Event.top, /hourglass
           WIDGET_CONTROL, SYN_MAIN_base, sensitive=0

           delvarx, spectrum

           IF photoexcitation THEN rphot1 = rphot ELSE delvarx, rphot1

           IF units(1) EQ 'photons' THEN photons = 1 ELSE IF $
              units(1) EQ 'ergs' THEN photons =0 ELSE photons =0

           if state.lookup_available then begin 
              widget_control,state.lookup_ev,get_value=lookup
           endif else begin
              lookup=0
           endelse 
           

           
           ch_synthetic, int_xrange[0], int_xrange[1], output=tran, $
                         err_msg=err_msg, msg=msg, $
                         pressure=pressure, density=density, model_file=model_file, $
                         all=theor_lines,sngl_ion=list_ions, photons=photons, $ 
                         LOGT_ISOTHERMAL=iso_logt, logem_isothermal=iso_logem, $
                         ioneq_name=ioneq_name, dem_name=dem_name, $
                         noprot=noprot, rphot=rphot1, radtemp=radtemp, /progress, lookup=lookup,$
                         ioneq_logt=ioneq_logt,advanced_model=advanced_model,ct=ct,$
                         atmosphere=atmosphere_file ;,he_abund=he_abund

           WIDGET_CONTROL, SYN_MAIN_base, sensitive=1


           IF err_msg[0]  NE '' THEN BEGIN 
              tran = ''
              widget_control, state.show_lines , $
                              set_val=err_msg
              xpopup, err_msg,  xsize=130, title='ERROR !'

           ENDIF  ELSE  BEGIN 

              widget_control, state.show_lines , $
                              set_val='Line intensities calculated. Now calculate spectrum !'

              xrange = int_xrange 
              yrange = [0.,0.]

;reset the units in the spectrum calculation box:
              WIDGET_CONTROL,state.wminw2,  set_value=strtrim(string(format='(f13.6)',xrange(0)),2)
              WIDGET_CONTROL,state.wmaxw2,  set_value=strtrim(string(format='(f13.6)',xrange(1)),2)
              widget_control,state.kev_info, set_val=string(197B)+' [ / keV ]'
              kev_yn=0
              units(0) = 'Angstroms'


              WIDGET_CONTROL,state.xmin_base,  set_value=strtrim(string(format='(f11.4)',xrange(0)),2)
              WIDGET_CONTROL,state.xmax_base,  set_value=strtrim(string(format='(f11.4)',xrange(1)),2)

              WIDGET_CONTROL,state.ymin_base,  set_value=strtrim(string(format='(e9.2)',yrange(0)),2)
              WIDGET_CONTROL,state.ymax_base,  set_value=strtrim(string(format='(e9.2)',yrange(1)),2)

           END 
           
           IF msg[0] NE '' THEN BEGIN 
              widget_control, state.show_lines , $
                              set_val=msg
;
;avoid popup widget:
;            xpopup, msg , xsize=130
           ENDIF 

        END 
     END  

     event.id EQ state.save_str OR event.id EQ state.gsave_str : BEGIN

        err_mess = ''

        IF  datatype(tran, /tname) NE 'STRUCT'  THEN  BEGIN 
           err_mess = [err_mess, 'Line intensities are not in memory !']
           dummy = WIDGET_MESSAGE(err_mess, /info)
        ENDIF ELSE BEGIN 

           if event.id EQ state.gsave_str then begin 
              ff = dialog_pickfile(file='ch_ss_int.genx', tit='Type genx file name ')
;            ff = BIGPICKFILE(file='ch_ss_int.genx', tit='Type genx file name ')

              IF ff NE  '' THEN BEGIN 
                 savegen, file=ff, struct=tran
                 widget_control, state.show_lines , /append, $
                                 set_val='Line intensities  saved in the IDL genx file '+ff
              ENDIF 

           ENDIF  else begin 

              ff = dialog_pickfile(file='ch_ss_int.fits', tit='Type FITS file name ') 
;            ff = BIGPICKFILE(file='ch_ss_int.fits', tit='Type FITS file name ') 

              IF ff NE  '' THEN BEGIN 
                 ch_write_fits, tran, ff
                 widget_control, state.show_lines , /append, $
                                 set_val='Line intensities  saved in the FITS file '+ff
              ENDIF 
           END 
        ENDELSE 
     END   


     event.id EQ state.plot_rat: BEGIN 

;print, event.type
        syn_cursor,event

     END 

     event.id eq state.restore_str OR event.id eq state.grestore_str: BEGIN

        delvarx,tran

        IF event.id eq state.grestore_str THEN begin 

           ff = dialog_pickfile(filter='ch_ss_int*.genx', $ 
                                tit='Select appropiate IDL genx save file to READ') 
           ;; ff = BIGPICKFILE(filter='ch_ss_int*.genx', $ 
           ;;                  tit='Select appropiate IDL genx save file to READ') 

           IF ff NE  '' THEN restgen, file=ff, struct=tran,/quiet

        ENDIF   ELSE begin 

           ff = dialog_pickfile(filter='ch_ss_int*.fits', $
                                tit='Select appropiate FITS file to READ')
           ;; ff = BIGPICKFILE(filter='ch_ss_int*.fits', $
           ;;                  tit='Select appropiate FITS file to READ')

           IF ff NE  '' THEN  BEGIN 
              ch_read_fits, ff, tran,  err_msg=err_msg
              IF err_msg[0]  NE '' THEN xpopup, err_msg,  xsize=130, title='Error'
           ENDIF 
        END     

;check we have the structure.

        IF n_elements(TRAN) NE  0 THEN BEGIN 

;check that we didn't pick up one that has the SPECTRUM:

           result = ch_check_str(TRAN, /int)

           IF result THEN BEGIN 
;we restore TRAN (COMMON)
              restore_line_int

           ENDIF ELSE   widget_control, state.show_lines , $
                                        set_val='Error, file did not have a valid line intensities structure !'

        ENDIF   ELSE widget_control, state.show_lines , $
                                     set_val='Error, no  file found !'
        
     END   


     (event.id EQ  state.restore_sp) OR (event.id EQ  state.grestore_sp): BEGIN

        delvarx,tran, spectrum

        IF event.id eq state.grestore_sp THEN begin 

           ff = dialog_pickfile(filter='ch_ss_sp*.genx', $ 
                                tit='Select appropiate IDL genx save file to READ') 
           ;; ff = BIGPICKFILE(filter='ch_ss_sp*.genx', $ 
           ;;                  tit='Select appropiate IDL genx save file to READ') 

           IF ff NE  '' THEN restgen, file=ff, struct=spectrum,/quiet

        ENDIF   ELSE IF (event.id EQ  state.restore_sp) THEN  BEGIN   

           ff = dialog_pickfile(filter='ch_ss_sp*.fits', $
                                tit='Select appropiate FITS file to READ')
           ;; ff = BIGPICKFILE(filter='ch_ss_sp*.fits', $
           ;;                  tit='Select appropiate FITS file to READ')

           IF ff NE  '' THEN BEGIN 
              ch_read_fits, ff, spectrum,  err_msg=err_msg
              IF err_msg[0]  NE '' THEN xpopup, err_msg,  xsize=130, title='Error'
           ENDIF 
        END     

;check we have the structure.

        IF n_elements(spectrum) NE  0 THEN BEGIN 

;check that we didn't pick up one that has the SPECTRUM:

           result = ch_check_str(spectrum, /sp)

           IF result THEN BEGIN 

;reconstruct the TRAN structure (see MAKE_CHIANTI_SPEC):

              tran = spectrum

;remove the abundance factor:
              tran.lines[*].int = tran.lines[*].int / tran.abund[tran.lines[*].iz-1]

              lines = tran.lines 
              lines = rem_tag(lines, 'peak')

              tran = rem_tag(tran, 'lines')
              tran = add_tag(tran, lines, 'lines')


              IF tag_exist(tran, 'EFFAREA' ) THEN BEGIN 
;get the nearest effective area value to find the 
; approximate new intensity values.

                 FOR i=0L, n_elements(tran.lines) -1 DO BEGIN 

                    dummy = min(abs(tran.lambda-tran.lines[i].wvl), ix)
                    ix = ix[0]
;GDZ: do not divide by 4*!pi
                    IF ix NE -1 THEN tran.lines[i].int = tran.lines[i].int /tran.effarea[ix] ELSE $
                       tran.lines[i].int =0.

                 ENDFOR 
                 tran.int_units ='photons cm-2 sr-1 s-1'

              ENDIF 

              tran = rem_tag(tran, 'LAMBDA')
              tran = rem_tag(tran, 'SPECTRUM')
              tran = rem_tag(tran, 'UNITS')
              tran = rem_tag(tran, 'INSTR_FWHM')
              tran = rem_tag(tran, 'BIN_SIZE' )
              tran = rem_tag(tran, 'ABUND_NAME')
              tran = rem_tag(tran, 'ABUND')
              tran = rem_tag(tran, 'MIN_ABUND')
              tran = rem_tag(tran, 'ABUND_REF')

              IF tag_exist(tran, 'CONTINUUM') THEN tran = rem_tag(tran,'CONTINUUM')
              IF tag_exist(tran, 'FILE_EFFAREA' ) THEN BEGIN 
                 FILE_EFFAREA = tran.FILE_EFFAREA
                 tran = rem_tag(tran,'FILE_EFFAREA' )
              END

              IF tag_exist(tran, 'EFFAREA' ) THEN tran = rem_tag(tran, 'EFFAREA')


              restore_line_int  ;will use TRAN   (COMMON)
              restore_spectrum  ; will use SPECTRUM (COMMON)

              index = WHERE(spectrum.lines[*].peak EQ max(spectrum.lines[*].peak))
              o_strength = spectrum.lines[index[0]].peak/10.
              WIDGET_CONTROL,state.strength_yn,set_value=trim(string(o_strength, format='(e9.2)'))
              plot_syn_spectrum
              widget_control, state.show_lines , $
                              set_val='Select lines with mouse'


           ENDIF ELSE   widget_control, state.show_lines , $
                                        set_val='Error, file did not have a valid spectrum structure !'

        ENDIF   ELSE widget_control, state.show_lines , $
                                     set_val='Error, no  file found !'
        
     END   


     event.id EQ state.calc_sp_butt: BEGIN 

        str=['This is a short HELP. For more details on single buttons please click the buttons where the " - HELP" sign is. ', $
             ' ', $
             'PLEASE READ THE CHIANTI USER GUIDES AND THE DOCUMENTATION IN THE HEADERS OF THE PROCEDURES to understand how to use this software. ', $
             '  ', $
             ' Send comments to Giulio Del Zanna  if the details in the documentation are not clear or if you have problems in running the software.', $
             ' ', $
             'This section of the program calculates, plots a synthetic spectrum and allows ', $
             'the creation of tables of line intensities and various outputs.', $
             'If you have already calculated a CHIANTI SPECTRUM STRUCTURE, you can restore it with the "RESTORE spectrum" button.', $
             'Otherwise, you can calculate the CHIANTI SPECTRUM STRUCTURE by pressing the button "Calculate and plot", which calls the subroutine MAKE_CHIANTI_SPEC.', $
             ' ', $
             ' The CHIANTI SPECTRUM STRUCTURE contains the X and Y values, plus all the details of the parameters used in the calculation ', $
             ' (see the CH_SS header for details), and all the details of the spectral lines.', $
             ' ', $
             'NOTE: a CHIANTI LINE INTENSITY STRUCTURE  MUST be in the program memory before this section of the program is run.', $
             'This can be done either by calculating the line intensities or by restoring a save file with previously calculated values. ', $
             ' ', $
             'Once the CHIANTI LINE INTENSITY STRUCTURE is in the program memory, you have to choose: ', $
             '-) Minimum and maximum wavelengths in Angstroms', $
             '-) the spectrum bin size in Angstroms (disallowed if an Effective area file is used).', $
             '-) the instrumental full-width-half-maximum (FWHM, in Angstroms). Setting this to a non-zero value broadens the lines with a Gaussian of FWHM=FWHM', $
             '-) if you want to add the continua (free-free, free-bound and two-photon (see  HELP button for details)', $
             '-) if you want to plot all lines or remove the "unobserved lines" (if they were orignally calculated - see  HELP button) ', $
             '-) an abundance file and a minimum abundance value (see  HELP button for details). ', $
             '-) (Eff. Area: Yes/No)', $
             '       If you want to fold the spectrum with an effective area.  ', $
             '       If set to Yes, you are requested to choose an input ascii file with two columns,  ', $
             '       the wavelength and the effective area values (cm^2).', $
             '       The wavelenghts in the file (that might not be linear)  are used  to create the spectrum,', $
             '       that is multiplied with the effective area values. The line intensities contributing to each ', $
             '       bin are summed, and  subsequently convolved with a gaussian of FWHM=FWHM, if FWHM is not set = 0. ', $
             '       Note that this option only works well if a sufficient number of bins is given.', $
             '       Also note that to have the correct output units  (counts s-1 bin-1) the appropiately scaled DEM (or EM) values must be provided. ', $
             ' ', $
             'Once the spectrum is displayed below, it is then possible to perform ', $
             'various operations by clicking ant typing on the buttons below. For example: ', $
             ' ', $
             '-) Add labels to the plot and select a cut (Minimum peak intensity -- see  HELP button for details). ', $
             '-) Click on a region of the plot and view in the bottom panel the details of the lines', $
             '-) ZOOM in and out with the buttons, or type in the X,Y ranges (in this case, whenever a new value is typed', $
             '     in the text widget, ENTER must be pressed to make the new value accepted. ', $
             '-) Create a postcript file of the plot, or send directly an hardcopy to the PRINTER. ', $
             '-) Save line details into latex or ascii files. NOTE that ONLY the details of the lines that have labels ', $
             '   shown in the plot are saved. If you want all of them, you have to Unzoom and set the Min. value to 0. ' , $
             ' ', $
             '-) Save the CHIANTI SPECTRUM  in three options: ', $
             '   1- simple ascii file with two columns, X and Y values.',$
             '   2- save (IDL file) the CHIANTI SPECTRUM STRUCTURE (this is done by calling SAVEGEN ). ', $
             '      NOTE: to restore the above file into an IDL stucture STRUCT, use (at the command line):  IDL> restgen,file=filename,struct=STRUCT ', $
             '   3- save the CHIANTI SPECTRUM STRUCTURE  as a FITS binary table.',$
             '      NOTE: to restore the above file into an IDL stucture STRUCT, use (at the command line):  IDL> ch_read_fits, filename, STRUCT ', $
             ' ', $
             'Giulio Del Zanna & Peter Young' ]
        xpopup, str, xsize=130
     END


     event.id EQ state.kev_info: BEGIN
;switch units:
        if kev_yn then begin 
           widget_control,state.kev_info, set_val=string(197B)+' [ / keV ]'
           kev_yn=0
           units(0) = 'Angstroms'

;switch the ranges from keV to Angstroms:
;----------------------------------------
           widget_control,state.wminw2, get_v=val 
           IF valid_num(val(0)) AND val(0) GT 0. THEN $ 
              xrange(1) =12.39854/float(val(0))   ELSE xrange(1) =tran.WVL_LIMITS[1]
           widget_control,state.wmaxw2, get_v=val 
           IF valid_num(val(0)) AND val(0) GT 0. THEN $
              xrange(0) =12.39854/float(val(0))   ELSE xrange(0) =tran.WVL_LIMITS[0]

           WIDGET_CONTROL,state.wminw2,  set_value=strtrim(string(xrange(0), format='(f13.6)'),2)
           WIDGET_CONTROL,state.wmaxw2,  set_value=strtrim(string(xrange(1), format='(f13.6)'),2)

        endif else if not kev_yn then begin 
           widget_control,state.kev_info, set_val='keV [ / '+string(197B)+' ]'
           kev_yn=1
           units(0) ='keV'

;switch the ranges from Angstroms to keV:
;---------------------------------------
           widget_control,state.wminw2, get_v=val 
           IF valid_num(val(0)) AND val(0) GT 0. THEN $ 
              xrange(1) =12.39854/float(val(0)) ELSE xrange(1) =12.39854/(tran.WVL_LIMITS[0]>1e-5)
           widget_control,state.wmaxw2, get_v=val 
           IF valid_num(val(0)) AND val(0) GT 0. THEN $
              xrange(0) =12.39854/float(val(0)) ELSE xrange(0) =12.39854/(tran.WVL_LIMITS[1]>1e-5)

           WIDGET_CONTROL,state.wminw2,  set_value=strtrim(string(xrange(0), format='(f13.6)'),2)
           WIDGET_CONTROL,state.wmaxw2,  set_value=strtrim(string(xrange(1), format='(f13.6)'),2)
        end 
     end 

     event.id EQ state.wminw2 :BEGIN 
        widget_control,state.wminw2, get_v=val 
        IF valid_num(val(0)) AND val(0) GT 0. THEN BEGIN 
           xrange(0) = float(val(0)) 
        ENDIF ELSE $
           WIDGET_CONTROL,state.wminw2,  set_value=strtrim(string(xrange(0), format='(f13.6)'),2)
        delvarx, spectrum
        WIDGET_CONTROL, state.plot_rat, GET_VALUE=plot_rat_id ,  sensitive=1 
        wset,plot_rat_id & wshow,plot_rat_id & tverase
        widget_control, state.show_lines , $
                        set_val='Please calculate the  spectrum !'
     END

     event.id EQ state.wmaxw2 :BEGIN 
        widget_control,state.wmaxw2, get_v=val 
        IF valid_num(val(0)) AND val(0) GT 0. THEN BEGIN 
           xrange(1) = float(val(0)) 
        ENDIF ELSE $
           WIDGET_CONTROL,state.wmaxw2,  set_value=strtrim(string(xrange(1), format='(f13.6)'),2)
        delvarx, spectrum
        WIDGET_CONTROL, state.plot_rat, GET_VALUE=plot_rat_id ,  sensitive=1 
        wset,plot_rat_id & wshow,plot_rat_id & tverase
        widget_control, state.show_lines , $
                        set_val='Please calculate the  spectrum !'
     END

     event.id EQ state.ang_butt: BEGIN
        str=['Bin size  in Angstroms of the spectrum', $
             'to be created. A linear spectrum is assumed.', $
             'In case an effective area file is used, the wavelenghts in the file', $
             '(that might not be linear) are used to create the spectrum, and ', $
             'this Bin size looses any meaning.']
        result = DIALOG_MESSAGE(str,/info)
     END


     event.id EQ  state.ang_read: BEGIN
        WIDGET_CONTROL,state.ang_read,get_value=bob
        ang=float(bob) & ang=ang[0]
        IF NOT valid_num(ang) THEN WIDGET_CONTROL,state.ang_read,set_value=trim(0.1)

        delvarx, spectrum
        WIDGET_CONTROL, state.plot_rat, GET_VALUE=plot_rat_id ,  sensitive=1 
        wset,plot_rat_id & wshow,plot_rat_id & tverase
        widget_control, state.show_lines , $
                        set_val='Please calculate the  spectrum !'

     END

     event.id EQ state.inst_butt: BEGIN
        str=['Instrumental broadening', $
             ' ', $
             'The spectral lines will have a Gaussian profile for which the width ', $
             'is given by a combination of thermal broadening and instrumental ', $
             'broadening. The thermal broadening is determined from the Tmax of ', $
             'the ion.', $
             '', $
             'The user can control the instrumental broadening by specifying a ', $
             'FWHM value in this box. The same value is applied to all spectral ', $
             'lines. A value of zero switches off the instrumental broadening, ', $
             'leaving only the thermal broadening.', $
             ' ', $
             'The units of FWHM are either angstroms or keV, depending on ', $
             'which unit is used for the spectrum.']
        result = DIALOG_MESSAGE(str,/info)
     END

     event.id eq state.inst_read: BEGIN
        WIDGET_CONTROL,state.inst_read,get_value=bob
        inst=float(bob) & inst=inst[0]
        IF NOT valid_num(inst) THEN  WIDGET_CONTROL,state.inst_read,set_value=trim(0.1)

        delvarx, spectrum
        WIDGET_CONTROL, state.plot_rat, GET_VALUE=plot_rat_id ,  sensitive=1 
        wset,plot_rat_id & wshow,plot_rat_id & tverase
        widget_control, state.show_lines , $
                        set_val='Please calculate the  spectrum !'

     END

     event.id EQ state.continua_butt: BEGIN
        str=['Add continua to the binned spectrum (free-free,  free-bound and two-photon)', $
             '',$
             'Please note that the continuum calculation takes some time and you may want ', $
             ' ',  'to define a minimum abundance value or choose a more coarse spectral resolution']
        result = DIALOG_MESSAGE(str,/info)
     END

     event.id EQ state.rem_theor2_butt: BEGIN
        str=['Add lines that have only theoretical energy levels ? ', $
             ' ', $
             'Note that it is possible to  remove these "unobserved lines" only if', $
             'they were originally calculated (i.e. they are present in the IDL structure).', $
             ' ', $
             'WARNING:  the wavelengths of these  "unobserved lines"', $
             'might not be very accurate']
        result = DIALOG_MESSAGE(str,/info)
     END


     event.id EQ state.rem_theor2: BEGIN

        IF n_elements(tran) GT 0 THEN BEGIN 

           WIDGET_CONTROL,state.rem_theor2,get_uvalue=bob,get_value=fred
           theor_lines =  bob[fred]

;add a check if there are any to add to the spectrum:

           IF theor_lines EQ 1 THEN BEGIN 

;flag = transitions.lines(*).flag  ;-1 for the theoretical lines, 0 for the
;observed ones..

              index1 = where(tran.lines[*].flag  EQ -1 AND $
                             tran.lines[*].wvl GE xrange[0] AND tran.lines[*].wvl LE xrange[1], nlines1)

              index2 = where(tran.lines[*].wvl GE xrange[0] AND $
                             tran.lines[*].wvl LE xrange[1], nlines2)

              IF nlines1 EQ 0 THEN BEGIN
                 theor_lines = 0
                 widget_control, state.show_lines , /append, $
                                 set_val='No "unobserved" lines present.....'
                 WIDGET_CONTROL,state.rem_theor2,  set_value=theor_lines
              ENDIF ELSE BEGIN 

                 delvarx, spectrum
                 WIDGET_CONTROL, state.plot_rat, GET_VALUE=plot_rat_id ,  sensitive=1 
                 wset,plot_rat_id & wshow,plot_rat_id
                 tverase

                 str=['There are '+trim(nlines1)+'  "unobserved" lines in the '+$
                      trim(xrange[0])+' - '+trim(xrange[1])+' range ', $
                      ' over a total of '+trim(nlines2)+' lines in the same wavelength range.', $
                      ' ', $
                      'WARNING:  the wavelengths of these  "unobserved" lines', $
                      'might not be very accurate']

                 result = DIALOG_MESSAGE(str,/info)

                 widget_control, state.show_lines , $
                                 set_val='Please calculate the  spectrum !'
              ENDELSE 
           ENDIF                ;theor_lines set to 1.
        END 
     END   

     event.id EQ state.continua: BEGIN

        WIDGET_CONTROL,state.continua, get_value=fred
        cont_ind= fix(fred)

        delvarx, spectrum
        WIDGET_CONTROL, state.plot_rat, GET_VALUE=plot_rat_id ,  sensitive=1 
        wset,plot_rat_id & wshow,plot_rat_id
        tverase

        widget_control, state.show_lines , $
                        set_val='Please calculate the  spectrum !'
        
     END

     event.id EQ state.ab_button: BEGIN
        str=[' ', $
             ' ', $
             'WARNING: only the lines of those elements which have an abundance defined', $
             'in the file  will be selected to create the spectrum.', $
             ' ', $
             ' ']
        result = DIALOG_MESSAGE(str,/info)
     END

     event.id EQ state.ab_pdmenu: BEGIN

        WIDGET_CONTROL,state.ab_pdmenu,  get_uvalue=bob

        abfile = ''
        abund_name = ''


        CASE  bob OF 

           '0': BEGIN 
              WIDGET_CONTROL,state.ab_show,set_value=''
           END 

           '1':BEGIN
              abund_path=concat_dir(!xuvtop,'abundance')
              abund_name = dialog_pickfile(filter='*abund', $
                                           title='Select appropiate ABUNDANCE file to READ', $
                                           path=abund_path)
;            abund_name = BIGPICKFILE(filter='*abund', tit='Select appropiate ABUNDANCE file to READ')
              IF abund_name NE '' THEN BEGIN
                 break_file, abund_name, disk, abdir, abfile, ext
                 abdir = concat_dir(disk, abdir)
                 abfile = abfile+ext
                 WIDGET_CONTROL,state.ab_show,set_value=abfile
              END 
           END 
           ELSE:BEGIN 
              abund_name= bob
              break_file, abund_name, disk, abdir, abfile, ext
              abdir = concat_dir(disk, abdir)
              abfile = abfile+ext

              WIDGET_CONTROL,state.ab_show,set_value=abfile

           END 
        ENDCASE 

        IF abund_name NE '' THEN BEGIN 

           read_abund,abund_name,abund,abund_ref

           index = where(abund GT 0.0)
           min_abund = min(abund(index))

           WIDGET_CONTROL, state.min_abund_ev,$
                           set_value=strtrim(string(MIN_ABUND, format='(e9.2)'),2)

           widget_control, state.show_lines,set_val= ''

           FOR  i=0,n_elements(abund_ref)-1 DO $
              widget_control, state.show_lines,/append , $
                              set_val=abund_ref[i]

; print the abundance values:

           IF index(0) GE 0 THEN BEGIN 

              el = elements(index)

              widget_control, state.show_lines , $
                              set_val='Abundance values: ', /append

              FOR  i=0,n_elements(el)-1 DO $
                 widget_control, state.show_lines,/append , $
                                 set_val=strpad(el(i), 3, /after)+'  '+ trim(abund(index[i]))

           ENDIF 

           delvarx, spectrum
           WIDGET_CONTROL, state.plot_rat, GET_VALUE=plot_rat_id ,  sensitive=1 
           wset,plot_rat_id & wshow,plot_rat_id
           tverase

           widget_control, state.show_lines , $
                           set_val='Please calculate the  spectrum !',/append
        END 

     END


     event.id EQ state.min_abund_button: BEGIN 
        str=['WARNING: If set not null, only the lines of those elements which have an abundance greater than', $
             'the value set here are selected and will be used to create the spectrum.', $
             ' ', $
             'Also, the continuum is calculated only for those elements which  have an abundance greater than', $
             'the value set. This can significantly speed up the calculations.', $
             ' ', $
             'By default, the minimum value in the selected abundance file is used. ', $
             ' ', $
             'To have an idea of what minimum abundance should be set, select an abundance file', $
             ' and read the bottom text message window.  ', $
             ' ', $
             ' As a reference, the abundances of Allen (1973) give:', $
             'abundance (H)  = 1. ', $
             'abundance (He) = 0.085 ', $
             'abundance (C)  = 3.3e-4 ', $
             'abundance (O)  = 6.6e-4 ', $
             'abundance (Si) = 3.3e-5 ', $
             'abundance (Fe) = 3.9e-5 ']

        result = DIALOG_MESSAGE(str,/info)
     END

     event.id EQ state.min_abund_ev: BEGIN 

        mina = min_abund
        widget_control,state.min_abund_ev, get_v=mina
        min_abund = float(mina(0))
        IF valid_num(min_abund) THEN BEGIN 
           IF  min_abund  LT 0. THEN BEGIN 
              min_abund =0.0
              WIDGET_CONTROL,state.min_abund_ev,set_value=trim(0.0) 
           ENDIF ELSE BEGIN 
              widget_control, state.show_lines , $
                              set_val=' Will be using minimum abundance: '+trim(min_abund)
           END 
        ENDIF  ELSE BEGIN 
           min_abund = 0.
           WIDGET_CONTROL,state.min_abund_ev,set_value=trim(min_abund)
        END 

        delvarx, spectrum
        WIDGET_CONTROL, state.plot_rat, GET_VALUE=plot_rat_id ,  sensitive=1 
        wset,plot_rat_id & wshow,plot_rat_id
        tverase
        widget_control, state.show_lines , $
                        set_val='Please calculate the  spectrum !'

     END 

     event.id EQ state.fold_info: BEGIN

        IF fold_yn THEN  BEGIN 

           widget_control,state.fold_info, set_val='Eff. Area: NO'
           widget_control,state.fold_show, set_val=' -- '
           WIDGET_CONTROL,state.ang_read,set_value=trim(0.1)
           fold_yn = 0
           delvarx, lambda, effarea, file_effarea

           widget_control, state.unit_info2, set_val='Units: ERGS'
           units(1) = 'ergs'

        endif else BEGIN

           widget_control,state.fold_info, set_val='Eff. Area: YES'
           fold_yn = 1
           widget_control, state.unit_info2, set_val='Units: PHOTONS'
           units(1) ='photons'

           widget_control,state.wminw2, get_v=val 
           IF valid_num(val(0)) AND val(0) GT 0. THEN BEGIN 
              xrange(0) = float(val(0)) 
           ENDIF ELSE $
              WIDGET_CONTROL,state.wminw2,  set_value=strtrim(string(xrange(0), format='(f13.6)'),2)

           widget_control,state.wmaxw2, get_v=val 
           IF valid_num(val(0)) AND val(0) GT 0. THEN BEGIN 
              xrange(1) = float(val(0)) 
           ENDIF ELSE $
              WIDGET_CONTROL,state.wmaxw2,  set_value=strtrim(string(xrange(1), format='(f13.6)'),2)

;read the effective area 

           dir = concat_dir(concat_dir(!xuvtop,'ancillary_data'), 'instrument_responses')

           IF DIR_EXIST(dir) THEN path = dir  ELSE BEGIN 
              cd, current=dir
              path = dir
           END 

           file_effarea =  dialog_pickfile(title=$
                                           'Select Effective Area File (TWO COLUMNS - wavelengths in Angstroms and cm^2 values )', $
                                           path= path, filter='*area' )
           ;; file_effarea =  bigpickfile(title=$
           ;;                             'Select Effective Area File (TWO COLUMNS - wavelengths in Angstroms and cm^2 values )', $
           ;;                             path= path, filter='*area' )


           IF file_exist(file_effarea) THEN BEGIN 

              data = read_ascii (file_effarea)

; lambda must be sorted in increasing order.

              lambda_effarea = reform( data.field1(0, *))
              lambda = lambda_effarea
              effarea = reform( data.field1(1, *))

; lambda must be sorted in increasing order.
              lambda=lambda[sort(lambda)]
              effarea =effarea(sort(lambda))

; convert to energy.
              if kev_yn then begin 
                 lambda=12.39854/lambda
                 lambda=lambda[sort(lambda)]
                 effarea =effarea(sort(lambda))
              end 

;check that we have at least say 10 points of overlap  in the
; wavelength range:

              in  = where(lambda GE xrange[0] AND lambda LE  xrange[1], ng)

              IF ng LT  10 THEN BEGIN 

                 result = DIALOG_MESSAGE('less than 10 points in the spectrum -- EXIT', /info)

                 widget_control,state.fold_info, set_val='Eff. Area: NO'
                 widget_control,state.fold_show, set_val=' -- '
                 WIDGET_CONTROL,state.ang_read,set_value=trim(0.1)
                 fold_yn = 0
                 delvarx, lambda, effarea, file_effarea
                 widget_control, state.unit_info2, set_val='Units: ERGS'
                 units(1) = 'ergs'

              ENDIF ELSE BEGIN  

;plot the values:

                 window,0,xs=600,ys=600
                 
                 break_file,file_effarea, disk, dir, f, ext
                 f = f+ext

                 IF  kev_yn THEN  xtit='keV' ELSE xtit = string(197B) 

                 plot, lambda[in], effarea[in], psym=-4,syms=1.2,$
                       tit='Effective area file: '+f,  xtit=xtit  ,$
                       ytit='cm!S!E-2 !N ',chars=1.2, /yno

                 wait, 2
                 wshow,!d.window,0

                 WIDGET_CONTROL,state.ang_read,set_value=' -- '
                 widget_control,state.fold_show, set_val=f+ext


              END
           ENDIF  ELSE BEGIN 
              widget_control,state.fold_info, set_val='Eff. Area: NO'
              fold_yn = 0
              delvarx, lambda, effarea, file_effarea

              widget_control,state.fold_show, set_val=' -- '
              WIDGET_CONTROL,state.ang_read,set_value=trim(0.1)
              widget_control, state.unit_info2, set_val='Units: ERGS'
              units(1) = 'ergs'

           ENDELSE 

        END  

        delvarx, spectrum
        WIDGET_CONTROL, state.plot_rat, GET_VALUE=plot_rat_id ,  sensitive=1 
        wset,plot_rat_id & wshow,plot_rat_id
        tverase
        widget_control, state.show_lines , $
                        set_val='Please calculate the  spectrum !'

     END    

     event.id EQ state.fold_button: BEGIN
        str=['  ', $
             'If you request an effective area file, you will have the option to select one of the  ', $
             'files stored in the !xuvtop/ancillary_data/instrument_responses/ directory, that contains   ', $
             'the effective areas of many instruments. These files are simple ascii files with two columns,   ', $
             'wavelengths in Ansgtroms and effective areas (cm-2).   ', $
             'These files have been created from the analogous files from PIMMS version 3.2c:  ', $
             'http://heasarc.gsfc.nasa.gov/docs/software/tools/pimms.html   ', $
             '   ', $
             'User-defined files can be used, but must have the same format.   ', $
             '   ', $
             'In case an effective area file is used, the wavelengths in the file', $
             '(that might not be linear) are used to create the spectrum. ', $
             '   ', $
             'WARNING: this option only works well if a sufficient number of bins is  ', $
             'given. The line intensities contributing to each bin are summed, and  ', $
             'subsequently convolved with a gaussian of full-width-half-maximum FWHM, ', $
             ' if FWHM is not set = 0. ', $
             '   ', $
             'WARNING: the convolution might not work if a small number of  bins is defined. ' ]

        result = DIALOG_MESSAGE(str,/info)
     END



     event.id EQ state.unit_info2: BEGIN

;
;  Change UNITS
;
;first remove any spectrum structure still around:
        delvarx, spectrum

        if units(1) EQ 'photons'  then begin
           widget_control, state.unit_info2, set_val='Units: ERGS'
           units(1) = 'ergs'

           WIDGET_CONTROL, state.plot_rat, GET_VALUE=plot_rat_id ,  sensitive=1 
           wset,plot_rat_id & wshow,plot_rat_id & tverase
           widget_control, state.show_lines , $
                           set_val='Please calculate the  spectrum !'

        endif else begin
           widget_control, state.unit_info2, set_val='Units: PHOTONS'
           units(1) ='photons'

           WIDGET_CONTROL, state.plot_rat, GET_VALUE=plot_rat_id ,  sensitive=1 
           wset,plot_rat_id & wshow,plot_rat_id
           tverase
           widget_control, state.show_lines , $
                           set_val='Please calculate the  spectrum !'

        endelse
     END 



     event.id EQ state.calc_sp : BEGIN

        WIDGET_CONTROL, SYN_MAIN_base, sensitive=0

;check that we have what we need:
;-------------------------------

        err_mess = ''


        IF NOT ch_check_str(TRAN, /int)  THEN $
           err_mess = [err_mess, 'Line intensities are not in memory !', $
                       'Either calculate them or restore a SAVE file! ']

;read the xrange:

        if kev_yn then valid_min_xrange=12.39854/int_xrange(1) else valid_min_xrange=int_xrange(0)
        if kev_yn then valid_max_xrange=12.39854/int_xrange(0) else valid_max_xrange=int_xrange(1)


; avoid  checking on the limits
        
        val = 0.
        widget_control,state.wminw2, get_v=val 
;AND float(val(0)) GE valid_min_xrange
        IF valid_num(val(0)) AND float(val(0)) GT 0.   THEN BEGIN 
           xrange(0) = float(val(0)) 
        ENDIF ELSE BEGIN 
           xrange(0) = valid_min_xrange
           WIDGET_CONTROL,state.wminw2,  set_value=strtrim(string(xrange(0), format='(f13.6)'),2)
           err_mess = [err_mess, 'No correct Min. X value defined -- resetting']
        END 

        val = 0.
        widget_control,state.wmaxw2, get_v=val 
;AND float(val(0)) LE valid_max_xrange
        IF valid_num(val(0)) AND ( float(val(0)) GT 0. AND float(val(0)) gt xrange(0) )  THEN BEGIN 
           xrange(1) = float(val(0)) 
        ENDIF ELSE BEGIN 
           xrange(1) = valid_max_xrange
           WIDGET_CONTROL,state.wmaxw2,  set_value=strtrim(string(xrange(1), format='(f13.6)'),2)
           err_mess = [err_mess, 'No correct Max. X value defined -- resetting']
        END 



        IF fold_yn EQ 0 THEN BEGIN 

           WIDGET_CONTROL,state.ang_read,get_value=bob
           ang=float(bob) & ang=ang[0]

           IF (NOT valid_num(ang))  OR (valid_num(ang) AND  (ang LE 0.)) $ 
              OR (valid_num(ang) AND ang gt (xrange(1)- xrange(0))/10.) THEN BEGIN 
              ang=(xrange(1)- xrange(0))/10.
              err_mess = [err_mess, 'Bin size not valid - resetting to 1/10 the range'] 
              WIDGET_CONTROL,state.ang_read,set_value=trim( ang )
           END 

           delvarx, lambda
        END 

        WIDGET_CONTROL,state.inst_read,get_value=bob
        inst=float(bob) & inst=inst[0]

        IF fold_yn EQ 0 THEN BEGIN 

           IF (NOT valid_num(inst)) OR (valid_num(inst) AND  (inst LT  0.)) or $
;             (valid_num(inst) AND inst lt ang) or $
              (valid_num(inst) AND inst gt (xrange(1)- xrange(0))/5.)   THEN BEGIN 
              inst=ang
              err_mess = [err_mess, 'FWHM not valid - reset to bin size'] 
              WIDGET_CONTROL,state.inst_read,set_value=trim(ang)
           ENDIF ELSE WIDGET_CONTROL,state.inst_read,set_value=trim(inst)

        endif  else begin 

;check the FWHM:

           data = read_ascii (file_effarea)

; lambda must be sorted in increasing order.

           lambda_effarea = reform( data.field1(0, *))
           lambda = lambda_effarea
; lambda must be sorted in increasing order.
           lambda=lambda[sort(lambda)]
; convert to energy.
           if kev_yn then begin 
              lambda=12.39854/lambda
              lambda=lambda[sort(lambda)]
           END 

           binsize = fltarr( n_elements(LAMBDA))
           binsize[0:n_elements(LAMBDA)-2]= $
              lambda[1+indgen(n_elements(LAMBDA)-1)]-lambda[indgen(n_elements(LAMBDA)-1)]
           binsize(n_elements(LAMBDA)-1) = binsize(n_elements(LAMBDA)-2)

           IF (NOT valid_num(inst))  OR (valid_num(inst) AND  (inst LT min(binsize) AND inst GT 0. )) THEN BEGIN 
              err_mess = [err_mess, 'FWHM less than bin-size - resetting '] 
              WIDGET_CONTROL,state.inst_read,set_value=trim(min(binsize))
           END 
        end 

;read continuum 

        WIDGET_CONTROL,state.continua, get_value=fred
        cont_ind= fix(fred)

;read the keyword all lines:

        WIDGET_CONTROL,state.rem_theor2,get_uvalue=bob,get_value=fred
        theor_lines = bob[fred]



        IF n_elements(abfile)  EQ 1 THEN BEGIN 
;check 
           dummy = ''
           WIDGET_CONTROL,state.ab_show, get_value=dummy
           dummy = dummy(0)
           IF dummy NE abfile THEN err_mess = [err_mess, 'Check ABUNDANCE file defined'] $
           ELSE BEGIN 
              abund_name = concat_dir(abdir, abfile)
              
              IF NOT  file_exist(abund_name) THEN  err_mess = [err_mess, 'Check ABUNDANCE file defined']
           END  
        ENDIF  ELSE err_mess = [err_mess, 'No ABUNDANCE file defined']


        mina = ''
        widget_control,state.min_abund_ev, get_v=mina
        min_abund = float(mina(0))
        IF valid_num(min_abund) THEN BEGIN 
           IF min_abund  LT 0. THEN BEGIN 
              min_abund =0.0
              err_mess = [err_mess, 'Min. Abund. value not valid. Reset to 0.']
              WIDGET_CONTROL,state.min_abund_ev,set_value=trim(0.0) 
           ENDIF 
        ENDIF ELSE BEGIN 
           min_abund = 0.
           err_mess = [err_mess, 'Min. Abund. value not valid. Reset to 0.']
           WIDGET_CONTROL, state.min_abund_ev,set_value=strtrim(string(MIN_ABUND, format='(e9.2)'),2)
        END 


        IF n_elements(err_mess) EQ 1 THEN BEGIN 

           read_abund,abund_name,abund,abund_ref

           line_abunds = abund[tran.lines.iz-1]

           index = where(line_abunds[*] GT min_abund , nl)
           IF nl EQ 0 THEN err_mess = [err_mess, 'No lines in the spectrum for the defined abundances']

           index = where(abund GT min_abund, nl)
           IF nl GT 0 THEN BEGIN 

              dummy = where_arr(tran.lines.iz-1, index , nn)
              IF nn NE n_elements(tran.lines.iz-1) THEN BEGIN 

                 elements=['H','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg',$
                           'Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr',$
                           'Mn','Fe','Co','Ni', 'Cu','Zn']

                 dummy2 = where_arr(tran.lines.iz-1, index ,nnn, /notequal)
                 missing = tran.lines[dummy2].iz-1
                 missing = missing(rem_dup(missing))

                 result = dialog_message(['WARNING: the following elements are present in the CHIANTI ', $
                                          'line intensity structure, but will *not* be included in the spectrum:', $
                                          elements[missing], '', 'If you want them to be included you should change ', $
                                          'the element abundances or the minimum abundance settings', ' ', 'CONTINUE ?'], /QUESTION)

                 IF result EQ 'No' THEN err_mess = [err_mess, 'Abundance settings to be changed']

              END  



           ENDIF 

        END 

        IF n_elements(err_mess) GT 1 THEN BEGIN 
           err_mess = ['Definition of the parameters for the calculation of the spectrum is wrong/incomplete: ', ' ', err_mess]
           dummy = WIDGET_MESSAGE(err_mess, /info)

           WIDGET_CONTROL, SYN_MAIN_base, sensitive=1

        ENDIF ELSE BEGIN 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;???  for now check the widget of the switch phot/ergs

           dummy = ''
           widget_control, state.unit_info2, get_val=dummy
           IF dummy EQ 'Units: ERGS' THEN units(1) = 'ergs' ELSE IF $
              dummy EQ 'Units: PHOTONS' THEN units(1) = 'photons'

;         dummy = ''
;         widget_control, state.unit_info3, get_val=dummy
;         IF dummy EQ 'Units: '+angstrom THEN units(0) = 'Angstroms' ELSE IF $
;           dummy EQ 'Units: keV' THEN units(0) = 'keV' 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;         dummy = ''
;         widget_control, state.proton_info, get_val=dummy
;         IF dummy EQ 'Protons: YES' THEN noprot =0 ELSE IF $
;           dummy EQ 'Protons: NO' THEN noprot =1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;reset no ZOOM:

           sty = 0

;do NOT reset the entire x-range available in the Structure:

           WIDGET_CONTROL,state.xmin_base,  set_value=strtrim(string(format='(f11.4)',xrange(0)),2)
           WIDGET_CONTROL,state.xmax_base,  set_value=strtrim(string(format='(f11.4)',xrange(1)),2)


;read the  abundance 

;reset to 

           WIDGET_CONTROL, state.plot_rat, GET_VALUE=plot_rat_id ,  sensitive=1 
           wset,plot_rat_id
           wshow,plot_rat_id
           tverase

           widget_control, state.show_lines , $
                           set_val='Please wait, calculating spectrum......'

           WIDGET_CONTROL,Event.top, /hourglass
           WIDGET_CONTROL, SYN_MAIN_base, sensitive=0

           calc_syn_spectrum, err_msg=err_msg

           WIDGET_CONTROL, SYN_MAIN_base, sensitive=1

           IF n_elements(err_msg) GT 0 THEN BEGIN 

              dummy = WIDGET_MESSAGE(err_msg, /info)

              WIDGET_CONTROL, state.plot_rat, GET_VALUE=plot_rat_id ,  sensitive=1 
              wset,plot_rat_id &  wshow,plot_rat_id &  tverase

              widget_control, state.show_lines , $
                              set_val='Please calculate the  spectrum !'


           ENDIF ELSE BEGIN 

;reset X,Y ranges:

              yrange = [0.,max(spectrum.spectrum)*1.2]

              WIDGET_CONTROL,state.ymin_base,  set_value=strtrim(string(format='(e9.2)',yrange(0)),2)
              WIDGET_CONTROL,state.ymax_base,  set_value=strtrim(string(format='(e9.2)',yrange(1)),2)

              o_lines = 1
              WIDGET_CONTROL,state.oplot_lines,set_value=0

;read oplot lines 
;            WIDGET_CONTROL,state.oplot_lines,get_uvalue=bob,get_value=fred
;            o_lines = bob[fred]

              index = WHERE(spectrum.lines[*].peak EQ max(spectrum.lines[*].peak))
              o_strength = spectrum.lines[index[0]].peak/10.
              WIDGET_CONTROL,state.strength_yn,set_value=trim(string(o_strength, format='(e9.2)'))

              plot_syn_spectrum

              widget_control, state.show_lines , $
                              set_val='Select lines with mouse'

           ENDELSE  
        END  
     END 


     event.id EQ state.oplot_butt: BEGIN
        str=['Overplot lines', $
             '', $
             'Setting this to yes plots a vertical line for each spectral', $
             'line in the spectrum, and also writes a label above the ', $
             'strongest lines indicating the ion from which the line', $
             'arises.']
        result = DIALOG_MESSAGE(str,/info)
     END


     event.id EQ state.oplot_lines: BEGIN
        WIDGET_CONTROL,state.oplot_lines,get_uvalue=bob,get_value=fred
        o_lines = bob[fred]
        IF n_elements(spectrum) NE  0 THEN BEGIN 
           IF o_lines EQ 1 THEN oplot_lines ELSE plot_syn_spectrum
        END 

     END


     event.id EQ state.strength_butt: BEGIN
        str=['Only list strong lines.', $
             '', $
             'When you click on a line in the spectrum in order to get the ', $
             'line ID, you will often receive a list of lines which include ', $
             'very weak lines. Only lines which have an intensity > the value ', $
             'set here will be listed and, if requested,  labelled and selected', $
             'for inclusion in the output. ', $
             'Setting the value=0.  will result ', $
             'in all lines being listed.']
        result = DIALOG_MESSAGE(str,/info)
     END

     event.id EQ state.strength_yn: BEGIN
        WIDGET_CONTROL,state.strength_yn,get_uvalue=aa,get_value=bb
        o_strength = float(bb) &  o_strength=o_strength[0]
        IF o_lines EQ 1 THEN plot_syn_spectrum
     END


     event.id EQ state.xmin_base: BEGIN
        IF n_elements(spectrum) NE  0 THEN BEGIN 
           dummy = xrange(0)
           WIDGET_CONTROL,state.xmin_base, get_v=dummy
           xrange(0) = float(dummy(0))
           sty = 1 
           IF xrange[0] LT MIN(lambda) THEN BEGIN 
              xrange[0] = MIN(lambda)
              WIDGET_CONTROL,state.xmin_base,  set_value=strtrim(string(format='(f11.4)',xrange(0)),2)
           END
           plot_syn_spectrum
        END 
     END 

     event.id EQ state.xmax_base: BEGIN
        IF n_elements(spectrum) NE  0 THEN BEGIN 
           dummy = xrange(1)
           WIDGET_CONTROL,state.xmax_base, get_v=dummy
           xrange(1) = float(dummy(0))
           sty = 1 
           IF xrange[1] GT MAX(lambda) THEN BEGIN 
              xrange[1] = MAX(lambda)
              WIDGET_CONTROL,state.xmax_base,  set_value=strtrim(string(format='(f11.4)',xrange(1)),2)
           END 
           plot_syn_spectrum
        END 
     END 

     event.id EQ state.ymin_base: BEGIN
        dummy = yrange(0)
        WIDGET_CONTROL,state.ymin_base, get_v=dummy
        yrange(0) = float(dummy(0))
        sty = 1 
        plot_syn_spectrum
     END 

     event.id EQ state.ymax_base: BEGIN
        dummy = yrange(1)
        WIDGET_CONTROL,state.ymax_base, get_v=dummy
        yrange(1) = float(dummy(0))
        sty = 1 
        plot_syn_spectrum
     END 

     event.id eq state.log_lin_info: BEGIN

        IF log EQ 0 THEN BEGIN 
           WIDGET_CONTROL,state.log_lin_info,set_value='Log [/Lin]' 
           log=1
        ENDIF ELSE IF log EQ 1 THEN BEGIN 
           WIDGET_CONTROL,state.log_lin_info,set_value='Lin [/Log]'
           log = 0
                                ;
                                ; PRY, 19-Nov-2020: the following resets yrange[0] to 0 when
                                ; switching from the log plot.
                                ;
           yrange[0]=0.
           widget_control,state.ymin_base,set_value=string(format='(e9.2)',yrange[0])
        END 
        plot_syn_spectrum
     END 


     event.id eq state.extras: BEGIN

        IF n_elements(spectrum) NE  0 THEN BEGIN 
           CASE event.value OF
              
              0: BEGIN

                 text=['Click-and-hold a mouse button on the plot to select a', $
                       'region for zooming.']
                 WIDGET_CONTROL,state.show_lines,set_value=text

                 box = ch_drawbox(/data, Color=0)

                 xrange = [box[0],box[2]]
                 xrange = xrange[sort(xrange)]
                 yrange = [box[1],box[3]]
                 yrange = yrange[sort(yrange)]

                 sty = 1

;            IF ABS(yrange[0]-0.) LE (!y.crange[1]-!y.crange[0])/20. THEN yrange[0] = 0.

                 IF xrange[0] LT MIN(lambda) THEN xrange[0] = MIN(lambda)
                 IF xrange[1] GT MAX(lambda) THEN xrange[1] = MAX(lambda)

                 WIDGET_CONTROL,state.xmin_base,  set_value=strtrim(string(format='(f11.4)',xrange(0)),2)
                 WIDGET_CONTROL,state.xmax_base,  set_value=strtrim(string(format='(f11.4)',xrange(1)),2)
                 WIDGET_CONTROL,state.ymin_base,  set_value=strtrim(string(format='(e9.2)',yrange(0)),2)
                 WIDGET_CONTROL,state.ymax_base,  set_value=strtrim(string(format='(e9.2)',yrange(1)),2)

                 plot_syn_spectrum
                 WIDGET_CONTROL,state.show_lines,set_value=''
              END

              1: BEGIN
                 xrange = [ MIN(spectrum.lambda),MAX(spectrum.lambda)]
                 yrange = [0.,max(spectrum.spectrum)*1.2]

                 WIDGET_CONTROL,state.xmin_base,  set_value=strtrim(string(format='(f11.4)',xrange(0)),2)
                 WIDGET_CONTROL,state.xmax_base,  set_value=strtrim(string(format='(f11.4)',xrange(1)),2)
                 WIDGET_CONTROL,state.ymin_base,  set_value=strtrim(string(format='(e9.2)',yrange(0)),2)
                 WIDGET_CONTROL,state.ymax_base,  set_value=strtrim(string(format='(e9.2)',yrange(1)),2)

                 sty = 0
                 plot_syn_spectrum
              END 
           ENDCASE   
        END  
     END 

     event.id eq state.extras1: BEGIN

        IF n_elements(spectrum) NE  0 THEN BEGIN 

           CASE event.value OF

              0: BEGIN

                 CD,CURRENT=curr_dir
                 ff = dialog_pickfile(file='idl.ps', tit='Type postscript file name ')
;               ff = BIGPICKFILE(file='idl.ps', tit='Type postscript file name ')
                 IF ff EQ '' THEN ff = concat_dir(curr_dir, 'idl.ps')

                 ps, ff, /landscape
                 plot_syn_spectrum
                 get_utc,utc
                 utc = anytim2cal(utc)
                 xyouts, 0.05,0.0,'Printed: '+utc,/norm,$
                         chars=0.7  
                 psclose
                 result=WIDGET_MESSAGE('Plot sent to the postscript file: '+ ff,$
                                       /info)

              END  

              1: BEGIN 

                 IF  n_elements(spectrum.spectrum) GT 0 THEN BEGIN 
                    ps 
                    plot_syn_spectrum
                    get_utc,utc
                    utc = anytim2cal(utc)
                    xyouts, 0.05,0.0,'Printed: '+utc,/norm,$
                            chars=0.7  

                    prin = getenv('PRINTER')
                    IF prin EQ '' THEN $
                       xsel_printer,prin,group=SYN_MAIN_base,$
                                    instruct='Select printer for hardcopy'
                    psplot,queu=prin
                    bell
                    result=WIDGET_MESSAGE('Plot sent to printer ' , /info)
                    widget_control, state.show_lines, /append,$
                                    set_val='Plot sent to printer '+prin
                 ENDIF ELSE BEGIN 
                    bell
                    widget_control, state.show_lines, /append,$
                                    set_val='Still nothing to plot!'

                 END 
                 
              END 

           END  
        END  
     END 

     event.id eq state.extras2: BEGIN

        IF n_elements(spectrum) NE  0 THEN BEGIN 

           CASE event.value OF

              0:BEGIN 

;create latex output

                 o_lines = 1
                 WIDGET_CONTROL,state.oplot_lines,  set_value=0
                 widget_control, state.show_lines , $
                                 set_val='The  details of the lines shown will be  saved in the latex file '

                 plot_syn_spectrum

                 good_pix = where(spectrum.lambda GE xrange[0] AND spectrum.lambda LE  xrange[1], ng1)

                 IF ng1 GT 0 THEN BEGIN 

                    good_lines = where((spectrum.lines[*].peak GE o_strength) , ng2 )

                    IF ng2 GT 0 THEN BEGIN 

                       minI = min(spectrum.lines[good_lines].int) > 0 

                       outname = dialog_pickfile(file='ch_ss.tex', tit='Type latex file name ') ;GROUP=
;                     outname = BIGPICKFILE(file='ch_ss.tex', tit='Type latex file name ') ;GROUP=

                       IF outname NE  '' THEN BEGIN 
                          
;                     IF units(1) EQ 'photons' THEN photons = 1 ELSE IF $
;                       units(1) EQ 'ergs' THEN photons =0 ELSE photons =0
;                     IF units(0) EQ 'Angstroms' THEN kev = 0 ELSE $ 
;                       IF units(0) EQ 'keV' THEN kev = 1 ELSE kev = 0

                          ch_line_list, spectrum, outname ,/spectrum,  /latex, $
                                        minI=minI, wmin=xrange[0],wmax=xrange[1],all=theor_lines

                          widget_control, state.show_lines , /append, $
                                          set_val='Line details saved in latex file: '+outname
                          widget_control, state.show_lines , /append, $
                                          set_val='Now latex three times the file: '+outname

                          str=['Line details saved in latex file: '+outname, $
                               'Now latex three times the file: '+outname]
                          result = DIALOG_MESSAGE(str,/info)

                       ENDIF   

                    ENDIF   ELSE  BEGIN 
;bell
                       widget_control, state.show_lines , /append, $
                                       set_val='No lines  to save !!'
                    ENDELSE  
                 ENDIF 
              END 
              
              1: BEGIN

;create ascii file

                 o_lines = 1
                 WIDGET_CONTROL,state.oplot_lines,  set_value=0
                 widget_control, state.show_lines , $
                                 set_val='The  details of the lines shown will be  saved in the ascii file '

                 plot_syn_spectrum

                 good_pix = where(spectrum.lambda GE xrange[0] AND spectrum.lambda LE  xrange[1], ng1)

                 IF ng1 GT 0 THEN BEGIN 

                    good_lines = where((spectrum.lines[*].peak GE o_strength) , ng2 )

                    IF ng2 GT 0 THEN BEGIN 

                       minI = min(spectrum.lines[good_lines].int) > 0 


                       ff = dialog_pickfile(file='ch_ss.txt', tit='Type ascii file name ') ;GROUP=
;                     ff = BIGPICKFILE(file='ch_ss.txt', tit='Type ascii file name ') ;GROUP=

                       IF ff NE  '' THEN BEGIN 

;                     IF units(1) EQ 'photons' THEN photons = 1 ELSE IF $
;                       units(1) EQ 'ergs' THEN photons =0 ELSE photons =0
;                     IF units(0) EQ 'Angstroms' THEN kev = 0 ELSE $
;                         IF units(0) EQ 'keV' THEN kev = 1 ELSE kev = 0

                          ch_line_list, spectrum, ff ,/spectrum,  /ascii, $
                                        minI=minI, wmin=xrange[0],wmax=xrange[1],all=theor_lines


                          widget_control, state.show_lines , /append, $
                                          set_val='Line intensities  saved in the ascii file '+ff

                          str=['Line details saved in the ascii file '+ff]
                          result = DIALOG_MESSAGE(str,/info)

                       ENDIF 
                    ENDIF  ELSE  BEGIN 
;bell
                       widget_control, state.show_lines , /append, $
                                       set_val='No lines  to save !!'
                    ENDELSE  
                 ENDIF 
              END  
           END 
        END 
     END 

     (event.id EQ  state.asave_sp) OR (event.id EQ  state.gsave_sp) OR $
        (event.id EQ  state.fsave_sp): BEGIN

;first re-check we have the right structure:

        IF  n_elements(spectrum) GT 0 THEN $ 
           result=ch_check_str(spectrum, /sp) ELSE result = 0


        IF result THEN BEGIN 

           IF  event.id EQ state.asave_sp THEN BEGIN 

;create an ascii file with two columns, lambdas and intensities.


              ff = dialog_pickfile(file='ch_ss_sp.ascii', tit='Type ascii file name ') ;GROUP=
;            ff = BIGPICKFILE(file='ch_ss_sp.ascii', tit='Type ascii file name ') ;GROUP=
              IF ff NE  '' THEN BEGIN 
                 openw, luo, ff, /get_lun

                 printf,luo,'% CH_SS: predicted spectrum'
                 printf,luo,'% CHIANTI database - Version '+spectrum.version

                 CASE tran.model_name OF
                    'Constant pressure': BEGIN 
                       dummy = 'Calculated with '+spectrum.model_name+ '= '+$
                               string(spectrum.model_pe,'(e9.2)')+ ' (cm-3 K) '
                    END 
                    'Constant density': BEGIN 
                       dummy = 'Calculated with '+spectrum.model_name+ '= '+$
                               string(spectrum.model_ne,'(e9.2)')+ ' (cm-3) '
                    END
                    'Function':BEGIN 
;         break_file, spectrum.model_file, disk, dir, file,ext
                       dummy = 'Calculated with the (Te,Ne) file: '+spectrum.model_file
                    END 
                 ENDCASE

                 printf,luo, '% '+dummy

                 w_min = min(spectrum.lambda)
                 w_max = max(spectrum.lambda)

                 IF spectrum.units(0) EQ 'Angstroms'  THEN  $
                    printf,luo,'% '+strtrim(string(w_min,'(f6.1)'),2)+' to '+strtrim(string(w_max,'(f6.1)'),2)+' Angstroms' ELSE $
                       printf,luo,'% '+strtrim(string(w_min,'(f6.3)'),2)+' to '+strtrim(string(w_max,'(f6.3)'),2)+' keV' 
                 

                 printf,luo,'% FWHM='+trim(spectrum.INSTR_FWHM)

                 IF tag_exist(spectrum, 'CONTINUUM') THEN printf,luo,'% CONTINUUM=YES' ELSE $
                    printf,luo,'% CONTINUUM=NO'

; 'Angstroms' 

                 printf,luo,'% Units are: '+ spectrum.units(0)+ ' vs. '+spectrum.units(1)


                 IF tag_exist(spectrum, 'FILE_EFFAREA' ) THEN $
                    printf,luo,'% Effective area file: '+spectrum.FILE_EFFAREA

                 printf,luo,'% Calculated: '+spectrum.DATE

                 printf,luo,'% Ionization Fractions file: '+ spectrum.ioneq_name
                 printf,luo,'% Elemental Abundance file: '+ spectrum.ABUND_NAME

                 printf,luo,'% Minimum abundance = '+$ 
                        strtrim(string(spectrum.MIN_ABUND, format='(e9.2)'),2)

                 IF spectrum.ADD_PROTONS EQ 1 THEN printf,luo,'% Proton excitation: YES' ELSE $
                    IF spectrum.ADD_PROTONS EQ 0 THEN printf,luo,'% Proton excitation: NO'

                 IF spectrum.PHOTOEXCITATION EQ 1 THEN printf,luo,'% Photo-excitation: YES' ELSE $
                    IF spectrum.PHOTOEXCITATION EQ 0 THEN printf,luo,'% Photo-excitation: NO' 


                 IF tag_exist(spectrum, 'DEM') THEN $
                    printf,luo,'% Differential Emission Measure file: '+spectrum.dem_name ELSE BEGIN 

                    printf,luo,'% Calculated with the isothermal approximation:'

                    dummy1 = spectrum.logt_isothermal
                    dummy2 = spectrum.logem_isothermal

                    dummy1 = trim(float(dummy1))
                    dummy2 = trim(float(dummy2))

                    printf,luo,'% log T='+arr2str(dummy1, ',',/trim)
                    printf,luo,'% log EM='+arr2str(dummy2, ',',/trim)

                 END 


                 FOR i=0L, n_elements(spectrum.lambda)-1 DO $
                    printf, luo, spectrum.lambda[i], spectrum.spectrum[i]

                 free_lun, luo
                 widget_control, state.show_lines , /append, $
                                 set_val='Spectrum saved in the ascii file '+ff

              END  

           ENDIF    ELSE IF (event.id EQ  state.gsave_sp) OR $
              (event.id EQ  state.fsave_sp)  THEN BEGIN 

;create a save files

              IF (event.id EQ  state.gsave_sp) THEN BEGIN 

;               ff = BIGPICKFILE(file='ch_ss_sp.genx', tit='IDL genx save file name ') ;GROUP=
                 ff = dialog_pickfile(file='ch_ss_sp.genx', tit='IDL genx save file name ') ;GROUP=

                 IF ff NE  '' THEN BEGIN 
                    savegen, struct=spectrum  ,file=ff 
;bell
                    widget_control, state.show_lines , /append, $
                                    set_val='Spectrum and line details saved in the IDL genx  file '+ff
                 END 

              ENDIF  ELSE IF (event.id EQ  state.fsave_sp) THEN BEGIN 

                 ff = dialog_pickfile(file='ch_ss_sp.fits', tit='FITS file name ') ;GROUP=
;               ff = BIGPICKFILE(file='ch_ss_sp.fits', tit='FITS file name ') ;GROUP=
                 IF ff NE  '' THEN BEGIN 

                    head1 = ["COMMENT  ", $
                             "COMMENT  PEAK  The peak intensity of the line in the spectrum (approx. value)", $
                             "COMMENT  "]


                    head2 = ["COMMENT  ", $
                             "COMMENT  This is additional information on the SPECTRUM: ", $
                             "COMMENT  ", $

                             "COMMENT  LAMBDA      The array of X-values", $
                             "COMMENT  SPECTRUM    The array of Y-values", $
                             "COMMENT  UNITS       The units of LAMBDA, SPECTRUM ", $
                             "COMMENT  INSTR_FWHM  The Instrumental FWHM", $
                             "COMMENT  BIN_SIZE    Width of the Bins  (fixed) in angstroms", $
                             "COMMENT  ABUND_NAME  The CHIANTI abundance file name", $
                             "COMMENT  ABUND       The abundance values", $
                             "COMMENT  MIN_ABUND   The minimum abundance value used", $
                             "COMMENT  ABUND_REF   The references", $
                             "COMMENT  CONTINUUM   The values of the continuum (if calculated)", $
                             "COMMENT  EFFAREA       The array of effective area values (optional)", $
                             "COMMENT  FILE_EFFAREA  The name of the effective area file used.", $
                             "COMMENT  "]

                    ch_write_fits, spectrum, ff,  head1=head1, head2=head2

;bell
                    widget_control, state.show_lines , /append, $
                                    set_val='Spectrum and line details saved in the FITS file '+ff

                 END   
              END 
           END 

        ENDIF     ELSE BEGIN 
;bell
           widget_control, state.show_lines , /append, $
                           set_val='Error, no spectrum to save !!!'
        END  
     END    

     event.id EQ state.quit:  WIDGET_CONTROL, event.top, /DESTROY ; quit

  ENDCASE     
END 


;-----------------------------------------------------------------------------
PRO syn_wid, GROUP=Group
;-----------------------------------------------------------------------------


COMMON base_com, syn_main_base, isothermal_base1
COMMON line_data, tran, state
COMMON plot_spec,kev_yn, ang, inst, o_lines, xrange, yrange, sty, theor_lines, $
  cont_ind, lambda, effarea, file_effarea, spectrum,fold_yn, log
COMMON wind_data, wxsiz, wysiz, o_strength
COMMON abundance, abfile, abdir, abund, abund_ref, min_abund
COMMON calc_int, int_xrange, const_names, const_nt, const_value, $
  ioneqfile, ioneqdir, demfile, demdir, iso_logt, iso_logem, isothermal_flag, $
  all_ions_yn, list_ions, units, noprot, photoexcitation, rphot, radtemp

defsysv,'!xuvtop', EXISTS = EXISTS 
IF NOT EXISTS THEN $
  message, 'system variable !xuvtop must be set '

xuvtop = !xuvtop


device, get_screen_size = sz
sz(0) = 1150 < 0.9*sz(0)
sz(1) = 950 < 0.9*sz(1)

wxsiz = 1100 < round(sz(0))     ;1100
wysiz = 400 < round(sz(1)/3.5)  ; 400

;print, wxsiz, wysiz

; Main base for everything
;
SYN_MAIN_BASE = WIDGET_BASE(col=1,/frame, MAP=1,x_scroll=sz(0),y_scroll=sz(1), $
                            UVALUE='SYN_MAIN',title='CHIANTI Spectral Synthesis Package', space=0)



;m1 = widget_base(syn_main,row=1)
;quit = widget_button(m1,value='Quit',uvalue='EXIT',/no_rel,font=font)

lin_base0 = WIDGET_BASE(syn_main_base, col=1, space=0) ;,/frame)

calc_int_butt = widget_BUTTON(lin_base0, value='Line intensities calculation - click here for a short HELP - Send comments to Giulio Del Zanna ', font=font)


lin_base = WIDGET_BASE(lin_base0,row=1, /frame, space=0)

angstrom = string(197B)
units = strarr(2)

ss = widget_base(lin_base, /frame)
ss5 =  WIDGET_BASE(ss, /column, /frame) 

;unit_info0 = widget_button(ss5,value='Units: '+angstrom)
unit_info0 = widget_label(ss5,value='Wavelength ['+angstrom+']')
;default:
units(0) = 'Angstroms'


ss1    = widget_base(ss5, /row)
wlss   = widget_label(ss1,value='Min.',font=font)
sss1   = widget_base(ss1,/row)
wminw1 = widget_text(sss1,/edit,xsize=8,ysize=1,uvalue='MINW1',$
                     value=strtrim(string(format='(f11.1)',xrange(0)),2), font=font)

ss2    = widget_base(ss5,/row)
wlss  = widget_label(ss2,value='Max.',font=font)
sss2   = widget_base(ss2,/row)
wmaxw1  = widget_text(sss2,/edit,xsize=8,ysize=1,uvalue='MAXW1',$
                      value=strtrim(string(format='(f11.1)',xrange(1)),2), font=font)


cc1 = WIDGET_BASE(ss5, /column, /frame)

c1 = WIDGET_BASE(cc1, /col)     ;, /frame)
pp=widget_label(c1,value='Model options')

const_names = ['Density (cm-3)', 'Pressure (cm-3 K)', 'Function (Ne,Te)']

desc = [{PSELECT_S, btext:'Const. Density', mtext:'Constant Density (cm-3)', uvalue:'0', flags:0}, $
        {PSELECT_S, btext:'Const. Pressure', mtext:'Constant Pressure (cm-3 K)', uvalue:'1', flags:0}, $
        {PSELECT_S, btext:'Function', mtext:'Function (Ne,Te)', uvalue:'2', flags:0} ]

const_widg = cw_pselect(c1, '', desc)


dummy = widget_base(c1, ysize=30, xsize=80, space=0)
const_read = WIDGET_TEXT(dummy,  value=strtrim(string(format='(e9.2)', const_value),2),$
                         /editable, /all, sensitive=1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;  selection of ionization fraction data
;


dummy=WIDGET_BASE(lin_base, /col, space=0,/frame)

ioneq_help_button=WIDGET_BUTTON(dummy,  value='Ion abund. HELP')

ibase2= WIDGET_BASE(dummy, /col, space=0)
ioneq_button=WIDGET_LABEL(ibase2,val='Calculate ion abund.?')

dummy_sel=WIDGET_BASE(ibase2, ysize=30)
ibase_sel=CW_BGROUP(dummy_sel,['Yes','No'],/row,$
                    set_value=0,/exclusive,/NO_RELEASE,$
                    uvalue=[1,0])

calc_ioneq_flag=1

output_ioneq_base=WIDGET_BASE(dummy, /col, space=0, map=1)
output_ioneq_label=widget_label(output_ioneq_base, value='Output file:')

if file_exist('adv.ioneq') then begin
              output_ioneq_file='adv_'+$
                  trim(string(SYSTIME(/JULIAN, /UTC), format='(f14.4)'))+'.ioneq'

endif else output_ioneq_file='adv.ioneq'

output_ioneq_show=WIDGET_TEXT(output_ioneq_base, value=output_ioneq_file, xsiz=14,/edit)


;----------------------------

ioneq_base =  WIDGET_BASE(dummy, /col, space=0, map=0) 

desc = [{PSELECT_S, btext:'Ioniz. Fraction', mtext:'Ioniz. Fraction selection:', uvalue:'0', flags:0}, $
        {PSELECT_S, btext:'Ioniz. Fraction', $
         mtext:'---CHIANTI DIRECTORY: *.ioneq---', uvalue:'0', flags:0}]

files = findfile(concat_dir(concat_dir(!xuvtop, 'ioneq'), '*.ioneq'))

FOR i=0,  N_ELEMENTS(files)-1 DO BEGIN 
   break_file,files[i], disk, dir, f, ext 
   desc = [desc, {PSELECT_S, btext:'Ioniz. Fraction', mtext:f[0]+ext, $
                  uvalue:files[i], flags:0}]
ENDFOR 

;get files in the working directory:

files = findfile('*.ioneq')

IF N_ELEMENTS(files) GT 0 AND files[0] NE '' THEN BEGIN 

   desc = [desc, {PSELECT_S, btext:'Ioniz. Fraction', $
                  mtext:'---WORKING DIRECTORY: *.ioneq---', uvalue:'0', flags:0}]

   FOR i=0,N_ELEMENTS(files)-1 DO BEGIN 
      break_file,files[i], disk, dir, f, ext 
      desc = [desc, {PSELECT_S, btext:'Ioniz. Fraction', mtext:f[0]+ext, $
                     uvalue:files[i], flags:0}]
   ENDFOR
ENDIF 

desc = [desc, {PSELECT_S, btext:'Ioniz. Fraction', $
               mtext:'SELECT FILES WITH WIDGET', uvalue:'1', flags:0}]

ioneq_pdmenu = cw_pselect(ioneq_base,'', desc)

dummy = widget_base(ioneq_base, ysize=30, xsize=120,space=0)

ioneq_show=WIDGET_TEXT(dummy, value=ioneqfile)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ibase1= WIDGET_BASE(lin_base, /col, space=0,/frame)
;-----------------------------------------------------------------
input_ibase1=WIDGET_BASE(ibase1, row=3, space=0, map=1)

input_ibase2 = widget_label(input_ibase1, value='min log T [K]: ')
input_ibase3 = WIDGET_BASE(input_ibase1, /row, ysize=30, xsize=40)

min_logt_ev = WIDGET_TEXT(input_ibase3,/edit, $
                          value=strtrim(string(format='(f3.1)', 4.0 )))
min_logt=4.0

input_ibase4 = widget_label(input_ibase1, value='max log T [K]: ')
input_ibase4 = WIDGET_BASE(input_ibase1, /row, ysize=30, xsize=35)

max_logt_ev = WIDGET_TEXT(input_ibase4,/edit, $
                          value=strtrim(string(format='(f3.1)', 8.0 )))
max_logt=8.0

input_e2 = widget_label(input_ibase1, value='d log T:')
input_e3 = WIDGET_BASE(input_ibase1, /row, ysize=30, xsize=50)
dlogt_ev = WIDGET_TEXT(input_e3,/edit,$
                           value=strtrim(string(format='(f4.1)', 0.1 )))
dlogt=0.1

input_ct_base=WIDGET_BASE(ibase1, ysize=30)
ct_base_sel=CW_BGROUP(input_ct_base,['CT on','CT off'],/row,$
                    set_value=1,/exclusive,/NO_RELEASE,$
                    uvalue=[1,0])
ct_flag=0
atmosphere_file=''

ct_base =widget_base(ibase1, column=1, /frame, map=0)
ct_show=WIDGET_TEXT(ct_base, value='', xsiz=14)


;-----------------------------------------------------------------


basic1_base = widget_base(lin_base, /col, /frame)

all_ions_base = WIDGET_BASE(basic1_base ,column=1, /frame)
all_ions_butt = WIDGET_BUTTON(all_ions_base,  value='All ions? - HELP')

dummy =  WIDGET_BASE(all_ions_base, ysize=30)
all_ions_ev = CW_BGROUP(dummy,['no', 'yes'],/row,$
                        set_value=all_ions_yn,/exclusive,/NO_RELEASE, $
                        uvalue=[0, 1])

rem_theor1_base = WIDGET_BASE(basic1_base ,column=1, /frame)
rem_theor1_butt = WIDGET_BUTTON(rem_theor1_base,  value='All lines? HELP')

dummy =  WIDGET_BASE(rem_theor1_base, ysize=30)
rem_theor1 = CW_BGROUP(dummy,['no','yes'],/row,$
                       set_value= 1 ,/exclusive,/NO_RELEASE,$
                       uvalue=[0,1]) ;, label_top='Add "unobserved" lines ?', /frame)

chck=getenv('CHIANTI_LOOKUP')
IF chck NE '' THEN BEGIN 
  lookup_base = WIDGET_BASE(basic1_base ,column=1, /frame)
  lookup_butt = WIDGET_BUTTON(lookup_base,  value='Lookup tables? HELP')

  dummy =  WIDGET_BASE(lookup_base, ysize=30)
  lookup_ev = CW_BGROUP(dummy,['no', 'yes'],/row,$
                        set_value=0,/exclusive,/NO_RELEASE, $
                        uvalue=[0, 1])
  lookup_available=1b
ENDIF ELSE BEGIN
  lookup_butt=0
  lookup_ev=0
  lookup_available=0b
ENDELSE 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


temp_base1 = widget_base(lin_base, column=1, /frame)

temp_base2 = WIDGET_BASE(temp_base1 ,column=1, space=0)
temp_butt = WIDGET_BUTTON(temp_base2,val=' HELP')


dummy_label=WIDGET_BUTTON(temp_base2,val='ISOTHERMAL?')
dummy = WIDGET_BASE(temp_base2, ysize=30)
temp_base=CW_BGROUP(dummy,['Yes','No (DEM)'],/row,$
                    set_value=1,/exclusive,/NO_RELEASE,$
                    uvalue=[1,0])
isothermal_flag = 0


isothermal_base1 = WIDGET_BASE(temp_base1,map=0 , row=2, space=0)

isothermal_base2 = widget_label(isothermal_base1, value='Log T [K]: ')
isothermal_base3 = WIDGET_BASE(isothermal_base1, /row, ysize=35, xsize=80)

iso_logt_ev = WIDGET_TEXT(isothermal_base3,/edit,/ALL_EVENTS, $
                          value=strtrim(string(format='(f3.1)', 6.0 )))
;String(13B) +; SCR_XSIZE=20,
isothermal_e2 = widget_label(isothermal_base1, value='Log EM' +  '[cm-5]:', YSIZE=35)
;isothermal_e2b = widget_label(isothermal_base1, value= '[cm-5]:', YSIZE=35)

isothermal_e3 = WIDGET_BASE(isothermal_base1, /row, ysize=35, xsize=60)
iso_logem_ev = WIDGET_TEXT(isothermal_e3,/edit,/ALL_EVENTS,$
                           value=strtrim(string(format='(f4.1)', 27.0 )))


basic2_base = widget_base(lin_base, /col, /frame)

;
;  selection of DEM data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

dem_base =widget_base(basic2_base, column=1, /frame, map=1)


desc = [{PSELECT_S, btext:'DEM', mtext:'DEM selection:', uvalue:'0', flags:0}, $
        {PSELECT_S, btext:'DEM', $
         mtext:'---CHIANTI DIRECTORY: *.dem---', uvalue:'0', flags:0}]

files = findfile(concat_dir(concat_dir(!xuvtop, 'dem'), '*.dem'))

FOR i=0,  N_ELEMENTS(files)-1 DO BEGIN 
   break_file,files[i], disk, dir, f, ext 
   desc = [desc, {PSELECT_S, btext:'DEM', mtext:f[0]+ext, $
                  uvalue:files[i], flags:0}]
ENDFOR 

;get files in the working directory:

files = findfile('*.dem')

IF N_ELEMENTS(files) GT 0 AND files[0] NE '' THEN BEGIN 

   desc = [desc, {PSELECT_S, btext:'DEM', $
                  mtext:'---WORKING DIRECTORY: *.dem---', uvalue:'0', flags:0}]

   FOR i=0,N_ELEMENTS(files)-1 DO BEGIN 
      break_file,files[i], disk, dir, f, ext
      desc = [desc, {PSELECT_S, btext:'DEM', mtext:f[0]+ext, $
                     uvalue:files[i], flags:0}]
   ENDFOR

ENDIF 

desc = [desc, {PSELECT_S, btext:'DEM', $
               mtext:'SELECT FILES WITH WIDGET', uvalue:'1', flags:0}]


dem_pdmenu = cw_pselect(dem_base, '', desc)

dem_show=WIDGET_TEXT(dem_base, value='', xsiz=8) 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pexc_base = widget_base(lin_base, /column, /frame, map=1)

;photoexcitation
pe_help_button = WIDGET_BUTTON(pexc_base,val='PE? - HELP')


photoexcitation_info = widget_button(pexc_base,value='PE: NO',$
                                     frame=2)
;define the default:

photoexcitation = 0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

photoexcitation_base=widget_base(pexc_base, row=2 , map=0)

dummy_base =  widget_label(photoexcitation_base, value='R/Ro:')
photoexcitation_rphot_ev = WIDGET_TEXT(photoexcitation_base,/edit, $
                                       xsize=6,value=strtrim(string(format='(f6.3)', 1.0 ), 2))


dummy_base =  widget_label(photoexcitation_base, value='T [K]:')
photoexcitation_radtemp_ev = WIDGET_TEXT(photoexcitation_base,/edit, $
                                         xsize=6,value=strtrim(string(format='(f6.1)', 6000.0 ), 2))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

dummy_base = widget_base(lin_base, /column)

unit_info = widget_button(dummy_base,value='Units: ERGS',uvalue='UNIT',$
                          frame=2)

;define the default:
units(1) = 'ergs'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

proton_info = widget_button(dummy_base,value='Protons: YES',uvalue='PROTON_RATES',$
                            frame=2)
;define the default:
noprot = 0


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

calc_lines = widget_button(dummy_base, value='Calculate intensities',uvalue='CALC_LINES') 
;,/no_rel,font=font)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

save_base = widget_button(dummy_base, value='Save/Restore ', menu=2)

gsave_str =  widget_button(save_base, value='SAVE line intensities (IDL)')
save_str =  widget_button(save_base, value='SAVE line intensities (FITS)')

grestore_str =  widget_button(save_base, value='RESTORE line intensities (IDL)')
restore_str =  widget_button(save_base, value='RESTORE line intensities (FITS)')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;dummy_base = widget_base(lin_base, /column)

quit = widget_button(dummy_base, $
                     value='Quit',uvalue='EXIT',/no_rel,font=font)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

angstrom = string(197B)
wavelength = '!7k!X'



;ang_base = WIDGET_BASE(syn_main, $
;      row=1, MAP=1,       UVALUE='ang_base')

sp_base0 = WIDGET_BASE(syn_main_base, col=1, space=0)

calc_sp_butt = widget_BUTTON(sp_base0, value='Calculate and plot a spectrum - click here for a short HELP - ', font=font)

sp_base =WIDGET_BASE(sp_base0,row=1,/frame)

ss5 =  WIDGET_BASE(sp_base, /col, /frame) 

kev_info= widget_button(ss5,value=angstrom+' [ / keV ]')
;default in Anstroms:
kev_yn=0

;unit_info3 = widget_button(ss5,value='Units: '+angstrom)
;unit_info3 = widget_label(ss5,value='Units: '+angstrom)
units(0) = 'Angstroms'

ss1    = widget_base(ss5, /row, ysize=30, space=0)
wlss   = widget_label(ss1,value='Min.',font=font)
sss1   = widget_base(ss1,/col)
wminw2 = widget_text(sss1,/edit,xsize=9,ysize=1,uvalue='MINW1',$
                     value=strtrim(string(format='(f13.6)',xrange(0)),2), font=font)

ss2    = widget_base(ss5,/row, ysize=30, space=0 )
wlss  = widget_label(ss2,value='Max.',font=font)
sss2   = widget_base(ss2,/col)
wmaxw2  = widget_text(sss2,/edit,xsize=9,ysize=1,uvalue='MAXW1',$
                      value=strtrim(string(format='(f13.6)',xrange(1)),2), font=font)


ang_base = WIDGET_BASE(sp_base ,column=1, /frame, UVALUE='ang_base')
dummy = WIDGET_LABEL(ang_base,value='Bin')
ang_butt=WIDGET_BUTTON(ang_base,value='HELP', ysize=20)
;
ang_READ=WIDGET_TEXT(ang_base,$
                     value=strtrim(string(format='(f6.3)',ang),2),$
                     xsiz=7,/editable)

;inst_base = WIDGET_BASE(syn_main, $
;      row=1, MAP=1, $
;      UVALUE='inst_base',/frame)

inst_base = WIDGET_BASE(sp_base ,column=1, /frame)
dummy = WIDGET_LABEL(inst_base,value='FWHM')
inst_butt=WIDGET_BUTTON(inst_base,value='HELP', ysize=20)
;
inst_READ=WIDGET_TEXT(inst_base,$
                      value=strtrim(string(format='(f6.3)',inst),2),$
                      xsiz=7,/editable)


;continua = CW_BGROUP(syn_main,/frame,/row,/nonexclusive,  $
;                     ['Free-bound','Free-free','None'], $
;                     uval=[0,1,2], $
;                     set_value=cont_ind, $
;                     label_left='Continua')


continua_base = WIDGET_BASE(sp_base ,column=1, /frame)
;continua = WIDGET_BUTTON(continua_base,value='Continuum: NO')

dummy = WIDGET_LABEL(continua_base,value='Continuum?')
continua_butt = WIDGET_BUTTON(continua_base, value='HELP', ysize=20)
continua = CW_BGROUP(continua_base,/row,/exclusive,/NO_RELEASE,  $
                     ['No', 'Yes'], $
;                     uval=[1,0], $
set_value=0) 

rem_theor2_base = WIDGET_BASE(sp_base ,column=1, /frame)
dummy = WIDGET_LABEL(rem_theor2_base,  value='All lines?')
rem_theor2_butt = WIDGET_BUTTON(rem_theor2_base,  value='HELP', ysize=20)

rem_theor2 = CW_BGROUP(rem_theor2_base,['no','yes'],/row,$
                       set_value= 0 ,/exclusive,/NO_RELEASE, $
                       uvalue=[0,1]) ;, $
;                        label_top='Add "unobserved" lines ?', /frame)


;
; Abundance base
;

ab_base = WIDGET_BASE(sp_base, column=1, /frame, space=0) ;/ALIGN_BOTTOM,


desc = [{PSELECT_S, btext:'Abundances', mtext:'Elemental abundance selection:', uvalue:'0', flags:0}, $
        {PSELECT_S, btext:'Abundances', $
         mtext:'---CHIANTI DIRECTORY: *.abund---', uvalue:'0', flags:0}]

files = FINDFILE(concat_dir(concat_dir(!xuvtop, 'abundance'), '*.abund'))

FOR i=0,  N_ELEMENTS(files)-1 DO BEGIN
   break_file,files[i], disk, dir, f, ext 
   desc = [desc, {PSELECT_S, btext:'Abundances', mtext:f[0]+ext, $
                  uvalue:files[i], flags:0}]
ENDFOR  


;get files in the working directory:

files = findfile('*.abund')

IF N_ELEMENTS(files) GT 0 AND files[0] NE '' THEN BEGIN 

   desc = [desc, {PSELECT_S, btext:'Abundances', $
                  mtext:'---WORKING DIRECTORY: *.abund---', uvalue:'0', flags:0}]

   FOR i=0,N_ELEMENTS(files)-1 DO BEGIN
      break_file,files[i], disk, dir, f, ext
      desc = [desc, {PSELECT_S, btext:'Abundances', mtext:f[0]+ext, $
                     uvalue:files[i], flags:0}]
   ENDFOR 
ENDIF 

desc = [desc, {PSELECT_S, btext:'Abundances', $
               mtext:'SELECT FILES WITH WIDGET', uvalue:'1', flags:0}]

dummy = widget_base(ab_base, xsize=80, ysize=30, space=0)

ab_pdmenu = CW_Pselect(dummy,'', desc) 

ab_button = WIDGET_BUTTON(ab_base,value='HELP', ysize=20)

dummy = widget_base(ab_base, xsize=100, ysize=30)
ab_show=WIDGET_TEXT(dummy, value=abfile)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

min_abund_base1 = WIDGET_BASE(sp_base, column=1, /frame) ;/ALIGN_BOTTOM,
dummy = WIDGET_LABEL(min_abund_base1,value='Min. Abund.')
min_abund_button = WIDGET_BUTTON(min_abund_base1,value='HELP', ysize=20)

min_abund_ev = widget_text(min_abund_base1,/edit,$
                           xsize=4,uvalue='MIN_ABUND',$
                           value=strtrim(string(MIN_ABUND, format='(e9.2)'),2), font=font)

;WIDGET_CONTROL, min_abund_ev,set_value=strtrim(string(MIN_ABUND, format='(e9.2)'),2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

fold_base = WIDGET_BASE(sp_base, column=1, /frame)
fold_info = widget_button(fold_base,value='Eff. Area: NO')
fold_yn = 0
fold_button = WIDGET_BUTTON(fold_base,value='HELP', ysize=20)

dummy = WIDGET_BASE(fold_base, xsize=100, ysize=30, space=0)
fold_show = widget_text(dummy, value='-',  sensitive=1) 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;calculation button:
;
;mm = WIDGET_BASE(sp_base , /ALIGN_BOTTOM, /frame)

dummy = WIDGET_BASE(sp_base ,column=1)

dummy2 = widget_button(dummy, value='RESTORE spectrum', menu=2)

grestore_sp =  widget_button(dummy2, value='RESTORE spectrum (IDL)')
restore_sp =  widget_button(dummy2, value='RESTORE spectrum (FITS)')


unit_info2 = widget_button(dummy,value='Units: ERGS',uvalue='UNIT',$
                           frame=2)

;define the default:
units(1) = 'ergs'

calc_sp = widget_button(dummy, value='Calculate and plot ',uvalue='CALC_SP') 
;,/no_rel,font=font)



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; Base for graphics window
;
WIND_BASE=WIDGET_BASE(SYN_MAIN_base, frame=2,$
                      col=1, $
                      MAP=1, $
                      xsiz=wxsiz,$
                      UVALUE='WIND_BASE', space=0)

;
; Plotting window
;
PLOT_RAT = WIDGET_DRAW( WIND_BASE, $
                        RETAIN=2, $
                        UVALUE='PLOT_RAT', $
                        XSIZE=wxsiz, $
                        YSIZE=wysiz, $
                        /motion_events, /button_events,  sensitive=0)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

plot_base = WIDGET_BASE(syn_main_base,row=1,/frame, space=0)

lop_base =  WIDGET_BASE(plot_base ,/row, space=0)

oplot_base = WIDGET_BASE(lop_base ,column=1)
oplot_butt = WIDGET_BUTTON(oplot_base,val='Labels? - HELP')
oplot_lines=CW_BGROUP(oplot_base,['yes','no'],/row,$
                      set_value=1,/exclusive,/NO_RELEASE,$
                      uvalue=[1,0])

;WIDGET_CONTROL,oplot_lines , SET_VALUE=1


strength_base = WIDGET_BASE(lop_base ,column=1, ysize=30, space=0)
strength_butt = WIDGET_BUTTON(strength_base,val='Min. - HELP')
strength_yn = WIDGET_TEXT(strength_base, value=strtrim(string(format='(e9.2)',o_strength)),$
                          xsiz=6,/editable)

xy_ranges_base =  WIDGET_BASE(plot_base, /col, space=0)

x1 = widget_base(xy_ranges_base,/row, ysize=30, space=0)
x1l = widget_label(x1, value='X:', font=font)
xx1 = widget_base(x1,/row)
xmin_base = widget_text(xx1,/edit,xsize=9,ysize=1,uvalue='MINX',$
                        value=strtrim(string(format='(f13.6)',xrange(0)), 2), font=font)
xmax_base = widget_text(xx1,/edit,xsize=9,ysize=1,uvalue='MAXX',$
                        value=strtrim(string(format='(f13.6)',xrange(1)), 2), font=font)

x2 = widget_base(xy_ranges_base,/row, ysize=30, space=0)
x2l = widget_label(x2, value='Y:', font=font)
xx2 = widget_base(x2,/row)
ymin_base = widget_text(xx2,/edit,xsize=9,ysize=1,uvalue='MINY',$
                        value=strtrim(string(format='(e9.2)',yrange(0)), 2), font=font)
ymax_base = widget_text(xx2,/edit,xsize=9,ysize=1,uvalue='MAXY',$
                        value=strtrim(string(format='(e9.2)',yrange(1)), 2), font=font)



;---------------------------------------0
;  buttons
;
;EXTRAS = CW_BGROUP( syn_main, $

dummy_base =  WIDGET_BASE(plot_base , /row)

EXTRAS = CW_BGROUP(dummy_base, $
                   ['Zoom', 'Unzoom'],$
                   /col)


EXTRAS1 = CW_BGROUP( dummy_base, $
                     ['Create PS file',  'Hardcopy'],$
                     /col)


EXTRAS2 = CW_BGROUP( dummy_base, $
                     [ 'Save line details (latex)','Save line details (ascii)'], $
                     /col)


dummy = WIDGET_BASE(plot_base , col=1)

log_lin_info =  WIDGET_button(dummy, value='Lin [/Log]')
log = 0

dummy2 = widget_button(dummy, value='SAVE spectrum', menu=2)
asave_sp = widget_button(dummy2, value='Save spectrum (ascii)')
gsave_sp = widget_button(dummy2, value='Save spectrum (IDL)' )
fsave_sp = widget_button(dummy2, value='Save spectrum (FITS)')


;EXTRAS3 = CW_BGROUP( dummy_base, $
;                     ['Save spectrum (ascii)', 'Save spectrum (IDL)', 'Save spectrum (FITS)' ],$
;                     /col ,/NO_RELEASE)


;---------------------------------------0


;
; Window for displaying information about the lines selected with mouse
;
show_lines = WIDGET_TEXT(syn_main_base,ysiz=8, $
                         value='Watch this space for information',/scroll,/wrap)

; if n_elements(atmosphere_file) eq 0 then atmosphere_file=''


state={calc_int_butt:calc_int_butt,unit_info:unit_info,unit_info2:unit_info2,$
       unit_info0:unit_info0,$  ;unit_info3:unit_info3,$
       proton_info:proton_info, $
       photoexcitation_info:photoexcitation_info, $
       calc_lines:calc_lines, calc_sp: calc_sp, quit:quit, $
       wminw1:wminw1, wmaxw1:wmaxw1,  wminw2:wminw2, wmaxw2:wmaxw2,$
       const_widg:const_widg, const_read:const_read, $
       ibase_sel:ibase_sel,input_ibase1:input_ibase1,ioneq_base:ioneq_base, ioneq_pdmenu:ioneq_pdmenu, ioneq_show:ioneq_show, $
       ioneq_help_button:ioneq_help_button,$
       output_ioneq_show:output_ioneq_show,output_ioneq_file:output_ioneq_file,$
       calc_ioneq_flag:calc_ioneq_flag,min_logt_ev:min_logt_ev, max_logt_ev:max_logt_ev,dlogt_ev:dlogt_ev ,$
       min_logt:min_logt, max_logt:max_logt,dlogt:dlogt,$
       input_ct_base:input_ct_base,ct_flag:ct_flag, ct_base_sel:ct_base_sel, ct_base:ct_base, ct_show:ct_show,atmosphere_file:atmosphere_file,$
;-----------       
       dem_pdmenu:dem_pdmenu, dem_show:dem_show, dem_base:dem_base, $
       rem_theor1:rem_theor1,rem_theor1_butt:rem_theor1_butt,$
       all_ions_ev:all_ions_ev, all_ions_yn:all_ions_yn, all_ions_butt:all_ions_butt, $
       xmin_base:xmin_base, xmax_base:xmax_base, ymin_base:ymin_base, ymax_base:ymax_base, $
       save_str:save_str, restore_str:restore_str, $ 
       gsave_str:gsave_str, grestore_str:grestore_str, $ 
; save_ascii:save_ascii,
temp_base:temp_base, temp_butt:temp_butt, $
  isothermal_base1:isothermal_base1, $
       iso_logt_ev:iso_logt_ev, iso_logem_ev:iso_logem_ev, $
       pe_help_button:pe_help_button,$
       pexc_base: pexc_base, $
  photoexcitation_base:photoexcitation_base, $
  photoexcitation_rphot_ev:photoexcitation_rphot_ev, $
  photoexcitation_radtemp_ev:photoexcitation_radtemp_ev, $
  log_lin_info:log_lin_info, $
  extras:extras,extras1:extras1, extras2:extras2,$
  asave_sp:asave_sp, gsave_sp:gsave_sp, fsave_sp:fsave_sp, $
; extras3:extras3, $
calc_sp_butt:calc_sp_butt, $
kev_info:kev_info,kev_yn:kev_yn,  ang_read:ang_read, ang_butt:ang_butt, $
  inst_read:inst_read, inst_butt:inst_butt, $
  show_lines:show_lines, plot_rat:plot_rat, $
  oplot_lines:oplot_lines, oplot_butt:oplot_butt, $
  rem_theor2:rem_theor2, rem_theor2_butt:rem_theor2_butt, $
  continua:continua,continua_butt:continua_butt, $
  ab_pdmenu:ab_pdmenu, ab_show:ab_show,ab_button:ab_button, $
  min_abund_ev:min_abund_ev,min_abund_button:min_abund_button, $
  strength_butt:strength_butt, strength_yn:strength_yn,$
       fold_info:fold_info, fold_button:fold_button, fold_show:fold_show, $
       lookup_butt: lookup_butt, lookup_ev: lookup_ev, lookup_available: lookup_available, $
  grestore_sp:grestore_sp, restore_sp:restore_sp}


WIDGET_CONTROL, SYN_MAIN_base, /REALIZE, set_uvalue=state


; Get drawable window index

WIDGET_CONTROL, plot_rat, GET_VALUE=plot_rat_id 
wset,plot_rat_id

;      WIDGET_CONTROL,state.const_widg,set_uvalue='0'

XMANAGER, 'SYN_MAIN', SYN_MAIN_base,  modal=modal, group_leader=group_leader

END


;-----------------------------------------------------------------------------
PRO CH_SS,  font=font
;-----------------------------------------------------------------------------

COMMON line_data, tran, state
COMMON plot_spec,kev_yn, ang, inst, o_lines, xrange, yrange, sty,theor_lines, $
  cont_ind, lambda, effarea, file_effarea, spectrum, fold_yn, log
COMMON wind_data, wxsiz, wysiz, o_strength
COMMON abundance, abfile, abdir, abund, abund_ref,min_abund 
COMMON calc_int, int_xrange, const_names, const_nt, const_value,$
  ioneqfile, ioneqdir, demfile, demdir, isothermal,iso_logem, isothermal_flag, $
  all_ions_yn, list_ions, units, noprot,photoexcitation, rphot, radtemp


defsysv,'!xuvtop', EXISTS = EXISTS 
IF NOT EXISTS THEN $
  message, 'system variable !xuvtop must be set  '
xuvtop = !xuvtop

;
;  initial checks
;

device,decomposed=0

;
;  check not already in use
;
if xregistered('CH_SS') then begin
   print,' CH_SS program is already registered.'
   return
endif

;
;  set only possible graphics device and load LUT
;set_plot,'X'

;save for later all the !p variable settings:

original_background = !p.background
original_color = !p.color
original_multi = !p.multi

tvlct,r,g,b,/get

loadct,0

!p.background = 255
!p.color = 0
!p.multi = 0


;widget dimensions:
;font_s='-adobe-helvetica-medium-r-normal--8-80-75-75-p-46-iso8859-1'
;font_n='-adobe-helvetica-bold-r-normal--12-120-75-75-p-70-iso8859-1'


IF n_elements(font) EQ  0 THEN BEGIN 
   font = get_dfont('-adobe-helvetica-bold-r-normal*120*75*')
;font = get_dfont('-adobe-helvetica-bold-r-normal*75*75*')
   IF font(0) NE '' THEN font = font(0) ELSE font = 'fixed'
END 

widget_control,default_font=font

;clear from memory:

delvarx, tran
delvarx, ang, inst, o_lines, xrange, yrange, sty,theor_lines
delvarx, cont_ind, lambda, effarea, file_effarea, spectrum, fold_yn, log
delvarx, wxsiz, wysiz, o_strength 
delvarx,  abfile, abdir, abund, abund_ref,min_abund 
delvarx, int_xrange, const_names, const_nt, const_value
delvarx,  ioneqfile, ioneqdir, demfile, demdir, isothermal,iso_logem, isothermal_flag
delvarx,   all_ions_yn, list_ions, units, noprot,photoexcitation, rphot, radtemp


;DEFAULT values of 

o_lines = 0                     ; don't display line labels by default
o_strength = 0.   
sty = 0                         ;do not zoom

;define as default not to show lines with only theoretical wavelengths:
theor_lines = 0


;define as default to calculate all ions

all_ions_yn = 1

;default: no continuum 
cont_ind = 0

;define a default abundance:

defsysv,'!abund_file', EXISTS = EXISTS 
IF  EXISTS THEN abund_name = !abund_file ELSE BEGIN 
   message, 'No default abundance file (!abund_file) defined ??'
END 

break_file, abund_name, disk, abdir, abfile, ext
abdir = disk+abdir
abfile = abfile+ext


ioneqfile=file_basename(!ioneq_file)
ioneqdir=file_dirname(!ioneq_file)


read_abund, abund_name ,abund,abund_ref

;define the minimum abundance 

index = where(abund GT 0.0)
min_abund = min(abund(index))

fold_yn = 0

;tran =''
ang = 0.1
const_nt = 0  
const_value=1e10
xrange =[150., 200.]
yrange = [0.,0.]

int_xrange = xrange 

const = 0.

;instrumental FWHM
inst = 0.5 

syn_wid


;RESTORE THE ORIGNAL SETTINGS 

!p.background = original_background 
!p.color = original_color
!p.multi = original_multi

;and  restore the color table

tvlct,r,g,b

END
