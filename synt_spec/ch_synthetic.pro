;+
; PROJECT:  CHIANTI
;
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the
;       University of Cambridge, Goddard Space Flight Center, and University of Michigan.
;
;        The recent update of this program was developed as part of CHIANTI-VIP. 
;        CHIANTI-VIP (Virtual IDL and Python) is  a member of the CHIANTI family
;        mantained by Giulio Del Zanna, to develop additional features and
;        provide them to the astrophysics community.
;
; NAME:
;	CH_SYNTHETIC
;
; PURPOSE:
;
;       to calculate CHIANTI line intensities or G(T) and output an IDL structure. 
;
; PROCEDURE:
;
;       This routine calculates as default line intensities for a user-specified 
;       differential emission measure and ionisation balance. The actual 
;       creation of a synthetic spectrum (i.e., wavelength vs. intensity) 
;       is performed by other routines - see CH_SS.PRO and 
;       MAKE_CHIANTI_SPEC.PRO.
;
;       Note that this routine does not include the element abundances 
;       in the line intensities, as this will be performed by 
;       make_chianti_spec. One of the reasons why  element abundances are not
;       included in the line intensities calculation is so that it is easier 
;       for the user to see how  modifying abundances affects their spectra in
;       e.g. CH_SS.PRO. 
;
;       The calculations are performed at constant pressure or 
;       at constant density, or with a user-defined grid of Te,Ne.   
;
;       The routine can also output line intensities calculated with an
;       isothermal approximation.
;
;       -----------------------------------------------------------------------------------       
;       ******* BEFORE CHIANTI VERSION 11, THE ROUTINE REQUIRED A
;               IONIZATION FRACTION FILE.
;       
;       ***** With the keyword /ADVANCED_MODEL,  WHICH IS THE DEFAULT,
;
;       extra effects are included. Currently, only
;       density-dependent  and charge transfer (CT) effects are included, see
;
;       Dufresne, R.P. and Del Zanna, G., 2019,  A&A, 626, A123
;       Modelling Ion Populations in Astrophysical Plasmas: Carbon in the Solar Transition Region
;
;       Dufresne, R.P. and Del Zanna, G., 2020, MNRAS, 497 (2), 1443
;       Effects of density on the oxygen ionisation equilibrium in collisional plasmas,
;       
;       Dufresne, R. P., Del Zanna, G., and Storey, P. J., 2021,MNRAS, 505, 3968
;       Modelling low charge ions in the solar atmosphere
;
;       Dufresne, R. P., Del Zanna, G., et al., 2024, ApJ (10.48550/arXiv.2403.16922)
;       -----------------------------------------------------------------------------------       
;
;       The routine then calculates the ion charge states on the fly.
;
;       To speed up the routine, an input array of log T can be input via the
;       keyword IONEQ_LOGT, otherwise the grid logT=4.0 -- 8.0 in 0.05 dex intervals is used.
;       Optionally, the atmospheric model and the helium abundance can be input via
;       ATMOSPHERE, HE_ABUND.
;       The routine calls CH_CALC_IONEQ to calculate the ion fractions, assuming several
;       approximations. If advanced models are not available for an element, the previous
;       zero density CHIANTI approximation is used.
;
;       Note: it is also possible to call CH_CALC_IONEQ first, write an ionization equilibrium
;             file, then use it within this routine by setting ADVANCED_MODEL=0 and
;             giving as input the ionization equilibrium file.
;
;       ****** If  ADVANCED_MODEL=0  then as before:
;       -----------------------------------------------------------------------------------       
;
;	If the isothermal approximation is not used, then the user will be asked
;	to select two  files, that can either be in the 
;	standard CHIANTI database or in the working directory. 
;
;       These files are: 
;       
;       - an ionization fraction file 
;       - a differential emission measure (DEM) file.
;
;       The routine can also output the contribution functions G(T) of the lines,
;       instead of the intensities, if the keyword GOFT is used. In this case,
;       only the ionization equilibrium file needs to be selected.
;       The G(T), or intensity per emission measure, is calculated as:
;
;        G=(hc/lambda_ij)*A_ji*(N_j(X^+m)/N(X^+m))*(N(X^+m)/N(X))/ N_e /(4.*!pi)
;
;       where A_ji is the A-value of the transition;
;             (N_j(X^+m)/N(X^+m)) is the population of the upper level,
;             calculated by solving the statistical equilibrium equations; 
;             (N(X^+m)/N(X)) is the ionization equilibrium
;             N_e is the electron density.
;
;       unless    /PHOTONS is set, in which case the  (hc/lambda_ij) factor
;       is not included.  
;
;       If not specified otherwise, with the use of the MASTERLIST or SNG_ION
;       keywords,  then the standard masterlist of the ions, which has 
;       all the ions in the current CHIANTI database, is used.
;
;       PROGRAMMING NOTES
;
;       The DEM is not assumed to be specified at 0.1 logT intervals (which 
;       is how the ion fraction are specified). Thus this routine reads 
;       in the DEM vs. logT information and then uses the IDL spline 
;       function to tabulate the DEM over 0.1 logT intervals. The minimum 
;       and maximum temperatures are those in the DEM file, rounded up to 
;       the nearest 0.1. The new DEM function tabulated over 0.1 logT 
;       intervals is contained in 'dem_int'.
;
;       For some of the dielectronic files, radiative decays that were in 
;       the standard .wgfa file will also be present in the dielectronic 
;       version of the .wgfa file. In these cases the line intensity 
;       produced from the latter file needs to be ignored and so we have a 
;       check in ch_synthetic to do this. An example is the 1-7 decay in 
;       the ca_19.wgfa and ca_19d.wgfa files. In the latter case, the 
;       model of the ion does not include electron excitation to level 7 
;       and so the model for the 1-7 decay is incorrect, hence we ignore 
;       it.
;
; CATEGORY:
;
;	spectral synthesis.
;
; CALLING SEQUENCE:
;
;       IDL> ch_synthetic,wmin,wmax, output=output, pressure=pressure,$
;            [MODEL_FILE=MODEL_FILE, err_msg=err_msg, msg=msg, $
;            density=density,all=all,sngl_ion=sngl_ion, $
;            photons=photons,  masterlist=masterlist, $
;            save_file=save_file , verbose=verbose, $
;            logt_isothermal=logt_isothermal,$
;            logem_isothermal=logem_isothermal,$
;            goft=goft, ioneq_name=ioneq_name, dem_name=dem_name,$
;            noprot=noprot, rphot=rphot, radtemp=radtemp, progress=progress ]
;
;
;
; INPUTS:
;
;	Wmin:  minimum of desired wavelength range in Angstroms
;	Wmax:  maximum of desired wavelength range in Angstroms
;
;
; OPTIONAL INPUTS :
;
;       PRESSURE:  pressure in emitting region (Pe,  cm^-3 K). 
;                  Only a single value is accepted, and the calculation is
;                  performed at constant pressure.
;
;       DENSITY:   density in emitting region (Ne, cm^-3). 
;                  Only a single value is accepted, and the calculation is
;                  performed at constant  density, unless LOGT_ISOTHERMAL is
;                  defined. In this case, DENSITY can be an array of values, but
;                  has to have the same number of elements as LOGT_ISOTHERMAL.
;
;
;       MODEL_FILE    Full path of the (Te,Ne) file if defined. 
;                     This file should have two columns, one with the Te (K)
;                     values, and one with the Ne (cm^-3) values. If these
;                     values are not sorted in ascending order of Te, the
;                     routine does sort them.
;              **** The  Ne values are interplotated over the grid of temperatures
;                   for the calculation, defined by either IONEQ_LOGT or
;                   LOGT_ISOTHERMAL or the temperatures in the IONEQ file.
;                     
;
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
;
;	IONEQ_NAME:  Name of the ionization equilization name to use.  If not
;		     passed, then the user is prompted for it, unless the
;                    ADVANCED_MODEL is used.
;
;                    **** Note: if (by default) the ADVANCED_MODEL is used,
;                       IONEQ_NAME will be the name of the file where the new
;                       ion charge states are written. 
;
;
;
;       IONEQ_LOGT: an array of log T [K] values, defining the grid for the
;                   calculation, unless the isothermal option is called, or
;                   an ion fraction file is used. 
;
;
;       ATMOSPHERE: A file with the H,He abundances as a function of temperature.
;                      By default, the file avrett_atmos.dat is read, with data from
;                      Avrett E.H., Loeser R., 2008, ApJ, 175, 229
;
;       ATMOS_PARAMS: This is a structure containing the atmospheric parameters, and is
;                     an alternative to giving ATMOSPHERE (see above). The tags are:
;
;                     .h_elec  ratio of hydrogen to electron number density (required)
;                     .h1_frac neutral hydrogen fraction (required)
;                     .he1_frac neutral helium fraction
;                     .he2_frac singly ionised helium fraction
;                     .temp    the temperatures (K) at which the above parameters are
;                              defined.
;
;                     The helium data are optional (helium ion fractions will be
;                     calculated from !ioneq_file if they are not specified). If the tags
;                     elec_dens, pressure and height are present, then they will be added
;                     to the output structure, but they are not essential for
;                     incorporating charge transfer.
;
;                     Special case: if h_elec and h1_frac are scalars, then they are
;                     applied to all of the input temperatures TEMP. This is specifically
;                     for creating lookup tables for a range of hydrogen (and helium)
;                     parameters. In this case the temp tag in atmos_params is ignored.
;
;       HE_ABUND:  The total helium abundance relative to hydrogen. 
;
;
;       DEM_NAME:  Name of the DEM file to used.  If not passed, then the user
;		   is prompted for it. Note that DEM is expected to be
;		   defined as N_e*N_H*dT/dh and *not* N_e^2*dT/dh.
;
;       LOGT_ISOTHERMAL
;                  Array of logarithmic temperatures.
;                  If defined, the emissivities are calculated with an
;                  isothermal approximation. The values are sorted in ascending
;                  order.
;
;       LOGEM_ISOTHERMAL
;                  Array of logarithmic emission measures.
;                  If defined, the emissivities are calculated with an
;                  isothermal approximation. The values are sorted in ascending
;                  order. If LOGT_ISOTHERMAL is specified without 
;                  LOGEM_ISOTHERMAL then the emission measures are set to 1 
;                  (logem_isothermal=0).
;
;       RADTEMP   The blackbody radiation field temperature (default 6000 K).
;
;       RPHOT    Distance from the centre of the star in stellar radius units.
;                I.e., RPHOT=1 corresponds to the star's surface. (Default is
;                infinity, i.e., no photoexcitation.)
;
;       RADFUNC         The name of a user-defined function that will generate
;                       a radiation spectrum as a function of temperature. 
;                       This radiation field will replace the black-body that
;                       is assumed when using the RADTEMP keyword in the call
;                       to pop_solver.
;
; OUTPUTS:
;
;       OUTPUT:    The name of the structure containing the line intensities and
;                  details.  
;                  The tags of the  structure are:
;
;       .lines     A structure containing information about the lines. 
;                  Its size is the number of lines in the spectrum. The 
;                  tags are:
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
;                  .tmax   The temperature of maximum emission of the line.
;
;                          If the G(T) are output, tmax is the maximum of G(T).
;
;                          If the isothermal approximation is used  tmax=0.
;
;                          If a DEM is used,  tmax is the maximum of the 
;                          emissivity that includes the product of the ion
;                          fraction and the DEM.
;                          Rounded to nearest 0.1
;
;                  .wvl    Wavelength of the transition, in Angstroms.
;
;                  .flag   A flag, =-1 if the line has only theoretical energy
;                          levels. Otherwise flag=0.
;
;                  .int    Intensity of line (erg/cm2/s/sr or phot/cm2/s/sr), 
;                          divided by the element abundance (exclusive with .goft). 
;
;                  .goft   The G(T) of the line (optional /exclusive with .int).
;
;
;       .ioneq_name     The ion balance file used (full path).
;       .ioneq_logt        The Log10 T values associated.
;       .ioneq_ref      The references.
;
;       .dem_name       The differential emission measure file eventually  used
;                       (full path).
;       .dem            The Log10 DEM values 
;       .dem_logt          The Log10 T values associated.
;       .dem_ref        The references.
;
;       .model_name    A string indicating the model used: 
;
;                    1- Constant density
;                    2- Constant pressure
;                    3- Function (Te,Ne)
;
;       .model_file    Full path of the (Te,Ne) file if defined. Null string otherwise.
;
;       .model_ne    the Ne value(s).
;
;                     - a scalar if 'Constant density' is selected.
;                     - an array if 'Function' is selected.
;                     - 0. if constant pressure is selected.
;
;       .model_pe    the Pe value.
;
;                     - a scalar if constant pressure is selected.
;                     - 0. if 'Constant density' is selected.
;                     - an array=density*temperature if 'Function' is selected.
;                          
;       .model_te    the Te values if 'Function' is selected. Otherwise 0.
;
;       .wvl_units  The wavelength units.
;
;       .wvl_limits    The wavelength limits specified by the user.
;
;       .int_units  The intensity units.
;
;                   1) If LOGT_ISOTHERMAL is defined, we have two cases:
;                      a) LOGEM_ISOTHERMAL is not defined, and is therefore
;                         assumed to be 0 (EM=1). In this case, units are
;                         'photons cm+3 sr-1 s-1' or 'erg cm+3 sr-1 s-1'.
;                      b)  LOGEM_ISOTHERMAL is defined. In this case, units are
;                         'photons cm-2 sr-1 s-1' or 'erg cm-2 sr-1 s-1'.
;
;                   2) If LOGT_ISOTHERMAL is not defined, we have two cases:
;                      a) intensities are calculated. In this case, units are
;                         'photons cm-2 sr-1 s-1' or 'erg cm-2 sr-1 s-1'.
;                      b) Contribution functions G(T) are calculated. In this
;                         case, units are 
;                         'photons cm+3 sr-1 s-1' or 'erg cm+3 sr-1 s-1'.
;
;       .logt_isothermal
;                  The Log10(T) values used. 
;
;       .logem_isothermal
;                  The Log10(EM) values used. 
;
;       .date      The date and time when the structure was created.
;
;       .version   The version number of the CHIANTI database used.
;
;       .lookup    Takes the value of keyword_set(lookup)
;
;       .add_protons 
;                  A flag (0/1) to indicate whether proton data were used (1)
;                  or not (0) to calculate the level population.
;
;       .photoexcitation
;                  A flag (0/1) to indicate if photoexcitation was included (1)
;                  or not (0).
;
;       .radtemp 
;                 The blackbody radiation field temperature used (if
;                 photoexcitation was included).
;
;       .rphot
;              Distance from the centre of the star in stellar radius units  
;              (if photoexcitation was included).
;
;
; OPTIONAL OUTPUTS:
;
;
;       SAVE_FILE: If defined, then an IDL save file is created, with the output
;                  structure. 
;
;
;       GOFT:      If set,  the G(T) of the lines are calculated, and put in
;                  the output structure, instead  of the line intensities.
;                  Units are 'photons cm+3 sr-1 s-1' or 'erg cm+3 sr-1 s-1'
;
;
;
; KEYWORDS:
;
;       PHOTONS:   The output intensities will be in photons instead of 
;                  ergs.
;
;       VERBOSE:   If set, the routine will list each ion it is looking at, 
;                  and how many lines from each ion it is including in the 
;                  spectrum.
;
;       GOFT:      If set,  the G(T) of the lines are calculated, and put in
;                  the output structure, instead  of the line intensities.
;                  Units are 'photons cm+3 sr-1 s-1' or 'erg cm+3 sr-1 s-1'
;
;       NOPROT     Switch off the inclusion of proton rates in the level 
;                  balance (default).
;
;      PROGRESS    If set, a widget appears, showing the progress of the
;                  calculation and allowing the user to halt the calculation.
;
;      NO_SUM_INT  Prevents the summing of intensities over temperature. 
;                  Only works in conjunction with the LOGT_ISOTHERMAL 
;                  option, and is implemented in order to work the 
;                  ISOTHERMAL routine. The .INT tag in OUT.LINES becomes 
;                  an array with the same number of elements as 
;                  LOGT_ISOTHERMAL, corresponding to the intensities at 
;                  each temperature.
;
;       FRAC_CUTOFF     The fraction of non-zero elements in the C matrix below
;                       which the sparse matrix solver is used. See the routine
;                       matrix_solver for more details.
;       NOIONREC: If set, then level-resolved ionization/recombination
;                rates are not read for the ion, even if they exist.
;
;       NO_RREC: If set, do not read the level-resolved radiative
;                recombination (RR) files.
;
;
;       LOOKUP:  If set, then routine will attempt to use the CHIANTI
;                lookup tables to obtain level populations, rather
;                than call pop_solver. This leads to a much quicker
;                calculation. 
;
;       REGULAR:  If set, then the regular inversion method for
;                 pop_solver will be used (i.e., the IDL invert
;                 routine). [This is the default.] 
;
;       SPARSE:   If set, then the sparse matrix inversion routine
;                 (linbcg) for pop_solver will be used.
;
;       LAPACK:   If set, then the LAPACK matrix inversion routine
;                 (la_invert) for pop_solver will be used. 
;
;       ADVANCED_MODEL: include density-dependent (and CT) effects.
;
;       CT: include charge transfer in advanced models
;
;       DR_SUPPRESSION: Switch on DR suppression from Nikolic et al (2018) for all ions 
;              not included in the advanced models. The comparison with Summers (1974) suppression
;              has not been checked for other elements when preparing the models.
;
;       NO_AUTO: If set, then the autoionization rates (contained in
;                the .auto file) are not read. The autoionization states are not
;           included in the calculations, i.e. a single ion rather than the
;           two-ion model  introduced in version 9 is calculated. This speeds
;           up the calculations without affecting the lines from the bound states.
;
; CALLS:  CH_GET_FILE
;          many CHIANTI standard routines,  including:
;          READ_IONEQ, READ_DEM, READ_MASTERLIST, ION2SPECTROSCOPIC,
;          ZION2FILENAME, READ_WGFA,READ_ELVLC,READ_SPLUPS,POP_SOLVER,
;          DESCALE_UPS, CONVERT_TERMS. 
;          CONVERT_TERMS uses some standard SolarSoft routines: 
;          REPSTR, STR_INDEX, DATATYPE, 
;          VALID_NUM, DELVARX, INFO_PROGRESS, SAVEGEN
;
; COMMON BLOCKS:
;          wgfa, wvl,gf,a_value
;          upsilon,splstr
;          elvlc,l1a,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref 
;          elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref
;          radiative, radt, dilute
;          proton, pstr
;          ionrec,rec_rate,ci_rate,temp_ionrec,luprec,lupci,status
;
; RESTRICTIONS:
;
; SIDE EFFECTS:
;
; CATEGORY:
;	spectral synthesis.
;	
; EXAMPLE:
;
;       This routine can be called in this way:
;
;       IDL> ch_synthetic,5.,10., output=structure, pressure=1.e+15
;
;       To make use of the output structure, use MAKE_CHIANTI_SPEC or CH_SS
;       
;
; PREV. HIST. :
;       Based on synthetic.pro, written by Ken Dere
;
;
; WRITTEN     : 
;       Ver.1, 22-Jun-00, Peter Young (PRY) and Giulio Del Zanna (GDZ)
;
; MODIFICATION HISTORY:
;
;       Ver.1, 22-Jun-00, Peter Young and Giulio Del Zanna
;
;       Ver.2, 25-Jul-00, PRY
;               Removed /all keyword; make_chianti_spec can be 
;                used to filter out negative wavelengths.
;               Added flabel tag to output in order to pick out 
;                dielectronic recombination lines.
;
;       Ver.3, 4-Oct-00, PRY
;               Replaced /all keyword.
;               Corrected bug when .wgfa files contain two A-values 
;                for the same transition.
;
;       Ver.4, 5-Oct-00, PRY
;               Corrected bug that gave rise to lines from the same 
;                transition when the dielectronic file existed.
;
;       V.5, 11-Oct-2000, GDZ
;            eliminate the abundance call; reinstate the /masterlist keyword;
;            added the tag  ident_latex to have the identification in 
;            late-style format; added a tag flag=-1 for the unobserved lines,
;            and =0 otherwise; reinstated all wavelengths > 0. ;
;            added the calculation of the G(T);
;            added  a few other tags in the output, and various checks and
;            comments. 
;       V.6 15-Oct-2000 ,GDZ
;             Replaced calls to solarsoft routines to  standard IDL ones. 
;             Corrected an error in the output creation, in relation to the
;             isothermal case. Added isothermal in the output. added checks to
;             the wavelengths. Default output name is TRANSITIONS. changed
;             const_net and added const_net_value + a few other things.
;
;       v.7, 27-Nov-2000, GDZ. Corrected an error in the calculation of the
;       G(T). 
;
;       Version 8, 5-Dec-2000, GDZ, DAMTP. Fixed a bug when checking the 
;       values in the .splups files.
;
;       V. 9, GDZ, 10-Apr-2001, corrected another error in the G(T) calc.
;
;       V. 10, GDZ, 30-Oct-2001 added CHIANTI Version number, changed isothermal
;            to logt_isothermal and added logem_isothermal to the output.
;            Removed the use of log T values, and the calculation. 
;            Added err_msg, a text string with an error message.
; 
;       Version 11, 8-Nov-01, GDZ
;
;            Changed the MASTERLIST keyword. Allowed double use, as a keyword 
;            and as a string.
;
;       Version 12, 18-Nov-01, Peter Young
;
;            Added /NOPROT, RPHOT and RADTEMP keywords; changed upsilon 
;            common block.
;
;        Version 13, 29-Apr-02, GDZ
;
;            Added no_protons, photoexcitation, rphot, radtemp 
;            tags into the output  structure. 
;            Revised Header. Added the PROGRESS widget.
;            Added a check if the ion is present in the Ion. Frac. file.
;            Added informative MSG keyword.
;            Now uses  savegen.pro to save the structure.
;
;        V. 14, 28-May-2002, GDZ: 
;                  generalize directory concatenation to work for Unix, Windows
;                  and VMS. 
;
;           modified tags: 
;                          limits -> wvl_limits
;                          ioneq_t -> ioneq_logt
;                          wvlunits -> wvl_units
;                          intunits -> int_units
;                          time --> date 
;                          no_protons -> add_protons
;                         dem_t -> dem_logt 
;                const_nte -> model_name
;                const_nte_value -> model_ne, model_pe, model_te
;                   removed from the main STR:   .ioneq  ctemp 
;                   removed from the LINES STR:  fwhm flabel
;
;            Added model_file  input for model Ne(T). Had to considerably
;            modify the routine.
;
;         V. 15, 16-Jul-2002, Peter Young
;                  Added keyword /NO_SUM_INT.
;
;         V. 16, 22-Jul-2002, Peter Young
;                  Corrected a bug related to /NO_SUM_INT; logt_isothermal 
;                  can now be specified without logem_isothermal.
;
;         V. 17, 23-July-2002, GDZ
;                  Modified a few checks on the input. Also, now it prints the
;                  error message whenever the program aborts
;
;         V.18, 2-Aug-02, GDZ
;                  Replaced all DBLARR and DOUBLE calls with floats.
;                  Added a comment at the end of the routine when it finishes.
;
;         V.19, 8-Aug-02, GDZ 
;                  Added more error info. Changed the use of the DENSITY
;                  keyword. It is possible to input an array of values if
;                  LOGT_ISOTHERMAL is defined.
;
;         V. 20, 17-Sep-02, GDZ
;                  Corrected a bug: the functional (T,N) form
;                  is now only accepted if DENSITY is an array with at least two
;                  values.
;         V. 21, 19-Sep-02, GDZ
;                  Corrected the definition of the UNITS in case LOGT_ISOTHERMAL
;                  is defined.
;
;         V. 22, 19-Aug-03, Peter Young
;                  when logem_isothermal is input, the derived EM is now a
;                  DOUBLE array rather than FLOAT, preventing infinities when
;                  logem_isothermal values are large.
;
;         V. 23,   4-Oct-2003, GDZ
;                  modified the input to POP_SOLVER, now it is dealt with an
;                  input structure.
;
;         V.24,    10-Oct-2003, K.Dere
;                  added modifications from K.Dere, regarding the satellite
;                  lines. 
;
;         V 25,    3-Nov-2003, GDZ
;                 Added GROUP keyword, and modified so the progress widget can
;                 be stopped within IDL Windows.
;
;         V 26,    17-Apr-2004, Enrico Landi (EL)
;                  Added the recombination/ionization population processes.
;
;         V.27,    13-Apr-2005, EL
;                  Replaced the main loop to calculate individual line intensities
;                  with operations among arrays, to speed the whole program in case
;                  of large numbers of lines.
;
;         v.28,  31-Aug-2005, GDZ
;                 Fixed bug concerning the case when multiple temperatures 
;                 (i.e. logt_isothermal) were defined as input. The program
;                 was, in some cases, returning null values.
;                 The problem was the use of nt, the number of good temperatures
;                 for each ion, for the definition of the arrays, instead of
;                 using the number of logt_isothermal values (and the t_index).
;
;         v.29,  22-Aug-2008, Peter Young
;                 Changed list_goft_new and this_goft to be double
;                 precision, and changed str.goft to be double precision.
;
;         v.30,  12-Jun-2009, Enrico Landi
;                 Changed the definition of the temperature array for ion fractions
;                 in the IONREC variable, now taken directly from the output of
;                 READ_IONEQ.PRO
;
;         v.31,  23-Sep-2010, Peter Young
;                 Changed name of LIST array to MLIST for
;                 compatibility with IDL 8.
;
;         v.32,  20-Oct-2011, Peter Young
;                 Corrected bug with intensity units (int_units tag)
;                 when /GOFT was set.
;
;         v.33,  20-Apr-2012, Peter Young
;                 Now changed units for logem_isothermal=0 case to
;                 cm^-2 s^-1 sr^-1 (were cm^3 before).
;
;         v.34,  10-May-2013, Peter Young
;                  Modified call to pop_solver so that it uses the new
;                  /pressure keyword.
;
;         v.35  1-May-2014  GDZ  
;               for v.8:  replaced read_splups with read_scups. added 
;                keyword frac_cutoff
;
;         v.36 4-July-2014  GDZ  
;               Modified the way splstr is used (make use of data
;               tag).
;
;         v.37 2-Sep-2014, Peter Young
;               Now checks to make sure all the temperatures in
;               logt_isothermal are unique.
;
;         v.38 17-Sep-2015, Peter Young
;               Modified the line "junk=temporary(np)" for the case np
;               is undefined, as this causes a crash in earlier
;               versions of IDL.
;
;         v.39 21-Mar-2016, Peter Young
;               Changed "pickfile" call to "dialog_pickfile".
;
;         v.40 5-Feb-2018, Peter Young
;               Major update! Removed common blocks and added call to
;               ch_setup_ion to load atomic data. I checked the output
;               against the previous version by searching for lines
;               with intensities > 0.01% different over the wavelength
;               range 1-2000 angstroms for AR. A small anomaly was found due
;               to double precision temperatures in the previous
;               version (e.g., 6.800002 is larger than 6.800), which
;               meant some temperatures were excluded from calculation.
;
;          v.41 GDZ, 4 Oct 2018 
;               Major  modifications for v.9 
;          v.42 GDZ, 18  Nov  2018 
;               Reinstated the 'dielectronic' files, as not all the
;               ions in v.9 have the new models.
;          v.44, 5 Dec 2018 GDZ, revised the handling of pop_solver.
;          v.45, 14-Dec-2018 GDZ, added NOIONREC, NO_RREC keywords.
;          v.46, 1-Jul-2020, Peter Young, added /LOOKUP keyword and
;                 added lookup tag to output structure.
;          v.47, 24-Aug-2020, Peter Young, fixed problem if the
;                 temperature-density model file does not have a
;                 unique set of temperature points.
;
;          v.48, 28-Sep-2020, GDZ, fixed a major bug, all the
;          satellite lines without an excitation rate to their upper
;          levels were removed from the calculation. This was a left
;          over from the earlier approximated 'dielectronic'
;          files.
;          v.49, 30-Sep-2020, Peter Young
;            Added some extra comments related to above bug. No change
;            to code.
;          v.50, 8-Oct-2020, GDZ, added radfunc which is passed to
;          pop_solver.
;          v.51, 11-Dec-2020, Peter Young
;            Changed intensity and DEM to double-precision; introduced
;            /regular, /sparse and /lapack keywords that get passed to
;            pop_solver.
;          v.52, 12-Apr-2021, Peter Young
;            Modified the check on the dielectronic ions to ensure
;            that it doesn't affect other ions.
;          v.53, 01-Nov-2021, Peter Young
;            Removed references to /all keyword as it doesn't
;            do anything (all lines are included by default). The
;            keyword is retained, however, for backwards
;            compatibility.
;          v.54, 12-May-2023, Peter Young
;            Added /lookup option to ch_setup_ion call in order to
;            speed up routine; added /no_auto option; introduced levmax
;            for lookup tables; print statements only in /verbose set;
;            passed /verbose to ch_setup_ion; removed /verbose from
;            pop_solver as information not useful.
;          v.55, 24-May-2023, Peter Young
;            Fixed bug when passing /lookup to ch_setup_ion.
;          v.56, 13-Jun-2023, Peter Young
;            Lines coming from autoionizing levels are now flagged with
;            "s" (for satellite) in the lines.snote tag. The dielectronic
;            ("d") ions are no longer flagged with d in snote.
;          v.57, 22-Jun-2023, Peter Young
;            Modified definition of t_index for the /goft case.
;
;          v.58, 17 May 2024, Giulio Del Zanna
;
;           Major rewrite, adding the advanced model option (the default).
;
;           Also added  NO_AUTO: If set, then the autoionization states are not
;           included in the calculations, i.e. a single ion rather than the
;           two-ion model I introduced in version 9 is calculated. This speeds
;           up the calculations without affecting the lines from the bound states.
;           Also fixed the VERBOSE option which did not work properly
;
;           Ask to choose ionization equilibrium file for transparency, in
;           case the NO ADVANCED MODEL option is chosen.
;
;          v.59, 14 June 2024, Giulio Del Zanna
;           clarified warning/error messages.
;
;          v.60, 3-Jul-2024, GDZ, fixed the logt_isothermal case for the
;                advanced model, which was previously commented.
;                Also added the dr_suppression keyword.
;
;          v.61, 05 Sep 2024, Roger Dufresne
;                Minor alteration to reading system time when creating ioneq name.
;
;          v.62, 29-Apr-2025, Graham Kerr & Peter Young
;                Added atmos_params= optional input.
;
;   VERSION 62
;-
PRO info_progress, pct,lastpct,pctage, pct_slider_id,$
           interrupt_id,halt,quiet, snote,  group=group

  IF pct GE lastpct+pctage THEN BEGIN
     WHILE pct GT lastpct DO lastpct = lastpct+pctage
;     IF NOT quiet THEN print,lastpct,format='($,i4,a)','%'+string(13b)
     IF n_elements(pct_slider_id) EQ 1 THEN $
        widget_control,pct_slider_id,set_value=lastpct,/show
     IF n_elements(interrupt_id) EQ 1 THEN BEGIN
     widget_control,interrupt_id,set_value=$
          snote+' -- Click here to STOP calculation',/show
        ev = widget_event(interrupt_id,/nowait)
        IF ev.id NE 0L THEN BEGIN
           halt = 1
        END
     END 
  ENDIF

END 

PRO ch_synthetic, wmin, wmax, output=output, err_msg=err_msg, msg=msg, $
                  pressure=pressure, density=density, $
                  model_file=model_file, all=all,sngl_ion=sngl_ion, $
                  photons=photons,  masterlist=masterlist,  $
                  save_file=save_file , verbose=verbose,$
                  logt_isothermal=logt_isothermal,  logem_isothermal=logem_isothermal,$
                  goft=goft, ioneq_name=ioneq_name, dem_name=dem_name, $
                  noprot=noprot, rphot=rphot, radtemp=radtemp,RADFUNC=RADFUNC, progress=progress, $
                  no_sum_int=no_sum_int,  group=group,frac_cutoff=frac_cutoff, $
                  noionrec=noionrec, no_rrec=no_rrec, lookup=lookup, $
                  regular=regular, sparse=sparse, lapack=lapack, $
                  no_auto=no_auto,ioneq_logt=ioneq_logt, advanced_model=advanced_model,ct=ct,$
                  atmosphere=atmosphere,he_abund=he_abund,dr_suppression=dr_suppression,$
                  atmos_params=atmos_params

;
  if n_params() lt 2 then begin
     print,' IDL> ch_synthetic, wmin, wmax, output= , pressure= ,$' 
     print,'         [err_msg= , msg= ,$'
     print,'         [model_file=,  density= , all=all, sngl_ion= , $'
     print,'          /photons,  masterlist= ,save_file=  , /verbose, $   '
     print,'          logt_isothermal= , logem_isothermal= , $ '
     print,'          /goft, ioneq_name=, dem_name= , /noprot, $'
     print,'           rphot= , radtemp= , /progress'
     print, ''
     print, 'e.g.: IDL> ch_synthetic, 5.,10., output=structure, pressure=1.e+15'
     print, ''
     return
  endif
;
  t1 = systime(1)

  err_msg = '' &  msg=''

  defsysv,'!xuvtop', EXISTS = EXISTS 
  IF NOT EXISTS THEN $
     message, 'system variable !xuvtop must be set  '
  xuvtop = !xuvtop

;check the wavelengths
  wmin = float(wmin) & wmax=float(wmax)
  IF (wmin GE wmax) OR (wmin LE 0.) OR (wmax LE 0.) THEN BEGIN 
     err_msg ='%CH_SYNTHETIC: Error in the wavelengths - EXIT'
     print, err_msg
     return
  END 

; GDZ: by default VERBOSE=0
  if n_elements(verbose) eq 0 then verbose=0

  if keyword_set(verbose) then quiet=0 else quiet=1 ; for CH_SETUP_ION

; ------ NEW ADVANCED MODEL:  by default calculate the ion fractions on the fly
  if  n_elements(advanced_model) eq 0 then advanced_model=1 

;  ------ NEW ADVANCED MODEL: by default do not calculate CT
  if  n_elements(ct) eq 0 then ct=0 



  IF n_elements(rphot) EQ 0 THEN BEGIN 
     photoexcitation =0
     dilute=0d0

     IF n_elements(RADFUNC) eq 1 THEN BEGIN
        err_msg ='%CH_SYNTHETIC: Error ! if you define RADFUNC you need to define RPHOT   - EXIT'
        print, err_msg
        return
     END
     
  ENDIF ELSE BEGIN
     
     photoexcitation = 1 
     dilute=r2w(rphot)

;by default choose a T=1e6 if not defined.
     IF n_elements(radtemp) EQ 0 and n_elements(RADFUNC) eq 0 THEN BEGIN 
        radtemp=6d3
                                ;this is used in the common block:
        radt=radtemp
        IF  keyword_set(verbose) THEN print,'setting a default radiation temperature of 6000 K'
     ENDIF  ELSE if n_elements(radtemp) EQ 1 and n_elements(RADFUNC) eq 0 THEN begin
        radtemp=float(radtemp)
                                ;this is used in the common block:
        radt=radtemp
     endif  else if n_elements(radtemp) EQ 1 and n_elements(RADFUNC) eq 1 THEN begin
        err_msg ='%CH_SYNTHETIC: Error ! you cannot define both RADTEMP and RADFUNC - EXIT'
        print, err_msg
        return
     ENDIF  
  ENDELSE 


  IF keyword_set(no_sum_int) AND n_elements(logt_isothermal) EQ 0 THEN BEGIN
     print,'Error, the keyword /NO_SUM_INT only works in conjunction with'+ $
           ' LOGT_ISOTHERMAL - EXIT'
     return
  ENDIF

  IF  keyword_set(goft) AND n_elements(logt_isothermal) NE 0 THEN BEGIN 
     err_msg = '%CH_SYNTHETIC: Error, you cannot have the G(T) calculated at constant T !! - EXIT'
     print,err_msg
     return
  END 

  IF KEYWORD_SET(photons) THEN intunits = 'photons' ELSE intunits = 'erg'

;
; PRY, 20-Oct-2011
; I've put units specification in this separate piece of code.
; (I'm sceptical that units for the logem_isothermal=0 case
; should be cm^+3; I think they should be cm^-2.)
;
; PRY, 20-Apr-2012
; I've now changed the logem_isothermal=0 case to cm^-2  :-)
;
  IF keyword_set(goft) THEN BEGIN
     intunits=intunits+' cm+3 sr-1 s-1'
  ENDIF ELSE BEGIN
;  IF n_elements(logt_isothermal) NE 0 AND n_elements(logem_isothermal) EQ 0 THEN BEGIN
;    intunits = intunits+' cm+3 sr-1 s-1'
;  ENDIF ELSE BEGIN
     intunits=intunits+' cm-2 sr-1 s-1'
;  ENDELSE 
  ENDELSE 




;
; PRY, 20-Oct-2011
; I've removed the specification of units from the following
; section of code.
;
  IF n_elements(logt_isothermal)  NE 0 THEN BEGIN
     logt_isothermal = float(logt_isothermal)

     IF n_elements(logem_isothermal) EQ 0 THEN BEGIN 
        logem_isothermal=fltarr(n_elements(logt_isothermal))
     ENDIF


     IF n_elements(logem_isothermal) NE n_elements(logt_isothermal) THEN BEGIN 
        err_msg = '%CH_SYNTHETIC: Error, inconsistent number of Log T,EM values '
        print, err_msg
        return
     ENDIF ELSE logem_isothermal = float(logem_isothermal)

;
; PRY, 2-Sep-2014
; Check to make sure there aren't any duplicate elements in
; logt_isothermal (this fixes a problem identified by Nic Labrosse). 
;
     n1=n_elements(logt_isothermal) 
     chck=uniq(logt_isothermal)
     n2=n_elements(chck)
     IF n1 NE n2 THEN BEGIN
        err_msg = '%CH_SYNTHETIC: Error, there are duplicate temperatures in the LOGT_ISOTHERMAL array.'
        print, err_msg
        return
     ENDIF 

     i_sort = sort(logt_isothermal)
     logt_isothermal = logt_isothermal[i_sort]
     logem_isothermal = logem_isothermal[i_sort]
     tmax = 0.

  ENDIF


  IF n_elements(logt_isothermal) NE 0 AND n_elements(dem_name) NE 0 THEN BEGIN 
     err_msg ='%CH_SYNTHETIC: Error, you cannot have both the isothermal and DEM defined !! - EXIT'
     print,err_msg
     return
  END 


  CASE 1  OF 

     n_elements(density) GT  0: BEGIN 
        
        IF n_elements(density) EQ  1 THEN BEGIN 
           IF  keyword_set(verbose) THEN  $
              print,'%CH_SYNTHETIC: using constant density = ',density
           model_name = 'Constant density' ;+string(density)
           model_ne = density
           model_pe = 0.
           model_te = 0.
;         model_file = ' '
        ENDIF ELSE BEGIN 

           IF n_elements(density) EQ  n_elements(logt_isothermal) THEN BEGIN 
              IF  keyword_set(verbose) THEN  $
                 print,'%CH_SYNTHETIC: using a functional (T,N) form' 
              model_name = 'Function'
              model_ne = density
              model_pe = density*10.^logt_isothermal
              temperature = 10.^logt_isothermal
              model_te = 10.^logt_isothermal
;            model_file = ' '
           ENDIF  ELSE BEGIN  
              err_msg = '%CH_SYNTHETIC: Array of density values not allowed -- EXIT '
              print, err_msg 
              return
           END  
        END   
     END   
     
     n_elements(pressure) EQ  1: BEGIN 
        IF  keyword_set(verbose) THEN    print,'%CH_SYNTHETIC: using constant pressure = ', pressure
        model_name = 'Constant pressure' ;+string(pressure)
        model_pe =pressure
        model_ne = 0.
        model_te = 0.
;      model_file = ' '
     END

     n_elements(model_file) GT 0: BEGIN 

        IF file_exist(model_file) THEN BEGIN 
           data = read_ascii (model_file)
           temperature = reform( data.field1(0, *))
           density = reform( data.field1(1, *))

;sort the temperatures: 
           i_sort = sort(temperature)
           temperature = temperature[i_sort]
           density =density[i_sort]

           IF  keyword_set(verbose) THEN  $
              print,'%CH_SYNTHETIC: using a functional (T,N) form' 

           model_name = 'Function'
           model_ne = density
           model_pe = density*temperature
           model_te = temperature

        ENDIF ELSE BEGIN 
           err_msg = '%CH_SYNTHETIC: No  model file found -- EXIT '
           print, err_msg 
           return
        END 

     END  
     ELSE: BEGIN 
        err_msg ='%CH_SYNTHETIC: Error, you have to  either define: '+$
                 ' a) density; b) pressure; c) model (T,N) -- EXIT'
        print,err_msg
        return
     END 
  ENDCASE 

;IF n_elements(output) EQ 0 THEN output = output


  wvlunits = 'Angstroms'
;
; wvlunits = 'keV' 

  if not keyword_set(advanced_model)  then begin 
;
; PRY, 31-Jan-2018
;  No longer ask user to choose ioneq file, but simply set to our
;  default file.
; IF n_elements(ioneq_name) EQ 0 THEN ioneq_name=!ioneq_file
;
; GDZ - revert to original way, ask to choose ionization equilibrium file for transparency.
     IF n_elements(ioneq_name) EQ 0 then begin
        dir=concat_dir(!xuvtop, 'ioneq')
        ioneq_name=ch_get_file(path=dir,filter='*.ioneq',title='Select ionization equilibrium file')
     endif
     
     ff = findfile(ioneq_name)
     IF  ff(0)  NE ''  THEN $
        read_ioneq,ioneq_name,ioneq_logt,ioneq,ioneq_ref ELSE BEGIN 
        err_msg = '%CH_SYNTHETIC: Error,  no ioneq file found ! -- EXIT'
        print,err_msg
        return
     END 

     IF KEYWORD_SET(verbose) THEN BEGIN 
        print, ''
        FOR i=0, n_elements(ioneq_ref)-1 DO print, ioneq_ref(i)
        print, ''
     END

     dlnt=ALOG(10.^(ioneq_logt[1]-ioneq_logt[0]))      

  endif     else begin

; GDZ- ADVANCED model:
;---------------------
     
     IF n_elements(ioneq_name) EQ 0 then begin
        print, '% CH_SYNTHETIC: ADVANCED ionization model WARNING: *** YOU HAVE NOT DEFINED THE NAME OF THE IONEQ FILE ...'

        pp=strsplit(anytim(!stime,/vms),/extract)
        ioneq_name='ch_adv_'+trim(pp[0])+'-'+strmid(pp[1],0,8)+'.ioneq'
        
        print, '% CH_SYNTHETIC: IONEQ FILENAME THAT WILL BE WRITTEN in the working directory is: '+ioneq_name
        print, '% CH_SYNTHETIC:  This file should be kept if creating a synthetic spectrum that includes the continuum.'
        print, '% CH_SYNTHETIC:  Otherwise it may be safely moved/deleted once the routine has completed.'
        print, ' '
        
       
     endif else begin
        IF file_exist(ioneq_name) THEN BEGIN 
           err_msg = '% CH_SYNTHETIC ERROR, ADVANCED ionization model option requested but the output ioneq file '+ioneq_name+' already exists, please give it a different name!  -- EXIT'
           print,err_msg
           return
        END 
     endelse

     IF n_elements(logt_isothermal) GT 0  THEN BEGIN
        
      if  n_elements(ioneq_logt) eq 0 then begin
               ioneq_logt = logt_isothermal 
            if verbose then print,'% CH_SYNTHETIC: calculating G(T) with the input logt_isothermal'
         endif else begin 
             err_msg = '% CH_SYNTHETIC ERROR: input Temperature via ioneq_logt '+$
              ' cannot be given if logt_isothermal is defined - EXIT '
            print,err_msg
            return
         endelse 
      endif else begin   
     
; we have asked for the advanced model but we have not supplied the input temperatures:     
     if  n_elements(ioneq_logt) eq 0 then begin
        ioneq_logt =findgen(81)/20.0+4.0
        if verbose then print,'% CH_SYNTHETIC: calculating G(T) in the log T=4-8 range'
     end
     
  endelse      
   endelse ; advanced ioneq model 
  
  n_ioneq_logt=n_elements(ioneq_logt)
  
; PRY, 31-Jan-2018
;  No longer compute proton-to-electron ratio as this will be done
;  from within pop_solver.
;
  IF NOT keyword_set(noprot) THEN add_protons=1 ELSE add_protons=0


  IF n_elements(logt_isothermal) EQ 0  THEN BEGIN

     CASE keyword_set(goft) OF 

        1: BEGIN 

;this is the counter for the number of G(T):
           n_goft = 0

           if not keyword_set(advanced_model) then $         
              gdt=WHERE(ioneq_logt GE 0) else begin

; ; GDZ- ADVANCED model: if input temperature array not defined:
              
              if  n_elements(ioneq_logt) eq 0 then begin
                 ioneq_logt =findgen(81)/20.0+4.0
                 n_ioneq_logt=n_elements(ioneq_logt)
                 if verbose then print,'% CH_SYNTHETIC: calculating G(T) in the log T=4-8 range'
              endif else n_ioneq_logt=n_elements(ioneq_logt)
              
           endelse 

        END 

        0: BEGIN 

;
; Choose DEM file
;
           IF n_elements(dem_name) EQ 0 THEN BEGIN
              dir=concat_dir(!xuvtop, 'dem')
              dem_name=ch_get_file(path=dir,filter='*.dem',title='Select DEM File')
           END

           ff = findfile(dem_name)

           IF  ff(0)  NE '' THEN $
              read_dem,dem_name,dem_logt,dem,dem_ref ELSE BEGIN 

              err_msg ='%CH_SYNTHETIC: No DEM file found -- EXIT'
              print,err_msg
              return
           END 

           IF KEYWORD_SET(verbose) THEN BEGIN 
              print, ''
              FOR i=0, n_elements(dem_ref)-1 DO print, dem_ref(i)
              print, ''
           END

; GDZ- ADVANCED model: if input temperature array not defined:
           if  n_elements(ioneq_logt) eq 0 then begin
              ioneq_logt =findgen(81)/20.0+4.0
              n_ioneq_logt=n_elements(ioneq_logt)
              if verbose then print,'% CH_SYNTHETIC: calculating G(T) in the log T=4-8 range'
           endif else n_ioneq_logt=n_elements(ioneq_logt)
           
           
;
; convert the DEM distribution to match the temperatures in ioneq_logt
;

           IF model_name EQ  'Function' THEN $
              gdt=WHERE((ioneq_logt GE MIN(dem_logt)) AND $
                        (ioneq_logt LE MAX(dem_logt)) AND (ioneq_logt GE min(alog10(temperature))) AND $
                        (ioneq_logt LE  max(alog10(temperature))) , nn) ELSE $ 
                           gdt=WHERE((ioneq_logt GE MIN(dem_logt)) AND $
                                     (ioneq_logt LE MAX(dem_logt)) , nn) 

;ask for at least three points:
           IF nn LT 3 THEN BEGIN 
              err_msg = '%CH_SYNTHETIC: No sufficient overlap in Temperature '+$
                        ' - please change the Temperature array -- EXIT '
              print,err_msg
              return
           END 


           dem_int1=10.^(SPLINE(dem_logt,dem,ioneq_logt[gdt]))
           dem_int=dblARR(n_ioneq_logt)
           ngt=n_elements(gdt)
           FOR igt=0,ngt-1 DO  dem_int[gdt[igt]]=dem_int1[igt] >  0.


           IF KEYWORD_SET(verbose) THEN BEGIN 

;plot the DEM

;      circle_sym,/fill
              window,0,xs=600,ys=600
              wset,0
              nb = strpos(dem_ref(0),':')
;plot_oo,10.^dem_logt,10.^dem,psym=8,syms=1.5,$
              plot, ioneq_logt, alog10(dem_int),psym=-5,syms=1.2,$
                    tit=strmid(dem_ref(0),nb+1,100),xtit='Log T',$
                    ytit='Log DEM',chars=1.4, /yno
              wait, 2
              wshow,!d.window,0
           END  
        END   
     ENDCASE   
  ENDIF ;  ELSE BEGIN              ; log T isothermal 

; now we have ioneq_logt defined, one way or the other. 
;------------------------------------------------------

  
;check that we have some overlap in T:
;---------------------------------------

  IF model_name EQ  'Function' THEN  BEGIN 
     IF n_elements(logt_isothermal) EQ 0 THEN BEGIN
        IF keyword_set(goft) THEN BEGIN 
           t_index=WHERE( (ioneq_logt GE min(alog10(temperature))) AND $
                          (ioneq_logt LE  max(alog10(temperature))))
        ENDIF ELSE BEGIN
           t_index=WHERE((dem_int NE 0.)  AND $
                         (ioneq_logt GE min(alog10(temperature))) AND $
                         (ioneq_logt LE  max(alog10(temperature))) )
        ENDELSE 

     ENDIF   ELSE $
        t_index = where((logt_isothermal GE MIN(ioneq_logt)) AND $
                        (logt_isothermal LE MAX(ioneq_logt)) AND $
                        (logt_isothermal GE min(alog10(temperature))) AND $
                        (logt_isothermal  LE  max(alog10(temperature))) )

     IF t_index[0] EQ  -1 THEN BEGIN
        err_msg = '%CH_SYNTHETIC: No Temperature overlap!! - please change parameters - EXIT '
        print,err_msg
        return
     END 
  END                           ; function  

; GDZ 
  IF KEYWORD_SET(verbose) AND   n_elements(model_file) gt 0 THEN BEGIN 

     IF n_elements(t_index) GT 1 THEN BEGIN 
        IF n_elements(logt_isothermal) EQ 0 THEN TEMP=10.^ioneq_logt[t_index] ELSE $
           TEMP = 10.^logt_isothermal[t_index]

        DENS = 10.^SPLINE(alog10(temperature),alog10(density), alog10(TEMP) )

        break_file, model_file, disk,dir,f,ext
        f = f+ext

        window,1,xs=600,ys=600
        plot_oo, temp, DENS , psym=-2, syms=.8,$
                 tit= 'File (Te,Ne): '+f+' curve:interpolated values', $
                 xtit='T (K)', ytit='Ne (cm-3 K)',chars=1.4, /yno
        oplot, temperature, density,  psym=5,  syms=1.5
        wait, 2
        wshow,!d.window,0
     END 
  END  

  IF keyword_set(no_sum_int) AND n_elements(logt_isothermal) NE 0 THEN BEGIN
     list_int=dblarr(n_elements(logt_isothermal),1)
  ENDIF ELSE list_int=0d0
  list_wvl=0.d
  list_flag = 0
  list_ident = ''
  list_ident_latex = ''
  list_snote = ''
  list_l1 = 0
  list_l2 = 0
  list_tmax = 0.
  list_ion = 0
  list_iz = 0
  list_goft = dblarr(n_ioneq_logt,1)

; GDZ- ADVANCED model:
;---------------------
; If one or a list of ions is given as input:
  IF  n_elements(sngl_ion) GT 0 THEN  BEGIN 
     mlist=sngl_ion
  ENDIF  ELSE  BEGIN

;  open the file that has the names of the ions
;----------------------------------------------
     mname=''

     IF n_elements(MASTERLIST) NE 0 THEN BEGIN 

        IF valid_num(MASTERLIST) THEN BEGIN 
           IF MASTERLIST EQ 1 THEN BEGIN 
              dir=concat_dir(!xuvtop,'masterlist')
                                ;
                                ; PRY, 21-Mar-2016: pickfile is obsolete so I've changed
                                ; to using diaglog_pickfile
              mname=dialog_pickfile(path=dir,filter='*', $
                                    title='Select list of IONS File')
           ENDIF ELSE $
              mname=concat_dir(concat_dir(!xuvtop,'masterlist'), 'masterlist.ions')

        ENDIF ELSE IF datatype(MASTERLIST, 0) EQ 'STR' THEN BEGIN 
           IF file_exist(MASTERLIST) THEN mname=MASTERLIST ELSE $ 
              mname=concat_dir(concat_dir(!xuvtop,'masterlist'), 'masterlist.ions')
        ENDIF ELSE mname=concat_dir(concat_dir(!xuvtop,'masterlist'), 'masterlist.ions')

     ENDIF ELSE mname=concat_dir(concat_dir(!xuvtop,'masterlist'), 'masterlist.ions')

;make sure it exists


     IF NOT file_exist(mname)  THEN BEGIN 
        err_msg ='%CH_SYNTHETIC: No masterlist file found -- EXIT!'
        print,err_msg
        return
     END 

;IF  keyword_set(sngl_ion) THEN  BEGIN 
;   mlist=sngl_ion
;ENDIF  ELSE  BEGIN
     read_masterlist,mname,mlist

  ENDELSE


;------ NEW ADVANCED MODEL HERE --- GDZ --------

; if the option ADVANCED MODEL is used (default), we should input the
; array ioneq_logt, otherwise gets defined here. 

; by default calculate the ion fractions
  if  keyword_set(advanced_model) then  begin 

; We need  first to find which elements need to be calculated. 
     requested_iz_ions=intarr(n_elements(mlist))
     
     for iion=0,n_elements(mlist)-1 do begin 
        convertname,mlist[iion],iz,ion
        requested_iz_ions[iion]=iz
     end
     requested_iz=requested_iz_ions[rem_dup(requested_iz_ions)]

                                ;   get elements list
                                ;
     zlabl=['H','He','Li','Be','B','C','N','O','F','Ne','Na',$
            'Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr', $
            'Mn','Fe','Co','Ni','Cu','Zn']

     elements= zlabl[requested_iz-1] 
     
     err_calc=''
;     IF KEYWORD_SET(verbose) THEN

     print,'% CH_SYNTHETIC: calculating ion fractions ...'

; calculate on the fly but also WRITE the new ionization equilibrium, for later use.
     
     ioneq= ch_calc_ioneq(10.^ioneq_logt, outname=ioneq_name, $
                          density=density,pressure=pressure,model_file=model_file,$
                          /advanced_model, ct=ct,$
                          elements=elements,$
                          atmosphere=atmosphere,he_abund=he_abund,verbose=verbose,$
                          err_msg =err_calc, warning_msg=warning_msg,dr_suppression=dr_suppression,$
                          atmos_params=atmos_params)

     
     if err_calc ne '' then begin
; ERROR checking?
        print,' % CH_SYNTHETIC ERROR:  ',err_calc
        return
     endif

     if KEYWORD_SET(verbose) and n_elements(warning_msg) gt 1 then begin
        for im=1, n_elements(warning_msg)-1 do print,'% CH_SYNTHETIC: '+warning_msg[im]
     endif
     
  endif                         ; else we have the ioneq data already


;-------- END OF ADVANCED MODEL ----------------   


  IF KEYWORD_SET(verbose) THEN print,'% CH_SYNTHETIC: getting emissivities '

;;
;; Create a widget to inform about progress if progress is set
;;

  IF keyword_set(progress) THEN BEGIN

     pct = 0
     halt = 0
     pctage = 0.1
     lastpct = -pctage

     base = widget_base(/column,title='Progress', group=group)
     pct_slider_id = widget_slider(base,xsize = 500, maximum=100,minimum=0,$
                                   title='Approximate % done')
     interrupt_id = widget_button(base,value='Click here to STOP calculation')
     xrealize,base,/center,group=group

; id = progmeter(/INIT,label='Progress Meter',button='Click here to STOP calculation')

  END

  nlist=n_elements(mlist)

;
; ---------------------------------------
; *** MAIN INPUT AND CALCULATION
; *** LOOP THROUGH EACH ION IN MASTERLIST
; ---------------------------------------
;
  FOR ilist=0,nlist-1 DO BEGIN

     gname=mlist[ilist]
     convertname,gname,iz,ion
     ion2spectroscopic,gname,snote, dielectronic=dielectronic
     zion2spectroscopic,iz,ion,snote
     
;convert z and ionisation stage to filename 
;(eg z=26, ion=24 > !xuvtop/fe/fe_24 ) :
;-------------------------------------------

     zion2filename,iz,ion,fname, diel=dielectronic

;add a check:

     IF NOT file_exist(fname+'.elvlc') THEN BEGIN 
        err_msg ='%CH_SYNTHETIC: Error, no files in the database correspondent to '+gname 
        print,err_msg
        return
     END 

     wname=fname+'.wgfa'
     elvlcname=fname+'.elvlc'
     upsname=fname+'.scups'
     pname=fname+'.psplups'

     this_ioneq=ioneq[*,iz-1,ion-1+dielectronic]
     
; first check if ion exists in the ion fraction file.
;---------------------------------------------------

     ind_gioneq = where(this_ioneq GT  0.)
     IF ind_gioneq[0] NE  -1 THEN BEGIN 

;-------------------------------------------------------------------------------
;
; 't_index' serves different purposes depending on whether the isothermal
;      approximation is used or not:
;
;  1) If isothermal is not specified, we have two options:
;
;    1-the G(T)'s are calculated. 
;     't_index' is simply the index of the
;     temperatures from the Ion Fraction file where the values are not null,
;     if constant Density or Pressure are used. Otherwise it also checks the 
;      minimum and maximum Te values given as input.
;
;    2-The line intensities are calculated:
;     then the DEM is being used, and t_index contains the indices of 
;     ioneq_logt where the DEM distribution and ion fraction overlap,
;     if constant Density or Pressure are used. Otherwise it also checks the 
;      minimum and maximum Te values given as input.
;
;  2)  If isothermal is specified, then t_index will 
;      contain the index of isothermal as long the temperature(s) lie 
;      within the ion's ion balance range,
;      if constant Density or Pressure are used. Otherwise it also checks the 
;      minimum and maximum Te values given as input.

        IF n_elements(logt_isothermal) EQ 0 THEN BEGIN

           IF keyword_set(goft) THEN BEGIN 

              IF model_name EQ  'Function' THEN BEGIN  
                 t_index=WHERE(this_ioneq NE 0. AND $
                               (ioneq_logt GE min(alog10(temperature))) AND $
                               (ioneq_logt LE  max(alog10(temperature))) )
              ENDIF ELSE BEGIN
                                ;
                                ; PRY, 22-Jun-2023:
                                ; Modified t_index due to problems with Li-sequence ions.
                                ;
;               t_index=WHERE(this_ioneq NE 0.)
                 t_index=WHERE(this_ioneq GE this_ioneq/1e6)
              ENDELSE 

           ENDIF ELSE BEGIN 

              IF model_name EQ  'Function' THEN BEGIN 
                 t_index=WHERE((dem_int NE 0.) AND (this_ioneq NE 0.) AND $
                               (ioneq_logt GE min(alog10(temperature))) AND $
                               (ioneq_logt LE  max(alog10(temperature))) )
              ENDIF ELSE BEGIN
                                ;
                                ; PRY 26-Oct-2020: modified t_index
                                ;
                 ff=dem_int*this_ioneq*10.^ioneq_logt
                 t_index=where(ff GE max(ff)*1e-4)
;              t_index=WHERE((dem_int NE 0.) AND (this_ioneq NE 0.))
              ENDELSE 

           ENDELSE   

        ENDIF   ELSE BEGIN

           IF model_name EQ  'Function' THEN $ 
              t_index = where((logt_isothermal GE MIN(ioneq_logt[ind_gioneq])) AND $
                              (logt_isothermal LE MAX(ioneq_logt[ind_gioneq])) AND $
                              (logt_isothermal GE min(alog10(temperature))) AND $
                              (logt_isothermal  LE  max(alog10(temperature))) ) ELSE $ 
                                 t_index = where((logt_isothermal GE MIN(ioneq_logt[ind_gioneq])) AND $
                                                 (logt_isothermal LE MAX(ioneq_logt[ind_gioneq])))

        ENDELSE
        
        IF t_index[0] NE -1 THEN BEGIN

; GDZ for v.9:

; This loads up the ion's atomic data and rates, which are then
; directly passed to pop_solver. Note: the routine calls ch_setup_ion
; to reead the basic atomic data, but then additionally creates the
; rate arrays and also take into account of the new v.9 of calculating
; the satellite lines, by including the autoionising levels. 
; Note: the number of levels is reduced, if requested, in this
; routine, by passing n_lev

; - Routine checks if there are lines in the range wmin to wmax
;   and also if these lines have non-zero A-values.
           
; define  temperature array (TEMP) to calculate level populations
;-------------------------------------------------------------

           IF n_elements(logt_isothermal) EQ 0 THEN BEGIN
;this works for both the G(T) case and the intensity case:
              TEMP=10.^ioneq_logt[t_index]

           ENDIF ELSE BEGIN

              log_temp = logt_isothermal[t_index]
              TEMP = 10.^log_temp

; GDZ - advanced model change: 
              if  keyword_set(advanced_model) then begin 
                 ion_frac = this_ioneq[t_index]
                 
              endif else begin 
;do a spline interpolation in the logs:
                 ion_frac = spline(ioneq_logt[ind_gioneq],alog10(this_ioneq[ind_gioneq]),log_temp)
                 ion_frac = 10.^ion_frac
              end 
           ENDELSE 

           nt = n_elements(t_index) ; number of temperatures

; define  density array (DENS) to calculate level populations
;-------------------------------------------------------------

           
           DENS = fltarr(N_ELEMENTS(TEMP))

           IF model_name EQ  'Function' THEN BEGIN 
              IF  n_elements(model_file) gt 0  THEN BEGIN ; GDZ 
                 ltemperature=alog10(temperature)
                 ldensity=alog10(density)
                 i_uniq=uniq(ltemperature,sort(ltemperature))
                 xx=ltemperature[i_uniq]
                 yy=ldensity[i_uniq]
                 DENS = 10.^SPLINE(xx,yy, alog10(TEMP) )
              ENDIF ELSE BEGIN 
;this is the case of isothermal + an array of densities:
                 dens=density[t_index]
              ENDELSE 
           ENDIF   ELSE BEGIN 
              IF n_elements(density) NE 0 THEN  dens(*)=density ELSE dens=pressure/temp
           ENDELSE 

           nd=n_elements(dens)  ; number of densities 

;  calculate level populations
;------------------------------
           
           IF n_elements(dens) NE n_elements(temp) THEN BEGIN
              print,'****TEMP AND DENS MUST HAVE SAME SIZE!!****'
              STOP
           ENDIF 

                                ;
                                ; PRY, 12-May-2023
                                ; I've moved the population lookup table here so that I can create
                                ; "levmax" which is then input to ch_setup_ion. This is just in
                                ; case the lookup table has less levels than the wgfa structure.
                                ;
           junk=temporary(levmax) ; make sure levmax does not exist
           no_lookup=1-keyword_set(lookup)
           IF keyword_set(lookup) THEN BEGIN
              p=ch_lookup_table_interp(gname,dens,temp,/quiet,/pad)
              IF n_tags(p) EQ 0 THEN BEGIN
                 no_lookup=1
                 print,'% CH_SYNTHETIC: lookup tables not found. Using standard method...'
              ENDIF ELSE BEGIN 
                 pops=rearrange(p.pop,[2,1,3])
                 s=size(pops,/dim)
                 levmax=s[2]
                 pops2=dblarr(s[0],s[2])
                 FOR i=0,s[0]-1 DO pops2[i,*]=pops[i,i,*]
                 pops=temporary(pops2)
              ENDELSE
           ENDIF

; GDZ : fixed a bug,   ioneq_file needs to be passed the file name        
           input=ch_setup_ion(gname,rphot=rphot,radtemp=radtemp,noprot=noprot, $
                              ioneq_file=ioneq_name,abund_file=abund_file,path=path, $
                              quiet=quiet, $
                              wvlmin=wmin,wvlmax=wmax, index_wgfa=anylines, $
                              noionrec=noionrec, no_auto=no_auto , no_rrec=no_rrec, $
                              opt_lookup=1b-no_lookup,n_levels=levmax)
           
           
           IF anylines[0] EQ -1 THEN BEGIN 
              IF  keyword_set(verbose) THEN  $
                 print, 'No lines in the selected range for Ion '+gname+' !!!'
              msg =  [msg, 'No lines in the selected range for Ion '+gname+' !!!']
           ENDIF
           
        ENDIF    ELSE BEGIN
           IF  keyword_set(verbose) THEN  $
              print, 'No Temperature overlap for Ion '+gname+' !!!'
           msg =  [msg, 'No Temperature overlap for Ion '+gname+' !!!']
           anylines = -1
        ENDELSE



; end of checks 
;---------------------------------------------------------------------------


        IF anylines[0] NE -1 THEN BEGIN ; do this ion

; number of lines:
           nn=n_elements(anylines)

           pct = (ilist+1)/float(nlist) *100.

           IF NOT keyword_set(verbose) THEN $
              print,format='($,"Progress:  %",i3,a)',pct, $
                    string(13b)

           IF keyword_set(progress) THEN BEGIN

              info_progress, pct,lastpct,pctage, pct_slider_id,$
                             interrupt_id,halt,quiet, 'Calculating  '+gname+' ('+trim(nn)+' lines)',  group=group
              IF halt THEN GOTO,halt

           END 

           IF KEYWORD_SET(verbose) THEN print,'% CH_SYNTHETIC: calculating  '+gname+' ('+trim(nn)+' lines)'

                                ;
                                ; Extract items from the elvlc structure
                                ;
           ecm=input.elvlcstr.data.energy
           mult=float(input.elvlcstr.data.mult)
           res_ascii=input.elvlcstr.data.full_level
           res_latex=input.elvlcstr.data.full_level_latex
;         convert_terms_all,res_ascii,res_latex

                                ;
                                ; Extract items from wgfa structure
                                ;
           lvl1=input.wgfastr.lvl1
           lvl2=input.wgfastr.lvl2
           wvl1=input.wgfastr.wvl
           a_value1=input.wgfastr.aval
           diel=input.wgfastr.diel

                                ;
                                ; Extract electron collision structure
                                ;
           splstr=input.splstr


                                ;
                                ; For dielectronic ions, the ground level is actually the
                                ; ground level of the recombining ion, so we need to read the
                                ; elvlc file of this ion and get the J value for this level
                                ; and insert it into INPUT. Note that INPUT.JJ is used by
                                ; pop_solver. 
                                ;
           IF dielectronic EQ 1 THEN BEGIN
              zion2filename,iz,ion+1,dname
              dlvlcname=dname+'.elvlc'
              read_elvlc,dlvlcname,elvlc=dlvlc
              input.jj[0]=dlvlc.data[0].j
           ENDIF



; GDZ - Oct 2018, the following has been modified for  v9:
; PRY, 11-Jun-2020: introduced /lookup option. Note I use the /pad
; keyword for ch_lookup_table_interp. Some ions have levels with zero
; population, which get omitted from the lookup tables, thus the
; population array ends up a different size to the pop_solver
; array. The /pad keyword makes them the same size.
; PRY, 12-May-2023: the lookup section has been moved earlier.
; GDZ, added verbose          
;---------------------------------------------------------

           IF keyword_set(no_lookup) THEN BEGIN
              pop_solver, input,temp,dens,pops,/pressure,radfunc=radfunc, frac_cutoff=frac_cutoff, $
                          regular=regular, sparse=sparse, lapack=lapack,verbose=verbose
           ENDIF 


; **** note that the pops  array has different dimensions depending on
; the input:  dblarr(nt,nd,nlev)  


; Removes lines with zero intensity and excludes lines that have already appeared in the
; standard CHIANTI file turning up in the dielectronic ion spectrum as well. Skips the
; ion if there are no lines that survive these tests.


           i_neg=where(total(pops[*, lvl2[anylines]],1) eq 0.) 

           IF (n_elements(i_neg) lt n_elements(anylines)) or $
              (n_elements(i_neg) eq n_elements(anylines) $
               and i_neg(0) lt 0) THEN BEGIN ; if there are lines left in the ion

              IF i_neg(0) ge 0 THEN begin
                 IF  keyword_set(verbose) THEN  $
                    print,gname+' removing lines from zero level population ...'
                 remove,i_neg,anylines
              end 


              
; PRY, 12-Apr-2021
; The following code is a leftover from the old dielectronic "d"
; ions. I've modified it so the code is only executed if
; dielectronic=1 as unexpected behavior occurred sometimes.
              
; Note that the removes any transitions for which the upper level is below
; the lowest level in the scups data-set. For normal ions no ions
; should  be flagged, but for the old "d" ions it will be all transitions
; below the ionization limit. I don't think this is correct as
; cascading from the AI levels are a legitimate population contributor
; to the bound levels. However, the code will be left as-is. The "d"
; ions are being phased out in favor of the new 2-ion models.
;
              IF keyword_set(dielectronic) THEN BEGIN 
                 
                 i_spl=where(lvl2[anylines] lt min(splstr.data.lvl2) or lvl2[anylines] gt max(splstr.data.lvl2))
                 IF (n_elements(i_spl) lt n_elements(anylines)) or $
                    (n_elements(i_spl) eq n_elements(anylines) and i_spl(0) lt 0) THEN BEGIN 
                    IF i_spl(0) ge 0  and dielectronic THEN begin
                       IF  keyword_set(verbose) THEN  $
                          print,gname+': removing '+trim(i_spl)+' lines... !!! '
                       remove,i_spl,anylines
                    ENDIF
                 ENDIF 
              ENDIF 

              nn=n_elements(anylines) ; new number of lines

; Calculates power, either in phot or in erg, for each line

              IF KEYWORD_SET(photons) THEN de = 1. $
              ELSE de = 1.986e-8/ABS(wvl1[anylines])

; GDZ: this should work for both constant density or pressure,
; isothermal, etc. : 

              de_nj_Aji=dblarr(nt,nn)
              FOR i=0,nt-1 DO de_nj_Aji[i,*]=pops[i, lvl2[anylines]-1]*de*a_value1[anylines]


; Calculates intensities or contribution functions, either isothermal or not

              IF n_elements(logt_isothermal) EQ 0 THEN BEGIN ; Non-isothermal case

                 IF keyword_set(goft) THEN BEGIN ; Contribution functions

                    tgt = 0.*de_nj_Aji
                    temp_temp = 0.*de_nj_Aji
                    FOR i=0,nt-1 DO BEGIN
                       tgt[i,*] = de_nj_Aji[i,*]*this_ioneq[t_index[i]]/dens[i]
                       temp_temp[i,*] = temp(i)
                    ENDFOR
                    tgt = tgt/4./!pi
                    get_tmax = MAX(tgt,tmax_ind,dim=1)
                    tmax=alog10(temp_temp[tmax_ind])

                    this_goft = dblarr(n_ioneq_logt,nn)  ;; PRY, 22-Aug-2008
                    this_goft[t_index,*] = tgt[*,*]

                 ENDIF ELSE BEGIN ; Intensities

; GDZ - ADVANCED model: the standard CHIANTI is a uniform grid of log T values but with the advanced
;  model the grid does not need to be uniform.                  
; dlnt=ALOG(10.^(ioneq_logt[1]-ioneq_logt[0]))      
                    
                    dlnt=fltarr(n_elements(ioneq_logt))
                    for kk=0, n_elements(ioneq_logt)-2 do $
                       dlnt[kk]=ALOG(10.^(ioneq_logt[kk+1]-ioneq_logt[kk]))
                                ; fix the last one:
                    dlnt[n_elements(ioneq_logt)-1]=dlnt[n_elements(ioneq_logt)-2]
                    
                    
                    dt=temp * dlnt

                    ieq_dem_dt=this_ioneq[t_index]*dem_int[t_index]*dt/dens
                    intensity=transpose(de_nj_aji)#ieq_dem_dt
                    intensity=intensity/4./!pi
                    tgt = 0.*de_nj_Aji
                    temp_temp = 0.*de_nj_Aji
                    FOR i=0,nt-1 DO BEGIN 
                       tgt[i,*] = de_nj_Aji[i,*]*ieq_dem_dt[i]
                       temp_temp[i,*] = temp(i)
                    ENDFOR
                    get_tmax = MAX(tgt,tmax_ind,dim=1)
                    tmax=alog10(temp_temp[tmax_ind])

                 END 

              ENDIF  ELSE BEGIN ; Isothermal case

                 ieq_dem_dt= ion_frac * $
                             double(10.)^(logem_isothermal[t_index]) /dens 

                 IF keyword_set(no_sum_int) THEN BEGIN


                    intensity = dblarr(n_elements(logt_isothermal),nn)

                    t_temp = 0.*intensity
                    FOR i=0,nt-1 DO BEGIN
;GDZ-corrected - need to keep t_index
                       intensity[t_index[i],*]=de_nj_Aji[i,*]*ieq_dem_dt[i]
                       t_temp[t_index[i],*]=logt_isothermal[i]
                    ENDFOR
                    get_tmax=max(intensity,tmax_ind,dim=1)
                    tmax=t_temp[tmax_ind]
                 ENDIF ELSE BEGIN
                    intensity=transpose(de_nj_Aji)#ieq_dem_dt
                    tmax=alog10(temp(0))+fltarr(n_elements(intensity))
                 ENDELSE

                 intensity=intensity/4./!pi


              ENDELSE 


;; isort=sort(wvl1[anylines[*]])
;; for ii=0, nn-1 do print, a_value1[anylines[isort[ii]]], wvl1[anylines[isort[ii]]],intensity[isort[ii]]

              lbl1: 
              
              IF keyword_set(goft) THEN BEGIN 

                 t_list_goft = this_goft

              ENDIF  ELSE BEGIN 
                 t_list_int=double(intensity)
              END 

;store the absolute value of the wavelengths

              t_list_wvl=FLOAT(abs(wvl1[anylines]))

;flag any lines with theoretical wavelength with a -1

              t_list_flag=fltarr(nn)
              negative_wvl=where(wvl1[anylines] lt 0)
              IF negative_wvl[0] ge 0 THEN t_list_flag(negative_wvl) = -1
              
              t_list_l1 = lvl1[anylines]
              t_list_l2 = lvl2[anylines]
              t_list_iz = iz+intarr(nn)
              t_list_ion = ion+intarr(nn)
              t_list_snote = replicate(snote,nn)
              t_list_tmax = tmax
                                ;
                                ; PRY, 13-Jun-2023
                                ; The following flags the satellite lines in snote by making
                                ; use of the diel tag of wgfastr.
                                ;
              t_list_diel=diel[anylines]
              k=where(t_list_diel EQ 1,nk)
              IF nk NE 0 THEN t_list_snote[k]=t_list_snote[k]+' s'

              IF NOT keyword_set(goft) THEN BEGIN
                 IF keyword_set(no_sum_int) AND n_elements(logt_isothermal) NE 0 THEN BEGIN
                    nn_exist=n_elements(reform(list_int[0,*]))
                    nn_new=nn+nn_exist

                    list_int_new=dblarr(n_elements(logt_isothermal),nn_new)
;GDZ-fixed bug.
                    FOR i=0,n_elements(logt_isothermal)-1 DO BEGIN
                       list_int_new[i,0:nn_exist-1]=list_int[i,*]    ; Intensity of previous lines stays the same
                       list_int_new[i,nn_exist:nn_new-1]=t_list_int[i,*] ; Stores the new ones

                    ENDFOR
                    
                    list_int=list_int_new

                 ENDIF  ELSE BEGIN
                    list_int=[list_int,t_list_int]
                 ENDELSE
              ENDIF ELSE BEGIN
                 nn_exist=n_elements(reform(list_goft[0,*]))
                 nn_new=nn+nn_exist
                 list_goft_new=dblarr(n_ioneq_logt,nn_new)  ;; PRY, 22-Aug-2008
                 FOR i=0,n_ioneq_logt-1 DO BEGIN 
                    list_goft_new[i,0:nn_exist-1]=list_goft[i,*]
                    list_goft_new[i,nn_exist:nn_new-1]=t_list_goft[i,*]
                 ENDFOR
                 list_goft=list_goft_new
              ENDELSE
              list_wvl=[list_wvl,t_list_wvl]
              list_flag=[list_flag,t_list_flag]
              list_l1=[list_l1,t_list_l1]
              list_l2=[list_l2,t_list_l2]
              list_iz=[list_iz,t_list_iz]
              list_ion=[list_ion,t_list_ion]
              list_snote=[list_snote,t_list_snote]
              list_tmax=[list_tmax,t_list_tmax]
              list_ident=[list_ident,res_ascii(t_list_l1-1)+' - '+res_ascii(t_list_l2-1)]
              list_ident_latex=[list_ident_latex,res_latex(t_list_l1-1)+' - '+res_latex(t_list_l2-1)]
                                ;
;          ENDIF                  ;  if there are lines left after zero population check 
           ENDIF                ;  if there are lines left after dielectronic check
        ENDIF                   ;  if block for anylines

     ENDIF  ELSE BEGIN          ;  ion present in ion fraction file                                ;
        IF  keyword_set(verbose) THEN  $
           print, 'Ion '+gname+' not present in the Ion. Frac. file!!!'
        msg =  [msg, 'Ion '+gname+' not present in the Ion. Frac. file!!!']
     END 

  ENDFOR                        ;  reading masterlist.ions


halt:

  IF keyword_set(progress) THEN BEGIN

     IF n_elements(interrupt_id) EQ 1 THEN $
        ev = widget_event(interrupt_id,/nowait,save_hourglass=0)

     xkill,base

;status = progmeter(id,/DESTROY)

  END 


  IF n_elements(list_iz) EQ 1 THEN BEGIN 
     err_msg = '%CH_SYNTHETIC: No lines in the selected wavelength range ! -- EXIT '
     print,err_msg
     return
  ENDIF ELSE msg = [msg, trim( n_elements(list_iz)-1)+' lines calculated and stored into memory !' ]

  list_wvl=list_wvl[1:*]
  list_flag = list_flag[1:*]
  list_ident = list_ident[1:*]
  list_ident_latex = list_ident_latex[1:*]
  list_snote =list_snote[1:*]
  list_l1 = list_l1[1:*]
  list_l2 = list_l2[1:*]
  list_tmax = list_tmax[1:*]
  list_iz = list_iz[1:*]
  list_ion = list_ion[1:*]


  IF NOT keyword_set(goft) THEN BEGIN
     IF keyword_set(no_sum_int) AND n_elements(logt_isothermal) NE 0 THEN BEGIN

        list_int=list_int[*,1:*]
     ENDIF ELSE BEGIN 
        list_int=list_int[1:*]
     ENDELSE
  ENDIF ELSE BEGIN
     list_goft=list_goft[*,1:*]
  ENDELSE

; SORT IN WAVELENGTH ???
; i_sort = sort(list_wvl)


;
; This rest of the routine sets up the output structure 
;

  n = N_ELEMENTS(list_wvl)

  IF NOT keyword_set(goft) THEN BEGIN
     IF NOT keyword_set(no_sum_int) THEN BEGIN
        str = {iz: 0     ,$
               ion: 0    ,$
               ident: '', $
               ident_latex: '', $
               snote:'', $
               lvl1: 0   ,$
               lvl2: 0   ,$
               tmax: 0.  ,$
               wvl:  0.d ,$
               flag: 0, $
               int:  0.d  } 
     ENDIF ELSE BEGIN
        str = {iz: 0     ,$
               ion: 0    ,$
               ident: '', $
               ident_latex: '', $
               snote:'', $
               lvl1: 0   ,$
               lvl2: 0   ,$
               tmax: 0.  ,$
               wvl:  0.d ,$
               flag: 0, $
               int:  dblarr(n_elements(logt_isothermal)) }
     ENDELSE
  ENDIF ELSE BEGIN
     str = {iz: 0     ,$
            ion: 0    ,$
            ident: '', $
            ident_latex: '', $
            snote:'', $
            lvl1: 0   ,$
            lvl2: 0   ,$
            tmax: 0.  ,$
            wvl:  0.d ,$
            flag: 0, $
            goft: dblarr(n_ioneq_logt) }     ;; PRY, 22-Aug-2008
  ENDELSE

  str2 = REPLICATE(str,n)

;logt_limits = 0.
;IF n_elements(logt_isothermal) EQ 0 THEN $
;logt_limits = [MIN(ioneq_logt[gdt]),MAX(ioneq_logt[gdt])]

  version = ' '
  ff = findfile(concat_dir(!xuvtop,'VERSION'))
  IF  ff(0)  NE ''  THEN BEGIN 
     openr, 1, ff(0)
     readf, 1, version
     close, 1
  ENDIF ELSE BEGIN 
     print,  'Please update your CHIANTI database version, which is older than 4.0'
     version = ' '
  END 

; reset for output
  if n_elements(model_file) eq 0 then model_file=' '

; ioneq_name must be defined by now. 
  IF n_elements(ioneq_name) EQ 0 THEN ioneq_name=''
  IF n_elements(ioneq_ref) EQ 0 THEN ioneq_ref=''


  OUTPUT = {lines:str2, $
            ioneq_logt:ioneq_logt,ioneq_name:ioneq_name, ioneq_ref:ioneq_ref, $
            wvl_limits:  [wmin,wmax], $
            model_file:model_file, model_name:model_name, model_ne:model_ne,$
            model_pe:model_pe, model_te:model_te, $
            wvl_units: wvlunits, $
            int_units: intunits, $
            add_protons:add_protons, $
            date: systime(), $
            version:version, $
            lookup: 1b-no_lookup, $
            photoexcitation:photoexcitation}

  IF photoexcitation THEN begin

     if n_elements(RADFUNC) eq 1 THEN $
        OUTPUT =CREATE_STRUCT(OUTPUT, 'rphot', rphot, 'radfunc', radfunc) $
     else OUTPUT =CREATE_STRUCT(OUTPUT, 'rphot', rphot, 'radtemp', radtemp)

  endif


  IF (NOT  keyword_set(goft)) THEN BEGIN 
     IF (n_elements(logt_isothermal) EQ   0)  THEN $
        OUTPUT =CREATE_STRUCT(OUTPUT, 'dem_name', dem_name , 'dem_ref', dem_ref, $
                              'dem_logt', dem_logt, 'dem', dem) ELSE BEGIN 
        OUTPUT =CREATE_STRUCT(OUTPUT, 'logt_isothermal', logt_isothermal)
        OUTPUT =CREATE_STRUCT(OUTPUT, 'logem_isothermal', logem_isothermal)
     END 
  END 

  OUTPUT.lines.iz = list_iz
  OUTPUT.lines.ion = list_ion
  OUTPUT.lines.ident = list_ident
  OUTPUT.lines.ident_latex = list_ident_latex
  OUTPUT.lines.snote = list_snote

  OUTPUT.lines.wvl = list_wvl
  OUTPUT.lines.flag = list_flag

  OUTPUT.lines.lvl1 = list_l1
  OUTPUT.lines.lvl2 = list_l2
  OUTPUT.lines.tmax = list_tmax

  IF NOT keyword_set(goft) THEN $ 
     OUTPUT.lines.int = list_int ELSE BEGIN 
     OUTPUT.lines.goft = list_goft
  END 

  IF n_elements(save_file) NE 0  THEN BEGIN 
     savegen, file=save_file, struct=OUTPUT
     IF KEYWORD_SET(verbose) THEN print,' Intensities stored in the file '+save_file
  END 

  if n_elements(msg) gt 1 then  msg=msg[1:*]

  t2 = systime(1)
  print,format='("% CH_SYNTHETIC: Line intensities computed in ",f8.1," seconds")',t2-t1

END
