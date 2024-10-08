
FUNCTION ch_lookup_gofnt, ionname, wmin=wmin, wmax=wmax, log_temp=log_temp, $
                          log_dens=log_dens, log_press=log_press, $
                          abund_file=abund_file, noabund=noabund, $
                          photons=photons, dir_lookup=dir_lookup, $
                          lower_levels=lower_levels, $
                          upper_levels=upper_levels, $
                          aval_struc=aval_struc, ioneq_file=ioneq_file, $
                          advanced_model=advanced_model

;+
; NAME:
;      CH_LOOKUP_GOFNT
;
; PURPOSE:
;      This routine takes a population lookup table and computes the
;      contribution function for the emission line. The definition of
;      the contribution is the same as that of the standard CHIANTI
;      routine gofnt.pro.
;
; INPUTS:
;      IonName: The name of an ion in CHIANTI format (e.g., 'o_6').
;
; OPTIONAL INPUTS:
;      Wmin:     If set, then only wavelengths above WMIN will be
;                listed in the widget.
;      Wmax:     If set, then only wavelengths below WMAX will be
;                listed in the widget.
;      Log_Temp: A 1D array of Log10 temperatures. If not specified,
;                then the temperature array in the lookup table will be used.
;      Log_Dens: Log10 of electron number density. If a scalar is given,
;                then the density will apply to all of the temperatures.
;                If an array of same size as log_temp is given, then the
;                population is calculated at each temperature with the
;                matching density. Either log_dens or log_press should be
;                defined. 
;      Log_Press: Log10 of electron pressure (N_e*T). If a scalar is given,
;                then the pressure will apply to all of the temperatures.
;                If an array of same size as log_temp is given, then the
;                population is calculated at each temperature with the
;                matching pressure. Either log_dens or log_press should be
;                defined. 
;      Abund_File: The name of  CHIANTI format element abundance
;                 file. If not specified then the user is asked to
;                 choose a file with a widget.
;      Dir_Lookup: The name of a directory containing population
;                 lookup table files. If set, then the routine will
;                 expect to find a filename of the form
;                 'pop_lookup_[ion_name].txt' in this directory.
;      Lower_Levels: Specify the lower level of the transition for
;                  which the contribution function is required. Can be
;                  an array. Upper_levels needs to be specified as
;                  well.
;      Upper_Levels: Specify the upper level of the transition for
;                  which the contribution function is required. Can be
;                  an array. Lower_levels needs to be specified as
;                  well.
;      Aval_Struc: A structure containing information about the
;                 transition(s) for which the contribution function is
;                 required. It must have the tags:
;                  .lvl1  Index of lower level
;                  .lvl2  Index of upper level
;                  .wvl   Wavelength of transition
;                  .aval  A-value for transition
;                 If there are multiple transitions, then AVAL_STRUC
;                 will be an array. Note that using AVAL_STRUC is the
;                 quickest way to calculate the contribution function
;                 as there is no need to read the wgfa file. The
;                 easiest way to create the structure is to use
;                 ch_wgfa_select. 
;      Ioneq_File: The name of an ionization equilibrium file. The
;                  default is to use the CHIANTI file (!ioneq_file).
;
; KEYWORD PARAMETERS:
;      NOABUND: If set, then the contribution function is not
;               multiplied by the element abundance.
;      PHOTONS: If set, then the contribution function is in photons
;               units rather than ergs.
;      ADVANCED_MODEL: If set uses the advanced models to compute the
;               ionization fraction.
;
; OUTPUTS:
;      An IDL structure with the tags:
;      .ltemp  Log temperatures at which contrib. fn. defined.
;      .gofnt  The contribution function in erg cm^3 s^-1 sr^-1.
;      .ldens  Log densities at which contrib. fn. defined. Either a
;              scalar or a 1D array (if log_press defined).
;      .lpress Log pressure at which contrib. fn. defined. Set to -1
;              if log_press not specified.
;      .abund_file  Name of the element abundance file.
;      .ioneq_file  Name of the ionization balance file; set to empty
;                   string if /advanced_model set.
;      .chianti_version  The version of CHIANTI that was used.
;      .ion    Name of the ion.
;      .advanced_model  Contains value of advanced_model keyword.
;      .time_stamp  The time at which the structure was created.
;
; CALLS:
;      CONVERTNAME, CH_WGFA_SELECT, ZION2FILENAME, READ_WGFA_STR
;      CH_READ_POP_LOOKUP_TABLE, READ_ABUND, READ_IONEQ, GET_IEQ,
;      CH_LOOKUP_TABLE_INTERP 
;
; EXAMPLES:
;      Choose the transition(s) using a widget:
;       IDL> g=ch_lookup_gofnt('o_6',log_dens=9.0)
;
;      Reduce the wavelength range for the widget display:
;       IDL> g=ch_lookup_gofnt('o_6',wmin=1030,wmax=1040,log_press=15.0)
;
;      Directly specify a transition using indices (no widget):
;       IDL> g=ch_lookup_gofnt('o_6',lower=1,upper=3,log_dens=9.0)
;
;      Use the aval_struc input (by first calling ch_wgfa_select):
;       IDL> aval_struc=ch_wgfa_select('o_6')
;       IDL> g=ch_lookup_gofnt('o_6',aval_struc=aval_struc,log_press=15.0)
; 
; MODIFICATION HISTORY:
;      Ver.1, 18-Dec-2019, Peter Young
;      Ver.2, 22-Jun-2022, Peter Young
;        The absolute value of the wavelength is used to compute the
;        transition energy to prevent negative G(T) values.
;      Ver.3, 27-Sep-2023, Peter Young
;        Modified to allow density and pressure to be arrays, but they
;        must match the size of the temperature array.
;      Ver.4, 26-Sep-2024, Peter Young
;        Added /advanced_model keyword.
;-


IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> g=ch_lookup_gofnt( ion_name ) '
  print,''
  print,' Keyword parameters:'
  print,'   wmin=, wmax=, log_temp=, log_dens=, log_press=, abund_file=, /noabund,'
  print,'   /photons, lower_levels=, upper_levels, aval_struc=, ioneq_file= dir_lookup='
  print,'   /advanced_model'
  print,''
  return,-1
ENDIF 


;
; Do checks on the level inputs
;
nl=n_elements(lower_levels)
nu=n_elements(upper_levels)
IF nl NE 0 AND nl NE nu THEN BEGIN
  print,'% CH_LOOKUP_GOFNT: the input LOWER_LEVELS must have the same size as UPPER_LEVELS. Returning...'
  return,-1
ENDIF

convertname,ionname,iz,iion

;
; There are three options for choosing the transition(s) for which the
; contribution function is required:
;   (1) using aval_struc (contains indices, wavelengths and A-values).
;   (2) using a widget (the default)
;   (3) specifying the indices of the transitions (lower_levels,
;       upper_levels)
; Each of these leads to arrays lower_levels, upper_levels, wvl and
; aval that are used later in the routine. The number of transitions
; is given by NU.
;
IF n_tags(aval_struc) EQ 0 THEN BEGIN
  IF nu EQ 0 THEN BEGIN
    wgfastr=ch_wgfa_select(ionname,wmin=wmin,wmax=wmax)
    IF n_tags(wgfastr) EQ 0 THEN BEGIN
      print,'% CH_LOOKUP_GOFNT: No lines selected. Returning...'
      return,-1
    ENDIF
  ENDIF ELSE BEGIN
    zion2filename,iz,iion,fname
    read_wgfa_str,fname+'.wgfa',wgfastr,ref
    index=lonarr(nu)
    FOR i=0,nu-1 DO BEGIN
      j=where(wgfastr.lvl1 EQ lower_levels[i] AND wgfastr.lvl2 EQ upper_levels[i])
      index[i]=j[0]
    ENDFOR
    k=where(index NE -1,nk)
    IF nk EQ 0 THEN BEGIN
      print,'% CH_LOOKUP_GOFNT: the levels specified by LOWER_LEVELS and UPPER_LEVELS do not exist! Returning...'
      return,-1
    ENDIF
    index=index[k]
    wgfastr=wgfastr[index]
  ENDELSE
 ;
  nu=n_elements(wgfastr)
  lower_levels=wgfastr.lvl1
  upper_levels=wgfastr.lvl2
  wvl=wgfastr.wvl
  aval=wgfastr.aval
ENDIF ELSE BEGIN
  nu=n_elements(aval_struc)
  lower_levels=aval_struc.lvl1
  upper_levels=aval_struc.lvl2
  wvl=aval_struc.wvl
  aval=aval_struc.aval
ENDELSE 


;
; If log_temp has not been specified, then I read the lookup table to
; find the temperature range at which it is defined, and then I use
; this. The lookup table is deleted.
;
IF n_elements(log_temp) EQ 0 THEN BEGIN
  lookup_filename=ch_lookup_filename(ionname,dir_lookup=dir_lookup,status=status)
  IF status EQ 0 THEN BEGIN
    print,'% CH_LOOKUP_GOFNT: the lookup table was not found. Returning...'
    return,-1
  ENDIF 
  a=ch_read_pop_lookup_table(lookup_filename)
  log_temp=a.ltemp
  junk=temporary(a)
ENDIF
nt=n_elements(log_temp)

nd=n_elements(log_dens) 
np=n_elements(log_press)

IF (nd NE 0 AND np NE 0) OR (nd EQ 0 AND np EQ 0) THEN BEGIN
  message,/cont,/info,'please specify EITHER log_dens OR log_press. Returning...'
  return,-1
ENDIF

  

;
; ch_lookup_table_interp takes only dens as an input, so convert pressure to dens.
;
IF np NE 0 THEN ldens=log_press-log_temp ELSE ldens=log_dens
nd=n_elements(log_dens) 

IF nd GT 1 AND nd NE nt THEN BEGIN
  message,/info,/cont,'If pressure or density are given as an array, they must have the same size as the temperature array. Returning...'
  return,-1
ENDIF 

IF np NE 0 OR nd GT 1 THEN diag_swtch=1b else diag_swtch=0b
 
;
; If /noabund has not been set, then read the abundance file. Note
; that is abund_file is an empty string then read_abund will ask the
; user to choose a file.
;
IF NOT keyword_set(noabund) THEN BEGIN
  IF n_elements(abund_file) EQ 0 THEN abund_file=''
  read_abund,abund_file,ab,ref
ENDIF ELSE BEGIN
  abund_file=''
ENDELSE 

;
; Define the output array.
;
contrib_data=dblarr(nt)

;
; Derive contribution functions from the lookup tables.
;
convertname,ionname,iz,iion

;
; Read the ioneq file and extract the ion fractions.
;
IF keyword_set(advanced_model) THEN BEGIN
  IF np GT 0 THEN press=10.^log_press
  IF nd GT 0 THEN dens=10.^log_dens
  ioneq_data=ch_calc_ioneq(10.^log_temp, dens=dens, press=press, /adv, ele=iz)
  frac=reform(ioneq_data[*,iz-1,iion-1])
  ioneq_file=''
ENDIF ELSE BEGIN
  IF n_elements(ioneq_file) EQ 0 THEN ioneq_file=!ioneq_file
  read_ioneq,ioneq_file,ioneq_t,ioneq_frac,ref
  frac=get_ieq(10.^log_temp,iz,iion,ioneq_logt=ioneq_t,ioneq_frac=ioneq_frac)
ENDELSE 


;
; Obtain the level populations for the specified density and
; temperatures. 
;
popstr=ch_lookup_table_interp(ionname,10.^ldens, $
                              10.^log_temp,/quiet,dir_lookup=dir_lookup)
IF n_tags(popstr) EQ 0 THEN BEGIN
  print,'% CH_LOOKUP_GOFNT: problem found with lookup table interpolation. Returning...'
  return,-1
ENDIF 


;
;
; Go through each transition for the ion and compute the contribution
; function and sum.
;
FOR j=0,nu-1 DO BEGIN
  l1=lower_levels[j]
  l2=upper_levels[j]
 ;
  ilev=where(popstr.levels EQ l2)
  IF ilev[0] NE -1 THEN BEGIN 
    pop=reform(popstr.pop[*,*,ilev])
   ;
    IF NOT keyword_set(photons) THEN energy=1.986e-8/abs(wvl[j]) ELSE energy=1.0
   ;
   ; If log_press was specified, then pop will be a 2D array and we need
   ; to take the diagonal of this array. Otherwise pop will be a 1D
   ; array. 
   ;
    IF diag_swtch THEN pop=diag_matrix(pop)
    contrib_data=contrib_data + energy*pop*frac*aval[j]/10.^ldens/4./!pi
  ENDIF ELSE BEGIN
    print,'% CH_LOOKUP_GOFNT: Upper level '+trim(l2)+' (wavelength: '+trim(string(wvl[j],format='(f15.2)'))+') not found in lookup table.'
    print,'                   Please check your lookup table.'
    return,-1
  ENDELSE 
ENDFOR
IF NOT keyword_set(noabund) THEN contrib_data=contrib_data*ab[iz-1]


IF n_elements(log_press) EQ 0 THEN lpress=-1. ELSE lpress=log_press

output={ ltemp: log_temp, $
         gofnt: contrib_data, $
         ldens: ldens, $
         lpress: lpress, $
         abund_file: abund_file, $
         ioneq_file: ioneq_file, $
         chianti_version: popstr.chianti_version, $
         ion: trim(ionname), $
         advanced_model: keyword_set(advanced_model), $
         time_stamp: systime() }


return, output

END
