;+
;
; PROJECT:  CHIANTI
;
;        This program was developed as part of CHIANTI-VIP. 
;        CHIANTI-VIP (Virtual IDL and Python) is  a member of the CHIANTI family
;        mantained by Giulio Del Zanna, to develop additional features and
;        provide them to the astrophysics community.
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the
;       University of Cambridge, Goddard Space Flight Center, and University of Michigan.
;
;
; NAME:
;	
;       CH_ADV_MODEL_SETUP_GSK
;
;
; PURPOSE:
;
;	      Sets up the parameters needed for the advanced models.
;
;
; EXPLANATION
;
;       Primarily it reads the list of ions included in the advanced models, and the number of
;       levels included in the level population solution.
;
;       It optionally obtains the model atmosphere parameters required for the charge
;       transfer (CT) calculation. ATMOSPHERE_FILE provides the model atmosphere file, which
;       contains the temperature, density, pressure, total hydrogen density, H and optionally He
;       ion fractions. If He fractions are not provided, it obtains these from the default
;       CHIANTI ion fraction file. If He abundance is not provided relative to H, it obtains
;       this from the default abundance file.
;  
;       The routine then interpolates the parameters over the temperature grid inputted, which
;       is the grid used for the ion balance calculation. Because the interpolation is over
;       temperature the model atmosphere parameters will be truncated below and above the
;       temperature minimum and maximum in the file.
;
;       Alternately, if ATMOS_PARAMS is provided, then user provided values of the neutral/ion
;       fractions are used instead of trying to read an atmosphere_file. This is mostly intended
;       as part of the process of building large lookup tables.
;
; CATEGORY
;
;       CHIANTI; model atmosphere data
;
;
; CALLING SEQUENCE:
;
;       temp=10.^(findgen(61)*0.05+3.5)
;       params=ch_adv_model_setup(temp)
;       params=ch_adv_model_setup(temp,/ct,atm=!xuvtop+ $
;           '/ancillary_data/advanced_models/model_atmospheres/avrett_atmosphere.dat')
;
;
; INPUTS:
;
;       TEMP:  1D array specifying the temperatures (in K) over which the ion balances
;              are needed.
;
;	
; KEYWORDS:
;
;       CT:  extract the model atmosphere parameters required for calculating CT rates
;
;
; OPTIONAL INPUTS:
;
;       ATMOSPHERE_FILE:  A string of the file location and name for the model atmosphere
;                    from which the parameters will be extracted for CT rates
;
;       HE_ABUND:  Floating point number for the He abundance relative to hydrogen
;
;       ATMOS_PARAMS: This is a structure containing the atmospheric parameters, and is
;                     an alternative to giving ATMOSPHERE_FILE. The tags are:
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
; OUTPUTS:
;
;       IONS_NLEVELS:  An array of integers indicating the maximum number of levels
;             to be included in the level population calculations for each ion in the
;             advanced models.
;
;       MODEL_IONS:  The list of ions to be included in the advanced models
;
;
; OPTIONAL OUTPUTS:
;
;       CT_MODEL:  Structure containing all the parameters from the model atmosphere file
;             needed for calculating CT rates
;
;
; CALLS:
;
;       ch_read_list_ions, rdfloat, read_ioneq, read_abund
;
;
; PREVIOUS HISTORY:
;
;       NONE
;
;
; WRITTEN:
;         
;       v.1 Roger Dufresne (RPD) and Giulio Del Zanna (GDZ)
;       DAMTP, University of Cambridge, 16 Sept 2023
;
;
; MODIFIED:
;
;       v.2, 17 Oct 2023 GDZ, added option to calculate populations only
;               with a number of levels.
;
;       v.3, 30 Oct 2023   RPD, interpolate linearly in the log.
;
;       v.4, 07 Dec 2023   RPD, remove coronal double-valued temperatures in atmosphere file
;
;       v.5, 21 Mar 2024   RPD, changed names of model list and recombination files
;
;       v.6, 01 May 2024,  RPD, change atmosphere file keyword
;
;       v.7, 28 Aug 2024,  RPD, corrected error in setting He ion fractions to zero
;                               outside atmosphere temperature limits
;
;       v.8, 29 Aug 2024,  RPD, no longer read and return the recombination fitting coefficients
;
;       v.9, 24 Apr 2025,  Graham Kerr & Peter Young
;                          added atmos_params functionality, providing required values    
;                          directly without needing to read an atmosphere (intended for 
;                          use in producing large lookup tables). The routine now calls
;                          ch_read_atmos and ch_interp_atmos
;
;       v.10, 30-Apr-2025, Peter Young
;                          The /spline option is used when interpolating the helium ion
;                          fraction data; force the ion fractions to be >= 0 (no negative
;                          values).
;
; VERSION:  10
;
;- 

function ch_adv_model_setup,temp,ct=ct,atmosphere_file=atmosphere_file,he_abund=he_abund, $
                                atmos_params=atmos_params

zlabl=['H','He','Li','Be','B','C','N','O','F','Ne','Na',$
       'Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr', $
       'Mn','Fe','Co','Ni','Cu','Zn']


; find ions to be included in advanced models and number of levels to be used in level
; population calculation

ion_list=concat_dir(concat_dir(concat_dir(!xuvtop,'ancillary_data'),'advanced_models'),$
  'advmodel_list.ions')
if not file_exist(ion_list) then $
  message,'% CH_ADV_MODEL_SETUP Error: no list of ions to be included in transition region models'

info=ch_read_list_ions(ion_list)
model_ions=info.list_ions
ions_nlevels=info.nlevels

ntemp=n_elements(temp)

;
; Check atmos_params. If the tags only have one element (swtch=1), then the H and He
; data apply to all of the temperatures TEMP. Otherwise (swtch=0), the H and He data
; are interpolated onto TEMP values.
;
swtch=0b
ntags=n_tags(atmos_params)
IF keyword_set(ct) AND ntags NE 0 THEN BEGIN
  IF n_elements(atmos_params.h_elec) NE n_elements(atmos_params.h1_frac) THEN BEGIN
    message,/info,/cont,'The h_elec and h1_frac tags of atmos_params must have the same number of elements. Returning...'
    return,-1
  ENDIF
  ;
  IF n_elements(atmos_params.h_elec) EQ 1 THEN swtch=1b
ENDIF 

;
; If the He abundance is not specified, then use the default CHIANTI
; abundance file.
if n_elements(he_abund) eq 0 then begin
  read_abund,!abund_file,abund_arr,ref,element=2
  he_abund=abund_arr[1]
  message,/info,/cont,'He adundance has not been set - using default file '+!abund_file
endif

IF keyword_set(ct) AND swtch EQ 0 THEN BEGIN
  ;
  IF n_tags(atmos_params) EQ 0 THEN BEGIN
    ;
    ; If atmosphere_file does not exist, the user will be prompted to select
    ; a file from the CHIANTI options.
    atmos_params=ch_read_atmos(atmosphere_file)
    IF n_tags(atmos_params) EQ 0 THEN return,-1.
  ENDIF 

  IF NOT tag_exist(atmos_params,'temp') THEN BEGIN
    message,/info,/cont,'The tag TEMP (temperature) does not exist for atmos_params. Returning...'
    return,-1
  ENDIF ELSE BEGIN 
    atmos_temp=atmos_params.temp
  ENDELSE 

  ;
  ; Note the use of index_good below. It prevents h2_frac (and he3_frac) from having
  ; non-zero values outside of the atmos parameter temperature range.
  ;
  h_elec=ch_interp_atmos(atmos_params.h_elec,atmos_temp,temp)
  h1_frac=ch_interp_atmos(atmos_params.h1_frac,atmos_temp,temp,index_good=index_good)
  h2_frac=dblarr(n_elements(temp))
  h2_frac[index_good]=1.0d0-h1_frac[index_good]

  IF tag_exist(atmos_params,'he1_frac') AND tag_exist(atmos_params,'he2_frac') THEN BEGIN
    he1_frac=ch_interp_atmos(atmos_params.he1_frac,atmos_temp,temp)
    he2_frac=ch_interp_atmos(atmos_params.he2_frac,atmos_temp,temp,index_good=index_good)
    he3_frac=dblarr(n_elements(temp))
    he3_frac[index_good]=1.0d0-he1_frac[index_good]-he2_frac[index_good]
  ENDIF ELSE BEGIN
    print,'He ion fractions have not been specified - using default file '+!ioneq_file
    read_ioneq,!ioneq_file,ioneqt,he_fracs,ref,element=2
    ;
    ; Limit ioneq temperature to the range of the atmosphere model.
    k=where(10.^ioneqt GE min(atmos_temp) AND 10.^ioneqt LE max(atmos_temp))
    ;
    he1_frac=ch_interp_atmos(he_fracs[k,1,0],10.^ioneqt[k],temp,/spline)
    he2_frac=ch_interp_atmos(he_fracs[k,1,1],10.^ioneqt[k],temp,index_good=index_good,/spline)
    he3_frac=dblarr(n_elements(temp))
    he3_frac[index_good]=1.0d0-he1_frac[index_good]-he2_frac[index_good]
  ENDELSE 

  ct_model={h_elec:h_elec,$
            h1_frac:h1_frac>0.,$
            h2_frac:h2_frac>0.,$
            he1_frac:he1_frac>0.,$
            he2_frac:he2_frac>0.,$
            he3_frac:he3_frac>0.,$
            he_abund:he_abund}
    
  IF tag_exist(atmos_params,'elec_dens') THEN BEGIN
    total_elec=ch_interp_atmos(atmos_params.elec_dens,atmos_temp,temp)
    ct_model=add_tag(ct_model,total_elec,'total_elec')
  ENDIF
  
  IF tag_exist(atmos_params,'height') THEN BEGIN
    model_height=ch_interp_atmos(atmos_params.height,atmos_temp,temp)
    ct_model=add_tag(ct_model,model_height,'model_height')
  ENDIF
  
  IF tag_exist(atmos_params,'pressure') THEN BEGIN
    total_pressure=ch_interp_atmos(atmos_params.pressure,atmos_temp,temp)
    ct_model=add_tag(ct_model,total_pressure,'total_pressure')
  ENDIF 

ENDIF 

;
; This is for the special case where the H and He data in atmos_params are
; scalars.
;
IF keyword_set(ct) AND swtch EQ 1 THEN BEGIN
  ;
  ; If the helium fractions are not specified in atmos_params, then set the
  ; He fractions to zero.
  ;
  IF tag_exist(atmos_params,'he1_frac') AND tag_exist(atmos_params,'he2_frac')  THEN BEGIN
    he1_frac=double(atmos_params.he1_frac[0])
    he2_frac=double(atmos_params.he2_frac[0])
    he3_frac=(1d0-he1_frac-he2_frac)>0.
  ENDIF ELSE BEGIN
    he1_frac=0d0
    he2_frac=0d0
    he3_frac=0d0
    he_abund=0d0
  ENDELSE 
  ;
  ct_model={h_elec: double(atmos_params.h_elec[0]), $
            h1_frac: double(atmos_params.h1_frac[0]), $
            h2_frac: (1d0-atmos_params.h1_frac[0])>0., $
            he1_frac:he1_frac,$
            he2_frac:he2_frac,$
            he3_frac:he3_frac,$
            he_abund:he_abund}
ENDIF



if keyword_set(ct) then $
  str={ions_nlevels:ions_nlevels,model_ions:model_ions,ct_model:ct_model} $
else str={ions_nlevels:ions_nlevels,model_ions:model_ions}

return,str


end
