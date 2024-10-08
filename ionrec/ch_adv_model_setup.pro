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
;       CH_ADV_MODEL_SETUP 
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
; VERSION:  8
;
;- 

function ch_adv_model_setup,temp,ct=ct,atmosphere_file=atmosphere_file,he_abund=he_abund



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


; retrieve model atmosphere data for charge transfer calculation
if keyword_set(ct) then begin

  atmos_data=strarr(file_lines(atmosphere_file))
  openr,lf,atmosphere_file,/get_lun
  readf,lf,atmos_data
  free_lun,lf
  
  comment_ind=where(atmos_data.startswith('-1'),nc)
  if nc gt 0 then nrecs=comment_ind[0] else nrecs=n_elements(atmos_data)
  
  rdfloat,atmosphere_file,atmos_temp,atmos_elec,atmos_height,$
    atmos_pressure,atmos_hyd,atmos_h1,atmos_he1,atmos_he2,numline=nrecs,/double

  ; all atmospheric parameters below the temperature minimum and above the temperature maximum
  ; have to be removed for temperature interpolation in ionization equilibrium
  atmos_min=min(atmos_temp,atmin)
  atmos_max=max(atmos_temp,atmax)
  if atmin ne 0 or atmax ne nrecs-1 then begin
    print,'Double-valued temperature array found in atmosphere file,'
    print,'parameters below temperature minimum and above temperature maximum will be removed'
    atmos_temp=atmos_temp[atmin:atmax]
    atmos_elec=atmos_elec[atmin:atmax]
    atmos_height=atmos_height[atmin:atmax]
    atmos_pressure=atmos_pressure[atmin:atmax]
    atmos_hyd=atmos_hyd[atmin:atmax]
    atmos_h1=atmos_h1[atmin:atmax]
    if n_elements(atmos_he2) gt 0 then begin
      atmos_he1=atmos_he1[atmin:atmax]
      atmos_he2=atmos_he2[atmin:atmax]
    endif
  endif
  
  log_temp=alog10(temp)
  ntemp=n_elements(temp)
  atmos_logt=alog10(atmos_temp)
    
  model_height=10.0d0^interpol(alog10(atmos_height),atmos_logt,log_temp)
  total_elec=10.0d0^interpol(alog10(atmos_elec),atmos_logt,log_temp)
  total_pressure=10.0d0^interpol(alog10(atmos_pressure),atmos_logt,log_temp)
  total_hyd=10.0d0^interpol(alog10(atmos_hyd),atmos_logt,log_temp)
  h1_frac=10.0d0^interpol(alog10(atmos_h1),atmos_logt,log_temp)

  if n_elements(atmos_he2) gt 0 then begin
    he1_frac=10.0d0^interpol(alog10(atmos_he1),atmos_logt,log_temp)
    he2_frac=10.0d0^interpol(alog10(atmos_he2),atmos_logt,log_temp)
    helim=where(log_temp lt atmos_logt[0] and log_temp gt atmos_logt[-1],natmhe)
  endif else begin
    read_ioneq,!ioneq_file,ioneqt,he_fracs,ref,element=2
    he1_frac=10.0d0^interpol(alog10(double(he_fracs[*,1,0])),ioneqt,log_temp)
    he2_frac=10.0d0^interpol(alog10(double(he_fracs[*,1,1])),ioneqt,log_temp)
    helim=where(log_temp lt ioneqt[0] or log_temp gt ioneqt[-1],natmhe)
    print,'He ion fractions have not been specified - using default file '+!ioneq_file
  endelse

  h_elec=total_hyd/total_elec
  h2_frac=1.0d0-h1_frac
  he3_frac=1.0d0-he1_frac-he2_frac

  ; set atmosphere parameters for CT to zero if requested temperature range lies outside
  ; range of model atmosphere file
  ctlim=where(log_temp lt atmos_logt[0] or log_temp gt atmos_logt[-1],natm)
  if natm gt 0 then begin
    print, 'Model atmosphere does not cover requested temperature range. Charge transfer rates will be set to zero outside model atmosphere limits'
    model_height[ctlim]=0.0d0 & total_elec[ctlim]=0.0d0 & total_pressure[ctlim]=0.0d0 & total_hyd[ctlim]=0.0d0
    h_elec[ctlim]=0.0d0 & h1_frac[ctlim]=0.0d0 & h2_frac[ctlim]=0.0d0
  endif

  if natmhe gt 0 then begin
    he1_frac[helim]=0.0d0 & he2_frac[helim]=0.0d0 & he3_frac[helim]=0.0d0
  endif
  he_err=where(he3_frac lt 0.0,nhe)
  if nhe gt 0 then he3_frac[he_err]=0.0d0
  ; test both h fracs and he fracs for values greater than 1

  if n_elements(he_abund) eq 0 then begin
    read_abund,!abund_file,abund_arr,ref,element=2
    he_abund=abund_arr[1]
    print,'He adundance has not been set - using default file '+!abund_file
  endif
  
  ct_model={h_elec:h_elec,h1_frac:h1_frac,h2_frac:h2_frac,he1_frac:he1_frac,$
    he2_frac:he2_frac,he3_frac:he3_frac,he_abund:he_abund,model_height:model_height,$
    total_elec:total_elec,total_pressure:total_pressure}

endif


if keyword_set(ct) then $
  str={ions_nlevels:ions_nlevels,model_ions:model_ions,ct_model:ct_model} $
else str={ions_nlevels:ions_nlevels,model_ions:model_ions}

return,str


end
