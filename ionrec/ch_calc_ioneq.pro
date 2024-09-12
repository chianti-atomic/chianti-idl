;+
;
; PROJECT:  CHIANTI
;
;       This program was developed as part of CHIANTI-VIP. 
;       CHIANTI-VIP (Virtual IDL and Python) is  a member of the CHIANTI family
;       mantained by Giulio Del Zanna, to develop additional features and
;       provide them to the astrophysics community.
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the
;       University of Cambridge, Goddard Space Flight Center, and University of Michigan.
;
;
; NAME:
;
;       CH_CALC_IONEQ
;
;
; PURPOSE:
;
;       Compute equilibrium ion fractions in CHIANTI format. 
;
;
; EXPLANATION:
; 
;       Adaption of make_ioneq_all routine, which calculates ion charge states at zero density,
;       reading ionization and recombination data from files present in each CHIANTI ion directory.
;
;       With the keyword /ADVANCED_MODEL, which is switched on by default, density-dependent
;       ionization and recombination from metastable levels and DR suppression are included.
;       Charge transfer (CT) effects can be included by using the keyword /CT.
;       For details of methods see
;
;       Dufresne, R.P. and Del Zanna, G., 2019, A&A, 626, A123
;       Modelling Ion Populations in Astrophysical Plasmas: Carbon in the Solar Transition Region
;
;       Dufresne, R.P. and Del Zanna, G., 2020, MNRAS, 497, 1443
;       Effects of density on the oxygen ionisation equilibrium in collisional plasmas,
;       
;       Dufresne R. P., Del Zanna G., Badnell N. R., 2021b, MNRAS, 503, 1976
;       The influence of photo-induced processes and charge transfer on carbon and oxygen
;       in the lower solar atmosphere
;
;       Dufresne, R. P., Del Zanna, G., and Storey, P. J., 2021, MNRAS, 505, 3968
;       Modelling low charge ions in the solar atmosphere
;
;       Dufresne, R. P., Del Zanna, G., et al., 2024, ApJ (10.48550/arXiv.2403.16922)
;       -----------------------------------------------------------------------------------
;
;       In the present models, densities corresponding to the temperatures need to be
;       defined. This can be done in three ways, similar to e.g. CH_SYNTHETIC:
;       1) input density, either as a constant scalar or an array corresponding to 
;            length of temperature array;
;       2) a constant pressure;
;       3) a model input file with the list of temperatures and densities.
;
;       The ion charge states are returned and can be output to a file.
;       The format of this array is the standard CHIANTI one:
;       ioneq_data=dblarr(ntemp,30,31) where ntemp is the number of temperatures. 
;
;       Optionally, only one or a list of elements can be calculated, via the keyword
;       ELEMENTS. In this case, the full array is still returned, but populated with zeros,
;       except the elements that have been calculated. This is to speed up calculation of
;       ion balance in on-the-fly routines, such as CH_SYNTHETIC.
;       
;
;       To include CT effects, the abundances of neutral H, He need to be known.
;       They can be input via a simple formatted file using ATMOSPHERE_FILE or via a.
;       widget. For the convenience of the user a few files have been prepared, from the
;       works of Avrett E.H., Loeser R., 2008, ApJ, 175, 229 and
;       Fontenla J. M., Landi E., Snow M., Woods T., 2014, Sol. Phys., 289, 515,
;       the latter including a variety of solar regions. The files are located in 
;       the ancillary_data/advanced_models directory of the CHIANTI database.
;
;       The routine first runs CH_ADV_MODEL_SETUP to read the list of ions
;       (advmodel_list.ions) for which advanced models are provided, and load into
;       memory the coefficients for the calculation of the rates.
;
;       After that, CH_ADV_MODEL_RATES is used to load up the rates.
;       Ionization rates for all levels are calculated from ab initio data
;       in files if available or using CHIANTI defaults for the ground level and using
;       Burgess, A., & Chidichimo, M. C. 1983, MNRAS, 203, 1269 approximation to estimate
;       metastable level rates. The independent atom model is used to solve the level
;       populations and then overall rates are calculated from these to give density-
;       dependent ionization and recombination. 
;
;       NOTE: to speed up the routine, level populations are calculated using a maximum
;       number of levels for some ions with large models. The higher states that are
;       neglected do not contribute significantly to the populations of the metastable states.
;
;       The density dependent effects are also affected by suppression of DR, which is
;       approximated using the Nikolic et al (2018) fits to the results of the hydrogenic models
;       developed by Burgess and Summers (1969) and with ion fractions published in Summers (1974).
;
;       Charge transfer rates are calculated if requested. The rate coefficients are extracted
;       from files using various sources. See the CHIANTI v.11 paper for the list of data 
;       sources. CT can only be calculated in the advanced models and only for ions for which
;       there are files; no approximations for rates are used in the absence of these files.
;
;       Instead of solving the full population arrays for an element, we use an approximation
;       where overall total rates from/to an ion are used to obtain the ion fraction
;       by taking the ratio of the ionization/recombination rates at each temperature,
;       as in the zero density case. The results can optionally be outputted to a CHIANTI .ioneq
;       file using the keyword OUTNAME.
;
;
; CATEGORY:
;
;       CHIANTI; ionization; recombination, ion fractions.
;
;
; CALLING SEQUENCE:
;
;	data=ch_calc_ioneq() 
;
;
; INPUTS:
;
;       TEMPERATURES: 1D array specifying the temperatures (in K) for which the ion fraction
;              curves are needed. If not specified, the range logT=4.0
;              to logT=8.0 in 0.05 dex intervals is used, unless MODEL_FILE is used.
;
;       DENSITY:  Specifies the electron number density (units: cm^-3)
;                 for which the density-dependent ion fractions are calculated.
;
;
; KEYWORDS:
;
;       ADVANCED_MODEL: include density-dependent ionization and recombination.
;              This can be switched off by using ADVANCED_MODEL=0.
;              Switching this off is only used for on-the-fly routines, such as CH_SYNTHETIC and
;              CHIANTI_DEM, otherwise the user should use default CHIANTI ion balance file.
;
;       CT: include charge transfer in advanced models
;
;       DR_SUPPRESSION: Switch on DR suppression from Nikolic et al (2018) for all ions 
;              not included in the advanced models. The comparison with Summers (1974) suppression
;              has not been checked for other elements when preparing the models.
;
;
; OPTIONAL INPUTS:
; 
;       OUTNAME:  The name of the output file. If not specified, then
;                 the ion fractions are NOT written to file. 
; 
;       PRESSURE: Specifies the pressure (N_e*T; units cm^-3 K) at
;                 which the density-dependent ion fractions are calculated.
;
;       MODEL_FILE: The name of the two-column model file with the list
;                   of temperatures (K) and densities (units: cm^-3). The model atmosphere
;                   files provided for CT data can also be used as a MODEL_FILE input;
;                   all the additional data in those files are ignored in this input.
;
;       ELEMENTS: The list of elements to be calculated. This is useful
;                 in combination with IONEQ_DATA, to calculate the ion fractions
;                 on-the-fly within other CHIANTI routines.
;
;       ATMOSPHERE_FILE: A file with the H, He ion fractions as a function of temperature.
;                 For the quiet Sun the file avrett_atmos.dat could be used, with data from
;                 Avrett and Loeser (2008), otherwise multiple regions are provided from
;                 Fontenla et al (2014). The files are located in the ancillary_data directory.
;                 If not given as an input and CT is switched on, a widget will request a file 
;                 to be selected. If the user prefers a self-consistent set of temperatures and
;                 densities for the ion balance as is used for the CT data, then the same
;                 atmosphere file here should also be given as an input in MODEL_FILE described above.
;
;       HE_ABUND:  The total helium abundance relative to hydrogen. The default CHIANTI 
;                  abundance for helium will be selected if this is not given and CT
;                  is switched on.
;
;       ERR_MSG: a variable to receive error messages
;
;       WARNING_MSG: a string array to receive relevant warning messages
;
;
; OUTPUTS:
;
;       IONEQ_DATA=dblarr(ntemp, 30,31) where ntemp is the number of
;       temperatures, containing ion fractions of all ions of all elements
;       up to and including zinc (if calculated, otherwise zero).
;
;       OUTNAME:  The file containing ion fractions. If not defined, a file is
;                 written in the working directory. 
;
; OPTIONAL OUTPUTS:
;
;
;       ERR_MSG: a string with an error message
;
;       WARNING_MSG: an array of strings with warning messages
;
;
; EXAMPLES:
;
;       iontemp=10.^(findgen(61)*0.05+3.5)
;       
;       data=ch_calc_ioneq(iontemp, dens=1.e10, /adv, ele='C')
;       
;       data=ch_calc_ioneq(iontemp, press=3e14, /adv, ele=['C','N','O'])
;
;       data=ch_calc_ioneq(iontemp, dens=1.e10, /adv, /ct, out='all_ions_ne=1e10.ioneq')
;
;       data=ch_calc_ioneq(iontemp, advanced_model=0, outname='all_ions_zero_density.ioneq')
;
;       ch_compare_ioneq, 'C', files=['all_ions_ne=1e10.ioneq','all_ions_zero_density.ioneq'],$
;       lab=['Ne=1e10','Ne=0'],/top,/right,psym=[6,5], lines=[0,2], ion=[1,4]
;
;
; CALLS:
; 
;      CH_ADV_MODEL_SETUP
;      CONVERTNAME, ZION2NAME, RDFLOAT
;      METASTABLE_LEVELS
;      CH_ADV_MODEL_RATES
;      IONIZ_RATE, RECOMB_RATE  for the zero-density case
;      
;
; PREVIOUS HISTORY:
;
;    The last part of this function is the same as the older MAKE_IONEQ_ALL
;    written by Ken Dere and modified by Peter Young.
;
;
; MODIFICATION HISTORY:
;
;    v.1, 17 Oct 2023  Roger Dufresne (RPD) and Giulio Del Zanna (GDZ)
;               DAMTP, University of Cambridge
; 
;         By default the routine works as before when called by make_ioneq_all,
;         calculating the abundances of the ion charge states at zero density.
;
;    v2, 24 Oct 2023, GDZ, removed a leftover keyword, ioneq_data
;                     Also, an error was converted to a warning, when the 
;                requested element is not available in the current advanced model;
;                in this way, the element is calculated with the zero-density approx.
;                This was necessary to run ch_synthetic. Fixed the return.
;
;    v.3, 16 Nov 2023, GDZ, modified dealing with atmosphere file.
;
;    v.4, 13 Dec 2023, RPD, allowed for comments at end of model file, to enable
;                use of the supplied atmosphere files for the ioneq temperature grid.
;                Also corrected retrieving metastable indices from recombination data.
;
;    v.5, 06 Feb 2024, RPD, added DR suppression for all coronal model ions unless
;                 keyword no_suppression selected. Altered reading of model file.
;
;    v.6, 27 Feb 2024, RPD, changed keyword for DR suppression for coronal ions
;                so that DR suppression is off unless requested. This is because 
;                of effect on Fe from Nikolic et al, which does not seem to match
;                Summers. Wrote header and finalised comments in the file.
;
;    v.7, 27 Mar 2024, RPD, Modified model file section of routine to fix certain issues.
;                Modified comments in output of ioneq file.
;                Allow a density array to be inputted as well as a single constant density.
;
;    v.8, 10 May 2024, RPD, partial preparation of density variable being output when
;                using a model file. Made temperature and density arrays double
;                precision in model file section for consistency with rest of routine.
;                Allowed elements to be entered by atomic number. Switched on
;                advanced models by default.
;
;    v.9, 17 May 2024, GDZ - added a reference to the paper and modified so that all elements
;                are now calculated. Only those requested (and available) are calculated with
;                the advanced model. All the others are calculated as before at zero density.
;                Also, by default a file is written in the working directory.
;                Completed density variable being output when using a model file.
;
;     v.10, 24 May 2024, GDZ - Completed routine to write output file.
;     v.11, 29 May 2024, RPD - Sundry changes to formatting, moved index_req check to element stage
;     v.12, 3-July-2024, GDZ - changed the definition of the output ionization equilibrium file.
;
;     v.13, 05 Sep 2025, RPD - altered to now use .rrcoeffs and .drcoeffs files for each ion
;                rather than loading the recombination coefficients in memory and passing to
;                the advanced model rate routine. Added comment to ioneq file if CT is included.
;                Included passing quiet keyword to other routines used. Altered reading system time
;                when creating ioneq name.
;
; VERSION    v.13
;-


function ch_calc_ioneq,temperatures,outname=outname,elements=elements,density=density,$
                       pressure=pressure,model_file=model_file,advanced_model=advanced_model,ct=ct,$
                       atmosphere_file=atmosphere_file,he_abund=he_abund,verbose=verbose,quiet=quiet,$
                       err_msg=err_msg,warning_msg=warning_msg,dr_suppression=dr_suppression

  t1=systime(1)
  err_msg=''
  warning_msg=''
  
  if n_params() lt 1 then begin
    err_msg='%CH_CALC_IONEQ: ERROR, please enter either a temperature range or a variable name to receive the temperature array'
    print,err_msg & return,-1
  endif
  
  if n_elements(verbose) eq 0 then verbose=0 
  if verbose then quiet=0 else quiet=1

  ; By default write an ion equilibrium file.
  if n_elements(outname) gt 0 then begin

     if file_exist(outname) then begin
         pp=strsplit(anytim(!stime,/vms),/extract)
       outname='ch_adv_'+trim(pp[0])+'-'+strmid(pp[1],0,8)+'.ioneq'
    
      print,'% CH_CALC_IONEQ: requested ionization file exists. Writing output to '+outname
    end
     
   endif else begin
 
      pp= strsplit(anytim(!stime, /vms),/extract)
       outname='ch_adv_'+pp[0]+'-'+strmid(pp[1],0,8)+'.ioneq'
      ;; outname='adv_'+$
      ;;   trim(string(SYSTIME(/JULIAN, /UTC), format='(f14.4)'))+'.ioneq'

     print,'% CH_CALC_IONEQ: Writing output to '+outname
 
  end
  
     
  ; By default for now do not include charge transfer, but make the advanced models the default
  if n_elements(advanced_model) eq 0 then advanced_model=1
  if n_elements(ct) eq 0 then ct=0 ; GDZ

  ioneqmin=1.e-20
  zlabl=['H','He','Li','Be','B','C','N','O','F','Ne','Na',$
        'Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr', $
        'Mn','Fe','Co','Ni','Cu','Zn']

     
  ; begin various tests on inputs

  if n_elements(elements) gt 0 then begin
    if isa(elements[0],/int) eq 1 then begin
        requested_iz=elements
        count=n_elements(elements)
    endif else requested_iz=1+where_arr(zlabl, elements, count)
    
    if count ne n_elements(elements) then begin 
        err_msg='% CH_CALC_IONEQ: ERROR in the input definition of the elements'
        print,err_msg & return,-1
    end
    
    if not keyword_set(advanced_model) then begin
        err_msg='% CH_CALC_IONEQ: ERROR, to calculate single elements you need to add density-dependent effects with the keyword /advanced_model'
        print,err_msg & return,-1
    end
    
  ; GDZ:  else request all the elements:      
  endif else requested_iz=1+indgen(30)

  
  if n_elements(model_file) gt 0 then begin ; GDZ

    if n_elements(density) gt 0 then begin 
      err_msg='% CH_CALC_IONEQ: ERROR, density cannot be defined if model is defined, it is an output !'
      print,err_msg & return,-1
    end
    
    if n_elements(temperatures) gt 0 then begin 
      err_msg='% CH_CALC_IONEQ: ERROR, temperatures cannot be defined if model is defined!'
      print,err_msg & return,-1
    end
    
    if not file_exist(model_file) then begin 
      err_msg='% CH_CALC_IONEQ: ERROR, input file '+model_file+' not found!'
      print,err_msg & return,-1
    end
    
    ; obtain the temperatures and densities from the model file
    mlines=file_lines(model_file)
    model_data=strarr(mlines)
    openr,lf,model_file,/get_lun
    readf,lf,model_data
    free_lun,lf
    
    comments=where(model_data.startswith('-1'),nc)
    if nc gt 0 then nrecs=comments[0] else nrecs=mlines

    rdfloat,model_file,temperatures,densities,numline=nrecs
    temperatures=double(temperatures) & densities=double(densities)
    
    ; all atmospheric parameters below the temperature minimum and above the temperature maximum
    ; have to be removed for temperature interpolation in ionization equilibrium

    atmos_min=min(temperatures,atmin)
    atmos_max=max(temperatures,atmax)
    if atmin ne 0 or atmax ne nrecs then begin
      print,'Double-valued temperature array found in atmosphere file,'
      print,'parameters below temperature minimum and above temperature maximum will be removed'
      temperatures=temperatures[atmin:atmax]
      densities=densities[atmin:atmax]
      log_temperatures=alog10(temperatures)
      ntemp=n_elements(temperatures)
    endif
    
    ne_pe_comment='Model used Te, Ne from file: '+model_file
    for ii=0,n_elements(temperatures)-1 do $
      ne_pe_comment=[ne_pe_comment,string(temperatures[ii])+' '+string(densities[ii])]

    ; now return density array
    density=densities
    
  endif else begin 

    ; the arrays of temperatures and density have to be defined if a model file is not used
    ; if temperature not set, then use log T= 4-8 

    if n_elements(temperatures) eq 0 then begin
      log_temperatures=findgen(81)/20.0+4.0
      temperatures=10.0d^log_temperatures
      print,'Temperatures array has not been specified: default CHIANTI values will be used'
    endif else temperatures=double(temperatures) 
    ntemp=n_elements(temperatures)

    if n_elements(density) gt 0 and n_elements(pressure) gt 0 then begin 
      err_msg='% CH_CALC_IONEQ: ERROR, either density or pressure must be specified'
      print,err_msg & return,-1
    end
    
    if keyword_set(advanced_model) then begin
      if n_elements(density) eq 0 and n_elements(pressure) eq 0 then begin 
        err_msg='% CH_CALC_IONEQ: ERROR, either density or pressure must be specified for the advanced model'
        print,err_msg & return,-1
      end
    endif

    if n_elements(density) eq 1 then begin 
      densities=dblarr(ntemp)+double(density)
      ne_pe_comment='Model used constant density='+string(density)
    endif else if n_elements(density) gt 1 then begin
      if n_elements(density) ne ntemp then begin 
        err_msg='% CH_CALC_IONEQ: ERROR, density array provided is not consistent with temperature array'
        print,err_msg & return,-1
      endif else densities=double(density)

    endif else if n_elements(density) eq 0 then begin 
      if n_elements(pressure) gt 1 then begin 
        err_msg='% CH_CALC_IONEQ: ERROR, pressure array provided but only scalar allowed'
        print,err_msg & return,-1
      endif else if n_elements(pressure) eq 1 then begin
        densities=double(pressure[0])/temperatures
        ne_pe_comment='Model used constant pressure='+string(pressure)
      endif else $
        if keyword_set(verbose) then print,'% CH_CALC_IONEQ: calculating the zero-density case...'
    endif 
    
  endelse    ; end of tests on inputs


  ; now retrieve the global parameters used throughout the advanced model
  if keyword_set(advanced_model) then begin

    ion_comments=''

    if keyword_set(ct) then begin
      if n_elements(atmosphere_file) eq 0 then begin
        atmosphere_file=ch_get_file(path=concat_dir(concat_dir(concat_dir(!xuvtop,'ancillary_data'),$
          'advanced_models'),'model_atmospheres'),filter='*.dat',tit='Select an atmosphere file for the charge transfer calculation') 
      endif else begin
        if not file_exist(atmosphere_file) then begin
          err_msg='% CH_CALC_IONEQ: ERROR, input atmosphere file does not exist'
          print,err_msg & return,-1
        endif
      endelse
      ion_comments=[ion_comments,'Charge transfer has been included in these ion balances']
      ion_comments=[ion_comments,' using the model atmosphere file - '+atmosphere_file]
    endif   

    params=ch_adv_model_setup(temperatures,ct=ct,atmosphere_file=atmosphere_file,he_abund=he_abund)

    model_ions=params.model_ions
    ions_nlevels=params.ions_nlevels ; GDZ

    ; calculate by default all elements
    calculate_iz=indgen(30)+1
    
    ; determine which elements are to be included in the advanced model calculation
    
    ; check which elements are present in the masterlist         
    available_iz_ions=intarr(n_elements(model_ions))
    for iion=0,n_elements(model_ions)-1 do begin 
      convertname,model_ions[iion],iz,ion
      available_iz_ions[iion]=iz
    end
    
    available_iz=available_iz_ions[rem_dup(available_iz_ions)]

    ; let user know whether all elements are part of advanced model
    for ii=0,n_elements(requested_iz)-1 do begin
      index=where(available_iz eq requested_iz[ii],nn)
      if nn eq 0 then begin
        warning_msg=[warning_msg,'% CH_CALC_IONEQ WARNING: requested element '+$
          zlabl[requested_iz[ii]-1]+' not available in current advanced model']
        if verbose then print,warning_msg 
      endif
    endfor 

    ; the final list of elements included in the ion balance calculation
    ; for which advanced models will be calculated
    requested_iz=requested_iz[sort(requested_iz)] 
    
  endif else begin           ; end of setup for advanced model
    
    ; setup for when advanced models are not requested
    if keyword_set(ct) then begin
      err_msg='% CH_CALC_IONEQ: ERROR, charge transfer can only be included in the advanced models'
      print,err_msg & return,-1
    endif
    
    calculate_iz=indgen(30)+1 ; calculate the zero density of all elements
    
  endelse

     
  ; initialise the standard CHIANTI array - main difference is that
  ; this is an array of double-precision numbers
  ioneq_data=dblarr(ntemp, 30,31)
  
  ; begin retrieving ionisation and recombination rates
  for iiz=0,n_elements(calculate_iz)-1 do begin

    z=calculate_iz[iiz]
    ion_rate=fltarr(ntemp,z+2)
    rec_rate=fltarr(ntemp,z+2)
    ; RPD: moved so that this is checked once for the element
    index_req=where(requested_iz eq z, nreq)

    for ion=0,z do begin

        zion2name,z,ion+1,gname
        if verbose then print,gname
        
        if keyword_set(advanced_model) then begin

          ; if ion is in advanced model list and we requested it,
          ; get initial-level resolved rates:
          ionmatch=where(params.model_ions eq gname)
          ionmatch=ionmatch[0] ; GDZ needed !
          
          if ionmatch gt -1 and nreq eq 1 then begin

            ; define the maximum number of levels to calculate the populations
            if ions_nlevels[ionmatch] gt 0 then $
              n_levels=ions_nlevels[ionmatch] else n_levels=999 ; GDZ
            
            ; call the sub-routine to calculate overall ionisation and recombination rates
            if keyword_set(ct) then $
              rates_tr=ch_adv_model_rates(gname,temperatures,densities,$
                model_atm=params.ct_model,verbose=verbose,n_levels=n_levels,quiet=quiet) $
            else rates_tr=ch_adv_model_rates(gname,temperatures,densities,$
              verbose=verbose,n_levels=n_levels,quiet=quiet)

            ion_rate(*,ion)=rates_tr.final_ioniz[*]
            rec_rate(*,ion)=rates_tr.final_recomb[*]

            ion_comments=[ion_comments,gname+' included in density effects model']
            if keyword_set(verbose) then print,gname+' included in density effects model'
              
          endif else begin

            ; for ions not in advanced model list use ground level ionisation and recombination rates
            if ion+1 le z then ion_rate(*,ion)=ioniz_rate(gname,temperatures,verbose=verbose)
            if ion gt 0 then $
              if keyword_set(dr_suppression) then $
                rec_rate(*,ion)=recomb_rate(gname,temperatures,verbose=verbose,density=densities) $
              else rec_rate(*,ion)=recomb_rate(gname,temperatures,verbose=verbose)
              
          endelse

        ; end of calculating ionisation and recombination rates for advanced models
        endif else begin

          ; use the coronal approximation for everything if advanced models not requested
          if ion+1 le z then ion_rate(*,ion)=ioniz_rate(gname,temperatures,verbose=verbose)
          if ion gt 0 then $
            if keyword_set(dr_suppression) then $
              rec_rate(*,ion)=recomb_rate(gname,temperatures,verbose=verbose,density=densities) $
            else rec_rate(*,ion)=recomb_rate(gname,temperatures,verbose=verbose)
          
        endelse

    endfor      ; end of obtaining ionisation and recombination rates

    
    ; begin calculating ion fractions - original part of Ken Dere's routine
    ioneq=dblarr(ntemp,z+2)

    for it=0,ntemp-1 do begin
      factor=fltarr(z+1)

      for ion=1,z-1 do begin
        rat=ion_rate(it,ion)/rec_rate(it,ion)
        factor(ion)=rat^2+rat^(-2)
      endfor

      factor(0)=max(factor)
      factor(z)=max(factor)

      idx=where(factor eq min(factor))

      most=idx(0)
      ioneq(it,most)=1.d

      ; Get ions above most
      for ion=most+1,z+1 do begin
        if rec_rate(it,ion) gt 0. then begin
          ioneq(it,ion)=ion_rate(it,ion-1)*ioneq(it,ion-1)/rec_rate(it,ion)
          ;if ioneq(it,ion) lt ioneqmin then ioneq(it,ion)=0.
        endif else ioneq(it,ion)=0.
      endfor
                          ;
      for ion=most-1,0,-1 do begin
        ioneq(it,ion)=rec_rate(it,ion+1)*ioneq(it,ion+1)/ion_rate(it,ion)
        if ioneq(it,ion) lt ioneqmin then ioneq(it,ion)=0.
      endfor

      ioneq(it,0)=ioneq(it,*)/total(ioneq(it,*))

    endfor                  ; it

    ioneq_data[*, z-1, 0:z]= ioneq[*, 0:z]
    
  endfor                     ; loop over z

     
  ; formats for the output      
  tformat='('+trim(ntemp)+'f6.2)'
  iformat='('+trim(ntemp)+'e10.3)'
  
  ; if an output file is requested
  if n_elements(outname) gt 0 then begin 
    openw,luw,outname,/get_lun
    printf,luw,ntemp,30,format='(2i3)'
    printf,luw,alog10(temperatures),format=tformat

    for z=1,30 do begin
      for ion=0,z do begin
        printf,luw,z,ion+1,ioneq_data[*,z-1,ion],format='(2i3,'+iformat+')'
      endfor
    endfor
    
    ion_string=' '

    printf,luw,' -1'
    printf,luw,' ionization equilibrium filename:  ',file_basename(outname)

    if keyword_set(advanced_model) then begin
      printf,luw,' the following ions have advanced models:'
      for ii=0,n_elements(ion_comments)-1 do printf,luw,ion_comments[ii]
      for ii=0,n_elements(ne_pe_comment)-1 do printf,luw,ne_pe_comment[ii]
    endif
    
    printf,luw,' Produced as part of the CHIANTI atomic data base collaboration'
    datetime=systime()
    printf,luw,' Created on ',datetime
    printf,luw,' -1'
    free_lun,luw
  
    if verbose then begin
      print,'% CH_CALC_IONEQ: the ion fractions have been written to the file'
      print,'                  '+outname
    endif
    
  endif ; end of writing the file
  
  if verbose then begin
    k=where(temperatures LT 1e4,nk)
    IF nk GT 0 THEN BEGIN
      print,'**WARNING**'
      print,'Ion fractions have been computed below 10^4 K. These values may not be accurate as charge transfer'
      IF keyword_set(advanced_model) THEN BEGIN
        IF keyword_set(ct) THEN $
          print,'can be significant at low temperatures and this process has not been included for all ions.' $
        ELSE BEGIN
          print,'can be significant at low temperatures and this process has not been included in these ion balances.'
          print,'It can be included for some ions by using the keyword /ct.'
        ENDELSE
      ENDIF ELSE BEGIN
        print,'can be significant at low temperatures and this process has not been included in these ion balances.'
        print,'It can be included for some ions by using the keywords /advanced_model and /ct.'
      ENDELSE
    ENDIF 
  endif
  
  t2 = systime(1)
  print,format='("% CH_CALC_IONEQ: ionization eq.  computed in ",f8.1," seconds")',t2-t1

  return, ioneq_data  
     
  END
