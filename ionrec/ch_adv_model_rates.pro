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
;       CH_ADV_MODEL_RATES 
;
;
; PURPOSE:
;
;	Calculate initial-level resolved or overall ionisation and recombination rates
;       out of an ion. For a given input ion, the ionisation rates are for ionisation out
;       of that charge state into the next higher charge state, and the recombination rates
;       are for recombination from the input ion into the next lower charge state.
;
;
; EXPLANATION:
;
;       Used for density-dependent advanced models where overall ionisation and recombination
;       rates depend on population of metastable levels, as opposed to the coronal
;       approximation where they only depend on the ground level rates.
;
;       Reads level-resolved direct ionisation rates in files ending .dilvl and indirect
;       ionisation (excitation--auto-ionisation) rates in files ending .ealvl, if available.
;       Uses the routine CH_IONIZ_RATE_LR to read the files and interpolate the rates over
;       the temperature grid given in the input.
;
;       If files are unavailable, uses CHIANTI rates for ground levels and approximates metastable
;       level rates by using Burgess and Chidichimo (1983) approximation. Their approximation
;       is improved by using it for ionisation rates from metastable and ground levels, taking
;       the ratio of these two rates and multiplying the CHIANTI ground rates by the ratio.
;       This is an equivalent way of estimating the constant C in the Burgess and Chidichimo 
;       formula. The current method could be improved because currently it uses the number of
;       equivalent electrons for the ground level, but some metastable levels have a different
;       number of equivalent electrons. A more accurate method is found in the thesis by
;       Dickson (1993, https://www.adas.ac.uk/theses/dickson_thesis.pdf).
;
;       Recombination rate coefficients are calculated from the fitting coefficients given
;       by N.R. Badnell (http://apap-network.org/DATA_AS/). These fit within a high degree of
;       accuracy the ab initio, total recombination data calculated by Badnell (2006) for 
;       radiative recombination (RR) and beginning with Badnell et al (2003) for dielectronic
;       recombination (DR). The total rates are resolved by initial level. DR rates are 
;       suppressed using the fitting formula of Nikolic et al (2013,2018), which mimic the
;       suppression of total DR rates into Rydberg levels as density increases. Recombination
;       into these levels can be followed by ionisation at higher densities before radiative decay
;       into lower levels. Nikolic et al fit the suppression shown in the Summers (1974) tables,
;       but in some cases the fitting formula does not match the results from Summers.
;
;       Charge transfer (CT) rates are calculated if files of rate coefficients are available.
;       They are found in files ending .ctilvl and .ctrlvl, corresponding to ionisation and
;       recombination data, respectively. The rates are collated from a variety of sources,
;       as detailed in the CHIANTI v.11 paper. They are mostly from quantum mechanical
;       calculations. Rate coefficients are read by CH_IONIZ_RATE_LR and interpolated over the
;       temperature grid. Rates are calculated by using the H, He ion fractions, H to electron
;       ratio and He abundance provided in the MODEL_ATM input. These are taken from model
;       atmosphere files. CT rates will not be calculated at temperature points above and
;       below the highest and lowest temperatures in the model atmosphere file. 
;
;       Overall rates are calculated by finding the level populations of an ion at each
;       point in the temperature and density grid and multiplying the initial-level resolved
;       rates by the fractional population of each ground and metastable level. These are 
;       passed back as the output. To speed up the calculation, level populations are 
;       calculated using a limited number of levels for ions where there is a large model
;       in CHIANTI. The number of levels are found next to the ion in the advmodel_list.ions file.
;
;
; CATEGORY:
;
;       CHIANTI; ionisation; recombination; overall rates
;
;
; CALLING SEQUENCE:
;
;       temp=10.^(findgen(61)*0.05+3.5)
;       dens=fltarr(61)+1.e11
;       metastable_levels,'c_1',metas,quiet=quiet
;       meta=where(metas eq 1)+1
;       recs=fltarr(24,n_elements(meta))
;
;       rates=ch_adv_model_rates('c_1',meta,temp,dens,recs)
;       rates=ch_adv_model_rates('c_1',meta,temp,dens,recs,/level_resolved)
;
;
; INPUTS:
;
;       THIS_ION:  The string in CHIANTI format for the ion.
;
;       META_INDEX:  Indexes of the metastable levels. These must match the level indexing
;              given in the ion energy file.
;
;       MODEL_TEMP:  1D array specifying the temperatures (in K) for which the rates
;              are needed.
; 
;       MODEL_DENSITY:  1D array specifying the electron number density (units: cm^-3)
;              for which the density-dependent rates are calculated. Array length
;              must match the temperature array length.
;
;       REC_COEFFS:  The coefficients required to calculate the recombination rates.
;
;	
; KEYWORDS:
;
;       VERBOSE:  Output various messages
;
;       LEVEL_RESOLVED: Instead of calculating overall ionisation and recombination rates,
;              output the rates resolved by initial level, i.e. for each ground and
;              metastable level in the ion
;
;
; OPTIONAL INPUTS:
;
;       MODEL_ATM:  Parameters from the model atmosphere for CT rates. These can be
;              drawn from CH_ADV_MODEL_SETUP. They contain the H, He ion fractions, H to electron
;              ratio, temperature and density grids and He abundance.
;
;       N_LEVELS:  The maximum number of levels to be included when solving the level populations.
;
;
; OUTPUTS:
;
;       FINAL_IONIZ: a 1D array of the overall ionisation rates for the ion at the specified
;               temperature points
;
;       FINAL_RECOMB: a 1D array of the overall recombination rates for the ion at the specified
;               temperature points
;
;
; OPTIONAL OUTPUTS:
;
;       NONE
;
;
; CALLS:
;
;       ch_ioniz_rate_lr, ch_burgchid_rate, ch_nikolic_dr_suppression
;        get_populations;
;
;
; PREVIOUS HISTORY:
;
;       NONE
;
;
; WRITTEN:
;         
;       Roger Dufresne (RPD) and Giulio Del Zanna (GDZ)
;       DAMTP, University of Cambridge, 16 Sept 2023 
;
;
; MODIFIED:
;
;       v.2, 17 Oct 2023 GDZ, added option to calculate populations
;            with a limited number of levels.
;
;       v.3, 10 Apr 2024  RPD, altered call to ch_ioniz_rate_lr because of changed
;              keywords. Created keyword for output of rates resolved by initial level
;              instead of calculating overall ionisation and recombination rates.
;
;       v.4, 18-May-2024 GDZ
;            moved the calc_recrates_fits function here. 
;
;
; VERSION:  4
;
;- 

function calc_recrates_fits,temp,rec_fits


ntemp=n_elements(temp)
nmeta=n_elements(rec_fits[0,*])
rr_rate=dblarr(ntemp,nmeta)
dr_rate=dblarr(ntemp,nmeta)


for im=0,nmeta-1 do begin

  rec_coeffs=double(rec_fits[*,im])

  rrb=rec_coeffs[1]+rec_coeffs[4]*exp(-rec_coeffs[5]/temp)
  rr_rate[*,im]=rec_coeffs[0]/(sqrt(temp/rec_coeffs[2])*$
    (1+sqrt(temp/rec_coeffs[2]))^(1.0d0-rrb)*$
    (1+sqrt(temp/rec_coeffs[3]))^(1.0d0+rrb))


  for it=0,8 do dr_rate[*,im]=dr_rate[*,im]+$
    temp^(-3.0d0/2.0d0)*rec_coeffs[6+it]*exp(-rec_coeffs[15+it]/temp)

endfor


rec_rates={rr_rate:rr_rate,dr_rate:dr_rate}

return,rec_rates


end


function ch_adv_model_rates,this_ion,meta_index,model_temp,model_density,$
            rec_coeffs,model_atm=model_atm,verbose=verbose,quiet=quiet,$
            n_levels=n_levels,level_resolved=level_resolved
            

ryd_ev=13.605698d0
equiv_elecs=[1,2,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,9,10,11,12]

convertname,this_ion,el,effchg
ion2filename,this_ion,ion_file
;print, this_ion,el,effchg
ntemp=n_elements(model_temp)
nmeta=n_elements(meta_index)

ci_rates=dblarr(ntemp,nmeta)
rec_rates=dblarr(ntemp,nmeta)
cti_rates=dblarr(ntemp,nmeta)
ctr_rates=dblarr(ntemp,nmeta)


; get ionisation rates

if effchg le el then begin

  di_file=ion_file+'.dilvl'

  if file_exist(di_file) then begin
    di_data=ch_ioniz_rate_lr(di_file,temp=model_temp,verbose=verbose)

    for im=0,nmeta-1 do begin
      dilvl=dblarr(ntemp)
      diind=where(di_data.lev_i eq meta_index[im],nr)
      for ir=0,nr-1 do $
        dilvl=dilvl+di_data.rates[*,diind[ir]]
      ci_rates[*,im]=ci_rates[*,im]+dilvl
    endfor
    
    ea_file=ion_file+'.ealvl'

    if file_exist(ea_file) then begin
      ea_data=ch_ioniz_rate_lr(ea_file,temp=model_temp,verbose=verbose)

      for im=0,nmeta-1 do begin
        ealvl=dblarr(ntemp)
        eaind=where(ea_data.lev_i eq meta_index[im],nr)
        for ir=0,nr-1 do $
          ealvl=ealvl+ea_data.rates[*,eaind[ir]]
        ci_rates[*,im]=ci_rates[*,im]+ealvl
      endfor
    endif
    
  endif else begin

    ; if ab initio rates are not available use CHIANTI rates for ground
    ; and Burgess and Chidichimo to estimate metastable rates
    ch_grd=ioniz_rate(this_ion,model_temp,verbose=verbose) ; GDZ
    ci_rates[*,0]=ch_grd

    if nmeta gt 1 then begin  
      grd_ip=ch_ip(this_ion)
      grd_burgchid=ch_burgchid_rate(effchg-1,equiv_elecs[el-effchg],grd_ip,model_temp)
      grdinf=where(grd_burgchid eq 0.0,nginf)

      read_elvlc,ion_file+'.elvlc',l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,ref

      for im=1,nmeta-1 do begin
        if eryd[meta_index[im]-1] gt 0.0 then meng=eryd[meta_index[im]-1] $
          else meng=erydth[meta_index[im]-1]
        meta_ip=grd_ip-meng*ryd_ev
        meta_burgchid=ch_burgchid_rate(effchg-1,equiv_elecs[el-effchg],meta_ip,model_temp)

        ionis_factor=meta_burgchid/grd_burgchid
        if nginf gt 0 then ionis_factor[grdinf]=0.0d0
        ci_rates[*,im]=ch_grd*ionis_factor
      endfor
    endif
      
  endelse

endif


; retrieve recombination data

if effchg gt 1 then begin

   rec_data=calc_recrates_fits(model_temp,rec_coeffs)
   sfactor=ch_nikolic_dr_suppression(this_ion,model_temp,density=model_density,quiet=quiet)
   
  for im=0,nmeta-1 do rec_rates[*,im]=rec_data.rr_rate[*,im]+rec_data.dr_rate[*,im]*sfactor

endif


; find charge transfer (CT) data if available

if n_elements(model_atm) gt 0 then begin

  if effchg le el then begin
    cti_file=ion_file+'.ctilvl'

    if file_exist(cti_file) then begin
      cti_data=ch_ioniz_rate_lr(cti_file,temp=model_temp,/ct)

      for im=0,nmeta-1 do begin
        ctilvl=dblarr(ntemp)
        ctiind=where(cti_data.lev_i eq meta_index[im],nr)
        for ir=0,nr-1 do begin
          if cti_data.perturber_z[ctiind[ir]] eq 1 then $
            ctilvl=ctilvl+cti_data.rates[*,ctiind[ir]]*model_atm.h2_frac $
          else if cti_data.perturber_z[ctiind[ir]] eq 2 then begin
            if cti_data.perturber_elecs[ctiind[ir]] eq 1 then $
              ctilvl=ctilvl+cti_data.rates[*,ctiind[ir]]*model_atm.he2_frac*model_atm.he_abund
            if cti_data.perturber_elecs[ctiind[ir]] eq 0 then $
              ctilvl=ctilvl+cti_data.rates[*,ctiind[ir]]*model_atm.he3_frac*model_atm.he_abund
          endif
        endfor

        cti_rates[*,im]=cti_rates[*,im]+ctilvl*model_atm.h_elec

      endfor
        
    endif
  endif
  
  if effchg gt 1 then begin
    ctr_file=ion_file+'.ctrlvl'

    if file_exist(ctr_file) then begin
      ctr_data=ch_ioniz_rate_lr(ctr_file,temp=model_temp,/ct)

      for im=0,nmeta-1 do begin
        ctrlvl=dblarr(ntemp)
        ctrind=where(ctr_data.lev_i eq meta_index[im],nr)
        for ir=0,nr-1 do begin
          if ctr_data.perturber_z[ctrind[ir]] eq 1 then $
            ctrlvl=ctrlvl+ctr_data.rates[*,ctrind[ir]]*model_atm.h1_frac $
          else if ctr_data.perturber_z[ctrind[ir]] eq 2 then begin
            if ctr_data.perturber_elecs[ctrind[ir]] eq 2 then $
              ctrlvl=ctrlvl+ctr_data.rates[*,ctrind[ir]]*model_atm.he1_frac*model_atm.he_abund
            if ctr_data.perturber_elecs[ctrind[ir]] eq 1 then $
              ctrlvl=ctrlvl+ctr_data.rates[*,ctrind[ir]]*model_atm.he2_frac*model_atm.he_abund
          endif
        endfor

        ctr_rates[*,im]=ctr_rates[*,im]+ctrlvl*model_atm.h_elec

      endfor
        
    endif
  endif
  
endif


; begin level population calculation and then determine overall rates,
; cut off for number of levels in ion will be used if given in advanced model ion list

if n_elements(level_resolved) gt 0 then begin

  final_ioniz=ci_rates+cti_rates
  final_recomb=rec_rates+ctr_rates

endif else begin
  if nmeta gt 1 then begin

    level_pops=dblarr(ntemp,nmeta)
    final_ioniz=dblarr(ntemp)
    final_recomb=dblarr(ntemp)

    get_populations,this_ion,model_temp,n_levels,pops=these_pops,densities=model_density,$ 
      /noionrec,/no_rrec,/no_auto,/pressure,verbose=verbose

    for im=0,nmeta-1 do begin
      level_pops[*,im]=reform(these_pops[meta_index[im]-1,*])

      final_ioniz=final_ioniz+level_pops[*,im]*(ci_rates[*,im]+cti_rates[*,im])
      final_recomb=final_recomb+level_pops[*,im]*(rec_rates[*,im]+ctr_rates[*,im])
    endfor

  endif else begin

    final_ioniz=ci_rates+cti_rates
    final_recomb=rec_rates+ctr_rates

  endelse
endelse


; density is not included in overall rates because some routines,
; such as make_ioneq_all, use rate coefficients rather than rates
final_rates={final_ioniz:final_ioniz,final_recomb:final_recomb}

return,final_rates


end
