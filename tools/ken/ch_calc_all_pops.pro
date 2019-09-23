
function ch_calc_all_pops, ionname, temp=temp, dens=dens, ltemp=ltemp, ldens=ldens, $
                      frac_cutoff=frac_cutoff, verbose=verbose

;+
;     This routine calculates level populations for a 2D array of
;     temperatures and densities.
;-

COMMON elvlc,l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref
COMMON wgfa, wvl,gf,a_value
COMMON upsilon, splstr
COMMON radiative, radt,dilute
COMMON proton, pstr, pe_ratio
COMMON elements,abund,abund_ref,ioneq,temp_all,ioneq_ref
COMMON ionrec,rec_rate,ci_rate,temp_ionrec,luprec,lupci,status

t0=systime(1)

IF n_elements(dens) NE 0. THEN ldens=alog10(dens)
IF n_elements(temp) NE 0. THEN ltemp=alog10(temp)

IF n_elements(ldens) EQ 0 THEN ldens=findgen(12)+1.
IF n_elements(ltemp) EQ 0 THEN ltemp=findgen(26)/5.+4.0


IF keyword_set(diel) THEN diel = 1 ELSE diel = 0

;
; The following extracts the names of the files to be read. I want to allow 
; a different path to be chosen, and I extract only the information I need 
; from filename.
;

convertname,ionname,iz,ion
zion2filename,iz,ion,filename,diel=diel,name=name
gname = trim(ionname)

IF N_ELEMENTS(path) NE 0 THEN filename=concat_dir(path, name)

chckfile=filename+'.splups'
chck=file_exist(chckfile)
IF chck NE 1 THEN BEGIN
  IF NOT keyword_set(quiet) THEN print,'% SHOW_POPS: no data exists for this ion in CHIANTI. Returning...'
  popstr=0
  return,-1
ENDIF


setup_ion,name,-1,-1,wvltst,lvl1,lvl2,wvl1,gf1,a_value1,path=path, $
     noprot=noprot

pe_ratio=proton_dens(ltemp)

IF n_elements(rphot) EQ 0 THEN dilute=0d0 ELSE dilute=r2w(rphot)
IF n_elements(radtemp) EQ 0 THEN radt=6d3 ELSE radt=double(radtemp)

status=-1

IF status lt 0 THEN ionrec = {rec:0., ci:0., temp:0., lev_up_rec:0., lev_up_ci:0., status:-1, ioneq:0., temp_ioneq:0.}

;----------X
;GDZ

input = {gname:name, jj:jj, ecm:ecm,ecmth:ecmth, $
 wvl:wvl, a_value:a_value, splstr:splstr, $
 pe_ratio:pe_ratio,prot_struc:pstr, dilute:dilute, radtemp:radt, ionrec:ionrec}


pop_solver,input, 10.^ltemp,10.^ldens,pop,n_levels=n_levels, $
     sum_mwl_coeffs=smc, radfunc=radfunc, frac_cutoff=frac_cutoff, verbose=verbose

t1=systime(1)
print,'% CH_CALC_ALL_POPS: time taken (s): ',t1-t0

return,pop

END
