PRO pop_solver, input, t, xne, pop,  data_str=data_str, $
                sum_mwl_coeffs=sum_mwl_coeffs, radfunc=radfunc, $
                frac_cutoff=frac_cutoff, pressure=pressure, $
                verbose=verbose, n_levels=n_levels, out_rates=out_rates, error=error,$
                noionrec=noionrec,no_rrec=no_rrec, all_levels=all_levels, $
                no_auto=no_auto

;+
; PROJECT:  CHIANTI
;
;       CHIANTI is an atomic database package for the calculation of
;       astrophysical emission line spectra.  
;
; NAME:  POP_SOLVER
;       
; PURPOSE:
;
;	To solve the level balance equations for Chianti ions.
;
; CATEGORY:
;       CHIANTI; level populations.
; 
; EXPLANATION:
;
;	This routine solves the level balance equations for the CHIANTI ions. 
;       The rates are passed as in input, and thge calculation is
;       carried out at each specific Te, Ne.
;
;       The matrix equation Ax=b is solved where A contains all the atomic 
;       data (electron rate coefficients, radiative decay rates, proton rate 
;       coefficients, photoexcitation rates), x are the level populations, 
;       and b a vector set to zeros except for the first element which is 1.
;
;       To solve the matrix equation, POP_SOLVER calls out to the CHIANTI
;       routine MATRIX_SOLVER.
;
;       From v.5 of CHIANTI the additional atomic processes of ionization
;       and recombination can be included when calculating the level
;       populations. These processes are not included in the matrix A.
;       Instead the level populations x are 'corrected' for ionization and
;       recombination afterwards. This correction is performed by the routine
;       CORRECT_POPS.
;      
;       NEW for v.9: 
;       The routines takes as input a structure with the data for the
;       ion, then runs  CH_LOAD_ION_RATES to load the rates of the
;       ion:
;
;       RATES   a structure created by CH_LOAD_ION_RATES with the  following tags:
;           n_levels  No. of levels in the model ion.
;           aa    A-values (2D array)
;           aax   Photoexcitation/stimulated emission (2D array).
;           qq    Electron rate coefficients (3D array:
;                 n_temperatures,n_levels,n_levels) 
;           ppr   Proton rate coefficients (3D array).
;
;           temp  Array of temperatures (the input TEMP).
;           ion_data
;                 The structure with the atomic data read by
;                 CH_SETUP_ION.
;
;       In the case of ions without autoionising levels,
;       the routine works in the same way as before, the only
;       difference is now the input. If an ion has autoionising
;       levels, then CH_LOAD_ION_RATES is called to load the rates of
;       the next higher ion, and CH_LOAD_2ION_RATES is used to combine
;       the rates of the two ions and between the ions into single
;       matrices. MATRIX_SOLVER is then called as usual. Note that the
;       total population is normalised to 1.
;
; CALLING SEQUENCE:
;
;	POP_SOLVER, RATES, T, XNE, POP
;
; INPUTS:
;
;       INPUT    The structure returned by ch_setup_ion. The 
;                structure has the tags:
;
;       .gname   Name of the ion in CHIANTI format.
;       .jj      Array containing J-values for all levels.
;       .ecm     Array containing level energies. Energies are
;                observed if available otherwise they are
;                theoretical. 
;       .ecmth   Array containing theoretical energies.
;       .wvl     2D array containing wavelengths. No negative values.
;       .a_value 2D array containing A-values.
;       .wgfastr Structure containing output from read_wgfa_str.
;       .splstr  Structure containing output from read_scups.
;       .ioneq_file The name of the ion balance file.
;       .ip      Ionization potential in cm^-1 units.
;
;       The following tags will be included if the relevant data-sets
;       exist:
;       .prot_struc  Structure containing proton rate data.
;       .ionrec  Structure containing level-resolved ionization and
;                recombination data.
;       .abund_file  The name of an element abundance file.
;       .dilute  The radiation dilution factor (derived from RPHOT).
;       .radtemp The value of RADTEMP.
;       .autostr Structure containing autoionization rates from the
;                .auto files (same format as returned by read_auto). 
;
;
;	T	Temperatures [K], e.g., 10.^6 (can be array)
;
;	XNE	Densities [cm-3], e.g., 10.^8 (can be array)
;
; OPTIONAL INPUTS:
;
;	N_LEVELS	This allows the number of levels in the model to 
;			be reduced. E.g., if the full model contains 100 
;			levels, one could set n_levels=50. This can be 
;			useful if one is interested in looking at the 
;			effects of cascading from higher levels
;
;       SUM_MWL_COEFFS  An array of coefficients of the same length as 
;                       the array of temperatures. Electron and proton rate 
;                       coefficients will be calculated at each temperature 
;                       and then a weighted sum of the coefficients is 
;                       performed using SUM_MWL_COEFFS. This allows 
;                       non-Maxwellian energy distributions to be 
;                       incorporated into the level balance equations.
;
;       RADFUNC         The name of a user-defined function that will generate
;                       a radiation spectrum as a function of temperature. 
;                       This radiation field will replace the black-body that
;                       is assumed when using the RADTEMP keyword in the call
;                       to pop_solver.
;
;       FRAC_CUTOFF     The fraction of non-zero elements in the C matrix below
;                       which the sparse matrix solver is used. See the routine
;                       matrix_solver for more details.
;
; KEYWORD PARAMETERS:
;
;       VERBOSE     If set, then a number of informational messages
;                   will be printed. (Intended to help with
;                   bug_checking.) 
;
;       NOIONREC: If set, then level-resolved ionization/recombination
;                rates are not read for the ion, even if they exist.
;
;       NO_RREC: If set, do not read the level-resolved radiative
;                recombination (RR) files.
;
;       ALL_LEVELS: If set, then the population array includes all
;                levels included in the calculation. Specifically, for
;                those ions for which recombination for the next
;                higher ionization stage is included, then /all_levels
;                means that POP will include the levels of this
;                additional ion. This is intended for use with the
;                DATA_STR optional output.
;
;       NO_AUTO:  If set, then the autoionization rates are not
;                 used. This means the 2-ion model will not be
;                 created.
;
;       PRESSURE: If set, then it is assumed that the populations are
;                 calculated for constant pressure. The number of
;                 temperatures must match the number of densities, and
;                 pressure=T*N_e. POP is returned as a 2D array of n_P
;                 x n_levels, rather than the usual 3D array.
;
; OUTPUTS:
;	POP	An array of level populations of size 
;		n_T x n_XNE x n_levels, unless the /PRESSURE keyword
;		is given. In this case a 2D array of size n_P x
;		n_levels is returned. If SUM_MWL_COEFFS has been
;		specified, then the output array is 1 x n_XNE x
;		n_levels. 
;
; OPTIONAL OUTPUTS:
;
;       DATA_STR If POP_SOLVER is called for just 1 temperature and density, 
;                then the individual data arrays for each of the physical 
;                processes can be output through DATA_STR. This allows the 
;                user to check for the dominant processes affecting the 
;                population of a given level. DATA_STR is a structure with 
;                the following tags:
;
;                .aa          A-values (2D array)
;                .aax         Photoexcitation/stimulated emission (2D array)
;                .cc          Electron rate (2D array)
;                .ccp         Proton rate (2D array)
;                .ai          Autoionization rate (2D array)
;                .dc          Dielectronic capture rate (2D array)
;                .rr          Radiative recombination rate (2D array)
;                .ion_rate    Ionization rate (1D array)
;                .rec_rate    Recombination rate (1D array)
;                .correction  Correction factor for level pop (1D array)
;                .frac_low    Ratio of N+1 ion fraction to N (scalar)
;                .frac_high   Ratio of N-1 ion fraction to N (scalar)
;
;                The 2D arrays are such that, e.g., aa[0,20] 
;                corresponds to an excitation, while aa[20,0] is a 
;                de-excitation. For the collisional quantities (cc,
;                ccp, etc.) the rates are given, not the rate
;                coefficients (i.e., they have been multiplied by
;                density). 
;
;                The 1D arrays are simply the rate coefficients into the
;                individual levels.
;
;       OUT_RATES: This is the structure returned by either
;                ch_load_ion_rates or ch_load_2ion_rates (depending if
;                the ion has autoionization rates), and contains the
;                rate and rate coefficient matrices for the different
;                processes. The format differs from DATA_STR.
;
; EXAMPLES:
;       IDL> dens=[1e9,1e10]
;       IDL> temp=[3e5,1e6]
;       IDL> input=ch_setup_ion('o_6')
;       IDL> pop_solver,input,temp,dens,pop
;
;       IDL> pop_solver,input,temp,dens,pop,/pressure
;
; CALLS:
;       ch_load_ion_rates,ch_load_2ion_rates,matrix_solver,correct_pops
;
; MODIFICATION HISTORY:
;      v.1  24-Aug-2018 Giulio Del Zanna (GDZ) 
;            New routine for version 9, based on the previous v.8
;            pop_solver,written by Peter Young.
;
;      v.2   GDZ  6-Oct-2018
;            fixed definition of status_ionrec; fixed default verbose=0
;            fixed if there are no data files returns an error and -1
;            Also, in the case of constant density calculations, speed
;            up the routine by just calculating once for each temperature.
;      v.3   GDZ 5 Dec 2018 
;            Significant rewrite, following more closely the original
;            pop_solver. Added verbose (3b)
;      v.4   14-Dec-2018 GDZ, added NOIONREC, NO_RREC keywords.
;      v.5,  5-Mar-2019, Peter Young
;            Added /all_levels keyword; added data_str= optional
;            output; added /no_auto keyword.
;      v.6,  6-Mar-2019, Peter Young
;            Fixed bug for /pressure case; updated header.
;
; VERSION     : v.6
;
;-

  if n_elements(verbose) eq 0 then verbose=0

xne = DOUBLE(xne)
t = DOUBLE(t)
; need the following to turn t into an array if it only has 1 element

IF n_elements(t) EQ 1 THEN BEGIN
  t0=t
  t=dblarr(1)
  t[0]=t0
ENDIF

nt=N_ELEMENTS(t)       ; no. of temperatures
nd=n_elements(xne)     ; no. of densities


rates1= ch_load_ion_rates(input, T, n_lev=n_levels, $
                            radfunc=radfunc, abund_file=abund_file, ioneq_file=ioneq_file, $ 
                            PATH=PATH, NOPROT=NOPROT, RPHOT=RPHOT, RADTEMP=RADTEMP,$
                            sum_mwl_coeffs=sum_mwl_coeffs, no_auto=no_auto, verbose=verbose,$
                          wvlmin=wmin,wvlmax=wmax, index_wgfa=anylines,obs_only=obs_only,$
                          noionrec=noionrec,no_rrec=no_rrec   )

convertname, rates1.ion_data.gname,iz,ion



IF keyword_set(pressure) THEN BEGIN
  IF nt NE nd THEN BEGIN
    print,'%POP_SOLVER: if /PRESSURE is set, then T and XNE must have the same size. Returning...'
error=1
    return
  ENDIF 
ENDIF



; if the selected ion has autoionizing states, then we add the next
; higher ion, combine the matrices and solve.

IF tag_exist(rates1.ion_data,'autostr') AND NOT keyword_set(no_auto) THEN BEGIN

  nlev1=rates1.n_levels

  zion2name,iz,ion+1,gname2

; GDZ: 
; We need to load the rates of the higher ion, but removing the 
; autoionizing levels, to avoid the complication of adding the next
; ion. 
; We also do not allow to restrict the number of levels in this case. 

  rates2=ch_load_ion_rates(gname2,t,sum_mwl_coeffs=sum_mwl_coeffs, /no_auto,$
                           radfunc=radfunc, abund_file=abund_file, ioneq_file=ioneq_file,$
                           NOPROT=NOPROT, RPHOT=RPHOT, RADTEMP=RADTEMP, verbose=verbose )

; if there are no data rates2=-1

  if is_number(rates2) then begin 
    print, '% POP_SOLVER: no data available for ion '+gname2
    error=1
    return
  endif 

  nlev2=rates2.n_levels
     
  rates=ch_load_2ion_rates(rates1,rates2, error=error, verbose=verbose, no_rrec=no_rrec)

  if error then return
  
ENDIF else rates=rates1


;
; Set up the POP output array, and handle SUM_MWL_COEFFS.
;
IF keyword_set(all_levels) THEN nlev=rates.n_levels ELSE nlev=rates1.n_levels
;
IF n_elements(sum_mwl_coeffs) EQ 0 THEN BEGIN
  sumtst=0
  IF keyword_set(pressure) THEN BEGIN
    pop=dblarr(nt,nlev)
  ENDIF ELSE BEGIN 
    pop=dblarr(nt,nd,nlev)
  ENDELSE 
ENDIF ELSE BEGIN
  sumtst=1
  pop=dblarr(1,nd,nlev)
  IF nt NE n_elements(sum_mwl_coeffs) THEN BEGIN
    print,'%POP_SOLVER: number of temperatures must match the size of '+$
         'SUM_MWL_COEFFS.'
    print,'             Populations not calculated.'
error=1
    return
  ENDIF
ENDELSE


  

IF tag_exist(rates1.ion_data,'ionrec')  then begin 
  if n_elements(sum_mwl_coeffs) NE 0 THEN  status_ionrec=0 else $
     status_ionrec=1
endif else status_ionrec=0 

IF sumtst EQ 0 THEN BEGIN
 ;
  IF NOT keyword_set(pressure) THEN BEGIN 

    FOR j=0,nt-1 DO BEGIN
      
      FOR i=0,nd-1 DO BEGIN

        if keyword_set(verbose) then $
           print,'calculating populations  for log T,N: ', alog10(t[j]),alog10(xne[i])

        p= matrix_solver(xne[i],rates,index_t=j, c=c, frac_cutoff=frac_cutoff,verbose=verbose)

        p=p/total(p)                                        ; renormalise

   ; Includes rec/ion corrections if ion is not He-like
   ; Note that this code works with the combined 2-ion models.
        IF status_ionrec GT 0 and iz-ion ne 1 THEN begin 
          p=correct_pops(p, t[j], xne[i],  rates1.ion_data.ionrec, C , $
                         rrate=rrate,crate=crate,correction=correction, $
                         frac_low=frac_low, frac_high=frac_high)
        ENDIF ; corrections
       ;
       ; Reduce  back to single ion model if /all_levels not set.
       ;
        IF NOT keyword_set(all_levels) THEN BEGIN
          p=p[0:nlev-1]
          p=p/total(p)
        ENDIF 
        pop[j,i,*]=p

      ENDFOR
    ENDFOR

  ENDIF   ELSE BEGIN             ; we have two arrays: T,N at constant pressure.

    FOR j=0,nt-1 DO BEGIN

      if keyword_set(verbose) then $               
         print,'calculating populations at constant pressure for log T,N: ',$
               alog10(t[j]), alog10(xne[j])

      p= matrix_solver(xne[j],rates,index_t=j, c=c, frac_cutoff=frac_cutoff,verbose=verbose)

;          IF NOT keyword_set(all_levels) THEN  p=p[0:nlev-1] ; remove the populaton of the higher ion:
      p=p/total(p)  ; renormalise

; Includes rec/ion corrections if ion is not He-like    
      IF status_ionrec GT 0 and iz-ion ne 1 THEN begin 
;              IF tag_exist(rates1.ion_data,'autostr')  THEN BEGIN
; the C matrix should only have the levels of the main ion, 
; so we remove the levels of the higher ion:

;                 C= C[0:nlev1-1,0:nlev1-1]
;              ENDIF  ; note: correct_pops renormalises to 1.
        p=correct_pops(p, t[j], xne[j],  rates1.ion_data.ionrec, C , $
                       rrate=rrate,crate=crate,correction=correction, $
                       frac_low=frac_low, frac_high=frac_high)
      ENDIF ; corrections 
       ;
       ; Reduce  back to single ion model if /all_levels not set.
       ;
      IF NOT keyword_set(all_levels) THEN BEGIN
        p=p[0:nlev-1]
        p=p/total(p)
      ENDIF 
      pop[j,*]=p
    ENDFOR 
  ENDELSE   

ENDIF    ELSE BEGIN   
  FOR ii=0,nd-1 DO BEGIN
    p=matrix_solver(xne[ii], rates,  /sumtst,c=c,frac_cutoff=frac_cutoff, $
                    verbose=verbose)
    IF NOT keyword_set(all_levels) THEN BEGIN
      p=p[0:nlev-1]
      p=p/total(p)
    ENDIF 
    pop[0,ii,*]=p
  ENDFOR
END 

bad=where(pop-pop NE 0,nbad)
IF nbad GT 0 THEN BEGIN
  print,'% POP_SOLVER: Warning - NaN values for ion '+gname
  pop[*]=0.
ENDIF




;
; DATA_STR is principally for sending arrays to the routine POP_PROCESSES,
; and since this routine will only call POP_SOLVER with 1 temperature and
; 1 density, then only bother filling DATA_STR if this is the case.
;
IF nt EQ 1 AND nd EQ 1 THEN BEGIN
  IF n_elements(rrate) EQ 0 THEN rrate=0.
  IF n_elements(crate) EQ 0 THEN crate=0.
  IF tag_exist(rates,'ai') THEN ai=rates.ai ELSE ai=0.
  IF tag_exist(rates,'rr') THEN rr=xne[0]*rates.rr ELSE rr=0.
  IF tag_exist(rates,'dc') THEN dc=xne[0]*rates.dc ELSE dc=0.
  IF n_elements(correction) EQ 0 THEN correction=0.
  IF n_elements(frac_low) EQ 0 THEN frac_low=0.
  IF n_elements(frac_high) EQ 0 THEN frac_high=0.
  data_str={nlev_ion: rates1.n_levels, $
            nlev_all: rates.n_levels, $
            aa: rates.aa, $
            aax: rates.aax, $
            cc: xne[0]*reform(rates.qq[0,*,*]),  $
            ccp: xne[0]*reform(rates.ppr[0,*,*]), $
            ai: ai, $
            rr: reform(rr), $
            dc: reform(dc), $
            rec_rate: rrate, $
            ion_rate: crate, $
            correction: correction, $
            frac_low: frac_low, $
            frac_high: frac_high}
ENDIF


out_rates=temporary(rates)


END 
