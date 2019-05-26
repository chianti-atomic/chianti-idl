FUNCTION ch_load_2ion_rates, rates1, rates2, error=error, verbose=verbose, no_rrec=no_rrec 

;+
; NAME:
;      CH_LOAD_2ION_RATES
;
; PURPOSE:
;      Combines rates for two neighboring ions into single
;      matrices. If the ion1 model has autoionization rates 
;
; CATEGORY:
;      CHIANTI; rates.
;
; CALLING SEQUENCE:
;      Result = CH_LOAD_2ION_RATES(Rates1, Rates2)
;
; INPUTS:
;      Rates1: The structure returned by CH_LOAD_ION_RATES for the
;              reference ion  (charge=z).
;      Rates2: The structure returned by CH_LOAD_ION_RATES for the
;              ionized ion  (charge=z-1).
;
; KEYWORD PARAMETERS:
;      NO_RREC:  If set, then level-resolved radiative recombination
;                rates are switched off.
;      VERBOSE:  If set, then information messages will be printed to
;                the IDL input window.
;
; OUTPUTS:
;      A structure with the following tags:
;         n_levels:   No. of levels in the model.
;         aa:        2D array containing A-values.
;         aax:       2D array containing stimulated emission and
;                    photoexcitation rates.
;         ppr:       3D array containing proton rate coefficients.
;         qq:        3D array containing electron rate coefficients.
;         rr:        3D array containing radiative recombination
;                    rates.
;         ai:        2D array containing autoionization rates.
;         dc:        3D array containing dielectronic capture rates.
;         dr:        3D array containing dielectronic recombination
;                    rate coefficients.
;
;      Note: for the 3D arrays, the first dimension is temperature.
;
;      If the SUM_MWL_COEFFS input has been passed through RATES1,
;      then the collision arrays will be returned as 2D arrays rather
;      than 3D arrays.
;
;      If a problem is found then the integer -1 will be returned.
;
; OPTIONAL OUTPUTS:
;      Error:  Returned as 1 if an error occurs, 0 otherwise. 
;
; CALLS:
;      CONVERTNAME
;
; EXAMPLE:
;
;      IDL> t=10.^(findgen(11)/10.+5.0)
;      IDL> rates1=ch_load_ion_rates('o_6',t)
;      IDL> rates2=ch_load_ion_rates('o_7',t)
;      IDL> rates=ch_load_2ion_rates(rates1,rates2)
;
; MODIFICATION HISTORY:
;
;
;      ver  1, 19 Jul 2018, Giulio Del Zanna (GDZ)
;          merged earlier codes written by GDZ and Peter Young (PRY)
;
;      v.2, 4-Oct-2018, GDZ, added the verbose keyword.
;
;      Ver.3, 2-Feb-2019, Peter Young
;       Changed interpolation of level-resolved recombination rates to
;       be performed on logT-logRate instead of T-Rate.
;
;      Ver.4, 5-Mar-2019, Peter Young
;       Added /NO_RREC keyword; implemented sum_mwl_coeffs keyword
;       (which is passed through RATES1 input).
;
;      Ver.5, 15-Mar-2019, Peter Young
;       Modified expression for tot_rr, removing the ground
;       transition. 
;-


  convertname, rates1.ion_data.gname,iz1,ion1
  convertname, rates2.ion_data.gname,iz2,ion2
  error=0

  IF iz2 NE iz1 OR ion2 NE ion1+1 THEN BEGIN
     print,'% CH_LOAD_2ION_RATES: ERROR ! ion2 must be one higher ionization state compared to ion1.'
     error=1
     return,-1
  ENDIF 

  if n_elements(verbose) eq 0 then verbose=0

;
; Get temperature array.
;
  t=rates1.temp

; GDZ: check that the rates of the second ion have the same 
; temperature array: 

  if n_elements(rates1.temp) ne n_elements(rates2.temp) then begin 
     print, '%CH_LOAD_2ION_RATES:  ERROR  in the temperatures '
     error=1
     return, -1
  endif 

  if max(abs(rates2.temp-rates1.temp)) gt 1. then begin 
     print, '%CH_LOAD_2ION_RATES:  ERROR, temperature arrays of the two ions differ ! '
     error=1
     return, -1
  endif 

 ;
 ; If sumtst=1, then we'll be summing the collision rates over
 ; the Maxwellians. In this case, ntm is the no. of temperatures
 ; involved in the sum, while nt=1. Note that nt is used to define the
 ; sizes of the rate arrays. For sumtst=0, we have ntm=nt.
 ;
  sumtst=rates1.sumtst
  sum_mwl_coeffs=rates1.sum_mwl_coeffs
  ntm=n_elements(rates1.temp)

;
; Get no. of temperatures from the size of the QQ array.
;
  s=size(rates1.qq)
  IF s[0] EQ 2 THEN nt=1 ELSE nt=s[1]
  


  
;
; Create the new rate arrays for the combined ion model.
;
; GDZ: *** in what follows we assume that the first level in each ion is
; the ground state ****. In other words, the ground state of the lower
; ion has an IDL index=0, the ground state of the higher ion has an
; IDL index=n1 
;

  n1=rates1.n_levels
  n2=rates2.n_levels
  nlev_matrix=n1 + n2
;
  aa=dblarr(nlev_matrix,nlev_matrix)
  aax=dblarr(nlev_matrix,nlev_matrix)
  qq=DBLARR(nt,nlev_matrix,nlev_matrix)
  ppr=DBLARR(nt,nlev_matrix,nlev_matrix)

  aa[0:n1-1,0:n1-1]=rates1.aa
  aax[0:n1-1,0:n1-1]=rates1.aax
  qq[*,0:n1-1,0:n1-1]=rates1.qq
  ppr[*,0:n1-1,0:n1-1]=rates1.ppr

  aa[n1:*,n1:*]=rates2.aa
  aax[n1:*,n1:*]=rates2.aax
  qq[*,n1:*,n1:*]=rates2.qq
  ppr[*,n1:*,n1:*]=rates2.ppr

 ;
 ; The following removes the temperature index from the qq and ppr
 ; arrays for the case that sum_mwl_coeffs has been set.
 ;
  IF sumtst THEN BEGIN
    qq=reform(qq)
    ppr=reform(ppr)
  ENDIF 

  
; Get the total radiative and dielectronic recombination rates.
;
;tot_rr=recomb_rate(ion_data.gname,t,/rad)
;tot_dr=recomb_rate(ion_data.gname,t,/diel)

;
; Ionization (IONIZ array). This array contains the ionization rates
; from the lower to the higher ionisation stage.
; ----------------------------------------------

; For v.9 we only include the ground-to-ground total ionization rate:
;
;  ioniz[*,0,n1]= ioniz_rate(rates1.ion_data.gname,t)

  IF sumtst THEN BEGIN
    ioniz=dblarr(nlev_matrix,nlev_matrix)
    FOR i=0,ntm-1 DO BEGIN
      ioniz[0,n1]=ioniz[0,n1]+sum_mwl_coeffs[i]*ioniz_rate(rates1.ion_data.gname,t[i])
    ENDFOR 
  ENDIF ELSE BEGIN
    ioniz=dblarr(nt,nlev_matrix,nlev_matrix)
    ioniz[*,0,n1]= ioniz_rate(rates1.ion_data.gname,t)
  ENDELSE 

;
; Radiative recombination (RR array)
; -----------------------------------
  rr=dblarr(ntm,nlev_matrix,nlev_matrix)

;  total RR rate from ground state of recombining ion:
  rr_rate_ground= recomb_rate(rates2.ion_data.gname,t,/rad) ;


  IF tag_exist(rates1.ion_data,'rrec') AND NOT keyword_set(no_rrec) THEN begin 
; We have level-resolved RR rates to include in the model.

rrec=rates1.ion_data.rrec

;  rrec={rate:rate,temp:temp,$
;    final_level:final_level, initial_level:initial_level,ref:rref}

; check if there are values in the T array outside of the temperature
; array in the RR rates:

     ind1=where(T lt min(rrec.temp), nt1) 
     ind2=where(T ge min(rrec.temp) and T le max(rrec.temp), nt2)
     ind3=where(T gt max(rrec.temp), nt3)

     if nt1 gt 0 then print, $
        '% WARNING: RR rates requested for temperatures below those available for '+rates1.ion_data.gname

     if nt3 gt 0 then print, $
        '% WARNING: RR rates requested for temperatures above those available for '+rates1.ion_data.gname


; for loop ?
     for ii=0, n_elements(rrec.final_level)-1 do begin 
        
        if nt1 gt 0 then begin 
           ind_min= where(rrec.temp eq min(rrec.temp))
           ind_min=ind_min[0]

           rr[ind1, n1+rrec.initial_level-1, rrec.final_level[ii]-1]=$
              rrec.rate[ii, ind_min]
           
        endif 

       ;
       ; PRY, 22-Feb-2019
       ;   changed interpolation to be on the logarithm of the
       ;   temperatures and rates.
        if nt2 gt 0 THEN BEGIN
          log_rr=interpol(alog10(rrec.rate[ii,*]), alog10(rrec.temp), alog10(t[ind2]) >0.)
          rr[ind2, n1+rrec.initial_level[ii]-1, rrec.final_level[ii]-1]=$
             10.^log_rr
        ENDIF 


        if nt3 gt 0 then begin 
           ind_max= where(rrec.temp eq max(rrec.temp))
           ind_max=ind_max[0]

           rr[ind3, n1+rrec.initial_level-1, rrec.final_level[ii]-1]=$
              rrec.rate[ii, ind_max]
           
        endif 

     endfor 


; Now we need to total the rates from the ground state of the
; recombining ion and remove them from the totals to
; avoid double counting.
;
; PRY, 15-Mar-2019
;  I've changed the sum in the 3rd dimension from 0:* to 1:* as
;  the ground transition shouldn't be included.
;
     tot_rr= total(rr[*, n1, 1:*],3)

;; if verbose then begin 
;; print,'Total RR                   :', rr_rate_ground
;; print,'Total RR from level-resolved:', tot_rr
;; end 

     if min(rr_rate_ground- tot_rr) lt 0 then begin 

        print, 'ERROR ! the total of the level-resolved RR rates is greater than the total !?? '

        print, 'total level-resolved RR rates: ',  arr2str(string(tot_rr,format='(e8.2)'),/trim)
        print, 'total  RR rates: ',arr2str(string(rr_rate_ground,format='(e8.2)'),/trim)


        rr[*, n1, 0]= ( rr_rate_ground- tot_rr) > 0.

        wait,3

      endif else  rr[*, n1, 0]=  rr_rate_ground- tot_rr

    ;
    ; If sum_mwl_coeffs specified, then need to sum over temperatures
    ;
     IF sumtst THEN BEGIN
       rrx=dblarr(nlev_matrix,nlev_matrix)
       FOR i=0,ntm-1 DO rrx=rrx+sum_mwl_coeffs[i]*reform(rr[i,*,*])
       rr=temporary(rrx)
     ENDIF 

     if verbose then print, 'Included  the level-resolved RR rate'

  endif  else begin 

;  we only include the ground-to-ground rates.

    rr[*, n1, 0]=  rr_rate_ground
    IF sumtst THEN BEGIN
      rrx=dblarr(nlev_matrix,nlev_matrix)
      FOR i=0,ntm-1 DO rrx=rrx+sum_mwl_coeffs[i]*reform(rr[i,*,*])
      rr=temporary(rrx)
    ENDIF 

     if verbose then print, 'Included only the total RR rate from the ground state'

  end 

;
; Now check if we have autoionization data. If not, then we just
; return OUTPUT. Otherwise, we keep going.
;

  IF NOT tag_exist(rates1.ion_data,'autostr') THEN begin 

     output={ n_levels: nlev_matrix, $
              aa: aa, $
              qq: qq, $
              aax: aax, $
              ppr: ppr, $
              ioniz:ioniz,$
              rr:rr }

  endif else begin 

; ------------------------------
; Create the autoionization (AI) and dielectronic capture (DC) rate
; arrays. 
;
     ai=dblarr(nlev_matrix,nlev_matrix)
     dc=dblarr(ntm,nlev_matrix,nlev_matrix)

;
; Fill the AI array. Note that the lower levels are in ION2.
; -----------------
     lvl_l=rates1.ion_data.autostr.lvl1
     lvl_s=rates1.ion_data.autostr.lvl2
     a_auto=rates1.ion_data.autostr.auto
;
; Only include levels below n_levels. By default this is all levels,
; but the keyword n_lev can change this.
;
     k=where(lvl_s LE rates1.n_levels,nk)
     IF nk GT 0 THEN ai[lvl_s[k]-1,lvl_l[k]+n1-1]=a_auto[k]


;
; Fill the DC array.
; -----------------
;
; Define physical constants
;
     planck = 6.6260693d-27     ; #erg s
     ev2Erg = 1.602176487d-12
     boltzmann = 1.3806504d-16  ;  # cgs, i.e. kT is in ergs
     emass = 9.10938215d-28     ; #  electron mass in gram
     invCm2Ev = 1./8.06554465e+3
     const= planck^3. /2.d/ (2.d*!pi*emass)^(3./2.)

;
; Only include levels below n_levels. By default this is all levels,
; but the keyword n_lev can change this.
;
     k=where(lvl_s LE rates1.n_levels,nk)
     IF nk GT 0 THEN BEGIN 
                                ;
                                ; Get statistical weights
                                ;

; gs is the statistical weight of the autoionizing state

        gs=1.+2.*rates1.ion_data.jj[lvl_s[k]-1]

; gl is the statistical weight of the level of the recombining ion 
; that produces the dielctronic capture. 
        gl=1.+2.*rates2.ion_data.jj[lvl_l[k]-1]
                                ;
                                ; Compute energy and exponential term. Recall that ecm  is either observed  or theoretical
                                ;
        des=(rates1.ion_data.ecm[lvl_s[k]-1] - (rates1.ion_data.ip + rates2.ion_data.ecm[lvl_l[k]-1]))*invCm2Ev* ev2Erg
        des_kt=(1/(boltzmann*T)) # des
                                ;
        FOR i=0,ntm-1 DO BEGIN 
           dc_t=reform(dc[i,*,*])
           dc_t[lvl_l[k]-1+n1,lvl_s[k]-1]=planck^3. /2.d/(2.d*!pi*emass*boltzmann)^1.5 * $
                                          exp(-des_kt[i,*]) * $
                                          ( t[i]^(-1.5) # (gs/gl*a_auto) )
           dc[i,*,*]=temporary(dc_t)
        ENDFOR
     ENDIF 

    ;
    ; Handle the sum over Maxwellians (sum_mwl_coeffs).
    ;
     IF sumtst THEN BEGIN
       dcx=dblarr(nlev_matrix,nlev_matrix)
       FOR i=0,ntm-1 DO dcx=dcx+sum_mwl_coeffs[i]*reform(dc[i,*,*])
       dc=temporary(dcx)
     ENDIF 

; GDZ **** the following part could be rewritten for speed removing
; the for loops ****

; Subtract individual DR rates from the total DR rate
; ---------------------------------------------------
; We need to consider only the autoionising levels, and the
; dielectronic capture rate from the ground state, as for the totals
; we only have the total from the ground state. 

     lvl_auto= get_uniq(lvl_s, count=n_auto)
     total_dr=fltarr(nt)

     branching_ratio1=fltarr(n_elements(lvl_auto))
     
     for ii=0, n_elements(lvl_auto)-1 do begin 
        
        inda=where(lvl_s eq lvl_auto[ii], nna) 
; 
        if nna gt 1 then begin 
           if verbose then print, 'found '+trim(nna)+' lower levels for autoionizing level: ',lvl_s[ii]
           a_auto_tot= total(a_auto[inda]) 
        endif else a_auto_tot=a_auto[ii]

; Do first the approximation. total the radiative decay down to bound levels

        a_tot=total(rates1.ion_data.a_value[indgen(min(lvl_auto)-1),lvl_auto[ii]-1 ])
        
        branching_ratio1[ii]= a_tot/(a_auto_tot+a_tot)
        
; gs is the statistical weight of the autoionizing state
        gs=rates1.mult[lvl_auto[ii]-1]
        
; now check how many levels of the recombining ion are connected to
; the ground state ******* 
        ind_this_auto= where(lvl_s eq lvl_auto[ii] and lvl_l eq 1, n_decays)
        
        if n_decays gt 0 then begin 

; ecms is the energy  of the autoionising level
           ecms=rates1.ion_data.ecm[lvl_auto[ii]-1] ; this is either observed  or theoretical

           for id=0,  n_decays-1 do begin 
              
              gl=rates2.mult[lvl_l[ind_this_auto[id]] -1] 
              
; des is the energy of the autoionizing state-energy of the level of
; the recombining ion in ergs
              des= (ecms - (rates1.ion_data.ip+ rates2.ion_data.ecm[ lvl_l[ind_this_auto[id]] -1])) *invCm2Ev* ev2Erg
              
; need to find the index 
              ind_auto= where(lvl_l eq lvl_l[ind_this_auto[id]] and lvl_s eq lvl_auto[ii],nna)
              if nna ne 1 then begin 
                 print,'error !'
                 return, -1
              endif 
; rate coefficient for dielectronic capture. T can be an array: 
              rate= planck^3. /2.d/(2.d*!pi*emass*boltzmann*T)^(3./2.) * $
                    exp(-des/(boltzmann*T) ) * gs / gl * a_auto[ind_auto] 
              
              total_dr=total_dr+ (rate * branching_ratio1[ii])
           endfor 
        endif 
     endfor  

     if keyword_set(verbose) then print, 'total DR related to the autoionising levels: ', total_dr

     total_dr_badnell=recomb_rate(rates2.ion_data.gname,t,/diel) 

     if keyword_set(verbose) then print,'total DR: ',total_dr_badnell

; Now compute the DR array. Although this has the same size as QQ, DC,
; etc., only one transition actually has a non-zero rate, namely the
; "ground-to-ground" transition. In the future this may change.
;
; Note that I force the subtracted DR rate to be > 0. In some cases
; the level-resolved DR may be greater than the total DR rate 
;
     dr=dblarr(ntm,nlev_matrix,nlev_matrix)

     dr[*,n1,0]=  (total_dr_badnell-total_dr )>0.

     IF sumtst THEN BEGIN
       drx=dblarr(nlev_matrix,nlev_matrix)
       FOR i=0,ntm-1 DO drx=drx+sum_mwl_coeffs[i]*reform(dr[i,*,*])
       dr=temporary(drx)
     ENDIF 

     output={ n_levels: nlev_matrix, $
              aa: aa, $
              qq: qq, $
              aax: aax, $
              ppr: ppr, $
              ioniz:ioniz, rr:rr,$
              ai:ai, dc:dc, dr:dr }


  endelse 


  return,output


END
