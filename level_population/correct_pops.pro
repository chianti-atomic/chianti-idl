
FUNCTION correct_pops, pp, t, xne, ionrec, cc, crate=crate, rrate=rrate, $
                       correction=correction, frac_low=frac_low, $
                       frac_high=frac_high

;+
; NAME
;
;    CORRECT_POPS()
;
; PROJECT
;
;    CHIANTI
;
; EXPLANATION
;
;    Corrects CHIANTI level populations with the ionization and recombination
;    rate coefficients
;
; INPUTS
;
;    PP      The level populations that need to be corrected.
;
;    T       Temperature at which calculation is performed. Units: K.
;
;    XNE     Electron density at which calculation is performed. Units: cm^-3
;
;    IONREC  Structure with the following tags
;            .rec        Effective recomb. rate coefficients
;            .ci         Effective ionization rate coefficients
;            .temp       Temperatures at which rates are tabulated
;            .lev_up_rec Levels to which recombination takes place
;            .lev_up_ci  Levels to which ionization takes place
;            .status     Either +1 (ion/rec data exists) or -1 (dosen't exist)
;            .ioneq      Ion fractions of the 3 ions
;
;    CC      2D matrix produced by MATRIX_SOLVER that contains the rate
;            coefficents from the standard CHIANTI processes.
;
; OPTIONAL OUTPUTS
;
;    CRATE   A 1D array of same size as POP containing the collisional
;            ionization rate coefficients (units: cm^3 s^-1).
;
;    RECRATE A 1D array of same size as POP containing the recombination
;            rate coefficients (units: cm^3 s^-1).
;
;    CORRECTION A 1D array of same size as POP containing the correction
;               factors for each level.
;
;    FRAC_LOW The ratio of the current ionization fraction to the fraction
;             of the one lower ion (i.e., less ionized).
;
;    FRAC_HIGH The ratio of the current ionization fraction to the fraction
;              of the one higher ion (i.e., more ionized).
;
; CALLS
;
;    ION_FRAC_INTERP(), CI_REC_INTERP()
;
; HISTORY
;
;    Ver.1, 10-Jun-2005, Peter Young
;        Taken original code of Enrico Landi and inserted it into a separate
;        routine.
;
;    Ver.2, 16-Aug-2005, Peter Young
;        Changed total_exc to be dblarr in order to prevent NaNs.
;
;    Ver.3, 1-Feb-2006, Peter Young
;        Corrected error: sum of populations is now renormalized to 1.
;
;    Ver.4, 12-Jun-2009, Enrico Landi
;        Changed the temperature array for ion fractions, now taken from
;        the IONREC variable
;
;    Ver.5,  6-Jul-2009, Enrico Landi
;        Corrected error in the definition of total_rate
;
;    Ver.6, 30-Oct-2009, Enrico Landi
;        Corrected error in the calculation of the correction factor
;
;-


siz=size(cc)
n_levels=siz[1]


;
; Create a cidata and recdata arrays for ion and rec rates respectively
;
cidata=dblarr(n_levels,n_elements(ionrec.temp))
recdata=cidata
;
IF n_elements(ionrec.ci) gt 1 THEN BEGIN
  cidata(ionrec.lev_up_ci-1,*)=ionrec.ci(*,*)
ENDIF
;
IF n_elements(ionrec.rec) gt 1 THEN BEGIN
  recdata(ionrec.lev_up_rec-1,*)=ionrec.rec(*,*)
ENDIF


;
; Note: ionrec.ioneq has dimensions [nt,3], where nt is the no. of
; temperatures in the ion balance file, and 3 is for the ions either side
; of the ion of interest.
;
nioneq=n_elements(reform(ionrec.ioneq(*,0)))
;temp_ioneq=4.0+0.1*findgen(nioneq)
temp_ioneq=reform(ionrec.temp_ioneq)


;
; Determines the ion fractions of the chosen ion, and of the two adjacent
; ions
;
ion_low=ion_frac_interp(t,temp_ioneq,reform(ionrec.ioneq[*,0]))
ion_middle=ion_frac_interp(t,temp_ioneq,reform(ionrec.ioneq[*,1]))
ion_high=ion_frac_interp(t,temp_ioneq,reform(ionrec.ioneq[*,2]))


;
; if ion fraction zero then just return the population array
;
IF ion_middle EQ 0. THEN return,pp

frac_low=ion_low/ion_middle
frac_high=ion_high/ion_middle

;
; Calculates the ionization and recombination rates at the temperature T
; Interpolation or extrapolation is performed by the routine ION_REC_INTERP.
; Note that extrapolation is only performed *below* the temperature range
; for ionization, and only *above* the temperature range for recombination.
;
; Note that I check if total(..datak)=0 in case there is no data for that
; transition.
;
rrate=dblarr(n_levels)
crate=dblarr(n_levels)
total_rate=dblarr(n_levels)
FOR k=0,n_levels-1 DO BEGIN
  cidatak=reform(cidata[k,*])
  IF total(cidatak) NE 0. THEN BEGIN
    cirate=ci_rec_interp(t,ionrec.temp,cidatak,/extrap_below)
  ENDIF ELSE BEGIN
    cirate=0d0
  ENDELSE
 ;
  recdatak=reform(recdata[k,*])
  IF total(recdatak) NE 0. THEN BEGIN
    recrate=ci_rec_interp(t,ionrec.temp,recdatak,/extrap_above)
  ENDIF ELSE BEGIN
    recrate=0d0
  ENDELSE
 ;
  total_rate[k]=(cirate*frac_low+recrate*frac_high)*xne
  rrate[k]=recrate
  crate[k]=cirate
ENDFOR


;
; Calculates the total excitation rate to each level, summing all the
; population-weighted excitations from lower levels and cascades from
; higher levels
; 
total_exc=dblarr(n_levels)
FOR k=1,n_levels-2 DO BEGIN
  total_exc(k)=total(cc(0:k-1,k)*pp(0:k-1)) + $
       total(cc(k+1:n_levels-1,k)*pp(k+1:n_levels-1))
ENDFOR
;
; No cascades for highest level!
;
total_exc(n_levels-1)=total(cc(0:n_levels-2,n_levels-1)*pp(0:n_levels-2))


;
; Calculates the correction to each level population, except the ground level
;
correction=dblarr(n_levels)+1.
IF ion_middle gt 0. THEN BEGIN
  correction(1:n_levels-1)=1. + $
       total_rate(1:n_levels-1)/(total_exc(1:n_levels-1))
ENDIF
; 
; Correct the populations of all levels, except the ground level, and
; renormalize to 1
;
new_pp=pp*correction
new_pp=new_pp/total(new_pp)
return,new_pp

END
