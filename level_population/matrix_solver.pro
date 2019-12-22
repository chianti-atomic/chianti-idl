
function matrix_solver, xne, rates, index_t=index_t, c=c, $
                        frac_cutoff=frac_cutoff, verbose=verbose, sumtst=sumtst

;+
; NAME
;
;    MATRIX_SOLVER()
;
; PROJECT
;
;    CHIANTI
;
; EXPLANATION
;
;    Takes the matrices for the various atomic processes employed in
;    CHIANTI and returns the level populations.
;
; INPUTS
;
;    XNE    Electron density (cm^-3). Scalar.
;
;    index_t 
;
;    RATES. A structure with the collected input rates:
;
;    AA     2-D array of transition probabilities.
;
;    AAX    2-D array of photoexcitation/stimulated emission rates.
;
;    QQ     3-D array (n_temperatures,n_levels,n_levels) of electron
;           excitation and de-excitation rates.
; 
;    PPR    3-D array (n_temperatures,n_levels,n_levels) of proton
;           excitation and de-excitation rates 
;
; OPTIONAL:  in the case that the matrix is composed of the rates for
;           two ions, these additional matrices are included in the
;           structure: 
;
;    RR     Radiative recombination rates from the higher ionization
;            stage.
;
;    IONIZ  Ionization rates from the lower to the higher ionisation
;           stage. 
;
;    AI     2-D array of autoionization rates from the lower to the
;           higher ionisation stage. 
;
;    DC     3-D array (n_temperatures,n_levels,n_levels) of
;           dielectronic capture rates 
;
;    DR     3-D array (n_temperatures,n_levels,n_levels) of
;           dielectronic recombination rates. For v.9, this includes
;           only the DR between the ground states, with the
;           contributions from the autoiomnizing states removed to
;           avoid double counting.
;
;
; OPTIONAL INPUT
;
;    FRAC_CUTOFF The fraction of non-zero elements in the C matrix below
;                which the sparse matrix solver is used. The default value
;                is zero (i.e., don't use sparse matrix solver).
;
; KEYWORDS
;
;    VERBOSE   If set, then some informational messages will be
;              printed to the IDL window.
;
;   SUMTST     If set, sum the rates for the non-Maxwellian calculation
;
; OPTIONAL OUTPUT
;
;
; OUTPUT
;
;    A 1-D array of level populations, scaled so that the sum of the 
;    populations is 1.
;
; HISTORY
;
;    Ver.1, 1-Jul-2005, Peter Young (PRY)
;
;    Ver.2, 10-Mar-2006, Peter Young
;        commented out warning message about status as the level populations
;        are still accurate to < 0.001% when this occurs and so it's not
;        neccessary to warn users.
;
;    Ver.3, 10-Jan-2014, Peter Young
;        Added /verbose keyword and some informational messages.
;    Ver.4, June 2018, PRY, modified input so it is given as a single structure.
;    Ver.5, 20-Jul-2018 Giulio Del Zanna (GDZ), minor modifications and
;            write-up of the header; added the matrix C as optional
;            output, used by correct_pops. 
;    ver.6, 5 Dec 2018, GDZ, added sumtst 
;
;-

if keyword_set(sumtst) then begin 
; in the case of non-Maxwellians, the rates have been summed in T already.

; add the electron and proton rates:
c_coll=rates.qq[*,*] + rates.ppr[*,*]

; add the radiative recombination:
IF tag_exist(rates,'rr') THEN c_coll=c_coll+rates.rr[*,*] 
; add the dielectronic recombination:
IF tag_exist(rates,'dr') THEN c_coll=c_coll+rates.dr[*,*]
; add dielectronic capture: 
IF tag_exist(rates,'dc') THEN c_coll=c_coll+rates.dc[*,*] 
; add the collisional ionization:
IF tag_exist(rates,'ioniz') THEN c_coll=c_coll+rates.ioniz[*,*] 


endif else begin 

IF n_elements(index_t) EQ 0 THEN index_t=0
;
; Add all the rate coefficients for the processes that are
; proportional to the electron density
;
; add the electron and proton rates:
c_coll=rates.qq[index_t,*,*] + rates.ppr[index_t,*,*]

; add the radiative recombination:
IF tag_exist(rates,'rr') THEN c_coll=c_coll+rates.rr[index_t,*,*] 
; add the dielectronic recombination:
IF tag_exist(rates,'dr') THEN c_coll=c_coll+rates.dr[index_t,*,*]
; add dielectronic capture: 
IF tag_exist(rates,'dc') THEN c_coll=c_coll+rates.dc[index_t,*,*] 
; add the collisional ionization:
IF tag_exist(rates,'ioniz') THEN c_coll=c_coll+rates.ioniz[index_t,*,*] 


end


;
; Add all the rate coefficients for the processes that are
; *not* proportional to the electron density
 
; Add A-values and photo-excitation, plus de-xcitation processes:
c_noncoll=rates.aa + rates.aax
; Add autoionization rates 
IF tag_exist(rates,'ai') THEN c_noncoll=c_noncoll+rates.ai

; Add the two stes of rates 
c=c_noncoll + xne*c_coll


;
; FRAC_CUTOFF determines when the sparse matrix solver is used. E.g., if 
; 0.3 then the C array must have less than 30% non-zero elements before the 
; sparse matrix routine is used.
;
; frac_cutoff=0 -> don't use sparse matrix routine
;
IF n_elements(frac_cutoff) EQ 0 THEN frac_cutoff=0.0

n_levels=rates.n_levels
b=dblarr(n_levels)

;c=aa + xne*(qq+ppr) + aax ; creates the c matrix

ind=where(c NE 0.)
frac_not_zero=float(n_elements(ind))/float(n_elements(c)) 

diag = -total(c,2)
c[findgen(n_levels),findgen(n_levels)] = diag

c(*,0)=0d0                      ; set this row to zeros...
c(0,0)=1d0                      ; ...except for ground level

b(0)=1d0                        ; b is zero except for first element


;
; Only use sparse matrix routine if we have more than 100 levels in ion
; model.
;
IF frac_not_zero LE frac_cutoff AND n_levels GT 100 THEN BEGIN
  pp1=LINBCG(SPRSIN(c,thresh=1d-30),b,b)
  IF keyword_set(verbose) THEN BEGIN
    print,'% MATRIX_SOLVER: sparse matrix method has been used.'
    print,'                 frac_not_zero=',string(format='(f5.2)',frac_NOT_zero)
    print,'                   frac_cutoff=',string(format='(f5.2)',frac_cutoff)
  ENDIF 
ENDIF ELSE BEGIN
  pp1=invert(transpose(c),status)#b
  IF keyword_set(verbose) THEN BEGIN
    print,'% MATRIX_SOLVER: standard matrix method has been used.'
    print,'                 frac_not_zero=',string(format='(f5.2)',frac_NOT_zero)
    print,'                   frac_cutoff=',string(format='(f5.2)',frac_cutoff)
    print,'                      n_levels=',n_levels
   ;    
    IF status NE 0 THEN BEGIN
      print,'%MATRIX_SOLVER: status='+trim(status)+' inversion error'
    ENDIF
  ENDIF 
ENDELSE

pp=pp1/total(pp1)               ; force total population to be 1

return,pp

END
