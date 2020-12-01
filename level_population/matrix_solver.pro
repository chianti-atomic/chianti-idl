
function matrix_solver, xne, rates, index_t=index_t, c=c, $
                        frac_cutoff=frac_cutoff, verbose=verbose, sumtst=sumtst, $
                        sparse=sparse, lapack=lapack, regular=regular

;+
; NAME:
;     MATRIX_SOLVER()
;
; PURPOSE:
;     Takes the matrices for the various atomic processes employed in
;     CHIANTI and returns the level populations. Intended to be called
;     from pop_solver.pro.
;
; CATEGORY:
;     CHIANTI; level populations.
;
; CALLING SEQUENCE:
;     Result = MATRIX_SOLVER( Xne, Rates )
;
; INPUTS:
;     Xne:   Electron density (cm^-3). Scalar.
;     Rates: A structure containing the rates for the various physical
;            processes. The tags are:
;             .aa  Radiative decay rates (2D).
;             .aax Photoexcitation/stimulated emission (2D).
;             .qq  Electron excitation rate coefficients (3D).
;             .ppr Proton excitation rate coefficients (3D).
;             .rr  Radiative recombination rate coefficients (3D).
;             .dr  Dielectronic recombination rate coefficients (3D).
;             .ioniz Ionization rate coefficients (3D).
;             .dc  Dielectronic capture rate coefficients (3D).
;             .ai  Autoionization rates (2D).
;           The 3D arrays are temperature*levels*levels.
;
; OPTIONAL INPUTS:
;     Index_T: Integer giving index of temperature array. If not set,
;              then 0 is used.
;     Frac_Cutoff:  The fraction of non-zero elements in the C matrix
;                   below which the sparse matrix solver is used. The
;                   default value is zero (i.e., don't use
;                   sparse matrix solver). [This is currently disabled.]
;	
; KEYWORD PARAMETERS:
;     VERBOSE:  If set, then some informational messages will be
;               printed to the IDL window.
;     SUMTST:   If set, sum the rates for the non-Maxwellian
;               calculation.
;     REGULAR:  If set, then the regular inversion method will be used
;               (i.e., the IDL invert routine). This takes precedence
;               over /sparse and /lapack if these are also set. [This
;               is the default.] 
;     SPARSE:   If set, then the sparse matrix inversion routine
;               (linbcg) will be used.
;     LAPACK:   If set, then the LAPACK matrix inversion routine
;               (la_invert) will be used. 
;
; OUTPUTS:
;     A 1-D array of level populations, scaled so that the sum of the 
;     populations is 1.
;
; MODIFICATION HISTORY:
;     Ver.1, 01-Jul-2005, Peter Young (PRY)
;     Ver.2, 10-Mar-2006, Peter Young
;        commented out warning message about status as the level populations
;        are still accurate to < 0.001% when this occurs and so it's not
;        neccessary to warn users.
;     Ver.3, 10-Jan-2014, Peter Young
;        Added /verbose keyword and some informational messages.
;     Ver.4, June 2018, PRY, modified input so it is given as a single structure.
;     Ver.5, 20-Jul-2018 Giulio Del Zanna (GDZ), minor modifications and
;            write-up of the header; added the matrix C as optional
;            output, used by correct_pops. 
;     Ver.6, 5 Dec 2018, GDZ, added sumtst 
;     ver.7, 12 Oct 2020, GDZ, using the LAPACK version by default for
;           large matrices.
;     Ver.8, 28-Oct-2020, Peter Young
;            Added /regular, /sparse and /lapack keywords for
;            switching to different inversion methods. Set the default
;            to /regular.
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
c_coll=reform(rates.qq[index_t,*,*]) + reform(rates.ppr[index_t,*,*])

; add the radiative recombination:
IF tag_exist(rates,'rr') THEN c_coll=c_coll+reform(rates.rr[index_t,*,*])
; add the dielectronic recombination:
IF tag_exist(rates,'dr') THEN c_coll=c_coll+reform(rates.dr[index_t,*,*])
; add dielectronic capture: 
IF tag_exist(rates,'dc') THEN c_coll=c_coll+reform(rates.dc[index_t,*,*])
; add the collisional ionization:
IF tag_exist(rates,'ioniz') THEN c_coll=c_coll+reform(rates.ioniz[index_t,*,*])


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


diag = -total(c,2)
c[findgen(n_levels),findgen(n_levels)] = diag

c(*,0)=0d0                      ; set this row to zeros...
c(0,0)=1d0                      ; ...except for ground level

b(0)=1d0                        ; b is zero except for first element


;
; Set the default option for solving the matrix.
;
dflt=0b   ; regular

;
; SWTCH takes the values
;   0 - use the invert routine (/regular)
;   1 - use the linbcg routine (/sparse)
;   2 - use the la_invert routine (/lapack)
;
; If the number of levels is 100 or less, then always use the regular
; method.
;
; If /regular is set, then it always takes precedence over /sparse and
; /lapack if these are also set.
;
swtch=dflt
CASE 1 OF
   keyword_set(sparse): swtch=1
   keyword_set(lapack): swtch=2
   ELSE: 
ENDCASE
IF keyword_set(regular) THEN swtch=0
IF n_levels LE 100 THEN swtch=0

;
; The following is relevant to the sparse inversion method.
; Note that if more than 30% of elements are not empty, then switch to
; the regular inversion method.
;
IF swtch EQ 1 THEN BEGIN
   sprs_thresh=1d-30
   ind=where(abs(c) GE sprs_thresh)
   frac_not_zero=float(n_elements(ind))/float(n_elements(c))
  ;
   IF frac_NOT_zero GE 0.30 THEN swtch=0
ENDIF ELSE BEGIN
   frac_NOT_zero=-1
ENDELSE 


;
; Only use LAPACK matrix routine if we have more than 100 levels in ion
; model.
;
status=0
CASE swtch OF
   0: BEGIN
      IF keyword_set(verbose) THEN print,'% MATRIX_SOLVER: standard matrix method has been used.'
      pp1=invert(transpose(c),status)#b
   END
  ;
   1: BEGIN
     ;
     ; For low temperatures all of the excitation rates may be below
     ; the threshold, giving rise to NaNs. I check this and then set
     ; all the populations to zero except the ground level.
     ;
      IF keyword_set(verbose) THEN print,'% MATRIX_SOLVER: sparse matrix method has been used.'
      pp1=LINBCG(SPRSIN(c,thresh=sprs_thresh),b,b)
      chck=finite(/nan,pp1)
      IF total(chck) EQ n_levels THEN BEGIN
         pp1=dblarr(n_levels)
         pp1[0]=1.0
      ENDIF 
   END
  ;
   2: BEGIN
      IF keyword_set(verbose) THEN print,'% MATRIX_SOLVER: LAPACK matrix method has been used.'
      pp1=la_invert(transpose(c),double=1,status=status)#b
   END
ENDCASE

IF swtch EQ 0 AND status NE 0 THEN BEGIN
   IF status EQ 2 AND keyword_set(verbose) THEN print,'% MATRIX_SOLVER: warning, small pivot element (invert.pro).'
   IF status EQ 1 THEN print,'% MATRIX_SOLVER: invalid inversion, singular array (invert.pro).'
ENDIF

IF swtch EQ 2 AND status NE 0 THEN BEGIN
   IF status EQ 1 THEN print,'% MATRIX_SOLVER: invalid inversion, singular array (la_invert.pro).'   
ENDIF


IF keyword_set(verbose) AND swtch EQ 1 THEN BEGIN
   print,'% MATRIX_SOLVER: sparse matrix, frac_not_zero = '+string(format='(f5.2)',frac_NOT_zero)
   ;; print,'                 frac_not_zero=',string(format='(f5.2)',frac_NOT_zero)
   ;; print,'                   frac_cutoff=',string(format='(f5.2)',frac_cutoff)
   ;; print,'                      n_levels=',n_levels
  ;    
ENDIF 

;; IF  n_levels GT 100 THEN BEGIN

;;    pp1=la_invert(transpose(c),double=1,status=status)#b

;;    IF status NE 0 THEN BEGIN
;;       print,'%MATRIX_SOLVER: status='+trim(status)+' inversion error !'

;; ; resort to sparse matrix solver      
;;       pp1=LINBCG(SPRSIN(c,thresh=1d-30),b,b)
      
;;    ENDIF

;;    if min(pp1) lt 0. then begin
;;       IF keyword_set(verbose) THEN $
;;          print,'%MATRIX_SOLVER: negative populations set to zero '
    
;;       pp1=pp1>0.
;;    end

 
;;    IF keyword_set(verbose) THEN BEGIN
;;       print,'% MATRIX_SOLVER: LAPACK matrix method has been used.'
;;       print,'                 frac_not_zero=',string(format='(f5.2)',frac_NOT_zero)
;;       print,'                   frac_cutoff=',string(format='(f5.2)',frac_cutoff)
;;    ENDIF 

; IF frac_not_zero LE frac_cutoff AND n_levels GT 100 THEN BEGIN

     
  ;;  pp1=LINBCG(SPRSIN(c,thresh=1d-30),b,b)
  ;; IF keyword_set(verbose) THEN BEGIN
  ;;   print,'% MATRIX_SOLVER: sparse matrix method has been used.'
  ;;   print,'                 frac_not_zero=',string(format='(f5.2)',frac_NOT_zero)
  ;;   print,'                   frac_cutoff=',string(format='(f5.2)',frac_cutoff)
  ;; ENDIF 


;; ENDIF ELSE BEGIN
;;    pp1=invert(transpose(c),status)#b
;;    IF keyword_set(verbose) THEN BEGIN
;;       print,'% MATRIX_SOLVER: standard matrix method has been used.'
;;       print,'                 frac_not_zero=',string(format='(f5.2)',frac_NOT_zero)
;;       print,'                   frac_cutoff=',string(format='(f5.2)',frac_cutoff)
;;       print,'                      n_levels=',n_levels
;;    ;    
;;       IF status NE 0 THEN BEGIN
;;          print,'%MATRIX_SOLVER: status='+trim(status)+' inversion error'
;;       ENDIF
;;    ENDIF 
;; ENDELSE 

pp=pp1/total(pp1)               ; force total population to be 1

return,pp

END
