
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
;     from pop_solver.pro. CHIANTI Technical Report 15 gives more
;     details about the procedure performed here.
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
;               is obsolete.] 
;     SPARSE:   If set, then the sparse matrix inversion routine
;               (linbcg) will be used. [This is obsolete.]
;     LAPACK:   If set, then the LAPACK matrix inversion routine
;               (la_invert) will be used. [This is obsolete.]
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
;     Ver.9, 22-May-2023, Peter Young
;            Now checks if there are any levels without a population mechanism,
;            and these are removed from the matrix equations. (They are reinstated
;            with zero populations after the equations are solved.) The invert
;            routine is no longer used for solving the equations. Instead, the
;            routines la_ludc and la_lusolve are used. The options /lapack,
;            /sparse, /regular and frac_cutoff= have been removed, although the
;            keywords are still present.
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

ENDELSE 


;
; Add all the rate coefficients for the processes that are
; *not* proportional to the electron density
 
; Add A-values and photo-excitation, plus de-xcitation processes:
c_noncoll=rates.aa + rates.aax
; Add autoionization rates 
IF tag_exist(rates,'ai') THEN c_noncoll=c_noncoll+rates.ai


; Add the two sets of rates 
c=c_noncoll + xne*c_coll

;
; n_levels is the number of levels for which the rates are defined.
;
n_levels=rates.n_levels
index=indgen(n_levels)

;
; The code below reduces the size of the c matrix if there are
; levels with no population mechanism. Such levels have
; total(c[*,i])=0.
;
; Iterations are required because, some levels may be have
; population mechanism rates (typically radiative decays from
; above), but the upper levels are actually zero-population levels.
; The iterations are performed until all of the zero-rate levels are removed.
;
; nlev is the number of levels of the final reduced c-matrix.
;
swtch=0
nlev_save=n_levels
count=1
WHILE swtch EQ 0 DO BEGIN
  c_chck=total(c,1)
  k=where(c_chck NE 0.,nlev)
  IF nlev EQ nlev_save OR nlev EQ 1 THEN BREAK
  IF keyword_set(verbose) THEN message,/info,/cont,'Removing zero-pop levels. Iteration '+trim(count)+': reducing levels from '+trim(nlev_save)+' to '+trim(nlev)+'.'
  karr=k#make_array(nlev,value=1.)
  karr_t=transpose(karr)
  c=c[karr,karr_t]
  index=index[k]
  nlev_save=nlev
  count=count+1
ENDWHILE

;
; It's possible that the above iteration gets rid of all the levels (!).
; This can happen for very low temperatures where the excitation rates
; become numerically zero. In this case I just set the ground level to
; 1 and exit.
;
IF nlev EQ 1 THEN BEGIN
  pp=dblarr(rates.n_levels)
  pp[0]=1.0
  return,pp
ENDIF 

b=dblarr(nlev)

;
; The array c has diagonal elements that are all zero. The array c_diag has
; non-zero diagonal values, and is used for the matrix inversion.
;
diag = -total(c,2)
c[findgen(nlev),findgen(nlev)] = diag

c[*,0]=0d0                      ; set this row to zeros...
c[0,0]=1d0                      ; ...except for ground level

b(0)=1d0                        ; b is zero except for first element


;
; Solve the linear equations (c.pp1 = b)
; 
c_ludc=c
la_ludc,c_ludc,ind,/double,status=status
IF status GT 0 AND keyword_set(verbose) THEN message,/info,/cont,'matrix has a zero diagonal element'
pp1=la_lusol(c_ludc,ind,b,/double)



;
; Insert populations into the full array and normalize so sum is 1.
;
pp=dblarr(n_levels)
pp[index]=pp1/total(pp1)


return,pp

END
