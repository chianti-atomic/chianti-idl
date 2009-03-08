
function matrix_solver, xne, aa, qq, ppr, aax, c=c, $
                        frac_cutoff=frac_cutoff

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
;    AA     2-D array of transition probabilities.
;
;    QQ     2-D array of electron excitation rates.
;
;    PPR    2-D array of proton rates
;
;    AAX    2-D array of photoexcitation/stimulated emission rates.
;
; OPTIONAL INPUT
;
;    FRAC_CUTOFF The fraction of non-zero elements in the C matrix below
;                which the sparse matrix solver is used. The default value
;                is zero (i.e., don't use sparse matrix solver).
;
; OPTIONAL OUTPUT
;
;    C      The 2D matrix containing all rate coefficients.
;
; OUTPUT
;
;    A 1-D array of level populations, scaled so that the sum of the 
;    populations is 1.
;
; HISTORY
;
;    Ver.1, 1-Jul-2005, Peter Young
;
;    Ver.2, 10-Mar-2006, Peter Young
;        commented out warning message about status as the level populations
;        are still accurate to < 0.001% when this occurs and so it's not
;        neccessary to warn users.
;-

;
; FRAC_CUTOFF determines when the sparse matrix solver is used. E.g., if 
; 0.3 then the C array must have less than 30% non-zero elements before the 
; sparse matrix routine is used.
;
; frac_cutoff=0 -> don't use sparse matrix routine
;
IF n_elements(frac_cutoff) EQ 0 THEN frac_cutoff=0.0

siz=size(aa)
n_levels=siz[1]
b=dblarr(n_levels)

c=aa + xne*(qq+ppr) + aax ; creates the c matrix

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
ENDIF ELSE BEGIN
  pp1=invert(transpose(c),status)#b
;  IF status NE 0 THEN BEGIN
;    print,'%MATRIX_SOLVER: status='+trim(status)+' inversion error'
;  ENDIF
ENDELSE

pp=pp1/total(pp1)               ; force total population to be 1

return,pp

END
