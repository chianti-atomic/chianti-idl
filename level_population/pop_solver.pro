

PRO POP_SOLVER, input, T, XNE, POP, N_LEVELS=N_LEVELS, data_str=data_str, $
                sum_mwl_coeffs=sum_mwl_coeffs, radfunc=radfunc, $
                frac_cutoff=frac_cutoff

;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an atomic database package for the calculation of
;       astrophysical emission line spectra.  It is a collaborative project
;       involving Ken Dere (Naval Research Laboratory, Washington DC), 
;       Brunella Monsignori-Fossi and Enrico Landi (Arcetri Observatory, 
;       Florence), and Helen Mason and Peter Young (DAMTP, Cambridge Univ.).
;
;
; NAME:  POP_SOLVER
;       
; PURPOSE:
;
;	To solve the level balance equations for Chianti ions.
;
; CATEGORY:
;
;       Scientific analysis
; 
; EXPLANATION:
;
;	This routine solves the level balance equations for the CHIANTI ions. 
;       Atomic data is pre-loaded into the COMMON blocks, and so POP_SOLVER 
;       can only be called indirectly through other routines.
;
;       The matrix equation Ax=b is solved where A contains all the atomic 
;       data (electron rate coefficients, radiative decay rates, proton rate 
;       coefficients, photoexcitation rates), x are the level populations, 
;       and b a vector set to zeros except for the first element which is 1.
;
;       To solve the matrix equation, pop_solver calls out to the CHIANTI
;       routine matrix_solver.
;
;       The matrix A is created from the atomic data in the COMMON blocks. 
;       In order to optimise POP_SOLVER, A is created where possible through 
;       array operations rather than FOR loops.
;
;       With v.5 of CHIANTI the additional atomic processes of ionization
;       and recombination can be included when calculating the level
;       populations. These processes are not included in the matrix A.
;       Instead the level populations x are 'corrected' for ionization and
;       recombination afterwards. This correction is performed by the routine
;       correct_pops. More details of this method are found in the CHIANTI
;       v.5 paper.
;
; CALLING SEQUENCE:
;
;	POP_SOLVER, T, XNE, POP, N_LEVELS=N_LEVELS
;
; INPUTS:
;
;	T	Temperatures, e.g., 10.^6
;
;	XNE	Densities, e.g., 10.^8
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
; OUTPUT:
;
;	POP	An array of level populations of size 
;		n_T x n_XNE x n_levels
;
; OPTIONAL OUTPUTS
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
;                .cc          Electron rate coefficients (2D array)
;                .ccp         Proton rate coefficients (2D array)
;                .ion_rate    Ionization rate (1D array)
;                .rec_rate    Recombination rate (1D array)
;                .correction  Correction factor for level pop (1D array)
;                .frac_low    Ratio of N+1 ion fraction to N (scalar)
;                .frac_high   Ratio of N-1 ion fraction to N (scalar)
;
;                The 2D arrays are such that, e.g., aa[0,20] 
;                corresponds to an excitation, while aa[20,0] is a 
;                de-excitation.
;
;                The 1D arrays are simply the rate coefficients into the
;                individual levels.
;
; PROGRAMMING NOTES:
;
;       PROTON RATES
;       ------------
;       To include the proton rates, it is necessary to have the 
;       proton-to-electron ratio. This needs to be calculated before the 
;       call to pop_solver, and the resulting ratio(s) passed through 
;       'pe_ratio' in the common block 'proton'.
;
;       Note that there is no keyword to switch off proton rates (i.e., 
;       no /NOPROT keyword). To switch off proton rates, it is necessary 
;       to set pstr=-1. This should be done by the calling routine.
;
;
; COMMON BLOCKS:
;
;	None.
;
; CALLS:
;
;	DESCALE_ALL, PROTON_DENS(), MATRIX_SOLVER(), CORRECT_POPS()
;
; HISTORY:
;
;	Ver 1, PRY 29-Mar-99
;	Ver 2, PRY 30-Apr-99, added call to get_prot_rates
;	Ver 3, PRY 15-Dec-99, added deu to upsilon common block in order 
;		to be consistent with the main Chianti routines.
;	Ver 4, PRY 9-May-00, corrected problem with threshold when dealing 
;		with sparse matrices. Basically values less than 1.e-30 in 
;		the c-matrix were being set to zero and giving rise to 
;		NaN's in certain circumstances.
;       Ver.5, PRY 14-Jul-00, changed elvl common block to the elvlc common 
;               block which is now the Chianti standard. Also, when 
;               descaling upsilons, the routine now uses the Delta-E from 
;               the .splups file.
;       Ver.6, PRY 9-Aug-00, changed routine to deal better with the 
;               dielectronic recombination files
;       Ver.7, PRY 17-Aug-00, routine does not call LINBCG now if radtemp 
;               is non-zero.
;       Ver.8, PRY 29-Aug-00, the sparse matrix section has been disabled.
;       Ver.9, PRY 12-Nov-01, calls routine proton_dens() to calculate the 
;               proton to electron ratio.
;       Ver.10, PRY, 6-Dec-01, corrected bug when there are more levels 
;               in .splups file than in .elvlc file (ZnXXV).
;       Ver.11, PRY, 11-Jul-02, removed ION keyword
;       Ver.12, PRY, 9-Aug-02, within the equation solving section, I've set 
;               the population of the ground level (rather the n_level level) 
;               to 1, and this seems to stop negative populations appearing 
;               in extreme conditions.
;       Ver.12, PRY, 21-Aug-02, changed exp(-1/1/a) to exp(-a) in electron
;               excitation section which caused a hang-up in some 
;               circumstances. Also, the routine now uses vector ECMC 
;               (combined experimental and theoretical energies) in 
;               determining if a level lies above or below another level. 
;               Previously only used the observed energy vector. Also, the 
;               exponential in the electron excitation section now uses the 
;               (accurate) .elvlc energy separation rather than the .splups 
;               energy separation, which can cause significant (~20-30%) 
;               differences in level populations of high-lying levels at 
;               low temperatures.
;       Ver.13, PRY, 10-Sep-02, corrected bug for proton rates. The excitation 
;               and de-excitation rates were being swapped.
;
;       V. 14  4-Oct-2003  Giulio Del Zanna (GDZ).
;               -removed all COMMON blocks (note that only proton_dens.pro has
;                one: COMMON elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref)
;               -only the essential information input is passed to the routine
;                via a new input structure.
;               -fixed a bug, that affected all the satellite lines, and was 
;                introduced in v.12,  included in  CHIANTI v.4.0.
;                basically the ionization potential was not subtracted when
;                calculating the Delta E in the exponential.
;
;       V. 15  7-Oct-2004  Enrico Landi (EL)
;               Included ionization and recombination as level population 
;               processes.
;
;       V. 16  6-Apr-2005  Enrico Landi (EL)
;               Included extrapolation of ionization and recombination rates
;               for temperatures beyond those provided in the .ci and .rec
;               files.
;
;       V. 17  10-Jun-2005  Peter Young
;               Tidied up code, introduced call to correct_pops for
;               ionization/recombination, and added radfunc= and
;               sum_mwl_coeffs= keywords
;
;       V. 18  12-Jul-2005, Peter Young
;               Improved implementation of RADFUNC keyword
;
;       V. 19  27-Jul-2005, Peter Young
;               Corrected bug when the ionrec structure does not exist.
;
;       V. 20  16-Aug-2005, Peter Young
;               Routine now catches any NaN values in the level populations
;               and prints a warning. All pops are set to zero in this case.
;
;       V. 21  1-Aug-2006, Enrico Landi
;               Changed the way recombination is handled for the He-like ions.
;               Collision rates are redefined including recombination rates
;               (plus cascades) for the levels for which they are available,
;               scaled by the ratio of the abundances of the recombining and
;               recombined ion.
;
;       V. 22  23-Apr-2013, Peter Young
;               For level-resolved recombination for He-like ions (see
;               v.21), the temperature array used for the ion
;               fractions (temp_ioneq) was hard-coded rather than
;               being obtained from the INPUT structure. This has been
;               fixed. 
;
; VERSION     : 22, 23-Apr-2013
;
;-

;these are required:
;-------------------

 gname = input.gname
 jj= input.jj
 ecm= input.ecm
 ecmth= input.ecmth
 wvl= input.wvl
 a_value= input.a_value
 splstr= input.splstr

;optional ones: 
;-------------------

IF tag_exist(input, 'radtemp') THEN  radtemp= input.radtemp
IF tag_exist(input, 'dilute') THEN  dilute= input.dilute
IF tag_exist(input, 'prot_struc') THEN  prot_struc= input.prot_struc
IF tag_exist(input, 'pe_ratio') THEN  pe_ratio= input.pe_ratio
IF tag_exist(input, 'ionrec') THEN BEGIN
  ionrec_struc= input.ionrec
  status=ionrec_struc.status
ENDIF ELSE BEGIN
  status=0
ENDELSE


IF n_params(0) LT 4 THEN BEGIN
   print,' use>  pop_solver,input, temperature,density,populations, $'
   print,'                   [n_levels=n_levels, data_str=data_str] '
   return
ENDIF

convertname,gname,iz,ion
ion2spectroscopic,gname,snote, dielectronic=dielectronic


IF dielectronic THEN BEGIN
  read_ip,concat_dir(concat_dir(!xuvtop, 'ip'), 'chianti.ip'),ionpot,ipref
  ip=ionpot(iz-1,ion-1)
ENDIF ELSE ip=0


xne = DOUBLE(xne)
t = DOUBLE(t)
;
; need the following to turn t into an array if it only has 1 element
;
IF n_elements(t) EQ 1 THEN BEGIN
  t0=t
  t=dblarr(1)
  t[0]=t0
ENDIF

ecmc=ecm
ind=where(ecm EQ 0.)
IF ind[0] NE -1 THEN ecmc[ind]=ecmth[ind]

mult=2.*jj+1.
;
hck=1.98648d-16/1.38062d-16
ryd2cm=109737.31534d
;
n_elvl=n_elements(ecm)
wsize=size(wvl)
n_wgfa1=wsize(1)
n_wgfa2=wsize(2)
usize=max([splstr.lvl1,splstr.lvl2])
;
IF N_ELEMENTS(n_levels) EQ 0 THEN n_levels=min([n_elvl,n_wgfa1,n_wgfa2,usize])
;
;
c=dblarr(n_levels,n_levels)
d=dblarr(n_levels,n_levels)
b=dblarr(n_levels)

diag=dblarr(n_levels)

nt=N_ELEMENTS(t)       ; no. of temperatures
nxne=N_ELEMENTS(xne)   ; no. of densities

IF n_elements(sum_mwl_coeffs) EQ 0 THEN BEGIN
  sum_mwl_coeffs=dblarr(nt)+1.
  sumtst=0
  pop=dblarr(nt,nxne,n_levels)
ENDIF ELSE BEGIN
  sumtst=1
  pop=dblarr(1,nxne,n_levels)
  IF nt NE n_elements(sum_mwl_coeffs) THEN BEGIN
    print,'%POP_SOLVER: number of temperatures must match the size of '+$
         'SUM_MWL_COEFFS.'
    print,'             Populations not calculated.'
    return
  ENDIF
ENDELSE
;pop=dblarr(nt,nxne,n_levels)

;;------------------------------[]
; The arrays
;
; e.g., aa(0,19) will be zero (no 0 -> 19 A value)
;       aa(19,0) will be non-zero
;
;       qq(0,19) electron excitation
;       qq(19,0) electron de-excitation
;;------------------------------[]

ident=make_array(n_levels,val=1.)  ; use for making de arrays
ecmn=ecmc(0:n_levels-1)        ; n_levels
den=ecmn#ident & dem=ident#ecmn

aat=a_value(0:n_levels-1,0:n_levels-1)        ; transpose of aa

aa=DOUBLE(TRANSPOSE(aat))

aax=c

IF N_ELEMENTS(dilute) EQ 0 THEN dilute=0.

;;------------------------------------------------------------------[-]
; The following loads up the photoexcitation (pexc) and stimulated 
; emission (stem) arrays)
;
IF dilute NE 0. THEN BEGIN
  stem=c & pexc=c & ede=c
 ;
  multn=mult[0:n_levels-1]      ; in case mult and ecm are bigger than 
 ;
  mm=TRANSPOSE(multn#(1/multn))      ; needed for photoexcitation
 ;
 ; Note below that there's a factor of 1d8^3. This is because the units
 ; of lambda are angstroms, but I need a cm^-3 in the units of the energy
 ; density.
 ;
  IF n_elements(radfunc) NE 0 THEN BEGIN
    en_rf=abs(den-dem)
    i=where(en_rf EQ 0.)   ; the index i prevents underflows
    IF i[0] NE -1 THEN en_rf[i]=1d50    ; set to arbitrarily large value
   ;
    lambda=1d8/en_rf
   ;
    bits=str_sep(radfunc,',')
    CASE n_elements(bits) OF
      1: result=call_FUNCTION(radfunc,lambda)
      2: BEGIN
        rfunc=bits[0]
        a1=double(bits[1])
        result=call_FUNCTION(rfunc,lambda,a1)
      END
      ELSE: BEGIN
        rfunc=bits[0]
        a1=double(bits[1])
        a2=double(bits[2])
        result=call_FUNCTION(rfunc,lambda,a1,a2)
      END
    ENDCASE
   ;
    IF i[0] NE -1 THEN result[i]=0d0
   ;
    result=result*lambda^5/(1d8)^3/8d0/!pi/1.986d-8
  ENDIF ELSE BEGIN
    dd=ABS(den-dem)*hck/radtemp
   ;
   ; the following lines are necessary to prevent infinities and floating
   ; underflow errors
   ;
    dd=dd < 150.
    i=where(dd EQ 0.)
    j=where(dd LE 1d-15 AND dd NE 0.)
    k=where(dd GE 1d-15)
   ;
    ede[k]=exp(dd[k]) - 1.
    IF j[0] NE -1 THEN ede[j]=dd[j]
    ede[i]=1d50     ; arbitrarily large value
   ;
    result=1d0/ede
    result[i]=0d0   ; set i's to zero since A-values are zero (thus don't
                    ; contribute to pexc and stem
  ENDELSE
 ;
  ind=WHERE( (aat NE 0.) AND (result NE 0.) )
  IF ind[0] NE -1 THEN pexc[ind]=aat[ind]*dilute*mm[ind]*result[ind]
 ;
  ind=where( (aa NE 0.) AND (result NE 0.) )
  IF ind[0] NE -1 THEN stem[ind]=aa[ind]*dilute*result[ind]
 ;
  aax=pexc+stem
ENDIF
;;------------------------------------------------------------------[-]

;________________
; create a ppr array for the proton rates
;
ppr=MAKE_ARRAY(n_levels,n_levels,nt,/double)
IF n_tags(prot_struc) NE 0 THEN BEGIN
 ;
  IF (n_elements(pe_ratio) NE nt) THEN BEGIN
    print,'%POP_SOLVER: WARNING, pe_ratio size does not match temp'
    print,n_elements(pe_ratio),nt
    pe_ratio=proton_dens(alog10(t))
  ENDIF
 ;
  FOR i=0,n_elements(prot_struc)-1 DO BEGIN
    l1=prot_struc[i].lvl1-1
    l2=prot_struc[i].lvl2-1
    de=ABS(prot_struc[i].de)
    descale_all,t,prot_struc,i,prate
    IF ecmc(l1) LT ecmc(l2) THEN BEGIN
      ppr[l1,l2,*]=prate*pe_ratio*sum_mwl_coeffs
      ppr[l2,l1,*]=prate*pe_ratio*mult[l1]/mult[l2]* $
           exp(de*13.61/8.617/10.^(-5)/t)*sum_mwl_coeffs
    ENDIF ELSE BEGIN
      ppr[l2,l1,*]=prate*pe_ratio*sum_mwl_coeffs
      ppr[l1,l2,*]=prate*pe_ratio*mult[l2]/mult[l1]* $
           exp(de*13.61/8.617/10.^(-5)/t)*sum_mwl_coeffs
    ENDELSE
  ENDFOR
ENDIF

;______________
; Create a qq array for electron rates
;
qq=MAKE_ARRAY(n_levels,n_levels,nt,/double)
;
l1=splstr.lvl1-1
l2=splstr.lvl2-1

;GDZ- added ip
kte=(hck*abs(ecmc[l1]-(ecmc[l2]-ip))) # (1d0/t)
;******************************************


; kte=(hck*ryd2cm*ABS(splstr.de)) # (1d0/t)


ind_pos=where(ecmc[l2] GT ecmc[l1])
ind_neg=where(ecmc[l2] LT ecmc[l1])
xx=dblarr(n_elements(splstr),nt)
yy=xx
;
; xx and yy contain all factors in the expression for the rate coefficient, 
; except for the upsilon. They can be generated using array operations - the 
; upsilons need a for loop.
;
IF ind_neg[0] NE -1 THEN BEGIN
  xx[ind_neg,*]=(8.63d-6/(mult[l1[ind_neg]])#(1./sqrt(t)))
  yy[ind_neg,*]=8.63d-6* exp(-kte[ind_neg,*]) * $
       1./( mult[l2[ind_neg]] # sqrt(t) )
ENDIF
IF ind_pos[0] NE -1 THEN BEGIN
  yy[ind_pos,*]=(8.63e-6/mult[l2[ind_pos]]) # (1./sqrt(t))
  xx[ind_pos,*]=8.63e-6* exp(-kte[ind_pos,*]) * $
       1./(mult[l1[ind_pos]] # sqrt(t))
ENDIF
;
; this is the for loop for the upsilons
;
FOR i=0,n_elements(splstr)-1 DO BEGIN
  IF (l1[i] LE n_levels-1) AND (l2[i] LE n_levels-1) THEN BEGIN
    descale_all,t,splstr,i,ups
; EL MODIFICATIONS START HERE
    IF status GT 0 and iz-ion eq 1 THEN BEGIN                ; Treats rec in He-like ions as excitation mechanism
      rate_lev_up=input.ionrec.lev_up_rec
      it_exists=where(l2[i] eq rate_lev_up-1 and l1[i] eq 0)
      IF it_exists(0) GE 0 THEN BEGIN 
        rate=input.ionrec.rec
        rate_temp=input.ionrec.temp
        ioneq_frac=input.ionrec.ioneq
        nioneq=n_elements(reform(ioneq_frac(*,0)))
        temp_ioneq=input.ionrec.temp_ioneq
;        temp_ioneq=4.0+0.1*findgen(nioneq)    ; this was changed in v.22
        rec_rate=dblarr(nt)
        ratio_ioneq=dblarr(nt)
        FOR k=0,nt-1 DO BEGIN 
          rec_rate(k)=ci_rec_interp(t(k),rate_temp,reform(rate(it_exists,*)),/extrap_above)
          ion_high=ion_frac_interp(t(k),temp_ioneq,reform(ioneq_frac(*,2)))
          ion_middle=ion_frac_interp(t(k),temp_ioneq,reform(ioneq_frac(*,1)))
          IF ion_middle gt 0 then ratio_ioneq(k)=ion_high/ion_middle ELSE ratio_ioneq(k)=0
          eff_rec_rate=ratio_ioneq*rec_rate
        ENDFOR
        qq[l1[i],l2[i],*]=xx[i,*]*ups*sum_mwl_coeffs+eff_rec_rate
        qq[l2[i],l1[i],*]=yy[i,*]*ups*sum_mwl_coeffs
      ENDIF ELSE BEGIN 
        qq[l1[i],l2[i],*]=xx[i,*]*ups*sum_mwl_coeffs
        qq[l2[i],l1[i],*]=yy[i,*]*ups*sum_mwl_coeffs
      ENDELSE
    ENDIF ELSE BEGIN
      qq[l1[i],l2[i],*]=xx[i,*]*ups*sum_mwl_coeffs
      qq[l2[i],l1[i],*]=yy[i,*]*ups*sum_mwl_coeffs
    ENDELSE
; EL MODIFICATIONS END HERE
  ENDIF
ENDFOR
;______________

IF sumtst EQ 0 THEN BEGIN
  FOR j=0,nt-1 DO BEGIN
    qqx=qq[*,*,j]
    FOR i=0,nxne-1 DO BEGIN
      pp=matrix_solver(xne[i],aa,qqx,ppr[*,*,j],aax,c=c)
; EL MODIFICATIONS START HERE
      IF status GT 0 and iz-ion ne 1 THEN BEGIN              ; Includes rec/ion corrections if ion is not He-like
; EL MODIFICATIONS END HERE
        pp=correct_pops(pp,t[j],xne[i],input.ionrec,c, $
                        rrate=rrate,crate=crate,correction=correction, $
                        frac_low=frac_low, frac_high=frac_high)
      ENDIF
      pop[j,i,*]=pp
    ENDFOR
  ENDFOR
ENDIF ELSE BEGIN          ; non-maxwellians
  qqx=total(qq,3)
  ppx=total(ppr,3)
  FOR i=0,nxne-1 DO BEGIN
    pp=matrix_solver(xne[i],aa,qqx,ppx,aax,c=c)
    pop[0,i,*]=pp
  ENDFOR
ENDELSE

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
IF nt EQ 1 AND nxne EQ 1 THEN BEGIN
  IF n_elements(rrate) EQ 0 THEN rrate=0.
  IF n_elements(crate) EQ 0 THEN crate=0.
  IF n_elements(correction) EQ 0 THEN correction=0.
  IF n_elements(frac_low) EQ 0 THEN frac_low=0.
  IF n_elements(frac_high) EQ 0 THEN frac_high=0.
  data_str={aa: aa, aax: aax, cc: xne[0]*qq[*,*,0],  $
            ccp: xne[0]*ppr[*,*,0], $
            rec_rate: rrate, ion_rate: crate, $
            correction: correction, $
           frac_low: frac_low, frac_high: frac_high}
ENDIF


END
