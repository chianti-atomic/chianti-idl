PRO chianti_ion::POPULATE, T, XNE, keystr = keystr, radtemp = radtemp, dilute = dilute, $   
	N_LEVELS = N_LEVELS, data_str = data_str, sum_mwl_coeffs = sum_mwl_coeffs, $
	radfunc = radfunc, frac_cutoff = frac_cutoff

;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an atomic database package for the calculation of
;       continuum and emission line spectra from astrophysical plasmas. It is a 
;       collaborative project involving George Mason University, the University
;		of Michigan, the Naval Research Laboratory (USA), Cambridge University (UK),
;		the Arcetri Observatory (Italy).
;
;
; NAME:  CHIANTI_ION::POPULATE
;       
; PURPOSE:
;
;	This a developmental copy of pop_solver to solve the level balance equations for Chianti
;   ions for testing with chianti idl objects.
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
;	DESCALE_ALL, PROTON_DENS(), MATRIX_SOLVER()
;
; HISTORY:
;
;
;       Ver.1, 6-jan-2011, Ken Dere
;			first release - provides the populate method for the chianti_ion object
;			derived from pop_solver V.21
;
;
;		V.1, 6-jan-2011, Ken Dere
;-
;
;
;-
IF n_params(0) LT 2 THEN BEGIN
   print,' use>  ->populate, temperature, density, keystr = keystr] '
   return
ENDIF

; radtemp=radtemp, dilute=dilute, N_LEVELS=N_LEVELS, $
; 		data_str=data_str, sum_mwl_coeffs=sum_mwl_coeffs, radfunc=radfunc, $
;                 frac_cutoff=frac_cutoff
;
;these are required:
;-------------------

 gname = self.ionS
 jj= (*self.elvlc).jj
 ecm= (*self.elvlc).ecm
 ecmth= (*self.elvlc).ecmth
 wvl= (*self.wgfa).wvl
 a_value= (*self.wgfa).a_value
 splstr= (*self.splups)

;
;-------------------


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
; zion2filename,z,ion,fname
; ciname = fname+'.cilvl'
; if file_test(ciname) then begin
;     ci=1
;     cidata = 
;  kpd
print,' n_levels = ',n_levels

if ion le iz then begin 
    if ptr_valid(self.reclvl) then begin
		print, ' ptr reclvl is valid'
        rec = 1 
        zion2name,iz,ion+1,uppername
        upper = obj_new('chianti_ion',uppername)
        upper->recombRate,t
        self->ionizRate,t
    endif else begin
		rec = 0
		print, ' ptr reclvl is NOT valid'
	endelse
endif else rec=0


if ion gt 1 then begin
    if ptr_valid(self.cilvl) then begin
		print, ' ptr cilvl is valid'
        ci = 1
        zion2name,iz,ion-1,lowername
        lower = obj_new('chianti_ion',lowername)
        lower->ionizRate,t
        self->recombRate,t
    endif else begin
		ci = 0
		print, ' ptr cilvl is NOT valid'
	endelse
endif else ci = 0

;  if rec == 1, then there are recombination contributions from the higher ionization stage and 
;  it will be included in the level population calculations
;
;  if ci == 1, then there are ionization contributions from the lower ionization stage and 
;  it will be included in the level population calculations


c=dblarr(ci+rec+n_levels,ci+rec+n_levels)
d=dblarr(ci+rec+n_levels,ci+rec+n_levels)
b=dblarr(ci+rec+n_levels)

diag=dblarr(ci+rec+n_levels)

nt=N_ELEMENTS(t)       ; no. of temperatures
nxne=N_ELEMENTS(xne)   ; no. of densities

IF n_elements(sum_mwl_coeffs) EQ 0 THEN BEGIN
  sum_mwl_coeffs=dblarr(nt)+1.
  sumtst=0
  pop=dblarr(nt,nxne,n_levels)
ENDIF ELSE BEGIN
  sumtst=1
  pop=dblarr(1,nxne,n_levels)
  IF nt NE n_element(sum_mwl_coeffs) THEN BEGIN
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

ident=make_array(ci+rec+n_levels,val=1.)  ; use for making de arrays
ecmn=ecmc(0:n_levels-1)        ; n_levels
ecmn = [0.,ecmn, 0.]
den=ecmn#ident & dem=ident#ecmn
aat = dblarr(ci+rec+n_levels,ci+rec+n_levels)
aat[ci,ci]=a_value(0:n_levels-1,0:n_levels-1)        ; transpose of aa

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
ppr=MAKE_ARRAY(ci+rec+n_levels,ci+rec+n_levels,nt,/double)

if ptr_valid(self.psplups) then begin
	pe_ratio = proton_dens(alog10(t))
	print, ' have proton data, pe_ratio = ',pe_ratio
	prot_struc = self.psplups
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
			ppr[ci+l1,ci+l2,*]=prate*pe_ratio*sum_mwl_coeffs
			ppr[ci+l2,ci+l1,*]=prate*pe_ratio*mult[l1]/mult[l2]* $
				exp(de*13.61/8.617/10.^(-5)/t)*sum_mwl_coeffs
			ENDIF ELSE BEGIN
			ppr[ci+l2,ci+l1,*]=prate*pe_ratio*sum_mwl_coeffs
			ppr[ci+l1,ci+l2,*]=prate*pe_ratio*mult[l2]/mult[l1]* $
				exp(de*13.61/8.617/10.^(-5)/t)*sum_mwl_coeffs
			ENDELSE
		ENDFOR
	ENDIF
endif else begin
	print, ' no proton data'
endelse
;______________
; Create a qq array for electron rates
;
qq=MAKE_ARRAY(ci+rec+n_levels,ci+rec+n_levels,nt,/double)
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
      qq[ci+l1[i],ci+l2[i],*]=xx[i,*]*ups*sum_mwl_coeffs
      qq[ci+l2[i],ci+l1[i],*]=yy[i,*]*ups*sum_mwl_coeffs
  ENDIF
ENDFOR
;
;______________________________________________
;
; now include ionization and recombination rates
;
if ci then begin
    cilvl = (*self.cilvl)
    cisize = size(cilvl.temp)
    ncitemp = cisize[1]
    ncirate = cisize[2]
    ;  the interpolated values should be saved somehow
    lowerIoniz = lower->getIonizRate()
    recomb = self->getRecombRate()
    for j=0,nt-1 do begin
        cisum = 0.
        for ici = 0,ncirate-1 do begin
            good = where(cilvl.temp[*,ici] gt 0.)
            lvl2 = cilvl.lvl2[ici]
;             print,' lvl2, temp, rate = ', lvl2, cilvl.temp[good,ici], cilvl.rate[good,ici]
            y2 = spl_init(cilvl.temp[good,ici],alog10(cilvl.rate[good,ici]))
            cirate = spl_interp(cilvl.temp[good,ici],alog10(cilvl.rate[good,ici]),y2,alog10(t[j]))
            cisum = cisum + 10.^cirate
;             print,' lvl2, temp, cirate = ', lvl2, t[j], 10.^cirate
            qq[0,ci+lvl2-1,j] = qq[0,ci,j] + 10.^cirate           
        endfor
        ;  make sure we aren't double-booking
        qq[0,ci,j] = qq[0,ci,j] + lowerIoniz.rate[j] - cisum
        qq[ci,0,j] = qq[ci,0,j] + recomb.rate[j]
;         print, 't, cirate, cisum = ',t[j], lowerIoniz.rate[j],  cisum
    endfor
endif
;
if rec then begin
    reclvl = (*self.reclvl)
    recsize = size(reclvl.temp)
    nrectemp = recsize[1]
    nrecrate = recsize[2]
    ;  the interpolated values should saved somehow
    upperRecomb = upper->getRecombRate()
    ioniz = self->getIonizRate()
    for j=0,nt-1 do begin
        recsum = 0.
        for ici = 0,nrecrate-1 do begin
            good = where(reclvl.temp[*,ici] gt 0.)
            lvl2 = reclvl.lvl2[ici]
;             print,' lvl2, temp, rate = ', lvl2, reclvl.temp[good,ici], reclvl.rate[good,ici]
            y2 = spl_init(reclvl.temp[good,ici],alog10(reclvl.rate[good,ici]))
            recrate = spl_interp(reclvl.temp[good,ici],alog10(reclvl.rate[good,ici]),y2,alog10(t[j]))
            recsum = recsum + 10.^recrate
;             print,' lvl2, temp, recrate = ', lvl2, t[j], 10.^recrate
            qq[ci+n_levels-1+rec,ci+lvl2-1,j] = qq[ci+n_levels-1+rec,ci+lvl2-1,j] + 10.^recrate     
        endfor
        ;  make sure we aren't double-booking
        qq[ci,ci+n_levels-1+rec,j] = qq[ci,ci+n_levels-1+rec,j] + ioniz.rate[j]
        qq[ci+n_levels-1+rec,ci,j] = qq[ci+n_levels-1+rec,ci,j] + upperRecomb.rate[j] - recsum
;         print, 't, recomb, recsum = ',t[j],upperRecomb.rate[j],  recsum
    endfor
endif

;______________

IF sumtst EQ 0 THEN BEGIN
  FOR j=0,nt-1 DO BEGIN
    qqx=qq[*,*,j]
    FOR i=0,nxne-1 DO BEGIN
      pp=matrix_solver(xne[i],aa,qqx,ppr[*,*,j],aax,c=c)
		;
		;  kpd - need to renormalize
		;
; 		print, ' ci,rec,n_levels, all = ',ci,rec,n_levels,ci+n_levels+rec
; 		help, pp
; 		help,pop
		tot = total(pp[ci:ci+n_levels-1])
		pp = pp/tot
		tot = total(pp[ci:ci+n_levels-1])
      pop[j,i,*]= pp[ci:ci+n_levels-1]
    ENDFOR
  ENDFOR
 ;  kpd :  non-maxwellians does not work properly now
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

popstr = {temperature:t, population:pop, density:xne, ci:ci, rec:rec, n_levels:n_levels}
ptr_pop = ptr_new(popstr)
self.population = ptr_pop

if ptr_valid(self.population) then begin
	print, 'population valid'
endif else print, 'population NOT valid'

END
