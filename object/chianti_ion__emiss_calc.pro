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
; NAME: CHIANTI_ION::EMISS_CALC
;       
; PURPOSE:
;
;       To compute the emissivities of all lines of a specified ion over
;	given ranges of temperature and density.
;
; CATEGORY:
;
;       Scientific analysis
;
; EXPLANATION:
;
; 	This routine calculates:
;
;		hc
;       	--  * N_j  * A_ji
;      	       lamb
;
; 	where hc = 1.986 * 10^-8 erg AA, lamb is in angstroms, N_j is the 
;	fraction of ions in the upper emitting level j, and A_ji is the 
;	radiative decay rate for the transition.
;
;	The emissivities are stored in a structure called EMISS that also
;	holds the wavelength of the transition, the level numbers i and j 
;	and also a 'flag', which is set to -1 if the wavelength is negative.
;
;	The temperature and density ranges can be specified directly using
;	the TEMP and DENS keywords. Setting TMAX to the log T_max of the 
;	ion, gives emissivities for 3 temperatures: log T_max +- 0.15. 
;	If DENS is not set, then it is set to 8 to 12 in 0.5 dex intervals. 
;	STDENS allows the start density (of 8) to be changed to some other 
;	value; ND allows the number of densities to be varied (default 9); 
;	DINT allows the density interval to be varied (default 0.5).
;
; CALLING SEQUENCE:
;
;       EMISS=EMISS_CALC (IZ, ION, [ TEMP=TEMP, DENS=DENS, RADT=RADT, $
;				DIL=DIL, PATH=PATH, /NO_DE, /PROTON, $
;				QUIET, PRESSURE=PRESSURE)
;
; EXAMPLES:
;
;	EMISS=EMISS_CALC(26,13)
;	EMISS=EMISS_CALC(26,13,temp=[6.2],dens=findgen(5)+8)
;	EMISS=EMISS_CALC(26,13,temp=findgen(11)/100.+5.5,press=10.^15)
;
; INPUTS:
;
;	IZ	The atomic number of the ion
;
;	ION	The spectroscopic number of the ion (e.g., 12 = XII)
;
; OPTIONAL INPUTS:
;
;	TEMP	Direct specification of the temperature range (log T)
;
;	DENS	Direct specification of the density range (log Ne)
;
;	RADTEMP	Specify background radiation temperature (default: 6000 K)
;
;	RPHOT   Distance from the centre of the star in stellar radius units.
;               I.e., RPHOT=1 corresponds to the star's surface. (Default is
;               infinity, i.e., no photoexcitation.)
;
;	PATH	If specified, the routine will look for the atomic data in 
;		the PATH directory, rather than in the CHIANTI database
;
;	PRESSURE	If a temperature array is given, and PRESSURE set, 
;			then the emissivities will be evaluated at the 
;			specified temperatures, but for densities = 
;			pressure/temperature. If DENS is set, then it will 
;			be ignored. The pressure is assumed to be in units 
;			K * cm^-3.
;
;       ABUND_FILE  The name of a CHIANTI abundance file. This is used for 
;               calculating the proton to electron ratio. Default is 
;               !abund_file.
;
;       IONEQ_FILE  The name of a CHIANTI ion balance file. This is used for 
;               calculating the proton to electron ratio and evaluating 
;               the T_max of the ion. Default is !ioneq_file.
;
;       SUM_MWL_COEFFS  An array of coefficients of the same length as 
;                       the array of temperatures. Electron and proton rate 
;                       coefficients will be calculated at each temperature 
;                       and then a weighted sum of the coefficients is 
;                       performed using SUM_MWL_COEFFS. This allows 
;                       non-Maxwellian energy distributions to be 
;                       incorporated into the level balance equations.
;                       This keyword is not compatible with the PRESSURE
;                       keyword.
;
;       RADFUNC         The name of a user-defined function that will generate
;                       a radiation spectrum as a function of temperature. 
;                       This radiation field will replace the black-body that
;                       is assumed when using the RADTEMP keyword in the call
;                       to pop_solver.
;
; KEYWORDS:
;
;	NO_DE	Drops the hc/lambda factor in the computation of the 
;		emissivities. Useful for emission measure analyses involving 
;		photon fluxes
;
;       NOPROT  If set, then the default setting will be NOT to use 
;               proton rates. This can be changed within the routine.
;
;	QUIET	If set, don't list the temperatures and densities at which 
;		the emissivities are caculated.
;
;       DIEL   If the dielectronic recombination files exist for the ion, 
;               then these are used to derive the emissivities.
;
;       NO_SETUP  If emiss_calc is called from a routine where the ELEMENTS
;                 common block has already been set up, then this keyword
;                 stops emiss_calc loading up the common block
;
; OUTPUT:
;
;	The structure that is output has the following tags:
;
;	.ion_name	string; contains ion name, e.g., 'Fe XIII'
;	.lambda		float; contains wavelength
;	.level1		integer; contains lower level of transition
;       .lvl1_desc      string; gives config. and term of lower level
;	.level2		integer; contains upper level of transition
;       .lvl2_desc      string; gives config. and term of upper level
;	.flag		integer; a flag to mark particular transitions
;	.em		fltarr(nt,nd); contains emissivities at particular 
;			  temperatures and densities.
;
; PROGRAMMING NOTES:
;
;	Transitions where only theoretical energies are available for at 
;	least one of the two levels are assigned negative wavelengths in 
;	the .wgfa file. With emiss_calc, the wavelength is set to be 
;	positive, but emiss.flag is set to -1 for that transition, so 
;	that other routines can keep track of the theoretical wavelengths.
;
; COMMON BLOCKS:  None
;
;
; CALLS:
;
;	SETUP_ION,
;       ZION2FILENAME, ZION2SPECTROSCOPIC
;
; HISTORY:
;
;
;		V.1, 5-jan-2011, Ken Dere
;			first release - provides the emiss_calc method for the chianti_ion object
;			derived from emiss_calc.pro V.30
;
; VERSION     :   1, 5-jan-2011
;
;-

Pro  CHIANTI_ION::EMISS_CALC, TEMP=TEMP, DENS=DENS, RADTEMP=RADTEMP, $
                      RPHOT=RPHOT, $
                      PATH=PATH, NO_DE=NO_DE, NOPROT=NOPROT, $
                      QUIET=QUIET, PRESSURE=PRESSURE, DIEL=DIEL, $
                      IONEQ_FILE=IONEQ_FILE, ABUND_FILE=ABUND_FILE, $
                      NO_SETUP=NO_SETUP, SUM_MWL_COEFFS=SUM_MWL_COEFFS, $
                      RADFUNC=RADFUNC



; COMMON elvlc,l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref
; COMMON wgfa, wvl,gf,a_value
; COMMON upsilon,splstr
; COMMON radiative, radt,dilute
; COMMON proton, pstr, pe_ratio
; COMMON elements,abund,abund_ref,ioneq,temp_all,ioneq_ref
; COMMON ionrec,rec_rate,ci_rate,temp_ionrec,luprec,lupci,status

; IF N_PARAMS() LT 2 THEN BEGIN
;   PRINT,'Use:  IDL> emiss=emiss_calc(iz,ion [, temp=temp, dens=dens, '
;   PRINT,'                              radtemp=radtemp, rphot=rphot, /no_de,'
;   PRINT,'                              path=path, /noprot, /quiet, $ '
;   PRINT,'                              pressure=pressure, ioneq_file=ioneq_file, $'
;   print,'                              abund_file=abund_file ] )'
;   RETURN
; ENDIF

IF n_elements(pressure) NE 0 AND n_elements(sum_mwl_coeffs) NE 0 THEN BEGIN
  print,'% EMISS_CALC: the PRESSURE and SUM_MWL_COEFFS keywords are incompatible with'
  print,'              each other. Please use one or the other.'
  return
ENDIF

;
; The following extracts the names of the files to be read. I want to allow 
; a different path to be chosen, and I extract only the information I need 
; from filename.
;
; IF keyword_set(diel) THEN diel = 1 ELSE diel = 0
; 
; gname = self.Ions
; ion2filename,gname,filename
; convertname,gname,iz,ion
iz = self.z
ion = self.stage

; setup_ion,name,-1,-1,wvltst,lvl1,lvl2,wvl1,gf1,a_value1,path=path, $
;      noprot=noprot

; IF NOT keyword_set(no_setup) OR n_elements(abund) EQ 0 OR $
;      n_elements(ioneq) EQ 0 THEN BEGIN
;   IF n_elements(ioneq_file) EQ 0 THEN ioneq_file=!ioneq_file
;   read_ioneq,ioneq_file,temp_all,ioneq,ioneq_ref
;   IF n_elements(abund_file) EQ 0 THEN abund_file=!abund_file
;   read_abund,abund_file,abund,abund_ref
; ENDIF

t_ioneq = (*self.ioneq).temperature 
temp_all = t_ioneq
; Defines the temperature array of ionization equilibrium

ioneq = (*self.ioneq).ioneq

;------------------------------------+
; If TEMP not specified, then set TEMP=[log(T_max),log(T_max)+-0.15]
; using the ion balance data in CHIANTI. 
;
;  do_pop determines whether we actually  need to calculate the populations

do_pop = 0

IF N_ELEMENTS(temp) EQ 0 THEN BEGIN
	if ptr_valid(self.population) then begin
		print, ' using population values'
		;  will use the temperature used to calculate the populations
		temp = (*self.population).temperature
		do_pop = 0
	endif else begin
		f_all=ioneq(*,iz-1,ion-1)
		ind=where(f_all EQ max(f_all))
		tmax=temp_all(ind) & tmax=tmax(0)
		temp=[tmax-0.15,tmax,tmax+0.15]
		do_pop = 1
	endelse
ENDIF
nt=N_ELEMENTS(temp)
;------------------------------------+

;;-----------------<>
; Default density range is 8 to 12 in .5 intervals:
;
IF n_elements(dens) EQ 0 THEN begin
	if ptr_valid(self.population) then begin
		;  will use the temperature used to calculate the populations
		dens = (*self.population).density
	endif else begin
		dens=findgen(9)/2. + 8.
		do_pop = 1
	endelse
endif
nd=n_elements(dens)
;
;-----------------<>

IF N_ELEMENTS(pressure) NE 0 THEN BEGIN
  nd=1
  dens=ALOG10(pressure/10.^temp)
  do_pop = 1
ENDIF

;---------------X
sumtst=0
IF n_elements(sum_mwl_coeffs) NE 0 THEN BEGIN
  IF n_elements(temp) NE n_elements(sum_mwl_coeffs) THEN BEGIN
    print,'% EMISS_CALC: SUM_MWL_COEFFS must have the same size as TEMP. Returning...'
    return
  ENDIF
 ;
 ; normalize sum_mwl_coeffs
 ;
  smc=double(sum_mwl_coeffs)
  IF total(smc) NE 1. THEN BEGIN
    print,'% EMISS_CALC: normalizing SUM_MWL_COEFFS'
    smc=smc/total(smc)
  ENDIF
  sumtst=1
ENDIF
;---------------X


;---------------O
zion2spectroscopic,iz,ion,spectname

; IF status GT 0 AND n_elements(sum_mwl_coeffs) NE 0 THEN BEGIN
;   print,'% EMISS_CALC: this ion ('+spectname+') has contributions from ionization'
;   print,'              and recombination that can not be included when the SUM_MWL_COEFFS'
;   print,'              keyword is set. The ionization and recombination processes'
;   print,'              will be switched off.'
;   status=-1   ; this switches off ioniz/rec
; ENDIF
;---------------O


IF n_elements(rphot) EQ 0 THEN dilute=0d0 ELSE dilute=r2w(rphot)
IF n_elements(radtemp) EQ 0 THEN radt=6d3 ELSE radt=double(radtemp)






;
; Create the structure with the ion/rec data, including ion fractions for
; the ion and the two adjacent ones
;

; IF status lt 0 THEN ionrec = {rec:0., ci:0., temp:0., lev_up_rec:0., $
;                               lev_up_ci:0., status:-1, ioneq:0., temp_ioneq:0.}
; IF status gt 0 THEN begin 
; 
;   IF ion gt 1 THEN ioneq_ionrec=reform(ioneq(*,iz-1,ion-2:ion))
;   IF ion eq 1 THEN ioneq_ionrec=reform(ioneq(*,iz-1,ion-1:ion)) ; No coll.ionization to neutral ions!
; 
; 
; ionrec = {rec:rec_rate, ci:ci_rate, temp:temp_ionrec, $
;                               lev_up_rec:luprec, lev_up_ci:lupci, $
;                               status:status, ioneq:ioneq_ionrec, $
;                               temp_ioneq:t_ioneq}
; endif 

;----------X
; Find the size of the population array
;
; pe_ratio=proton_dens(temp[0])
; 
; input = {gname:gname, jj:jj, ecm:ecm,ecmth:ecmth, $
;  wvl:wvl, a_value:a_value, splstr:splstr, $
;  pe_ratio:pe_ratio,prot_struc:pstr, dilute:dilute, radtemp:radt, ionrec:ionrec}

term = (*self.elvlc).data.term
jj = (*self.elvlc).data.j
ecm = (*self.elvlc).data.obs_energy
ecmth = (*self.elvlc).data.theory_energy

lvl1a = (*self.wgfa).lvl1
lvl2a = (*self.wgfa).lvl2
wvl1a = (*self.wgfa).wvl
gf1a = (*self.wgfa).gf
a_value1a = (*self.wgfa).a_value


ntrans = n_elements(lvl1a)
; for ii = 0,ntrans-1 do begin
; 	print,ii,lvl1a[ii],lvl2a[ii],wvl1a[lvl1a[ii]-1,lvl2a[ii]-1],a_value1a[lvl1a[ii]-1,lvl2a[ii]-1]
; endfor
;
lvl1 = intarr(ntrans)
lvl2 = intarr(ntrans)
wvl1 = fltarr(ntrans)
gf1 = fltarr(ntrans)
a_value1 = fltarr(ntrans)
;
for itr = 0,ntrans-1 do begin
	l1 = lvl1a[itr]
	lvl1[itr] = l1
	l2 = lvl2a[itr]
	lvl2[itr] = l2
	wvl1[itr] = wvl1a[l1-1,l2-1]
	gf1[itr] = gf1a[l1-1,l2-1]
	a_value1[itr] = a_value1a[l1-1,l2-1]
endfor

; thision->populate,10.^temp[0],10.^9,pop,radtemp=radtemp,dilute=dilute
if do_pop eq 1 then begin
	print, ' calculating pop'
	self->populate, 10.^temp[0],10.^9
endif else print, ' using existing population'

pop1=self->getPopulation()


nlvls2 = pop1.n_levels
;----------X

;
; remove lines with zero wavelengths
;
ind=where(wvl1 NE 0.)
IF ind[0] NE -1 THEN BEGIN
  lvl1=lvl1[ind]
  lvl2=lvl2[ind]
  wvl1=wvl1[ind]
  gf1=gf1[ind]
  a_value1=a_value1[ind]
ENDIF

;
; remove lines involving levels not in pop_solver model
;
ind=where((lvl1 LE nlvls2) AND (lvl2 LE nlvls2))
IF ind[0] NE -1 THEN BEGIN
  lvl1=lvl1[ind]
  lvl2=lvl2[ind]
  wvl1=wvl1[ind]
  gf1=gf1[ind]
  a_value1=a_value1[ind]
ENDIF

ntrans2=n_elements(wvl1)

;----------------------------------o
; Create the EMISS structure, using ntrans2 to set size
;
fname=concat_dir(!xuvtop,'VERSION')
str1=' '
IF file_exist(fname) THEN BEGIN
  openr,lun,fname,/get_lun
  readf,lun,str1
  free_lun,lun
;  str1='CHIANTI '+str1
ENDIF ELSE BEGIN
  str1=' '
print, 'Please update your version of CHIANTI to v.4 or later'
ENDELSE
;
zion2spectroscopic,iz,ion,species, diel=diel
;
IF keyword_set(sumtst) THEN em=dblarr(1,nd) ELSE em=dblarr(nt,nd)
;
str={lambda: 0., $
	level1: 0, $
        lvl1_desc: '', $
	level2: 0, $
        lvl2_desc: '', $
	avalue: 0., $
	flag:   0, $
	em: em $
        }
;
emiss=replicate(str,ntrans2)
;----------------------------------o

emiss.lambda=abs(wvl1)
emiss.level1=lvl1
emiss.level2=lvl2
emiss.lvl1_desc=term[lvl1-1]
emiss.lvl2_desc=term[lvl2-1]
emiss.avalue = a_value1
ind=where(wvl1 LT 0.)

IF ind[0] NE -1 THEN emiss[ind].flag=-1


IF NOT KEYWORD_SET(quiet) THEN BEGIN
  PRINT,''
  PRINT,'Log_10 temperatures...'
  PRINT,FORMAT='("   ",9f6.2)',alog10(temp)
  PRINT,''
  PRINT,'Log_10 densities...'
  PRINT,FORMAT='("   ",9f6.2)',alog10(dens)
ENDIF

IF NOT keyword_set(no_de) THEN begin
	mult_factor=(1.986d-8/emiss.lambda)*a_value1
endif ELSE mult_factor=a_value1
;
IF N_ELEMENTS(pressure) NE 0 THEN BEGIN
	print, ' number of pressure = ',N_ELEMENTS(pressure)
;   FOR j=0,nt-1 DO BEGIN
;     pe_ratio=proton_dens(temp[j])
;    ;
;     input = {gname:gname, jj:jj, ecm:ecm,ecmth:ecmth, $
;              wvl:wvl, a_value:a_value, splstr:splstr, $
;              pe_ratio:pe_ratio,prot_struc:pstr, dilute:dilute, $
;              radtemp:radt, ionrec:ionrec}
   ;
;   ENDFOR
    self->populate, 10.^temp,pressure/10.^temp, radfunc=radfunc
    pop=self->getPopulation()
   ;
    emiss.em=mult_factor*reform(pop.population[emiss.level2-1])
ENDIF ELSE BEGIN
;   pe_ratio=proton_dens(temp)
 ;
;   input = {gname:gname, jj:jj, ecm:ecm,ecmth:ecmth, $
;            wvl:wvl, a_value:a_value, splstr:splstr, $
;            pe_ratio:pe_ratio,prot_struc:pstr, dilute:dilute, $
;            radtemp:radt, ionrec:ionrec}
 ;
	if do_pop eq 1 then begin
		self->populate, 10.^temp,10.^dens, radfunc=radfunc
	endif
    pop=self->getPopulation()
 ;
  IF sumtst EQ 1 THEN BEGIN
   ;
    IF nd NE 1 THEN emiss.em[0,*]=pop.population[*,*,emiss.level2-1]
    IF nd EQ 1 THEN emiss.em[0,*]=reform(pop.population[*,*,emiss.level2-1],1, $
                                         n_elements(emiss.level2))
   ;
  ENDIF ELSE BEGIN
   ;
    IF nt EQ 1 THEN BEGIN
		IF nd NE 1 THEN emiss.em=pop.population[*,*,emiss.level2-1] ELSE $
			emiss.em=reform(pop.population[*,*,emiss.level2-1],1,n_elements(emiss.level2))
		ENDIF ELSE begin
			emiss.em=reform(pop.population[*,*,emiss.level2-1])
		endelse
   ;
	ENDELSE
 ;
  FOR i=0L,ntrans2-1 DO emiss[i].em=emiss[i].em*mult_factor[i]
ENDELSE

ind=sort(emiss.lambda)
emiss=emiss[ind]
topstr = {density:dens, temperature:temp, ion_name:species, version:str1, emstr:emiss}
self.emiss = ptr_new(topstr)

END
