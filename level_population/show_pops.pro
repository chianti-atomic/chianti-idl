;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory (USA), the University of Florence (Italy), the
;       University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
;
; NAME: SHOW_POPS
;       
; PURPOSE:
;
;	To display populations of significant levels in a CHIANTI ion 
;	model
;
; CATEGORY:
;
;       Scientific analysis
;
; EXPLANATION:
;
;	Displays percentage level populations and level IDs for all levels 
;	in the specified ion with populations greater than 0.01%. If the 
;	temperature is not specified, then it is taken to be where the 
;	maximum of the ionisation fraction is.
;
;       If the keyword /ALL is set, the output populations are relative 
;       values and not percentages. All levels are shown.
;
; CALLING SEQUENCE:
;
;	SHOW_POPS, IZ, ION [, POPSTR, DENS= , TEMP= , /NOPROT,
;                           RADTEMP= , RPHOT= , /ALL, /DIEL, PATH=,
;                           IONEQ_FILE= , ABUND_FILE= ]
;
; EXAMPLES:
;
;	show_pops,26,13,popstr
;	show_pops,26,13,dens=7.5,temp=6.0,rphot=1.2
;
; INPUTS:
;
;	IZ	The atomic number of the ion
;
;	ION	The spectroscopic number of the ion (e.g., 12 = XII)
;
; OPTIONAL INPUTS:
;
;	DENS	Logarithm of electron density
;
;	TEMP	Logarithm of electron temperature. If not specified, then
;               T_max of the ion is used
;
;	RADTEMP	Specify background radiation temperature (default: 6000 K)
;
;	RPHOT   Distance from the centre of the star in stellar radius units.
;               I.e., RPHOT=1 corresponds to the star's surface. (Default is
;               infinity, i.e., no photoexcitation.)
;
;	PATH	Directly specify the path where the ion data is contained, 
;		e.g., path='/home/other_data/fe_13'
;
;       N_LEVELS  The size of the ion model is automatically determined 
;                 from the information in the CHIANTI data files. Setting 
;                 this keyword allows the number of levels in the model to 
;                 be reduced. E.g., N_LEVELS=14 reduces the model to the 
;                 first 14 levels given in the .ELVLC file.
;
;       IONEQ_FILE To include proton rates in the level balance equations 
;                  requires the number of density of protons to be known, 
;                  and this requires an ion balance file and an abundance 
;                  file to be specified, which are done through the 
;                  IONEQ_FILE and ABUND_FILE keywords. If they are not set 
;                  then the default files specified by !ioneq_file and 
;                  !abund_file are used.
;
;       ABUND_FILE See IONEQ_FILE.
;
;       SUM_MWL_COEFFS  An array of coefficients of the same length as 
;                       the array of temperatures. Electron and proton rate 
;                       coefficients will be calculated at each temperature 
;                       and then a weighted sum of the coefficients is 
;                       performed using SUM_MWL_COEFFS. This allows 
;                       non-Maxwellian energy distributions to be 
;                       incorporated into the level balance equations.
;                       If this keyword is set for an ion that has ionization
;                       and recombination included in the level balance, then
;                       these processes will be switched off for the
;                       calculation since their rates implicitly assume a
;                       single Maxwellian to describe the ion fractions of
;                       the neighbouring ions.
;
;       LEVEL   Allows the control of which level populations are displayed
;               to the screen. If set to a positive scalar or array, then
;               only those levels are printed. If it is set to a negative
;               scalar (-n), then all level populations up to level n are
;               printed. E.g., LEVEL=20 (only level 20); LEVEL=[5,7,20]
;               (levels 5, 7 and 20); LEVEL=-20 (all levels up to level 20).
;
;       RADFUNC         The name of a user-defined function that will generate
;                       a radiation spectrum as a function of temperature. 
;                       This radiation field will replace the black-body that
;                       is assumed when using the RADTEMP keyword in the call
;                       to pop_solver.
;
;; OPTIONAL OUTPUTS
;
;       POPSTR  Send level population information to a structure. POPSTR has 
;               the tags:
;               .dens    Density (cm^-3)
;               .temp    Temperature (K). Can be an array if SUM_MWL_COEFFS=
;                        is used.
;               .radtemp RADTEMP. Set to -1 if RADTEMP not 
;                        set.
;               .rphot   RPHOT value. Set to -1 if RPHOT not set.
;               .proton  String set to 'yes' if proton rates included, 'no' 
;                        otherwise
;               .version CHIANTI version used to derive populations.
;               .date    Date and time at which structure created.
;               .level   Structure containing level information. Tags are:
;                  .index   CHIANTI level index
;                  .term    String containing level identifier.
;                  .pop     Population of level
;               .sumtst  Set to 1 if the SUM_MWL_COEFFS keyword has been used.
;                        Set to 0 otherwise.
;               .sum_mwl_coeffs Contains SUM_MWL_COEFFS. Set to -1 if sumtst=0.
;
; KEYWORDS:
;
;	ALL	Show populations for all levels.
;               the output populations are relative 
;               values and not percentages.
;
;       NOPROT  If set, then the default setting will be NOT to use 
;               proton rates. This can be changed within the routine.
;
;       DIEL    Use the dielectronic recombination files. E.g., for Fe XXII, 
;               the routine will read in the fe_22d.* files instead of the 
;               fe_22.* files.
;
; COMMON BLOCKS:
;
;	ELVLC, WGFA, UPSILON, RADIATIVE, PROTON, ELEMENTS
;
; CALLS:
;
;	ZION2FILENAME, POP_SOLVER, SETUP_ION, READ_IONEQ,
;       READ_ABUND, PROTON_DENS, CONCAT_DIR
;
; HISTORY:
;
;	Ver 1, PRY 22-Sep-97
;	Ver.2, PRY 5-Sep-98  - added call to choose_ioneq
;	Ver.3, PRY 23-Apr-99 - calls pop_solver now; added DENS keyword
;	Ver.4, PRY 18-Dec-99 - added deu to upsilon common block to be 
;			consistent with main Chianti routines.
;       Ver.5, PRY 7-Aug-00  - added /DIEL keyword to allow populations of 
;                       the dielectronic recombination files to be studied.
;                       Also changed elvlc common block to match new 
;                       version of pop_solver.
;       Ver 6, PRY 10-Oct-00 - now calls setup_ion to read ion data
;       Ver 7, PRY 12-Nov-01 - modified for proton rates, photoexcitation,
;                       and 9 pt splines.
;       Ver 8, PRY 9-Dec-01  - completed modifications for v4 of CHIANTI.
;
;       V.  9, 25-May-2002, Giulio Del Zanna (GDZ) 
;                   generalized directory concatenation to work for
;                   Unix, Windows  and VMS. 
;                   Now we only call zion2filename
;
;       V. 10, 9-Aug-2002, Peter Young
;                   corrected !ioneq_file problem
;
;       V. 11, 12-Aug-2002, Peter Young
;                   added POPSTR output, and tidied up header.
;
;       V. 12, 4-Oct-2003, GDZ
;                  modified the input to POP_SOLVER, now it is dealt with an
;                  input structure.
;
;       V 13, 4-May-2005, Enrico Landi (EL)
;                  Modified in order to include ionization and recombination
;                  data in the input to POP_SOLVER
;
;       V.14, 26-May-2005, Peter Young (implemented by GDZ)
;                  added SUM_MWL_COEFFS optional input for allowing
;                  non-Maxwellian distributions to be considered.
;
;                  added LEVEL= optional input to only print the populations of
;                  a few levels.
;
;       V.15, 5-Jul-2005, Peter Young
;                  added RADFUNC= and /QUIET keywords
;
;       V.16, 14-Sept-2005 Giulio Del Zanna 
;              modified header.
;
;       v.17, 12-Dec-2005, Peter Young
;              routine now performs check on existence of .splups file
;              instead of .elvlc since with v.5 some ions have a .elvlc
;              file but nothing else (fe_3, fe_4)
;
;       v.18, 12-Jun-2009, Enrico Landi
;              Changed the definition of the temperature array for ion fractions
;              in the IONREC variable, now taken directly from the output of
;              READ_IONEQ.PRO
;
; VERSION     :   18, 12-Jun-2009
;-

PRO  SHOW_POPS, IZ, ION, popstr, DENS=DENS, TEMP=TEMP, RPHOT=RPHOT, $
                ALL=ALL, PATH=PATH, NOPROT=NOPROT, $
                N_LEVELS=N_LEVELS, RADTEMP=RADTEMP, DIEL=DIEL, $
                IONEQ_FILE=IONEQ_FILE,ABUND_FILE=ABUND_FILE, $
                SUM_MWL_COEFFS=SUM_MWL_COEFFS, RADFUNC=RADFUNC, $
                LEVEL=LEVEL, QUIET=QUIET


COMMON elvlc,l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref
COMMON wgfa, wvl,gf,a_value
COMMON upsilon, splstr
COMMON radiative, radt,dilute
COMMON proton, pstr, pe_ratio
COMMON elements,abund,abund_ref,ioneq,temp_all,ioneq_ref
COMMON ionrec,rec_rate,ci_rate,temp_ionrec,luprec,lupci,status

IF N_PARAMS() LT 1 THEN BEGIN
  PRINT,'Use: IDL> show_pops, iz, ion, [popstr, dens= , temp= , /all, '
  PRINT,'                           radtemp= , rphot= , path= , /noprot, '
  PRINT,'                           n_levels= , /diel, ioneq_file= , '
  PRINT,'                           abund_file= , level= ]'
  RETURN
ENDIF


IF keyword_set(diel) THEN diel = 1 ELSE diel = 0

;
; The following extracts the names of the files to be read. I want to allow 
; a different path to be chosen, and I extract only the information I need 
; from filename.
;

zion2filename,iz,ion,filename,name=name,diel=diel
gname = name

IF N_ELEMENTS(path) NE 0 THEN filename=concat_dir(path, name)

chckfile=filename+'.splups'
chck=file_exist(chckfile)
IF chck NE 1 THEN BEGIN
  IF NOT keyword_set(quiet) THEN print,'% SHOW_POPS: no data exists for this ion in CHIANTI. Returning...'
  popstr=0
  return
ENDIF


setup_ion,name,-1,-1,wvltst,lvl1,lvl2,wvl1,gf1,a_value1,path=path, $
     noprot=noprot

IF N_ELEMENTS(temp) EQ 0 THEN BEGIN
;------------------------------------+
; Work out the T_max for the ion using the ion balance data in CHIANTI. 
;

dir=concat_dir(!xuvtop,'ioneq')
IF n_elements(ioneq_file) NE 0 THEN ioneq_name=ioneq_file $
     ELSE ioneq_name=!ioneq_file

read_ioneq,ioneq_name,temp_all,ioneq,ioneq_ref
;
f_all=ioneq(*,iz-1,ion-1)
ind=where(f_all EQ max(f_all))
temp=temp_all(ind) 
temp=temp(0)
;------------------------------------+
ENDIF

IF n_elements(abund_file) EQ 0 THEN abund_file=!abund_file
  read_abund,abund_file,abund,abund_ref

pe_ratio=proton_dens(temp)

IF n_elements(rphot) EQ 0 THEN dilute=0d0 ELSE dilute=r2w(rphot)
IF n_elements(radtemp) EQ 0 THEN radt=6d3 ELSE radt=double(radtemp)
IF N_ELEMENTS(dens) EQ 0 THEN dens=10d0 ELSE dens=FLOAT(dens)

IF n_elements(sum_mwl_coeffs) NE 0 THEN BEGIN
  smc=double(sum_mwl_coeffs)
  IF total(smc) NE 1. THEN BEGIN
    print,'%SHOW_POPS: normalizing SUM_MWL_COEFFS'
    smc=smc/total(smc)
  ENDIF
  sumtst=1
ENDIF ELSE BEGIN
  IF n_elements(temp) NE 1 THEN BEGIN
    print,'%SHOW_POP: only 1 temperature should be input if SUM_MWL_COEFFS not specified'
    temp=temp[0]
    pe_ratio=pe_ratio[0]
  ENDIF
  sumtst=0
ENDELSE


;
; Ionization and recombination need the ion fractions of the neighbouring
; ions in order to be included correctly in the level balance equations.
; Equilibrium ion fractions are not compatible with electron distributions
; described by multiple Maxwellians, so the lines below switch ionization
; and recombination when SUM_MWL_COEFFS is set.
;
IF sumtst EQ 1 AND status GT 0 THEN BEGIN
  status=-1   ; this switches off ioniz/recomb
  print,''
  print,'%SHOW_POPS: ionization and recombination will be switched off for this ion'
  print,'            since their implementation in CHIANTI is incompatible with the'
  print,'            SUM_MWL_COEFFS keyword'
ENDIF


;-------------------------------------


IF status lt 0 THEN ionrec = {rec:0., ci:0., temp:0., lev_up_rec:0., lev_up_ci:0., status:-1, ioneq:0., temp_ioneq:0.}
IF status GT 0 THEN BEGIN
    IF ion gt 1 THEN ioneq_ionrec=reform(ioneq(*,iz-1,ion-2:ion))
    IF ion eq 1 THEN ioneq_ionrec=reform(ioneq(*,iz-1,ion-1:ion))  ; No coll.ionization to neutral ions!
    t_ioneq=temp_all                                       ; Defines temperature array for ion abundances for ionrec variable
    ionrec = {rec:rec_rate, ci:ci_rate, temp:temp_ionrec, lev_up_rec:luprec, $
    lev_up_ci:lupci, status:status, ioneq:ioneq_ionrec, temp_ioneq:t_ioneq}
ENDIF

;----------X
;GDZ

input = {gname:gname, jj:jj, ecm:ecm,ecmth:ecmth, $
 wvl:wvl, a_value:a_value, splstr:splstr, $
 pe_ratio:pe_ratio,prot_struc:pstr, dilute:dilute, radtemp:radt, ionrec:ionrec}


pop_solver,input, 10.^temp,10.^dens,pop,n_levels=n_levels, $
     sum_mwl_coeffs=smc, radfunc=radfunc

pop=pop(0,0,*)
nlvls2=n_elements(pop)
;----------X

IF NOT keyword_set(quiet) THEN PRINT,''
IF NOT keyword_set(quiet) THEN PRINT,FORMAT='("Log10 density:     ",f5.1)',dens
IF n_elements(sum_mwl_coeffs) GT 0 THEN BEGIN
  print,'Using a sum of Maxwellians'
ENDIF ELSE BEGIN
  IF NOT keyword_set(quiet) THEN PRINT,FORMAT='("Log10 temperature: ",f5.1)',temp
ENDELSE
IF NOT keyword_set(quiet) THEN PRINT,''

str={index: 0, term: '', pop: 0d0}

lev=replicate(str,nlvls2)
lev.index=findgen(nlvls2)+1
lev.term=term[0:nlvls2-1]
lev.pop=reform(pop)
IF sumtst EQ 0 THEN smc=-1.
popstr={dens: 10.^dens, temp: 10.^temp, level: lev, radtemp: float(radt), $
        rphot: 0., proton: 'yes', version: '', date: systime(), $
       sum_mwl: sumtst, sum_mwl_coeffs: smc}

IF NOT keyword_set(rphot) THEN BEGIN
  popstr.rphot=-1
  popstr.radtemp=-1
ENDIF
;
fname=concat_dir(!xuvtop,'VERSION')
str1=''
IF file_exist(fname) THEN BEGIN
  openr,lun,fname,/get_lun
  readf,lun,str1
  free_lun,lun
  str1='CHIANTI '+str1
ENDIF ELSE BEGIN
  str1='Please update your version of CHIANTI to v.4 or later'
ENDELSE
popstr.version=str1
IF keyword_set(noprot) THEN popstr.proton='no'

IF keyword_set(quiet) THEN return

IF n_elements(level) NE 0 THEN BEGIN
  IF level[0] LT 0 THEN BEGIN
    FOR i=0,-level[0]-1 DO PRINT,FORMAT='(i3,a23,e10.2)',i+1,term(i),pop(i)
  ENDIF ELSE BEGIN
    n=n_elements(level)
    FOR i=0,n-1 DO BEGIN
      IF level[i] GT 0 THEN BEGIN
        print,format='(i3,a23,e10.2)',level[i],term(level[i]-1),pop(level[i]-1)
      ENDIF
    ENDFOR
  ENDELSE
  return
ENDIF

IF KEYWORD_SET(all) THEN BEGIN
  FOR i=0,nlvls2-1 DO $
       PRINT,FORMAT='(i3,a23,e10.2)',i+1,term(i),pop(i)
ENDIF ELSE BEGIN
  FOR i=0,nlvls2-1 DO BEGIN
    IF pop(i) GT 0.0001 THEN PRINT,FORMAT='(i3,a23,f9.2)', $
                                     i+1,term(i),100.*pop(i)
  ENDFOR
ENDELSE

END
