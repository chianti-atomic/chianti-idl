
PRO pop_processes, IONNAME, DENS=DENS, TEMP=TEMP, LEVEL=LEVEL, $
                   PATH=PATH, NOPROT=NOPROT, RPHOT=RPHOT, RADTEMP=RADTEMP, $
                   DIEL=DIEL, data_str=data_str, n_levels=n_levels, $
                   radfunc=radfunc, quiet=quiet, pop=pop, verbose=verbose

;+
; NAME
;
;      POP_PROCESSES
;
; PROJECT
;
;      CHIANTI
;
; PURPOSE:
;
;      Outputs to the screen the contributions of the different physical 
;      processes to the population of the specified level within the ion. 
;
;
;      E.g., for Fe XIII, level 4, the output is:
;
;      Population leaving level 4
;        rad. decay:     1.51e+01     39.17%
;        e de-exc:       3.56e-01      0.92%
;        e exc:          2.28e+01     59.12%
;        p de-exc:       2.63e-01      0.68%
;        p exc:          4.05e-02      0.11%
;        stim. emiss:    0.00e+00      0.00%
;        photoexc:       0.00e+00      0.00%
;                        --------
;                 TOTAL  3.85e+01
;      
;      Population entering level 4
;        rad. decay:     3.59e+01     93.19%
;        e de-exc:       3.81e-02      0.10%
;        e exc:          1.46e+00      3.79%
;        p de-exc:       3.18e-03      0.01%
;        p exc:          1.12e+00      2.91%
;        stim. emiss:    0.00e+00      0.00%
;        photoexc:       0.00e+00      0.00%
;                        --------
;                 TOTAL  3.85e+01
;
;      which shows that the level population is dominated by electron
;      excitation and cascading into the level, and by radiative decay
;      out of the level.
;
;      Note that the rates for each physical process are multiplied by the 
;      population of originating level (this results in the totals for 
;      entering and leaving the level to balance).
;
;      For some ions, ionization and recombination are additional processes
;      included when working out the level balance. Because of the way these
;      processes are included (see the v.5 paper for details), the populations
;      entering and leaving some levels *will not balance*.
;
; INPUTS
;
;      IZ      The atomic number of the ion
;
;      ION     The spectroscopic number of the ion (e.g., 12 = XII)
;
; OPTIONAL INPUTS
;
;      LEVEL   The ion level for which information is required.
;
;      DENS    Electron density at which rates calculated (units: cm^-3).
;              If not specified, a value of 10^10 is assumed.
;
;      TEMP    Temperature at which rates calculated (units: K). If not set,
;              then T_max of the ion is used
;
;      PATH    If the ion data-files are not in the CHIANTI directories, 
;              then PATH allows you to choose an alternative location.
;
;      RPHOT   Distance from the centre of the star in stellar radius units.
;              I.e., RPHOT=1 corresponds to the star's surface. (Default is
;              infinity, i.e., no photoexcitation.)
;
;      RADTEMP Specify background radiation temperature (default: 6000 K)
;
;      N_LEVELS Restrict the ion model to this number of levels. E.g., if
;               the CHIANTI model contains 40 levels for the ion, then
;               setting N_LEVELS=12 reduces the model to 12 levels.
;
;      RADFUNC  The name of a user-defined function that will generate
;               a radiation spectrum as a function of temperature. 
;               This radiation field will replace the black-body that
;               is assumed when using the RADTEMP keyword in the call
;               to pop_solver.
;
; KEYWORDS
;
;      NOPROT  If set, then the default setting will be NOT to use 
;              proton rates. This can be changed within the routine.
;
;      QUIET   If set, then do not print information to the
;              screen. (This is useful if you are simply using
;              pop_processes to access the DATA structure.)
;
;      VERBOSE If set, then additional information about the level is
;              given, specifically the dominant channels through which
;              the level is populated.
;
; OPTIONAL OUTPUTS
;
;      POP     A 1D array giving the ion level populations for the
;              specified temperature and density
;
;      DATA_STR This is a structure containing the atomic data that
;               goes into the level population calculation. It is
;               identical to the structure returned by
;               POP_SOLVER. Please see this routine for more details. 
;
; OUTPUTS
;
;      Information about the population processes is written to the
;      IDL window.
;
; CALLS
;
;      R2W, ZION2FILENAME, PROTON_DENS, POP_SOLVER, SETUP_ION
;
; EXAMPLES
;
;      One can compare the effect of cascading on a level population by using
;      the N_LEVELS keyword. Consider the case of Fe XIV:
;        IDL> pop_processes,'fe_14',lev=5
;        IDL> pop_processes,'fe_14',lev=5,n_levels=12
;      With the first call there are two dominant terms to the population
;      entering level 5: approximately 47% for radiative decays (cascading)
;      and 53% for electron excitation. Setting n_levels=12, one finds that
;      the cascading contribution disappears as there are no longer any
;      high-lying levels that cascade into level 5. The cascading provides
;      a strong contribution to the population of this level.
;
; HISTORY
;
;      Ver.1, 11-Sep-2002, Peter Young
;      Ver.2, 15-Jan-2004, Peter Young
;          modified call to pop_solver following recent revision to
;          pop_solver; changed input from IZ,ION to IONNAME to match other
;          CHIANTI routines
;      Ver.3, 26-May-2005, Peter Young
;          changed TEMP and DENS keywords
;      Ver.4, 10-Jun-2005, Peter Young
;          added common block for ionization/recombination data and modified
;          INPUT structure.
;      Ver.5, 14-Jun-2005, Peter Young
;          routine now prints the percentage contribution of each process;
;          added N_LEVELS= keyword
;      Ver.6, 1-Jul-2005, Peter Young
;          added warning for ions with ionization/recombination
;      Ver.7, 12-Jun-2009, Enrico Landi
;          Changed the definition of the temperature array for ion fractions
;          in the IONREC variable, now taken directly from the output of
;          READ_IONEQ.PRO.
;      Ver.8, 9-May-2013, Peter Young
;          Added /QUIET keyword.
;      Ver.9, 29-Jan-2014, Peter Young
;          I've added the keyword /VERBOSE.
;      Ver.10, 21-May-2014, Peter Young
;          Fixed bug when level is the last level in the model.
;-

COMMON elvlc,l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref
COMMON wgfa, wvl,gf,a_value
COMMON upsilon, splstr
COMMON radiative, radt,dilute
COMMON proton, pstr, pe_ratio
COMMON elements,abund,abund_ref,ioneq,temp_all,ioneq_ref
COMMON ionrec,rec_rate,ci_rate,temp_ionrec,luprec,lupci,status

IF N_PARAMS() LT 1 THEN BEGIN
  PRINT,'Use: IDL> pop_processes, ionname, [level= , dens= , temp= , '
  PRINT,'                           radtemp= , rphot= , path= , /noprot, '
  PRINT,'                           /diel, /quiet, /verbose, data_str=, pop='
  print,'                           radfunc=, n_levels= ]'
  RETURN
ENDIF

convertname,ionname,iz,ion
zion2filename,iz,ion,filename,name=name,diel=diel

IF n_elements(radtemp) EQ 0 THEN radt=6000. ELSE radt=radtemp
IF n_elements(rphot) EQ 0 THEN dilute=0d0 ELSE dilute=r2w(rphot)
IF n_elements(temp) EQ 0 THEN BEGIN
  read_ioneq,!ioneq_file,temp_all,ioneq,ioneq_ref
  f_all=ioneq(*,iz-1,ion-1)
  ind=where(f_all EQ max(f_all))
  temp=temp_all(ind) & temp=temp(0)
  temp=10.^temp
  t_ioneq=temp_all
ENDIF
IF n_elements(dens) EQ 0 THEN dens=10.^10
IF n_elements(level) EQ 0 THEN level=1

IF keyword_set(diel) THEN diel = 1 ELSE diel = 0

IF N_ELEMENTS(path) NE 0 THEN filename=concat_dir(path, name) 

setup_ion,name,-1,-1,wvltst,lvl1,lvl2,wvl1,gf1,a_value1,path=path, $
     noprot=noprot

pe_ratio=proton_dens(alog10(temp))

;
; ionization/recombination data
;
IF status lt 0 THEN ionrec = {rec:0., ci:0., temp:0., lev_up_rec:0., $
                              lev_up_ci:0., status:-1, ioneq:0., temp_ioneq:0.}
IF status gt 0 THEN BEGIN
  IF ion gt 1 THEN ioneq_ionrec=reform(ioneq(*,iz-1,ion-2:ion))
  IF ion eq 1 THEN ioneq_ionrec=reform(ioneq(*,iz-1,ion-1:ion)) ; No coll.ionization to neutral ions!
  ionrec = {rec:rec_rate, ci:ci_rate, temp:temp_ionrec, lev_up_rec:luprec, $
            lev_up_ci:lupci, status:status, ioneq:ioneq_ionrec, temp_ioneq:t_ioneq}
ENDIF

input = {gname:ionname, jj:jj, ecm:ecm,ecmth:ecmth, $
         wvl:wvl, a_value:a_value, splstr:splstr, $
         pe_ratio:pe_ratio,prot_struc:pstr, $
         dilute:dilute, radtemp:radt, $
         ionrec: ionrec}

;----------X
; Find the size of the population array
;
pop_solver,input,temp,dens,pop,n_levels=n_levels, data_str=data_str, $
     radfunc=radfunc
pop=reform(pop(0,0,*))
nlvls2=n_elements(pop)
;----------X

levstr=strpad(term[level-1],30,/after)

IF NOT keyword_set(quiet) THEN BEGIN 
print,''
print,format='("Level:   ",a30)',levstr
print,''
print,format='("Log10 Temperature: ",f5.1)',alog10(temp)
print,format='("Log10 Density:     ",f5.1)',alog10(dens)

print,''
print,'Population leaving level '+strtrim(string(level),2)
ENDIF 

sum=0.

;
; The following lines go through each process and extract the rate out of
; the level. The rates are stored in variables _OUT, and the total rate is
; stored in SUM.
;
IF level GT 1 THEN BEGIN
  a_leave=total(data_str.aa[level-1,0:level-2])
  a_out=a_leave*pop[level-1]
  sum=sum+a_leave*pop[level-1]
ENDIF

IF level GT 1 THEN BEGIN
  c_leave=total(data_str.cc[level-1,0:level-2])
  e_deexc_out=c_leave*pop[level-1]
  sum=sum+c_leave*pop[level-1]
ENDIF

IF level LT nlvls2 THEN BEGIN 
  c_leave=total(data_str.cc[level-1,level:*])
  e_exc_out=c_leave*pop[level-1]
  sum=sum+c_leave*pop[level-1]
ENDIF ELSE BEGIN
  e_exc_out=0.
ENDELSE 

IF level GT 1 THEN BEGIN
  p_leave=total(data_str.ccp[level-1,0:level-2])
  p_deexc_out=p_leave*pop[level-1]
  sum=sum+p_leave*pop[level-1]
ENDIF

IF level LT nlvls2 THEN BEGIN 
  p_leave=total(data_str.ccp[level-1,level:*])
  p_exc_out=p_leave*pop[level-1]
  sum=sum+p_leave*pop[level-1]
ENDIF  ELSE BEGIN
  p_exc_out=0.
ENDELSE 


IF level GT 1 THEN BEGIN
  ph_leave=total(data_str.aax[level-1,0:level-2])
  stem_out=ph_leave*pop[level-1]
  sum=sum+ph_leave*pop[level-1]
ENDIF

IF level LT nlvls2 THEN BEGIN 
  ph_leave=total(data_str.aax[level-1,level:*])
  ph_out=ph_leave*pop[level-1]
  sum=sum+ph_leave*pop[level-1]
ENDIF ELSE BEGIN
  ph_out=0.
ENDELSE 


;
; The following lines print the various rates to the screen.
; The ground level does not have any decays, so the if statements filter
; out these cases below.
;
IF NOT keyword_set(quiet) THEN BEGIN 
IF level GT 1 THEN print,format='("  rad. decay:   ",e10.2,f10.2,"%")', $
     a_out,a_out/sum*100.
IF level GT 1 THEN print,format='("  e de-exc:     ",e10.2,f10.2,"%")', $
     e_deexc_out,e_deexc_out/sum*100.
print,format='("  e exc:        ",e10.2,f10.2,"%")',$
     e_exc_out, e_exc_out/sum*100.
IF level GT 1 THEN print,format='("  p de-exc:     ",e10.2,f10.2,"%")', $
     p_deexc_out,p_deexc_out/sum*100.
print,format='("  p exc:        ",e10.2,f10.2,"%")', $
     p_exc_out, p_exc_out/sum*100.
IF level GT 1 THEN print,format='("  stim. emiss:  ",e10.2,f10.2,"%")', $
     stem_out,stem_out/sum*100.
print,format='("  photoexc:     ",e10.2,f10.2,"%")', $
     ph_out,ph_out/sum*100.

print,'                  --------'
print,format='("           TOTAL",e10.2)',sum

print,''
ENDIF


;------------
IF NOT keyword_set(quiet) THEN BEGIN 
print,'Population entering level '+strtrim(string(level),2)
ENDIF 

sum=0.

IF level LT nlvls2 THEN BEGIN 
  a_enter=total(data_str.aa[level:*,level-1]*pop[level:*])
  sum=sum+a_enter
ENDIF ELSE BEGIN
  a_enter=0.
ENDELSE

IF level GT 1 THEN BEGIN
  e_exc_enter=total(data_str.cc[0:level-2,level-1]*pop[0:level-2])
  sum=sum+e_exc_enter
ENDIF

IF level LT nlvls2 THEN BEGIN 
  e_deexc_enter=total(data_str.cc[level:*,level-1]*pop[level:*])
  sum=sum+e_deexc_enter
ENDIF ELSE BEGIN
  e_deexc_enter=0.
ENDELSE 

IF level GT 1 THEN BEGIN
  p_exc_enter=total(data_str.ccp[0:level-2,level-1]*pop[0:level-2])
  sum=sum+p_exc_enter
ENDIF

IF level LT nlvls2 THEN BEGIN 
  p_deexc_enter=total(data_str.ccp[level:*,level-1]*pop[level:*])
  sum=sum+p_deexc_enter
ENDIF ELSE BEGIN
  p_deexc_enter=0.
ENDELSE 

IF level LT nlvls2 THEN BEGIN 
  stem_enter=total(data_str.aax[level:*,level-1]*pop[level:*])
  sum=sum+stem_enter
ENDIF ELSE BEGIN
  stem_enter=0.
ENDELSE 

IF level GT 1 THEN BEGIN
  ph_enter=total(data_str.aax[0:level-2,level-1]*pop[0:level-2])
  sum=sum+ph_enter
ENDIF

;
; It's possible that status=1, but the ioniz/recomb coefficients aren't calculated because the ion fraction of the ion is zero at the input temperature. The check on data_str.correction below accounts for this case
;
IF status GT 0 AND n_elements(data_str.correction) GT 1 THEN BEGIN
  ci_enter=data_str.ion_rate[level-1]*dens*data_str.frac_low
  rec_enter=data_str.rec_rate[level-1]*dens*data_str.frac_high
  sum=sum+ci_enter+rec_enter
ENDIF


;
; The following lines print the various rates to the screen.
; The ground level does not have any decays, so the if statements filter
; out these cases below.
;
IF NOT keyword_set(quiet) THEN BEGIN 
print,format='("  rad. decay:   ",e10.2,f10.2,"%")', $
     a_enter,a_enter/sum*100.
print,format='("  e de-exc:     ",e10.2,f10.2,"%")', $
     e_deexc_enter,e_deexc_enter/sum*100.
IF level GT 1 THEN print,format='("  e exc:        ",e10.2,f10.2,"%")',$
     e_exc_enter, e_exc_enter/sum*100.
print,format='("  p de-exc:     ",e10.2,f10.2,"%")', $
     p_deexc_enter,p_deexc_enter/sum*100.
IF level GT 1 THEN print,format='("  p exc:        ",e10.2,f10.2,"%")', $
     p_exc_enter, p_exc_enter/sum*100.
print,format='("  stim. emiss:  ",e10.2,f10.2,"%")', $
     stem_enter,stem_enter/sum*100.
IF level GT 1 THEN print,format='("  photoexc:     ",e10.2,f10.2,"%")', $
     ph_enter,ph_enter/sum*100.
IF status GT 0 AND n_elements(data_str.correction) GT 1  THEN BEGIN
  print,format='("  ionization:   ",e10.2,f10.2,"%")', $
       ci_enter,ci_enter/sum*100.
  print,format='("  recomb.:      ",e10.2,f10.2,"%")', $
       rec_enter,rec_enter/sum*100.
ENDIF

print,'                  --------'
print,format='("           TOTAL",e10.2)',sum
ENDIF 


;
; In this section I print out the principal channels for populating
; the specified level. At present I've only coded radiative
; decay and electron excitation.
;
IF keyword_set(verbose) THEN BEGIN
  n=n_elements(pop)
  s={lev: 0, levstr: '', rate: 0., type: ''}
  rstr=replicate(s,n*2)
 ;
  rstr[0:n-1].rate=pop*data_str.aa[*,level-1]
  rstr[0:n-1].lev=indgen(n)+1
  rstr[0:n-1].levstr=strpad(term,30,/after)
  rstr[0:n-1].type='cascade'
 ;
  rstr[n:2*n-1].rate=pop*data_str.cc[*,level-1]
  rstr[n:2*n-1].lev=indgen(n)+1
  rstr[n:2*n-1].levstr=strpad(term,30,/after)
  rstr[n:2*n-1].type='e- excitation'
 ;
  k=reverse(sort(rstr.rate))
  tot_rate=total(rstr.rate)
  print,''
  print,'The channels through which the level is populated are given below, organized according'
  print,'to their importance.'
  FOR i=0,9 DO BEGIN
    j=k[i]
    print,format='(i10,2x,a30,f9.1,"%",2x,a15)',rstr[j].lev,rstr[j].levstr, $
          rstr[j].rate/tot_rate*100.,rstr[j].type
  ENDFOR 
  print,'Warning: the above calculation only considers radiative decay and electron excitation.'
  print,'         Other processes have yet to be added.'
ENDIF 

IF status GT 0 AND NOT keyword_set(quiet) THEN BEGIN
  print,''
  print,'** For this ion ionization and recombination processes are included when'
  print,'** working out the level populations. Due to the way these processes are'
  print,'** included (see CHIANTI v.5 paper for details) the populations entering and'
  print,'** leaving the above level may not balance.'
  print,''
ENDIF

END
