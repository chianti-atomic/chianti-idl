
PRO pop_processes, IONNAME, DENS=DENS, TEMP=TEMP, LEVEL=LEVEL, $
                   PATH=PATH, NOPROT=NOPROT, RPHOT=RPHOT, RADTEMP=RADTEMP, $
                   DIEL=DIEL, data_str=data_str, n_levels=n_levels, $
                   radfunc=radfunc, quiet=quiet, pop=pop, verbose=verbose, $
                   abund_file=abund_file, ioneq_file=ioneq_file, $
                   noionrec=noionrec, output=output

;+
; NAME:
;      POP_PROCESSES
;
; CATEGORY:
;      CHIANTI; level populations.
;
; PURPOSE:
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
; CALLING SEQUENCE:
;      IDL> POP_PROCESSES, ION_NAME
; 
; INPUTS:
;      Ion_Name:  The name of the ioN in CHIANTI format, e.g., 'fe_13'.
;
; OPTIONAL INPUTS:
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
;      Abund_File: The name of an element abundance file. If not
;               specified, then !abund_file will be used.
;
;      Ioneq_File: The name of an ionization equilibrium file. If not
;               specified, then !ioneq_file will be used.
;
; KEYWORD PARAMETERS:
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
;      NOIONREC If set, then level-resolved ionization and
;              recombination will not be included in the population
;              calculation. 
;
;
; OPTIONAL OUTPUTS:
;      POP     A 1D array giving the ion level populations for the
;              specified temperature and density
;
;      DATA_STR This is a structure containing the atomic data that
;               goes into the level population calculation. It is
;               identical to the structure returned by
;               POP_SOLVER. Please see this routine for more details.
;
;      Output: This contains the data that is printed to the
;              screen. It is a structure with the tags:
;
;              .ionname   Ion name (CHIANTI format).
;              .lvl       Level index
;              .lvlname   The description of the level.
;              .density   Electron number density (cm^-3).
;              .temperature  Temperature (K).
;              .pop       Population of the level.
;              .out       Structure containing the rates out of the
;                         level for each process. 
;              .in        Structure containing the rates into the
;                         level for each process.
;
;              The tags for out and in can be matched up with the
;              processes printed to the screen.
;
; OUTPUTS:
;      Information about the population processes is written to the
;      IDL window. Optionally they can be output to IDL to the
;      OUTPUT structure (see above).
;
; CALLS:
;      ZION2FILENAME, POP_SOLVER, CH_SETUP_ION
;
; EXAMPLES:
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
; MODIFICATION HISTORY:
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
;      Ver.11, 11-Aug-2016, Peter Young
;          Fixed bug: t_ioneq wasn't defined if temp was input;
;          only affected ions with ionrec data.
;      Ver.12, 29-Jan-2018, Peter Young
;          Added call to ch_setup_ion and removed common blocks;
;          tidied up header.
;      Ver.13, 20-Feb-2019, Peter Young
;          Added /NOIONREC keyword.
;      Ver.14, 5-Mar-2019, Peter Young
;          Modified to work correctly with the new CHIANTI 9 models.
;      Ver.15, 2-May-2019, Peter Young
;          Added the optional output OUTPUT.
;      Ver.16, 05-Nov-2020, Peter Young
;          For 2-ion models, normalize populations to only the
;          recombined ion populations (don't include
;          recombining ion populations). 
;-


IF N_PARAMS() LT 1 THEN BEGIN
  PRINT,'Use: IDL> pop_processes, ionname, [level= , dens= , temp= , '
  PRINT,'                           radtemp= , rphot= , path= , /noprot, '
  PRINT,'                           /diel, /quiet, /verbose, data_str=, pop='
  print,'                           radfunc=, n_levels=, /noionrec, output= ]'
  RETURN
ENDIF

convertname,ionname,iz,ion
zion2filename,iz,ion,filename,name=name,diel=diel

IF n_elements(temp) EQ 0 THEN temp=ch_tmax(name,ioneqname=!ioneq_file)
IF n_elements(dens) EQ 0 THEN dens=10.^10
IF n_elements(level) EQ 0 THEN level=1

IF keyword_set(diel) THEN diel = 1 ELSE diel = 0

IF N_ELEMENTS(path) NE 0 THEN filename=concat_dir(path, name) 

;
; Load the atomic data into INPUT.
;
input=ch_setup_ion(name,rphot=rphot,radtemp=radtemp,noprot=noprot, $
                   ioneq_file=ioneq_file,abund_file=abund_file,path=path, $
                   noionrec=noionrec)
term=input.elvlcstr.data.full_level
n1=n_elements(input.ecm)

;----------X
; Calculate level populations and output rate coefficients to
; DATA_STR. 
;
pop_solver,input,temp,dens,pop,n_levels=n_levels, data_str=data_str, $
     radfunc=radfunc, /all_levels
pop=reform(pop(0,0,*))
n2=n_elements(pop)
nlvls2=data_str.nlev_ion
;----------X

;
; For 2-ion models, n2>n1. In this case we have to scale the
; populations (which are computed for all levels in the 2-ion model),
; such that the populations of the recombined ion sum to 1.
;
IF n1 LT n2 THEN BEGIN
   scale_factor=total(pop[0:n1-1])
   pop=pop/scale_factor
ENDIF 


IF level GT nlvls2 THEN BEGIN
  print,'% POP_PROCESSES: The input LEVEL must be less than '+trim(nlvls2+1)+'. Returning...'
  return
ENDIF 

;
; The new processes introduced in CHIANTI 9 are autoionization (ai),
; dielectronic capture (dc) and radiative recombination (rr). If they
; exist, then there should be arrays containing the rates.
;
IF n_elements(data_str.ai) EQ 1 THEN no_ai=1 ELSE no_ai=0
IF n_elements(data_str.dc) EQ 1 THEN no_dc=1 ELSE no_dc=0
IF n_elements(data_str.rr) EQ 1 THEN no_rr=1 ELSE no_rr=0

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

out={ rad_decay: 0., $
      e_exc: 0., $
      e_deexc: 0., $
      p_exc: 0., $
      p_deexc: 0., $
      ph_exc: 0., $
      ph_deexc: 0., $
      ai: 0. }
in={ rad_decay: 0., $
     e_exc: 0., $
     e_deexc: 0., $
     p_exc: 0., $
     p_deexc: 0., $
     ph_exc: 0., $
     ph_deexc: 0., $
     rr: 0., $
     dc: 0.}

;
; The following lines go through each process and extract the rate OUT of
; the level. The rates are stored in variables _OUT, and the total rate is
; stored in SUM.
;
IF level GT 1 THEN BEGIN
  a_leave=total(data_str.aa[level-1,0:level-2])
  a_out=a_leave*pop[level-1]
  sum=sum+a_leave*pop[level-1]
  out.rad_decay=a_out
ENDIF

IF level GT 1 THEN BEGIN
  c_leave=total(data_str.cc[level-1,0:level-2])
  e_deexc_out=c_leave*pop[level-1]
  sum=sum+c_leave*pop[level-1]
  out.e_deexc=e_deexc_out
ENDIF

IF level LT nlvls2 THEN BEGIN 
  c_leave=total(data_str.cc[level-1,level:*])
  e_exc_out=c_leave*pop[level-1]
  sum=sum+c_leave*pop[level-1]
  out.e_exc=e_exc_out
ENDIF ELSE BEGIN
  e_exc_out=0.
ENDELSE 

IF level GT 1 THEN BEGIN
  p_leave=total(data_str.ccp[level-1,0:level-2])
  p_deexc_out=p_leave*pop[level-1]
  sum=sum+p_leave*pop[level-1]
  out.p_deexc=p_deexc_out
ENDIF

IF level LT nlvls2 THEN BEGIN 
  p_leave=total(data_str.ccp[level-1,level:*])
  p_exc_out=p_leave*pop[level-1]
  sum=sum+p_leave*pop[level-1]
  out.p_exc=p_exc_out
ENDIF  ELSE BEGIN
  p_exc_out=0.
ENDELSE 


IF level GT 1 THEN BEGIN
  ph_leave=total(data_str.aax[level-1,0:level-2])
  stem_out=ph_leave*pop[level-1]
  sum=sum+ph_leave*pop[level-1]
  out.ph_exc=stem_out
ENDIF

IF level LT nlvls2 THEN BEGIN 
  ph_leave=total(data_str.aax[level-1,level:*])
  ph_out=ph_leave*pop[level-1]
  sum=sum+ph_leave*pop[level-1]
  out.ph_deexc=ph_out
ENDIF ELSE BEGIN
  ph_out=0.
ENDELSE

;
; This handles population lost through autoionization (introduced in
; CHIANTI 9).
;
IF NOT keyword_set(no_ai) THEN BEGIN
  ai_out=total(data_str.ai[level-1,*])*pop[level-1]
  sum=sum+ai_out
  out.ai=ai_out
ENDIF


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
IF NOT keyword_set(no_ai) THEN print,format='("  autoioniz.:   ",e10.2,f10.2,"%")', $
     ai_out,ai_out/sum*100.

print,'                  --------'
print,format='("           TOTAL",e10.2)',sum

print,''
ENDIF


;------------
;Now extract the individual rates coming INTO the level
;
IF NOT keyword_set(quiet) THEN BEGIN 
print,'Population entering level '+strtrim(string(level),2)
ENDIF 

sum=0.

IF level LT nlvls2 THEN BEGIN 
  a_enter=total(data_str.aa[level:*,level-1]*pop[level:*])
  sum=sum+a_enter
  in.rad_decay=a_enter
ENDIF ELSE BEGIN
  a_enter=0.
ENDELSE

IF level GT 1 THEN BEGIN
  e_exc_enter=total(data_str.cc[0:level-2,level-1]*pop[0:level-2])
  sum=sum+e_exc_enter
  in.e_exc=e_exc_enter
ENDIF

IF level LT nlvls2 THEN BEGIN 
  e_deexc_enter=total(data_str.cc[level:*,level-1]*pop[level:*])
  sum=sum+e_deexc_enter
  in.e_deexc=e_deexc_enter
ENDIF ELSE BEGIN
  e_deexc_enter=0.
ENDELSE 

IF level GT 1 THEN BEGIN
  p_exc_enter=total(data_str.ccp[0:level-2,level-1]*pop[0:level-2])
  sum=sum+p_exc_enter
  in.p_exc=p_exc_enter
ENDIF

IF level LT nlvls2 THEN BEGIN 
  p_deexc_enter=total(data_str.ccp[level:*,level-1]*pop[level:*])
  sum=sum+p_deexc_enter
  in.p_deexc=p_deexc_enter
ENDIF ELSE BEGIN
  p_deexc_enter=0.
ENDELSE 

IF level LT nlvls2 THEN BEGIN 
  stem_enter=total(data_str.aax[level:*,level-1]*pop[level:*])
  sum=sum+stem_enter
  in.ph_exc=stem_enter
ENDIF ELSE BEGIN
  stem_enter=0.
ENDELSE 

IF level GT 1 THEN BEGIN
  ph_enter=total(data_str.aax[0:level-2,level-1]*pop[0:level-2])
  sum=sum+ph_enter
  in.ph_deexc=ph_enter
ENDIF

;
; This handles population entering the level due to radiative
; recombination (RRLVL files).
;
IF  NOT keyword_set(no_rr) THEN BEGIN
  rr_enter=total(data_str.rr[level:*,level-1]*pop[level:*])
  sum=sum+rr_enter
  in.rr=rr_enter
ENDIF ELSE BEGIN
  rr_enter=0.
ENDELSE

;
; This handles population entering the level due to dielectronic
; capture. 
;
IF NOT keyword_set(no_dc) THEN BEGIN
  dc_enter=total(data_str.dc[level:*,level-1]*pop[level:*])
  sum=sum+dc_enter
  in.dc=dc_enter
ENDIF ELSE BEGIN
  dc_enter=0.
ENDELSE


;
; It's possible that status=1, but the ioniz/recomb coefficients aren't calculated because the ion fraction of the ion is zero at the input temperature. The check on data_str.correction below accounts for this case
;
IF tag_exist(input,'ionrec') eq 1 AND n_elements(data_str.correction) GT 1 THEN BEGIN
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
IF tag_exist(data_str,'rr') THEN BEGIN
  print,format='("  rad recomb:   ",e10.2,f10.2,"%")', $
     rr_enter, rr_enter/sum*100.
ENDIF 
IF tag_exist(data_str,'dc') THEN BEGIN
  print,format='("  diel capture: ",e10.2,f10.2,"%")', $
     dc_enter, dc_enter/sum*100.
ENDIF 
IF tag_exist(input,'ionrec') AND n_elements(data_str.correction) GT 1  THEN BEGIN
  print,format='("  ionization:   ",e10.2,f10.2,"%")', $
       ci_enter,ci_enter/sum*100.
  print,format='("  recomb.:      ",e10.2,f10.2,"%")', $
       rec_enter,rec_enter/sum*100.
ENDIF

print,'                  --------'
print,format='("           TOTAL",e10.2)',sum
ENDIF 

output={ ionname: ionname, $
         lvl: level, $
         lvl_name: trim(levstr), $
         density: dens[0], $
         temperature: temp[0], $
         pop: pop[level-1], $
         out: out, $
         in: in }


;
; In this section I print out the principal channels for populating
; the specified level. At present I've only coded radiative
; decay and electron excitation.
;
; 5-Mar-2019, PRY: I've had to introduce n2 below for the 2-ion
; models of CHIANTI 9. TERM is only for the single ion (n2 levels),
; whereas the other arrays are for the two-ion system (n levels). 
;
IF keyword_set(verbose) THEN BEGIN
  n=n_elements(pop)
  n2=n_elements(term)
  s={lev: 0, levstr: '', rate: 0., type: ''}
  rstr=replicate(s,n*2)
 ;
  rstr[0:n-1].rate=pop*data_str.aa[*,level-1]
  rstr[0:n-1].lev=indgen(n)+1
  rstr[0:n2-1].levstr=strpad(term,30,/after)
  rstr[0:n-1].type='cascade'
 ;
  rstr[n:2*n-1].rate=pop*data_str.cc[*,level-1]
  rstr[n:2*n-1].lev=indgen(n)+1
  rstr[n:n+n2-1].levstr=strpad(term,30,/after)
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

IF tag_exist(input,'ionrec') AND NOT keyword_set(quiet) THEN BEGIN
  print,''
  print,'** For this ion ionization and recombination processes are included when'
  print,'** working out the level populations. Due to the way these processes are'
  print,'** included (see CHIANTI v.5 paper for details) the populations entering and'
  print,'** leaving the above level may not balance.'
  print,''
ENDIF

END
