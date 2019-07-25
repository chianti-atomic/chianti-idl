
function ch_critical_density, ion_name, level, temp=temp, name=name, quiet=quiet

;+
; NAME
;
;      CH_CRITICAL_DENSITY()
;
; PROJECT
;
;      CHIANTI
;
; EXPLANATION
;
;      This routine computes the critical density for an atomic level
;      within an ion. The critical density is the electron number
;      density at which the rate at which electron collisions
;      de-populate a level exactly balances the rate at which
;      radiative decays de-populate the level.
;
;      The definition for the critical density comes from Osterbrock &
;      Ferland ("Astrophysics of Gaseous Nebulae and Active Galactic
;      Nuclei", 2nd Ed., University Science Books, 2006).
;
; INPUTS
;
;      ION_NAME   The ion name in CHIANTI format. E.g., 'fe_13' for Fe
;                 XIII. 
;
;      LEVEL      The CHIANTI index number for the level. See the
;                 .elvlc data files to get the index numbers.
;
; OPTIONAL INPUTS
;
;      TEMP        The temperature for which the critical density is
;                  required. Units: K.  Note that if TEMP is not
;                  specified then the T_max of the ion will be used
;                  (as given by the routine ch_tmax.pro).
;
; KEYWORDS
;
;      QUIET       If set, then no information will be printed to the
;                  screen.
;
; OPTIONAL OUTPUTS
;
;      NAME        The name of the level. E.g., '3s2.3p2 3P1'.
;
; OUTPUTS
;
;      The critical density. This is the electron number density in
;      units of cm^-3.
;
; CALLS
;
;      POP_PROCESSES, CONVERTNAME, ZION2FILENAME, READ_ELVLC, CH_TMAX
;
; HISTORY
;
;      Ver.1, 9-May-2013, Peter Young
;      Ver.2, 14-May-2013, Peter Young
;         changed get_tmax to ch_tmax
;-

IF n_elements(temp) EQ 0 THEN BEGIN
  temp=ch_tmax(ion_name)
ENDIF 

convertname,ion_name,iz,ion
zion2filename,iz,ion,fname
elvlcname=fname+'.elvlc'

read_elvlc,elvlcname,l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,ref

name=trim(term[level-1])

IF NOT keyword_set(quiet) THEN BEGIN
  print,format='("Log10 (T/K):   ",f6.2)',alog10(temp)
  print,'Selected level:  ',name
ENDIF 

dens=1e10
pop_processes,ion_name,temp=temp,data=data,dens=dens,/quiet

;
; Note that cc = N_e * q, where q is the electron excitation rate
; coefficient. 
;
;  aa[i,j]  corresponds to decay from level i to level j
;  cc[i,j]  corresponds to excitation from level i to level j
;
aa=data.aa   ; radiative decay rates
cc=data.cc   ; electron excitation rates

nc=total(aa[level-1,*])/total(cc[level-1,*])*dens

return,nc

END

