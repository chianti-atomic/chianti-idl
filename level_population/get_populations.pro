;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. 
;
;
; NAME:
;	GET_POPULATIONS
;
; PURPOSE:
;	get the population of a number of the lowest levels as a function of 
;       electron density for a specific temperature
;
; CATEGORY:
;
;	science.
;
; CALLING SEQUENCE:
;
;       GET_POPULATIONS,Ion,T,Nlevels
;
;
; INPUTS:
;
;       gname:  CHIANTI style name for the ion, i.e., 'c_6' for C VI
;       T:  electron temperature (K)
;
; OPTIONAL INPUT:
;
;       NLEVELS: the maximum number of levels displayed. If not
;                passed, all the levels are displayed.
;
;       DENSITIES: the array  of electron densities at which the
;                  populations are calculated. 
;
;       Path:    This directly specifies the path where the
;                ion's data files are stored. If not set, then
;                the files are taken from the user's CHIANTI
;                distribution. 
;
;       Radtemp: If photon excitation is included (by defining RPHOT),
;                then this input specifies the blackbody radiation
;                temperature in K. If not specified, then it is set to
;                6000 K.
;       Rphot:   Distance from the centre of the star in stellar radius units.
;                That is, RPHOT=1 corresponds to the star's
;                surface. If RPHOT is not specified, then photon
;                excitation will be switched off when pop_solver is
;                called.
;
; KEYWORDS:
;
;       NOIONREC: If set, then level-resolved ionization/recombination
;                rates are not read for the ion, even if they exist.
;
;       NO_RREC: If set, do not read the level-resolved radiative
;                recombination (RR) files.
;       NOPROT:  If set, then proton rates are not read for the ion,
;                even if they exist.
;
;       v.1, 23 Nov 2018, Giulio Del Zanna (GDZ)
;       v.2, 14 Dec 2018, GDZ, added keywords       
;       v.3  16 Jan 2019, GDZ added double
;
; VERSION     : 3 
;
;-
pro get_populations, gname, t, nlevels, densities=densities, $
                      outfile=outfile, noionrec=noionrec, no_rrec=no_rrec, $
                     noprot=noprot,radtemp=radtemp,rphot=rphot, $
                    PATH=path, verbose=verbose

t=double(t)
densities=double(densities)

if n_elements(verbose) eq 0 then verbose=1 

if verbose then quiet=0 else quiet=1

;  This loads up the ion's atomic data  

input=ch_setup_ion(gname,rphot=rphot,radtemp=radtemp,noprot=noprot, $
                        ioneq_file=ioneq_file,abund_file=abund_file,path=path, $
                        quiet=quiet,   noionrec=noionrec, no_rrec=no_rrec )


; either the nlevels or levels should be defined.

IF n_elements(nlevels) EQ 0 THEN $ 
nlevels = n_elements(input.ecm)

level_numbers = indgen(nlevels)+1 

term=input.elvlcstr.data.conf+' '+input.elvlcstr.data.level


print,' ' 
print,' Levels: '
FOR  i=0,n_elements(level_numbers)-1 DO $
  print,level_numbers[i],strpad(term(level_numbers[i]-1),30,/after),format='(i5,4x,a30)'


; define the densities

IF n_elements(densities) eq 0 then densities=10.^(indgen(10)+4)

nd=n_elements(densities)

;  calculate level populations
;------------------------------

     pop_solver, input, t , densities ,pop, radfunc=radfunc,verbose=verbose
                                ;

pop = REFORM(pop)
pops = TRANSPOSE(pop[*,0:nlevels-1])

;
print,' '
print,'Density       Populations'
print,' '

FOR i=0,nd-1 DO begin 

     print,' '
     print, 'Density: '+trim(densities[i]) 

      FOR j=0, nlevels-1 DO $
        print, string(trim(level_numbers[j]), pops[level_numbers[j]-1,i], $
                          format='(" ",i4,": ",500e10.3)')

endfor 

; save the values:

if n_elements(outfile) eq 0 then outfile='pops.save'

save, file=outfile, t , densities ,pops


end
