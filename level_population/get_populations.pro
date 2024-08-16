;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. 
;
;        This program was developed as part of CHIANTI-VIP. 
;        CHIANTI-VIP (Virtual IDL and Python) is  a member of the CHIANTI family
;        mantained by Giulio Del Zanna, to develop additional features and
;        provide them to the astrophysics community.
;
;
;
; NAME:
;
;	GET_POPULATIONS
;
;
; PURPOSE:
;
;	get the population of a number of the lowest levels as a function of 
;       electron density for a specific temperature
;
;
; CATEGORY:
;
;	Atomic processes; level populations
;
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
;
; OUTPUTS:
;
;       The populations of the specified number of energy levels at the
;       input temperature over a range of densities
;
;
; OPTIONAL OUTPUTS:
;
;       NONE
;
;
; CALLS:
;
;       CH_SETUP_ION
;       POP_SOLVER
;
;
; WRITTEN:
;         
;       v.1, 23 Nov 2018, Giulio Del Zanna (GDZ)
;
;
; MODIFIED:
;
;       v.2, 14 Dec 2018, GDZ, added keywords       
;       v.3  16 Jan 2019, GDZ added double
;       v.4  24 Aug 2023, Roger Dufresne, output population array if no
;               outfile set. Uses the keyword pressure when pop_solver
;               is to solve a density array that is connected to
;               the temperature array.
;       v.5 14 Sept 2023, GDZ
;           modified the way the number of levels is dealt with. 
;       v.6  04 Jan 2024, RPD
;               stopped level population output to screen in verbose mode
;
;
; VERSION     : 6
;
;-

pro get_populations, gname, t, n_levels, pops=pops, densities=densities, $
                      outfile=outfile, noionrec=noionrec, no_rrec=no_rrec, $
                     noprot=noprot,radtemp=radtemp,rphot=rphot, $
                    PATH=path, verbose=verbose, no_auto=no_auto, pressure=pressure

t=double(t)

; define the densities

IF n_elements(densities) eq 0 then densities=10.d^(indgen(10)+4) else densities=double(densities)
nd=n_elements(densities)


if n_elements(verbose) eq 0 then verbose=0 
if verbose then quiet=0 else quiet=1

;  This loads up the ion's atomic data  
; if no_auto is set, n_levels will include only the bound states; n_levels gets modified by 
; ch_setup_ion

input=ch_setup_ion(gname,rphot=rphot,radtemp=radtemp,noprot=noprot, $
                        ioneq_file=ioneq_file,abund_file=abund_file,path=path, $
                   quiet=quiet,  no_auto=no_auto,  noionrec=noionrec, no_rrec=no_rrec,$
                   n_levels=n_levels) ; GDZ 


if keyword_set(verbose) then begin 

level_numbers = indgen(n_levels)+1 
term=input.elvlcstr.data.conf+' '+input.elvlcstr.data.level

   print,' ' 
print,' Levels: '
FOR  i=0,n_elements(level_numbers)-1 DO $
  print,level_numbers[i],strpad(term(level_numbers[i]-1),30,/after),format='(i5,4x,a30)'

endif 

;  calculate level populations
;------------------------------

     pop_solver, input, t , densities ,pop, n_levels=n_levels, radfunc=radfunc,verbose=verbose,$
                    noionrec=noionrec, no_rrec=norrec, no_auto=no_auto, pressure=pressure
                                ;
         
pop = REFORM(pop)
if size(pop,/n_dim) eq 1 then pops=pop $
  else if size(pop,/n_dim) eq 2 then pops=TRANSPOSE(pop[*,0:n_levels-1]) $
  else if size(pop,/n_dim) eq 3 then pops=TRANSPOSE(pop[*,*,0:n_levels-1])

;
;if keyword_set(verbose) then begin
;  print,' '
 
;  FOR i=0,nd-1 DO begin 

;      print,' '
;      print, gname+ ' Temperature='+string(t[i],format='(e12.4)')+$
;             '  Density='+string(densities[i],format='(e12.4)') 
;      print, 'Level:  relative population '
      
;        FOR j=0, n_levels-1 DO $
;          print, string(trim(level_numbers[j]), pops[level_numbers[j]-1,i], $
;                            format='(" ",i4,": ",500e10.3)')

;  endfor 
;endif


; save the values:

if n_elements(outfile) eq 1 then $
  save, file=outfile, t , densities ,pops


end
