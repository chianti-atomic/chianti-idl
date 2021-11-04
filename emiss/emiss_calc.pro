;+
; NAME:
;      EMISS_CALC()
;       
; PURPOSE:
;      To compute the emissivities of all lines of a specified ion over
;      given ranges of temperature and density.
;
; CATEGORY:
;      CHIANTI; emissivity.
;
; EXPLANATION:
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
;       EMISS=EMISS_CALC (ION_NAME)  or
;       EMISS=EMISS_CALC (IZ, ION)
;
; EXAMPLES:
;	EMISS=EMISS_CALC(26,13)
;	EMISS=EMISS_CALC('fe_13')
;	EMISS=EMISS_CALC(26,13,temp=[6.2],dens=findgen(5)+8)
;	EMISS=EMISS_CALC('fe_13',temp=findgen(11)/100.+5.5,press=10.^15)
;
; INPUTS:
;	IZ	Either the atomic number of the ion, or the name of
;               the ion in CHIANTI format (e.g., 'fe_13'). If the
;               latter, then ION will be ignored.
;
;	ION	The spectroscopic number of the ion (e.g., 12 = XII)
;
; OPTIONAL INPUTS:
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
; KEYWORDS PARAMETERS:
;	NO_DE	Drops the hc/lambda factor in the computation of the 
;		emissivities. Useful for emission measure analyses involving 
;		photon fluxes
;
;       NOPROT  If set, then the default setting will be NOT to use 
;               proton rates. This can be changed within the routine.
;
;	QUIET	This keyword has been disabled as the routine is now
;       	quiet by default. Please use /VERBOSE to print
;       	information messages.
;
;       DIEL    If the dielectronic recombination files exist for the ion, 
;               then these are used to derive the emissivities.
;
;       NO_SETUP  If emiss_calc is called from a routine where the ELEMENTS
;                 common block has already been set up, then this keyword
;                 stops emiss_calc loading up the common block.
;
;       NO_CALC  If set, then the emissivities are not calculated but
;                everything else in the output structure is populated.
;
;       NOIONREC: If set, then level-resolved ionization/recombination
;                rates are not read for the ion, even if they exist.
;       NO_RREC: If set, do not read the level-resolved radiative
;                recombination (RR) files.
;       VERBOSE: If set, then print information messages to the screen.
;
;
; OUTPUT:
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
;	Transitions where only theoretical energies are available for at 
;	least one of the two levels are assigned negative wavelengths in 
;	the .wgfa file. With emiss_calc, the wavelength is set to be 
;	positive, but emiss.flag is set to -1 for that transition, so 
;	that other routines can keep track of the theoretical wavelengths.
;
; CALLS:
;	CH_SETUP_ION, POP_SOLVER, CONVERTNAME,
;       ZION2FILENAME, ZION2SPECTROSCOPIC, CH_TMAX
;
; HISTORY:
;	Ver 1, PRY 28-Jun-97
;	Ver 2, PRY 26-Jul-97  - corrected problem with size of emiss
;	Ver 3, PRY 22-Sep-97  - allowed photo-excitation to be included
;	Ver 4, PRY 6-Jul-98   - added PATH
;	Ver 5, PRY 5-Sep-98   - added call to choose_ioneq
;       Ver 6, PRY 3-Dec-98   - dosen't crash if no params given
;	Ver 7, PRY 9-Jan-99   - allowed proton rates to be added through 
;				/PROTON keyword.
;	Ver 8, PRY 10-Feb-99  - added /QUIET keyword
;	Ver 9, PRY 8-Oct-99   - for H-like ions, there's a 2-photon 
;				transition with a non-zero A-value that is 
;				assigned a zero wavelength (as it does not 
;				produce an emission line). This caused 
;				problems for dens_plotter, so emiss_calc 
;				now removes this transition from the emiss 
;				structure.
;	Ver 10, PRY 15-Dec-99 -	added deu to the upsilon common block in 
;				order to be consistent with the main Chianti 
;				routines.
;	Ver 11, PRY 8-May-00  - added PRESSURE
;       Ver 12, PRY 17-Aug-00 - changed elvlc common block to match new 
;                               version of pop_solver
;       Ver 13, PRY 10-Oct-00 - now calls setup_ion to read ion data
;       Ver 14, PRY 1-Jun-00  - removed call to choose_ioneq; now makes use 
;                               of the !ioneq_file system variable.
;       Ver 15, PRY 25-Sep-01 - modified for 9-point splines and proton rates
;       Ver 16, PRY 9-Dec-01  - completed changes for v.4
;
;       V. 17, 29-May-2002, Giulio Del Zanna (GDZ) 
;                   generalized directory concatenation to work for
;                   Unix, Windows  and VMS. 
;                   Now we only call zion2filename, corrected the call to
;                   zion2spectroscopic (dielectronic cases were not handled
;                   correctly).
;
;       V.18, 6-Aug-2002, Peter Young
;                   Theoretical wavelengths weren't being flagged so this 
;                   has now been corrected.
;
;       V.19, 7-Aug-2002, Peter Young
;                   Stopped "lines" with zero wavelength being included in 
;                   the structure.
;                   Changed the use of CHIANTI system variables.
;
;       V.20 14-Aug-2002, GDZ 
;                   Modified the use of Version, to make it compatible with the
;                   other CHIANTI v.4 routines.
;
;       V.21 10-Sep-2002, Peter Young
;                   Allowed a density of 1 cm-3 to be input.
;
;       V.22 7-Aug-2003, Peter Young
;                   Added keyword /NO_SETUP
;
;       V.23  4-Oct-2003, GDZ
;                  modified the input to POP_SOLVER, now it is dealt with an
;                  input structure.
;
;       V.24  5-Mar-2004, Enrico Landi (EL)
;                  included ionization and recombination as level population
;                  processes
;
;       V.25  6-May-2005, EL
;                  corrected a minor incompatibility with IDL 5.6
;
;       V.26  5-Jul-2005, Peter Young
;                  added RADFUNC= and SUM_MWL_COEFFS= keywords
;
;       V.27  1-Aug-2005, Peter Young
;                  re-ordered code for setup to prevent crash
;
;       V.28  3-Aug-2005,GDZ
;                  fixed bug, only define ioneq_ionrec when files are
;                  there (it was failing with neutrals)
;
;       V.29 , 7-Dec-07, GDZ changed a for loop to long to allow long
;              lists of lines.
;
;       V.30 , 12-Jun-2009, Enrico Landi
;               Changed the definition of the temperature array for ion fractions
;               in the IONREC variable, now taken directly from the output of
;               READ_IONEQ.PRO
;
;       V.31, 28-Apr-2013, Peter Young
;               The first call to pop_solver (which is only used to
;               get the size of the population array) is now made with
;               the /no_calc keyword in order to speed it up.
;
;       V.32, 25-Feb-2014, Peter Young
;               Added /no_calc to emiss_calc.
;
;       v.33, 26-Jan-2018, Peter Young
;               Removed common blocks and added call to ch_setup_ion;
;               IZ can now be the ion name
;
;       v.34   12-May-2018 Giulio Del Zanna
;              major revision for version 9. 
;       v.35   24 Jul 2018 GDZ minor revisions.
;       v.36   10 Oct 2018 GDZ, replaced call to ch_load_ion_rates
;              with a call to ch_setup_ion
;
;       v.37   16 Nov 2018 GDZ, fixed a  bug when T is a single value
;       and N is an array. Reinstated the 'dielectronic' option, which
;       was the default for all ions in v.8 (and previous versions).
;
;       v.38   5 Dec 2018 GDZ revised handling of pop_solver.
;       v.39   14-Dec-2018 GDZ, added NOIONREC, NO_RREC keywords.
;       v.40, 19-Nov-2020, Peter Young
;           Added keyword /VERBOSE and disabled /QUIET so that routine
;           no longer prints information to the IDL window (unless you
;           ask for it).
;       v.41, 04-Nov-2021, Peter Young
;           Restored original behavior, whereby transitions with
;           negative wavelength in the .wgfa file become positive in
;           the emiss structure but have flag set to -1.
;
; VERSION     :   41
;
;-

FUNCTION  EMISS_CALC, IZ, ION, TEMP=TEMP, DENS=DENS, RADTEMP=RADTEMP, $
                      RPHOT=RPHOT, $
                      PATH=PATH, NO_DE=NO_DE, NOPROT=NOPROT, $
                      QUIET=QUIET, PRESSURE=PRESSURE, DIEL=DIEL, $
                      IONEQ_FILE=IONEQ_FILE, ABUND_FILE=ABUND_FILE, $
                      NO_SETUP=NO_SETUP, SUM_MWL_COEFFS=SUM_MWL_COEFFS, $
                      RADFUNC=RADFUNC, NO_CALC=NO_CALC, noionrec=noionrec, $
                      no_rrec=no_rrec, verbose=verbose




IF N_PARAMS() LT 1 THEN BEGIN
   PRINT,'Use:  IDL> emiss=emiss_calc(iz,ion [, temp=temp, dens=dens, '
   PRINT,'                              radtemp=radtemp, rphot=rphot, /no_de,'
   PRINT,'                              path=path, /noprot, $ '
   PRINT,'                              pressure=pressure, ioneq_file=ioneq_file, $'
   print,'                              abund_file=abund_file, /no_calc, /verbose ] )'
   RETURN,0
ENDIF


quiet=1b-keyword_set(verbose)
;IF  KEYWORD_SET(quiet) THEN verbose=0 else verbose=1

;
; PRY, 26-Jan-2018
; If the 1st parameter is a string, then it is interpreted as the
; CHIANTI ion format, e.g., 'o_6'.
;
IF datatype(iz) EQ 'STR' THEN BEGIN
   iz_save=iz
   convertname,iz_save,iz,ion
ENDIF


IF n_elements(pressure) NE 0 AND n_elements(sum_mwl_coeffs) NE 0 THEN BEGIN
   print,'% EMISS_CALC: the PRESSURE and SUM_MWL_COEFFS keywords are incompatible with'
   print,'              each other. Please use one or the other.'
   return,0
ENDIF

;
; The following extracts the names of the files to be read. I want to allow 
; a different path to be chosen, and I extract only the information I need 
; from filename.
IF keyword_set(diel) THEN diel = 1 ELSE diel = 0

zion2filename,iz,ion,filename, name=name, diel=diel
gname = name


;
; Set ioneq_file.
; Note the file is used to get Tmax (if necessary) and for the proton
; density.  
;
IF n_elements(ioneq_file) EQ 0 THEN ioneq_file=!ioneq_file


;
; Set abundance file
;  The file is used for the proton density.
;
IF n_elements(abund_file) EQ 0 THEN abund_file=!abund_file

;------------------------------------+

; ****** NOTE THAT TEMP IS LOG T ************

; If TEMP not specified, then set TEMP=[log(T_max),log(T_max)+-0.15]
;
IF N_ELEMENTS(temp) EQ 0 THEN BEGIN
   tmax=ch_tmax(gname,/log,ioneqname=ioneq_file)
   temp=[tmax-0.15,tmax,tmax+0.15]
ENDIF
nt=N_ELEMENTS(temp)
;------------------------------------+

;---------------X
sumtst=0
IF n_elements(sum_mwl_coeffs) NE 0 THEN BEGIN
   IF n_elements(temp) NE n_elements(sum_mwl_coeffs) THEN BEGIN
      print,'% EMISS_CALC: SUM_MWL_COEFFS must have the same size as TEMP. Returning...'
      return,0
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

;;-----------------<>
; Default density range is 8 to 12 in .5 intervals:
;
IF n_elements(dens) EQ 0 THEN dens=findgen(9)/2. + 8.
nd=n_elements(dens)
;;-----------------<>


IF N_ELEMENTS(pressure) NE 0 THEN BEGIN
   nd=1
   dens=ALOG10(pressure/10.^temp)
ENDIF

;------------------------------------------------------------------------------------
;  This loads up the ion's atomic data  

    input=ch_setup_ion(gname,rphot=rphot,radtemp=radtemp,noprot=noprot, $
                        ioneq_file=ioneq_file,abund_file=abund_file,path=path, $
                        quiet=quiet, noionrec=noionrec, no_rrec=no_rrec)


zion2spectroscopic,iz,ion,spectname

IF tag_exist(input,'ionrec')  AND n_elements(sum_mwl_coeffs) NE 0 THEN BEGIN
   print,'% EMISS_CALC: this ion ('+spectname+') has contributions from ionization'
   print,'              and recombination that can not be included when the SUM_MWL_COEFFS'
   print,'              keyword is set. The ionization and recombination processes'
   print,'              will be switched off.'

;   status=-1                    ; v.8  switch. In V.9 the switch is
;   done by checking for the existence of the IONREC tag:

  input=rem_tag(input,'ionrec')

ENDIF
;---------------O


; GDZ - calculate  the size of the population array.
  n_elvl=n_elements(input.ecm)
  n_wgfa=max([input.wgfastr.lvl1,input.wgfastr.lvl2])
  usize=max([input.splstr.data.lvl1,input.splstr.data.lvl2])

  IF tag_exist(input,'autostr') THEN asize=max(input.autostr.lvl2) ELSE asize=usize
  ausize=max([asize,usize])

  nlvls2=min([n_elvl,ausize,n_wgfa])

;
; Reduce wgfa down to those transitions that have a non-zero
; wavelength (i.e., remove two-photon  transitions),
; and any transitions beyond the pop_solver model (in case wgfa file
; has more transitions than needed).
; I've also added a check on the A-value being non-zero as we
; seem to have quite a few of these transitions in CHIANTI.
;
wgfa=input.wgfastr
ind=where(wgfa.wvl NE 0. AND wgfa.lvl1 LE nlvls2 AND wgfa.lvl2 LE nlvls2 AND wgfa.aval NE 0.)
wgfa=wgfa[ind]
ntrans2=n_elements(wgfa)


;----------------------------------o
; Create the EMISS structure, using ntrans2 to set size
;
verstr=ch_get_version()
;
zion2spectroscopic,iz,ion,species, diel=diel
;
IF keyword_set(sumtst) THEN em=dblarr(1,nd) ELSE em=dblarr(nt,nd)
;
str={   ion_name: species, $
        lambda: 0., $
        level1: 0, $
        lvl1_desc: '', $
        level2: 0, $
        lvl2_desc: '', $
        flag:   0, $
        em: em, $
        version: verstr $
    }
;
emiss=replicate(str,ntrans2)
;----------------------------------o

;
; The following sets all wavelengths to be positive, but transitions
; with negative wavelengths have flag=-1.
;
k=where(wgfa.wvl LT 0.,nk)
IF nk NE 0 THEN emiss[k].flag=-1
emiss.lambda=abs(wgfa.wvl)
;
emiss.level1=wgfa.lvl1
emiss.level2=wgfa.lvl2
emiss.lvl1_desc=input.elvlcstr.data[wgfa.lvl1-1].full_level
emiss.lvl2_desc=input.elvlcstr.data[wgfa.lvl2-1].full_level


IF NOT KEYWORD_SET(quiet) THEN BEGIN
   PRINT,''
   PRINT,'Log_10 temperatures...'
   PRINT,FORMAT='("   ",9f6.2)',temp
   PRINT,''
   PRINT,'Log_10 densities...'
   PRINT,FORMAT='("   ",9f6.2)',dens
   
   verbose=1
ENDIF else verbose=0

IF NOT keyword_set(no_de) THEN mult_factor=1.986d-8/emiss.lambda*wgfa.aval $
ELSE mult_factor=wgfa.aval

;
; The following code gets the level populations and computes the
; emissivities (if /no_calc not set).
;
IF NOT keyword_set(no_calc) THEN BEGIN 
   
   IF N_ELEMENTS(pressure) NE 0 THEN BEGIN
                                ;
      pop_solver, input, 10.^temp,10.^dens,pop, radfunc=radfunc, /pressure,verbose=verbose
                                ;
      emiss.em=pop[*,emiss.level2-1]
      FOR i=0,nt-1 DO emiss.em[i]=emiss.em[i]*mult_factor
   ENDIF ELSE BEGIN
                                ;
      pop_solver, input, 10.^temp,10.^dens,pop, sum_mwl_coeffs=smc, $
                  radfunc=radfunc,verbose=verbose

           ;
      IF sumtst EQ 1 THEN BEGIN
                                ;
         IF nd NE 1 THEN emiss.em[0,*]=pop[*,*,emiss.level2-1]
         IF nd EQ 1 THEN emiss.em[0,*]=reform(pop[*,*,emiss.level2-1],1, $
                                              n_elements(emiss.level2))
                                ;
      ENDIF ELSE BEGIN

    IF nt EQ 1 THEN BEGIN
      IF nd NE 1 THEN emiss.em=pop[*,*,emiss.level2-1] ELSE $
           emiss.em=reform(pop[*,*,emiss.level2-1],1,n_elements(emiss.level2))
    ENDIF ELSE emiss.em=reform(pop[*,*,emiss.level2-1])
                                ;
     ENDELSE 
                                ;
      FOR i=0L,ntrans2-1 DO emiss[i].em=emiss[i].em*mult_factor[i]
   ENDELSE
ENDIF 


ind=sort(emiss.lambda)
emiss=emiss[ind]

;
; This restores the input value of IZ if it was specified as a string.
;
IF n_elements(iz_save) NE 0 THEN iz=iz_save

;
; Tidy up
;
junk=temporary(input)
junk=temporary(wgfa)

RETURN,emiss

END
