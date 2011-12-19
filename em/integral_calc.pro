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
; NAME: INTEGRAL_CALC
;       
; PURPOSE:
;
;       To compute the atomic data integral for use in column or volume
;	emission measure work.
;
; CATEGORY:
;
;       Scientific analysis
;
; EXPLANATION:
;
; 	Defining
;
;	G(T) = 0.83 * Fr(T) * N_j * A_ji / N_e
;
;	where Fr(T) is the ionisation fraction (e.g., from Arnaud & 
;	Rothenflug), N_j the relative population of level j, A_ji the 
;	A-value for the j->i transition and N_e the electron density. The 
;	0.83 is the ratio of hydrogen to free electrons, which is constant 
;	above around 10^4 K. This function is sharply-peaked at a 
;	temperature T_mem (the temperature of maximum emission, which can 
;	be different from the temperature of maximum ionisation, T_max) 
;	and so an oft-used approximation is to take G(T) constant in the 
;	range log T_mem - 0.15 to log T_mem + 0.15 and zero outside. The 
;	value of the constant, which we call C_lambda here, is then given 
;	by
;
;	C_lambda =    integral { G(T) dT }
;                     --------------------
;                  T_mem * (10^0.15 - 10^-0.15)
;
;	If EM(s) is the column emission measure, F the flux (erg cm-2 s-1)
;	in a line lambda, Ab the abundance of the element and DE (erg) the
;	energy for the transition, then:
;
;	F = DE * Ab * C_lambda * EM(s)
;
;	If we are dealing with intensities I (erg cm-2 s-1 sr-1) then:
;
;	4pi * I = DE * Ab * C_lambda * EM(s)
;
;	This program extracts the ionisation balance and emissivities from 
;	the CHIANTI database and calculates C_lambda for all lines in the 
;	specified wavelength interval WRANGE by integrating over 
;	0.02 dex temperature intervals.
;
;	The C_lambda functions for all the lines in the selected wavelength 
;	range WRANGE are displayed as well as the temperature of maximum 
;	emission (T_mem), DE*C_lambda and 4pi/(DE*C_lambda). These latter 
;	two quantities are useful for the emission measure analysis.
;
;	Any combination of the displayed lines can then be blended and the 
;	corresponding quantities for the blend will be displayed.
;
;	The function Fr(T) * N_j * A_ji can also be plotted at this stage.
;
; CALLING SEQUENCE:
;
;       INTEGRAL_CALC, IZ, ION, [WRANGE=WRANGE, /CHOOSE, DENS=DENS]
;
; EXAMPLES:
;
;	INTEGRAL_CALC, 26, 13, WRANGE=[200,205], /CHOOSE
;
;	INTEGRAL_CALC, 14, 10, WRANGE=[250,270], DENS=10.
;
; INPUTS:
;
;	IZ:	The atomic number of the ion
;	ION:	The spectroscopic number of the ion (e.g., 12 = XII)
;
; OPTIONAL INPUTS:
;
;	DENS:	The density at which the emissivities are calculated 
;		(default=10.)
;	WRANGE: Wavelength range for which C_lambda functions are 
;		calculated. If not given, then the 10 strongest lines 
;		are printed.
;	INDEX:	Particular elements in the emissivity structure can be 
;		selected with INDEX. This allows integral_calc to be run 
;		'silently'. The output is contained in the OUTSTR structure. 
;		If index is given as, e.g., [7,8], then the C_lambda 
;		functions for these two lines are summed and output.
;	PATH:	Directly specify the directory path where the Chianti data 
;		for the ion is found
;
;       ABUND_FILE  The name of a CHIANTI abundance file. This is used for 
;               calculating the proton to electron ratio. Default is 
;               !abund_file.
;
;       IONEQ_FILE  The name of a CHIANTI ion balance file. This is used for 
;               calculating the proton to electron ratio and evaluating 
;               the T_max of the ion. Default is !ioneq_file.
;
;       RADTEMP  The blackbody radiation field temperature (default 
;                6000 K).
;
;       RPHOT    Distance from the centre of the star in stellar radius 
;                units. I.e., RPHOT=1 corresponds to the star's surface. 
;                (Default is infinity, i.e., no photoexcitation.)
;
; KEYWORDS:
;
;	CHOOSE:	Allow ion balance calculations to be selected manually 
;		(see choose_ioneq.pro routine).
;
; OPTIONAL OUTPUTS:
;
;	OUTSTR:	A structure with the following tags
;
;		.tmem	- the T_mem for the line(s)
;		.dec	- total( de * c_lambda )
;		.pidec	- 4 * pi / total( de * c_lambda )
;
;	Only output when INDEX is specified.
;
; COMMON BLOCKS:
;
;	None.
;
; CALLS:
;
;	EMISS_CALC, CH_GET_FILE, READ_IONEQ, GET_IEQ
;
; HISTORY:
;
;	Ver.1: PRY, 28-JUN-97.
;	Ver.2: PRY, 7-OCT-97. Added TEMPI and GOFT, for plotting.
;	Ver.3: PRY, 31-JUL-98. Added PATH.
;	Ver.4: PRY, 6-APR-99. Added INDEX, OUTSTR. Removed TEMPI and GOFT 
;		(these can be got from the g_of_t routine).
;       Ver.5: PRY, 9-Dec-01. Modified for v.4 of CHIANTI.
;
;       V.6, 21-May-2002, Giulio Del Zanna (GDZ) 
;                   generalized directory concatenation to work for
;                   Unix, Windows  and VMS. 
;
;       V.7, 06-Aug-02 GDZ
;              Changed the use of CHIANTI system variables. 
;
;       V.8, 6-Feb-2006, Peter Young
;              Integration is now performed on d(logT) intervals, rather than
;              dT intervals following problems identified by Luca Teriaca.
;       
;
; VERSION     :   8, 06-Feb-06

;-

PRO INTEGRAL_CALC, IZ, ION, WRANGE=WRANGE, DENS=DENS, CHOOSE=CHOOSE, $
                   PATH=PATH, OUTSTR=OUTSTR, INDEX=INDEX, $
                   IONEQ_FILE=IONEQ_FILE, RPHOT=RPHOT, RADTEMP=RADTEMP, $
                   NOPROT=NOPROT, ABUND_FILE=ABUND_FILE


IF N_PARAMS() LT 2 THEN BEGIN
  PRINT,'Use: IDL> integral_calc, iz, ion [, wrange= ,/choose, /noprot, '
  PRINT,'                           dens= , index= , radtemp= , '
  PRINT,'                           rphot= , outstr= , path= , '
  PRINT,'                           ioneq_file= , abund_file= ]'
  RETURN
ENDIF

;---------------------------------------<0>

dir=concat_dir(!xuvtop,'ioneq')

CASE 1 OF
 ;
  n_elements(ioneq_file) NE 0: ioneq_name=ioneq_file
 ;
  keyword_set(choose): ioneq_name=ch_get_file(path=dir,filter='*.ioneq', $
                           title='Select Ionization Equilibrium File')
 ;
  ELSE: ioneq_name= !ioneq_file
 ;
ENDCASE

read_ioneq,ioneq_name,temp_all,ioneq,ref
ioneq=double(ioneq)

IF N_ELEMENTS(index) EQ 0 THEN BEGIN
  PRINT,''
  PRINT,'Ion balance file:   ',ioneq_name
ENDIF

f_all=ioneq[*,iz-1,ion-1]
ind=WHERE(f_all NE 0.)
f=f_all(ind)
temp=temp_all(ind)

; use a finer T scale for working out the integral
;
nt=n_elements(temp)
tempi=dindgen( (nt-1)*5 +1 )/50. + min(temp)

; calculate ion fraction on this scale
;
frac=get_ieq(10.^tempi,iz,ion,ioneq_logt=temp_all,ioneq_frac=ioneq)


;--------<
; Pick density
;
IF NOT KEYWORD_SET(dens) THEN dens=10.
IF N_ELEMENTS(index) EQ 0 THEN BEGIN
  PRINT,FORMAT='("Log10 density: ",f7.1)',dens
ENDIF
;--------<

;-----------------+
; Create the emissivity structure. Note that, to save time, emissivities
; are calculated at 0.1 dex logT intervals. Also note that the returned
; emissivities do not contain Delta-E.
;
emiss=emiss_calc(iz,ion,dens=dens,temp=temp,/no_de,path=path,/quiet, $
                 noprot=noprot, ioneq_file=ioneq_name, $
                 abund_file=abund_file, rphot=rphot, radtemp=radtemp)
;-----------------+

;---------------------------[
; Get info from emiss
;
IF N_ELEMENTS(index) EQ 0 THEN BEGIN
  IF N_ELEMENTS(wrange) NE 0 THEN BEGIN
    ind=WHERE( (emiss.lambda GE wrange(0)) AND (emiss.lambda LE wrange(1)))
    IF ind[0] EQ -1 THEN BEGIN
      print,''
      print,'% INTEGRAL_CALC: no lines in specified range. Returning...'
      return
    ENDIF
  ENDIF ELSE BEGIN
    n_em=N_ELEMENTS(emiss(0).em)
    ind=REVERSE( SORT (emiss.em(FIX(n_em/2))) )
    ind=ind(0:9)
  ENDELSE
ENDIF ELSE BEGIN
  ind=index
ENDELSE
waves=emiss(ind).lambda
;
n=N_ELEMENTS(waves)
c_lambda=fltarr(n)
t_mem=fltarr(n)
;---------------------------[


;-------------------0
; Work out C_lambda for each line in the specified wavelength range
;
FOR i=0,n-1 DO BEGIN
  em=emiss[ind[i]].em
 ;
  y2=spl_init(temp,em)
  emi=spl_interp(temp,em,y2,tempi)
 ;
 ; dT = alog(10) * T * d(logT), and d(logT)=0.02
 ;
  integrand=int_tabulated(tempi,emi*frac*10.^tempi*alog(10d0),/double)
  get_imax=MAX(emi*frac,imax)
  t_mem(i)=tempi(imax)
  c_lambda(i)=0.83*integrand/(0.7046*10.^(tempi(imax))*10.^(dens))
ENDFOR
;-------------------0

nindex=N_ELEMENTS(index)
IF nindex NE 0 THEN BEGIN
  outstr={tmem: 0., dec: double(0.), pidec: double(0.)} 
  de_c=DOUBLE(0.)
  de=1.986d-8/waves
  FOR i=0,nindex-1 DO BEGIN
    de_c=de_c+c_lambda(i)*de(i)
  ENDFOR
  pide_c=4d0*!pi/de_c
  outstr.tmem=t_mem(0)
  outstr.dec=de_c
  outstr.pidec=pide_c
  RETURN
ENDIF

;--------------------------------<
; Print out the C_lambda information for each ion
;
ind2=SORT(waves)
;
PRINT,''
PRINT,'        Lambda    T_mem   C_lambda  DE*C_lambda 4pi/DE*C_lambda'
PRINT,''
FOR i=0,n-1 DO BEGIN
  PRINT,FORMAT='(i3,f12.3,f8.2,e12.4,e12.4,e12.4)',i+1,waves(ind2(i)), $
	t_mem(ind2(i)),c_lambda(ind2(i)),$
        c_lambda(ind2(i))*1.986d-8/waves(ind2(i)),$
        4d0*!pi/(c_lambda(ind2(i))*1.986d-8/waves(ind2(i)))
ENDFOR
;--------------------------------<

;-----------------+
; Use displayed information?
;
ans=''
PRINT,''
READ,'Blend lines or plot C_lambda function? (e.g., 1+3+7 or p2)  ',ans
;
ans=STRCOMPRESS(ans,/REMOVE_ALL)         ; Tidy up ans
;-----------------+

;-----------0
; If there's a `p' it indicates that a plot is required. A plot of a blend is
; allowed
;
test=STR_SEP(ans,'p')
IF N_ELEMENTS(test) GT 1 THEN BEGIN
  ans=test(1)
  plot=1
ENDIF
;-----------0

;------------------------{
; Lines linked by `+' indicates a blend, and the C_lambda info for the
; blend is displayed.
;
indices=STR_SEP(ans,'+')
n=n_elements(indices)
IF n GT 1 THEN BEGIN
  indices=fix(indices)-1
  PRINT,''
  PRINT,'Chosen wavelengths:  ', $
	string(format='(6f10.3)',waves(ind2(indices)))

  de_c=0.
  FOR i=0,n-1 DO BEGIN
    de_c=de_c+1.986d-8/waves(ind2(indices(i))) * c_lambda(ind2(indices(i)))
  ENDFOR
  PRINT,''
  PRINT,' TOTAL(DE*C_lambda)   4pi/TOTAL(DE*C_lambda)'
  PRINT,''
  PRINT,FORMAT='(7x,e12.4,13x,e12.4)',de_c,4d0*!pi/de_c
  PRINT,''
ENDIF
;------------------------{

;----------------------<
; The following requires the emissivity array in EMISS to be 1 dimensional,
; i.e., no temp dependence.
;
IF KEYWORD_SET(plot) THEN BEGIN
  n=N_ELEMENTS(indices)
  plot_em=emiss(0).em-emiss(0).em
  FOR i=0,n-1 DO plot_em=emiss(ind(indices)).em+plot_em
  y2=spl_init(temp,plot_em)
  plot_emi=spl_interp(temp,plot_em,y2,tempi)
  PLOT,tempi,plot_emi*frac,charsiz=1.4, $
        xtitle='!3Log!d10!n ( Temperature [K] )',$
	ytitle='!3n!dj!n A!dji!n Fr(T)!3'
  goft=plot_emi*frac
ENDIF
;----------------------<

END
