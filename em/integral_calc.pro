
PRO integral_calc, ionname, spect_num, ioneq_file=ioneq_file, quick=quick, $
                   wrange=wrange, index=index, outstr=outstr, all=all, $
                   emiss=emiss, choose=choose, abund_file=abund_file, $
                   radtemp=radtemp, rphot=rphot, quiet=quiet, noprot=noprot, $
                   path=path, dens=dens

;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. 
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
;       INTEGRAL_CALC, ION_NAME, [WRANGE=WRANGE, INDEX=INDEX, DENS=DENS]
;
; EXAMPLES:
;
;	INTEGRAL_CALC, 'fe_13'
;
;	INTEGRAL_CALC, 'si_10', WRANGE=[250,270], DENS=10.
;
; INPUTS:
;
;	IONNAME:  The name of the ion in CHIANTI format. E.g., 'fe_13'
;                 for Fe XIII. Alternatively, the atomic number of the
;                 ion can be specified as an integer, but then the
;                 spectroscopic number should be specified through the
;                 ION keyword. 
;
; OPTIONAL INPUTS:
;
;	ION:	The spectroscopic number of the ion (e.g., 12 =
;               XII). (This is only used if ionname was specified as
;               an integer.
;
;	DENS:	The logarithm (base 10) of the density at which the emissivities are calculated 
;		(default=10.)
;
;	WRANGE: Wavelength range for which C_lambda functions are 
;		calculated. If not given, then the 10 strongest lines 
;		are printed.
;
;	INDEX:	Particular elements in the emissivity structure can be 
;		selected with INDEX. This allows integral_calc to be run 
;		'silently'. The output is contained in the OUTSTR structure. 
;		If index is given as, e.g., [7,8], then the C_lambda 
;		functions for these two lines are summed and output.
;
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
;	EMISS_CALC, CH_GET_FILE, READ_IONEQ, GET_IEQ, PROTON_DENS,
;       READ_ABUND 
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
;              dT intervals following problems identified by Luca
;              Teriaca.
;
;       V.9, 9-Nov-2012, Peter Young
;              I forced the interpolation of the ionization fraction
;              to be on intervals of 0.01 dex, and to restrict it to
;              between the min and max of the ion fraction. The output
;              from integral_calc will not be affected.
;
;       V.10, 25-Mar-2013, Peter Young
;              This routine is now very slow for the large iron ion
;              models so a new keyword /QUICK has been implemented to
;              speed up calculations. In addition, the emissivity
;              structure can be output, so if the routine is called
;              again for the same ion the structure can be re-used
;              rather than calculated again.
;
;       V.11, 27-Mar-2013, Peter Young
;              Modified to ensure compatibility with IDL versions 7
;              and earlier. I now correctly compute the H-to-e ratio
;              rather than just assume 0.83.
;
;       V.12, 07-Dec-2020, Peter Young
;              Added ioneq_file and abund_file inputs to proton_dens,
;              allowing me to remove the common block.
;-



IF N_PARAMS() LT 1 THEN BEGIN
  PRINT,'Use: IDL> integral_calc, ion_name [, wrange= ,/choose, /noprot, dens=, index= '
  PRINT,'                          radtemp=, rphot= , outstr= , path= , ioneq_file= , '
  PRINT,'                          abund_file= ]'
  print,'E.g.,'
  print,"   IDL> integral_calc,'fe_13'  "
  print,"               - for 10 strongest lines"
  print,"   IDL> integral_calc,'fe_13',wrange=[200,205] "
  print,"               - for 5 strongest lines in wavelength range"
  print,"   IDL> integral_calc,'fe_13',wrange=[200,205],/all "
  print,"               - for all lines in wavelength range"
  print,"   IDL> integral_calc,'fe_13',index=1000,outstr=outstr"
  print,"               - send output for transition with emiss_calc index of 1000 to a structure"
  RETURN
ENDIF


IF keyword_set(choose) THEN BEGIN
  print,'%INTEGRAL_CALC: this keyword is now obsolete. Please use the input IONEQ_FILE to specify'
  print,'                 an alternative ion balance file.'
  return
ENDIF 


print_mess=0
CASE 1 OF 
  datatype(ionname) EQ 'STR': convertname,ionname,iz,ion
  datatype(ionname) EQ 'INT' OR datatype(ionname) EQ 'LON': BEGIN
    iz=ionname
    IF n_params() LT 2 THEN print_mess=1 ELSE ion=spect_num
  END 
  ELSE: print_mess=1
ENDCASE

IF print_mess EQ 1 THEN BEGIN
  print,''
  print,'Please check your inputs. The routine should be called as:'
  print,"     IDL> integral_calc,'fe_13'"
  print,"or   IDL> integral_calc,26,13"
  print,''
  return
ENDIF

IF n_elements(ioneq_file) EQ 0 THEN ioneq_file=!ioneq_file
IF n_elements(dens) EQ 0 THEN dens=10.0
;
IF NOT keyword_set(quiet) THEN BEGIN
  print,''
  print,'Ion balance file: ',trim(ioneq_file)
  print,'Log (density/cm^-3): ',trim(string(format='(f10.1)',dens))
ENDIF 


IF n_elements(abund_file) EQ 0 THEN abund_file=!abund_file
read_abund,abund_file,abund,abund_ref


;
; Temperature arrays used in this routine:
;
;  ion_t   - (log) temperatures over which ion fraction defined
;  em_t    - (log) temperatures over which emissivity defined
;  tempi   - ion_t interpolated onto 0.01 dex intervals (used for integration).
;


;
; Obtain ionization fraction for selected ion.
;
read_ioneq,ioneq_file,temp_all,ioneq,ioneq_ref
;
ion_f=ioneq[*,iz-1,ion-1]
k=where(ion_f NE 0.)
ion_t=temp_all[k]       ; note that ion_t is log(T)
ion_f=ion_f[k]

IF keyword_set(quick) THEN BEGIN
  k=where(ion_f GE max(ion_f)/1000.)
  ion_f=ion_f[k]
  ion_t=ion_t[k]
ENDIF


; 
; PRY, 9-Nov-2012
; Interpolate ion fraction onto a temperature scale with dlogT=0.01
;
num_t=round((max(ion_t)-min(ion_t))/0.01)+1
tempi=dindgen(num_t)*0.01+min(ion_t)

; calculate ion fraction on this scale
;
frac=get_ieq(10.^tempi,iz,ion,ioneq_logt=temp_all,ioneq_frac=ioneq)


;
; For large ion models, the number of temperatures used for the
; emissivity calculation can really slow down the routine. To
; alleviate this, the default option is to calculate the emissivities
; at 0.10 dex intervals in logT (the ion fraction is usually tabulated
; at 0.05 dex intervals). For further speed improvements, the /QUICK
; keyword means that only n_em temperatures are used. I find for the
; coronal iron ions that there's a minimal loss of accuracy (1%
; or less) when doing this. For other ions, e.g., H I or Li-like ions
; the loss of accuracy can be higher.
;
IF keyword_set(quick) THEN BEGIN 
  n_em=6
  em_dt=(max(ion_t)-min(ion_t))/float(n_em-1)
  em_t=findgen(n_em)*em_dt+min(ion_t)
ENDIF ELSE BEGIN 
  em_t=findgen(151)/10.
  t_min=floor(min(ion_t)*10.)/10.
  t_max=ceil(max(ion_t)*10.)/10.
  k=where(em_t GE t_min AND em_t LE t_max,n_em)
  em_t=em_t[k]
ENDELSE 


IF n_tags(emiss) NE 0 THEN BEGIN
  zion2spectroscopic,iz,ion,name
  IF name NE emiss[0].ion_name THEN BEGIN
    print,'%INTEGRAL_CALC: the input EMISS is not compatible with the specified ion. Returning...'
    return
  ENDIF 
  IF n_em NE n_elements(emiss[0].em) THEN BEGIN
    print,'%INTEGRAL_CALC: the input EMISS is not compatible with other input options. Returning...'
    return
  ENDIF 
ENDIF ELSE BEGIN 
  emiss=emiss_calc(iz,ion,dens=dens,temp=em_t,/no_de,path=path,/quiet, $
                   noprot=noprot, ioneq_file=ioneq_file, $
                   abund_file=abund_file, rphot=rphot, radtemp=radtemp)
ENDELSE


;
; The code below loads up the index (of EMISS) for the lines that are
; going to be printed. There are 3 cases depending if INDEX, WRANGE or
; neither are input. The index is called IND.
;
swtch=0
CASE 1 OF
 ;
 ; Use line(s) specified by INDEX. 
 ;
  n_elements(index) NE 0: BEGIN
    ind=index
    n_ind=n_elements(ind)
    IF NOT keyword_set(quiet) THEN BEGIN
      print,'Lines selected by INDEX are:'
      FOR i=0,n_ind-1 DO BEGIN 
        trans=trim(emiss[ind[i]].lvl1_desc)+' - '+trim(emiss[ind[i]].lvl2_desc)
        trans=strpad(trans,60,fill=' ',/after)
        print,format='(i7,f15.3,2x,a60)',ind[i],emiss[ind[i]].lambda,trans
      ENDFOR 
    ENDIF 
  END 
 ;
 ; All lines in specified wavelength range
 ;
  n_elements(wrange) NE 0: BEGIN
    swtch=1
    ind=where(emiss.lambda GE wrange[0] AND emiss.lambda LE wrange[1],n_ind)
    IF n_ind EQ 0 THEN BEGIN
      print,'%INTEGRAL_CALC: no transitions were found in the specified wavelength range. Returning...'
      return 
    ENDIF 
   ;
   ; The following retains only the five strongest lines in the
   ; wavelength range, unles /all was set.
   ;
    IF NOT keyword_set(all) THEN BEGIN
      k=reverse(sort(emiss[ind].em[n_em/2]))
      n_ind=min([n_elements(emiss[ind]),5])
      ind=ind[k[0:n_ind-1]]
    ENDIF 
  END
 ;
 ; Select the "top 10" lines
 ;
  ELSE: BEGIN
    swtch=1
    ind=reverse(sort(emiss.em(fix(n_em/2))))
    n_ind=min([n_elements(emiss),10])
    ind=ind[0:n_ind-1]
  END
ENDCASE 

;
; compute H to electron ratio
;
h_rat=proton_dens(ion_t,/hydrogen, abund_file=abund_file, ioneq_file=ioneq_file)
y2=spl_init(ion_t,h_rat)
h_to_e_ratio=spl_interp(ion_t,h_rat,y2,tempi)


t_mem=fltarr(n_ind)
c_lambda=dblarr(n_ind)
waves=emiss[ind].lambda
warning=bytarr(n_ind)

;-------------------0
; Work out C_lambda for each line in the specified wavelength range
;
FOR i=0,n_ind-1 DO BEGIN
  em=emiss[ind[i]].em
 ;
  y2=spl_init(em_t,em)
  emi=spl_interp(em_t,em,y2,tempi)
 ;
 ; dT = alog(10) * T * d(logT), and d(logT)=0.02
 ;
  integrand=int_tabulated(tempi,h_to_e_ratio*emi*frac*10.^tempi*alog(10d0),/double)
  get_imax=MAX(emi*frac,imax)
  t_mem(i)=tempi(imax)
  c_lambda(i)=integrand/(0.7046*10.^(tempi(imax))*10.^(dens))
 ;
  IF keyword_set(quick) THEN BEGIN
    func=emi*frac
    nf=n_elements(func)
    IF func[0]/max(func) GE 0.01 OR func[nf-1]/max(func) GE 0.01 THEN warning[i]=1b
  ENDIF 
ENDFOR
;-------------------0

;
; Return numbers 
;
IF swtch EQ 0 THEN BEGIN
  outstr={tmem: 0., dec: double(0.), pidec: double(0.)} 
  de_c=DOUBLE(0.)
  de=1.986d-8/waves
  FOR i=0,n_ind-1 DO BEGIN
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
PRINT,'      Emiss ind    Lambda    T_mem   C_lambda  DE*C_lambda 4pi/DE*C_lambda'
PRINT,''
FOR i=0,n_ind-1 DO BEGIN
  IF warning[i] EQ 1 THEN wstr='**' ELSE wstr=''
  PRINT,FORMAT='(i3,i11,f12.3,f8.2,e12.4,e12.4,e12.4,a5)',i+1, ind[ind2[i]], $
        waves(ind2(i)), $
	t_mem(ind2(i)),c_lambda(ind2(i)),$
        c_lambda(ind2(i))*1.986d-8/waves(ind2(i)),$
        4d0*!pi/(c_lambda(ind2(i))*1.986d-8/waves(ind2(i))),wstr
ENDFOR
;--------------------------------<

IF total(warning) GT 0 THEN BEGIN
  print,''
  print,'**WARNING** lines marked with "**" are inaccurate. Please do not use /QUICK.'
ENDIF 

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
  y2=spl_init(em_t,plot_em)
  plot_emi=spl_interp(em_t,plot_em,y2,tempi)
  PLOT,tempi,plot_emi*frac,charsiz=1.4, $
        xtitle='!3Log!d10!n ( Temperature [K] )',$
	ytitle='!3n!dj!n A!dji!n Fr(T)!3',/xsty
  goft=plot_emi*frac
ENDIF
;----------------------<


END
