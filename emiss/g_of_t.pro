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
; NAME: G_OF_T()
;       
; PURPOSE:
;
;	To compute DE * F(T) * n_j * A_ji / N_e for selected emission lines. 
;	Option to also multiply by abundance.
;
; CATEGORY:
;
;       Atomic data analysis
;
; EXPLANATION:
;
;	The G-of-T function has a number of different definitions in the 
;	literature. In the most basic form it contains only the temperature 
;	dependent parts (i.e., 0.83*n_j*A_ji*F(T)/N_e), but often a Delta-E 
;	and Ab(X) are added as well. Here, the _default_ form is:
;
;	Delta-E * 0.83 * n_j * A_ji * F(T) / N_e
;
;	By using the NO_DE keyword, the Delta-E can be omitted, while the 
;	ABUND keyword allows the abundance to be added.
;
;	The function that is output is tabulated over 4.0 <= logT <= 8.0 
;	in 0.1 dex intervals. If you want the function tabulated over 
;	smaller intervals, run the ION_INTERP routine afterwards.
;
; CALLING SEQUENCE:
;
;
; EXAMPLES:
;
;	RESULT=G_OF_T(26,13)
;
;	RESULT=G_OF_T(26,13,DENS=7)
;
;	RESULT=G_OF_T(26,13,GOFTFILE='my_gofts.dat')
;
;       RESULT=G_OF_T(26,13,/ABUND)
;
;       RESULT=G_OF_T(26,13,ABUND_FILE=ABUND_FILE, IONEQ_FILE=IONEQ_FILE)
;
; INPUTS:
;
;	IZ:	The atomic number of the ion (e.g., 26 = Fe)
;
;	ION:	The spectroscopic number of the ion (e.g., 12 = XII)
;
; OPTIONAL INPUTS:
;
;	DENS:	The logarithm (to base 10) of the density at which the 
;		emissivities are calculated (default=10.)
;
;	WRANGE: Wavelength range from which lines are required. If not 
;		specified, then the 10 strongest lines are displayed.
;
;	PATH:	If specified, the routine will look for the atomic data in 
;		the PATH directory, rather than in the CHIANTI database
;
;	GOFTFILE:	By specifying GOFTFILE as a filename, the G-of-T 
;			function can be stored in this file. It is stored 
;		in the form a structure (called goft_list) with the following 
;		labels:
;
;	  goft_list.label: user-specified string, e.g., 'Si XII  520.7'
;	  goft_list.func:	 fltarr(41), the G-of-T function
;
;	If the same GOFTFILE is specified for another ion, then the 
;	G-of-T function is added to the existing structure. The GOFTFILE 
;	option only works when the ABUND keyword is set. The GOFTFILE is 
;	meant to be read by another routine called GOFT_PLOTTER.
;
;	INDEX:	Allows the direct specification of indices within the 
;		emiss structure. This allows the g_of_t routine to be 
;		run without the menu widget appearing. If the /quiet 
;		keyword is also used, then the routine will run "silently".
;
;	RADTEMP	Specify background radiation temperature (default: 6000 K)
;
;	RPHOT   Distance from the centre of the star in stellar radius units.
;               I.e., RPHOT=1 corresponds to the star's surface. (Default is
;               infinity, i.e., no photoexcitation.)
;
;       IONEQ_FILE  Directly specify the name of the ion balance file 
;               (including directory path). If not set, then a widget will 
;               pop up allowing you to select a file.
;
;       ABUND_FILE  Directly specify the name of the abundance file 
;               (including directory path). One can also use /ABUND_FILE 
;               to include the abundances in the G(T) function, but allow 
;               the abundance file to be picked through a widget.
;
; KEYWORDS:
;
;       NOPROT  If set, then the default setting will be NOT to use 
;               proton rates. This can be changed within the routine.
;
;	NO_DE:	If set, then the output function will *not* contain the 
;		Delta-E. Be careful using this if you are using blends 
;		(as Delta-E is different for different wavelengths).
;
;	QUIET	If set, then don't plot the G(T) function or print out 
;               information to the screen.
;
; CALLS:
;
;	EMISS_CALC, READ_IONEQ, READ_ABUND, EMISS_SELECT, CH_GET_FILE
;
; RESTRICTIONS:
;
;
; HISTORY:
;
;	Ver.1, PRY 28-Jul-97.
;	Ver.2, PRY 22-Jun-98, added CHOOSE keyword and removed RAY
;	Ver.3, PRY 4-Aug-98, cosmetic changes following comments of Giulio 
;			Del Zanna
;	Ver.4, PRY 5-Sep-98, added call to choose_ioneq
;	Ver.5, PRY 23-Sep-98, added GOFTFILE input
;	Ver.6, PRY 3-Mar-99, now calls EMISS_SELECT
;       Ver.7, PRY 6-Nov-00, removed the /CHOOSE keyword; also changed 
;                       PICKFILE call to DIALOG_PICKFILE and removed call 
;                       to the CHIANTI routine ADD\_SUBDIR
;       Ver.8, PRY 18-Oct-01, adjusted for proton rates, and 
;                       photoexcitation.
;       Ver.9, PRY 9-Dec-01, completed modifications for v.4 of CHIANTI.
;
;       V.  10, 21-May-2002, Giulio Del Zanna (GDZ) 
;                   generalized directory concatenation to work for
;                   Unix, Windows  and VMS.
;       Ver. 11, 9-Feb-2005, Peter Young
;                   changed keyword_set(abund) to keyword_set(abund_file)
;
; VERSION     :   11, 9-Feb-2005
;
;-

FUNCTION G_OF_T, IZ, ION, WRANGE=WRANGE, DENS=DENS, $
                 NO_DE=NO_DE, PATH=PATH, GOFTFILE=GOFTFILE, $
                 INDEX=INDEX, QUIET=QUIET, $
                 RADTEMP=RADTEMP, RPHOT=RPHOT, NOPROT=NOPROT, $
                 IONEQ_FILE=IONEQ_FILE, ABUND_FILE=ABUND_FILE


IF N_PARAMS() LT 2 THEN BEGIN
  PRINT,'Use:  IDL> result=g_of_t( ion, iz [,wrange=, dens=, path=,'
  PRINT,'                            goftfile=, /no_de, /noprot,'
  PRINT,'                            /quiet, index=, radtemp=, rphot=,'
  print,'                            ioneq_file=, abund_file= ] )'
  RETURN,0.
ENDIF


IF NOT KEYWORD_SET(dens) THEN dens=10.
IF dens GT 100. THEN BEGIN
  PRINT,'** Please specify the logarithm of the density **'
  RETURN,0.
ENDIF
IF NOT KEYWORD_SET(quiet) THEN PRINT,FORMAT='("Log10 density:   ",f5.1)',dens

;---------------------------------------<0>
IF n_elements(ioneq_file) NE 0 THEN BEGIN
  ioneq_name=ioneq_file
  print,'Using file: ',ioneq_name
ENDIF ELSE BEGIN

  dir=concat_dir(!xuvtop,'ioneq')
  ioneq_name=ch_get_file(path=dir,filter='*.ioneq', $
                         title='Select Ionization Equilibrium File')
ENDELSE

read_ioneq,ioneq_name,temp_all,ioneq,ref

;IF NOT KEYWORD_SET(quiet) THEN PRINT,$
;        'Ion balance file:   ',ioneq_name
;
; Pick out non-zero elements of ioneq
;
ioneq=REFORM(ioneq(*,iz-1,ion-1))
ind=WHERE(ioneq NE 0.)
;---------------------------------------<0>


;----------------------------------[]
IF keyword_set(abund_file) THEN BEGIN
  res=findfile(string(abund_file))
  IF res[0] EQ '' THEN BEGIN
    dir=concat_dir(!xuvtop,'abundance')
    abund_name=ch_get_file(path=dir,filter='*.abund', $
                           title='Select Abundance File')
  ENDIF ELSE BEGIN
    abund_name=abund_file
    print,'Using file: ',abund_name
  ENDELSE
  read_abund,abund_name,ab,ref
  ab=ab(iz-1)
;  IF NOT KEYWORD_SET(quiet) THEN PRINT,$
;        'Abundance file:     ',abund_name
ENDIF ELSE BEGIN
  ab=1.
ENDELSE
;----------------------------------[]

emiss=emiss_calc(iz,ion,temp=temp_all(ind),dens=dens,no_de=no_de, $
                 path=path, /quiet, rphot=rphot, radtemp=radtemp, $
                 noprot=noprot, abund_file=abund_name, $
                 ioneq_file=ioneq_name)

IF N_ELEMENTS(index) EQ 0 THEN BEGIN
  IF KEYWORD_SET(wrange) THEN BEGIN
    cal_emiss=emiss_select(emiss,wrange=wrange,sel_ind=sel_ind)
  ENDIF ELSE BEGIN
    n_em=N_ELEMENTS(ind)
    index=REVERSE( SORT (emiss.em(FIX(n_em/2))) )
    index=index(0:9)
    cal_emiss=emiss_select(emiss,index,sel_ind=sel_ind)
  ENDELSE
ENDIF ELSE BEGIN
  sel_ind=index
  cal_emiss=emiss(index).em
  IF N_ELEMENTS(index) GT 1 THEN cal_emiss=TOTAL(cal_emiss,2) $
               ELSE cal_emiss=REFORM(emiss(index).em)
ENDELSE

chosen_wavels=emiss(sel_ind).lambda
IF NOT KEYWORD_SET(quiet) THEN PRINT,'Chosen wavelengths: ',chosen_wavels

;--------------[]
; Calculate the G(T) function over the non-zero temp range. Note: dividing
; by the density. Once calculated, work out G(T) over entire temp range.
;
func_short=0.83*cal_emiss*ioneq(ind)/10.^(dens) * ab
;
func=ioneq-ioneq
func(ind)=func_short
;--------------[]

IF NOT KEYWORD_SET(quiet) THEN BEGIN
;---------------------------------------<O>
; Create the title for the plot
;
n=N_ELEMENTS(chosen_wavels)
label=''
FOR i=0,n-1 DO BEGIN
  label=label+'+'+strtrim(string(format='(f10.3)',chosen_wavels(i)),2)
ENDFOR
len=strlen(label)
label=strmid(label,1,len-1)
;
zion2spectroscopic,iz,ion,name
title=name+'  '+label
;
ytit='n!dj!n A!dji!n F(T) / N!de!n'
IF KEYWORD_SET(abund_file) THEN ytit='Ab(X) '+ytit
IF NOT KEYWORD_SET(no_de) THEN ytit='!4D!3E '+ytit
ytit='0.83 '+ytit
;
plot,temp_all,func,charsiz=1.3,xmarg=[15,3], $
    xtit='Log!d10!n ( Temperature [K] )!3', $
    ytit=ytit,title=title,xticklen=-0.015
;---------------------------------------<O>
ENDIF


;--------------------------------------[I]
IF (N_ELEMENTS(goftfile) NE 0) AND KEYWORD_SET(abund_file) THEN BEGIN
  ans=''
  READ,'Add to the G(T) list? (y/n) ',ans
  IF ans EQ 'y' THEN BEGIN
    result=FINDFILE(EXPAND_PATH(goftfile))
    IF result(0) EQ '' THEN BEGIN
      PRINT,'The G-of-T file does not exist, so one will be created...'
      str={label:'', func:FLTARR(41)}
      goft_list=REPLICATE(str,100)
    ENDIF ELSE BEGIN
      RESTORE,EXPAND_PATH(goftfile)
    ENDELSE
    ;
    ind=WHERE(goft_list.label NE '')
    n=N_ELEMENTS(ind)
    ;
    IF ind(0) EQ -1 THEN BEGIN
      n=0 
    ENDIF ELSE BEGIN
      PRINT,''
      PRINT,'Current list of lines:'
      FOR i=0,n-1 DO PRINT,goft_list(ind(i)).label
      PRINT,''
    ENDELSE
    ;
    READ,'Give a name for the label: ',ans
    goft_list(n).label=ans
    goft_list(n).func=func
    SAVE,file=EXPAND_PATH(goftfile),goft_list
    PRINT,'New entry has been added, and structure has been saved'
  ENDIF
ENDIF
;--------------------------------------[I]

RETURN,func

END
