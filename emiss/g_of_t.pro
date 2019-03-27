;+
; NAME: G_OF_T()
;       
; PURPOSE:
;	To compute the contribution function, G(T), for selected
;	emission lines. See the section 'OUTPUT' for the defintion of
;	G(T) used in this code.
;
;	Note that another CHIANTI routine, gofnt.pro, also calculates
;	the G(T) function and we generally recommend users to use
;	gofnt.pro. 
;
; CATEGORY:
;       CHIANTI; contribution function.
;
; CALLING SEQUENCE:
;       Result = G_of_T ( IZ )
;
; INPUTS:
;	IZ:	Either, the atomic number of the ion (e.g., 26 = Fe),
;	        or, the name of the ion in CHIANTI format (e.g.,
;	        'o_6'= O VI). We recommend the former to be
;	        consistent with other CHIANTI routines.
;
;	ION:	The spectroscopic number of the ion (e.g., 12 =
;        	XII). Not needed if IZ was specified as an ion string.
;
; OPTIONAL INPUTS:
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
; OUTPUTS:
;       The G(T) function for the selected emission line. The default
;       definition is:
;
;       G(T) = DE * F(T) * n_j * A_ji N_H/ N_e^2
;
;       where DE is the energy for the transition (erg), F(T) is the
;       ionization fraction, n_j is the population of the upper level
;       of the transition relative to the ion as a whole, A_ji is the
;       radiative decay rate, N_H is hydrogen number density, and N_e
;       is the electron number density. The units are erg cm^-3 s^-1.
;
;       If the input ABUND_FILE is given, then the function will be
;       multiplied by the element abundance.
;
;       One major difference with the output from gofnt.pro is that
;       the hydrogen-to-electron (N_H/N_e) is included for g_of_t, but
;       not for gofnt. Note that this ratio is calculated with
;       proton_dens.pro. 
;
; OPTIONAL OUTPUTS:
;       LogT:   The log temperature (units: K) array for which the
;               output function is tabulated.
;
; CALLS:
;	EMISS_CALC, READ_IONEQ, READ_ABUND, EMISS_SELECT, CH_GET_FILE,
;	PROTON_DENS. 
;
; EXAMPLES:
;       Note that the ion can be specified as 'fe_13' or 26,13.
;
;	RESULT=G_OF_T('fe_13')
;	RESULT=G_OF_T(26,13,DENS=7)
;	RESULT=G_OF_T('fe_13',GOFTFILE='my_gofts.dat')
;       RESULT=G_OF_T(26,13,/ABUND)
;       RESULT=G_OF_T('fe_13',ABUND_FILE=ABUND_FILE, IONEQ_FILE=IONEQ_FILE)
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
;                   changed keyword_set(abund) to
;                   keyword_set(abund_file)
;       Ver.12, 23-Jan-2018, Peter Young
;          Added optional output for temperature (logt); now calculate NH/Ne
;          ratio using proton_dens rather than use 0.83; allow the ion
;          name to be specified as an input; updated header.
;
; VERSION     :   12, 23-Jan-2018
;
;-

FUNCTION G_OF_T, IZ, ION, WRANGE=WRANGE, DENS=DENS, $
                 NO_DE=NO_DE, PATH=PATH, GOFTFILE=GOFTFILE, $
                 INDEX=INDEX, QUIET=QUIET, $
                 RADTEMP=RADTEMP, RPHOT=RPHOT, NOPROT=NOPROT, $
                 IONEQ_FILE=IONEQ_FILE, ABUND_FILE=ABUND_FILE, LogT=LogT


IF N_PARAMS() LT 1 THEN BEGIN
  PRINT,'Use:  IDL> result=g_of_t( ion_name [,wrange=, dens=, path=,'
  PRINT,'                            goftfile=, /no_de, /noprot,'
  PRINT,'                            /quiet, index=, radtemp=, rphot=,'
  print,'                            ioneq_file=, abund_file=, logt=logt ] )'
  print,''
  print,"  The ion can be specified either with 'o_6' or 8,6."
  RETURN,0.
ENDIF

;
; PRY, 23-Jan-2018
; If the 1st parameter is a string, then it is interpreted as the
; CHIANTI ion format, e.g., 'o_6'.
;
IF datatype(iz) EQ 'STR' THEN BEGIN
  iz_save=iz
  convertname,iz_save,iz,ion
ENDIF


IF NOT KEYWORD_SET(dens) THEN dens=10.
IF dens GT 100. THEN BEGIN
  PRINT,'** Please specify the logarithm of the density **'
  RETURN,0.
ENDIF
IF NOT KEYWORD_SET(quiet) THEN PRINT,FORMAT='("Log10 density:   ",f5.1)',dens

;---------------------------------------<0>
; PRY, 23-Jan-2018
; Changed this to automatically select !ioneq_file without widget. 
;
ioneq_name=!ioneq_file
IF n_elements(ioneq_file) NE 0 THEN BEGIN
  ioneq_name=ioneq_file
  chck=file_search(ioneq_file,count=count)
  IF count GT 0 THEN ioneq_name=ioneq_file
ENDIF 
read_ioneq,ioneq_name,temp_all,ioneq,ref
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

;
; PRY, 23-Jan-2018
; Get hydrogen-to-electron ratio.
; Note that if abund_file has not been set, then the CHIANTI default
; will be selected by proton_dens. 
;
nh_NE=proton_dens(temp_all[ind],ioneq_file=ioneq_file,abund_file=abund_file)


;--------------[]
; Calculate the G(T) function over the non-zero temp range. Note: dividing
; by the density. Once calculated, work out G(T) over entire temp range.
;
func_short=nh_ne*cal_emiss*ioneq(ind)/10.^(dens) * ab
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
ytit='N!dH!n/N!de!n '+ytit
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

;
; This restores the input value of IZ if it was specified as a string.
;
IF n_elements(iz_save) NE 0 THEN iz=iz_save

logt=temp_all

RETURN,func

END
