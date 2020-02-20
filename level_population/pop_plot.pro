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
; NAME: POP_PLOT
;       
; PURPOSE:
;
;	To compute n_j A_ji / N_e for a selected transition(s) and plot it
;	against N_e. If it is insensitive to N_e, then the line(s) is 
;	suitable for emission measure analysis.
;
; CATEGORY:
;
;       Atomic data analysis
;
; EXPLANATION:
;
;	The routine calls EMISS_CALC to give values of DE*n_j*A_ji at the 
;	temperature TEMP and densities from 10^8 to 10^12. You are then 
;	asked to select which transition(s) you are interested in. (If 
;	more than one line is selected, the lines are blended.) 
;	DE*n_j*A_ji/N_e is then plotted against density.
;
;	If TEMP is not specified, then the temperature at which the 
;	ionisation fraction has its maximum is calculated. For the iron 
;	ions, the ion balance calcs of Arnaud & Raymond are used, 
;	otherwise Arnaud & Rothenflug are used. If TEMP is specified, 
;	and the value is less than 20, then it is assumed that the log 
;	of the temperature has been specified.
;
;	In emission measure work it is important to isolate lines for 
;	which DE*n_j*A_ji/N_e is insensitive to density. If only such lines 
;	are used, then the derived emission measure curve is independent 
;	of density.
;
; CALLING SEQUENCE:
;
;       POP_PLOT, IZ, ION, WRANGE=WRANGE, [TEMP=TEMP, /QUICK, DATA=DATA, $
;				DENS_RANGE=DENS_RANGE, DILUTE=DILUTE]
;
; EXAMPLES:
;
;	POP_PLOT, 26, 14, WRANGE=[330,360]
;
;	- 3 lines should appear in the widget. Selecting 334.17 should show 
;	a curve that falls off with density. Choosing 353.83 shows a curve 
;	that rises with density. By selecting a blend of the two lines, 
;	the curve will be insensitive to density, telling us that only a 
;	blend of 334.17 and 353.83 is suitable for emission measure work.
;
;	POP_PLOT, 8, 4, WRANGE=[550,560], TEMP=6.0, /QUICK, DENS_RANGE=[6,10]
;
;	- O IV is a member of the boron sequence, and so calculations take a 
;	lot longer. Giving the QUICK keyword speeds things up. The 
;	temperature is well away from the T_max of the ion
;
; INPUTS:
;
;	IZ	The atomic number of the ion
;
;	ION	The spectroscopic number of the ion (e.g., 12 = XII)
;
; OPTIONAL INPUTS:
;
;	DILUTE	Used to set radiative dilution factor. (Default: 0.0)
;
;	TEMP	The temperature at which calculations are required. Usually 
;		this will be the Tmax of the ion.
;
;	DENS_RANGE  The default density range is log Ne = 8 to 12. By 
;		    inputting two integers, a different range can be chosen.
;
;	WRANGE  Wavelength range from which lines are required. If not 
;		given, then the user is allowed to choose from the complete 
;		set of lines for the ion.
;
;       ABUND_FILE  The name of a CHIANTI abundance file. This is used for 
;               calculating the proton to electron ratio. Default is 
;               !abund_file.
;
;       IONEQ_FILE  The name of a CHIANTI ion balance file. This is used for 
;               calculating the proton to electron ratio and evaluating 
;               the T_max of the ion. Default is !ioneq_file.
;
;       RADTEMP The blackbody radiation field temperature (default 
;               6000 K).
;
;       RPHOT   Distance from the centre of the star in stellar radius 
;               units. I.e., RPHOT=1 corresponds to the star's surface. 
;               (Default is infinity, i.e., no photoexcitation.)
;
; OPTIONAL OUTPUTS:
;
;	DATA:	An array that contains the data that is plotted: data(*,0) 
;		contains
;		the densities, while data(*,1) contains the Y-axis values.
;
; KEYWORDS:
;
;       NOPROT  If set, then the default setting will be NOT to use 
;               proton rates. This can be changed within the routine.
;
;	QUICK:	The density range over which the calculations are done is 
;		8 to 12 in 0.2 intervals. This keyword forces the 
;		calculations to be done at 0.5 intervals.
;
; CALLS:
;
;	EMISS_CALC, EMISS_SELECT, READ_IONEQ
;
; HISTORY:
;
;	Ver.1, PRY 28-Jul-97.
;	Ver.2, PRY 23-Sep-97 - added DILUTE keyword for photo-excitation
;	Ver.3, PRY 30-May-98 - added DENS_RANGE
;	Ver.4, PRY 5-Sep-98  - added call to choose_ioneq
;	Ver.5, PRY 7-Apr-99  - tidied up, and introduced call to emiss_select
;       Ver.6, PRY 7-Dec-01  - modified for v.4 of CHIANTI
;
;       V. 7, 21-May-2002, Giulio Del Zanna (GDZ) 
;                   generalized directory concatenation to work for
;                   Unix, Windows  and VMS. 
;
;       V.8, 06-Aug-02 GDZ
;              added ABUND_FILE to the call to emiss_calc (was missing).
;              Changed the use of CHIANTI system variables. 
;
;
; VERSION     : 8,  06-Aug-02
;-

PRO POP_PLOT, IZ, ION, WRANGE=WRANGE, TEMP=TEMP, QUICK=QUICK, DATA=DATA, $
              RPHOT=RPHOT, RADTEMP=RADTEMP, DENS_RANGE=DENS_RANGE, $
              IONEQ_FILE=IONEQ_FILE, ABUND_FILE=ABUND_FILE, NOPROT=NOPROT


IF N_PARAMS() LT 2 THEN BEGIN
  PRINT,'Use: IDL> pop_plot, iz, ion, wrange= , [data= , temp= , /noprot,'
  PRINT,'                      dens_range= , /quick, radtemp= ,'
  print,'                      rphot= , abund_file= , ioneq_file= ]'
  RETURN
ENDIF

IF N_ELEMENTS(temp) EQ 0 THEN BEGIN
;
;------------------------------------+

;dir=concat_dir(!xuvtop,'ioneq')
IF n_elements(ioneq_file) NE 0 THEN ioneq_name=ioneq_file $
     ELSE ioneq_name= !ioneq_file

read_ioneq,ioneq_name,temp_all,ioneq,ref

f_all=ioneq(*,iz-1,ion-1)
ind=WHERE(f_all EQ max(f_all))
ind=ind(0)      ; in case ind has more than one element (e.g., Ar IX)
temp=temp_all(ind)
;------------------------------------+
;
ENDIF ELSE BEGIN
  IF temp GT 20. THEN temp=alog10(temp)
ENDELSE

PRINT,format='("Temperature: ",f7.1)',temp

IF N_ELEMENTS(dens_range) EQ 0 THEN dens_range=[8,12]

;----------------<
; The density range is set at [8:12]. The default is 0.2 dex intervals,
; but with the QUICK keyword 0.5 dex intervals can be used.
;
dens=FINDGEN( (dens_range(1)-dens_range(0))*5+1)/5. + dens_range(0)
IF KEYWORD_SET(QUICK) THEN dens=FINDGEN((dens_range(1)-dens_range(0))*2+1)/2. + dens_range(0)
;----------------<

IF NOT KEYWORD_SET(dilute) THEN dilute=0.

emiss=emiss_calc(iz,ion,temp=temp,dens=dens,/quiet,rphot=rphot, $
                ioneq_file=ioneq_file,ABUND_FILE=ABUND_FILE, $
                radtemp=radtemp, $
                noprot=noprot)

IF N_ELEMENTS(wrange) EQ 0 THEN wrange=0 ELSE wrange=FLOAT(wrange)

cal_emiss=emiss_select(emiss,wrange=wrange)

data_plot=cal_emiss/10.^(dens) * 10.^(20)

PLOT, dens, data_plot, th=2, charsiz=1.3, $
	xtit='Log!d10!n N!de!n',$
	ytit='!4D!3E n!dj!n A!dji!n x 10!u20!n / N!de!n'

data=fltarr(n_elements(data_plot),2)
data(*,0)=dens
data(*,1)=data_plot

END
