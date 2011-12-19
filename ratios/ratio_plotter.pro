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
; NAME: RATIO_PLOTTER
;       
; PURPOSE:
;	A widget-based routine to allow the analysis of density or 
;       temperature sensitive ratios.
;
; EXPLANATION:
;
;	RATIO_PLOTTER is designed to study temperature and density 
;       sensitive line ratios. For temperature ratios, just use the 
;       keyword /temp; for density ratios just use /dens.
;
;       The routine allows line ratios relative to a 
;	particular line (called the 'denominator line') to be viewed. One 
;	can either choose a particular numerator, or one can simply view 
;	all ratios within a specified wavelength range (see the button widget 
;	just above the main plot window).
;
;	In the former case, the numerator line can be chosen either from the 
;	strongest lines, or from all lines within a specified wavelength 
;       range.
;
;	A wavelength range is selected via the sliders at the lower left of 
;	the main widget.
;
;       A default set of emissivity parameters is used initially to create 
;       the first plot that you see. To change the parameters, click on the 
;       'CHANGE PARAMETERS' button near the top of the GUI. This activates 
;       a new set of widgets allowing the user to change whether energy or 
;       photon units are needed, whether proton rates are included, the 
;       temperature/density range, etc. Note that while the parameters are 
;       being changed, the user can not modify the line selection or plot 
;       parameters.
;
;	When a particular numerator line has been selected, observed 
;	intensities (and error bars) for both the denominator and numerator 
;	can be input, and the derived density (plus error bars) can be seen 
;	on the plot (click on 'PLOT ERROR BARS'), or in a text widget (see 
;       the 'SHOW DERIVED DENSITIES' button).
;
;	The accuracy of the derived density or temperature is limited by the 
;	intervals in the emissivity calculation. The intervals can be 
;       changed on the widget (there are 4 choices given) and if a new 
;       interval is selected, new emissivities should be calculated. 
;       Using the smallest intervals will give the highest accuracy for 
;       the computed density or temperature.
;
;
; CALLING SEQUENCE:
;
;       RATIO_PLOTTER, IZ, ION [, /TEMP, /DENS, EM, PATH=, /NOPROT, $
;                      IONEQ_FILE=, ABUND_FILE= ]
;
; EXAMPLES:
;
;	RATIO_PLOTTER, 26, 13, /TEMP        ; Fe XIII
;
;	RATIO_PLOTTER, 10, 6, /DENS	    ; Ne VI
;
;	RATIO_PLOTTER, 26, 13, /DENS, PATH='/home/new_chianti_data'
;
; INPUTS:
;
;	IZ:	The atomic number of the ion
;	ION:	The spectroscopic number of the ion (e.g., 12 = XII)
;
; OPT. INPUTS:
;
;	EM:	Save the displayed emissivities to structure EM. This 
;		structure is simply the structure EMISS_SELECT used in 
;		the routine, with some extra tags. This structure has the 
;               tags:
;
;               .lines.emiss  Emissivities of line at requested densities
;               .lines.wavel  Wavelength(s) of line (blend)
;               .dens         Log10 electron density/densities
;               .temp         Log10 electron temperature(s)
;               .rphot        Photoexcitation radius
;               .radt         Radiation temperature
;               .proton       String indicating whether proton rates are 
;                             included.
;               .date         The date and time at which structure was 
;                             created.
;               .version      The version of CHIANTI that created the 
;                             structure.
;
;	PATH:	Data in the CHIANTI format that is not in the CHIANTI 
;		database can be read by specifying the directory in which 
;		it lies through PATH.
;
;       ABUND_FILE  The name of a CHIANTI abundance file. This is used for 
;               calculating the proton to electron ratio. Default is 
;               !abund_file.
;
;       IONEQ_FILE  The name of a CHIANTI ion balance file. This is used for 
;               calculating the proton to electron ratio and evaluating 
;               the T_max of the ion. Default is !ioneq_file.
;
; KEYWORDS:
;
;       DENSITY If set then ratios are plotted as a function of 
;               density rather than temperature.
;
;       TEMPERATURE  If set then the ratios are plotted as a function of 
;               temperature rather than density.
;
;       NOPROT  If set, then the default setting will be NOT to use 
;               proton rates. This can be changed within the routine.
;
; PROGRAMMING NOTES:
;
;	RATIO_PLOTTER is actually a collection of several routines....
;
;	DENS_FINDER	- converts observed ratios into densities
;	INDEX_EXTRACTOR	- a device for extracting indices out of arrays
;       WAVEL_PLOT      - plots the little window in bottom left of GUI
;	DENS_PLOT	- updates the plot window
;	DENS_MAIN_EVENT	- responds to widget actions
;	SAMPLE_WID 	- sets up the widgets
;	RATIO_PLOTTER	- loads up various initial parameters
;
;	The emissivities of all the ion's lines are stored in the 
;	structure "emiss", while those of the selected lines are stored 
;	in "emiss_sel". This latter can be output into IDL by giving the 
;	"EM" optional input on the command line.
;
;	The routine is set up to allow more than one numerator to be 
;	displayed at once, but I've never actually found a need to 
;	implement this yet. (Note that you'll see "Numerator 1" on 
;	the widget.)
;
; CALLS:
;
;	EMISS_CALC, READ_IONEQ, EMISS_SELECT, ZION2NAME,
;	ZION2FILENAME, ZION2SPECTROSCOPIC, PS, PSCLOSE, IS_NUMBER,
;       CHIANTI_FONT, TRIM
;
; INTERNAL PROCEDURES AND FUNCTIONS
;
;       DENS_FINDER, INDEX_EXTRACTOR, DENS_PLOT, DENS_MAIN_EVENT,
;       SAMPLE_WID, WAVEL_PLOT
;
; COMMON BLOCKS
;
;       EMISS_DATA, SELECT, SLIDERS, PLOTTING, RAD_DATA, PROTON_RATES,
;       EXTRA, ELVLC, FILES, TEMP
;
; HISTORY:
;
;	Ver.1: PRY, 15-SEP-97.
;	Ver.2: PRY, 6-JUL-98,  added EM and PATH
;	Ver.3: PRY, 5-SEP-98,  added call to choose_ioneq
;	Ver.4: PRY, 11-JAN-99, added PROTON keyword and widget to allow 
;			the addition of proton rates
;	Ver.5: PRY, 3-FEB-99,  added a title to the widget, and the name 
;			of the ion to the ratio plot title
;	Ver.6: PRY, 25-MAR-99, corrected bug for wavelengths > 1000 
;			angstroms
;	Ver.7: PRY, 1-JUN-99,  routine now works in CDS environment 
;			without use_dere_chianti
;	Ver.8: PRY, 22-FEB-00, corrected colour problem
;       Ver.9: PRY, 17-AUG-00, added windows identifying the numerator 
;                       and denominator transitions. This required an 
;                       extra routine (MAKE_STRINGS) and a common block 
;                       (ELVLC).
;       Ver.10: PRY, 25-AUG-00, added buttons to allow a choice between 
;                       ratios in energy or photon units.
;       Ver.11: PRY, 25-AUG-00, can now specify a wavelength range by 
;                       directly typing in the numbers
;       Ver.12: PRY, 30-Nov-00, changed call to emiss_select
;       Ver.13: PRY, 21-Dec-00, removed set_plot,'x' following help from 
;                       Bill Thompson
;       Ver.14, PRY, 27-Dec-00, changed switch to tst1 for IDL v5.4
;       Ver.15, PRY,  7-Dec-01, changed /prot keyword to /noprot to be 
;                       compatible with other CHIANTI routines.
;                       Added /temperature keyword.
;       Ver.16, PRY, 28-May-02, removed SPLINE calls, changing them to 
;                       SPL_INIT and SPL_INTERP; changed density labels to 
;                       temperature labels where appropriate.
;       Ver.17, PRY, 29-May-02, made treatment of photoexcitation consistent 
;                       with other CHIANTI routines.
;
;       V.  18, 29-May-2002, Giulio Del Zanna (GDZ)
;                        generalized directory concatenation to work for  Unix,
;                        Windows  and VMS.  
;                        Now we only call zion2filename.
;                        When creating the ps file, ps and psclose are used.
;
;       v.19, 12-Jun-02 GDZ 
;             fixed a small bug when finding the names of the files when the
;             keyword PATH is given.
;
;       v.20, 2-Aug-02, Peter Young
;             made some cosmetic changes following Giulio's suggestions.
;
;       v.21, 5-Aug-02, Peter Young
;             made various changes:
;               - the emissivity parameters (at top of widget) can only be 
;                 changed independently of line selection now
;               - the text widgets now check to make sure the user inputs 
;                 are numbers.
;               - extended the number of tags on the output EM structure
;
;       v.22, 6-Aug-02, Peter Young
;             corrected plot problem when viewing all temperature ratios 
;             in a fixed wavelength range.
;
;       v.23, 8-Aug-02, Peter Young
;             a number of cosmetic changes to make the GUI more 
;             user-friendly
;
;       v.24, 13-Aug-02, Peter Young
;             changed the temperature/density text widget so that numbers 
;             are registered even if the return key has not been hit.
;
;       v.25, 14-Aug-02, Peter Young
;             photoexcitation button now works again; also some cosmetic 
;             changes
;       V.26, 17-Sept-02, GDZ 
;             added scrolling in main frame
;
;       V.27, 10-Feb-03, Peter Young
;             corrected bug related to fonts on Windows PCs
;
;       V.28, 13-Feb-03, Peter Young
;             corrected problem when ion balance data dosen't exist for 
;             an ion
;
;       V.29, 18-Feb-03, Peter Young
;             added call to routine chianti_font to deal with fonts.
;
;       V.30, 28-Aug-03, Peter Young
;             corrected bug when two lines have the same wavelength by
;             introducing .ind tag to emiss_sel
;
;       V.31, 3-Nov-03, Giulio Del Zanna 
;             Added printout of Ne(Te) vs. ratio values
;             Modified format e8.2 to e9.2 for Windows compatibility.
;
;       V.32, 6-Nov-03, Peter Young
;             Corrected bug when new emissivities are calculated that
;             prevented the emissivites in emiss_sel from being updated.
;             Also, the Ne(Te), ratio values are now printed to a pop-up
;             window through the new 'SHOW RATIO VALUES' button rather
;             than printed to the text window (see V.31).
;
;       V.33, 7-Nov-03, Giulio Del Zanna (GDZ)
;             Modified the popup widget by calling xpopup and adding 
;             labels. Now it is possible to paste and copy the RATIO 
;             values into a file.
;
;       V.34, 12-Dec-03, Peter Young
;             Lowered the position of the postscript plot on the paper so
;             that the title isn't chopped off on US letter-size paper.
;
;       V.35, 2-Aug-2005, GDZ
;
;             Various modifications. Now the routine handles the
;             dielectronic case. Added printing of line lists in the
;             line ratio plots. Added the option to print directly to
;             the PRINTER. Added the option to SAVE a structure
;             containing the emissivities into an IDL save file (using
;             SAVEGEN). Added retain=2 to the window. Added
;             logarithmic or linear plot switch. Added cleaning of the
;             variables at the start of the routine. Various minor
;             cosmetic changes.
;
;       V.36, 12-Aug-2005, GDZ
;             Reinstated some previous cosmetic features. Also adjusted
;             some sizes, added some checks (e.g. minimum is set to 1e-10 if
;             plot is in log scale -to avoid log(0)-; automatic scaling cannot
;             be set if log scale is on).
;
;       V.37, 15-Aug-2005, Peter Young
;             Corrected bug with derived densities when the theoretical ratio
;             is double-valued.
;
;       V.38, 15-Aug-2005, Peter Young
;             Corrected bug introduced by modification above.
;
;       V.39, 17-Aug-2005, Peter Young
;             Adjusted layout of widget, modified printing of wavelengths
;             in the plot window, and adjusted plot limits for log scale.
;
;       V.40, 18-Aug-2005, Peter Young
;             Corrected bug for temperature version of routine that led to
;             mismatch between stated parameters and plotted parameters
;
;       V.41, 11-Oct-2007, Peter Young
;             Corrected bug when PATH= is specified.
;
;       V.42, 14-Jun-2010, Peter Young
;             For ions with a large number of transitions (e.g., Fe
;             XVIII) it was impossible to select lines with large
;             wavelengths as the integers weren't defined as
;             longs. This has been fixed now. Also fixed bug when
;             selecting strongest lines for ions with less than 15
;             transitions. 
;
; VERSION     :   42, 14-Jun-2010
;
;-


;------------------------------------------------------------------------------
PRO MAKE_STRINGS, I, STR
;------------------------------------------------------------------------------
;GDZ - changed emiss_data COMMON
COMMON emiss_data,iz,ion,diel,emiss,lo_dens,hi_dens,dint,dens,temp,units
COMMON elvlc, l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,ref


index_extractor,i,index=index

str = strarr(n_elements(index))

;GDZ
FOR j=0L,n_elements(index)-1 DO BEGIN
   l1 = emiss[index[j]].level1
   l2 = emiss[index[j]].level2
   wvl = strtrim(string(format='(f12.3)',emiss[index[j]].lambda),2)+ $
        ' '+string(197b)
   term1 = strtrim(term[l1-1],2)
   term2 = strtrim(term[l2-1],2)
   l1 = strtrim(string(format='(i4)',l1),2)
   l2 = strtrim(string(format='(i4)',l2),2)
   str[j] = wvl+'  '+l1+' - '+l2+'   '+term1+' - '+term2
ENDFOR

END


;------------------------------------------------------------------------------
PRO DENS_FINDER, INDEX, PS=PS
;------------------------------------------------------------------------------
; Uses the observed ratio (and 1-sigma error) to estimate the density.
;
; GDZ - changed plotting and emiss_data COMMONS
; PRY - corrected error for doubled-valued ratios
;
COMMON emiss_data,iz,ion,diel,emiss,lo_dens,hi_dens,dint,dens,temp,units
COMMON select,emiss_sel
COMMON sliders,lo_w,hi_w,min_w,max_w
COMMON plotting,lims,set_scale,all, loglin, line_labels

IF keyword_set(ps) THEN BEGIN
  th=2 
  usersym,[-1,1,1,-1,-1],[-1,-1,1,1,-1],th=2
  psym=8
ENDIF ELSE BEGIN
  th=1
  psym=6
ENDELSE


;----
; OBS_RATIO is the observed ratio value
; SIG_RATIO is the 1-sigma error on the observed ratio
; RATIO contains the theoretical values for the ratio
;
I_num=emiss_sel(index).obs & I_den=emiss_sel(0).obs
sig_num=emiss_sel(index).sig & sig_den=emiss_sel(0).sig
;
obs_ratio=I_num/I_den
sig_ratio=SQRT((I_num*sig_den)^2 + (I_den*sig_num)^2)/(I_den^2)
;
ratio=emiss_sel(index).em/emiss_sel(0).em


; -----
; Check if observed ratio lies within theoretical ratio range
;
IF (obs_ratio GT max(ratio)) OR (obs_ratio LT min(ratio)) THEN BEGIN
  result=DIALOG_MESSAGE(['The observed ratio lies outside of the plotted theoretical values.', $
                         'Try increasing the density range.'],/info)
  RETURN
ENDIF


;-----
; Need to be at least one point away from the endpoints of the density range
; in order to run a spline through the data-points (see code below)
;
result=MIN(abs(ratio-obs_ratio),ind)
;
IF (ind LT 1) OR (ind GT N_ELEMENTS(dens)-2) THEN BEGIN
  result=DIALOG_MESSAGE(['The ratio lies close to the edges of the plot.', $
                         'Try a different density range.'],/info)
  RETURN
ENDIF


;-------
; Define dens_start
;
CASE ind OF
  1: BEGIN 
       dens_start=dens(ind-1)
       ind=ind+1
     END
  N_ELEMENTS(dens)-2: BEGIN
       dens_start=dens(ind-3)
       ind=ind-1
     END
  ELSE: dens_start=dens(ind-2)
ENDCASE


;
; Fit a spline through the RATIO values on a x10 smaller density scale to
; derive the density value.
; ------
densi=FINDGEN(41)*dint/10. + dens_start
y2=spl_init(dens(ind-2:ind+2),ratio(ind-2:ind+2))
rati=spl_interp(dens(ind-2:ind+2),ratio(ind-2:ind+2),y2,densi)
;
result=MIN(ABS(rati-obs_ratio),ind)
;
OPLOT,[densi(ind),densi(ind)],[rati(ind),rati(ind)],psym=psym,symsiz=2
OPLOT,[densi(ind),densi(ind)],[rati(ind)+sig_ratio,rati(ind)-sig_ratio],th=th
;
PRINT,format='("Observed ratio:    ",f8.3,"+/-",f8.3)',obs_ratio,sig_ratio
PRINT,format='("Predicted density: ",f5.2)',densi(ind)
emiss_sel(index).dens=densi(ind)


;--------------------[O]
; Now derive the density (DENS1) corresponding to OBS_RATIO+SIG_RATIO
;
result=MIN(ABS(ratio-(obs_ratio+sig_ratio)),ind)
IF (ind LT 2) OR (ind GT N_ELEMENTS(dens)-3) THEN BEGIN
  dens1=-1
ENDIF ELSE BEGIN
  densi=FINDGEN(41)*dint/10. + dens(ind-2)
  y2=spl_init(dens(ind-2:ind+2),ratio(ind-2:ind+2))
  rati=spl_interp(dens(ind-2:ind+2),ratio(ind-2:ind+2),y2,densi)
 ;
  result=MIN(abs(rati-(obs_ratio+sig_ratio)),ind)
 ;
  dens1=densi(ind)
ENDELSE
;--------------------[O]

;--------------------(O)
; Now derive the density (DENS2) corresponding to OBS_RATIO-SIG_RATIO
;
result=MIN(abs(ratio-(obs_ratio-sig_ratio)),ind)
IF (ind LT 2) OR (ind GT N_ELEMENTS(dens)-3) THEN BEGIN
  dens2=-1
ENDIF ELSE BEGIN
  densi=FINDGEN(41)*dint/10. + dens(ind-2)
  y2=spl_init(dens(ind-2:ind+2),ratio(ind-2:ind+2))
  rati=spl_interp(dens(ind-2:ind+2),ratio(ind-2:ind+2),y2,densi)
 ;
  result=MIN(ABS(rati-(obs_ratio-sig_ratio)),ind)
 ;
  dens2=densi(ind)
ENDELSE
;--------------------(O)

;--------
; Deal with the case when DENS1 or DENS2 are outside the range
;
tst1=((dens2 EQ -1) AND (dens1 EQ -1))+(dens2 EQ -1)
;
CASE tst1 OF
 ;
  0: BEGIN                                  ;---dens2 not equal to -1
       IF (dens2 GE emiss_sel(index).dens) THEN BEGIN
         emiss_sel(index).densup=dens2
         emiss_sel(index).denslo=dens1
       ENDIF ELSE BEGIN
         emiss_sel(index).densup=dens1
         emiss_sel(index).denslo=dens2
       ENDELSE

     END
 ;
  1: BEGIN                                  ;---only dens2 equal to -1
       IF (dens1 GE emiss_sel(index).dens) THEN BEGIN
         emiss_sel(index).densup=dens1
         emiss_sel(index).denslo=dens2
       ENDIF ELSE BEGIN
         emiss_sel(index).densup=dens2
         emiss_sel(index).denslo=dens1
       ENDELSE
     END
 ;
  2: BEGIN                                  ;---dens1 and dens2 equal to -1
       emiss_sel(index).densup=-1
       emiss_sel(index).denslo=-1
     END
 ;
ENDCASE


IF ((obs_ratio+sig_ratio) GT MAX(ratio)) THEN BEGIN
  ind=WHERE(ratio EQ MAX(ratio))
  IF dens(ind(0)) GT emiss_sel(index).dens $
     THEN emiss_sel(index).densup=-1 $
     ELSE emiss_sel(index).denslo=-1
ENDIF
IF ((obs_ratio-sig_ratio) LT MIN(ratio)) THEN BEGIN
  ind=WHERE(ratio EQ MIN(ratio))
  IF dens(ind(0)) GT emiss_sel(index).dens $
     THEN emiss_sel(index).densup=-1 $
     ELSE emiss_sel(index).denslo=-1
ENDIF


END

;;-----------------------------------------------------------------------------
PRO INDEX_EXTRACTOR, I, INDEX=INDEX, PLOT_LABEL=PLOT_LABEL, PS=PS
;------------------------------------------------------------------------------

;+
; This procedure has a dual purpose: (i) to extract `emiss' indices from
; `emiss_sel' and, (ii) create labels that will be displayed on the plot
; window.
; The lines that comprise the blend are contained in emiss_sel.label 
; as, e.g., '203.79+203.80'. This routine separates this string into 
; the separate wavelengths and looks in emiss.lambda to see which 
; indices they correspond to. The set of indices are output through INDEX.
;
; INPUT
;
;   I    Index of component in emiss_sel. E.g., for the denominator I=0,
;        for the numerator I=1.
;
; OPTIONAL OUTPUT
;
;   INDEX Integer array of same length as the number of lines in the
;         numerator or denominator (based on value of I). The elements
;         give the indices of the lines as they appear in the EMISS
;         structure.
;
;   PLOT_LABEL A string containing the index to be used when labelling
;              the line ratio plots.
;
; KEYWORD
;
;   PS    Setting this adds an angstrom to the wavelength labels. (Only
;         works for the postscript device.)
;-
;GDZ - changed emiss_data COMMON

COMMON emiss_data,iz,ion,diel,emiss,lo_dens,hi_dens,dint,dens,temp,units
COMMON select, emiss_sel

IF keyword_set(ps) THEN ang=string(197b) ELSE ang=''

n=emiss_sel(i).n_ind - 1
index=lonarr(n+1)

lresult=STR_SEP(emiss_sel(i).label,'+')
result=str_sep(emiss_sel[i].ind,'+')
index=long(result)

;
; In the next bit, if the ion had a large number of lines (>1024), IDL 
; does not like it when you convert emiss.lambda to a string array. 
; Thus instead I've used 'ind' to pick out only a sub-set of the 
; wavelengths.
;
; IF N_ELEMENTS(index) NE 0 THEN BEGIN
;   FOR j=0,n DO begin
;     ind=WHERE( ABS( FLOAT(result(j)) - emiss.lambda ) LT 2. )
;     index(j)=WHERE(result(j) EQ STRTRIM(STRING(FORMAT='(f10.3)', $
;                          emiss(ind).lambda),2))
;     index(j)=ind(index(j))
;   ENDFOR
; ENDIF


IF N_ELEMENTS(plot_label) THEN BEGIN
  plot_label=lresult(0)+ang
  IF n GT 0 THEN FOR j=1,n DO plot_label=plot_label+'!c'+lresult(j)+ang
ENDIF

END

;------------------------------------------------------------------------------
PRO DENS_PLOT, state=state, ps=ps
;------------------------------------------------------------------------------

; This is where the line ratios get plotted. The common block 'plotting' is
; used to denote what type of scaling is being used (set_scale) and whether
; all lines are plotted or just selected lines (all)

;GDZ - changed plotting and emiss_data COMMONS
COMMON emiss_data,iz,ion,diel,emiss,lo_dens,hi_dens,dint,dens,temp,units
COMMON select,emiss_sel
COMMON sliders,lo_w,hi_w,min_w,max_w
COMMON plotting,lims,set_scale,all, loglin, line_labels
COMMON temp, t_switch
COMMON proton_rates, prot_nincl,prot_switch,prot_mess,pname
COMMON proton, pstr, pe_ratio

;help, lims,set_scale,all, loglin
;stop

;-------------------------------<
; The following definition creates space at the side of the plot for
; displaying labels. To change space, change 0.8 to a different value.
;
nd=n_elements(dens)
width=(max(dens)-min(dens))/0.8

;GDZ-
xrange=[min(dens), width + min(dens) + 0.2]
;xrange=[min(dens), width + min(dens)]
;-------------------------------<

;GDZ-
  ang=' '+string(197b)
IF keyword_set(ps) THEN BEGIN
  th=2
charsize = 0.6
ENDIF ELSE BEGIN
  th=1
charsize = 1.1
ENDELSE

;GDZ
zion2spectroscopic,iz,ion,title, diel=diel

IF keyword_set(all) THEN rat_txt='ratios' ELSE rat_txt='ratio'
title=title+' line '+rat_txt+' relative to '+emiss_sel(0).label+ang

IF t_switch EQ 1 THEN xtit='Log!d10!n (Electron temperature [K])' $
ELSE xtit='Log!d10!n (Electron density [cm!u-3!n])'

IF units EQ 0 THEN YTIT='Ratio (Energy)'  ELSE IF $
    units EQ 1 THEN YTIT='Ratio (Photons)'

;
; adjust upper limit of plot depending on whether it's log or linear
;
IF loglin EQ 1 THEN upfact=2.0 ELSE upfact=1.15

IF all eq 0 THEN BEGIN         ;;----------------; PLOT SELECTED LINES
 ;
  IF set_scale EQ 1 THEN BEGIN        ;---get y-scale for plot
    yrange=lims 
  ENDIF ELSE BEGIN
    max_rat=0. & min_rat=0.01
    FOR i=1,4 DO BEGIN
      max_i=MAX(emiss_sel(i).em/emiss_sel(0).em)
      IF max_i GT max_rat THEN max_rat=max_i
      min_i = min(emiss_sel(i).em/emiss_sel(0).em)
      IF min_i GT min_rat THEN min_rat = min_i
    ENDFOR
    IF set_scale EQ 2 THEN yrange=[min_rat*0.95,(upfact-0.1)*max_rat] $
         ELSE yrange=[0.,upfact*max_rat]
    lims=yrange
    WIDGET_CONTROL,state.lo_read,set_value=string(format='(e10.2)',lims[0])
    WIDGET_CONTROL,state.up_read,set_value=string(format='(e10.2)',lims[1])
  ENDELSE
 ;
  line_label=strarr(5)          ;---construct the line labels
  FOR i=0,4 DO BEGIN
    IF emiss_sel(i).n_ind NE -1 THEN BEGIN
      result=''
      index_extractor,i,plot_label=result,ps=ps
      line_label(i)=result
    ENDIF
  ENDFOR

; GDZ
  PLOT,dens,emiss_sel(1).em/emiss_sel(0).em,XSTY=1, $     ;---plot first ratio
        XRANGE=xrange,YTICKLEN=-0.015, $
        YRANGE=yrange,YSTY=1,XTICKLEN=-0.015,TITLE=title, $
        XTITLE=xtit, YTITLE=ytit,th=th,xth=th,yth=th, ylog=loglin

 ;
   XYOUTS,min(dens)+width*0.83, $                          ;---add first label
            emiss_sel(1).em(nd-1)/emiss_sel(0).em(nd-1),$
            line_label(1)

;   dummy0 = 'CHIANTI Version '
;   XYOUTS,0.5, 0.9, align=0.5, /norm, dummy0

;   IF  PROT_SWITCH THEN dummy1 = 'Proton rates included, p/e='+trim(pe_ratio) $
;        ELSE dummy1 = 'Proton rates not included' 
;   XYOUTS,0.5, 0.85, align=0.5, /norm, dummy1

;   dummy2 = 'Calculated at Log T='+trim(temp)
;   XYOUTS,0.5, 0.8, align=0.5, /norm, dummy2

 ;
  IF emiss_sel(1).obs NE 0. THEN dens_finder,1,ps=ps    ;--add observed ratio
 ;
  IF emiss_sel(2).label NE '' THEN oplot,dens,emiss_sel(2).em/emiss_sel(0).em,$
       th=2
 ;
  IF emiss_sel(3).label NE '' THEN oplot,dens,emiss_sel(2).em/emiss_sel(0).em,$
       th=2


 show_densities, ratio_txt, dens_txt, disp_string, full_string, err=err

 IF err EQ '' THEN BEGIN 

;   XYOUTS,0.5, 0.75, align=0.5, /norm,'Observed Ratio:  '+ratio_txt
;                  IF t_switch EQ 1 THEN  dummy = 'Log T [K]:  '+dens_txt ELSE $
;             dummy = ' Log N_e [cm-3]: '+dens_txt

;   XYOUTS,0.5, 0.7, align=0.5, /norm,dummy

;stop

 END 
;
;
;
;
ENDIF ELSE BEGIN              ;-------;  PLOT ALL LINES IN WAVELENGTH RANGE
 ;
  index=WHERE(emiss.lambda GE lo_w AND emiss.lambda LE hi_w)
 ;
  IF index(0) EQ -1 THEN BEGIN         ;---no lines in specified range
    result=WIDGET_MESSAGE('No lines in specified range',/info)
  ENDIF ELSE BEGIN
    IF set_scale EQ 1 THEN BEGIN         ;---get y-scale for plot
      yrange=lims 
    ENDIF ELSE BEGIN
      max_rat=0. & min_rat=1e10
      FOR i=0,nd-1 DO BEGIN
        max_i=MAX(emiss(index).em(i)/emiss_sel(0).em(i),m_ind)
        min_i=min(emiss(index).em(i)/emiss_sel(0).em(i),m_ind)
        IF max_i GT max_rat THEN max_rat=max_i
        IF min_i GT min_rat THEN min_rat=min_i
      ENDFOR
      IF set_scale EQ 2 THEN yrange=[0.95*min_i,(upfact-0.1)*max_rat] $
           ELSE yrange=[0.,upfact*max_rat]
      lims=yrange
      WIDGET_CONTROL,state.lo_read,set_value=string(format='(e10.2)',lims[0])
      WIDGET_CONTROL,state.up_read,set_value=string(format='(e10.2)',lims[1])
  ENDELSE

; GDZ- added log option

    PLOT,xrange,yrange,/NODATA,YSTY=1,XSTY=1, $       ;---show plot axes
         TITLE=title,XTICKLEN=-0.015,YTICKLEN=-0.015, 	$
         XTITLE=xtit, $
         YTITLE='Ratio',th=th,xth=th,yth=th, ylog=loglin
   ;

;GDZ - added definition of labels and points where to plot

line_labels = ' '
ypoints = -1

    FOR i=1L,N_ELEMENTS(index) DO BEGIN

  ypoint = emiss(index(i-1)).em(nd-1)/emiss_sel(0).em(nd-1)

IF ypoint LE yrange[1] AND ypoint GE yrange[0] THEN BEGIN 

      OPLOT,dens,emiss(index(i-1)).em/emiss_sel(0).em,th=th  ;---plot lines

     ;
     ; adjust no. of decimal places of wavelength for printing
     ;
      lambda=emiss(index(i-1)).lambda
      CASE 1 OF
        lambda LT 50.: wvl_label=trim(string(format='(f12.4)',lambda))
        lambda GT 10000.: wvl_label=trim(string(format='(f12.1)',lambda))
        ELSE: wvl_label=trim(string(format='(f12.3)',lambda))
      ENDCASE
     ;
      line_label=trim(emiss(index(i-1)).level1)+'-'+$
                 trim(emiss(index(i-1)).level2)+' '+$
           wvl_label
;      strtrim( string(format='(f12.4)',emiss(index(i-1)).lambda),2)
      
      XYOUTS,min(dens)+width*0.83, ypoint, $                          ;---add labels
              line_label+ang, charsize=charsize

ypoints =[ypoints, ypoint]

line_labels = [line_labels, $
    trim(emiss(index(i-1)).level1)+'-'+$
                 trim(emiss(index(i-1)).level2)+' '+$
      strtrim(string(format='(f12.4)',emiss(index(i-1)).lambda),2)+ang+$
       ' '+ emiss(index(i-1)).LVL1_DESC+' - '+ emiss(index(i-1)).LVL2_DESC]

ENDIF  
    ENDFOR 

IF n_elements(line_labels) GT 1 THEN BEGIN 
line_labels = line_labels[1:*]
ypoints = ypoints[1:*]

isort = reverse(sort(ypoints))
line_labels = line_labels(isort)
endif 

  ENDELSE
 ;
ENDELSE

IF t_switch EQ 1 THEN xtit='Log Electron temperature [K] :' $
ELSE xtit='Log Electron density [cm-3] :'

IF units EQ 0 THEN YTIT='Ratio (Energy) values: '  ELSE IF $
    units EQ 1 THEN YTIT='Ratio (Photons) values: '

dummy1 = XTIT+arr2str(trim(float(dens)), ',',/trim)
dummy2 = YTIT+arr2str(trim(float(emiss_sel(1).em/emiss_sel(0).em)), ',',/trim)

make_currplot_label,label
;GDZ
;widget_control,state.info_txt,set_value=[dummy1, dummy2,label]
widget_control,state.info_txt,set_value=label

make_dens_label,label

;GDZ
;IF n_elements(line_labels) GT 0 THEN $
;  widget_control,state.dens_txt,set_value=line_labels, /append

widget_control,state.dens_txt,set_value=label

END


;------------------------------------------------------------------------------
PRO wavel_plot
;------------------------------------------------------------------------------

; creates the little plot at the bottom left of the widget which shows the 
; distribution of the ion's lines with wavelength.

COMMON emiss_data,iz,ion,diel,emiss,lo_dens,hi_dens,dint,dens,temp,units
COMMON sliders,lo_w,hi_w,min_w,max_w
COMMON Comm, plot_rat_id, plot_spec_id

wset,plot_spec_id
device,decomposed=0
;plot,[min_w,max_w],[0,1],/nodata,xticklen=0.00001,yticklen=0.000001,$
;     xsty=1,yticks=1,charsiz=.01
plot,[min_w,max_w],[0,1],/nodata,pos=[0,0,1,1],xsty=5,ysty=5
n_el=long(n_elements(emiss))
for i=0l,n_el-1 do $
     oplot,alog10([emiss(i).lambda,emiss(i).lambda]),[0,1]
oplot,alog10([lo_w,lo_w]),[0,1],th=3,col=50
oplot,alog10([hi_w,hi_w]),[0,1],th=3,col=50
wset,plot_rat_id
device,decomposed=0

END

;GDZ - added struct_info to be passed when saving the output 
; structure

PRO make_currplot_label, label, struct_info= struct_info

COMMON emiss_data,iz,ion,diel,emiss,lo_dens,hi_dens,dint,dens,temp,units
COMMON select,emiss_sel
COMMON sliders,lo_w,hi_w,min_w,max_w
;GDZ
COMMON plotting,lims,set_scale,all, loglin, line_labels
COMMON Comm, plot_rat_id, plot_spec_id
COMMON rad_data, radtemp,rphot,r_tst
COMMON proton_rates, prot_nincl,prot_switch,prot_mess,pname
COMMON proton, pstr, pe_ratio
COMMON extra, fpath
COMMON elvlc, l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,ref
COMMON files, ioneq_file, abund_file
COMMON temp, t_switch

IF t_SWITCH EQ 1 THEN BEGIN
  td_label='Log Dens [cm-3]: '+string(format='(f5.2)',temp)
  int_label='Temperature intervals: '+strtrim(string(format='(f5.1)',dint),2)

struct_info = {temperature:10.^dens, density:10.^temp, $
COMMENT1:'Calculated with CHIANTI at constant '+td_label, $
COMMENT2:int_label}

ENDIF ELSE BEGIN
  td_label='Log Temp [K]: '+string(format='(f5.2)',temp)
  int_label='Density intervals: '+strtrim(string(format='(f5.1)',dint),2)

struct_info = {density:10.^dens, temperature:10.^temp, $
COMMENT1:'Calculated with CHIANTI at constant '+td_label, $
COMMENT2:int_label}

ENDELSE


IF r_tst EQ 1 THEN BEGIN
  IF rphot GE 1000. THEN rphot_str=string(format='(e9.2)',rphot) ELSE $
       rphot_str=strtrim(string(format='(f6.2)',rphot),2)
  IF radtemp GT 1d5 THEN rt_str=strtrim(string(format='(e9.2)',radtemp),2) $
       ELSE rt_str=strtrim(string(long(radtemp)),2)
  pexc_label='Photoexc: rphot='+rphot_str+', rtemp='+rt_str+'K'
ENDIF ELSE BEGIN
  pexc_label='Photoexc: not included'
ENDELSE

prot_label='Protons: '+strlowcase(prot_mess[prot_SWITCH])

struct_info = join_struct(struct_info,$
  {COMMENT3:pexc_label,COMMENT4:prot_label})

IF units EQ 1 THEN BEGIN 
units_label='Units: photons' 
struct_info = join_struct(struct_info,{Units:'photons'})
ENDIF  ELSE BEGIN 
units_label='Units: energy'
struct_info = join_struct(struct_info,{Units:'energy'})
END 

;label=['CURRENT PLOT PARAMETERS', $
;       td_label,int_label,pexc_label,prot_label,units_label]

;GDZ
label=[td_label , int_label, pexc_label, prot_label, units_label]

END


;------------------------------------------------------------------------------
PRO make_dens_label, label
;------------------------------------------------------------------------------

;
; This creates the label that goes into the dens_txt widget, showing the 
; derived value of the density/temperature and error bars
;

COMMON emiss_data,iz,ion,diel,emiss,lo_dens,hi_dens,dint,dens,temp,units
COMMON select,emiss_sel
COMMON temp, t_switch
;GDZ
COMMON plotting,lims,set_scale,all, loglin, line_labels


IF t_SWITCH EQ 1 THEN BEGIN
  head_label='DERIVED TEMPERATURE [K]:'
ENDIF ELSE BEGIN
  head_label='DERIVED DENSITY [cm-3]:'
ENDELSE

IF all EQ 1 THEN label=[head_label,'Not applicable']


density=emiss_sel[1].dens
denslo=emiss_sel[1].denslo
densup=emiss_sel[1].densup
;
ratio=emiss_sel[1].obs/emiss_sel[0].obs
sigma=sqrt((emiss_sel[1].obs*emiss_sel[0].sig)^2+ $
     (emiss_sel[0].obs*emiss_sel[1].sig)^2)/emiss_sel[0].obs^2

rat_label='Ratio: '+$
     strtrim(string(format='(f6.2)',ratio),2)
IF sigma NE 0. THEN rat_label=rat_label+' '+string(177b)+' '+$
     strtrim(string(format='(f6.2)',sigma),2)

ld_str=strtrim(string(format='(f8.2)',density),2)
d_str=strtrim(string(format='(e9.2)',10.^density),2)
IF denslo EQ -1. THEN BEGIN
  ldlo_str='LOW' 
  dlo_str='LOW'
ENDIF ELSE BEGIN
  ldlo_str=strtrim(string(format='(f8.2)',denslo),2)
  dlo_str=strtrim(string(format='(e9.2)',10.^(denslo)),2)
ENDELSE
IF densup EQ -1. THEN BEGIN
  ldhi_str='HIGH' 
  dhi_str='HIGH'
ENDIF ELSE BEGIN
  ldhi_str=strtrim(string(format='(f8.2)',densup),2)
  dhi_str=strtrim(string(format='(e9.2)',10.^(densup)),2)
ENDELSE


CASE 1 OF
  ratio EQ 0.: BEGIN
    label=[head_label,'Input line intensities']
  END
 ;
  (ratio NE 0.) AND (sigma EQ 0.): BEGIN
    label=[head_label,'     '+d_str+'   [log='+ld_str+']', $
          '','','',rat_label]
  END
 ;
  ELSE: BEGIN
    label=[head_label,'     '+d_str+'   [log='+ld_str+']', $
           '  + '+dhi_str+'   [log='+ldhi_str+']',$
           '  - '+dlo_str+'   [log='+ldlo_str+']','', $
          rat_label]
  ENDELSE
ENDCASE


END


PRO show_densities, ratio_txt, dens_txt, disp_string,full_string, err=err

;GDZ
COMMON emiss_data,iz,ion,diel,emiss,lo_dens,hi_dens,dint,dens,temp,units
COMMON select,emiss_sel
COMMON temp, t_switch

err = ''

          n_lines=intarr(5)
           FOR i=1,4 DO BEGIN          ; for each numerator
             IF emiss_sel(i).obs NE 0. THEN n_lines(i)=1
           ENDFOR
           ind=WHERE(n_lines EQ 1)
           IF ind(0) NE -1 THEN BEGIN
             line_labels=emiss_sel(where(n_lines EQ 1)).label+'/'$
                             +emiss_sel(0).label
             line_text=STRARR(N_ELEMENTS(ind))

             FOR i=0,N_ELEMENTS(ind)-1 DO BEGIN
              ;
              ; Create string containing ratio plus errors
              ;
               sig1=emiss_sel(ind(i)).sig
               sig2=emiss_sel(0).sig
               i1=emiss_sel(ind(i)).obs
               i2=emiss_sel(0).obs
               ratio=i1/i2
               sig=sqrt( (i1*sig2)^2 + (i2*sig1)^2 )/i2^2
               ratio_txt=strtrim(string(format='(f10.3)',ratio),2)+'+/-'+$
                         strtrim(string(format='(f10.3)',sig),2)
              ;
              ; Create string containing density plus errors
              ;
               density=emiss_sel(ind(i)).dens
               densup=emiss_sel(ind(i)).densup
               IF densup EQ -1 THEN BEGIN
                 densup='HI'
               ENDIF ELSE BEGIN
                 densup=strtrim(string(format='(f10.2)',densup-density),2)
               ENDELSE
               denslo=emiss_sel(ind(i)).denslo
               IF denslo EQ -1 THEN BEGIN
                 denslo='LO'
               ENDIF ELSE BEGIN
                 denslo=strtrim(string(format='(f10.2)',density-denslo),2)
               ENDELSE
               density=strtrim(string(format='(f10.2)',density),2)
               dens_txt=density+'+'+densup+'-'+denslo
              ;
               IF i EQ 0 THEN BEGIN
                 disp_string=[line_labels(i),$
                    '    Observed Ratio:  '+ratio_txt]
                 IF t_switch EQ 1 THEN  $
                      disp_string=[disp_string,'    Log T [K]:  '+dens_txt] $
                 ELSE disp_string=[disp_string,'    Log N_e [cm-3]: '+dens_txt]
                 disp_string=[disp_string,'']

                 IF t_SWITCH EQ 1 THEN tede_str='temperature' ELSE $
                      tede_str='density'

                 full_string=[disp_string,'The accuracy of the derived '+$
               tede_str,'depends on the '+tede_str+$
               ' intervals','you are using. The smaller the intervals,',$
               'the higher the accuracy.']
               ENDIF
;???
               IF i NE 0 THEN BEGIN
                 disp_string=[disp_string,line_labels(i),$
                              'Observed Ratio:  '+ratio_txt]
                 IF t_switch EQ 1 THEN  $
                      disp_string=[disp_string,'Temperature:  '] ELSE $
                      disp_string=[disp_string,'Density:  ']
               ENDIF
             ENDFOR
ENDIF ELSE err = 'no lines'
END 


;------------------------------------------------------------------------------
PRO DENS_MAIN_Event, Event
;------------------------------------------------------------------------------

;GDZ - changed plotting and emiss_data COMMONS
COMMON emiss_data,iz,ion,diel,emiss,lo_dens,hi_dens,dint,dens,temp,units
COMMON select,emiss_sel
COMMON sliders,lo_w,hi_w,min_w,max_w
COMMON plotting,lims,set_scale,all, loglin, line_labels
COMMON Comm, plot_rat_id, plot_spec_id
COMMON rad_data, radtemp,rphot,r_tst
COMMON proton_rates, prot_nincl,prot_switch,prot_mess,pname
COMMON proton, pstr, pe_ratio
COMMON extra, fpath
COMMON elvlc, l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,ref
COMMON files, ioneq_file, abund_file
COMMON temp, t_switch


WIDGET_CONTROL,Event.top, get_uvalue=state ;,/no_copy

txt_mess='Please type a valid number'

CASE 1 OF 

    event.id EQ state.ch_params: BEGIN
        widget_control,state.ch_params,get_uval=chck

        IF chck EQ 0 THEN BEGIN
            widget_control,state.ch_params,set_val='CALCULATE NEW EMISSIVITIES'
            widget_control,state.ch_params,set_uval=1
                                ;
            widget_control,state.base_2,sensitive=0
            widget_control,state.base_sub1,sensitive=1
        ENDIF

        IF chck EQ 1 THEN BEGIN

            WIDGET_CONTROL,state.temp_read,get_value=aa
            IF is_number(aa) NE 0 THEN tchck=float(aa[0]) ELSE tchck=-1.
            IF tchck GT 0. THEN BEGIN
                temp=tchck
            ENDIF ELSE BEGIN
                str1=['Temperature','Density']
                txt=str1[t_switch]+': please type a valid number (> 0)'
                result=dialog_message(txt)
                WIDGET_CONTROL,state.temp_read, $
                  set_value=strtrim(string(format='(f5.2)',temp),2)
                GOTO,lbl2
            ENDELSE

            WIDGET_CONTROL,state.rt_read,get_value=aa
            IF is_number(aa) NE 0 THEN rchck=float(aa[0]) ELSE rchck=-1.
            IF rchck GT 0. THEN BEGIN
                radtemp=rchck
            ENDIF ELSE BEGIN
                txt='PHOTOEXCITATION: please type a valid number (> 0) for the'+$
                  ' radiation temperature'
                result=dialog_message(txt)
                WIDGET_CONTROL,state.rt_read,set_value=string(format='(f7.1)',radtemp)
                GOTO,lbl2
            ENDELSE

            WIDGET_CONTROL,state.dil_read,get_value=aa
            IF is_number(aa) NE 0 THEN rchck=float(aa[0]) ELSE rchck=-1.
            IF rchck GE 1. THEN BEGIN
                rphot=rchck
            ENDIF ELSE BEGIN
                txt='PHOTOEXCITATION: please type a valid number (> or = 1) for the'+$
                  ' distance'
                result=dialog_message(txt)
                rphot_str=strtrim(string(format='(f7.2)',rphot),2)
                WIDGET_CONTROL,state.dil_read,set_value=rphot_str
                GOTO,lbl2
            ENDELSE

            nd=ROUND((hi_dens-lo_dens)/dint +1)
            dens=FINDGEN(nd)*dint + lo_dens

            WIDGET_CONTROL,/hourglass

            IF t_switch EQ 1 THEN BEGIN
                IF r_tst EQ 1 THEN BEGIN
                    emiss=emiss_calc(iz,ion,temp=dens,dens=temp,rphot=rphot, $
                                     path=fpath,noprot=noprot,radt=radtemp, $
                                     /quiet, ioneq_file=ioneq_file,DIEL=DIEL, $
                                     abund_file=abund_file,no_de=units)
                ENDIF ELSE BEGIN
                    emiss=emiss_calc(iz,ion,temp=dens,dens=temp, $
                                     path=fpath,noprot=noprot, $
                                     /quiet, ioneq_file=ioneq_file,DIEL=DIEL, $
                                     abund_file=abund_file,no_de=units)
                ENDELSE
            ENDIF ELSE BEGIN
                IF r_tst EQ 1 THEN BEGIN
                    emiss=emiss_calc(iz,ion,dens=dens,temp=temp,rphot=rphot, $
                                     radt=radtemp,path=fpath,noprot=prot_nincl,/quiet, $
                                     no_de=units, ioneq_file=ioneq_file,DIEL=DIEL, $
                                     abund_file=abund_file)
                ENDIF ELSE BEGIN
                    emiss=emiss_calc(iz,ion,dens=dens,temp=temp, $
                                     path=fpath,noprot=prot_nincl,/quiet, $
                                     no_de=units, ioneq_file=ioneq_file,DIEL=DIEL, $
                                     abund_file=abund_file)
                ENDELSE
            ENDELSE

            labels=emiss_sel.label & n_inds=emiss_sel.n_ind ; save info before 
            obs=emiss_sel.obs & sig=emiss_sel.sig ; destroying emiss_sel
            inds=emiss_sel.ind

            str={label:'', em:fltarr(nd), n_ind:-1, obs:0., sig:0., $
                 dens:0., densup:0., denslo:0., ind:''}
            emiss_sel=replicate(str,5) ; recreate emiss_sel

            IF prot_nincl EQ 0 THEN BEGIN
                result=FINDFILE(EXPAND_PATH(pname))
                IF result(0) NE '' THEN BEGIN
                    prot_switch=1
                ENDIF ELSE BEGIN
                    prot_switch=2
                    prot_nincl=1
                    WIDGET_CONTROL,state.prot_bgrp,set_value=0
                ENDELSE
            ENDIF ELSE BEGIN
                prot_switch=0
            ENDELSE
            WIDGET_CONTROL,state.prot_disp,set_value=prot_mess(prot_switch)

            FOR i=0,4 DO BEGIN
                emiss_sel(i).label=labels(i) & emiss_sel(i).n_ind=n_inds(i)
                emiss_sel(i).obs=obs(i) & emiss_sel(i).sig=sig(i)
                emiss_sel[i].ind=inds[i]
                IF emiss_sel(i).n_ind NE -1 THEN BEGIN
                    index_extractor,i,index=index
                    IF emiss_sel(i).n_ind EQ 1 THEN BEGIN
                        emiss_sel[i].em=reform(emiss[index[0]].em)
                    ENDIF ELSE BEGIN
                        em_all=reform(emiss[index].em)
                        emiss_sel[i].em=reform(total(em_all,2))
                    ENDELSE
                ENDIF
            ENDFOR
            dens_plot,state=state

            widget_control,state.ch_params,set_val='CHANGE PARAMETERS'
            widget_control,state.ch_params,set_uval=0
                                ;
            widget_control,state.base_2,sensitive=1
            widget_control,state.base_sub1,sensitive=0
        ENDIF
        lbl2:
    END


    event.id EQ state.units_buts: widget_control,state.units_buts, $
      get_value=units

    event.id EQ state.prot_bgrp: BEGIN
        WIDGET_CONTROL,state.prot_bgrp,get_value=aa
        prot_nincl=aa
    END

    event.id EQ state.lo_dens_slid: BEGIN ;---dens slider
        IF t_switch EQ 1 THEN BEGIN
            lo_dens=(event.value)/5. + 2.0
            widget_control,state.ld_label, $
              set_value=string(format='(f4.1)',lo_dens)
        ENDIF ELSE BEGIN
            lo_dens=event.value
        ENDELSE
    END

    event.id EQ state.hi_dens_slid: BEGIN ;---dens slider
        IF t_switch EQ 1 THEN BEGIN
            hi_dens=(event.value)/5. + 2.0
            widget_control,state.hd_label, $
              set_value=string(format='(f4.1)',hi_dens)
        ENDIF ELSE BEGIN
            hi_dens=event.value
        ENDELSE
    END

    event.id EQ state.dens_bgrp: BEGIN ;---density interval
                                ;                                                  ;   buttons
        CASE event.value OF
                                ;
            0: BEGIN
                WIDGET_CONTROL,state.dens_bgrp,get_value=aa,get_uvalue=bb
                dint=float(bb(aa))
            END
                                ;
            1: BEGIN
                WIDGET_CONTROL,state.dens_bgrp,get_value=aa,get_uvalue=bb
                dint=float(bb(aa))
            END
                                ;
            2: BEGIN
                WIDGET_CONTROL,state.dens_bgrp,get_value=aa,get_uvalue=bb
                dint=float(bb(aa))
            END
                                ;
            3: BEGIN
                WIDGET_CONTROL,state.dens_bgrp,get_value=aa,get_uvalue=bb
                dint=float(bb(aa))
            END
                                ;
        ENDCASE
    END


    event.id EQ state.pexc_bgrp: BEGIN
        widget_control,state.pexc_bgrp,get_value=getval
        IF getval EQ 0 THEN BEGIN
            widget_control,state.dil_set,map=1
            widget_control,state.rt_set,map=1
            r_tst=1
        ENDIF ELSE BEGIN
            widget_control,state.dil_set,map=0
            widget_control,state.rt_set,map=0
            r_tst=0
        ENDELSE
    END

;
; note: I used to allow the dilution factor to be input, but now I allow the
;  distance to be set. I didn't change the name of the widgets, though.
;
;   event.id eq state.dil_read: BEGIN                  ;---read distance
;      WIDGET_CONTROL,state.dil_read,get_value=rphot
;      rphot=float(rphot[0])
;      IF rphot LT 1.0 THEN BEGIN
;        result=dialog_message([ $
;           'Distance is the radial distance of the plasma from the', $
;           'centre of the star, in stellar radii units, and must', $
;           'be great-or-equal to 1, corresponding to the plasma', $
;           'being either on or above the stellar surface.', $
;           '', $
;           'The value you have specified is less than 1 and will', $
;           'be changed to 1.'])
;        rphot=1.0
;        widget_control,state.dil_read,set_value='1.00'
;      ENDIF
;   END

;   event.id eq state.rt_read: BEGIN                   ;---read rad temp
;      WIDGET_CONTROL,state.rt_read,get_value=bob
;      radtemp=float(bob) & radtemp=radtemp(0)
;   END

;;-------------------------------------------------; LINE SELECTION
;;The following bit allows new lines to be inserted 
;;into emiss_sel, and also allows blends to be analysed.
;;`identity' determines which of the denominator or
;;numerators is being altered.
;;
;; To add more numerators simply add an extra OR to the line below
;
    (event.id eq state.d_lines) OR (event.id eq state.n1_lines): BEGIN
                                ;
        WIDGET_CONTROL,event.id,get_uvalue=identity ; identity is an integer
                                ;                                             ; that identifies what is
                                ;                                             ; being altered
        CASE event.value of
            1: BEGIN
                n=n_elements(dens)
                index=reverse( sort (emiss.em(fix(n/2))) )
                index=index(0:min([14,n_elements(emiss.lambda)-1])) ; 15 strongest lines.
                index2=sort(emiss(index).lambda) ; Sort the strongest lines into
                index=index(index2) ; wavelength order.
                result=emiss_select(emiss,index,sel_ind=ind,group=event.top)
                IF result(0) NE -1 THEN BEGIN
                    emiss_sel(identity).em=result
                    emiss_sel(identity).n_ind=n_elements(ind)
                    emiss_sel(identity).sig=0.
                    IF identity EQ 0 THEN emiss_sel(identity).obs=1. $
                    ELSE emiss_sel(identity).obs=0.
                    WIDGET_CONTROL,state.ints(identity), $
                      set_value=strtrim(string(format='(f5.1)', $
                                               emiss_sel(identity).obs),2)
                    WIDGET_CONTROL,state.sigs(identity), $
                      set_value=strtrim(string(format='(f5.1)', $
                                               emiss_sel(identity).sig),2)
                    label=''
                    indstr=''
                    for i=1,n_elements(ind) do BEGIN
                        label=label+'+'+strtrim(string(format='(f10.3)',$
                                                       emiss(ind(i-1)).lambda),2)
                        IF indstr EQ '' THEN indstr=trim(ind[i-1]) ELSE $
                          indstr=indstr+'+'+trim(ind[i-1])
                    ENDFOR
                    len=strlen(label)
                    label=strmid(label,1,len-1)
                    emiss_sel(identity).label=label
                    emiss_sel(identity).ind=indstr
                    WIDGET_CONTROL,state.labels(identity),set_value=label
                    make_strings,identity,str
                    widget_control,state.windows[identity],set_val=str
                    dens_plot,state=state
                ENDIF
            END
                                ;
            2: BEGIN
                result=emiss_select(emiss,wrange=[lo_w,hi_w],sel_ind=ind,$
                                    group=event.top)
                IF result(0) NE -1 THEN BEGIN
                    emiss_sel(identity).em=result
                    emiss_sel(identity).n_ind=N_ELEMENTS(ind)
                    IF identity EQ 0 THEN emiss_sel(identity).obs=1. $
                    ELSE emiss_sel(identity).obs=0.
                    emiss_sel(identity).sig=0.
                    WIDGET_CONTROL,state.ints(identity), $
                      set_value=STRTRIM(STRING(FORMAT='(f5.1)', $
                                               emiss_sel(identity).obs),2)
                    WIDGET_CONTROL,state.sigs(identity), $
                      set_value=STRTRIM(STRING(FORMAT='(f5.1)', $
                                               emiss_sel(identity).sig),2)
                    label=''
                    indstr=''
                    FOR i=1l,N_ELEMENTS(ind) DO BEGIN
                        label=label+'+'+STRTRIM(STRING(FORMAT='(f10.3)',$
                                                       emiss(ind(i-1)).lambda),2)
                        IF indstr EQ '' THEN indstr=trim(ind[i-1]) ELSE $
                          indstr=indstr+'+'+trim(ind[i-1])
                    ENDFOR

                    len=STRLEN(label)
                    label=STRMID(label,1,len-1)
                    emiss_sel(identity).label=label
                    emiss_sel(identity).ind=indstr
                    WIDGET_CONTROL,state.labels(identity),SET_VALUE=label
                    make_strings,identity,str
                    widget_control,state.windows[identity],set_val=str
                    dens_plot,state=state
                ENDIF
            END
                                ;
            3: BEGIN
                index_extractor,identity,index=index
                ind_main=-1
                FOR i=0,N_ELEMENTS(index)-1 DO BEGIN
                    wavel=emiss(index(0)).lambda
                    ind=WHERE(emiss.lambda LT wavel+2. AND emiss.lambda GT wavel-2.)
                    print,format='(f12.3)',emiss(ind).lambda
                ENDFOR
                print,''
            END
                                ;
        ENDCASE
                                ;
    END
;;-------------------------------------------------; LINE SELECTION

;--------------------------------------------; INPUT OBSERVED INTENSITIES
;
; Note that the text widgets are only read when the PLOT ERROR BARS 
; button on the widget is pressed. This allows intensities and sigmas 
; to be registered even if the enter key hasn't been pressed.
;
    (event.id EQ state.dn_plot[0]) OR (event.id EQ state.dn_plot[1]): BEGIN

        text_arr=[['denominator intensity','numerator intensity'], $
                  ['denominator sigma','numerator sigma']]
        tst1=0
        FOR i=0,1 DO BEGIN
            widget_control,state.ints[i],get_value=getval
            chk=valid_num(getval[0])
            IF chk EQ 0 THEN BEGIN
                tst1=tst1+1
                result=dialog_message('The '+text_arr[i,0]+' is not a valid number')
            ENDIF ELSE BEGIN
                emiss_sel[i].obs=getval[0]
            ENDELSE

            widget_control,state.sigs[i],get_value=getval
            chk=valid_num(getval[0])
            IF chk EQ 0 THEN BEGIN
                tst1=tst1+1
                result=dialog_message('The '+text_arr[i,1]+' is not a valid number')
            ENDIF ELSE BEGIN
                emiss_sel[i].sig=getval[0]
            ENDELSE
        ENDFOR 

        IF tst1 EQ 0 THEN BEGIN
            dens_plot,state=state
        ENDIF ELSE BEGIN
            emiss_sel.sig=0.
            emiss_sel.obs=0.
            dens_plot,state=state
        ENDELSE

    END


;--------------------------------------------; INPUT OBSERVED INTENSITIES

;--------------------------------------------; WAVELENGTH SLIDER
    event.id eq state.lo_slid: BEGIN
        IF event.value GE alog10(hi_w) THEN BEGIN
            lo_w=hi_w-1 
            WIDGET_CONTROL,state.lo_slid,set_value=alog10(lo_w)
        ENDIF ELSE BEGIN
            lo_w=10.^(event.value)
        ENDELSE
        IF lo_w GE 1.e6 THEN show_val = strtrim(string(format='(e9.2)',lo_w),2) $
        ELSE show_val = strtrim(string(lo_w),2)
        WIDGET_CONTROL,state.lo_txt,set_value=show_val
        wavel_plot
        dens_plot,state=state
    END

    event.id EQ state.lo_txt: BEGIN
        widget_control,state.lo_txt,get_value=value
        value=value[0]
                                ;
        IF valid_num(value) NE 0 THEN BEGIN
            value = double(value[0])
            IF value GE hi_w THEN BEGIN
                lo_w = hi_w-1 
                IF lo_w GE 1.e6 THEN show_val = $
                  strtrim(string(format='(e9.2)',lo_w),2) $
                ELSE show_val = strtrim(string(lo_w),2)
                WIDGET_CONTROL,state.lo_txt,set_value=show_val
            ENDIF ELSE BEGIN
                lo_w = value
            ENDELSE
            widget_control,state.lo_slid,set_value=alog10(lo_w)
            wavel_plot
            dens_plot,state=state
                                ;
        ENDIF ELSE BEGIN
            result=dialog_message(txt_mess)
            IF lo_w GE 1.e6 THEN show_val = $
              strtrim(string(format='(e9.2)',lo_w),2) $
            ELSE show_val = strtrim(string(lo_w),2)
            widget_control,state.lo_txt,set_value=show_val
        ENDELSE
    END

    event.id EQ state.hi_slid: BEGIN
        IF event.value LE ALOG10(lo_w) THEN BEGIN
            hi_w=lo_w+1 
            WIDGET_CONTROL,state.hi_slid,SET_VALUE=ALOG10(hi_w)
        ENDIF ELSE BEGIN 
            hi_w=10.^(event.value)
        ENDELSE
        IF hi_w GE 1.e6 THEN show_val = strtrim(string(format='(e9.2)',hi_w),2) $
        ELSE show_val = strtrim(string(hi_w),2)
        WIDGET_CONTROL,state.hi_txt,set_value=show_val
                                ;
        wset,plot_spec_id
        wavel_plot
        dens_plot,state=state
    END

    event.id EQ state.hi_txt: BEGIN
        widget_control,state.hi_txt,get_value=value
        value=value[0]
                                ;
        IF valid_num(value) NE 0 THEN BEGIN
            value = double(value[0])
            IF value LE lo_w THEN BEGIN
                hi_w = lo_w+1 
                IF hi_w GE 1.e6 THEN show_val = $
                  strtrim(string(format='(e9.2)',hi_w),2) $
                ELSE show_val = strtrim(string(hi_w),2)
                WIDGET_CONTROL,state.hi_txt,set_value=show_val
            ENDIF ELSE BEGIN
                hi_w = value
            ENDELSE
            widget_control,state.hi_slid,set_value=alog10(hi_w)
            wset,plot_spec_id
            wavel_plot
            dens_plot,state=state
                                ;
        ENDIF ELSE BEGIN
            result=dialog_message(txt_mess)
            IF hi_w GE 1e6 THEN show_val = $
              strtrim(string(format='(e9.2)',hi_w),2) $
            ELSE show_val = strtrim(string(hi_w),2)
            widget_control,state.hi_txt,set_value=show_val
        ENDELSE
    END

;--------------------------------------------; WAVELENGTH SLIDER


    event.id eq state.up_read: BEGIN ; upper limit for manual plots
        WIDGET_CONTROL,state.up_read,get_value=result
        IF is_number(result) NE 0 THEN BEGIN
            lims[1]=result
            IF set_scale EQ 0 THEN BEGIN
                set_scale=1
                widget_control,state.man_aut,set_value=1
            ENDIF
            dens_plot,state=state
        ENDIF ELSE BEGIN
            result=dialog_message(txt_mess)
            WIDGET_CONTROL,state.up_read,set_value=string(format='(e10.2)',lims[1])
        ENDELSE
    END

    event.id eq state.LO_read: BEGIN ; lower limit for manual plots
        WIDGET_CONTROL,state.lo_read,get_value=result
        IF is_number(result) NE 0 THEN BEGIN
            lims(0)=result
            IF set_scale EQ 0 THEN BEGIN
                set_scale=1
                widget_control,state.man_aut,set_value=1
            ENDIF
            dens_plot,state=state
        ENDIF ELSE BEGIN
            result=dialog_message(txt_mess)
            WIDGET_CONTROL,state.lo_read,set_value=string(format='(e10.2)',lims[0])
        ENDELSE
    END

;GDZ-
    event.id eq state.log_plot_info: BEGIN ; switch between log or linear
        CASE event.value OF     ;  scaling
                                ;
            0: BEGIN
                loglin = 1
;in the log case, check that we do not have automatic scaling
; with zero.
  if  set_scale eq 0 then begin 
;switch to auto /ynozero 
 set_scale=2
widget_control,state.man_aut,set_value=1
               widget_control,state.lim_base,map=0

if lims[0] le 0 then begin 
lims[0]=1e-10
WIDGET_CONTROL,state.lo_read,set_value=string(format='(e10.2)',lims[0])
end

  endif else if set_scale eq 2 then   begin 

if lims[0] le 0 then begin 
lims[0]=1e-10
WIDGET_CONTROL,state.lo_read,set_value=string(format='(e10.2)',lims[0])
end 

  endif else if set_scale eq 1 then   begin   ;manual
;widget_control,state.man_aut,set_value=1
;avoid zeros:

if lims[0] le 0 then begin 
lims[0]=1e-10
WIDGET_CONTROL,state.lo_read,set_value=string(format='(e10.2)',lims[0])
endif 

  end 
            END
                                ;
            1: BEGIN
                loglin = 0
            END
        ENDCASE  
;refresh:
        dens_plot,state=state
    END  

    event.id eq state.man_aut: BEGIN ; switch between manual and
        CASE event.value OF     ; automatic scaling
                                ;
            0: BEGIN
;GDZ - do not allow if log option is set.
              if not loglin then begin 
                set_scale=0   ;automatic scaling
                widget_control,state.lim_base,map=0
            endif else begin 
set_scale=2   ;automatic scaling , /ynozero
widget_control,state.man_aut,set_value=1
widget_control,state.lim_base,map=0

            end 
            END
                                ;
            1: BEGIN
                set_scale=2   ;automatic scaling , /ynozero
                widget_control,state.lim_base,map=0
            END
                                ;
            2: BEGIN
                set_scale=1  ;manual scaling
                widget_control,state.lim_base,map=1
            END
                                ;
        ENDCASE
        dens_plot,state=state
    END

    event.id EQ state.refresh_butt: dens_plot,state=state

    event.id EQ state.pr_butt: BEGIN

        nd=ROUND((hi_dens-lo_dens)/dint +1)
        pr_string=strarr(nd)
        ratio=emiss_sel[1].em/emiss_sel[0].em
        FOR i=0,nd-1 DO BEGIN
            format='(f12.3)'
            IF ratio[i] LT 0.01 THEN format='(e12.3)'
            IF ratio[i] GE 1d7 THEN format='(e12.3)'
            pr_string[i]=string(format='(f7.2)',dens[i])+ $
              string(format=format,ratio[i])
        ENDFOR

;GDZ
;        zion2filename,iz,ion,filename, name=name, diel=diel
    zion2spectroscopic,iz,ion,ion_name, diel=diel
        ion_name=strcompress(ion_name)

        IF t_switch EQ 1 THEN xtit='Log Electron temperature [K] ' $
        ELSE xtit='Log Electron density [cm-3] '

        IF units EQ 0 THEN YTIT='   Ratio (Energy) values'  ELSE IF $
          units EQ 1 THEN YTIT='   Ratio (Photons) values'

        pr_string=['Ratio values for '+ion_name, $
                   emiss_sel[1].label+'/'+emiss_sel[0].label, $
                   '', $
                   xtit+ytit, $
                   '', $
                   pr_string, $
                   '', $
                   'Mean: '+string(format=format,total(ratio)/n_elements(ratio)), $
                   'Std Dev: '+string(format=format,stdev(ratio))]
;    result=dialog_message(pr_string,/info)
        xpopup, pr_string

    END

    event.id EQ state.b_all: BEGIN

          CASE event.value OF
                                  ;
              0: BEGIN
                  all=0
                  widget_control,state.pr_butt,sens=1
              END
              1: BEGIN
                  all=1
                  widget_control,state.pr_butt,sens=0
              END
                                  ;
          ENDCASE

      dens_plot,state=state

    END

;---------------------------------------buttons at bottom of widget
    event.id eq state.extras: BEGIN
        CASE 1 OF
                                ;
;GDZ
            event.value EQ 1:BEGIN 

;save the values in a file
;-------------------------

                index=WHERE(emiss.lambda GE lo_w AND emiss.lambda LE hi_w)
                IF index(0) EQ -1 THEN BEGIN ;---no lines in specified range
                    result=WIDGET_MESSAGE('No lines in specified range',/info)
                ENDIF ELSE BEGIN 

                    cd,curr=curr
                    result=dialog_pickfile(file=concat_dir(curr,'_emiss.genx'),filter='*.genx', $
                                           tit='Choose savegen file ')
                    IF result[0] EQ '' THEN GOTO,lbl1

                    make_currplot_label, label, struct_info= struct_info

                    struct_save = join_struct({created:systime(), Directory:fpath, emiss:emiss(index)}, $
                                              struct_info)

                    savegen, file=result[0], struct=struct_save

                    result=WIDGET_MESSAGE('Values saved in IDL genx file '+result[0]+$
                                          ' - to restore, use restgen ' , /info) 

;help, struct_save, /st


                END 
            END 

            event.value EQ 3 OR event.value EQ 4 OR event.value EQ 5: BEGIN
                cd,curr=curr
                IF event.value EQ 3 OR event.value EQ 5 THEN BEGIN
                    filename=concat_dir(curr,'idl.ps')
                ENDIF ELSE IF event.value EQ 4  THEN BEGIN
                    result=dialog_pickfile(file=concat_dir(curr,'idl.ps'),filter='*.ps')
                    IF result[0] EQ '' THEN GOTO,lbl1
                    filename=result[0]
                END 

                                ;
                dname=!d.name
                !p.font=0
                SET_PLOT,'ps'
                device,xsiz=7,ysiz=6,/inches,/isolatin1,file=filename,/portrait, $
                  yoff=4.0
                dens_plot,state=state,ps=1
                                ;
                ypos=lims[0]-(lims[1]-lims[0])*0.20
                xpos=!x.crange[0]
                make_currplot_label,label
                label=strjoin(label,'!c')

;        xyouts,xpos,ypos,label
                                ;
                IF emiss_sel[1].obs NE 0. THEN BEGIN
                    xpos=0.6*!x.crange[1]+0.4*!x.crange[0]
                    make_dens_label,label
                    label=strjoin(label,'!c')
;          xyouts,xpos,ypos,label
                ENDIF
                                ;
                fname=concat_dir(!xuvtop,'VERSION')
                str1=''
                IF file_exist(fname) THEN BEGIN
                    openr,lun,fname,/get_lun
                    readf,lun,str1
                    free_lun,lun
                ENDIF ELSE BEGIN
                    str1='Please update your version of CHIANTI to v.4 or later'
                ENDELSE

;        xyouts,!x.crange[0],lims[0]-(lims[1]-lims[0])*0.40, $
;             'Plot created on '+systime()+'!c'+$
;             'CHIANTI version: '+str1+'!c'+$
;             'File: '+filename

;GDZ- create a list of lines to be printed in the lower part of the
;     plot. 

                label=label+'!c  !c'+'Plot created on '+systime()+'!c'+$
                  'CHIANTI version: '+str1+'!c'+$
                  'File: '+filename+'!c'+$
                  'Directory: '+fpath+' !c '

IF all eq 0 THEN BEGIN         ;;----------------; PLOT SELECTED LINES

;? 


endif else begin 

                FOR i=0L, n_elements(line_labels) -1 DO $
                  label=label+'!c'+line_labels(i)
            end 
                
                xyouts,0.05, -0.08,/norm, label, chars=0.9
                                ;
                DEVICE,/CLOSE
                !p.font=-1
                SET_PLOT,dname

;          ps,/portrait
;          dens_plot,state=state
                                ;
                                ; add text info below plot
                                ;
;          psclose
;          !p.font=-1

                IF  event.value EQ 3 OR event.value EQ 4 THEN $
                  result=WIDGET_MESSAGE('Plot sent to '+$
                                        filename,$
                                        /info) ELSE  IF event.value EQ 5 THEN BEGIN 

                    prin = getenv('PRINTER')
                    IF prin EQ '' THEN $
                      xsel_printer,prin,group=DENS_MAIN,$
                      instruct='Select printer for hardcopy'

                    psplot,filename,queu=prin
                    result=WIDGET_MESSAGE('Plot sent to idl.ps file and to PRINTER '+prin, /info) 
                END 
                lbl1:

            END  

            event.value EQ 6: WIDGET_CONTROL, event.top, /DESTROY ; quit

        ENDCASE 
    END

    ELSE:

ENDCASE

    END 




;-----------------------------------------------------------------------------
PRO sample_wid, GROUP=Group
;-----------------------------------------------------------------------------
;GDZ - changed plotting and emiss_data COMMONS
COMMON emiss_data,iz,ion,diel,emiss,lo_dens,hi_dens,dint,dens,temp,units
COMMON select,emiss_sel
COMMON sliders,lo_w,hi_w,min_w,max_w
COMMON plotting,lims,set_scale,all, loglin, line_labels
COMMON Comm, plot_rat_id, plot_spec_id
COMMON proton_rates, prot_nincl,prot_switch,prot_mess,pname
COMMON proton, pstr, pe_ratio
COMMON rad_data, radtemp,rphot,r_tst
COMMON extra, fpath
COMMON elvlc, l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,ref
COMMON files, ioneq_file, abund_file
COMMON temp, t_switch

IF N_ELEMENTS(Group) EQ 0 THEN GROUP=0


;--------------------------------O
; Create a title for the widget
;

;GDZ - allow dielectronic files
zion2spectroscopic,iz,ion,title, diel=diel

;fpath must be defined.
direc=fpath 

;IF KEYWORD_SET(fpath) THEN direc=fpath ELSE $
; zion2filename,iz,ion, filename, name=name, diel=diel

; zion2filename,iz,ion,direc

title='ION: '+title+' --- DIRECTORY: '+direc
;--------------------------------O


; Main base for everything

;GDZ- added scrolling

device, get_screen_size = sz
IF sz[0] LT 1200 THEN x_scroll=sz[0]*0.9 ELSE x_scroll=0
IF sz[1] LT 1000 THEN y_scroll=sz[1]*0.9 ELSE y_scroll=0

;sz(0) = 1100 < 0.9*sz(0)
;sz(1) = 900 < 0.9*sz(1)

;chianti_font,font

DENS_MAIN = WIDGET_BASE(GROUP_LEADER=Group, $
      col=1, MAP=1, x_scroll=x_scroll,y_scroll=y_scroll, $
;      col=1, MAP=1,  $
      UVALUE='DENS_MAIN',title=title)

fname=concat_dir(!xuvtop,'VERSION')
str1=''
IF file_exist(fname) THEN BEGIN
  openr,lun,fname,/get_lun
  readf,lun,str1
  free_lun,lun
  str1='CHIANTI version '+str1
ENDIF ELSE BEGIN
  str1='Please update your version of CHIANTI to v.4 or later'
ENDELSE
IF t_SWITCH EQ 1 THEN rout='TEMP_PLOTTER' ELSE rout='DENS_PLOTTER'
label=str1+' --- ROUTINE: '+rout

chianti_font,bfont,/big

chianti_txt1=widget_label(dens_main,val=label,/align_left, $
                         font=bfont)


BASE_1=WIDGET_BASE(DENS_MAIN,col=1,map=1,uvalue='BASE_1')


;---------------------------------------------------[]
; base_1 will contain emissivity stuff. Everything else is on base_2
;
BASE_SUB1=WIDGET_BASE(BASE_1,row=1,map=1,uvalue='BASE_SUB1', $
                     /align_center,sensitive=0)

base_emiss = widget_base(base_sub1,col=1,/frame)
units_lbl=widget_label(base_emiss,/align_left,value='Emissivity units')
units_buts = cw_bgroup(base_emiss,['Energy','Photons'],/row, $
                           set_value=0,/exclusive)


PROT_BASE=WIDGET_BASE(BASE_SUB1,col=1,map=1,uvalue='PROT_BASE',frame=2)
;
PROT_LBL=WIDGET_LABEL(PROT_BASE,/align_left,value='Proton rates')
PROT_BGRP=CW_BGROUP(PROT_BASE,['yes','no'],/row,$
                        set_value=prot_nincl,/exclusive,$
                        uvalue=[0,1])
PROT_DISP=WIDGET_TEXT(PROT_BASE,value=prot_mess(prot_switch),xsiz=11)


B_D_SLID=WIDGET_BASE(BASE_SUB1,col=1,map=1,uvalue='B_D_SLID', frame=2)
;
IF t_switch EQ 1 THEN str1='Temperature range, log T [K]' $
     ELSE str1='Density range, log Ne [cm-3]'
D_SLID_TXT=WIDGET_LABEL(B_D_SLID,/align_left,value=str1)
IF t_switch EQ 1 THEN BEGIN
  ld_i=fix((lo_dens-2.0)*5.)
  hd_i=fix((hi_dens-2.0)*5.)
  LD_LABEL=widget_label(b_d_slid,xsiz=40, $
                        val=string(format='(f4.1)',lo_dens))
  LO_DENS_SLID=WIDGET_SLIDER(B_D_SLID,min=0,max=40,val=ld_i,xsiz=200, $
                             /suppress_value)
  HD_LABEL=widget_label(b_d_slid,xsiz=40, $
                        val=string(format='(f4.1)',hi_dens))
  HI_DENS_SLID=WIDGET_SLIDER(B_D_SLID,min=0,max=40,val=hd_i,xsiz=200, $
                             /suppress_value)
ENDIF ELSE BEGIN
  ld_label=0
  hd_label=0
  LO_DENS_SLID=WIDGET_SLIDER(B_D_SLID,min=1,max=16,val=lo_dens,xsiz=200)
  HI_DENS_SLID=WIDGET_SLIDER(B_D_SLID,min=1,max=16,val=hi_dens,xsiz=200)
ENDELSE

TEMP_DENS=WIDGET_BASE(BASE_SUB1,col=1,map=1,uvalue='TEMP_DENS')
;
; Note: the set_value below selects the second button in the group by default
;
DENS_INT=WIDGET_BASE(TEMP_DENS,col=1,map=1,uvalue='DENS_INT',frame=2)
IF t_switch EQ 1 THEN str1='Temperature intervals' $
ELSE str1='Density intervals'
DENS_INT_TXT=WIDGET_LABEL(DENS_INT,/align_left,value=str1)
IF t_switch EQ 1 THEN str1=['0.2','0.1','0.05','0.01'] $
     ELSE str1=['1.0','0.5','0.2','0.1']
IF t_switch EQ 1 THEN setval=1 ELSE setval=2
DENS_BGRP=CW_BGROUP(DENS_INT,str1,/row,$
                        set_value=setval,/exclusive,$
                        uvalue=str1)

TEMP_SET=WIDGET_BASE(TEMP_DENS,row=1,map=1,uvalue='TEMP_SET',frame=2)
;
IF t_switch EQ 1 THEN str1='Dens, log N_e [cm-3]:' $
ELSE str1='Temp, log T [K]:'
TEMP_LBL=WIDGET_LABEL(TEMP_SET,/align_left, $
               value=str1)
TEMP_READ=WIDGET_TEXT(TEMP_SET,$
                          value=strtrim(string(format='(f5.2)',temp),2),$
                          xsiz=5,/editable)


RAD_SET=WIDGET_BASE(BASE_SUB1,col=1,map=1,frame=2)
;;

pexc_lbl=widget_label(rad_set,/align_left,value='Include photoexcitation?')
pexc_BGRP=CW_BGROUP(rad_set,['yes','no'],/row,$
                        set_value=1,/exclusive,$
                        uvalue=[0,1])

DIL_SET=WIDGET_BASE(RAD_SET,row=1,map=0)
;
DIL_LBL=WIDGET_LABEL(DIL_SET,/align_left,value='Distance (R_* units):')
DIL_READ=WIDGET_TEXT(DIL_SET,$
                     value=strtrim(string(format='(f4.2)',rphot),2),$
                     xsiz=5,/editable)

RT_SET=WIDGET_BASE(RAD_SET,row=1,map=0)
;
RT_LBL=WIDGET_LABEL(RT_SET,/align_left,value='Radiation Temp (K):')
RT_READ=WIDGET_TEXT(RT_SET,$
                          value=strtrim(string(format='(f7.1)',radtemp),2),$
                          xsiz=8,/editable)

;
; the following is the button which makes the parameter change widgets 
; appear
;
ch_params=WIDGET_BUTTON(BASE_1,val='CHANGE PARAMETERS',xsiz=250,uval=0, $
                       /align_center)


;---------------------------------------------------[]
; BASE_2 contains all other widgets
;
BASE_2=WIDGET_BASE(DENS_MAIN,row=1,map=1,uvalue='BASE_2')
;
;
; Base to contain denominator and numerator bases
;
BASE_ND = WIDGET_BASE(BASE_2, COLUMN=1, MAP=1, xsiz=450, UVALUE='BASE_ND')


junk   = { CW_PDMENU_S, flags:0, name:'' }
;
desc=[ { CW_PDMENU_S, 1, 'Choose a new line' }, $
         { CW_PDMENU_S, 0, '...from strongest lines' }, $
         { CW_PDMENU_S, 2, '...from specified wavelength range' } ];, $
;         { CW_PDMENU_S, 2, 'Show blends' } ]
;
;---------------------------------------<>
; Label and buttons for numerator
;
N1_BASE=WIDGET_BASE(BASE_ND, FRAME=2, $
      COLUMN=1, MAP=1, xsiz=500, UVALUE='N1_BASE')

N1_BASE_1=WIDGET_BASE(N1_BASE, row=1, MAP=1, xsiz=500, UVALUE='D_BASE_1')
;
N1_TXT2=WIDGET_LABEL(N1_BASE_1,/align_left,value='NUMERATOR:  ')
N1_TXT =WIDGET_LABEL(N1_BASE_1,/align_left,xsiz=400,VALUE=emiss_sel(1).label)

N1_LINES=CW_PDMENU(N1_BASE,desc,uvalue=1,font=font)

N1_BASE_2=WIDGET_BASE(N1_BASE, row=1, MAP=1, xsiz=500, UVALUE='N1_BASE_1')
;
L_INT_N1=WIDGET_LABEL(N1_BASE_2,/align_center,value='INTENSITY:')
INT_N1=WIDGET_TEXT(N1_BASE_2,/EDITABLE, xsiz=7, uvalue=1, $
                    value=strtrim(string(format='(f6.2)', $
                                         emiss_sel(1).obs),2))
L_SIG_N1=WIDGET_LABEL(N1_BASE_2,/align_center,value='SIGMA:')
SIG_N1=WIDGET_TEXT(N1_BASE_2,/EDITABLE, xsiz=7, uvalue=1, $
                    value=strtrim(string(format='(f6.2)', $
                                         emiss_sel(1).obs),2))
n1_plot=widget_button(n1_base_2,val='PLOT ERROR BARS',uvalue=1)
;
make_strings,1,str
;GDZ
dummy_base=WIDGET_BASE(n1_base, /col, xsize=430)
n1_window = WIDGET_TEXT(dummy_base,ysiz=3, $
                         value=str,/scroll,/wrap)
;---------------------------------------<>

;---------------------------------------[]
; Label and buttons for denominator
;
D_BASE=WIDGET_BASE(BASE_ND, FRAME=2, $
                       COLUMN=1, MAP=1, xsiz=500, UVALUE='D_BASE')

D_BASE_1=WIDGET_BASE(D_BASE, row=1, MAP=1, xsiz=500, UVALUE='D_BASE_1')
;
D_TXT2=WIDGET_LABEL(D_BASE_1,/align_left,value='DENOMINATOR:  ')
D_TXT =WIDGET_LABEL(D_BASE_1,/align_left,xsiz=400,VALUE=emiss_sel[0].label)

D_LINES=CW_PDMENU(D_BASE,desc,uvalue=0,font=font)

D_BASE_2=WIDGET_BASE(D_BASE, row=1, MAP=1, xsiz=500, UVALUE='D_BASE_1')
;
L_INT_D=WIDGET_LABEL(D_BASE_2,/align_center,value='INTENSITY:')
INT_D=WIDGET_TEXT(D_BASE_2,/EDITABLE, xsiz=7, uvalue=0, $
                     value=strtrim(string(format='(f6.2)', $
                                          emiss_sel(0).obs),2), $
                  /all_events)
L_SIG_D=WIDGET_LABEL(D_BASE_2,/align_center,value='SIGMA:')
SIG_D=WIDGET_TEXT(D_BASE_2,/EDITABLE, xsiz=7, uvalue=0, $
                     value=strtrim(string(format='(f6.2)', $
                                          emiss_sel(0).sig),2), $
                  /all_events)
d_plot=widget_button(d_base_2,val='PLOT ERROR BARS',uval=0)
make_strings,0,str
;GDZ
dummy_base=WIDGET_BASE(d_base, /col, xsize=430)
den_window = WIDGET_TEXT(dummy_base,ysiz=3, $
                         value=str,/scroll,/wrap)
;---------------------------------------[]

;----------------------------------(I)
slider_base=widget_base(base_nd,/col,/frame)
;
LO_SLID_B=WIDGET_BASE(slider_base, ROW=1, MAP=1, xsiz=400,$
                         /align_center,UVALUE='LO_SLID_B')
;
LO_SLID=cw_fslider(LO_SLID_B,min=min_w,max=max_w,xsiz=300,$
      tit='Low wavelength limit ('+string(197b)+')', $
                   val=alog10(lo_w),/suppress_value)
IF lo_w GE 1.e6 THEN show_val = strtrim(string(format='(e9.2)',lo_w),2) $
                ELSE show_val = strtrim(string(lo_w),2)
LO_TXT=WIDGET_TEXT(LO_SLID_B,xsiz=9,$
      /align_left,VALUE=show_val,/edit)

HI_SLID_B=WIDGET_BASE(slider_base,ROW=1,MAP=1, xsiz=400,/align_center, $
                          UVALUE='HI_SLID_B')
;
HI_SLID=cw_fslider(HI_SLID_B,min=min_w,max=max_w,xsiz=300,$
      tit='High wavelength limit ('+string(197b)+')', val=alog10(hi_w),/suppress_value)
IF hi_w GE 1.e6 THEN show_val = strtrim(string(format='(e9.2)',hi_w),2) $
                ELSE show_val = strtrim(string(hi_w),2)
HI_TXT=WIDGET_TEXT(HI_SLID_B,xsiz=9,$
      /align_left,VALUE=show_val,/edit)
;----------------------------------(I)

;GDZ retain=2

;P_S_BASE=WIDGET_BASE(slider_base,row=1,map=1,xsiz=380,/align_center)
PLOT_SPEC=WIDGET_DRAW(slider_base, $
      RETAIN=2, xsiz=350,ysiz=30, $
      UVALUE='PLOT_RAT',/align_center)

;---------------------------------------0
; Quit button
;
IF t_SWITCH EQ 1 THEN tede_str='temperatures' ELSE tede_str='densities'
; EXTRAS = CW_BGROUP( BASE_ND, $
;                       ['Refresh plot','Show derived '+tede_str,'Hardcopy',$
;                        'QUIT'],$
;                       /row)
;junk   = { CW_PDMENU_S, flags:0, name:'' }
desc=[ { CW_PDMENU_S, 1, 'SAVE' }, $
         { CW_PDMENU_S, 2, 'emiss. of selected lines' }, $
{ CW_PDMENU_S, 1, 'HARDCOPY' }, $
         { CW_PDMENU_S, 0, 'send to idl.ps' }, $
         { CW_PDMENU_S, 0, 'specify filename' }, $
         { CW_PDMENU_S, 2, 'send to printer' }, $
         { CW_PDMENU_S, 2, 'QUIT' } ]
extras=CW_PDMENU(base_nd,desc,uvalue=1,font=bfont)
;EXTRAS = CW_BGROUP( BASE_ND, ['Hardcopy','QUIT'],/row,font=bfont)
;---------------------------------------0

help_txt1=widget_label(base_nd,/align_left,val=' ')
val='Send comments to chianti_help@halcyon.nrl.navy.mil'
help_txt2=widget_label(base_nd,/align_left,val=val,font=font)


;---------------------------------------<>
; Base for plotting window
;
WIND_BASE=WIDGET_BASE(BASE_2,col=1,MAP=1,UVALUE='WIND_BASE')

;
; create widgets containing info 
;
info_base=widget_base(wind_base,/row)

;GDZ- added various things:

IF t_switch EQ 1 THEN xtit='Log Electron temperature [K] :' $
ELSE xtit='Log Electron density [cm-3] :'

IF units EQ 0 THEN YTIT='Ratio (Energy) values: '  ELSE IF $
    units EQ 1 THEN YTIT='Ratio (Photons) values: '

dummy1 = XTIT+arr2str(trim(float(dens)), ',',/trim)
dummy2 = YTIT+arr2str(trim(float(emiss_sel(1).em/emiss_sel(0).em)), ',',/trim)

make_currplot_label,label
;GDZ
;info_txt=widget_text(info_base,ysiz=6,xsiz=20,value=[dummy1, dummy2, label], /scroll)
info_txt=widget_text(info_base,ysiz=6,xsiz=25,value=label,/wrap)

;info_txt=widget_text(info_base,ysiz=6,xsiz=32,value=label)
;
make_dens_label,label
dens_txt=widget_text(info_base,ysiz=6,xsiz=32,value=label,/wrap)
                     
;
; Plotting window
;
PLOT_RAT = WIDGET_DRAW( WIND_BASE, $
                        RETAIN=2, $
                        UVALUE='PLOT_RAT', $
                        XSIZE=520, $
                        YSIZE=400)


;----------------------X
; Base for plotting widgets
;
PLOT_BASE= WIDGET_BASE(WIND_BASE,$
                       ROW=1, MAP=1,  UVALUE='PLOT_BASE')


; Choose between selected lines and all lines. Default=0
;
plot_b_base = widget_base(plot_base,/col,/base_align_center)
B_ALL = CW_BGROUP( plot_b_BASE,['Selected lines','All lines'],/col,$
                       set_value=0,/exclusive,/frame)
pr_butt=widget_button(plot_b_base,val='SHOW RATIO VALUES',sens=1)
refresh_butt=widget_button(plot_b_base,val='REFRESH PLOT')

;
; Choose between selected lines and all lines. Default=0
;
yaxis_base=widget_base(plot_base,/row,/frame)
;
man_aut_base=widget_base(yaxis_base,/col)
man_aut_txt=widget_label(man_aut_base,/align_left,val='Y-AXIS SCALING')
man_aut_base2=widget_base(man_aut_base,/row)
;
log_plot_info = CW_BGROUP(man_aut_base2,['Log', 'Linear'], $
                  /col,set_value=1,/exclusive)
loglin=0
;
MAN_AUT=CW_BGROUP(man_aut_base2,['Automatic','Automatic (ynozero)','Manual'], $
                  /col,set_value=0,/exclusive)

;
; Base for containing plot upper and lower limits
;
LIM_BASE=WIDGET_BASE(yaxis_base,$
                     col=1, MAP=0, UVALUE='LIM_BUT_B')
;
UP_READ_L=WIDGET_LABEL(LIM_BASE,value='Upper limit:',/align_left)
UP_READ=WIDGET_TEXT(LIM_BASE,$
                    value=strtrim(string(format='(e10.2)',lims(1)),2),$
                    xsiz=8,/editable)
;
LO_READ_L=WIDGET_LABEL(LIM_BASE,value='Lower limit:',/align_left)
LO_READ=WIDGET_TEXT(LIM_BASE,$
                    value=strtrim(string(format='(e10.2)',lims(0)),2),$
                    xsiz=8,/editable)

; pr_base=widget_base(plot_base,/frame)
; pr_butt=widget_button(pr_base,val='SHOW RATIO VALUES')


;---------------------------------------<>

;; To add more numerators, add an extra n?_bgroup and n?_txt to `state'
;
state={ch_params: ch_params, base_sub1: base_sub1, base_2: base_2, $
       prot_bgrp:prot_bgrp, prot_disp:prot_disp, $
       lo_dens_slid:lo_dens_slid, hi_dens_slid:hi_dens_slid, $
       dens_bgrp:dens_bgrp, temp_read:temp_read, $
       pexc_bgrp:pexc_bgrp, dil_set:dil_set, rt_set:rt_set, $
       dil_read:dil_read, rt_read:rt_read, $
       d_lines:d_lines, n1_lines:n1_lines, $
       labels:[d_txt,n1_txt], $
       windows:[den_window,n1_window], $
       dn_plot: [d_plot,n1_plot], $
       ints:[int_d,int_n1], sigs:[sig_d,sig_n1], $
       lo_txt:lo_txt, lo_slid:lo_slid, $
       hi_txt:hi_txt, hi_slid:hi_slid, lim_base:lim_base, $
       up_read:up_read, lo_read:lo_read, man_aut:man_aut, $
       b_all:b_all, units_buts:units_buts, $
       extras:extras, ld_label: ld_label, hd_label: hd_label, $
       info_txt:info_txt, dens_txt:dens_txt, $
       refresh_butt:refresh_butt, pr_butt: pr_butt,log_plot_info:log_plot_info}

  WIDGET_CONTROL, DENS_MAIN, /REALIZE, set_uvalue=state

  ; Get drawable window index

  WIDGET_CONTROL, plot_rat, GET_VALUE=plot_rat_id
  WIDGET_CONTROL, plot_spec, GET_VALUE=plot_spec_id

  wavel_plot
  dens_plot,state=state

  XMANAGER, 'DENS_MAIN', DENS_MAIN, group=group

END


;-----------------------------------------------------------------------------
PRO RATIO_PLOTTER, ION_Z, ION_SP, EM, PATH=PATH, NOPROT=NOPROT, $
                  IONEQ_FILE=IONEQ_FILE, ABUND_FILE=ABUND_FILE, $
                  TEMPERATURE=TEMPERATURE, DENSITY=DENSITY,DIELECTRONIC=DIELECTRONIC

;-----------------------------------------------------------------------------
;GDZ - changed plotting and emiss_data COMMONS
COMMON emiss_data,iz,ion,diel,emiss,lo_dens,hi_dens,dint,dens,temp,units
COMMON select,emiss_sel
COMMON sliders,lo_w,hi_w,min_w,max_w
COMMON plotting,lims,set_scale,all, loglin, line_labels
COMMON rad_data, radtemp,rphot,r_tst
COMMON proton_rates, prot_nincl,prot_switch,prot_mess,pname
COMMON proton, pstr, pe_ratio
COMMON extra, fpath
COMMON elvlc, l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,ref
COMMON files, ioneq_name, abund_name
COMMON temp, t_switch


IF N_PARAMS() LT 2 THEN BEGIN
  PRINT,'Use: IDL> ratio_plotter, iz, ion [, em, path=path, /noprot, $'
  print,'               ioneq_file=ioneq_file, abund_file=abund_file, $'
  print,'               /temperature, /density ] '
  print,''
  print,'** to look at density sensitive ratios, use /density'
  print,'** to look at temperature sensitive ratios, use /temperature'
  RETURN
ENDIF

;GDZ - clean up memory in COMMON blocks:

delvarx, iz, ion, emiss,lo_dens,hi_dens,dint,dens,temp,units
delvarx, emiss_sel, lo_w,hi_w,min_w,max_w, lims,set_scale,all, loglin
delvarx, radtemp,rphot,r_tst, prot_nincl,prot_switch,prot_mess,pname, pstr, pe_ratio
delvarx, fpath, l1,term,conf,ss,ll,jj,ecm,eryd
delvarx, ecmth,erydth,ref,t_switch 
delvarx, ioneq_name, abund_name

IF n_elements(DIELECTRONIC) gt 0 then diel=DIELECTRONIC else diel=0

IF keyword_set(temperature) THEN t_switch=1 ELSE t_switch=0

; The yellow bars on the bottom left display did not work on my 24-bit 
; display, so the following corrects for this. Specifying col=150 gives 
; yellow.
;
device,decomposed=0

tvlct,r,g,b,/get
loadct,0
tvlct,[[255],[255],[0]],50


chianti_font,font

widget_control,default_font=font

IF n_elements(abund_file) NE 0 THEN abund_name=abund_file

;------------------<
; If, in a previous call to dens_plotter, you have set path, then it 
; is `remembered' in the common block. The following command removes it 
; from the common block
;
IF N_ELEMENTS(fpath) NE 0 THEN junk=temporary(fpath)
;------------------<

IF N_ELEMENTS(path) NE 0 THEN begin 
  fpath=path 
  zion2filename,ion_z,ion_sp, filename,name=name, diel=diel
 ;
 ; pname is the complete path specification of where the .psplups file
 ; should be 
 ;
  pname=concat_dir(path, name)+'.psplups'
ENDIF ELSE BEGIN 
 zion2filename,ion_z,ion_sp,pname,name=name, diel=diel
 pname=pname+'.psplups'
 pp = str_sep(name, '_')
 fpath=concat_dir(concat_dir(!xuvtop, pp[0]), name)
END 

iz=ion_z & ion=ion_sp

zion2spectroscopic,iz,ion,ion_name,diel=diel


;-------------------------------<I>
; Proton rates.  Note the difference between prot_nincl and prot_switch. 
; The former is either 0 or 1 and shows whether proton rates are included 
; in the level balance (=0) or not (=1). The latter is meant for use with 
; prot_mess 
; and can be either 0, 1 or 2. In the case of 2, it means that the routine 
; has tried to find the .psplups file, but failed.
;
prot_mess=['Not included','Included','No .psplups file']

IF NOT KEYWORD_SET(noprot) THEN BEGIN
  result=FINDFILE(EXPAND_PATH(pname))
  IF result(0) NE '' THEN BEGIN
    prot_nincl=0
    prot_switch=1
  ENDIF ELSE BEGIN
    prot_nincl=1
    prot_switch=2
  ENDELSE
ENDIF ELSE BEGIN
  prot_nincl=1
  prot_switch=0
ENDELSE
;-------------------------------<I>


units=0                          ; ratios in energy units
lims=[0.,0.] & set_scale=0       ; default: automatic scaling
all=0                            ; default: display selected lines


radtemp=6000. & rphot=1.0 & r_tst=0   ; default: radiation data

;------------------------------------+
; Work out the T_max for the ion using the ion balance data in CHIANTI. 
;
dir=concat_dir(!xuvtop,'ioneq')        
IF n_elements(ioneq_file) NE 0 THEN ioneq_name=ioneq_file $
     ELSE ioneq_name=!ioneq_file

read_ioneq,ioneq_name,temp_all,ioneq,ioneq_ref

; GDZ

f_all=ioneq(*,iz-1,ion-1+diel)

IF total(f_all) EQ 0. THEN BEGIN
  print,''
  print,'** Ion fraction data is not available for this ion. Please select'
  print,'** another ion balance file with the IONEQ_FILE keyword.'
  print,''
  print,'Current ion balance file: ',ioneq_name
  return
ENDIF
;
ind=WHERE(f_all EQ max(f_all))
ind=ind(0)      ; in case ind has more than one element (e.g., Ar IX)
temp=temp_all(ind)
;------------------------------------+

IF t_switch EQ 1 THEN BEGIN
; GDZ
;  lo_dens=temp-1.0 & hi_dens=temp+1.0 & dint=0.2
  lo_dens=temp-0.5 & hi_dens=temp+0.5 & dint=0.1
  dens=findgen(11)*dint + lo_dens
  temp=10.0
ENDIF ELSE BEGIN
  lo_dens=8 & hi_dens=12 & dint=0.2
  n=fix((hi_dens-lo_dens)/dint)+1
  dens=findgen(n)*dint + lo_dens     ; default: density range
ENDELSE


PRINT,''
PRINT,' -- Please wait while emissivities are calculated --'
PRINT,''

IF t_switch EQ 1 THEN BEGIN
  emiss=emiss_calc(iz,ion,temp=dens,dens=temp, $
                   path=fpath,noprot=noprot, $
                   /quiet, ioneq_file=ioneq_name,DIEL=DIEL, $
                   abund_file=abund_file)
ENDIF ELSE BEGIN
  emiss=emiss_calc(iz,ion,temp=temp,dens=dens,path=fpath,noprot=noprot, $
                   /quiet, ioneq_file=ioneq_name,DIEL=DIEL, abund_file=abund_file)
ENDELSE


lo_w=MIN(emiss.lambda-1) & hi_w=MAX(emiss.lambda+1)
min_w=ALOG10(lo_w) & max_w=ALOG10(hi_w)

;;--------------------------------------------------<
;; Create the `emiss_sel' structure from emiss
;
; emiss_sel is a structure similar to emiss that contains the
; emissivities of the lines being plotted. emiss_sel(0) is always the
; denominator.
; It is not possible to directly store the indices of blends in `emiss_sel'
; and so I'm going to store the number of indices and then use a little
; routine to extract the indices from emiss_sel.label (the routine is 
; index_extractor)
;
nd=n_elements(dens)
str={label:'', em:FLTARR(nd), n_ind:-1, obs:0., sig:0., $
          dens:0. , densup:0. , denslo:0., ind: ''}
emiss_sel=REPLICATE(str,5)       ; allow for up to 4 numerators
;
ind=WHERE(emiss.em(0) EQ MAX(emiss.em(0)))
ind=ind[0]
emiss_sel(0).label=STRTRIM(STRING(FORMAT='(f10.3)',emiss(ind).lambda),2)
emiss_sel(0).em=emiss(ind).em
emiss_sel(0).n_ind=1
emiss_sel(0).obs=1.0 & emiss_sel(0).sig=0.
emiss_sel[0].ind=trim(ind)
;
ind=WHERE(reform(emiss.em[nd-1]) EQ MAX(reform(emiss.em[nd-1])))
ind=ind[0]
emiss_sel(1).label=STRTRIM(STRING(FORMAT='(f10.3)',emiss(ind).lambda),2)
emiss_sel(1).em=emiss(ind).em
emiss_sel(1).n_ind=1
emiss_sel(1).obs=0. & emiss_sel(1).sig=0.
emiss_sel[1].ind=trim(ind)
;;----------------------------------------------------<


;;----------------------------[O]
; To display transition information, I need to read the .elvlc file
; 
IF KEYWORD_SET(path) THEN BEGIN
zion2filename,iz,ion,filename,name=name, diel=diel
ename=concat_dir(path, name)+'.elvlc'  
ENDIF ELSE BEGIN                         ; .prot file should be
  zion2filename,ion_z,ion_sp,ename, diel=diel
  ename=ename+'.elvlc'
ENDELSE
read_elvlc,ename,l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,ref
;;----------------------------[O]


sample_wid

ind = where(emiss_sel.label NE '')
nind = n_elements(ind)

str = {emiss: fltarr(n_elements(dens)), $
       wavel: ''}
str2 = replicate(str,nind)

FOR i=1,nind DO BEGIN
  str2[i-1].emiss = emiss_sel[i-1].em
  str2[i-1].wavel = emiss_sel[i-1].label
ENDFOR

IF r_tst EQ 0 THEN BEGIN
  rphot=-1 
  radtemp=-1
ENDIF

;
IF keyword_set(temperature) THEN BEGIN
  temp=dens
  dens=temp[0]
ENDIF ELSE BEGIN
  temp=temp[0]
ENDELSE
;
fname=concat_dir(!xuvtop,'VERSION')
str1=''
IF file_exist(fname) THEN BEGIN
  openr,lun,fname,/get_lun
  readf,lun,str1
  free_lun,lun
  str1='CHIANTI '+str1
ENDIF ELSE BEGIN
  str1='Please update your version of CHIANTI to v.4 or later'
ENDELSE
;
em = {ion: ion_name, $
      lines: str2, $
      dens: float(dens), $
      temp: float(temp), $
      rphot: rphot, $
      radt: radtemp, $
      proton: prot_mess[prot_SWITCH], $
      date: systime(), $
      version: str1 }

tvlct,r,g,b

END

