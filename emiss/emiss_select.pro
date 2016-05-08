
;+
; PROJECT:  CHIANTI
;
;       CHIANTI is an atomic database package for the calculation of
;       astrophysical emission line spectra.  It is a collaborative project
;       involving Ken Dere (Naval Research Laboratory, Washington DC), 
;       Brunella Monsignori-Fossi and Enrico Landi (Arcetri Observatory, 
;       Florence), and Helen Mason and Peter Young (DAMTP, Cambridge Univ.).
;
; NAME
;
;	EMISS_SELECT()
;       
; PURPOSE
;
;	Allows emissivity arrays within EMISS structures to be selected 
;	via a widget.
;
; CATEGORY:
;
;       Atomic data analysis
;
; CALLING SEQUENCE
;
;	em=emiss_select(emiss,index,sel_ind=sel_ind)
;
; EXAMPLES
;
;	emiss=emiss_calc(26,13)       ; calculate emiss for FeXIII
;	em=emiss_select(emiss)
;	em=emiss_select(emiss,indgen(10))
;	em=emiss_select(emiss,wra=[200,300],sel_ind=sel_ind)
;
; INPUTS
;
;	EMISS	A structure that is output by the routine EMISS_CALC
;
; OPTIONAL INPUTS
;
;	INDEX	A list of indices that allow a sub-set of wavelengths to 
;		be displayed in the widget. If not specified, then 
;		complete list of wavelengths is shown.
;
;	WRANGE	Only show wavelengths lying in the specified wavelength 
;		range, e.g., wrange=[200,300]
;
;	SEL_IND	Contains the emiss index/indices of the selected 
;		wavelengths.
;
;       GROUP   If emiss_select is being called from another widget-based 
;               routine, then one needs to set GROUP to the widget ID of 
;               the top-level widget.
;
; OUTPUTS
;
;	An array that contains the emissivities. If two or more lines were 
;	selected, then the arrays are added together.
;
; CALLS
;
;
; HISTORY
;
;	Ver.1, PRY 1-Sep-97.
;	Ver.2, PRY 23-Feb-99,  added WRANGE keyword
;       Ver.3, PRY 13-Sep-00,  changed format of wavelengths in the 
;                  widget list
;       Ver.4, PRY 20-Nov-00,  removed call to xselect_s - the widget 
;                  is now created internally.
;       Ver.5, PRY 27-Dec-00,  changed switch to tst1 for IDL v5.4
;       Ver.6, PRY 6-Aug-02, added transition information to the widget, 
;                  and tidied up some of the code.
;       
;       V.7, 15-Aug-2002 Giulio Del Zanna (GDZ) 
;
;             Replaced the entire  widget with a call to ch_xmenu_sel.
;             This works fine with very long lists of lines, while previously
;             the widget was very slow and sometimes would hang.
;    
;       V.8, 14-Jun-2010, Peter Young
;             Converted some arrays/integers to long.
;
; VERSION     :  8, 14-Jun-2010
;
;-

;;------------------------
;PRO lbase_event, event

;COMMON list, data, dataind

;WIDGET_CONTROL,Event.top, get_uvalue=state

;CASE event.id of

;  state.list: BEGIN
;    widget_control,state.list,get_val=aaa
;    dataind=where(aaa EQ 1)
;  END
 
;  state.extras: BEGIN
;    CASE event.value OF
;      0: BEGIN
;        dataind=-1
;        WIDGET_CONTROL, event.top, /DESTROY ; quit
;      END
;      1: WIDGET_CONTROL, event.top, /DESTROY ; quit
;    ENDCASE
;  END

;ENDCASE

;END




;;------------------------
FUNCTION emiss_select, EMISS, INDEX, SEL_IND=SEL_IND, WRANGE=WRANGE, $
                GROUP=GROUP

COMMON list,data,dataind


n_elt=N_ELEMENTS(emiss.lambda)
 
IF N_ELEMENTS(index) EQ 0 THEN index=lindgen(n_elt)
 
IF N_ELEMENTS(wrange) EQ 2 THEN BEGIN
  wavels=emiss(index).lambda
  index=WHERE( (wavels GE wrange(0)) AND (wavels LE wrange(1)) )
  IF index[0] EQ -1 THEN BEGIN
    IF n_elements(group) EQ 0 THEN BEGIN
      PRINT,' ** There are no lines in the requested wavelength range! **'
    ENDIF ELSE BEGIN
      result=WIDGET_MESSAGE(['There are no lines in the', $
                             'specified wavelength range.'],/INFORMATION)
    ENDELSE
    RETURN,-1
  ENDIF
ENDIF

IF n_elements(index) EQ 0 THEN BEGIN
  emiss2=emiss
ENDIF ELSE BEGIN
  emiss2=emiss[index]
  n_elt=n_elements(index)
ENDELSE
 
;------------------------o
; Wavelengths that are theoretical are shown with a "*" next to them in 
; the widget. The position of the theoretical wavelengths in the emiss 
; structure is denoted by the "flag" tag.
; 
ind_add=WHERE(emiss2.flag EQ -1)
text_add=make_array(n_elt,val='  ')
IF ind_add(0) NE -1 THEN text_add(ind_add)=' *'
;------------------------o


;
; 'options' contains the list of wavelengths to be displayed in the widget
;
wavels=strtrim(STRING(format='(f15.3)',emiss2.lambda),2)
mlen=max(strlen(wavels))
wavels=strpad(wavels,mlen)+' '+string(197b)+text_add
;
lvl_desc=emiss2.lvl1_desc+' - '+emiss2.lvl2_desc
mlen=max(strlen(lvl_desc))
lvl_desc=strpad(lvl_desc,mlen+2,/after)
;
options=wavels+'  '+lvl_desc+'  ['+$
     strtrim(string(emiss2.level1),2)+'-'+$
     strtrim(string(emiss2.level2),2)+']'
mlen=max(strlen(options))
xsiz=mlen*8


;
; The following sets up the widget for choosing the wavelength(s).
;
data=options
dataind=long(-1)

;
; use Giulio's font from ch_ss
;
;font = get_dfont('-adobe-helvetica-bold-r-normal*120*75*')
;font=get_dfont('-misc-fixed-bold-r-normal--13-100-100-100-c-70-iso8859-1')
font=get_dfont('-misc-fixed-bold-r-normal--13-120-75-75-c-70-iso8859-1')
IF font(0) NE '' THEN font = font(0) ELSE font = 'fixed'
widget_control,default_font=font

;-------------------------<X>
; If emiss_select is being called from another widget (e.g., dens_plotter) 
; then I need to set modal to 1 which means that dens_plotter has to wait 
; until a selection has been made with emiss_select.
;
; IF n_elements(group) EQ 0 THEN modal=0 ELSE modal=1


;title='Choose one or more wavelengths and click Continue'
;lbase=widget_base(col=1, MAP=1, $
;        UVALUE='lbase',title=title,modal=modal,group_leader=group)
;list=cw_bgroup(lbase, options, $
;               x_scroll_siz=xsiz, $
;               y_scroll_siz=500,/scroll, $
;               /nonexclusive)
;info=widget_label(lbase,val='* indicates that the wavelength is theoretical', $
;                 /align_left)
;EXTRAS = CW_BGROUP( LBASE, ['Quit','Continue'], /row)

;state={list:list, extras:extras}

;WIDGET_CONTROL, lbase, /REALIZE, set_uvalue=state
;XMANAGER, 'lbase', lbase
;;-------------------------<X>


;
; dataind contains the indices of those lines which have been selected
;
 dataind = ch_xmenu_sel(options, group=group, $
  tit='Choose one or more wavelengths', $
   text=['If more than one line is selected, ',$
         ' the emissivities of the lines will be summed.', $
         'A  * indicates that the wavelength is theoretical' ])


;----------------------------------x
; Take the selected wavelengths and store the total emissivity
; into cal_emiss
;
IF dataind[0] NE -1 THEN BEGIN
  cal_emiss=emiss2[dataind[0]].em
  n=n_elements(dataind)
  IF n GE 2 THEN BEGIN
    FOR i=1,n-1 DO cal_emiss=cal_emiss+emiss2[dataind[i]].em
  ENDIF
  sel_ind=index[dataind]
ENDIF ELSE cal_emiss=-1
;----------------------------------x
 
RETURN,cal_emiss
 
END
