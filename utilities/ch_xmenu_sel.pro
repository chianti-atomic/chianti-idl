;
; NAME: 
; 	CH_XMENU_SEL_EV
;
; PURPOSE:
;       The event handler for CH_XMENU_SEL
;
PRO  ch_xmenu_sel_ev, event
;
;
;
common ch_xmenu_sel_blk, all

widget_control, event.top, get_uvalue=all
n = n_elements(all.status)

qrefresh = 0
case (strtrim(event_name(event),2)) of
   "BUTTON" : begin
      case strupcase(get_wvalue(event.id)) of
         'OK - EXIT': begin
            widget_control, event.top,/destroy
         endcase
         'SELECT ALL': begin
            if (all.qone) then begin
               tbeep, 3
               print, 'Sorry... You are only allowed to select a single item'
            end else begin
               all.status = 1   ;set them all as active
               all.ordered_index = indgen(n)
            end
            qrefresh = 1
         endcase
         'DE-SELECT ALL': begin
            all.status = 0      ;reset them all
            all.ordered_index[*] = -1
            qrefresh = 1
         endcase
         'SELECT ALL BETWEEN LAST 2 CHOICES': begin
            if (min(all.last_two) eq -1) then begin
               tbeep, 3
               print, 'You must select two items before pressing ALL BETWEEN LAST 2' 
            end else BEGIN 
               i1 = all.last_two(0)<all.last_two(1) 
               i2 = all.last_two(1)>all.last_two(0)
               all.status(i1:i2) = 1
               nsel = i2-i1+1
               all.ordered_index[*] = -1
               all.ordered_index[0:nsel-1] = indgen(nsel)+i1
               qrefresh = 1
            end
         ENDCASE 
         ELSE : print, 'Cannot recognize button ', strupcase(get_wvalue(event.id))
      endcase
   end
   "LIST": begin
      i1 = event.index
      all.status(i1) = abs(all.status(i1)-1) ;toggle the value

      IF all.status(i1) EQ 1 THEN BEGIN 

;find how many we have selected already and add 

         istart = max(where(all.ordered_index NE -1))+1
         all.ordered_index[istart[0]] = i1

      ENDIF ELSE BEGIN 
;remove if present 

         irem = where(all.ordered_index EQ i1)

         IF irem[0] NE -1 THEN BEGIN 

            inum = max(where(all.ordered_index NE -1))

            new_list =  all.ordered_index[0:inum]
            remove, irem, new_list

;now that  we have removed the value(s), stack them back in

            all.ordered_index[*] =-1
            all.ordered_index[0:n_elements(new_list)-1] = new_list

         ENDIF  

      END   

      all.last_two(0) = all.last_two(1)
      all.last_two(1) = i1
      qrefresh = 1
   end
   ;;for i = frst, secnd do widget_control, butt(i), /set_button 
   ELSE: print, 'Cannot  recognize ', (strtrim(event_name(event),2))
ENDCASE 

if (qrefresh) then begin
   ch_xmenu_sel_lab, all, lab
   ii = widget_info(all.id_list, /list_top)
   widget_control, all.id_list, set_value=lab, set_list_top=ii
end
;
widget_control, event.top, set_uvalue=all, bad_id=destroyed ;update structure holding all of the info
;
if (all.qone) and (strtrim(event_name(event),2) eq 'LIST') then begin
   widget_control, event.top,/destroy ;all done - got the single one that was needed
end
;
END 
;=======================================================================

PRO  ch_xmenu_sel_lab, all, lab
;
;
;
n = n_elements(all.status)
ss = where(all.status)
lead = '   '+strarr(n)  
if (ss(0) ne -1) then lead(ss) = ' + '
;
lab = lead + all.list
end
;=======================================================================
;+
;
; NAME: 
;	CH_XMENU_SEL
;
; PURPOSE:
;     	Allow user to select a set of items from a list.  Widget equivalent
;	of WMENU_SEL
;
; CALLING SEQUENCE:
;       ss = ch_xmenu_sel(array)
;       ss = ch_xmenu_sel(array, /one)
;       ss = ch_xmenu_sel(array, /fixed) - use fixed font (keep column alignement)
;
; INPUTS:
;       ARRAY       A string or string array of values to be displayed for
;		    selection
; OPTIONAL KEYWORD INPUT:
;	ONE	    If set then only one button may be turned on at a time.
;	TIT	    The title of the widget
;	GROUP	    The parent widget id (so that if the parent widget exits,
;		    this widget is destroyed too)
;       FIXED_FONT  If set, use fixed font (keep columns aligned)
;       SIZE_FONT   Size of (fixed) font to use - default=15 (implies /FIXED)
;       NLINES      How many lines to display.  Default is 20
;
; OUTPUTS:
;	The result returns the select indices of the array ARRAY.
;
; RESTRICTIONS:
;	Must have widgets available.
;
; HISTORY:
;	Written 30-Jan-95 by M.Morrison using Elaine Einfalt
;	YO_TAPE_WIDG as a starting point
;	10-Jul-96 (MDM) - Ajustment to make the output scaler if it
;			  is a single item
;       11-nov-96 (SLF) - add FIXED_FONT and SIZE_FONT keywords
;	15-Apr-97 (MDM) - Re-added the 9-Jan-97 (MDM) modification
;			  to merge with 11-Nov-96 version
;			      9-Jan-97 (MDM) - Added NLINES option
;	22-Jul-97 (MDM) - Added call to WMENU_SEL if the device is
;			  not X (so that terminal prompting is
;			  enabled.
;
;       9-may-2001 Giulio Del Zanna (GDZ) 
;        added keywords  text and include, renamed ch_xmenu_sel
;
;       V.7, 15-Aug-2002 GDZ 
;         Modified a few cosmetics, and add the option to have a TEXT string
;         array to add as a comment in the top widget.
;    
;       V.8  30 Jan 2002, GDZ
;          Now it returns an ordered list, according to how the lines where
;          selected. It works only the first time the routine is called.
;
;       V.9, 30-Jun-2003, Peter Young
;          Added call to CHIANTI_FONT to correct problem with fixed font in 
;          Windows.
;       V.10 3-Oct-2003, GDZ 
;          inserted the ,/MODAL keyword within the main widget, for 
;          compatibility with IDL v.5.4

; VERSION     :  10 3-Oct-2003
;-
;
function ch_xmenu_sel, array, one=one, tit=tit, group=group, $
           fixed_font=fixed_font, size_font=size_font, nlines=nlines,$
           include = include,  text=text

;
; if (!d.name ne 'X') then begin  ;MDM added 22-Jul-97
;    out = wmenu_sel(array, one=one)
;    return, out
; end
;
;
common ch_xmenu_sel_blk, all	;to pass results back and forth from event driver

if (n_elements(nlines) eq 0) then nlines = 20

n_el   = n_elements(array)
if n_el eq 0 then return, -1


all = {status: bytarr(n_el), $
       list: string(array), $
       last_two: intarr(2)-1, $
       qone: 0B, $
       id_list: 0L, $
       base: 0L, $
       ordered_index: intarr(n_elements(array))-1}


IF n_elements(include) NE 0 THEN BEGIN 
   IF  n_elements(include) NE n_el THEN return, -1

   all.status = include

   index = where(include, nn)
   IF nn GT 0 THEN all.ordered_index[0:nn-1] = index 

ENDIF  


if (n_elements(tit) eq 0) then tit = ''
if (n_elements(group) eq 0) then group = 0 ;needed to avoid crash in widget_base call

base = widget_base(title='CH_XMENU_SEL'+tit, /column, xpad=20, ypad=20, group=group)

xmenu, ['OK - EXIT', 'Select all', 'De-select all', 'Select all between last 2 choices'], $
  base, /row

instruc = widget_base(base, /row, space=10)
labs = widget_base(instruc, /column)
lab = widget_label(labs, value='Selection is made by pointing and clicking')

if keyword_set(one) then begin
   lab = widget_label(labs, value='You can only select a single item')
end else begin
   lab = widget_label(labs, value='A "+" preceding the item indicates it is selected, while')
   lab = widget_label(labs, value='a second click de-selects. Select one or as many as desired')
end
lab = widget_label(labs, value='Click on OK  to  accept selection and exit')

IF n_elements(text) NE 0 THEN $
  FOR i=0, n_elements(text)-1 DO lab = widget_label(labs, value=text[i])

ch_xmenu_sel_lab, all, lab
ba_base = widget_base(base, /column)

fixed=keyword_set(fixed_font) or keyword_set(size_font)
fixed=1

chianti_font,fontfix,/fixed

if fixed then begin
   if n_elements(size_font) eq 0 then size_font=10
   listit = widget_list(ba_base, value=lab, ysize=nlines<n_el, /frame, $
                        font=fontfix)
;                        font=get_xfont(/only_one,/fixed,closest=size_font))
endif else listit = widget_list(ba_base, value=lab, ysize=nlines<n_el, /frame)

all.id_list = listit
all.qone    = keyword_set(one)
all.base    = base

widget_control,set_uvalue=all, base
widget_control, base, /realize
xmanager, 'ch_xmenu_sel', base, $
event_handler='ch_xmenu_sel_ev', /modal

;=(group ne 0), just_reg=(group ne 0)


include = all.status

index = where(all.ordered_index NE -1, nn)
index2 = where(all.status EQ 1b, nn2)
IF nn NE nn2 THEN message, 'Error in CH_XMENU_SEL ?! '


;out = where(all.status)
IF nn EQ 0 THEN out = -1 ELSE out = all.ordered_index[index]

if (n_elements(out) eq 1) then out = out(0) ;make it a scalar
return, out
end
