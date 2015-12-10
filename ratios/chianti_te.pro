;+
; Project     : SOHO - CDS     
;                   
; Name        : CHIANTI_TE
;               
; Purpose     : Calculate and plot CHIANTI temperature sensitive line ratios.
;               
; Explanation : CHIANTI_TE (temperature ratios)
;               calculates and plots temp. sensitive line ratios based on 
;               the CHIANTI atomic database of Dere et. al. 
;               
; Use         : IDL> chianti_te 
;    
; Inputs      : None
;               
; Opt. Inputs : None
;               
; Outputs     : None
;               
; Opt. Outputs: none. 
;               
; Keywords    : None
;
; Calls       : CALC_DMM_TR, make_ion_list
;
; Common      : dmm_tr_com
;               
; Restrictions: None
;               
; Side effects: None
;               
; Category    : Spectral
;               
; Prev. Hist. : Started life as 'temperature_ratios' by Ken Dere
;
; Written     : H.E. Mason, 3 Oct 1996 
;
; Modified    : Update for IDL v5.2.  CDP, 20-Apr-99
;               V.3.  Update list of elements.  CDP, 18-Jun-99
;
;               V.4 , Giulio Del Zanna (DAMTP), 10 Oct-2000
;               Rewritten completely, adding possibility to select the density
;               at which the intensities are calculated,  making this routine
;               compatible ith CHIANTI v.3, and with the same characteristics as
;               CHIANTI_TE. Added ratio plots and hardcopies in linear scale, added various
;               checks and  minor things. Added a few tags to the 
;               ratio output structure (density, units, comment). 
;               Removed optional output structure.
;		Version 5, 21-Dec-2000, William Thompson, GSFC
;			Modified for better cross-platform graphics capability
;
;               V. 6, 1-May-02, GDZ, commented out a few refs.
;               V.7,  21-May-2002, GDZ, fixed a bug with the fonts.
;               V.8, 17-Sep-2002 (GDZ) 
;                    added !p.multi = 0 upon exit and added X-label
;
; Version     : Version 8, 17-Sep-2002
;-            

pro chianti_te_event,ev
;
;  event processor for temperature ratio program (chianti_te)
;

;
;  pick up common definitions
;

common dmm_tr_com,   dmm_tr_base, image_area,$
  mess_win, plot_title,      $
  xuvtop, wavmin_id, wavmax_id, tempmin_id, tempmax_id,$
  wslider_min, wslider_max, dslider_min, dslider_max,$
  el_info, ion_info, density_info, inten_info, unit_info,$
  ion_base, els, ions, phot_mod, wind_keep

common dmm_lines_te, list_wvl, list_int, list_descr1, list_descr2, species,$
  temperature, ntem, ratio, description, nlist, savetext, $
  dens, intens

common wavmm, wminw1, wminw2, wminw3, wminw4, $
  wmaxw1, wmaxw2, wmaxw3, wmaxw4, $
  minw1, minw2, minw3, minw4, $
  maxw1, maxw2, maxw3, maxw4

common dmm_refs_te, ioneq_ref, wgfaref, eref, upsref

;
;  get the values, user values and type of widget involved
;
widget_control,ev.id,get_uvalue=uvalue
name = strmid(tag_names(ev,/str),7,1000)
if n_elements(uvalue) eq 0 then uvalue = ''
if n_elements(value)  eq 0 then  value = ''

;print,'name: ',name,' value: ',value,' uvalue: ',uvalue
;help,ev,/str

;
;  act according to type of widget
;
;print,'**',name,'**'
case name of
;
;  Button will be the type of most of the user input from various wwidgets
;
   'BUTTON': begin
      case strmid(uvalue,0,4) of
;
;  Exit from the system (read text widgets first)
;
         'EXIT': begin
            DEVICE, WINDOW_STATE=STATE
            IF N_ELEMENTS(STATE) GT 32 THEN FOR I_WIN=32,N_ELEMENTS(STATE)-1 DO $
              IF STATE(I_WIN) THEN WDELETE, I_WIN
            widget_control,/destroy,ev.top
!p.multi = 0
         end
;               
;  HELP will put the help file to the screen (not implemented yet..)
;
         'HELP': begin
            widg_help,'dmm_te.hlp'
         end

;
;  Change UNITS
;
         'UNIT': begin
            if phot_mod then begin
               widget_control, unit_info, set_val='Units for the ratio plot: ERGS'
               phot_mod = 0
            endif else begin
               widget_control, unit_info, set_val='Units for the ratio plot: PHOTONS'
               phot_mod = 1
            endelse
         end

;
;  Set the calculation going
;
         'CALC': begin
            widget_control,wminw1, get_v=minw1 
            minw1 = float(minw1(0))
            widget_control,wminw2, get_v=minw2 
            minw2 = float(minw2(0))
            widget_control,wminw3, get_v=minw3 
            minw3 = float(minw3(0))
            widget_control,wminw4, get_v=minw4 
            minw4 = float(minw4(0))
            widget_control,wmaxw1, get_v=maxw1 
            maxw1 = float(maxw1(0))
            widget_control,wmaxw2, get_v=maxw2 
            maxw2 = float(maxw2(0))
            widget_control,wmaxw3, get_v=maxw3 
            maxw3 = float(maxw3(0))
            widget_control,wmaxw4, get_v=maxw4 
            maxw4 = float(maxw4(0))

            widget_control, el_info,get_val=el
            el = strtrim(el(0),2)
            widget_control, ion_info,get_val=ion
            ion = strtrim(ion(0),2)

            widget_control, density_info,get_val=dens
            dens= float(strtrim(dens(0),2))
            widget_control, inten_info,get_val=intens
            intens = float(strtrim(intens(0),2))

            widget_control, image_area, get_v=wind

            wset, wind
            !p.multi=[0,2,1]
;            !p.multi=0
            widget_control, dmm_tr_base,/hour

            calc_dmm_tr,el,ion,wslider_min,wslider_max,$
              dslider_min,dslider_max,mess_win, $
              density=dens

            !p.multi=[0,2,1]
;            !p.multi=0
            if nlist gt 0 then plot_dmm_tr_fig
            DEVICE, WINDOW_STATE=STATE
            wind_keep = where(state eq 1)
            wind_keep = wind_keep(0)
         end


;
;  Delete plot windows to tidy up
;
         'DELP': begin
            DEVICE, WINDOW_STATE=STATE
            if wind_keep ge 0 then begin
               state(wind_keep) = 0
               IF N_ELEMENTS(STATE) GT 32 THEN FOR I_WIN=32,N_ELEMENTS(STATE)-1 DO $
                 IF STATE(I_WIN) THEN WDELETE, I_WIN
            endif
         end

;
;  plot ratio after user supplies lines involved
;
         'PRAT': begin

            CASE  strmid(uvalue,5,3) OF 
               'LIN': plot_log = 0
               'LOG': plot_log = 1
            END 

            xinput,text,['Enter lines involved in ratio',$
                         'eg. enter  2,3 4 / 1+5',$
                         'to calculate ratio from sum of lines',$
                         '2, 3 and 4 all divided by lines 1+5'],$
              group=dmm_tr_base,xoff=400,yoff=50
            stext = text
            if strpos(text,'/') gt 0 then begin
               text = repchar(text,',',' ')
               text = repchar(text,'+',' ')
               text = repchar(text,'.',' ')
               text = strtrim(text,2)

               text = str2arr(text,del='/')
               text(0) = strtrim(text(0),2)
               text(1) = strtrim(text(1),2)

               numer = fix(str2arr(text(0),del=' '))
               denom = fix(str2arr(text(1),del=' '))
               numval = 0.0
               denval = 0.0
               nlines = indgen(n_elements(list_wvl))
               for kk=0,n_elements(numer)-1 do begin
                  hit = where(nlines eq numer(kk),count)
                  if count eq 0 then begin
                     bell
                     widget_control, mess_win, /append,$
                       set_val='Invalid line number in numerator --> ' + stext
                     return
                  ENDIF

                 if phot_mod then begin
                     numval = numval + list_int(numer(kk),*)
                  endif else begin
                     numval = numval + list_int(numer(kk),*)/$
                       list_wvl(numer(kk))
                  endelse
               endfor
               for kk=0,n_elements(denom)-1 do begin
                  hit = where(nlines eq denom(kk),count)
                  if count eq 0 then begin
                     bell
                     widget_control, mess_win, /append,$
                       set_val='Invalid line number in denominator --> ' + stext
                     return
                  endif
                  if phot_mod then begin
                     denval = denval + list_int(denom(kk),*)
                  endif else begin
                     denval = denval + list_int(denom(kk),*)/$
                       list_wvl(denom(kk))
                  endelse
               endfor
               ratio = numval/denval
               ratio=reform(ratio)

               if n_elements(numer) gt 1 then begin
                  numdesc = fmt_vect(list_wvl(numer),$
                                     form='(f10.4)',delim='+')
               endif else numdesc = string(list_wvl(numer(0)),'(f10.4)')
               numdesc = strcompress(numdesc,/rem)

               if n_elements(denom) gt 1 then begin
                  dendesc = fmt_vect(list_wvl(denom),$
                                     form='(f10.4)',delim='+')
               endif else dendesc = string(list_wvl(denom(0)),'(f10.4)')
               dendesc = strcompress(dendesc,/rem)
               description=species+'   '+numdesc+'/'+dendesc
               window,/free
               wset,!d.window
               !p.multi=0

               if phot_mod then ytit = 'Intensity (photons) Ratio' else $
                 ytit = 'Intensity (ergs) Ratio'

               IF plot_log THEN $
                 plot_oo,temperature,ratio,xr=[temperature(0),temperature(ntem-1)],$
                 title=description,xtit='Log!d10!n (Electron temperature [K])',$
                  ytitle=ytit ,chars=1.3 ELSE $ 
                 plot,temperature,ratio,xr=[temperature(0),temperature(ntem-1)],$
                 title=description,xtit='Electron temperature [K]', $
                 ytitle=ytit ,chars=1.3

            endif else begin
               bell
               widget_control, mess_win, /append,$
                 set_val='Invalid ratio specification -->' + stext
            endelse
         end

;
;  SAVE the calculated ratio
;
         'DUMP': BEGIN

;check that we have defined the ratio:
            IF n_elements(description) EQ 0 THEN BEGIN
               bell 
               widget_control, mess_win, /append,$
                 set_val='Still no ratio to SAVE!'
            ENDIF ELSE BEGIN 

               xinput,text,'Enter file name to store ratio data (a .CH_TE will be appended)',$
                 group=dmm_tr_base,/modal,xoff=400,yoff=50

               IF  text ne '' THEN  BEGIN

                  if phot_mod then ytit = 'Intensity (photons) Ratio' else $
                    ytit = 'Intensity (ergs) Ratio'

                  te_ratio = {density:dens, $
                              units:ytit, $
                              comment:'Calculated with CHIANTI at contant density', $
                              desc:description,temp:temperature,$
                              ratio:ratio}
                  save,te_ratio,file=text+'.CH_TE',/xdr
                  bell
                  widget_control, mess_win, /append,$
                    set_val='Ratio saved in file '+text+'.CH_TE'
               ENDIF
            endelse
         END 

;
;  HARDCOPY requested
;
         'HARD': BEGIN

            case strmid(uvalue,5,3) of

               'RAT': begin

;check that we have a ratio to plot:

                  IF n_elements(description) EQ 0 THEN BEGIN
                     bell 
                     widget_control, mess_win, /append,$
                       set_val='Still no ratio defined  to plot!'
                  ENDIF ELSE BEGIN 

                     IF strmid(uvalue,8,1) EQ '1' THEN plot_log = 0 ELSE $
                       IF strmid(uvalue,8,1) EQ '2' THEN plot_log = 1 
                     ps
                     !p.multi=0

                     if phot_mod then ytit = 'Intensity (photons) Ratio' else $
                       ytit = 'Intensity (ergs) Ratio'

                     IF plot_log THEN $
                       plot_oo,temperature,ratio,$
                       xr=[temperature(0),temperature(ntem-1)],$
                       title='CHIANTI: '+description,$
xtit='Log!d10!n (Electron temperature [K])',$
                       ytitle=ytit,chars=1.3 ELSE $
                       plot,temperature,ratio,$
                       xr=[temperature(0),temperature(ntem-1)],$
                       title='CHIANTI: '+description,$
xtit='Electron temperature [K]', $
                       ytitle=ytit,chars=1.3

                     get_utc,utc
                     utc = anytim2cal(utc)
                     xyouts, 0.05,0.0,'Printed: '+utc,/norm,$
                       chars=0.7  
                     psplot
                  ENDELSE 
               END  

               'INT': BEGIN
                  
                  if nlist gt 0 then begin
                     ps
                     plot_dmm_tr_fig
                     prin = getenv('PRINTER')
                     IF prin EQ '' THEN $
                       xsel_printer,prin,group=dmm_tr_base,$
                       instruct='Select printer for hardcopy'
                     psplot,queu=prin
                     bell
                     widget_control, mess_win, /append,$
                       set_val='Plot sent to printer '+prin
                  endif else BEGIN
                     bell
                     widget_control, mess_win, /append,$
                       set_val='Still nothing to plot!'
                  endelse
               END 
               'LIN': begin
                  get_utc,utc
                  utc = anytim2cal(utc)
                  text = 'CHIANTI '+species+'    Printed: '+$
                    utc
                  savetext = [text,' ',$
                              'Calculation at constant Ne='+$
                              string(dens,format='(e10.3)')+' (cm-3)', $
                              ' ',$
                              '#    Wavelength  Intensity    Transition ', $
                              '        (A)   ', $
                              savetext]
;                  savetext = [savetext,' ','Ion. Eq. refs:',$
;                              ioneq_ref]
;                  savetext = [savetext,' ','Wgfa refs:',$
;                              wgfaref]
;                  savetext = [savetext,' ','Elvl refs:',$
;                              eref]
;                  savetext = [savetext,' ','Upsilon refs:',$
;                              upsref]
                  print_str,savetext,/hard,/quiet
               end
               else:
            endcase
         end
         else:  print,'Button not recognised'

      endcase
   end

;
;  wavelength sliders
;
   'SLIDER'   : begin
      case strmid(uvalue,0,4) of
         'WAVE': begin
            widget_control,ev.id,get_uvalue=uvalue,get_value=value
            if strmid(uvalue,4,3) eq 'MIN' then begin
               if value(0) ge wslider_max then value(0) = value(0)-2 
               wslider_min = value(0)
               if value(0) gt wslider_max then begin
                  widget_control, wavmax_id,set_val=value(0)+2
                  wslider_max = value(0)+2
               endif
            endif else begin
               if value(0) le wslider_min then value(0) = value(0)+2 
               wslider_max = value(0)
               if value(0) lt wslider_min then begin
                  widget_control, wavmin_id,set_val=value(0)-2
                  wslider_min = value(0)-2
               endif
            endelse  
            widget_control,wavmin_id,set_val=string(wslider_min)
            widget_control,wavmax_id,set_val=string(wslider_max)
         end
         'TEMP': begin
            widget_control,ev.id,get_uvalue=uvalue,get_value=value
            if strmid(uvalue,4,3) eq 'MIN' then begin
               if value(0) ge dslider_max then value(0) = value(0)-1 
               dslider_min = value(0)
               if value(0) ge dslider_max then begin
                  widget_control, tempmax_id,set_val=value(0)+1
                  dslider_max = value(0)+1
               endif
            endif else begin
               if value(0) le dslider_min then value(0) = value(0)+1 
               dslider_max = value(0)
               if value(0) le dslider_min then begin
                  widget_control, tempmin_id,set_val=value(0)-1
                  dslider_min = value(0)-1
               endif
            endelse  
            widget_control,tempmin_id,set_val=string(dslider_min)
            widget_control,tempmax_id,set_val=string(dslider_max)
         end
         else:
      endcase
   end
;
;  must have been from cw_pdmenu
;
   '': begin
      case uvalue of
         'ELEMENT': begin
            widget_control, el_info ,set_val=ev.value
            nit = where(strupcase(els) eq strupcase(ev.value))
            for i=0,n_elements(els) do begin
               widget_control, ion_base(i), map=0
            endfor
            widget_control, ion_base(nit(0)+1), map=1
            widget_control, ion_info ,set_val=''    
         end
         'ION'    : widget_control, ion_info ,set_val=ev.value     
         else:
      endcase
   end

   else: print,'name not recognised'

endcase
end

;
; --------------------------------
;  The main routine
; --------------------------------
;

pro chianti_te,  modal=modal, group_leader=group_leader,$
               font=font


;
;  pick up common definitions
;

common dmm_tr_com,   dmm_tr_base,  image_area,$
                     mess_win, plot_title,      $
                     xuvtop, wavmin_id, wavmax_id, tempmin_id, tempmax_id,$
                     wslider_min, wslider_max, dslider_min, dslider_max,$
                     el_info, ion_info, density_info, inten_info, unit_info,$
                     ion_base, els, ions, phot_mod, wind_keep
  
   
common dmm_lines_te, list_wvl, list_int, list_descr1, list_descr2, species,$
                  temperature, ntem, ratio, description, nlist, savetext, $
                  dens, intens

common wavmm, wminw1, wminw2, wminw3, wminw4, $
              wmaxw1, wmaxw2, wmaxw3, wmaxw4, $
              minw1, minw2, minw3, minw4, $
              maxw1, maxw2, maxw3, maxw4


;
;  pick up top level directory where the atomic data and instrument files exist
;
defsysv,'!xuvtop', EXISTS = EXISTS 
IF NOT EXISTS THEN $
message, 'system variable !xuvtop must be set '
xuvtop = !xuvtop


nlist = 0
wind_keep = -1

;
;  initial checks
;
if (!d.flags and 65536) eq 0 then message,'Widgets are NOT available.'

;
;  check not already in use
;
if xregistered('chianti_te') then begin
   print,'CDS chianti_te program is already registered.'
   return
endif

;clear all ratio arrays
;
delvarx, list_wvl, list_int, list_descr1, list_descr2, species,$
                  temperature, nden, ratio, description

;
;  set only possible graphics device and load LUT
;
set_x
loadct,0

;
;  the base widget is all-screen and column-orientated
;
if not keyword_set(font) then begin
  font = get_dfont('-adobe-helvetica-bold-r-normal*120*')
  IF font(0) NE '' THEN font = font(0) ELSE font = 'fixed'
endif
device, get_screen_size = sz
sz = 0.9*sz

if (sz(0) ge 1280) and (sz(1) ge 1024) then sz(*) = 0
sz = sz < [1280,1024]
dmm_tr_base = widget_base(title='CHIANTI Temperature Dependent Ratios',$
                   /column,/frame,x_scroll=sz(0),y_scroll=sz(1))

;
;  the TOP widget is where all the user/numerical/text input takes place
;
top = widget_base(dmm_tr_base,/row,space=10)

;
;  the left hand side of the top contains...
;
top_left = widget_base(top,/column,space=5)


wltit  = widget_label(top_left,$
             value='Select up to 4 wavelength ranges.',font=font)
ss1    = widget_base(top_left,/row)
wlss   = widget_label(ss1,value='Minimum Wavelength.',font=font)
sss1   = widget_base(ss1,/row)
wminw1 = widget_text(sss1,/edit,xsize=5,ysize=1,uvalue='MINW1',$
                     value='1', font=font)
wminw2  = widget_text(sss1,/edit,xsize=5,ysize=1,uvalue='MINW2',$
                     value='0', font=font)
wminw3  = widget_text(sss1,/edit,xsize=5,ysize=1,uvalue='MINW3',$
                     value='0', font=font)
wminw4  = widget_text(sss1,/edit,xsize=5,ysize=1,uvalue='MINW4',$
                     value='0', font=font)

ss2    = widget_base(top_left,/row)
wlss  = widget_label(ss2,value='Maximum Wavelength.',font=font)
sss2   = widget_base(ss2,/row)
wmaxw1  = widget_text(sss2,/edit,xsize=5,ysize=1,uvalue='MAXW1',$
                     value='1700', font=font)
wmaxw2  = widget_text(sss2,/edit,xsize=5,ysize=1,uvalue='MAXW2',$
                     value='0', font=font)
wmaxw3  = widget_text(sss2,/edit,xsize=5,ysize=1,uvalue='MAXW3',$
                     value='0', font=font)
wmaxw4  = widget_text(sss2,/edit,xsize=5,ysize=1,uvalue='MAXW4',$
                     value='0', font=font)

;
;  temperature range sliders
;
tempmin_id = widget_slider(top_left,title='Log minimum temperature',$
                           uvalue='TEMPMIN',xsize=500,min=4,max=7)
dslider_min = 4

tempmax_id = widget_slider(top_left,title='Log maximum temperature',$
                           uvalue='TEMPMAX',xsize=500,min=5,max=8)
widget_control, tempmax_id,set_val='8'
dslider_max = 8


;
; the top middle widget contains...
;
top_middle= widget_base(top,/column,space=10)

;
;  selection of element
;
cs3 = widget_base(top_middle,/row)
junk = {CW_PDMENU_S,FLAGS:0,NAME:''}


;element=['H','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na',$
;         'Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti',$
;         'V ','Cr','Mn','Fe','Co','Ni','Cu','Zn']

;
;  large font
;             

bfont = get_dfont('-adobe-helvetica-bold-r-normal*180*')
;  bfont = '-adobe-helvetica-bold-r-normal--0-120-75-75-p-0-iso8859-1'
IF bfont(0) NE '' THEN bfont = bfont(0) ELSE bfont = 'fixed'

make_ion_list,els,ions

desc = [{CW_PDMENU_S,FLAGS:1,NAME:'Select element'}]
for i=0,n_elements(els)-1 do begin
  desc = [desc,{CW_PDMENU_S,FLAGS:0,NAME:els(i)}]
endfor
desc(n_elements(desc)-1).flags=2
menu = cw_pdmenu(cs3,desc,/return_name,uvalue='ELEMENT',font=font)
cs4 = widget_base(cs3,/column)
el_info = cw_field(cs4,title='Current element:',value='    ',$
            /row,xsize=4,fieldfont=bfont)

;
;  selection of ion
;

junk = {CW_PDMENU_S,FLAGS:0,NAME:''}
ionstage=['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII',$
          'XIII','XIV','XV','XVI','XVII','XVIII','XIX','XX','XXI','XXII',$
          'XXIII','XXIV','XXV','XXVI','XXVII','XXVIII']
  

cs51 = widget_base(top_middle,/row)
cs5 = widget_base(cs51)

nel = n_elements(els)
ion_base = lonarr(nel+1)

ion_base(0) = widget_base(cs5,map=0)

for k=1,nel do begin
   desc = [{CW_PDMENU_S,FLAGS:1,NAME:'Select ionization'}]
   ionstage = ions(k-1,*)
   ionstage = ionstage(where(ionstage ne ''))
   for i=0,n_elements(ionstage)-1 do begin
     desc = [desc,{CW_PDMENU_S,FLAGS:0,NAME:ionstage(i)}]
   endfor
   desc(n_elements(desc)-1).flags=2
   ion_base(k) = cw_pdmenu(cs5,desc,/return_name,uvalue='ION',font=font)
   widget_control, ion_base(k), map=0
endfor
cs6 = widget_base(cs51,/column)
ion_info = cw_field(cs6,title='Current ion:',value='    ',$
                    /row,xsize=5,fieldfont=bfont)


;
unit_info = widget_button(top_middle,value='Units for the ratio plot: ERGS',uvalue='UNIT',$
                    frame=2)

;define the default:
phot_mod = 0
widget_control,unit_info,set_val='Units for the ratio plot: ERGS'


;
;  density selection
;
eb2 = widget_base(top_middle,/column)
density_info = cw_field(eb2,title='Density (cm-3):',value='1.e8',$
             /row,xsize=8,/return,uvalue='USERN')

;default:
dens = 1.e8

;
;  minimun line intensity selection
;
eb3 = widget_base(top_middle,/column)
inten_info = cw_field(eb3,title='Minimum intensity ratio (for the calculation):',value='0.0001',$
             /row,xsize=10,/return,uvalue='USERR')


widget_control, ion_base(1), map=1

;
;  MIDDLE "command" line has the options...
;
middle = widget_base(dmm_tr_base,/row,space=10,/frame)

;
;  ...Quit without any further ado
;
m1 = widget_base(middle,/column)
dummy = widget_button(m1,value='Quit',uvalue='EXIT',/no_rel,font=font)


;
;  ...calculate new data with current instrument settings
;
m2 = widget_base(middle,/column)
dummy = widget_button(m2,value='CALCULATE line intensities',uvalue='CALC',$
                      /no_rel,font=font)


;
;  replot graph with revised wavelength range
;
m9 = widget_base(middle,/column)
XPDmenu, ['"PLOT ratio" {','"in linear scale"   PRAT_LIN',$
                         '" in log scale"    PRAT_LOG','}'], m9, font=font

;
;  request HARDCOPY
;
m4 = widget_base(middle,/column)
XPDmenu, ['"HARDCOPY" {','"Of ratio plot (linear scale)"   HARD_RAT1',$
                        '"Of ratio plot (log scale)"   HARD_RAT2',$
                         '"Line details (+ refs)"    HARD_LIN',$
                         '"Of intensity plots"  HARD_INT','}'], m4, font=font
 
;
;  dump details
;
m6 = widget_base(middle,/column)
XPDmenu, ['"SAVE" {','"Ratio values"     DUMP_TR','}'],m6, font=font


;dummy = widget_button(m9,value='PLOT ratio',uvalue='PRAT',$
;                      /no_rel,font=font)

;
;  delete plot windows to tidy up
;
m10 = widget_base(middle,/column)
dummy = widget_button(m10,value='Delete plot windows',uvalue='DELP',$
                      /no_rel,font=font)



;
;  BOTTOM display area
;   
bottom = widget_base(dmm_tr_base,/frame,/column)

;
;  ...plotting area
;
image_area = widget_draw(bottom,xsize=800,ysize=225,/button_events)
have_vert = 0
widget_control,image_area,sensitive=0


;
;  messages
;
mess_win = widget_text(dmm_tr_base,ysize=8,xsize=120,/scroll)



;
;  make the whole thing happen
;
widget_control,dmm_tr_base,/realize

;
;  find out window id for plotting
;
widget_control,get_value=window,image_area
wset,window 

;
; make the boss do some work
;
xmanager,'chianti_te',dmm_tr_base, modal=modal, group_leader=group_leader


end


