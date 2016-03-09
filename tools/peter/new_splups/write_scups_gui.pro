
;+
; NAME
;
;     WRITE_SCUPS_GUI
;
; PROJECT
;
;     CHIANTI
;
; EXPLANATION
;
;     This routine provides a way of manually "fitting" the upsilon
;     data stored in the .UPS files. The scaled upsilons are
;     plotted, and the user has the options to:
;      1. remove some of the points
;      2. choose a different transition type
;      3. change the scaling parameter
;
; INPUTS
;
;     ION_NAME  The ion name in CHIANTI format (e.g., 'o_6' for O
;               VI). The routine assumes that the upsilons are stored
;               in file the file ION_NAME.ups.
;
; KEYWORDS
;
;     NO_FLAG   The routine assumes that the only transitions it needs
;               to fit are identified in the ION_NAME_flags.txt
;               file. By setting /NO_FLAG the flag file is ignored
;
; INTERNAL ROUTINES
;
;     SCUPS_WRITE_SPLSTR
;
; HISTORY
;
;     Ver.1, 28-May-2013, Peter Young
;-


;-------------------
PRO gui_write_splstr, state
;
; This writes the spline structure to the .scups file.
;

splfile=state.wid_data.splfile
upsstr=state.upsstr
ind=state.wid_data.index

;
; Make a backup of the current SCUPS file (if it exists)
;
chck=file_search(splfile)
IF chck[0] NE '' THEN BEGIN 
  backup_file=splfile+'_backup2'
  file_copy,splfile,backup_file,/overwrite
ENDIF 


read_scups,splfile,splstr,/add_empty

i=splstr.info.ntrans-1

splstr.data[i].lvl1=upsstr.data[ind].lvl1
splstr.data[i].lvl2=upsstr.data[ind].lvl2
splstr.data[i].de=upsstr.data[ind].de
splstr.data[i].gf=upsstr.data[ind].gf
splstr.data[i].lim=upsstr.data[ind].lim
splstr.data[i].c_ups=state.wid_data.c_ups
splstr.data[i].t_type=state.wid_data.t_type

st=state.wid_data.st
sups=state.wid_data.sups
missing=state.wid_data.missing
k=where(st NE missing,nk)
splstr.data[i].stemp[0:nk-1]=st[0:nk-1]
splstr.data[i].spl[0:nk-1]=sups[0:nk-1]
splstr.data[i].nspl=nk

write_splstr,splstr,/overwrite

END




;-----------------------
function gui_info_string, state
;
; This creates the information string displayed in the GUI.
;

ind=state.wid_data.index
ntrans=state.upsstr.info.ntrans
n_remain=ntrans-ind
str=[ 'Ion: '+state.upsstr.info.ion_roman, $
      '', $
      'No. of transitions to fit: '+trim(ntrans), $
      'No. of transitions remaining: '+trim(n_remain), $
      'No. of transitions completed: '+trim(state.wid_data.fit_count), $
      'No. of transitions skipped: '+trim(state.wid_data.skip_count), $
      '', $
      'Current transition: ', $
      'Indices: '+trim(state.upsstr.data[ind].lvl1)+'-'+trim(state.upsstr.data[ind].lvl2) ]
      
return,str
END


;--------------
PRO gui_plot_spl, state
;
; This plots the scaled upsilons.
;

wset,state.wid_data.spl_id

k=where(state.wid_data.st NE state.wid_data.missing)
st=state.wid_data.st[k]
sups=state.wid_data.sups[k]

ind=state.wid_data.index
l1=state.upsstr.data[ind].lvl1
l2=state.upsstr.data[ind].lvl2

title=trim(l1)+'-'+trim(l2)

plot,st,sups,xra=[-0.05,1.05],psym=6,title=title,/xsty, $
     yra=[0,max(sups)*1.10],/ysty
oplot,[0,0],[-1e5,1e5],line=2
oplot,[1,1],[-1e5,1e5],line=2

y2=spl_init(st,sups)
sti=findgen(51)/50.
supsi=spl_interp(st,sups,y2,sti)
oplot,sti,supsi

END



;------------
PRO gui_scale_ups, state
;
; This subroutine loads the st and sups data. Note that the following
; wid_data parameters need to have already been set:
;   .index
;   .t_type
;   .c_ups
;

;
; reset the stored st and sups arrays
;
state.wid_data.st=state.wid_data.missing
state.wid_data.sups=state.wid_data.missing

;
; Get the transition type
;
t_type=state.wid_data.t_type

;
; Get the C-value
;
c_ups=state.wid_data.c_ups


udata=state.upsstr.data[state.wid_data.index]


t_ind0=state.wid_data.t_ind0
t_ind1=state.wid_data.t_ind1
temp=udata.temp[t_ind0:t_ind1]
upsilon=udata.ups[t_ind0:t_ind1]

de_in=udata.de
c_ups_curr=state.wid_data.c_ups
;
kt=1.38062d-16*temp/2.1799d-11   ; in Ryd
;
CASE t_type OF

1: begin
  xt=kt/de_in
  st=1.-alog(c_ups_curr)/alog(xt+c_ups_curr)
  sups=upsilon/alog(xt+exp(1.))
END
;
2: BEGIN
  xt=kt/de_in
  st=xt/(xt+c_ups_curr)
  sups=upsilon
END
;
3: BEGIN
  xt=kt/de_in
  st=xt/(xt+c_ups_curr)
  sups=(xt+1.)*upsilon
END
;
4: BEGIN
  xt=kt/de_in
  st=1.-alog(c_ups_curr)/alog(xt+c_ups_curr)
  sups=upsilon/alog(xt+c_ups_curr)
END
;
5: BEGIN
  xt=kt/de_in
  st=xt/(xt+c_ups_curr)
  sups=upsilon*kte
END
;
6: BEGIN
  xt=kt/de_in
  st=xt/(xt+c_ups_curr)
  sups=alog10(upsilon)
END
;
ELSE:  print,' t_type_ups ne 1,2,3,4=',t_type_ups
;
ENDCASE


;
; extrapolate to zero
;
cc=linfit(st[0:1],sups[0:1])
st=[0.,st]
sups=[cc[0],sups]

;
; extrapolate to one (if necessary)
;
IF udata.lim EQ state.upsstr.info.missing THEN BEGIN
  n=n_elements(sups)
  cc=linfit(st[n-2:n-1],sups[n-2:n-1])
  st=[st,1.]
  sups=[sups,cc[0]+cc[1]*1.0]
ENDIF ELSE BEGIN
  st=[st,1.]
  sups=[sups,udata.lim]
ENDELSE 

n=n_elements(st)

state.wid_data.st[0:n-1]=st
state.wid_data.sups[0:n-1]=sups

widget_control,state.scups_base,set_uval=state

END


;---------------------
PRO scups_load_trans, state
;
; Given wid_data.index, this routine computes t_type and c_ups and
; stores them in wid_data.
;

ind=state.wid_data.index

de=state.upsstr.data[ind].de
gf=state.upsstr.data[ind].gf
t_type=burly_get_ttype(gf,de)
state.wid_data.t_type=t_type

nt=state.upsstr.data[ind].nt
temp=state.upsstr.data[ind].temp[0:nt-1]
c_ups=burly_find_mid_c(temp,de,t_type)
state.wid_data.c_ups=c_ups

state.wid_data.t_ind0=0
state.wid_data.t_ind1=nt-1

widget_control,state.lo_slider,set_val=0
widget_control,state.lo_slider,set_slider_max=nt-1
widget_control,state.hi_slider,set_val=nt-1
widget_control,state.hi_slider,set_slider_max=nt-1
widget_control,state.type_buts,set_val=t_type-1
widget_control,state.c_over_val,set_val=string(format='(f10.4)',c_ups)

;IF state.upsstr.data[ind].lim NE -1 THEN state.wid_data.use_limit=1 ELSE state.wid_data.use_limit=0

widget_control,state.scups_base,set_uval=state

END



;---------------------
PRO scups_base_event, event
;
; This is the event handling routine.
;

WIDGET_CONTROL,Event.top, get_uvalue=state

CASE 1 OF
  event.id EQ state.exit: BEGIN
    CASE event.value OF 
      0: BEGIN                   ; --- RE-FIT TRANSITION
      END
     ;
      1: BEGIN                   ; --- SKIP TRANSITION
        state.wid_data.index=state.wid_data.index+1
        IF state.wid_data.index EQ state.upsstr.info.ntrans THEN BEGIN
          result=dialog_message('All transitions have been processed. Please exit the program.')
        ENDIF ELSE BEGIN
          state.wid_data.skip_count=state.wid_data.skip_count+1
          scups_load_trans,state
          gui_scale_ups, state
          gui_plot_spl, state
          str=gui_info_string(state)
          widget_control,state.info_txt,set_val=str
        ENDELSE 
      END 
     ;
      2: BEGIN                   ; --- NEXT TRANSITION
        state.wid_data.index=state.wid_data.index+1
        IF state.wid_data.index EQ state.upsstr.info.ntrans THEN BEGIN
          result=dialog_message('All transitions have been processed. Please exit the program.')
        ENDIF ELSE BEGIN
          widget_control,/hourglass
          state.wid_data.fit_count=state.wid_data.fit_count+1
          gui_write_splstr, state
          scups_load_trans,state
          gui_scale_ups, state
          gui_plot_spl, state
          str=gui_info_string(state)
          widget_control,state.info_txt,set_val=str
        ENDELSE 
      END
     ;
      3: BEGIN                   ; --- EXIT
       ;
       ; If the user exits without having fitted a transition, then don't
       ; want to save the spline values. Check if l1=0.
       ;
;        IF n_tags(state.splstr) NE 0 THEN  $
;             write_splups,state.data.splupsname,state.data
        WIDGET_CONTROL, event.top, /DESTROY
      END
      ELSE:
    ENDCASE
  END

  event.id EQ state.lo_slider: BEGIN
    widget_control,state.lo_slider,get_val=val
    state.wid_data.t_ind0=val
   ;
   ; Update scaling parameter
   ;
    t_ind0=state.wid_data.t_ind0
    t_ind1=state.wid_data.t_ind1
    ind=state.wid_data.index
    de=state.upsstr.data[ind].de
    nt=state.upsstr.data[ind].nt
    temp=state.upsstr.data[ind].temp[t_ind0:t_ind1]
    t_type=state.wid_data.t_type
    c_ups=burly_find_mid_c(temp,de,t_type)
    state.wid_data.c_ups=c_ups
    widget_control,state.c_over_val,set_val=string(format='(f10.4)',c_ups)
   ;
    widget_control,state.scups_base,set_uval=state
    gui_scale_ups, state
    gui_plot_spl, state
  END 

  event.id EQ state.hi_slider: BEGIN
    widget_control,state.hi_slider,get_val=val
    state.wid_data.t_ind1=val
   ;
   ; Update scaling parameter
   ;
    t_ind0=state.wid_data.t_ind0
    t_ind1=state.wid_data.t_ind1
    ind=state.wid_data.index
    de=state.upsstr.data[ind].de
    nt=state.upsstr.data[ind].nt
    temp=state.upsstr.data[ind].temp[t_ind0:t_ind1]
    t_type=state.wid_data.t_type
    c_ups=burly_find_mid_c(temp,de,t_type)
    state.wid_data.c_ups=c_ups
    widget_control,state.c_over_val,set_val=string(format='(f10.4)',c_ups)
   ;
    widget_control,state.scups_base,set_uval=state
    gui_scale_ups, state
    gui_plot_spl, state
  END 

  event.id EQ state.type_buts: BEGIN
    widget_control,state.type_buts,get_val=val
    t_type=val+1
    state.wid_data.t_type=t_type
   ;
   ; Update scaling parameter
   ;
    t_ind0=state.wid_data.t_ind0
    t_ind1=state.wid_data.t_ind1
    ind=state.wid_data.index
    de=state.upsstr.data[ind].de
    nt=state.upsstr.data[ind].nt
    temp=state.upsstr.data[ind].temp[t_ind0:t_ind1]
    c_ups=burly_find_mid_c(temp,de,t_type)
    state.wid_data.c_ups=c_ups
    widget_control,state.c_over_val,set_val=string(format='(f10.4)',c_ups)
   ;
    widget_control,state.scups_base,set_uval=state
    gui_scale_ups, state
    gui_plot_spl, state
  END 

  event.id EQ state.c_over_val: BEGIN
    widget_control,state.c_over_val,get_val=val
    state.wid_data.c_ups=val
    widget_control,state.scups_base,set_uval=state
    gui_scale_ups, state
    gui_plot_spl, state
  END

ENDCASE 

END 


;----------------------
PRO manu_scups_widget, group=group, upsstr=upsstr

missing_val=-1

nmax=max(upsstr.data.nt)+2

splfile=trim(upsstr.info.ion_name)+'.scups'

wid_data={ spl_id: 0, $
           splfile: splfile, $
           st: fltarr(nmax), $
           sups: fltarr(nmax), $
           t_type: 0, $
           c_ups: 0., $
           index: 0, $
           t_ind0: 0, $
           t_ind1: 0, $
           skip_count: 0, $
           fit_count: 0, $
           missing: missing_val}
          


; Set main base
;
scups_base=widget_base(/col,map=1,title='WRITE_SCUPS_GUI')

chianti_font,font
chianti_font,bigfont,/big
chianti_font,fixfont,/fix

base_top=widget_base(scups_base,/row,map=1)
base_left=widget_base(base_top,/col,map=1,/frame)
base_right=widget_base(base_top,/col,map=1,/frame)

exit=cw_bgroup(base_left, $
               ['RE-FIT','SKIP','WRITE and NEXT','EXIT'], $
               /row,font=bigfont)

sub_base1=widget_base(base_left,/row,map=1)


sub_base5=widget_base(sub_base1,/col,map=1)
info_txt=widget_text(sub_base5,val='',ysiz=10,xsiz=40,font=font) 


;
; Sliders for temperature range
;
n=upsstr.data[0].nt
t_ind0=0
t_ind1=n-1
wid_data.t_ind0=t_ind0
wid_data.t_ind1=t_ind1
sub_base3=widget_base(sub_base1,/col,map=1,/frame)
lohi_lbl=widget_label(sub_base3,val='Temperature range',font=font)
lo_slider=widget_slider(sub_base3,min=0,max=n-1,value=t_ind0)
hi_slider=widget_slider(sub_base3,min=0,max=n-1,value=t_ind1)
;; but_lbl=widget_label(sub_base3,val='Save range?',font=font)
;; but_slider=cw_bgroup(sub_base3,['No','Yes'],font=font, $
;;                      set_val=0,/row,/exclusive)



;
; Buttons for selecting transition type. Note that the number of transition
; types is hardcoded here.
;
gf=upsstr.data[0].gf
de=upsstr.data[0].de
t_type=burly_get_ttype(gf,de)
wid_data.t_type=t_type
;
sub_base4=widget_base(sub_base1,/col,map=1)
type_buts=cw_bgroup(sub_base4,['Type 1 (allowed)', $
                               'Type 2 (forbidden)', $
                               'Type 3 (intercombination)', $
                               'Type 4 (allowed, small gf)', $
                               'Type 5 (dielectronic)', $
                               'Type 6 (log forbidden)'], $
                    set_val=t_type-1,/col,/exclusive,font=font,/frame)

de=upsstr.data[0].de
nt=upsstr.data[0].nt
temp=upsstr.data[0].temp[0:nt-1]
c_ups=burly_find_mid_c(temp,de,t_type)
wid_data.c_ups=c_ups
sub_base5=widget_base(sub_base1,/row,map=1)
c_over_val_lbl=widget_label(sub_base5,val='C-value:',font=font)
c_over_val=widget_text(sub_base5,value=trim(string(format='(f10.4)',c_ups)),/editable,font=font, $
                       xsiz=10,sens=1)


;
; Plot widgets
;
plot_base=widget_base(scups_base,/row,map=1)
;
spl_plot=WIDGET_DRAW(plot_base, $
                       RETAIN=1, $
                       uval='splups_plot', $
                       XSIZE=600, $
                       YSIZE=400)



;
; Set up structure that contains all information used by the program.
;
state={wid_data: wid_data, $
       upsstr: upsstr, $
       exit: exit, $
       spl_plot: spl_plot, $
       info_txt: info_txt, $
       scups_base: scups_base, $
       lo_slider: lo_slider, $
       hi_slider: hi_slider, $
       type_buts: type_buts, $
       c_over_val: c_over_val}


WIDGET_CONTROL, scups_base, /REALIZE, set_uvalue=state

WIDGET_CONTROL, spl_plot, GET_VALUE=spl_id

state.wid_data.spl_id=spl_id

WIDGET_CONTROL, scups_base, /REALIZE, set_uvalue=state

str=gui_info_string(state)
widget_control,state.info_txt,set_val=str

gui_scale_ups,state
gui_plot_spl,state

XMANAGER, 'scups_base', scups_base, group=group


END



PRO write_scups_gui, ionname, no_flag=no_flag

upsfile=ionname+'.ups'
chck=file_search(upsfile)
IF chck[0] EQ '' THEN BEGIN
  print,'%WRITE_SCUPS_GUI: The .ups file ('+upsfile+') was not found. Returning...'
  return
ENDIF 
read_ups,upsfile,upsstr

ntrans=upsstr.info.ntrans
k=lonarr(ntrans)+1    ; initially select all transitions in upsstr.data


; 
; The following checks if .scups already exists. If yes, then it is
; copied to .scups_backup1.  If no, then check if the .scups_auto
; file exists and copy it to .scups (unless /no_flag is set).
;
; Note that there's no need to backup the auto file.
;
splfile=ionname+'.scups'
autofile=splfile+'_auto'
chck=file_search(splfile)
IF chck[0] NE '' THEN BEGIN 
  backup_file=splfile+'_backup1'
  file_copy,splfile,backup_file,/overwrite
  print,''
  print,'** A .scups file already exists, this has been copied to '+backup_file
  print,''
ENDIF ELSE BEGIN
  chck2=file_search(autofile)
  IF chck2[0] NE '' AND NOT keyword_set(no_flag) THEN BEGIN
    file_copy,autofile,splfile
    print,''
    print,'**The '+autofile+' file has been copied to '+splfile
    print,''
  ENDIF 
ENDELSE 


;
; The following code creates a new upsstr (new_ups) that contains only
; those transitions that are going to be fit with write_splupx_gui.
;
; Firstly, if the flag file exists, then this is used to filter
; upsstr.
;
; Secondly, if the .scups file (not .scups_auto) exists, then the
; transitions are further filtered using this file.
;


;
; The logic below is that the scups_auto file will only be read if
; the flag file exists. If the flag file doesn't exist, then it
; is assumed that the user wants to manually fit all the transitions
; (this if scups_auto exists, then it will be over-written).
;
IF NOT keyword_set(no_flag) THEN BEGIN 
  flagfile=ionname+'_flags.txt'
  chck2=file_search(flagfile)
  IF chck2[0] NE '' THEN BEGIN

    k=lonarr(ntrans)   ; de-select all transitions

    openr,lin,flagfile,/get_lun
    WHILE eof(lin) NE 1 DO BEGIN
      readf,format='(3i7)',lin,a,b,c
      i=where(upsstr.data.lvl1 EQ a AND upsstr.data.lvl2 EQ b,ni)
      IF ni NE 0 THEN k[i[0]]=1
    ENDWHILE 
    free_lun,lin

  ENDIF 
ENDIF 


;
; The following checks if a SCUPS file already exists. If yes, then
; a backup copy of this file is made, and any transitions that are
; already in the SCUPS file are flagged as zeros in the K array.
;
chck=file_search(splfile)
IF chck[0] NE '' THEN BEGIN 
  read_scups,splfile,splstr
 ;
 ; Flag as 0 any transitions that are already in the scups file.
 ;
  FOR i=0,ntrans-1 DO BEGIN
    j=where(upsstr.data[i].lvl1 EQ splstr.data.lvl1 AND upsstr.data[i].lvl2 EQ splstr.data.lvl2,nj)
    IF nj GT 0 AND k[i] EQ 1 THEN k[i]=0
  ENDFOR 
ENDIF 


;
;  Check which transitions have K=1, and create a new upsilon
;  structure that contains only these transitions. These are the
;  transitions that will be fit with WRITE_SCUPS_GUI
;
j=where(k EQ 1,nj)
IF nj GT 0 THEN BEGIN 
  new_ups={info: upsstr.info, $
           data: upsstr.data[j] }
  new_ups.info.ntrans=nj
  new_ups.info.time_stamp=systime(1)
ENDIF ELSE BEGIN
  print,'%WRITE_SCUPS_GUI:  There are no transitions left to fit for this ion!'
  return
ENDELSE 


;
; Start up the GUI.
;
manu_scups_widget, upsstr=new_ups


END
