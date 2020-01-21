
;+
; NAME:
;     WRITE_SCUPS_GUI
;
; CATEGORY:
;     CHIANTI; data assessment.
;
; PURPOSE:
;     This routine provides a way of manually "fitting" the upsilon
;     data stored in the .UPS files. The scaled upsilons are
;     plotted, and the user has the options to:
;      1. remove some of the points
;      2. choose a different transition type
;      3. change the scaling parameter
;     It is recommended that the user run the automatic fitting
;     routine WRITE_SCUPS first, prior to calling WRITE_SCUPS_GUI.
;
; INPUTS:
;     Ion_Name: The ion name in CHIANTI format (e.g., 'o_6' for O
;               VI). The routine assumes that the upsilons are stored
;               in file the file ION_NAME.ups.
;
; OPTIONAL INPUTS:
;     Elvlcfile: The name of a CHIANTI elvlc file. If not set, then
;                the routine searches for [ion_name].elvlc.
;
; KEYWORDS PARAMETERS:
;     NO_FLAG:  The routine assumes that the only transitions it needs
;               to fit are identified in the ION_NAME_flags.txt
;               file. By setting /NO_FLAG the flag file is ignored
;
; OUTPUTS:
;     Creates the file ION_NAME.scups.
;
; INTERNAL ROUTINES:
;     GUI_SET_TEMP_LABELS, GUI_WRITE_SPLSTR, GUI_INFO_STRING,
;     GUI_PLOT_SPL, GUI_SCALE_UPS, SCUPS_LOAD_TRANS, SCUPS_BASE_EVENT,
;     MANU_SCUPS_WIDGET.
;
; EXAMPLES:
;     IDL> write_scups_gui,'mg_5'
;     IDL> write_scups_gui,'mg_5', elvlcfile='mg_5.elvlc_new'
;
; MODIFICATION HISTORY:
;     Ver.1, 28-May-2013, Peter Young
;     Ver.2, 27-Jan-2017, Peter Young
;         Fixed bug when SCUPS file does not already exist.
;     Ver.3, 3-Feb-2017, Peter Young
;         Many changes in preparation for first release.
;     Ver.4, 4-Aug-2017, Peter Young
;         Now reads .elvlc file in order to add transition information
;         to GUI; increase plot font size; added parameter input
;         check; added buttons to choose method for computing
;         C-value.
;     Ver.5, 17-Aug-2017, Peter Young
;         Fixed bug where the level labels were not synced to
;         transition; added 'Use previous value' button for C-value;
;         added slider and two buttons for manually adjusting
;         C-value.
;     Ver.6, 28-Jul-2020, Peter Young
;         Fixed crash when skipping the final transition; added the
;         widget 'Use high temp limit point?' which gives the user the
;         option of not using the high temperature limit point for
;         allowed transitions; default method for finding C has
;         changed from 0 to 2 ("hybrid maxmin,grad") method, as this
;         generally gives a better spread of upsilons.
;-


PRO gui_set_temp_labels, state, lo_label=lo_label, hi_label=hi_label
;
; This returns the strings that display the minimum and maximum
; temperatures of scaled upsilon plot.
;
t_ind0=state.wid_data.t_ind0
t_ind1=state.wid_data.t_ind1
temp=state.upsstr.data[state.wid_data.index].temp

lo_label=string(format='("Low temp: ",e9.2)',temp[t_ind0])
hi_label=string(format='("High temp: ",e9.2)',temp[t_ind1])

END


;-------------------
PRO gui_write_splstr, state
;
; This writes the spline structure to the .scups file.
;

splfile=state.wid_data.splfile
upsstr=state.upsstr
ind=state.wid_data.index


chck=file_search(splfile,count=count)
IF count NE 0 THEN BEGIN
 ;
 ; SCUPS file already exists.
 ; --------------------------
 ; First I make a backup of the existing file
  backup_file=splfile+'_backup2'
  file_copy,splfile,backup_file,/overwrite
 ;
 ; Here I read scups data from the most recent scups file, and then add
 ; an extra entry (using the /add_empty) so that the new transition
 ; will go in the empty row (index i).
 ;
  read_scups,splfile,splstr,/add_empty
  i=splstr.info.ntrans-1
 ;
 ; Transfer the uppstr data into the empty row of SPLSTR.
 ;
  splstr.data[i].lvl1=upsstr.data[ind].lvl1
  splstr.data[i].lvl2=upsstr.data[ind].lvl2
  splstr.data[i].de=upsstr.data[ind].de
  splstr.data[i].gf=upsstr.data[ind].gf
  splstr.data[i].lim=upsstr.data[ind].lim
  splstr.data[i].c_ups=state.wid_data.c_ups
  splstr.data[i].t_type=state.wid_data.t_type
 ;
 ; Now load the scaled temperatures and upsilons, taking care to deal
 ; with missing data.
 ;
  st=state.wid_data.st
  sups=state.wid_data.sups
  missing=state.wid_data.missing
  k=where(st NE missing,nk)
  splstr.data[i].stemp[0:nk-1]=st[0:nk-1]
  splstr.data[i].spl[0:nk-1]=sups[0:nk-1]
  splstr.data[i].nspl=nk
 ;
ENDIF ELSE BEGIN
 ;
 ; SCUPS file does not exist.
 ; --------------------------
 ; Need to create a fresh SPLSTR using the current transition.
 ;
  info={ion_name: state.wid_data.ionname, $
        ntrans: 1l, $
        comments: upsstr.info.comments}
 ;
 ; The following is to make sure I get the size of STEMP and SPL
 ; correct, based on the maximum possible size (given information in
 ; UPSSTR). 
  max_nspl=max(upsstr.data.nt)+2
  st=state.wid_data.st
  sups=state.wid_data.sups
  missing=state.wid_data.missing
  k=where(st NE missing,nspl)
  
  data={ lvl1: upsstr.data[ind].lvl1, $
         lvl2: upsstr.data[ind].lvl2, $
         de: upsstr.data[ind].de, $
         gf: upsstr.data[ind].gf, $
         lim: upsstr.data[ind].lim, $
         c_ups: state.wid_data.c_ups, $
         t_type: state.wid_data.t_type, $
         nspl: nspl, $
         stemp: fltarr(max_nspl), $
         spl: fltarr(max_nspl)}
  data.stemp[0:nspl-1]=st[0:nspl-1]
  data.spl[0:nspl-1]=sups[0:nspl-1]
  splstr={info: info, data: data}
ENDELSE 

;
; Write out new SCUPS file.
;
write_scups_struc,splstr,/overwrite

END




;-----------------------
function gui_info_string, state
;
; This creates the information string displayed in the GUI.
;

ind=state.wid_data.index
ntrans=state.upsstr.info.ntrans
n_remain=ntrans-ind
elvlcstr=state.elvlcstr

str=[ 'Ion: '+state.upsstr.info.ion_roman, $
      '', $
      'No. of transitions to fit: '+trim(ntrans), $
      'No. of transitions remaining: '+trim(n_remain), $
      'No. of transitions completed: '+trim(state.wid_data.fit_count), $
      'No. of transitions skipped: '+trim(state.wid_data.skip_count) ]


IF ind LT ntrans THEN BEGIN
  i=state.upsstr.data[ind].lvl1
  j=state.upsstr.data[ind].lvl2
  str=[str, $
       '', $
       'Current transition: '+trim(i)+'-'+trim(j)]
 ;
  IF n_tags(elvlcstr) NE 0 THEN BEGIN
    k=where(elvlcstr.data.index EQ i)
    str=[str,'  Level '+trim(i)+': '+elvlcstr.data[k[0]].full_level ]
   ;
    k=where(elvlcstr.data.index EQ j)
    str=[str,'  Level '+trim(j)+': '+elvlcstr.data[k[0]].full_level ]
  ENDIF
ENDIF 

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
IF ind EQ state.upsstr.info.ntrans THEN BEGIN
  plot,/nodata,[0,1],[0,1],xticklen=1e-6,yticklen=1e-6, $
       charsize=1e-6
  xyouts,align=0.5,0.5,0.5,'All transitions completed',charsize=2.0
ENDIF ELSE BEGIN 
  l1=state.upsstr.data[ind].lvl1
  l2=state.upsstr.data[ind].lvl2

  title=trim(l1)+'-'+trim(l2)

  plot,st,sups,xra=[-0.05,1.05],psym=6,title=title,/xsty, $
       yra=[0,max(sups)*1.10],/ysty, $
       xtitle='Scaled temperature', $
       ytitle='Scaled upsilon',charsize=1.4
  oplot,[0,0],[-1e5,1e5],line=2
  oplot,[1,1],[-1e5,1e5],line=2

  y2=spl_init(st,sups)
  sti=findgen(51)/50.
  supsi=spl_interp(st,sups,y2,sti)
  oplot,sti,supsi
ENDELSE 

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

widget_control,state.droplim_but,get_val=use_lim

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
sups=[max([cc[0],0.]),sups]

;
; extrapolate to one (if necessary)
;
IF udata.lim EQ state.upsstr.info.missing OR use_lim EQ 0 THEN BEGIN
  n=n_elements(sups)
  cc=linfit(st[n-2:n-1],sups[n-2:n-1])
  st=[st,1.]
  sups=[sups,max([cc[0]+cc[1]*1.0,0])]
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

;
; This deals with the option for including the high-temp limit for
; allowed transitions.
;
IF t_type EQ 1 OR t_type EQ 4 THEN BEGIN
  sens=1
  use_lim=1
ENDIF ELSE BEGIN
  sens=0
ENDELSE 
widget_control,state.droplim_base,sens=sens
widget_control,state.droplim_but,set_val=use_lim


nt=state.upsstr.data[ind].nt
temp=state.upsstr.data[ind].temp[0:nt-1]
ups=state.upsstr.data[ind].ups[0:nt-1]
c_method=state.wid_data.c_method
c_ups=burly_optimize_c(temp,ups,de,t_type,method=c_method)
state.wid_data.c_ups=c_ups

state.wid_data.t_ind0=0
state.wid_data.t_ind1=nt-1

widget_control,state.lo_slider,set_val=0
widget_control,state.lo_slider,set_slider_max=nt-1
widget_control,state.hi_slider,set_val=nt-1
widget_control,state.hi_slider,set_slider_max=nt-1
widget_control,state.type_buts,set_val=t_type-1
widget_control,state.c_over_val,set_val=trim(string(format='(f10.4)',c_ups))

IF t_type EQ 1 OR t_type EQ 4 THEN widget_control,state.droplim_but,set_val=1

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
     ;
      0: BEGIN                   ; --- SKIP TRANSITION
        IF state.wid_data.index LT state.upsstr.info.ntrans THEN BEGIN
          state.wid_data.skip_count=state.wid_data.skip_count+1
          state.wid_data.index=state.wid_data.index+1
          str=gui_info_string(state)
          widget_control,state.info_txt,set_val=str
          widget_control,state.scups_base,set_uval=state
        ENDIF
       ;
        IF state.wid_data.index EQ state.upsstr.info.ntrans THEN BEGIN
          widget_control,state.lo_slider,sensitive=0
          widget_control,state.hi_slider,sensitive=0
          widget_control,state.type_buts,sensitive=0
          widget_control,state.c_over_val,sensitive=0
          gui_plot_spl,state
          result=dialog_message('All transitions have been processed. Please exit the program.')
        ENDIF ELSE BEGIN
          scups_load_trans,state
          gui_scale_ups, state
          gui_plot_spl, state
          widget_control,state.info_txt,set_val=str
          widget_control,state.scups_base,set_uval=state
        ENDELSE 
      END 
     ;
      1: BEGIN                   ; --- NEXT TRANSITION
       ;
       ; Write the current transition to the SCUPS file.
       ;
        IF state.wid_data.index LT state.upsstr.info.ntrans THEN BEGIN
          widget_control,/hourglass
          state.wid_data.fit_count=state.wid_data.fit_count+1
          gui_write_splstr, state
          state.wid_data.index=state.wid_data.index+1
          str=gui_info_string(state)
          widget_control,state.scups_base,set_uval=state
        ENDIF 
       ;
       ; Check if we've done all the transitions
       ;
        IF state.wid_data.index EQ state.upsstr.info.ntrans THEN BEGIN
          widget_control,state.lo_slider,sensitive=0
          widget_control,state.hi_slider,sensitive=0
          widget_control,state.type_buts,sensitive=0
          widget_control,state.c_over_val,sensitive=0
          gui_plot_spl,state
          result=dialog_message('All transitions have been processed. Please exit the program.')
        ENDIF ELSE BEGIN
         ;
         ; Load the next transition
         ;
          state.wid_data.prev_c_ups=state.wid_data.c_ups
          scups_load_trans,state
          gui_scale_ups, state
          gui_plot_spl, state
          widget_control,state.scups_base,set_uval=state
          str=gui_info_string(state)
          widget_control,state.info_txt,set_val=str
        ENDELSE 
      END
     ;
      2: BEGIN                   ; --- EXIT
       ;
       ; If the user exits without having fitted a transition, then don't
       ; want to save the spline values. Check if l1=0.
       ;
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
    ups=state.upsstr.data[ind].ups[t_ind0:t_ind1]
    t_type=state.wid_data.t_type
    c_method=state.wid_data.c_method
    c_ups=burly_optimize_c(temp,ups,de,t_type,method=c_method)
    state.wid_data.c_ups=c_ups
    widget_control,state.c_over_val,set_val=trim(string(format='(f10.4)',c_ups))
   ;
    gui_set_temp_labels,state,lo_label=lo_label,hi_label=hi_label
    widget_control,state.lo_label,set_val=lo_label
    widget_control,state.hi_label,set_val=hi_label
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
    ups=state.upsstr.data[ind].ups[t_ind0:t_ind1]
    t_type=state.wid_data.t_type
    c_method=state.wid_data.c_method
    c_ups=burly_optimize_c(temp,ups,de,t_type,method=c_method)
    state.wid_data.c_ups=c_ups
    widget_control,state.c_over_val,set_val=trim(string(format='(f10.4)',c_ups))
   ;
    gui_set_temp_labels,state,lo_label=lo_label,hi_label=hi_label
    widget_control,state.lo_label,set_val=lo_label
    widget_control,state.hi_label,set_val=hi_label
   ;
    widget_control,state.scups_base,set_uval=state
    gui_scale_ups, state
    gui_plot_spl, state
  END 

  event.id EQ state.type_buts: BEGIN
    widget_control,state.type_buts,get_val=val
    t_type=val+1
    state.wid_data.t_type=t_type
    IF t_type NE 1 AND t_type NE 4 THEN BEGIN
      widget_control,state.droplim_base,sens=0
    ENDIF ELSE BEGIN 
      widget_control,state.droplim_base,sens=1
    ENDELSE 
   ;
   ; Update scaling parameter
   ;
    t_ind0=state.wid_data.t_ind0
    t_ind1=state.wid_data.t_ind1
    ind=state.wid_data.index
    de=state.upsstr.data[ind].de
    nt=state.upsstr.data[ind].nt
    temp=state.upsstr.data[ind].temp[t_ind0:t_ind1]
    ups=state.upsstr.data[ind].ups[t_ind0:t_ind1]
    c_method=state.wid_data.c_method
    c_ups=burly_optimize_c(temp,ups,de,t_type,method=c_method)
    state.wid_data.c_ups=c_ups
    widget_control,state.c_over_val,set_val=trim(string(format='(f10.4)',c_ups))
   ;
    widget_control,state.scups_base,set_uval=state
    gui_scale_ups, state
    gui_plot_spl, state
  END 

  event.id EQ state.droplim_but: BEGIN
    widget_control,state.droplim_but,get_val=use_lim
    gui_scale_ups,state
    gui_plot_spl,state
    widget_control,state.scups_base,set_uval=state    
  END 
  
  event.id EQ state.c_over_val: BEGIN
    widget_control,state.c_over_val,get_val=val
    state.wid_data.c_ups=val
    widget_control,state.scups_base,set_uval=state
    gui_scale_ups, state
    gui_plot_spl, state
  END

  event.id EQ state.c_options: BEGIN
    widget_control,state.c_options,get_val=val
    state.wid_data.c_method=val
   ;
    t_type=state.wid_data.t_type
    t_ind0=state.wid_data.t_ind0
    t_ind1=state.wid_data.t_ind1
    ind=state.wid_data.index
    de=state.upsstr.data[ind].de
    nt=state.upsstr.data[ind].nt
    temp=state.upsstr.data[ind].temp[t_ind0:t_ind1]
    ups=state.upsstr.data[ind].ups[t_ind0:t_ind1]

    c_ups=burly_optimize_c(temp,ups,de,t_type,method=val)

    state.wid_data.c_ups=c_ups
    widget_control,state.c_over_val,set_val=trim(string(format='(f10.4)',c_ups))
    c_range=state.wid_data.c_range
    getmin=min(abs(c_range-state.wid_data.c_ups),imin)
    widget_control,state.c_slider,set_value=imin
   ;
    widget_control,state.scups_base,set_uval=state
    gui_scale_ups, state
    gui_plot_spl, state
  END

 ;
 ; Button to select the C-value from the previous transition.
 ;
  event.id EQ state.c_prev: BEGIN
    prev_c_ups=state.wid_data.prev_c_ups
    IF prev_c_ups NE -1. THEN BEGIN 
      state.wid_data.c_ups=prev_c_ups
      widget_control,state.c_over_val,set_val=trim(string(format='(f10.4)',state.wid_data.c_ups))
      c_range=state.wid_data.c_range
      getmin=min(abs(c_range-state.wid_data.c_ups),imin)
      widget_control,state.c_slider,set_value=imin
      widget_control,state.scups_base,set_uval=state    
      gui_scale_ups, state
      gui_plot_spl, state
    ENDIF 
  END 

  event.id EQ state.c_slider: BEGIN
    widget_control,state.c_slider,get_value=i
    state.wid_data.c_ups=state.wid_data.c_range[i]
    widget_control,state.c_over_val,set_val=trim(string(format='(f10.4)',state.wid_data.c_ups))
    widget_control,state.scups_base,set_uval=state    
    gui_scale_ups, state
    gui_plot_spl, state
  END 

 ;
 ; This button moves the C-value slider backwards by two steps.
 ;
  event.id EQ state.c_slider_butt1: BEGIN
    c_range=state.wid_data.c_range
    c_ups=state.wid_data.c_ups
    getmin=min(abs(c_range-c_ups),imin)
    IF imin GT 1 THEN i=imin-2 ELSE i=0
    state.wid_data.c_ups=c_range[i]
    widget_control,state.c_over_val,set_val=trim(string(format='(f10.4)',state.wid_data.c_ups))
    widget_control,state.c_slider,set_val=i
    widget_control,state.scups_base,set_uval=state    
    gui_scale_ups, state
    gui_plot_spl, state    
  END 

 ;
 ; This button moves the C-value slider forward by two steps.
 ;
  event.id EQ state.c_slider_butt2: BEGIN
    c_range=state.wid_data.c_range
    nc=n_elements(c_range)
    c_ups=state.wid_data.c_ups
    getmin=min(abs(c_range-c_ups),imin)
    IF imin LE nc-3 THEN i=imin+2 ELSE i=nc-1
    state.wid_data.c_ups=c_range[i]
    widget_control,state.c_over_val,set_val=trim(string(format='(f10.4)',state.wid_data.c_ups))
    widget_control,state.c_slider,set_val=i
    widget_control,state.scups_base,set_uval=state    
    gui_scale_ups, state
    gui_plot_spl, state    
  END 

ENDCASE 

END 


;----------------------
PRO manu_scups_widget, group=group, upsstr=upsstr, ionname=ionname, $
                       elvlcstr=elvlcstr,  prev_c_ups=prev_c_ups
;
; UPSSTR is the upsilon structure containing data only for those
;        transitions that are not already in the SCUPS file.
;
missing_val=-1

nmax=max(upsstr.data.nt)+2

splfile=trim(upsstr.info.ion_name)+'.scups'

wid_data={ spl_id: 0, $
           splfile: splfile, $
           ionname: ionname, $
           st: fltarr(nmax), $
           sups: fltarr(nmax), $
           t_type: 0, $
           c_ups: 0., $
           c_range: fltarr(101), $
           prev_c_ups: prev_c_ups, $
           c_method:  2, $   ; set to hybrid gradient-maxmin method
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
               ['SKIP','WRITE and NEXT','EXIT'], $
               /row,font=bigfont)

sub_base1=widget_base(base_left,/row,map=1)


sub_base5=widget_base(sub_base1,/col,map=1)
info_txt=widget_text(sub_base5,val='',ysiz=12,xsiz=40,font=font) 


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
lo_label=widget_label(sub_base3,val='',font=font,/align_left,/dynamic_resize)
hi_label=widget_label(sub_base3,val='',font=font,/align_left,/dynamic_resize)



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

IF t_type EQ 1 OR t_type EQ 4 THEN sens=1 ELSE sens=0
droplim_base=widget_base(sub_base4,/col,map=1,/frame,sens=sens)
droplim_txt=widget_label(droplim_base,value='Use high temp limit point?',font=font)
droplim_but=cw_bgroup(droplim_base,/row,['No','Yes'],set_val=1,font=font,/exclusive)


de=upsstr.data[0].de
nt=upsstr.data[0].nt
temp=upsstr.data[0].temp[0:nt-1]
ups=upsstr.data[0].ups[0:nt-1]
c_ups=burly_optimize_c(temp,ups,de,t_type,method=wid_data.c_method,c_range=c_range)
wid_data.c_range=c_range
nc=n_elements(c_range)
getmin=min(abs(c_range-c_ups),imin)
wid_data.c_ups=c_ups
sub_base_c=widget_base(sub_base1,/col,map=1,frame=1)
c_slider_base=widget_base(sub_base_c,/row,map=1)
c_slider=widget_slider(c_slider_base,min=0,max=nc-1,value=imin,font=font,xsize=150)
c_slider_butt1=widget_button(c_slider_base,value='<',font=bigfont)
c_slider_butt2=widget_button(c_slider_base,value='>',font=bigfont)
sub_base_cval=widget_base(sub_base_c,/row,map=1,frame=1)
c_over_val_lbl=widget_label(sub_base_cval,val='C-value:',font=bigfont)
c_over_val=widget_text(sub_base_cval,value=trim(string(format='(f10.4)',c_ups)), $
                       /editable,font=bigfont, $
                       xsiz=10,sens=1)
c_prev=widget_button(sub_base_c,value='Use previous value',font=font)
c_options=cw_bgroup(sub_base_c,['Minimize gradients', $
                                'Max (min step)', $
                                'Hybrid (maxmin,grad)', $
                                'Min (max step)', $
                                'Hybrid (minmax,grad)'], $
                    /col,font=font,set_val=wid_data.c_method, $
                    /exclusive )


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
       elvlcstr: elvlcstr, $
       exit: exit, $
       spl_plot: spl_plot, $
       info_txt: info_txt, $
       scups_base: scups_base, $
       lo_slider: lo_slider, $
       hi_slider: hi_slider, $
       type_buts: type_buts, $
       c_over_val: c_over_val, $
       lo_label: lo_label, $
       hi_label: hi_label, $
       c_options: c_options, $
       c_prev: c_prev, $
       c_slider: c_slider, $
       c_slider_butt1: c_slider_butt1, $
       c_slider_butt2: c_slider_butt2, $
       droplim_but: droplim_but, $
       droplim_base: droplim_base }


WIDGET_CONTROL, scups_base, /REALIZE, set_uvalue=state

WIDGET_CONTROL, spl_plot, GET_VALUE=spl_id

state.wid_data.spl_id=spl_id

WIDGET_CONTROL, scups_base, /REALIZE, set_uvalue=state

str=gui_info_string(state)
widget_control,state.info_txt,set_val=str

gui_scale_ups,state
gui_plot_spl,state

;
; Set the min and max temperature labels.
;
gui_set_temp_labels,state,lo_label=lo_label,hi_label=hi_label
WIDGET_CONTROL, state.lo_label, set_value=lo_label
WIDGET_CONTROL, state.hi_label, set_value=hi_label

XMANAGER, 'scups_base', scups_base, group=group


END


;------------------
PRO write_scups_gui, ionname, no_flag=no_flag, elvlcfile=elvlcfile
                    


IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> write_scups_gui, ion_name [ ,/no_flag ]'
  return
ENDIF 

upsfile=ionname+'.ups'
chck=file_search(upsfile,count=count)
IF count EQ 0 THEN BEGIN
  print,'%WRITE_SCUPS_GUI: The .ups file ('+upsfile+') was not found. Returning...'
  return
ENDIF 
read_ups,upsfile,upsstr

ntrans=upsstr.info.ntrans
k=lonarr(ntrans)+1    ; initially select all transitions in upsstr.data


;
; Now read elvlc file
;
IF n_elements(elvlcfile) EQ 0 THEN elvlcfile=ionname+'.elvlc'
chck=file_search(elvlcfile,count=count)
IF count NE 0 THEN BEGIN
  read_elvlc,elvlcfile,elvlcstr=elvlcstr
ENDIF ELSE BEGIN
  print,'%WRITE_SCUPS_GUI: The .elvlc file ('+elvlcfile+') was not found. Transition information will not be displayed.'
  elvlcstr=0
ENDELSE 


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
; Note that I extract the previous scaling parameter (prev_c_ups)
; which is needed by the button 'Use previous value'.
;
chck=file_search(splfile)
prev_c_ups=-1.
IF chck[0] NE '' THEN BEGIN 
  read_scups,splfile,splstr
 ;
 ; Flag as 0 any transitions that are already in the scups file.
 ;
  FOR i=0,ntrans-1 DO BEGIN
    j=where(upsstr.data[i].lvl1 EQ splstr.data.lvl1 AND upsstr.data[i].lvl2 EQ splstr.data.lvl2,nj)
    IF nj GT 0 AND k[i] EQ 1 THEN k[i]=0
  ENDFOR
  n=splstr.info.ntrans
  prev_c_ups=splstr.data[n-1].c_ups
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
manu_scups_widget, ionname=ionname, upsstr=new_ups, elvlcstr=elvlcstr, prev_c_ups=prev_c_ups


END
