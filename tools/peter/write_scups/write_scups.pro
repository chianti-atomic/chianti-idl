
;+
; NAME:
;     WRITE_SCUPS
;
; CATEGORY:
;     CHIANTI; data assessement.
;
; PURPOSE:
;     This routine reads the specified CHIANTI UPS file format and
;     writes out a CHIANTI SCUPS file containing spline fits to the
;     effective collision strength data. The routine checks the
;     upsilons for various problems. If any are found, then those
;     transitions are not written to the output SCUPS file.
;
;     The routine optionally creates a plot file that can display all
;     problem transitions, and all good successful transitions.
;
;     A separate document (CHIANTI technical report #4) gives more
;     detail about how WRITE_SCUPS works.
;
; INPUTS:
;     Ion_Name: The name of the ion in the standard CHIANTI format
;               (e.g., o_6 for O VI). Note that the routine expects to
;               find the file ION_NAME.ups in the user's
;               working directory. 
;
; OPTIONAL INPUTS:
;     Plot_prob_type:  Takes an integer between 1 and 6 corresponding
;                 to the different problem types (see Technical Report
;                 No. 4). All problem transitions of this type will be
;                 plotted to mg_5_plots.ps. If no problems are found
;                 then the plot will be empty.
;     Plot_trans_type: Takes an integer between 1 and 6 corresponding
;                 to the transition type with the usual CHIANTI
;                 indexing. All transitions of the specified type will
;                 be plotted whether they are problem transitons or
;                 not.
;     Params:     A four element array that controls how strict the
;                 routine is when checking for bad points. Please
;                 check Technical Report No. 4 for more details. 
;
; KEYWORDS PARAMETERS:
;     IGNORE_LIMIT:  If set, then if the program encounters a
;                    transition that is not tending to its
;                    high-temperature limit point (i.e., problem type
;                    6), then it will ignore the limit point.
;     STRICT:  If set, then the routine sets a very strict set of
;              PARAMS values when checking for bad points.
;
; OUTPUTS:
;     Creates the file [ionname].scups_auto containing the scaled upsilon
;     data. ([ionname] is the CHIANTI format for an ion name, e.g.,
;     'o_6' is O VI.)
;
;     The text files [ionname]_flags.txt and [ionname]_log.txt are
;     also created.
;
;     A postscript file containing plots of various transition can be
;     produced if one of the plot keywords is set.
;
; EXAMPLES:
;     IDL> write_scups, 'mg_5.ups'
;     IDL> write_scups, 'mg_5.ups', plot_prob_type=6
;
; CALLS:
;     READ_UPS, CONVERTNAME, BURLY_GET_TTYPE, BURLY_OPTIMIZE_C,
;     BURLY_SCALE_UPS
;
; INTERNAL ROUTINES:
;     UPDATE_FLAGSTR.
;
; MODIFICATION HISTORY:
;     Ver.1, 21-May-2013, Peter Young
;     Ver.2, 24-May-2013, Peter Young
;        Tidied up some of the plots.
;     Ver.3, 24-May-2013, Peter Young
;        Modified output format to be more consistent with upsdatx
;        file. 
;     Ver.4, 31-May-2013, Peter Young
;        Changed routine name to reflect new names for the .ups and
;        .scups files.
;     Ver.5, 31-May-2013, Peter Young
;        Tweaked the parameters for type 4 and type 6 problems in
;        order to reduce the number of transitions that are flagged.
;     Ver.6, 30-Jan-2017, Peter Young
;        Now adds '-1' to the end of the file, consistent with updated
;        file format.
;-



;----------------------
FUNCTION update_flagstr, flagstr, l1, l2, flag=flag

IF n_elements(flag) EQ 0 THEN flag=1

IF n_tags(flagstr) EQ 0 THEN BEGIN
  flagstr={l1: l1, l2: l2, flag: flag}
ENDIF ELSE BEGIN
  i=where(flagstr.l1 EQ l1 AND flagstr.l2 EQ l2,nf)
  IF nf EQ 0 THEN BEGIN
    fstr={l1: l1, l2:l2, flag: flag}
    flagstr=[flagstr,fstr]
  ENDIF 
ENDELSE

return,flagstr

END


;----------------
PRO write_scups, ion_name, plot_file=plot_file, plot_type=plot_type, $
                 problem=problem, all=all, limit=limit, no_limit=no_limit, $
                 plot_prob_type=plot_prob_type, plot_trans_type=plot_trans_type, $
                 ignore_limit=ignore_limit, strict=strict, params=params


IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> write_scups, ion_name [, plot_prob_type=, plot_trans_type, /ignore_limit'
  print,'                          params=, /strict ]'
  print,"  ion_name - standard CHIANTI format, e.g., 'fe_13' for Fe XIII"
  return
ENDIF 

;
; The array PARAMS controls how strict the routine is when working out
; whether there are any irregular points. 
;
IF n_elements(params) EQ 0 THEN BEGIN
  IF keyword_set(strict) THEN BEGIN
    params=[ 1.0, 1.0, 1.0, 0.0, 0.1]
  ENDIF ELSE BEGIN
    params=[ 2.5, 2.5, 1.5, 0.1, 0.3]
  ENDELSE 
ENDIF ELSE BEGIN
  np=n_elements(params)
  IF np NE 4 THEN BEGIN
    print,'%WRITE_SCUPS: the input PARAMS must have 4 elements.'
    return
  ENDIF 
ENDELSE

;
; Read the UPS file and extract information.
;
upsfile=ion_name+'.ups'
chck=file_search(upsfile,count=count)
IF count EQ 0 THEN BEGIN
  print,'%WRITE_SCUPS: the file '+upsfile+' was not found. Returning...'
  return
ENDIF 
read_ups, upsfile, upsstr
n=upsstr.info.ntrans
ionname=upsstr.info.ion_name

;
; Set up output files
;
IF n_elements(outfile) EQ 0 THEN outfile=ionname+'.scups_auto'
IF n_elements(plot_file) EQ 0 THEN plot_file=ionname+'_plots.ps'
logfile=ionname+'_log.txt'

IF n_elements(plot_type) EQ 0 THEN BEGIN
  plot_type=0
  IF keyword_set(problem) THEN plot_type=1
  IF keyword_set(all) THEN plot_type=2
  IF keyword_set(limit) THEN plot_type=3
  IF keyword_set(no_limit) THEN plot_type=4
ENDIF 

;
; These two keywords control what gets plotted. They are not
; compatible with each other. A -1 indicates the keyword will be
; ignored. 
;
IF n_elements(plot_trans_type) EQ 0 THEN plot_trans_type=-1
IF n_elements(plot_prob_type) EQ 0 THEN plot_prob_type=-1
;
IF plot_trans_type GT 0 AND plot_prob_type GT 0 THEN BEGIN
  print,'%WRITE_SCUPS: Please specify only one of PLOT_PROB_TYPE or PLOT_TRANS_TYPE, not both.'
  return
ENDIF


;
; If necessary, open the plot file.
;
IF plot_trans_type GT 0 OR plot_prob_type GT 0 THEN BEGIN 
  dname=!d.name
  set_plot,'ps'
  device,file=plot_file,xsiz=8,ysiz=10,/inches,yoff=0.5,xoff=0.25
  !p.multi=[0,2,4]
ENDIF


;
; Open the output files.
;
openw,lout,outfile,/get_lun
openw,llog,logfile,/get_lun

;
; The flag structure (FSTR) keeps track of problem transitions. If
; fstr.flag is 0, then there's no problem. Otherwise, it takes
; a value from 1 to 6.
;
fstr={l1: 0, l2: 0, flag: 0}
flagstr=0

ttype_count=intarr(10)
prob_cnt=intarr(6)

;
; This loop goes through all transitions in UPSSTR.
;
FOR i=0,n-1 DO BEGIN
  prob=0b  
 ;
  l1=upsstr.data[i].lvl1
  l2=upsstr.data[i].lvl2
  nt=upsstr.data[i].nt
  temp=upsstr.data[i].temp[0:nt-1]
  ups=upsstr.data[i].ups[0:nt-1]
  de=upsstr.data[i].de
  gf=upsstr.data[i].gf
  lim=upsstr.data[i].lim

  IF keyword_set(no_lim) THEN lim=-1
 ;
 ; Get transition type and c-value
 ;
  ttype=burly_get_ttype(gf,de)
  ttype_count[ttype-1]=ttype_count[ttype-1]+1
  c_val=burly_optimize_c(temp,ups,de,ttype)
;  c_val=burly_find_mid_c(temp,de,ttype)
 ;
 ; Check for problem type 1 (zero energy).
 ; ------------------------
  IF de EQ 0. THEN BEGIN 
    prob=1b
    prob_cnt[0]=prob_cnt[0]+1
    printf,llog,'Transition: '+trim(l1)+'-'+trim(l2)+'; ups index: '+trim(i)+' -- *Zero energy*'
    flagstr=update_flagstr(flagstr,l1,l2,flag=1)
  ENDIF 

 ;
 ; Check for problem type 2 (negative upsilons).
 ; ------------------------
  k=where(ups LT 0.,nk)
  IF nk GT 0 THEN BEGIN
    prob=1b
    prob_cnt[1]=prob_cnt[1]+1
    printf,llog,'Transition: '+trim(l1)+'-'+trim(l2)+'; ups index: '+trim(i)+' -- *Negative upsilon*'
    flagstr=update_flagstr(flagstr,l1,l2,flag=2)
  ENDIF 
 ;
 ; Scale upsilons
 ;
  sclups=burly_scale_ups(temp, ups, de, ttype, c_val)  
  st=sclups.st
  sups=sclups.sups
  nt=n_elements(st)

 ;
 ; Check for problem type 3
 ; ------------------------
 ;   - the 2nd and 3rd points are extrapolated back to the 1st point
 ;     and we look to see if the 1st point is consistent with the
 ;     extrapolated point. Technical Report 4 has more details.
 ;
  cc=linfit(st[1:2],sups[1:2])

  angle=atan(max(sups))*params[0]
 ;
 ; angle_plus can not be greater than +90 degrees (1.57 rads)
 ; angle_minus can not be less than -90 degrees (-1.57 rads)
 ;
  angle_plus=max([-1.57,atan(cc[1])-angle])
  angle_minus=min([1.57,atan(cc[1])+angle])
 ;
  m_plus=tan(angle_plus)
  m_minus=tan(angle_minus)
 ;
 ; lim_range specifies the range of values within which sups[nt-1]
 ; should lie.
 ;
  lim_range=[m_minus*(st[0]-st[1])+sups[1], $
             m_plus*(st[0]-st[1])+sups[1] ]
  lim_range=lim_range[sort(lim_range)]
 ;
 ; Increase lim_range by amount based on max(sups) [see TR4]
  lim_range[0]=lim_range[0]-max(sups)*params[3]
  lim_range[1]=lim_range[1]+max(sups)*params[3]
  lim_range=lim_range>0
 
  IF sups[0] LT lim_range[0] OR sups[0] GT lim_range[1] THEN BEGIN
    prob=1b
    prob_cnt[2]=prob_cnt[2]+1
    printf,llog,'Transition: '+trim(l1)+'-'+trim(l2)+'; ups index: '+trim(i)+' -- *low-T problem*'
    flagstr=update_flagstr(flagstr,l1,l2,flag=3)
    IF (plot_prob_type EQ 3 AND prob EQ 1) OR ttype EQ plot_trans_type THEN BEGIN
      title='*PT=3 (low-T)* Trans:'+trim(l1)+'-'+trim(l2)+' (TT='+trim(ttype)+')'
      plot,st,sups,xra=[-0.05,1.05],/xsty,psym=1,title=title,charsiz=1.5, $
           yra=[0,max([sups,lim_range])*1.10],/ysty
;      plots,psym=6,st[0],chck
      plots,psym=6,st[0],lim_range[0]
      plots,psym=6,st[0],lim_range[1]
      xr=!x.crange
      yr=!y.crange
      xyouts,xr[0]*0.95+xr[1]*0.05,yr[0]*0.95+yr[1]*0.05,'upsstr index='+trim(i)
      oplot,[0,0],[-1e3,1e3],line=2
      oplot,[1,1],[-1e3,1e3],line=2
    ENDIF 
  ENDIF 

 ;
 ; Check for problem type 4
 ; ------------------------
 ;   - the (n-2)th and (n-1)th points are extrapolated forwards to the
 ;     nth point and we look to see if the nth point is consistent with the
 ;     extrapolated point. Technical Report 4 has more details.
 ;
  cc=linfit(st[nt-3:nt-2],sups[nt-3:nt-2])

  angle=atan(max(sups))*params[1]
 ;
 ; angle_plus can not be greater than +90 degrees (1.57 rads)
 ; angle_minus can not be less than -90 degrees (-1.57 rads)
 ;
  angle_plus=min([1.57,atan(cc[1])+angle])
  angle_minus=max([-1.57,atan(cc[1])-angle])
 ;
  m_plus=tan(angle_plus)
  m_minus=tan(angle_minus)
 ;
 ; lim_range specifies the range of values within which sups[nt-1]
 ; should lie.
 ;
  lim_range=[m_minus*(st[nt-1]-st[nt-2])+sups[nt-2], $
             m_plus*(st[nt-1]-st[nt-2])+sups[nt-2] ]
  lim_range=lim_range[sort(lim_range)]
 ;
 ; Increase lim_range by amount based on max(sups) [see TR4]
  lim_range[0]=lim_range[0]-max(sups)*params[3]
  lim_range[1]=lim_range[1]+max(sups)*params[3]
  lim_range=lim_range>0
 ;
 ; Now check if we have a problem and optionally plot the result.
 ;
  IF sups[nt-1] LT lim_range[0] OR sups[nt-1] GT lim_range[1] THEN BEGIN 
   ;
    prob=1b
    prob_cnt[3]=prob_cnt[3]+1
    printf,llog,'Transition: '+trim(l1)+'-'+trim(l2)+'; ups index: '+trim(i)+' -- *high-T problem*'
    flagstr=update_flagstr(flagstr,l1,l2,flag=4)
    IF (plot_prob_type EQ 4 AND prob EQ 1) OR ttype EQ plot_trans_type THEN BEGIN 
      title='*PT=4 (high-T)* Trans:'+trim(l1)+'-'+trim(l2)+' (TT='+trim(ttype)+')'
      plot,st,sups,xra=[-0.05,1.05],/xsty,psym=1,title=title,charsiz=1.5, $
           yra=[0,max([sups,lim_range])*1.10],/ysty
      plots,psym=6,st[nt-1],lim_range[1]
      plots,psym=6,st[nt-1],lim_range[0]
      xr=!x.crange
      yr=!y.crange
      xyouts,xr[0]*0.95+xr[1]*0.05,yr[0]*0.95+yr[1]*0.05,'upsstr index='+trim(i)
      oplot,[0,0],[-1e3,1e3],line=2
      oplot,[1,1],[-1e3,1e3],line=2
    ENDIF 
  ENDIF 

 ;
 ; Check for problem type 5
 ; ------------------------
 ;    - For each of the "intermediate" points (i.e., not the first or
 ;      last), a spline is fit through all points but the intermediate
 ;      point.
 ;    - If the difference is > (params[5]*100) % then the transition
 ;      is flagged. 
 ;
  uu=sups
  tt=st
  u_spl=fltarr(nt)
  nu=nt-2
  FOR j=0,nu-1 DO BEGIN
    ux=[uu[0:j],uu[j+2:*]]
    tx=[tt[0:j],tt[j+2:*]]
    y2=spl_init(tx,ux)
    yi=spl_interp(tx,ux,y2,tt[j+1])
    IF abs(yi-uu[j+1])/uu[j+1] GE params[4] THEN u_spl[j+1]=yi
  ENDFOR
  k=where(u_spl NE 0,nk)
  IF nk GT 0 THEN BEGIN
    prob=1b
    prob_cnt[4]=prob_cnt[4]+1
    printf,llog,'Transition: '+trim(l1)+'-'+trim(l2)+'; ups index: '+trim(i)+' -- *spline anomaly*'
    flagstr=update_flagstr(flagstr,l1,l2,flag=5)
    IF plot_prob_type EQ 5 OR plot_trans_type EQ ttype THEN BEGIN
      title='*PT=5 (spline)* Trans:'+trim(l1)+'-'+trim(l2)+' (TT='+trim(ttype)+')'
      plot,st,sups,xra=[-0.05,1.05],/xsty,psym=1,title=title,charsiz=1.5, $
           yra=[0,max([sups,lim])*1.10],/ysty
      oplot,psym=6,tt[k],u_spl[k]
      xr=!x.crange
      yr=!y.crange
      xyouts,xr[0]*0.95+xr[1]*0.05,yr[0]*0.95+yr[1]*0.05,'upsstr index='+trim(i)
      oplot,[0,0],[-1e3,1e3],line=2
      oplot,[1,1],[-1e3,1e3],line=2
   ENDIF 
  ENDIF 


 ;
 ; Extrapolate to zero
 ; -------------------
 ; Draw a line between the lowest two points to determine the scaled
 ; upsilon value at st=0. If the extrapolated value is less than zero,
 ; then set it to zero.
 ;
  cc=linfit(st[0:1],sups[0:1])
  IF cc[0] LT 0. THEN cc[0]=0.
  st=[0,st]
  sups=[cc[0],sups]
  nt=n_elements(st)

 ; Extrapolate to one
 ; ------------------
 ; Draw a line between the highest two points to determine the scaled
 ; upsilon value at st=1. If the extrapolated value is less than zero,
 ; then set it to zero. Only do extrapolation if the limit value has
 ; not been defined.
 ;
  IF lim EQ -1 THEN BEGIN
  ;
  ; This is the case where there is no limit point.
  ;
    cc=linfit(st[nt-2:nt-1],sups[nt-2:nt-1])
    extrap1=cc[1]+cc[0]
    IF extrap1 LT 0. THEN extrap1=0.
    st=[st,1.]
    sups=[sups,extrap1]
   ;
  ENDIF ELSE BEGIN
   ;
   ; This is the case where there is a limit point.
   ;
    cc=linfit(st[nt-2:nt-1],sups[nt-2:nt-1])
    extrap1=cc[1]+cc[0]
    extrap1=extrap1>0
   ;
    angle=atan(max(sups))*params[2]
   ;
   ; angle_plus can not be greater than +90 degrees (1.57 rads)
   ; angle_minus can not be less than -90 degrees (-1.57 rads)
   ;
    angle_plus=min([1.57,atan(cc[1])+angle])
    angle_minus=max([-1.57,atan(cc[1])-angle])
   ;
    m_plus=tan(angle_plus)
    m_minus=tan(angle_minus)
   ;
   ; lim_range specifies the range of values within which sups[nt-1]
   ; should lie.
   ;
    lim_range=[m_minus*(1.0-st[nt-1])+sups[nt-1], $
               m_plus*(1.0-st[nt-1])+sups[nt-1] ]
    lim_range=lim_range[sort(lim_range)]
    lim_range[0]=lim_range[0]-max(sups)*params[3]
    lim_range[1]=lim_range[1]+max(sups)*params[3]
    lim_range=lim_range>0

    IF lim LT lim_range[0] OR lim GT lim_range[1] AND NOT keyword_set(ignore_limit) THEN BEGIN
      prob=1b
      prob_cnt[5]=prob_cnt[5]+1
      printf,llog,'Transition: '+trim(l1)+'-'+trim(l2)+'; ups index: '+trim(i)+' -- *high-T limit problem*'
      flagstr=update_flagstr(flagstr,l1,l2,flag=6)
      IF plot_prob_type EQ 6 OR plot_trans_type EQ ttype THEN BEGIN
        title='*Limit point* ('+trim(l1)+' - '+trim(l2)+')'
        yra=[0,max([sups,lim])*1.10]
        title='*PT=6 (limit)* Trans:'+trim(l1)+'-'+trim(l2)+' (TT='+trim(ttype)+')'
        plot,st,sups,xra=[-0.05,1.05],/xsty,psym=1,title=title,charsiz=1.5,yra=yra,/ysty
        oplot,[1,1],lim_range,psym=2
        plots,psym=6,1,lim
        xr=!x.crange
        yr=!y.crange
        xyouts,xr[0]*0.95+xr[1]*0.05,yr[0]*0.95+yr[1]*0.05,'upsstr index='+trim(i)
        oplot,[0,0],[-1e3,1e3],line=2
        oplot,[1,1],[-1e3,1e3],line=2
      ENDIF
    ENDIF ELSE BEGIN   ; no problem
      st=[st,1.]
      sups=[sups,lim]
      nt=n_elements(st)
    ENDELSE 
  ENDELSE 

  IF prob EQ 0 THEN BEGIN 
    nt=n_elements(st)
   ;
   ; NOTE: if you adjust the format here, make sure to check against
   ; the routine WRITE_SCUPS_STRUC, as they should be the same.
   ;
   ;  Write out first line of data
   ;
    IF c_val GE 1e4 OR c_val LT 0.1 THEN cform='e12.4' ELSE cform='f12.5'
    IF lim GE 0 THEN limform='(e12.3)' ELSE limform='(i12)'
    format_str='(2i7,2e12.3,'+limform+',2i5,'+cform+')'

;    IF lim GE 0 THEN format_str='(2i7,3e12.3,2i5,e12.5)' ELSE format_str='(2i7,2e12.3,i12,2i5,e12.5)'
    printf,lout,format=format_str, $
           l1,l2,de,gf,lim,nt,ttype,c_val
   ;
   ;  Write out second line of data (scaled temperatures)
   ;
    format_str='('+trim(nt)+'e12.3)'
    printf,lout,format=format_str,st
   ;
   ;  Write out third line of data (scaled upsilons)
   ;
    printf,lout,format=format_str,sups
   ;
   ;  Plot the "good" transitions.
   ;
    swtch=0
    IF plot_trans_type GT 0 THEN begin
    ;; IF plot_type EQ 2 THEN swtch=1
    ;; IF plot_type EQ 3 AND lim GE 0 THEN swtch=1
    ;; IF plot_type EQ 4 AND lim LT 0 THEN swtch=1
      IF ttype EQ plot_trans_type THEN BEGIN 
        title='Trans: '+trim(l1)+'-'+trim(l2)+', Type='+trim(ttype)+ $
              ', C='+trim(string(format='(f10.2)',c_val))
        plot,st,sups,xra=[-0.05,1.05], $
             yra=[0,max([sups,lim])*1.1],psym=1,charsiz=1.5,xsty=1, $
             title=title
        IF lim GE 0 THEN oplot,[1,1],lim_range,psym=2
        oplot,[0,0],[-1e3,1e3],line=2
        oplot,[1,1],[-1e3,1e3],line=2
       ;
       ; Overplot the spline fit to the data as solid line
       ;
        xi=findgen(51)/50.
        y2=spl_init(st,sups)
        yi=spl_interp(st,sups,y2,xi)
        oplot,xi,yi
      ENDIF
    ENDIF 
  ENDIF 

ENDFOR

;
; Add comments to the scups file.
;
printf,lout,' -1'
printf,lout,'File:  '+outfile
printf,lout,'File created: '+systime()
refs=upsstr.info.comments
nrefs=n_elements(refs)
FOR i=0,nrefs-1 DO printf,lout,refs[i]
printf,lout,' -1'


;
; Close the postscript plot.
;
IF plot_prob_type GT 0 OR plot_trans_type GT 0 THEN BEGIN 
  !p.multi=0
  device,/close
  set_plot,dname
ENDIF 

free_lun,lout,llog

;
; Write out the problem transitions to the flag file.
;
flagfile=ionname+'_flags.txt'
IF n_tags(flagstr) NE 0 THEN BEGIN
  openw,lout,flagfile,/get_lun
  n=n_elements(flagstr)
  FOR i=0,n-1 DO BEGIN
    printf,lout,format='(3i7)',flagstr[i].l1,flagstr[i].l2,flagstr[i].flag
  ENDFOR 
  free_lun,lout
ENDIF 

;
; Give a summary of the different transition types.
;
print,''
print,'Transition type summary: '
FOR i=0,9 DO BEGIN
  IF ttype_count[i] NE 0 THEN print,'  - type '+trim(i+1)+': ',trim(ttype_count[i])
ENDFOR 


;
; Print information about any problem transitions.
;
IF total(prob_cnt) NE 0 THEN BEGIN
  print,''
  print,'The following problems were found. These transitions were NOT output to the scups file.'
  print,'PT refers to problem type. See Technical Report No. 4 for more details.'
  prob_str=['PT 1 (zero energy): ', $
            'PT 2 (negative ups): ', $
            'PT 3 (low temp): ', $
            'PT 4 (high temp): ', $
            'PT 5 (spline): ', $
            'PT 6 (allowed limit): ']
  FOR i=0,5 DO BEGIN
    IF prob_cnt[i] NE 0 THEN print,'  - '+prob_str[i]+trim(prob_cnt[i])
  ENDFOR 
  print,''
  print,'Please check the file '+logfile+' for more details.'
  print,'To make plots of the problem transitions, use the PLOT_PROB_TYPE keyword.'
ENDIF


;
; Summarize the files that have been created.
;
print,''
print,'The following files have been created:'
chck=file_search(outfile)
IF chck[0] NE '' THEN print,'  - '+outfile
chck=file_search(flagfile)
IF chck[0] NE '' THEN print,'  - '+flagfile
IF plot_type GT 0 THEN print,'  - '+plot_file
chck=file_search(logfile)
IF chck[0] NE '' THEN print,'  - '+logfile


END
