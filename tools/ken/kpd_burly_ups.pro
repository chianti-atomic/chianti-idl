;+
; NAME
;
;    KPD_BURLY_UPS
;
; PROJECT
;
;    CHIANTI
;
; EXPLANATION
;
;    Widget-based routine for fitting electron excitation data for input
;    to the CHIANTI database. It is an update of the original routine
;    burly_ups.pro. New features include:
;
;      - automatic minimization to determine the best fit scaling
;        parameter, C
;      - routine determines the metastable levels for which fitting is
;        required
;      - transitions are automatically loaded consecutively by level number
;      - all fitted transitions are stored internally in an IDL structure,
;        and only printed when the program is exited. Thus it is not
;        possible for the file to contain two entries for the same transition.
;
;    **ADDING EXTRA SPLINE POINTS
;    This routine can be modified in a fairly straightforward manner to
;    perform spline fits with a greater number of nodes. Note that some
;    code is left over from my experiments with 15 point splines, so this
;    can be used (do a search on "15" to find this code).
;    - check for "xs=" and change this (defines spline values)
;    - change definition of trans.spl to include extra points
;    - change NPT_BUT widget to include extra button (also change uval)
;    - make sure to modify the separate routines burly_fit_fn.pro,
;      descale_all.pro and read_splups.pro
;
;   This version has been renamed from new_burly_ups since it allows the use of 5 through 9
;   spline nodes
;
; INPUTS
;
;    DIRECTORY Location of the CHIANTI input files.
;
;    IZ        Atomic number of ion, e.g., IZ=26 for iron.
;
;    ION       Spectroscopic number of ion, e.g., ION=13 for Fe XIII.
;
; OPTIONAL INPUTS
;
;    CUTOFF    A level is considered to be metastable if it has no
;              A-values greater than CUTOFF. The default is 1e5.
;
;    META      An array of same size as the number of levels in the
;              ion, and with values of either 0 or 1. A 1 indicates
;              that the level is metastable and should be included for
;              spline fits. The best way to specify META is to run the
;              routine metastable_levels to yield a default META array
;              and then edit the entries by hand.
;
; KEYWORDS
;
;    PROTON    Fit proton rates rather than electron rates. NOT IMPLEMENTED
;              YET!
;
;    DIEL      Fit dielectronic rates. NOT IMPLEMENTED YET!
;
; PROGRAMMING NOTES
;
;    All data used by the program is stored in the structure STATE (defined
;    in burly_ups_widget). The tags of STATE are
;
;    .DATA   Contains the CHIANTI data, including the upsilons, energy
;            levels and final spline fits. The tags of DATA include
;            .STRUPS  contents of .upsdat file
;            .SPLUPS  contents of .splups file
;            .ESTR    energy levels
;
;    .TRANS  This contains the information for the particular transition
;            that is being fitted. The fields of .TRANS are filled by the
;            routine NEXT_TRANS.
;    --
;    OUTPUT FILES
;    The spline values are only written to the .splups file when the 'EXIT'
;    button is clicked. However, after each transition is fitted the set
;    of spline fits are saved to the file .splups_bck1. The previous
;    contents of the _bck1 file are sent to .splups_bck2. Thus if
;    new_burly_ups crashes, there should always be a backup of the latest
;    spline fits in one of these files.
;    --
;    METASTABLES
;    Metastable levels are identified through the external routine
;    METASTABLE_LEVELS. This routine has an in-built cut-off for A-values
;    to determine which levels are metastable. This cutoff can be varied
;    (keyword CUTOFF). Thus if a metastable level is not being picked up
;    new_burly_ups, then the user should vary this parameter.
;    --
;    EXTRAPOLATION
;    Extrapolation is performed within the subroutine spline_values, and
;    not in scale_ups as was done in the original burly_ups. Thus the
;    extrapolations are done 'on-the-fly' during the fitting rather than
;    being added to the scaled temperature and upsilon arrays that are
;    passed around in state.trans.
;
; CALLS
;
;    METASTABLE_LEVELS, BURLY_FIT_FN, DESCALE_ALL, READ_WGFA_STR,
;    READ_UPSDAT_STR, READ_ELVLC, READ_SPLUPS
;
; HISTORY
;
;    Beta 1, Jun-2005, Peter Young
;    Beta 2, 3-Jun-2006, Peter Young
;       added widget showing number of transitions that have been fitted.
;    Beta 3, 7-Jun-2006, Peter Young
;       added C over-ride button
;    Beta 4, 13-Jun-2006, Peter Young
;       done some tidying up of the widget
;    Beta 5, 11-Jul-2006, Peter Young
;       I've experimented with 15 point splines, but decided to scrap it.
;       However I've left some my code in, so this can be used
;       in future.
;    Beta 6, 2-Feb-2009, Peter Young
;       I've modified the extrapolation widget to allow the
;       user to manually specify the extrapolation value at 0.
;    Beta 7, 4-Feb-2009, Peter Young
;       Made cosmetic changes to plot displays.
;    Beta 8, 23-Mar-2009, Peter Young
;       Sometimes metastable_levels does not guess the metastable
;       levels will so I've added the keyword META= so that the
;       metastables can be directly specified by the user.
;    Beta 0, 3-Jan-2012, Ken Dere
;        Provides for the number of spline nodes between 5 and 9
;        renamed kpd_burly_ups for the time being.
;-
;
; --------------------------------------------
;
PRO number_transitions, data, show_str
;
; A widget in the top-left shows the number of transitions that have been
; fitted for each metastable. This routines updates the text for this widget
; based on the information in DATA.
;
i=where(data.estr.meta EQ 1,n_i)
metastr=data.estr.id[i]
;
show_str=strarr(n_i)
FOR j=0,n_i-1 DO BEGIN
  lev_n=i[j]+1
  IF n_tags(data.splstr) EQ 0 THEN BEGIN
    show_str[j]=metastr[j]+'   0'
  ENDIF ELSE BEGIN
    k=where(data.splstr.lvl1 EQ lev_n OR data.splstr.lvl2 EQ lev_n,nk)
    show_str[j]=metastr[j]+' '+strpad(trim(nk),3)
  ENDELSE
 ;
  k=where(data.strups.lvl1 EQ lev_n OR data.strups.lvl2 EQ lev_n,nk)
  show_str[j]=show_str[j]+' '+strpad(trim(nk),3)
ENDFOR

END

;-----------
PRO save_spl, state
;+
; Updates the STATE.DATA.SPLSTR structure by adding the transition that
; has just been fitted. Note that the ADD_TAG and REM_TAG routines are used
; to remove old data tags and add the new tags.
;-

splstr=state.data.splstr
trans=state.trans

str={lvl1: trans.l1, lvl2: trans.l2, t_type: trans.t_type, $
     gf: trans.gf, de: trans.de, c_ups: trans.c_ups, $
     nspl: trans.nspl, spl: trans.spl}

IF trans.l1 EQ 0 THEN print,'*** SAVING ZERO VALUES !! **'


IF n_tags(splstr) EQ 0 THEN BEGIN
  splstr=str
ENDIF ELSE BEGIN
  help,splstr.spl,str.spl
  splstr=[splstr,str]
ENDELSE

data=state.data

data=rem_tag(data,'splstr')
data=add_tag(data,splstr,'splstr')

state=rem_tag(state,'data')
state=add_tag(state,data,'data')

widget_control,state.burly_base,set_uval=state

END




;------------------
FUNCTION lin_extrap, x, y, x0
;+
; Fits a straight line through the points X,Y, and returns the fit value
; at X0.
;-

a=linfit(x,y)
return,a[0]+a[1]*x0

END


;---------------
PRO reset_extrap, state, allowed=allowed
;+
; When the next transition is loaded in, the extrapolation widget needs to
; be set back to default values.
;-

state.trans.extrapstr.lo=0.
state.trans.extrapstr.loset=0
IF keyword_set(allowed) THEN BEGIN
  state.trans.extrapstr.hi=0.
  state.trans.extrapstr.hiset=0
ENDIF ELSE BEGIN
  state.trans.extrapstr.hi=0.
  state.trans.extrapstr.hiset=0
ENDELSE

widget_control,state.e0_user_but,set_val=state.trans.extrapstr.loset
widget_control,state.e0_auto_but,set_val=state.trans.extrapstr.loset
widget_control,state.e1_but,set_val=state.trans.extrapstr.hiset
widget_control,state.e0_user_val,set_val=trim(state.trans.extrapstr.lo)
widget_control,state.e0_user_val,sens=0
widget_control,state.e0_auto_val,set_val=trim(state.trans.extrapstr.lo)
widget_control,state.e1_val,set_val=trim(state.trans.extrapstr.hi)
widget_control,state.burly_base,set_uval=state


END


;----------------
PRO spline_values, state, verbose=verbose, plot=plot
;+
; Works out the spline values using the C-parameter and scaled upsilon
; values. Currently applies the same method for 5 and 9 point splines -
; this is different from burly_ups_pry2, where I retained Ken's method
; for 5 points and used my method for 9 points. This uses only my method.
;-

st=state.trans.st
sups=state.trans.sups
mask=state.trans.mask
nspl=state.trans.nspl
extrap=state.trans.extrapstr

;
; Use MASK to restrict temperatures that are being considered
;
i=where(mask EQ 1)
st=st[i]
sups=sups[i]
;  have correct # of values here - kpd
; print, ' in spline_values, st = ',st
; print,' in spline_values, sups = ',sups

CASE extrap.loset OF
  1: BEGIN
    x0=0.
    y0=lin_extrap(st[0:1],sups[0:1],x0)
    st=[x0,st]
    sups=[y0,sups]
    widget_control,state.e0_auto_val,set_val=trim(y0)
  END
 ;
  2: BEGIN
    x0=0.
    y0=extrap.lo
    st=[x0,st]
    sups=[y0,sups]
  END
  ELSE:
ENDCASE

IF extrap.hiset EQ 1 THEN BEGIN
  x1=1.
  n=n_elements(st)
  y1=lin_extrap(st[n-2:n-1],sups[n-2:n-1],x1)
  st=[st,x1]
  sups=[sups,y1]
  widget_control,state.e1_val,set_val=trim(y1)
ENDIF ELSE BEGIN
  IF state.trans.hi NE 0. THEN BEGIN
    st=[st,1.0]
    sups=[sups,state.trans.hi]
  ENDIF
ENDELSE

; print, ' in spline_values, nspl = ',nspl
; print, '  in spline_values, st = ',st
; print,'  in spline_values, sups = ',sups

; IF nspl EQ 15 THEN BEGIN
;   xs1=dindgen(5)*0.05
;   xs2=dindgen(5)*0.125 + 0.25
;   xs=[xs1,xs2,xs1+0.8]
; ENDIF ELSE BEGIN
  xs=dindgen(nspl)/(nspl-1)
; ENDELSE
; print,'  in spline_values, xs = ',n_elements(xs),xs

y2=spl_init(st,sups)
ys=spl_interp(st,sups,y2,xs)
;
; curvefit dosen't work if the yy values are too small, so I
; need to scale up low numbers
;
getmin=min(abs(sups(where(sups NE 0d0))))
IF getmin LT 1.d-7 THEN scale=getmin*1d3 ELSE scale=1d0
;
; get initial Chi-square and print to screen
;
y2=spl_init(xs,ys)
init_sups=spl_interp(xs,ys,y2,st)
IF keyword_set(verbose) THEN BEGIN
  print,''
  print,format='("Initial chisq: ",e10.3)', $
       total(((init_sups-sups)/sups)^2)/ $
       (n_elements(sups)-n_elements(ys))
ENDIF
yy=sups/scale
weight=1d0/yy^2*100d0
aa=ys/scale
fit=curvefit(st,yy,weight,aa,sigmaa, $
             funct='burly_fit_fn', chisq=chisq,/noderiv, $
             iter=iter,tol=1d-10)
IF keyword_set(verbose) THEN print,format='("Final chisq: ",e10.3)',chisq
IF keyword_set(verbose) THEN print,format='("No. iterations: ",i3)',iter
;
nfree=n_elements(yy)-n_elements(aa)
IF keyword_set(verbose) THEN print,format='("Check: ",e10.3)',total(abs(yy))/1d7/nfree
;
; The YS are the spline values, so store them in state.trans.
; Make sure to set previous spline values to all -1's to prevent problems
; when changing from 9 points to 5 points.
;
;
;  only 9 values here - kpd
; print, ' in spline_values [344],  spl  = ',state.trans.spl
ys=aa*scale
state.trans.spl=-1.
; print, ' in spline_values [346], nspl = ',state.trans.nspl
; print, ' in spline_values [346],  ys  = ',ys
; print, ' in spline_values [346],  spl  = ',state.trans.spl
; help,/str,state.trans.spl
;  kpd
state.trans.spl[0:state.trans.nspl-1] = ys
; CASE state.trans.nspl OF
;   5: state.trans.spl[0:4]=ys
;   7: state.trans.spl[0:state.trans.nspl-1] = ys
;   9: state.trans.spl[0:8]=ys
;   11: state.trans.spl[0:10] = ys
;   15: state.trans.spl=ys
; ENDCASE

; print, ' in spline_values [353], spl = ',state.trans.spl
state.trans.spl[0:state.trans.nspl-1] = ys
;
; print warning if spline values fall below 0
;
i=where(ys LT 0.)
IF i[0] NE -1 AND keyword_set(verbose) THEN BEGIN
  istr=trim(i[0])
  IF n_elements(i) GT 1 THEN BEGIN
    FOR j=1,n_elements(i)-1 DO istr=istr+','+trim(i[j])
  ENDIF
  txt='**WARNING: Spline value ['+istr+'] is less than zero'
  widget_control,state.mess_text,set_val=txt
ENDIF ELSE BEGIN
  widget_control,state.mess_text,set_val=''
ENDELSE
;
widget_control,state.burly_base,set_uval=state
;
xplot=findgen(41)*0.025
y2=spl_init(xs,ys)
yplot=spl_interp(xs,ys,y2,xplot)
sups_fit=spl_interp(xs,ys,y2,st)


; Now use the spline values (stored in trans.spl) to derive upsilons over
; the temperature range. These values are stored in trans.upsfit.
; Note the temperature range is the original one and not that selected by
; MASK.
; TRANS has been set up so that it has exactly the form required by the
; standard CHIANTI routine descale_all. As there's only 1 transition in
; TRANS, then I select index 0.
;
descale_all,state.trans.temp,state.trans,0,upsfit
state.trans.upsfit=upsfit
widget_control,state.burly_base,set_uval=state

END


;------------
PRO scale_ups, state

trans=state.trans

;
; Note: when transition data is first loaded into burly_ups, the
; subroutine select_data sets sups0 and sups1 (the scaled upsilons
; values at 0 and 1) to be -1. In the code below, checks are made to
; see if sups0 and sups1 are greater than 0 in order that these points
; aren't included.
;
; For proton rates I'm allowing for negative extrapolations at 0. For
; this reason, I check in this case for sups0 to be equal to -1
;

extrapstr=state.trans.extrapstr
;
trans=state.trans
de_in=trans.de
c_ups_curr=trans.c_ups
;
kt=1.38062d-16*trans.temp/2.1799d-11   ; in Ryd
upsilon=trans.ups

; print, ' in scale_ups, de_in = ', de_in
; print, ' in scale_ups, c_ups = ',c_ups_curr
; print, ' in scale_ups, kt = ',kt
; print, ' in scale_ups, upsilon = ',upsilon
;
CASE trans.t_type OF

1: begin
  xt=kt/de_in
  st=1.-alog(c_ups_curr)/alog(xt+c_ups_curr)
  sups=upsilon/alog(xt+exp(1.))
;   IF extrapstr.loset EQ 1 THEN BEGIN
;     st=[0.,st,1.]
;     sups=upsilon/alog(xt+exp(1.))
;     sups=[extrapstr.lo,sups,trans.hi]
;   ENDIF ELSE BEGIN
;     st=[st,1.]
;     sups=upsilon/alog(xt+exp(1.))
;     sups=[sups,trans.hi]
;   ENDELSE
END
;
2: BEGIN
  xt=kt/de_in
  st=xt/(xt+c_ups_curr)
  sups=upsilon
;   IF extrapstr.hiset EQ 1 THEN BEGIN
;     st=[st,1.]
;     sups=[sups,extrapstr.hi]
;   endif
;   IF extrapstr.loset EQ 1 THEN BEGIN
;     st=[0.,st]
;     sups=[extrapstr.lo,sups]
;   ENDIF
END
;
3: BEGIN
  xt=kt/de_in
  st=xt/(xt+c_ups_curr)
  sups=(xt+1.)*upsilon
;   IF extrapstr.hiset EQ 1 THEN BEGIN
;     st=[st,1.]
;     sups=[sups,extrapstr.hi]
;   endif
;   IF extrapstr.loset EQ 1 THEN BEGIN
;     st=[0.,st]
;     sups=[extrapstr.lo,sups]
;   endif
end
;
4:  begin
  xt=kt/de_in
  st=1.-alog(c_ups_curr)/alog(xt+c_ups_curr)
  sups=upsilon/alog(xt+c_ups_curr)
;   IF extrapstr.loset EQ 1 THEN BEGIN
;     st=[0.,st]
;     sups=[extrapstr.lo,sups]
;   ENDIF
;   st=[st,1.]
;   sups=[sups,trans.hi]
end
;
5: BEGIN
  xt=kt/de_in
  st=xt/(xt+c_ups_curr)
  sups=upsilon*kte
END
;
6:  BEGIN
  xt=kt/de_in
  st=xt/(xt+c_ups_curr)
  sups=alog10(upsilon)
;   IF extrapstr.hiset EQ 1 THEN BEGIN
;     st=[st,1.]
;     sups=[sups,extrapstr.hi]
;   endif
;   IF extrapstr.loset EQ 1 THEN BEGIN
;     st=[0.,st]
;     sups=[extrapstr.lo,sups]
;   endif
END
;
ELSE:  print,' t_type_ups ne 1,2,3,4=',t_type_ups
;
ENDCASE

; print, ' in scale_ups, xt = ',xt
; print, ' in scale_ups, st = ',st
; print, ' in scale_ups, sups = ',sups

state.trans.st=st
; print, ' in scale_ups, state.trans.sups = ',state.trans.sups
state.trans.sups=sups
widget_control,state.burly_base,set_uval=state

END


;------------------
FUNCTION find_mid_c, trans, st=st
;+
; Works out the scaling parameter (C) that puts the middle temperature
; value at a scaled temperature value of 0.5.
;-

IF n_elements(st) EQ 0 THEN st=0.5

t=trans.temp
mask=trans.mask
i=where(mask EQ 1)
t=t[i]

n=n_elements(t)/2
t=t[n]
kt_e=t*8.61735d-5/(trans.de*13.606)

CASE 1 OF
  trans.t_type EQ 1 OR trans.t_type EQ 4: BEGIN
    IF st EQ 0.5 THEN BEGIN
      return,sqrt(0.25+kt_e)+0.5
    ENDIF ELSE BEGIN
     ;
     ; The lower bound in c1 is 1.0, the upper bound is (1+kt_e)
     ;
      c1=findgen(101)/100.*alog(kt_e+1)*(1.-st)/st
      c1[0]=1e-8
      c1=exp(c1)
      getmin=min(abs(c1^(1./(1.-st))-c1-kt_e),i)
      print,'find_mid_c: ',st,c1[i]
      return,c1[i]
;       IF kt_e GE 5. THEN BEGIN
;         c1=(findgen(101)/200.-0.25)*kt_e^(1./st) + kt_e^(1./st)
;         getmin=min(abs(c1^(1./st)-c1-kt_e),i)
;         return,c1[i]
;       ENDIF ELSE BEGIN
;         c1=exp(findgen(101)*0.01609+0.001)
;         getmin=min(abs(c1^(1./st)-c1-kt_e),i)
;         return,c1[i]
;       ENDELSE
    ENDELSE
  END
 ;
  ELSE: return,kt_e*(1./st -1.)
ENDCASE

END


;----------------
PRO do_transition, state
;+
; Contains the commands to perform the minimization for C, do the spline
; fit and then plot the results.
;
; The range over which C is varied is determined by finding the C values
; which place the central temperature point at scaled temperatures of 0.10
; and 0.75. The 0.10 value is used because for some allowed transitions
; it is necessary to shove the points over to small scaled values in order
; to get a decent fit.
;-

plot_ups,state
low_st=0.10
hi_st = 0.8  ;  was = 0.75 (kpd - but does not make much difference)
c1=find_mid_c(state.trans,st=low_st)   ; get bounds for C
c2=find_mid_c(state.trans,st=hi_st)   ;
max_perc=100.
c_all=0.
FOR i=0,100 DO BEGIN
  diff=abs(c2-c1)/100.
  state.trans.c_ups=diff*i+min([c1,c2])
  c_all=[c_all,state.trans.c_ups]
  widget_control,state.burly_base,set_uval=state
  scale_ups,state
  spline_values,state
 ;
  mask=state.trans.mask
  j=where(mask EQ 1)
  ups=state.trans.ups[j]
  upsfit=state.trans.upsfit[j]
 ;
  max_perc=[max_perc,max(abs((ups-upsfit)/ups))]
ENDFOR
max_perc=max_perc[1:*]
c_all=c_all[1:*]
print,'Min max perc: ',min(max_perc)
yra=[0,min(max_perc)*1000.]
;
wset,state.plot_dat.c_id
plot,c_all,max_perc*100.,xmarg=[6,1],ymarg=[3,2],yra=yra,/xsty, $
     tit='C minimization', $
     xtit='C value', $
     ytit='Max % error in upsilon fit'
;
getmin=min(max_perc,i)
;
print,'C-range: ',c1,c2
oplot,c_all[i]*[1,1],[-1d6,1d6],line=1
;
state.trans.c_ups=c_all[i]
widget_control,state.burly_base,set_uval=state
scale_ups,state
spline_values,state,/verbose
plot_splups,state,/spline
plot_ups,state

END


;--------------
PRO plot_splups, state, spline=spline
;+
; Plot the scaled upsilons in the bottom right window, and optionally
; (with /SPLINE) overplot the spline fit.
;-

wset,state.plot_dat.spl_id

cstr=string(state.trans.c_ups,format='(f15.2)')
tstr=trim(state.trans.t_type)

st=state.trans.st
sups=state.trans.sups
mask=state.trans.mask

i=where(mask EQ 1)
st=st[i]
sups=sups[i]


yplot=-1
IF keyword_set(spline) THEN BEGIN
  ys=state.trans.spl
  n=state.trans.nspl
  ys=ys[0:n-1]       ; picks out 5 points if set
  IF n EQ 15 THEN BEGIN
    xs1=dindgen(5)*0.05
    xs2=dindgen(5)*0.125 + 0.25
    xs=[xs1,xs2,xs1+0.8]
  ENDIF ELSE BEGIN
    xs=dindgen(n)/(n-1)
  ENDELSE
;  xs=dindgen(n)/(n-1)
  xplot=findgen(41)*0.025
  y2=spl_init(xs,ys)
  yplot=spl_interp(xs,ys,y2,xplot)
ENDIF

yrange=[0,max([sups,state.trans.hi,yplot])]

plot,st,sups,psym=5,xra=[0,1], $
     charsiz=1.5, $
     xtit='Scaled Temperature', $
     ytit='Scaled Upsilon', $
     tit='C='+trim(cstr)+', Type='+tstr, $
     yticklen=-0.015, yrange=yrange

oplot,xplot,yplot


;
; Draw high-T limit as filled circle
;
i=findgen(32)/31.*2.*!pi
x0=2.*cos(i)
y0=2.*sin(i)
usersym,x0,y0,/fill
;
IF state.trans.hi NE 0. AND state.trans.extrapstr.hiset EQ 0 THEN BEGIN
  plots,1.0,state.trans.hi,psym=8
ENDIF


END


;-----------
PRO plot_ups, state
;
; Plots upsilons in bottom right window (triangles). The derived upsilons
; are overplotted (stars) after the fit has been performed, and an additional
; dotted line is plotted showing the error in the fit for each point.
;
wset,state.plot_dat.ups_id

t=state.trans.temp
ups=state.trans.ups
mask=state.trans.mask

i=where(mask EQ 1)
t=t[i]
ups=ups[i]

title='Transition: '+trim(state.trans.l1)+'-'+trim(state.trans.l2)
plot,t,ups,psym=5,/xlog, $
     charsiz=1.5, $
     xtit='Temperature / K', $
     ytit='Upsilon', $
     title=title

ylims=!y.crange
xlims=!x.crange
ypos=ylims[1]*0.10+ylims[0]*0.9
xpos=xlims[1]*0.05+xlims[0]*0.95
;xyouts,10.^xpos,ypos,'Transition: '+trim(state.trans.l1)+'-'+trim(state.trans.l2),charsiz=1.2


IF total(state.trans.upsfit) NE 0. THEN BEGIN
  upsfit=state.trans.upsfit
  upsfit=upsfit[i]
  oplot,t,upsfit,psym=2
 ;
  grad= ( ylims(1)-ylims(0) ) / 0.20
  coeff= ( ylims(1) + ylims(0) ) * 0.50
  oplot,[min(t),max(t)],[1.,1.]*coeff,line=2   ; dashed line
  oplot,t,(upsfit/ups-1.0)*grad+coeff,line=1
 ;
  max_perc_diff=max(abs(upsfit/ups-1.0))*100.
  ypos=ylims[1]*0.06+ylims[0]*0.94
  xyouts,10.^xpos,ypos,'Max diff: '+trim(string(format='(f10.1)',max_perc_diff))+'%',charsiz=1.2
ENDIF


END


;---------------
PRO write_splupsk, fname, data
;+
; Writes out a .splups file
;-
n=n_elements(data.splstr)
splstr=data.splstr
iz=data.iz
ion=data.ion

form1='(5i3,8e10.3)'
form2='(5i3,12e10.3)'
form2='(5i3,18e10.3)'

openw,lout,fname,/get_lun
FOR i=0,n-1 DO BEGIN
;     print, ' i, nspl = ',i,splstr[i].nspl
    form1 = '(5i3,' + string(splstr[i].nspl+3) + 'e10.3)'
    printf, lout, format=form1,iz,ion, $
               splstr[i].lvl1,splstr[i].lvl2,splstr[i].t_type, $
               splstr[i].gf,splstr[i].de,splstr[i].c_ups, $
               splstr[i].spl[0:splstr[i].nspl-1]
;       CASE splstr[i].nspl OF
;     5: BEGIN
;       printf,lout,format=form1,iz,ion, $
;            splstr[i].lvl1,splstr[i].lvl2,splstr[i].t_type, $
;            splstr[i].gf,splstr[i].de,splstr[i].c_ups, $
;            splstr[i].spl[0:4]
;     END
;     9: BEGIN
;       printf,lout,format=form2,iz,ion, $
;            splstr[i].lvl1,splstr[i].lvl2,splstr[i].t_type, $
;            splstr[i].gf,splstr[i].de,splstr[i].c_ups, $
;            splstr[i].spl
;     END
;     15: BEGIN
;       printf,lout,format=form3,iz,ion, $
;            splstr[i].lvl1,splstr[i].lvl2,splstr[i].t_type, $
;            splstr[i].gf,splstr[i].de,splstr[i].c_ups, $
;            splstr[i].spl
;     END
;   ENDCASE

ENDFOR

printf,lout,' -1'
n=n_elements(data.upsref)
FOR i=0,n-1 DO printf,lout,data.upsref[i]
printf,lout,' -1'

free_lun,lout

END


;---------------
PRO write_splbck_k, data
;+
; Writes out a backup file (.splups_bck1) containing the latest splups data.
; The previous file is sent to .splups_bck2.
;-

file1=data.splupsname+'_bck1'
file2=data.splupsname+'_bck2'
chck1=file_search(file1)
IF chck1 NE '' THEN file_move,file1,file2,/overwrite
write_splupsk,file1,data

END


;-------------
PRO next_trans, state, exit=exit
;+
; This routine works out the next transition to be fitted. It starts from
; the lowest of the metastable levels and works upwards until it finds the
; first transition that hasn't been fitted.
;
; The routine also does the job of setting the structure TRANS that
; contains all the data for the transition to be fitted.
;
; EXIT is used to identify if all transitions have been fitted. (EXIT=1
; implies yes.)
;-

;
; reset the warning message text
;
widget_control,state.mess_text,set_val=''

exit=0

splstr=state.data.splstr
strups=state.data.strups
trans=state.trans
;
IF n_tags(splstr) EQ 0 THEN BEGIN
  l1=1
  i=where(strups.lvl1 EQ 1)
  l2=min(strups[i].lvl2)
ENDIF ELSE BEGIN
  i=where(state.data.estr.meta EQ 1)  ; positions of all metastables
  meta=i+1                            ; list of metastable levels
  nm=n_elements(meta)
  FOR i=0,nm-1 DO BEGIN
    chcks=where(splstr.lvl1 EQ meta[i] OR splstr.lvl2 EQ meta[i],ns)
    chcku=where(strups.lvl1 EQ meta[i] OR strups.lvl2 EQ meta[i],nu)
    print,'**** ',i,ns,nu
    IF ns LT nu THEN BEGIN
      FOR j=0,nu-1 DO BEGIN
        l1=strups[chcku[j]].lvl1
        l2=strups[chcku[j]].lvl2
        k=where( (splstr.lvl1 EQ l1 AND splstr.lvl2 EQ l2) OR $
                 (splstr.lvl1 EQ l2 AND splstr.lvl2 EQ l1))
        IF k[0] EQ -1 AND l1 NE l2 THEN GOTO,lbl1
      ENDFOR
    ENDIF
  ENDFOR
  txt=['All transitions have been fitted.', $
       'Please click EXIT to exit and write out the .splups file.']
  result=dialog_message(txt)
  exit=1
  return
ENDELSE

;
; Swap levels if l1 > l2
;
IF l1 GT l2 THEN BEGIN
  lsave=l1
  l1=l2
  l2=lsave
ENDIF


lbl1:

i=where(strups.lvl1 EQ l1 AND strups.lvl2 EQ l2)
;
hi=4.*strups[i].gf/strups[i].de
IF hi NE 0. THEN BEGIN
  IF strups[i].gf GE 1d-3 THEN type=1 ELSE type=4
ENDIF ELSE BEGIN
  type=2
ENDELSE
widget_control,state.type_buts,set_val=type-1
;
n=n_elements(strups[i].ups)
;
trans.l1=l1
trans.l2=l2
trans.gf=strups[i].gf
trans.de=strups[i].de
trans.temp=strups[i].temp
trans.ups=strups[i].ups
trans.t_type=type
trans.hi=hi
trans.upsfit=0.
trans.st=0.
trans.sups=0.
;
; only reset MASK if MASKSAVE is set
;
IF trans.masksave EQ 0 THEN BEGIN
  trans.mask=intarr(n)+1
  widget_control,state.lo_slider,set_val=0
  widget_control,state.hi_slider,set_val=n-1
ENDIF
;
trans.c_ups=0.
trans.spl=0.

;
; reset the C over-ride widgets
;
widget_control,state.c_over_but,set_val=0
widget_control,state.c_over_val,sensitive=0


state.trans=trans
IF hi NE 0. THEN allowed=1 ELSE allowed=0
reset_extrap,state,allowed=allowed
widget_control,state.burly_base,set_uval=state

;
; The following lines create the text (displayed in trans_text) that
; identify the transition being fitted.
;
zion2spectroscopic,state.data.iz,state.data.ion,name
;
CASE 1 OF
  trans.gf LT 0.01 AND trans.gf NE 0.: gf_text=trim(string(format='(e10.2)',trans.gf))
  trans.gf EQ 0.: gf_text='0.0'
  trans.gf GE 0.01: gf_text=trim(string(format='(f6.3)',trans.gf))
ENDCASE
gf_text='gf='+gf_text
;
de_text='DE='+trim(string(format='(f8.3)',trans.de))+' Ryd'
trans_text=name+' *** '+trim(state.data.estr.id[l1-1])+' - '+ $
     trim(state.data.estr.id[l2-1])+' *** '+gf_text+' *** '+de_text
widget_control,state.trans_text,set_val=trans_text

END


;*******************************
PRO burly_base_event, event
;*******************************

WIDGET_CONTROL,Event.top, get_uvalue=state

CASE 1 OF
  event.id EQ state.exit: BEGIN
    CASE event.value OF
      0: BEGIN                   ; --- RE-FIT TRANSITION
        widget_control,state.c_over_but,get_val=val
        IF val EQ 1 THEN BEGIN
          widget_control,state.c_over_val,get_val=val
          state.trans.c_ups=val
          widget_control,state.burly_base,set_uval=state
         ;
         ; perform spline fit and plot results
         ;
          scale_ups,state
          spline_values,state
          plot_splups,state,/spline
          plot_ups,state
        ENDIF ELSE BEGIN
          do_transition,state
        ENDELSE
      END
     ;
      1: BEGIN                   ; --- NEXT TRANSITION
        IF state.trans.l1 NE 0 THEN BEGIN
          save_spl,state
          write_splbck_k,state.data
        ENDIF
       ;
        number_transitions,state.data,show_str
        widget_control,state.status_label,set_val=show_str
       ;
        next_trans,state, exit=exit
        IF exit NE 1 THEN do_transition,state
      END
     ;
      2: BEGIN                   ; --- EXIT
       ;
       ; If the user exits without having fitted a transition, then don't
       ; want to save the spline values. Check if l1=0.
       ;
;        IF state.trans.l1 NE 0 THEN save_spl,state
        IF n_tags(state.data.splstr) NE 0 THEN  $
             write_splupsk,state.data.splupsname,state.data
        WIDGET_CONTROL, event.top, /DESTROY
      END
      ELSE:
    ENDCASE
  END
 ;
 ; Event handler for no. of spline points (5/9)
 ; --------------------------------------------
  event.id EQ state.npt_but: BEGIN
    widget_control,state.npt_but,get_uvalue=uval
    widget_control,state.npt_but,get_value=val
    nspl=uval[val]
    IF nspl NE state.trans.nspl THEN BEGIN
      state.trans.nspl=nspl
      widget_control,state.burly_base,set_uval=state
   ;  kpd added this may 22, 2012
   ; perform spline fit and plot results
   ;
		scale_ups,state
		spline_values,state
		plot_splups,state,/spline
		plot_ups,state
    ENDIF

  END
 ;
 ; Event handler for extrapolation widgets
 ; ---------------------------------------
  event.id EQ state.e0_auto_but: BEGIN
    widget_control,state.e0_auto_but,get_val=value
    state.trans.extrapstr.loset=value
    state.trans.extrapstr.lo=0
   ;
   ; Make sure the user button is switched off
   ;
    IF state.trans.extrapstr.loset EQ 1 THEN BEGIN
      widget_control,state.e0_user_but,set_val=0
      widget_control,state.e0_user_val,sens=0
    ENDIF
   ;
    widget_control,state.burly_base,set_uval=state
  END
 ;
  event.id EQ state.e0_user_but: BEGIN
    widget_control,state.e0_user_but,get_val=value
    IF value GT 0 THEN BEGIN
      state.trans.extrapstr.loset=value+1 ; note: set to 2 rather than 1
      widget_control,state.e0_auto_val,get_val=value
      widget_control,state.e0_user_val,sens=1,set_val=value
      state.trans.extrapstr.lo=value
     ;
     ; Make sure the auto button is switched off
     ;
      widget_control,state.e0_auto_but,set_val=0
    ENDIF ELSE BEGIN
      state.trans.extrapstr.loset=value
    ENDELSE
    widget_control,state.burly_base,set_uval=state
  END
 ;
  event.id EQ state.e1_but: BEGIN
    widget_control,state.e1_but,get_val=value
    state.trans.extrapstr.hiset=value
    state.trans.extrapstr.hi=0
    widget_control,state.burly_base,set_uval=state
  END
 ;
  event.id EQ state.e0_user_val: BEGIN
    widget_control,state.e0_user_val,get_val=value
    state.trans.extrapstr.lo=value
    widget_control,state.burly_base,set_uval=state
  END
 ;
  event.id EQ state.e1_val: BEGIN
    widget_control,state.e1_val,get_val=value
    state.trans.extrapstr.hiset=value
    state.trans.extrapstr.hi=0
    widget_control,state.burly_base,set_uval=state
  END
 ;
 ;
 ; Event handler for temperature range sliders
 ; -------------------------------------------
  event.id EQ state.lo_slider: BEGIN
    mask=state.trans.mask*0     ; set to all 0's
    widget_control,state.lo_slider,get_val=val1
    widget_control,state.hi_slider,get_val=val2
    nspl=state.trans.nspl
    IF val1 GT val2-nspl+1 THEN BEGIN
      val1=val2-nspl+1
      widget_control,state.lo_slider,set_val=val1
    ENDIF
    mask[val1:val2]=1
    print,'Lo slider: ',val1,val2
    state.trans.mask=mask
    widget_control,state.burly_base,set_uval=state
  END
 ;
  event.id EQ state.hi_slider: BEGIN
    mask=state.trans.mask*0     ; set to all 0's
    widget_control,state.lo_slider,get_val=val1
    widget_control,state.hi_slider,get_val=val2
    nspl=state.trans.nspl
    IF val2 LT val1+nspl-1 THEN BEGIN
      val2=val1+nspl-1
      widget_control,state.hi_slider,set_val=val2
    ENDIF
    mask[val1:val2]=1
    print,'Hi slider: ',val1,val2
    state.trans.mask=mask
    widget_control,state.burly_base,set_uval=state
  END
 ;
  event.id EQ state.but_slider: BEGIN
    widget_control,state.but_slider,get_val=val
    state.trans.masksave=val
    widget_control,state.burly_base,set_uval=state
  END
 ;
 ; Event handler for transition types
 ; ----------------------------------
  event.id EQ state.type_buts: BEGIN
    widget_control,state.type_buts,get_val=val
    state.trans.t_type=val+1
    widget_control,state.burly_base,set_uval=state
  END
 ;
 ; Event handler for C over-ride button
 ; ------------------------------------
  event.id EQ state.c_over_but: BEGIN
    widget_control,state.c_over_but,get_val=val
    state.trans.c_over_val=val
    widget_control,state.c_over_val,sensitive=val
   ;
   ; Set c_over_val to be current value of c_ups
   ;
    c_ups_str=trim(string(format='(f10.3)',state.trans.c_ups))
    IF val EQ 1 THEN  widget_control,state.c_over_val, $
         set_val=c_ups_str
    widget_control,state.burly_base,set_uval=state
  END
 ;
 ; Event handler for C over-ride value
 ; -----------------------------------
  event.id EQ state.c_over_val: BEGIN
    widget_control,state.c_over_val,get_val=val
    state.trans.c_ups=val
    widget_control,state.burly_base,set_uval=state
   ;
   ; perform spline fit and plot results
   ;
    scale_ups,state
    spline_values,state
    plot_splups,state,/spline
    plot_ups,state
  END
;
ENDCASE

END


;******************************************
PRO burly_ups_widget, group=group, data=data, trans=trans
;******************************************
;+
; Creates the widget (base name BURLY_BASE).
; All data is stored in the structure STATE which is set to the UVALUE
; of BURLY_BASE.
;-


; Set main base
;
burly_base=widget_base(/col,map=1,title='BURLY_UPS')

chianti_font,font
chianti_font,bigfont,/big
chianti_font,fixfont,/fix

base_top=widget_base(burly_base,/row,map=1)
base_left=widget_base(base_top,/col,map=1,/frame)
base_right=widget_base(base_top,/col,map=1,/frame)

exit=cw_bgroup(base_left, $
               ['RE-FIT TRANSITION','NEXT TRANSITION','EXIT'], $
               /row,font=bigfont)

sub_base1=widget_base(base_left,/row,map=1)

;
; Text widget for displaying no. of transitions that have been fitted.
;
number_transitions,data,show_str
;
meta_base=widget_base(sub_base1,/col,map=1,/frame)
status_label=widget_text(meta_base,val=show_str,$
                         font=fixfont,/align_left, $
                         ysiz=n_elements(show_str))


;
; Sliders for temperature range
;
n=n_elements(data.strups[0].temp)
sub_base3=widget_base(sub_base1,/col,map=1,/frame)
lohi_lbl=widget_label(sub_base3,val='Temperature range',font=font)
lo_slider=widget_slider(sub_base3,min=0,max=n-1,value=0)
hi_slider=widget_slider(sub_base3,min=0,max=n-1,value=n-1)
but_lbl=widget_label(sub_base3,val='Save range?',font=font)
but_slider=cw_bgroup(sub_base3,['No','Yes'],font=font, $
                     set_val=trans.masksave,/row,/exclusive)


;
; Buttons for selecting transition type. Note that the number of transition
; types is hardcoded here.
;
sub_base4=widget_base(sub_base1,/col,map=1)
type_buts=cw_bgroup(sub_base4,['Type 1 (allowed)', $
                               'Type 2 (forbidden)', $
                               'Type 3 (intercombination)', $
                               'Type 4 (allowed, small gf)', $
                               'Type 5 (dielectronic)', $
                               'Type 6 (log forbidden)'], $
                    set_val=0,/col,/exclusive,font=font,/frame)


;
; Widgets for extrapolation
;
extrap_base=widget_base(sub_base1,/col,map=1,/frame)
extrap_label=widget_label(extrap_base,val='EXTRAPOLATION',$
                          font=font,/align_left)
;
extrap0_base=widget_base(extrap_base,/col,map=1)
;
extrap0_base_auto=widget_base(extrap0_base,/row,map=1)
e0_auto_but=cw_bgroup(extrap0_base_auto,['0 (auto)'],font=font,/nonexclusive)
e0_auto_label=widget_label(extrap0_base_auto,val='Value:',font=font)
e0_auto_val=widget_text(extrap0_base_auto,value='',font=font,xsiz=10)
;
extrap0_base_user=widget_base(extrap0_base,/row,map=1)
e0_user_but=cw_bgroup(extrap0_base_user,['0 (user)'],font=font,/nonexclusive)
e0_user_label=widget_label(extrap0_base_user,val='Value:',font=font)
e0_user_val=widget_text(extrap0_base_user,value='',/editable,font=font,xsiz=10,sens=0)

;
extrap1_base=widget_base(extrap_base,/row,map=1)
e1_but=cw_bgroup(extrap1_base,['1'],font=font,/nonexclusive)
e1_label=widget_label(extrap1_base,val='Value:',font=font)
e1_val=widget_text(extrap1_base,value='',/editable,font=font,xsiz=10)


;
; 5/9 point splines
;
;
;sub_base2=widget_base(sub_base1,/col,map=1)
uval=[5,6,7,8,9]
i=where(trans.nspl EQ uval)
if i LE 0 then begin
    if trans.nspl lt min(uval) then begin
        i = 0
    endif else i = n_elements(uval) - 1
endif
npt_but=cw_bgroup(base_right,['5 points','6 points','7 points','8 points','9 points'],font=font,$
                  /exclusive,/row,set_val=i,uval=uval,/frame)
;
c_over_base=widget_base(base_right,/row,map=1)
c_over_but=cw_bgroup(c_over_base,['C over-ride'],font=font, $
                     /nonexclusive)
c_over_val=widget_text(c_over_base,value='',/editable,font=font, $
                       xsiz=10,sens=0)
;
c_plot=widget_draw(base_right,retain=1,uval='c_plot',$
                   xsiz=250,ysiz=200)

txt_base1=widget_base(base_left,/row,map=1)
trans_label=widget_label(txt_base1,val='TRANSITION:',font=bigfont)
trans_text=widget_text(txt_base1,val='',font=bigfont,xsiz=60)
;
txt_base2=widget_base(base_left,/row,map=1)
mess_label=widget_label(txt_base2,val='MESSAGES:',font=font)
mess_text=widget_text(txt_base2,val='',font=font,xsiz=50)

;
; Plot widgets
;
plot_base=widget_base(burly_base,/row,map=1)
;
ups_plot= WIDGET_DRAW(plot_base, $
                       RETAIN=1, $
                       uval='ups_plot', $
                       XSIZE=500, $
                       YSIZE=400)
;
spl_plot=WIDGET_DRAW(plot_base, $
                       RETAIN=1, $
                       uval='splups_plot', $
                       XSIZE=500, $
                       YSIZE=400)
;
plot_dat={ups_id: 0, spl_id: 0, c_id: 0}


;
; Set up structure that contains all information used by the program.
;
state={trans: trans, data: data, exit: exit, plot_dat: plot_dat, $
       ups_plot: ups_plot, spl_plot: spl_plot, c_plot: c_plot, $
       burly_base: burly_base, npt_but: npt_but, $
       e0_auto_but: e0_auto_but, e0_user_but: e0_user_but, $
       e1_but: e1_but, $
       e0_auto_val: e0_auto_val, e0_user_val: e0_user_val, $
       e1_val: e1_val, $
       trans_text: trans_text, $
       lo_slider: lo_slider, hi_slider: hi_slider, but_slider: but_slider, $
       type_buts: type_buts, mess_text: mess_text, $
       status_label: status_label, $
       c_over_but: c_over_but, c_over_val: c_over_val}


WIDGET_CONTROL, burly_base, /REALIZE, set_uvalue=state

WIDGET_CONTROL, ups_plot, GET_VALUE=ups_id
WIDGET_CONTROL, spl_plot, GET_VALUE=spl_id
WIDGET_CONTROL, c_plot, GET_VALUE=c_id

state.plot_dat.ups_id=ups_id
state.plot_dat.spl_id=spl_id
state.plot_dat.c_id=c_id

WIDGET_CONTROL, burly_base, /REALIZE, set_uvalue=state


XMANAGER, 'burly_base', burly_base, group=group

END


;----------------
PRO kpd_burly_ups, directory,iz,ionn,proton=proton, diel=diel, cutoff=cutoff, $
                   meta=meta

IF n_params() LT 3 THEN BEGIN
  print,'Use  IDL> new_burly_ups, directory, iz, ion, proton=, diel='
  return
ENDIF

zion2name,iz,ionn,name
root=concat_dir(directory,name)
fname=root+'.upsdat'
read_upsdat_str,fname,strups,upsref

IF keyword_set(proton) THEN BEGIN
  splupsname=root+'.psplups'
ENDIF ELSE BEGIN
  splupsname=root+'.splups'
ENDELSE
;
chck=file_search(splupsname)
IF chck NE '' THEN BEGIN
  print,'%BURLY_UPS: a .splups file already exists. Fitting will continue from the last'
  print,'            entry in this file.'
  read_splups,splupsname,splstr,splref,prot=proton
ENDIF ELSE BEGIN
  splstr=0
ENDELSE

wgfaname=concat_dir(directory,name+'.wgfa')
read_wgfa_str,wgfaname,wgfastr,ref

IF n_elements(cutoff) EQ 0 THEN cutoff=1e5
IF n_elements(meta) EQ 0 THEN metastable_levels,name,meta,cutoff=cutoff,path=directory


elvlcname=concat_dir(directory,name+'.elvlc')
read_elvlc,elvlcname,l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,ref
n1=n_elements(l1)
n2=n_elements(meta)
IF n2 GT n1 THEN meta=meta[0:n1-1]
IF n2 LT n1 THEN BEGIN
  print,'%KPD_BURLY_UPS: nlevels, nmeta = ',n1,n2
  print,'%KPD_BURLY_UPS: too many levels in .wgfa file. kpd says to continue anyway'
;  print,'%KPD_BURLY_UPS: too many levels in .wgfa file. Returning...'
;  return
ENDIF
estr={i: l1, id: term, meta: meta}

data={splstr:splstr, strups: strups, estr: estr, splupsname: splupsname, $
     iz: iz, ion: ionn, upsref: upsref, cutoff: cutoff}

;
; Set up the TRANS structure that stores information on the transition
; to be fitted. Note the initial type of spline fit (5 or 9 points) is
; set here
;
;  kpd size of nspl set here???
n=n_elements(strups[0].ups)
IF n LT 9 THEN nspl=5 ELSE nspl=11
;
; The extrapolation structure works as follows:  .loset and .hiset
; indicate whether extrapolation to low and high temperatures will be
; performed (0-no; 1-yes).  The tags .lo and .hi give the values of
; the extrapolated points, but only if they've been typed in
; manually by the user in the text boxes
extrapstr={lo: 0., loset: 0, hi: 0., hiset: 0}
;
;  !!!  nspl set here - kpd
;
trans={l1: 0, l2: 0, gf: 0., de: 0., $
       temp: fltarr(n), ups: fltarr(n), $
       t_type: 0, $
       upsfit:fltarr(n), $
       st: fltarr(n), sups: fltarr(n), mask: intarr(n)+1, $
       masksave: 0, $
       c_ups: 0., hi: 0., $
       extrapstr: extrapstr, nspl: nspl, spl: fltarr(11), $
       c_over_but: 0, c_over_val: 0.0}

burly_ups_widget, data=data, trans=trans

END
