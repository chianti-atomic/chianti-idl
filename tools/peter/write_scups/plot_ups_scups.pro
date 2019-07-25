
FUNCTION plot_ups_scups, input, ups_ind=ups_ind, transition=transition

;+
; NAME:
;     PLOT_UPS_SCUPS
;
; PURPOSE:
;     Allows the scaled upsilons for a specific transition to be
;     directly plotted to the screen.
;
; CATEGORY:
;     CHIANTI; data assessment, plotting.
;
; CALLING SEQUENCE:
;     Result = PLOT_UPS_SCUPS( Input )
;
; INPUTS:
;     Input:   Either the name of a CHIANTI .ups file, or a structure
;              returned by read_ups.pro. 
;
; OPTIONAL INPUTS:
;     Ups_Ind: The index of a transition within UPSSTR.DATA.
;     Transition:  A 2-element array specifying the lower and upper
;                  levels of a transition. For example, [1,20].
;
; OUTPUTS:
;     Returns a plot object containing a plot of scaled temperature
;     vs. scaled upsilons for the specified transition.
;
; MODIFICATION HISTORY:
;     Beta 1, 8-Feb-2017, Peter Young
;-

IF n_tags(input) EQ 0 THEN BEGIN
  chck=file_search(input,count=count)
  IF count EQ 0 THEN BEGIN
    print,'%PLOT_UPS_SCUPS:  please check your input parameters. Returning...'
    return,-1
  ENDIF ELSE BEGIN
    read_ups,input,upsstr
  ENDELSE
ENDIF ELSE BEGIN
  upsstr=input
ENDELSE 

IF n_elements(ups_ind) EQ 0 AND n_elements(transition) EQ 0 THEN BEGIN
  print,'%PLOT_UPS_SCUPS:  please specify a transition with either UPS_IND= or TRANSITION=. Returning...'
  return,-1
ENDIF

IF n_elements(transition) EQ 2 THEN BEGIN
  l1=transition[0]
  l2=transition[1]
  k=where(upsstr.data.lvl1 EQ l1 AND upsstr.data.lvl2 EQ l2,nk)
  IF nk NE 0 THEN BEGIN
    ups_ind=k[0]
  ENDIF ELSE BEGIN
    print,'%PLOT_UPS_SCUPS:  TRANSITION should be a two element array. Returning...'
    return,-1
  ENDELSE 
ENDIF

temp=upsstr.data[ups_ind].temp
ups=upsstr.data[ups_ind].ups
de=upsstr.data[ups_ind].de
gf=upsstr.data[ups_ind].gf

t_type=burly_get_ttype(gf,de)

c=burly_optimize_c(temp,ups,de,t_type)

scupstr=burly_scale_ups(temp,ups,de,t_type,c)

st=scupstr.st
scups=scupstr.sups

p=plot(st,scups,symbol='+', $
      xrange=[0,1],/xsty)

return,p

END
