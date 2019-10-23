
FUNCTION burly_optimize_c, temp, ups, de, t_type, plot=plot, method=method, $
                           c_range=c_range

;+
; NAME:
;     BURLY_OPTIMIZE_C
;
; PURPOSE:
;     Takes a set of upsilon data and computes an optimum scaling
;     parameter (C). The method for determining the optimum C-value is
;     set by the input METHOD.
;
; CATEGORY:
;     CHIANTI; data assessment.
;
; CALLING SEQUENCE:
;     Result = BURLY_OPTIMIZE_C( Temp, Ups, De, T_Type )
;
; INPUTS:
;     Temp:   An array of temperatures (units: K).
;     Ups:    An array of upsilons (same size as TEMP).
;     De:     Energy for transition (Rydbergs).
;     T_Type: CHIANTI transition type.
;
; OPTIONAL INPUTS:
;     Method: Specifies the method used to compute the optimum
;             C-value. The options are:
;              0: minimize the gradients between points (i.e., makes a
;                 "flat" curve). This is the default.
;              1: maximize the minimum separation between two points
;              2: hybrid: takes the midway C-value between
;                 method 0 and method 1.
;              3: minimize the maximum separation between two points.
;              4: hybrid: takes the midway C-value between
;                 method 0 and method 3.
;
; KEYWORD PARAMETERS:
;     PLOT:   If set, then a plot is made showing the scaled upsilons
;             using the C-value derived from BURLY_FIND_MID_C
;             (crosses), and the scaled upsilons using the optimized
;             C-value (squares). The right panel shows all C-values
;             that were tried (X-axis) plotted with the maximum
;             gradient (Y-axis). The dashed vertical line shows the
;             MID_C value, and the vertical line shows the optimized
;             C-value.
;
; OUTPUTS:
;     Returns the optimum C-value for the transition. If a problem is
;     found, then -1 is returned.
;
; OPTIONAL OUTPUTs:
;     C_Range:  The range of C-values used for the optimization
;               process. 
;
; EXAMPLE:
;     IDL> read_ups,'mg_5.ups',upsstr
;     IDL> temp=upsstr.data[304].temp
;     IDL> ups=upsstr.data[304].ups
;     IDL> de=upsstr.data[304].de
;     IDL> c=burly_optimize_c(temp,ups,de,4)  ; transition type 4
;
; PROCEDURE:
;     The old manual way to scaling upsilons was to adjust C to
;     minimize steep gradients in the scaled data points. I mimic this
;     here by computing the gradients from point-to-point, and picking
;     out the maximum gradient (stored in GRAD_ARR).
;
;     My initial guess for C is obtained from BURLY_FIND_MID_C, and I
;     then compute an array of C-values around this value.
;
;     Note that when computing the gradients, I normalize sups by
;     dividing by max(sups). This is because some transition types
;     (e.g. 4) scale the absolute value of sups with the scaling
;     parameter. 
;
; MODIFICATION HISTORY:
;     Ver.1, 8-Feb-2017, Peter Young
;     Ver.2, 4-Aug-2017, Peter Young
;        Added the METHOD= input for selecting the method used to
;        compute the optimum C-value.
;     Ver.3, 17-Aug-2017, Peter Young
;        Added optional output C_RANGE.
;-


IF n_params() LT 4 THEN BEGIN
  print,'Use:  IDL> c = burly_optimize_c(temp,ups,de,t_type [,/plot, method=, c_range= ])'
  return,-1
ENDIF 


;
; This specifies the default method for computing the optimum
; C-value. 
;
IF n_elements(method) EQ 0 THEN method=0


;
; Compute my starting C-value (mid_c)
;
mid_c=burly_find_mid_c(temp, de, t_type)

;
; Apply mid_c to scale the input upsilon data.
;
scupstr=burly_scale_ups(temp, ups, de, t_type, mid_c)
sups=scupstr.sups
st=scupstr.st


;
; This is the maximum separation between two scaled temperature values
; allowed by the routine.
;
st_step=1.0/(n_elements(st)+1)
st_step_limit=st_step*2.5
;IF n_elements(st_step_limit) EQ 0 THEN st_step_limit=0.30


IF keyword_set(plot) THEN BEGIN
  !p.multi=[0,2,1]
  plot,st,sups/max(sups),psym=1,charsize=1.5, $
       xtitle='Scaled temperature', $
       ytitle='Normalized scaled upsilon',symsize=1.5, $
       yrange=[0,1.10],/ysty

ENDIF 
 
nt=n_elements(temp)
i=indgen(nt-1)
j=indgen(nt-1)+1

;
; Compute gradients for mid_c
;
grad=(sups[j]-sups[i])/(st[j]-st[i])

;
; I consider 'ref_grad' to be an acceptable gradient, and only do the
; optimization procedure if the maximum gradient in the mid_c scaled
; data-set is larger than this.
;
ref_grad=max(sups)/0.5
maxgrad=max(abs(grad))

;
; I've switched off ref_grad (by setting it to -1) so the
; optimization procedure will always proceed.
;
ref_grad=-1

IF maxgrad GT ref_grad THEN BEGIN
 ;
 ; Compute an array of C-values for the optimization procedure.
 ;
  nc=101
  c=(findgen(nc)/float(nc-1)-0.9)*mid_c+mid_c
 ;
 ; For type 1 and 4, C must be greater than 1, so modify C.
 ;
  IF t_type EQ 1 OR t_type EQ 4 THEN c=c+(1-min(c))+0.001
 ;
 ; Run through all C-values and store the maximum gradient in the
 ; array GRAD_ARR.
 ;
  grad_arr=fltarr(nc)
  max_st_step=fltarr(nc)
  min_st_step=fltarr(nc)
  FOR k=0,nc-1 DO BEGIN
    scupstr=burly_scale_ups(temp, ups, de, t_type, c[k])
    sups=scupstr.sups
   ;
   ; IMPORTANT: I normalize the upsilons using max(sups)
    sups=sups/max(sups)
    st=scupstr.st
    grad=(sups[j]-sups[i])/(st[j]-st[i])
    grad_arr[k]=max(abs(grad))
   ;
    stx=[0.,st,1.]
    nstx=n_elements(stx)
    st2=stx[1:*]
    st1=stx[0:nstx-2]
    st_step=st2-st1
    max_st_step[k]=max(st_step)
    min_st_step[k]=min(st_step)
  ENDFOR
 ;
 ; These are the methods for computing the optimimum C-value. The
 ; value of METHOD determines which one is used. 
 ;
  getmin=min(grad_arr,imin0)
  getmax=max(min_st_step,imin1)
  getmin=min(max_st_step,imin2)
  hy1_imin=round(mean([imin0,imin1]))
  hy2_imin=round(mean([imin0,imin2]))
 ;
  CASE method OF
    1: imin=imin1
    2: imin=hy1_imin
    3: imin=imin2
    4: imin=hy2_imin
    ELSE: imin=imin0
  ENDCASE 
 ;
 ; Now do scaling of data using the optimum C-value.
 ;
  scupstr=burly_scale_ups(temp, ups, de, t_type, c[imin])
 ;
  IF keyword_set(plot) THEN BEGIN 
    sups=scupstr.sups
    st=scupstr.st
    oplot,st,sups/max(sups),psym=6,symsize=1.5
    plot,c,grad_arr,charsize=1.5, $
         xtitle='Scaling parameter', $
         ytitle='Maximum gradient'
    yr=!y.crange
    oplot,mid_c*[1,1],[0,1e5],line=2
    xyouts,mid_c,0.90*yr[1]+0.10*yr[0],'Mid-C',charsize=1.5
    oplot,c[imin]*[1,1],[0,1e5]
    xyouts,c[imin],0.95*yr[1]+0.05*yr[0],'Optimum-C',charsize=1.5
  ENDIF 
 ;
  c_range=c
  return,c[imin]
ENDIF ELSE BEGIN
  return,mid_c
ENDELSE 

END
