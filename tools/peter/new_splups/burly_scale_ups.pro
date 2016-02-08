function burly_scale_ups, temp, ups, de, ttype, c_val, lim=lim, $
                          high_t=high_t

;+
; NAME
;
;    BURLY_SCALE_UPS
;
; EXPLANATION
;
;    Scales a set of temperatures and effective collision strengths
;    (upsilons) using the Burgess & Tully (1992) method.
;
; HISTORY
;
;    Beta 1, 22-Jul-2010
;    Beta 2, 4-Aug-2010
;      added HIGH_T and LIM optional inputs
;-

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

de_in=de
c_ups_curr=c_val

IF n_elements(lim) EQ 0 THEN lim=-1

IF lim EQ -1 THEN lim_set=0 ELSE lim_set=1

;
; The following reduces the temperature range based on the HIGH_T
; optional input.
;
IF n_elements(high_t) EQ 0 THEN high_t=max(temp)
k=where(temp LE high_t)
kt=1.38062d-16*temp[k]/2.1799d-11   ; in Ryd

upsilon=ups
;
CASE ttype OF

1: begin
  xt=kt/de_in
  st=1.-alog(c_ups_curr)/alog(xt+c_ups_curr)
  sups=upsilon/alog(xt+exp(1.))
  IF lim_set EQ 1 THEN BEGIN
    st=[st,1.0]
    sups=[sups,lim]
  ENDIF 
 ;
;  st=[st,1.0]
;  sups=sups
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
4:  BEGIN
  xt=kt/de_in
  st=1.-alog(c_ups_curr)/alog(xt+c_ups_curr)
  sups=upsilon/alog(xt+c_ups_curr)
  IF lim_set EQ 1 THEN BEGIN
    st=[st,1.0]
    sups=[sups,lim]
  ENDIF 
END 
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
END
;
ELSE:  print,' Transition type must be between 1 and 6',t_type_ups
;
ENDCASE

sclups={st: st, sups: sups}

return,sclups

END

