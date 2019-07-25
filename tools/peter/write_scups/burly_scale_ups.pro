function burly_scale_ups, temp, ups, de, ttype, c_val, lim=lim, $
                          high_t=high_t

;+
; NAME:
;    BURLY_SCALE_UPS
;
; PURPOSE:
;    Scales a set of temperatures and effective collision strengths
;    (upsilons) using the Burgess & Tully (1992) method.
;
; CATEGORY:
;    CHIANTI; data.
;
; CALLING SEQUENCE:
; 
;    Result = BURLY_SCALE_UPS( Temp, Ups, De, Ttype, C_val )
;
; INPUTS:
;    Temp:   An array of temperatures (units: K).
;    Ups:    An array of same size as TEMP containing effective
;            collision strengths.
;    De:     The energy (in Rydbergs) for the atomic transition.
;    Ttype:  The transition type specified as an integer.
;    C_val:  The scaling parameter for the transition.
;
; OPTIONAL INPUTS:
;    Lim:    The high temperature limit for a transition. This is
;            usually only available for allowed transitions.
;    High_T: A temperature value above which any temperatures in TEMP
;            are ignored.
;
; OUTPUTS:
;     A structure with the following tags is returned:
;       .ST    An array of scaled temperatures.
;       .SUPS  An array of scaled upsilons.
;     The arrays will have the same size as the input TEMP, except
;     when LIM= has been specified. Then ST and SUPS have one
;     additional element.
;
;     If a problem is found, then -1 is returned.
;
; MODIFICATION HISTORY:
;     Ver.1, 30-Jan-2017
;        Expanded header, tidied code and added additional
;        comments. No change to functionality, though.
;-


IF n_params() LT 5 THEN BEGIN
  print,'Use:  IDL> str = burly_scale_ups( temp, ups, de, ttype, c_val [, high_t=, lim=] )'
  return,-1
ENDIF


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
; The following scales the temperatures and upsilons based on the
; transition type (TTYPE).
;
; Transition types 1 to 4 are taken from Burgess & Tully
; (1992). Transition types 5 and 6 were introduced in CHIANTI 3 (Dere
; et al. 2001) and CHIANTI 4 (Young et al. 2003), respectively.
;
CASE ttype OF
 ;
  1: BEGIN
    xt=kt/de_in
    st=1.-alog(c_ups_curr)/alog(xt+c_ups_curr)
    sups=upsilon/alog(xt+exp(1.))
    IF lim_set EQ 1 THEN BEGIN
      st=[st,1.0]
      sups=[sups,lim]
    ENDIF 
   ;
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
  6: BEGIN
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

