PRO DESCALE_ALL, TEMP, SPLSTR, INDEX, UPS

;+
; EXPLANATION
;
;      This routine descales all types of spline fits into upsilons or
;      rates, i.e., it does both electron upsilons and proton rates,
;      and both 5-point and 10-point splines. In addition it can
;      simultaneously descale several temperatures at once.
;
; INPUTS
;
;      TEMP     Temperature(s), K.
;
;      SPLSTR   Structure output by read_splups.
;
;      INDEX    Index of structure.
;
; OUTPUTS
;
;      UPS      Upsilon value(s) at temperature(s) TEMP.
;
; EXAMPLES
;
;      read_splups,splupsfile,splstr
;      descale_all,[1.e6,2.e6],splstr,5,ups
;      print,ups
;
; HISTORY
;
;      Ver.1, 15-Mar-01, Peter Young
;               adapted from Ken Dere's descale_ups.pro.
;
;      Ver.2, 12-Nov-01, Peter Young
;               added type 6 transitions (for protons)
;
;      Ver.3, 3-Jan-2012, Ken Dere
;               allowed for the number of spline points from 5 through 9
;-


de=splstr[index].de
cc=splstr[index].c_ups
tt=splstr[index].t_type
spl=splstr[index].spl
st=splstr[index].nspl

kte=temp/de/1.57888d5

CASE 1 OF
  (tt EQ 1) OR (tt EQ 4): xt=1 - alog(cc)/(alog(kte + cc))
  (tt EQ 2) OR (tt EQ 3) OR (tt EQ 5) OR (tt EQ 6): xt=kte / (kte +cc)
ENDCASE
;
; this is from Ver. 2
; IF st EQ 5 THEN BEGIN     ; 5-point spline
;   spl=spl[0:4]
;   xs=0.25*findgen(5)
; ENDIF ELSE BEGIN          ; 9-point spline
;   spl=spl[0:8]
;   xs=0.125*findgen(9)
; ENDELSE

;  kpd - want to allow any number of splines
spl=spl[0:st-1]
xs = findgen(st)/float(st-1.)

y2=spl_init(xs,spl)
sups=spl_interp(xs,spl,y2,xt)

CASE tt OF

  1: ups=sups*alog(kte + exp(1.))

  2: ups=sups

  3: ups=sups/(kte+1.)

  4: ups=sups*alog(kte+cc)

  5: ups=sups/(kte)

  6: ups=10.^sups

  ELSE:  print,' t_type ne 1,2,3,4,5,6=',tt,l1,l2

ENDCASE

ups=ups>0.

END
