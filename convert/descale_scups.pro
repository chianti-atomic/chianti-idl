

PRO descale_scups, temp, splstr, index, ups

;+
; NAME:
;	DESCALE_SCUPS
;
; PURPOSE:
;	This routine takes a structure containing data from a CHIANTI
;	'scups' file and de-scales the spline fit for a
;	transition's collision strengths to yield the effective
;	collision strengths (upsilons) for the requested temperatures.
;
;       Note that if the de-scaled upsilon is negative, then it is set
;       to zero.
;
; CATEGORY:
;	CHIANTI.
;
; CALLING SEQUENCE:
;       DESCALE_SCUPS, Temp, Splstr, Index, Ups
;
; INPUTS:
;	Temp:	The temperature(s) for which the upsilons are
;       	required. Units: K.
;       Splstr: The structure returned by the routine read_scups.pro.
;       Index:  The index within SPLSTR that identifies the transition
;               for which upsilons are required.
;
; OUTPUTS:
;       Ups     The upsilons for the transition identified by INDEX. 
;
; EXAMPLE:
;       IDL> read_scups, 'fe_12.scups', splstr
;       IDL> descale_scups, [1e6, 2e6], splstr, 0, ups
;
; MODIFICATION HISTORY:
;       Ver.1, 8-Jan-2014, Peter Young.
;          - adapted from descale_all.pro.
;       Ver.2, 20-Jun-2014, Peter Young
;          - fixed structure problem.
;-


de=splstr.data[index].de
stemp=splstr.data[index].stemp
spl=splstr.data[index].spl
cc=splstr.data[index].c_ups
tt=splstr.data[index].t_type
nspl=splstr.data[index].nspl

spl=spl[0:nspl-1]
stemp=stemp[0:nspl-1]

kte=temp/de/1.57888d5

CASE 1 OF
  (tt EQ 1) OR (tt EQ 4): xt=1 - alog(cc)/(alog(kte + cc))
  (tt EQ 2) OR (tt EQ 3) OR (tt EQ 5) OR (tt EQ 6): xt=kte / (kte +cc)
ENDCASE

y2=spl_init(stemp,spl)
sups=spl_interp(stemp,spl,y2,xt)

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

