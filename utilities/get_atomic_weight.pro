
FUNCTION GET_ATOMIC_WEIGHT, IZ
;+
; PROJECT:  CHIANTI
;
;       CHIANTI is an atomic database package for the calculation of
;       continuum and emission line spectra from astrophysical plasmas. It is a 
;       collaborative project involving the Naval Research Laboratory
;       (Washington DC, USA), the Arcetri Observatory (Firenze, Italy), and the
;       Cambridge University (United Kingdom).
;
;
; NAME:
;	GET_ATOMIC_WEIGHT
;
;
; PURPOSE:
;
;       to produce  an atomic weight. 
;
; PROCEDURE:
;
;
; This routine converts the ion atomic number into an atomic weight. Where 
; an element has more than one isotope, the weight is that of the most 
; common.
;
; I've written this routine in order to apply a thermal broadening to 
; lines in a CHIANTI synthetic spectrum (see routine make_chianti_spec.pro).
;
;
; CALLING SEQUENCE:
;
;       IDL>
;
; EXAMPLES:
;
; INPUTS:
;        IZ 
;
; OPT. INPUTS :
;
; OUTPUTS:
;       the atomic weight
;
; OPTIONAL OUTPUTS:
;
; KEYWORDS:
;
; CALLS: 
;
; COMMON:
;
; RESTRICTIONS:
;            I go up to zinc (iz=30)
;
;
; SIDE EFFECTS:
;
; CATEGORY:
;	
;	synthetic spectra
;
;
; PREV. HIST. :
;
; WRITTEN     : 
;       Ver.1, 22-Jun-00, Peter Young (PRY)
;
; MODIFICATION HISTORY:
;
;  VERSION     : V 1, 22-Jun-00, Peter Young (PRY)
;
;-

weights = [1,4,7,9,11,12,14,16,19,20,23,24,27,28,31,32,35,40,39,40,$
           45,48,51,52,55,56,59,58,63,64] 

RETURN, weights[iz-1]

END
