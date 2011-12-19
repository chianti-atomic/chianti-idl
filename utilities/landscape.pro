;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an atomic database package for the calculation of
;       continuum and emission line spectra from astrophysical plasmas. It is a 
;       collaborative project involving the Naval Research Laboratory
;       (Washington DC, USA), the Arcetri Observatory (Firenze, Italy), and the
;       Cambridge University (United Kingdom).
;
; NAME:
;	LANDSCAPE
;
; PURPOSE:
;
;	:
;
;
; CALLING SEQUENCE:
;
;       LANDSCAPE
;
;
; INPUTS:
;
;	None	
;
;	
; KEYWORD PARAMETERS:
;
;	None
;
; OUTPUTS:
;
;	None
;
;
;
; COMMON BLOCKS:
;
;	None
;
;
; EXAMPLE:
;
;    to make a postscript file in landscape orientation
;
;             > set_plot,'ps'
;             > landscape
;             > plot,x,y
;             > device,/close
;             > set_plot,'x'
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	March 1996:     Version 2.0
;
;-
Pro landscape
device,xoffset=2.5, yoffset=25.5, /landscape
device,xsize=23.5,ysize=16.5
end
