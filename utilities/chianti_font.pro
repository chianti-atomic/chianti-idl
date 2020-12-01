
PRO chianti_font, font, big=big, fixed=fixed

;+
; NAME:
;
;      CHIANTI_FONT
;
; PURPOSE:
;
;      Generates standard fonts for CHIANTI GUIs suitable for both Unix and
;      Windows operating systems.
;
; CATEGORY:
;
;      Widgets, fonts
;
; CALLING SEQUENCE:
;
;      CHIANTI_FONT, FONT [, /BIG, /FIXED ]
;
; INPUTS:
;
;      None.
;
; OPTIONAL INPUTS:
;
;      None.
;
; KEYWORD PARAMETERS:
;
;      BIG    Output a descriptor for a large font.
;
;      FIXED  Output a descriptor for a fixed-width font.
;
; OUTPUTS:
;
;      FONT   A descriptor for a font suitable for passing to IDL widget
;             routines.
;
; OPTIONAL OUTPUTS:
;
;      None.
;
; COMMON BLOCKS:
;
;      None.
;
; SIDE EFFECTS:
;
;      None.
;
; RESTRICTIONS:
;
;      Has not been tried with a MAC OS.
;
; PROCEDURE:
;
;      CHIANTI_FONT, FONT [, /BIG, /FIXED ]
;
; EXAMPLE:
;
;      IDL> chianti_font,font
;      IDL> print,font
;      Arial*bold*16
;
; MODIFICATION HISTORY:
;
;      Ver.1, 6-Aug-2003, Peter Young
;-


CASE !version.os_family OF

  'unix': BEGIN
    IF keyword_set(fixed) THEN fstr='-*-courier-' ELSE $
         fstr='-adobe-helvetica-'
    IF keyword_set(big) THEN str='18' ELSE str='12'
    font=fstr+'bold-r-*-*-'+str+'-*'
  END

  ELSE: BEGIN
    IF keyword_set(fixed) THEN fstr='Courier' ELSE $
         fstr='Arial'
    IF keyword_set(big) THEN str='20' ELSE str='16'
    font=fstr+'*bold*'+str
  END

ENDCASE

END

