;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory (USA), the University of Florence (Italy), the
;       University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
;
; NAME:
;	SPECTROSCOPIC2ION
;
; PURPOSE:
;	provide identification string
;
; CATEGORY:
;
;	database.
;
; CALLING SEQUENCE:
;
;       SPECTROSCOPIC2ION,  snote, Ion, dielectronic=dielectronic
;
; INPUTS:
;
;	snote:  the spectroscopic notation for the ion, i.e. 'C II'
;
; OPTIONAL INPUTS: none
;
; KEYWORD PARAMETERS:
;              
; OUTPUTS:
;       Ion:   CHIANTI notation for an ion, i.e., 'c_2' for C II	
;
; OPTIONAL OUTPUTS:  dielectronic (0/1)
;
; EXAMPLE:
;             > spectroscopic2ion, 'O VI d', ion,dielectronic=dielectronic
;             > help, ion, dielectronic
;             > ION             STRING    = 'o_6d'
;             > DIELECTRONIC    INT       =        1
;
; WRITTEN     : 
;       Version 1, Written by: Giulio Del Zanna (GDZ) 25-May-2002
;
; MODIFICATION HISTORY:
;
;       V.2, 12-Aug-02, GDZ 
;           Corrected a typo concerning XXII.
;
; VERSION     :   2, 12-Aug-02
;
;-

PRO spectroscopic2ion, snote, ion, dielectronic=dielectronic

ionstage=['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII', $
          'XIII','XIV','XV','XVI','XVII','XVIII','XIX','XX','XXI','XXII', $
          'XXIII','XXIV','XXV','XXVI','XXVII','XXVIII','XXIX','XXX','XXXI']

dielectronic=0

bits=str_sep(snote, ' ')

IF n_elements(bits) LT  2 THEN BEGIN
   ion=''
   return 
ENDIF ELSE BEGIN
   index = where(ionstage EQ bits[1], nn)
   IF nn NE 1 THEN BEGIN 
      ion = ''
      return
   END 
   ion = trim(strlowcase(bits[0]))+'_'+trim(index+1)
END 

IF n_elements(bits) EQ 3 THEN BEGIN 
   IF bits[2] NE 'd' THEN BEGIN 
      dielectronic=0
      return
   ENDIF  ELSE BEGIN 
      dielectronic=1
      ion =ion+'d'
   END 
END  

END  
