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
;	ION2SPECTROSCOPIC
;
; PURPOSE:
;
;	provide identification string
;
; CATEGORY:
;
;	database.
;
; CALLING SEQUENCE:
;
;       ION2SPECTROSCOPIC, Ion, Spectroscopic
;
; INPUTS:
;
;       Ion:   CHIANTI notation for an ion, i.e., 'c_2' for C II	
;
; OUTPUTS:
;
;	Name:  the spectroscopic notation for the ion, i.e. 'C II'
;
; EXAMPLE:
;             > ion2spectroscopic,'fe_13',name
;             > print,name
;             > Fe XIII   
;
;             > ion2spectroscopic,'o_6d',name
;             > print,name
;             > O VI d   
;
; MODIFICATION HISTORY:
; 	Written by:	    Ken Dere
;	September 1999:     Version 3.0
;
;       Ver.2, 11-Dec-2001, Peter Young
;           Revised routine; removed call to repstr
;
;       V.3, 29-May-2002, Giulio Del Zanna (GDZ) 
;              Added output keyword dielectronic.
;
;       V.4, 12-Aug-02, GDZ 
;           Corrected a typo concerning XXII.
;
; VERSION     :   4, 12-Aug-02
;
;-
pro ion2spectroscopic,ion,snote, dielectronic=dielectronic
;

if n_params(0) lt 2 then begin
   print,' '
   print,'    use> ion2spectroscopic,ion,snote'
   print,'    use> ion2filename,''c_6'',snote'
   print,'           giving snote =C VI, for example'
   print,' '
   return
endif

element=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si', $
         'P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co', $
         'Ni','Cu','Zn']

ionstage=['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII', $
          'XIII','XIV','XV','XVI','XVII','XVIII','XIX','XX','XXI','XXII', $
          'XXIII','XXIV','XXV','XXVI','XXVII','XXVIII','XXIX','XXX','XXXI']

bits=str_sep(ion,'_')

 dielectronic=0

IF n_elements(bits) NE 2 THEN BEGIN
   iz=0
   ion=0
ENDIF ELSE BEGIN
   first=strupcase(bits[0])
   ind=where(first EQ strupcase(element))
   IF ind[0] EQ -1 THEN BEGIN
      print,'% Error, ion not in database  ',name
      iz=0 &  ion=0
      return
   ENDIF
   snote=element[ind[0]]
                                ;
   bits2=str_sep(bits[1],'d')
   ind=fix(bits2[0])
   snote=snote+' '+ionstage[ind[0]-1]
   IF n_elements(bits2) EQ 2 THEN BEGIN 
      dielectronic=1
      snote=snote+' d'
   ENDIF 
ENDELSE 

END
