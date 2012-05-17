;+
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory (USA), the University of Florence (Italy), the
;       University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
;
; NAME: ch_check_str
;
; PURPOSE:
;       To check that an input structure is of the right type
;
; PROCEDURE:
;       This function checks that the input structure has at least the  basic
;       tags that should be present. Two types of basci structures are checked:
;       1) the standard CHIANTI structure, output of the synthetic program
;       CH_SYNTHETIC, that contains line INTENSITIES
;       2) The standard  CHIANTI structure output of MAKE_CHIANTI_SPEC, that contains a
;       synthetic SPECTRUM.
;        
; CATEGORY:
;
;	spectral synthesis.
;
; CALLING SEQUENCE:
;
;       IDL> result=ch_check_str (tran, [/int , /sp])
;
;
; INPUTS: the IDL structure TRAN
;
;
;
; OPTIONAL INPUTS :
;
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; KEYWORDS:
;        intensities  
;        spectrum 
;
; CALLS: required_tags
;
; COMMON BLOCKS: none
;
; RESTRICTIONS:
;
; SIDE EFFECTS:
;
; EXAMPLE:
;
; PREV. HIST. :
;
; WRITTEN     : 
;       Ver.1, 22-May-02, Giulio Del Zanna (GDZ), DAMTP 
;
; MODIFICATION HISTORY:
;
; VERSION     : 1, 22-May-02, GDZ
;
;-


function ch_check_str, tran, intensities=intensities, spectrum=spectrum, $
         error=error

;DEM ???

;check the common TAGS of the two structures:

 taglist = 'IONEQ_LOGT,IONEQ_NAME,IONEQ_REF,WVL_LIMITS,'+$
  'model_name,WVL_UNITS,INT_UNITS,ADD_PROTONS,'+$
  'DATE,VERSION,PHOTOEXCITATION,LINES'

result = required_tags(tran, taglist,  missing_tags=missing_tags)

;Exit
IF NOT result THEN return, result


IF keyword_set(intensities) THEN BEGIN 

taglist = 'LAMBDA,SPECTRUM,UNITS,INSTR_FWHM,BIN_SIZE,'+$
  'ABUND_NAME,ABUND,MIN_ABUND,ABUND_REF'

dummy = required_tags(tran, taglist)
IF dummy THEN return, 0 ELSE return,1


ENDIF ELSE IF keyword_set(spectrum) THEN BEGIN 
 taglist = 'IONEQ_LOGT,IONEQ_NAME,IONEQ_REF,WVL_LIMITS,'+$
  'model_name,WVL_UNITS,INT_UNITS,ADD_PROTONS,'+$
  'DATE,VERSION,PHOTOEXCITATION,LINES,'+$
  'LAMBDA,SPECTRUM,UNITS,INSTR_FWHM,BIN_SIZE,'+$
  'ABUND_NAME,ABUND,MIN_ABUND,ABUND_REF'

result = required_tags(tran, taglist,  missing_tags=missing_tags)
return, result

ENDIF ELSE BEGIN 
return, result
END 



END
