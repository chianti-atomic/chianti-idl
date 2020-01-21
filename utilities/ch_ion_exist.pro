
FUNCTION ch_ion_exist, ionname

;+
; NAME:
;	CH_ION_EXIST
;
; PURPOSE:
;       Checks if the specified ion exists in CHIANTI, i.e., that the
;       collisional/radiative data exist for population modeling.
;
; CATEGORY:
;	CHIANTI; file management.
;
; CALLING SEQUENCE:
;	Result = CH_ION_EXIST( IonName )
;
; INPUTS:
;	IonName: The name of an ion in CHIANTI format, e.g., 'fe_13'
;                for Fe XIII.
;
; OUTPUTS:
;       Returns 1 if the ion exists or 0 otherwise.
;
; EXAMPLE:
;       IDL> result=ch_ion_exist('fe_13')
;       IDL> result=ch_ion_exist(
;
; MODIFICATION HISTORY:
; 	Ver.1, 15-Aug-2016, Peter Young
;-

IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> result = ch_ion_exist ( ionname )'
  return,0
ENDIF 

convertname,ionname,iz,ion
IF iz EQ 0 THEN return,0
zion2filename,iz,ion,basename

chck1=file_search(basename+'.wgfa',count=count1)
chck2=file_search(basename+'.elvlc',count=count2)
chck3=file_search(basename+'.scups',count=count3)

count=count1+count2+count3


IF count LT 3 THEN return,0 ELSE return,1

END
