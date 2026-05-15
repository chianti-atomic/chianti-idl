
FUNCTION ch_check_ion_advanced_model, ion_name, adv_file=adv_file

;+
; NAME:
;     CH_CHECK_ION_ADVANCED_MODEL
;
; PURPOSE:
;     Checks if the specified ion has an advanced model.
;
; CATEGORY:
;     CHIANTI; advanced model.
;
; CALLING SEQUENCE:
;     Result = CH_CHECK_ION_ADVANCED_MODEL( Ion_Name )
;
; INPUTS:
;     Ion_Name:  Name of an ion in CHIANTI format (e.g., 'o_6').
;
; OPTIONAL INPUTS:
;     Adv_File:  The name of the advanced model ion list file. If not
;                specified then the default file in the CHIANTI distribution
;                is used.
;	
; OUTPUTS:
;     Returns 1 if the ion has an advanced model, 0 otherwise.
;
; EXAMPLE:
;     IDL> output=ch_check_ion_advanced_model('o_6')
;
; MODIFICATION HISTORY:
;     Ver.1, 15-May-2026, Peter Young
;-


output=0b

IF n_elements(adv_file) EQ 0 THEN BEGIN
  adv_file=concat_dir(concat_dir(concat_dir(!xuvtop,'ancillary_data'),'advanced_models'),$
                      'advmodel_list.ions')
ENDIF

chck=file_info(adv_file)
IF chck.exists EQ 1 THEN BEGIN 
  info=ch_read_list_ions(adv_file)
  model_ions=info.list_ions
  k=where(trim(ion_name) EQ model_ions,nk)
  IF nk GT 0 THEN output=1b
ENDIF

return,output

END
