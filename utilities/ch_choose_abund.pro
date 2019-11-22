
FUNCTION ch_choose_abund

;+
; NAME:
;     CH_CHOOSE_ABUND
;
; PURPOSE:
;     Allows the user to choose an element abundance file from the
;     selection offered in CHIANTI through a widget interface.
;
; CATEGORY:
;     CHIANTI; abundances; file selection.
;
; CALLING SEQUENCE:
;     Result = CH_CHOOSE_ABUND()
;
; INPUTS:
;     None.
;
; OUTPUTS:
;     The name of CHIANTI element abundance file, as selected by the
;     user from the widget. If no selection was made, then an empty
;     string is returned. 
;
; MODIFICATION HISTORY:
;     Ver.1, 26-Apr-2019, Peter Young
;-

dir=concat_dir(!xuvtop,'abundance')
dir=expand_path(dir)
output=dialog_pickfile(path=dir,filter='*.abund',title='Select Abundance File')

return,output

END

