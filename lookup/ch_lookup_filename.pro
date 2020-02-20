
FUNCTION ch_lookup_filename, ion_name, dir_lookup=dir_lookup, verbose=verbose, status=status

;+
; NAME:
;     CH_LOOKUP_FILENAME
;
; PURPOSE:
;     Creates the population lookup table filename for the specified
;     ion. 
;
; CATEGORY:
;     CHIANTI; lookup tables.
;
; CALLING SEQUENCE:
;     Result = CH_LOOKUP_FILENAME( Ion_Name )
;
; INPUTS:
;     Ion_Name: The name of an ion in CHIANTI format. For example,
;               'fe_13'.
;
; OPTIONAL INPUTS:
;     Dir_Lookup: Specifies the directory where the lookup table is
;                 expected to be
;	
; KEYWORD PARAMETERS:
;     VERBOSE:  Prints information messages to the IDL window.
;
; OUTPUTS:
;     Returns the name of the population lookup table for the specified
;     ion. The filename will be 'pop_lookup_[ion_name].txt'. The
;     routine checks if the environment variable $CHIANTI_LOOKUP is
;     defined. If yes, then the filename is concatenated with it. If
;     DIR_LOOKUP is specified, then the filename is concatenated
;     with it instead.
;
; OPTIONAL OUTPUTS:
;     Status:  The routine checks if the table exists. If yes, then
;              status=1 otherwise status=0. 
;
; EXAMPLE:
;     IDL> fname=ch_lookup_filename('o_6')
;     IDL> fname=ch_lookup_filename('o_6',dir_lookup='~/lookup',status=status)
;
; MODIFICATION HISTORY:
;     Ver.1, 18-Dec-2019, Peter Young
;-

IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> filename=ch_lookup_filename(ion_name [, dir_lookup= ])'
  return,''
ENDIF 

;
; Note that the lookup file name has to be in a standard format.
;
filename='pop_lookup_'+trim(ion_name)+'.txt'
chianti_lookup=getenv('CHIANTI_LOOKUP')
IF chianti_lookup NE '' AND keyword_set(verbose) THEN BEGIN
  print,'% CH_LOOKUP_FILENAME: the environment variable $CHIANTI_LOOKUP is defined, with value'
  print,'                      '+chianti_lookup
ENDIF
;
IF n_elements(dir_lookup) NE 0 THEN BEGIN
  chck=file_info(dir_lookup)
  IF chck.exists EQ 0 THEN BEGIN
    print,'% CH_LOOKUP_FILENAME: The input DIR_LOOKUP was specified, but the directory was not found. Returning...'
    return,-1
  ENDIF 
ENDIF ELSE BEGIN
  IF chianti_lookup NE '' THEN dir_lookup=chianti_lookup
ENDELSE
;
IF n_elements(dir_lookup) NE 0 THEN filename=concat_dir(dir_lookup,filename)

chck=file_info(filename)
IF chck.exists EQ 0 THEN status=0 ELSE status=1

return,filename

END
