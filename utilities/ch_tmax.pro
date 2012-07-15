FUNCTION ch_tmax, ionname, ioneqname=ioneqname, log=log

;+
; NAME
;
;     CH_TMAX()
;
; PROJECT
;
;     CHIANTI
;
; EXPLANATION
;
;     Returns the temperature of maximum ionization (T_max) for the
;     specified ion.
;
; INPUTS
;
;     IONNAME   The name of the ion in CHIANTI format. E.g., 'fe_13'
;               for Fe XIII.
;
; OPTIONAL INPUTS
;
;     IONEQNAME  The routine reads the default CHIANTI ion balance
;                file (!ioneq_file). To use a different ion balance
;                file, specify the full pathname with this keyword.
;
; KEYWORDS
;
;     LOG     If set, then the logarithm (base 10) of T_max is
;             returned. 
;
; OUTPUTS
;
;     Returns the temperature of maximum ionization in K. If the ion
;     balance file is not found, then a value of -1 is returned.
;
; CALLS
;
;     READ_IONEQ, CONVERTNAME
;
; HISTORY
;
;     Ver.1, 9-May-2013, Peter Young
;         This routine is the same as a previous routine named
;         get_tmax.pro, only I've added the /LOG keyword.
;-

IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> tmax=ch_tmax( ion_name [, /log, ioneqname=] )'
  return,-1.
ENDIF 

IF n_elements(ioneqname) EQ 0 THEN ioneqname=!ioneq_file

IF NOT file_exist(ioneqname) THEN BEGIN
  print,'%GET_TMAX: the specified ionization balance file does not exist. Returning...'
  return,-1.
ENDIF 

read_ioneq,ioneqname,tt,ii,ref

convertname,ionname,iz,ion

ii=reform(ii[*,iz-1,ion-1])

getmax=max(ii,index)


IF keyword_set(log) THEN return,tt[index] ELSE return,10.^tt[index]

END
