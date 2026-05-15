
FUNCTION ch_compare_ioneq, ioneq_file1, ioneq_file2, verbose=verbose

;+
; NAME:
;     CH_COMPARE_IONEQ
;
; PURPOSE:
;     Compares to CHIANTI ionization equilibrium files (ioneq) and gives
;     information about the differences.
;
; CATEGORY:
;     CHIANTI; ioneq; validation.
;
; CALLING SEQUENCE:
;     Result = CH_COMPARE_IONEQ( File1, File2 )
;
; INPUTS:
;     Ioneq_File1:  String giving the name of the reference file.
;     Ioneq_File2:  String giving the name of the comparison file.
;	
; KEYWORD PARAMETERS:
;     VERBOSE:  If set, then the lists of ions for which there are diffences
;               are printed to the IDL window.
;
; OUTPUTS:
;     An IDL structure array (an entry for each ion) with the tags:
;      .ion   Ion name in CHIANTI format (e.g., 'o_6').
;      .log_tmax1 Log Tmax value for the ion in the reference file.
;      .log_tmax2 Log Tmax value for the ion in the comparison file.
;      .max_val1  The ion fraction at the Tmax value for reference file.
;      .max_val2  The ion fraction at the Tmax value for comparison file.
;      .temp_flag Takes value 1 if there is a discrepancy for Tmax.
;      .frac_flag Takes value 1 if there is a discrepancy for peak fraction.
;
; EXAMPLE:
;     Compare with the current default ioneq file:
;
;     IDL> r=ch_compare_ioneq(!ioneq_file,'new.ioneq')
;
; MODIFICATION HISTORY:
;     Ver.1, 15-May-2026, Peter Young
;-

IF n_params() LT 2 THEN BEGIN
  print,'Use:  IDL> r=ch_compare_ioneq( file1, file2 [, /verbose] )'
  return,-1
ENDIF

chck=file_info(ioneq_file1)
IF chck.exists EQ 0 THEN BEGIN
  message,/info,/cont,'ioneq_file1 does not exist! Please check your inputs.'
  return,-1.
ENDIF 

chck=file_info(ioneq_file2)
IF chck.exists EQ 0 THEN BEGIN
  message,/info,/cont,'ioneq_file2 does not exist! Please check your inputs.'
  return,-1.
ENDIF 

read_ioneq,ioneq_file1,tt1,ii1,ref
read_ioneq,ioneq_file2,tt2,ii2,ref

nt1=n_elements(tt1)
nt2=n_elements(tt2)

IF nt1 NE nt2 THEN BEGIN
  message,/info,/cont,'The two files are tabulated for two different temperature ranges. Returning...'
  return,-1
ENDIF 

s=size(ii1,/dim)
n_elt=s[1]
n_ion=s[0]

str={ ion: '', $
      log_tmax1: 0., $
      log_tmax2: 0., $
      max_val1: 0., $
      max_val2: 0., $
      temp_flag: 0b, $
      frac_flag: 0b}

output=0

FOR i=0,n_elt-1 DO BEGIN
  FOR j=0,i+1 DO BEGIN
    zion2name,i+1,j+1,name
   ;
    str.ion=name
    getmax=max(ii1[*,i,j],imax)
    str.log_tmax1=tt1[imax]
    str.max_val1=getmax
    getmax=max(ii2[*,i,j],imax)
    str.log_tmax2=tt2[imax]
    str.max_val2=getmax
   ;
    IF str.log_tmax1 NE str.log_tmax2 THEN str.temp_flag=1b ELSE str.temp_flag=0b
    perc=(str.max_val2-str.max_val1)/str.max_val1*100.
    IF abs(perc) GT 10. THEN str.frac_flag=1b ELSE str.frac_flag=0b
   ;
    IF n_tags(output) EQ 0 THEN output=str ELSE output=[output,str]
  ENDFOR 
ENDFOR

print,format='("No. of ions: ",i3)',n_elements(output)

k=where(output.temp_flag EQ 1,nk)
print,format='("No. of ions for which T_max has shifted: ",i3)',nk
IF keyword_set(verbose) THEN BEGIN
  FOR i=0,nk-1 DO BEGIN
    j=k[i]
    print,format='(5x,a7,2f10.2)',strpad(output[j].ion,7,fill=' '), $
          output[j].log_tmax1,output[j].log_tmax2
  ENDFOR 
ENDIF 

k=where(output.frac_flag EQ 1,nk)
print,format='("No. of ions for which peak ion fraction has changed by more than 10%: ",i3)',nk
IF keyword_set(verbose) THEN BEGIN
  FOR i=0,nk-1 DO BEGIN
    j=k[i]
    print,format='(5x,a7,2f10.2)',strpad(output[j].ion,7,fill=' '), $
          output[j].max_val1,output[j].max_val2
  ENDFOR 
ENDIF 

k=where(abs(output.max_val2-output.max_val1)/output.max_val1 GE 0.01,nk)
print,format='("No. of ions for which peak ion fraction has changed by more than 1%: ",i3)',nk
IF keyword_set(verbose) THEN BEGIN
  FOR i=0,nk-1 DO BEGIN
    j=k[i]
    print,format='(5x,a7,2f10.2)',strpad(output[j].ion,7,fill=' '), $
          output[j].max_val1,output[j].max_val2
  ENDFOR 
ENDIF 

return,output

END
