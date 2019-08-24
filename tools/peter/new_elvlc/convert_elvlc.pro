
PRO convert_elvlc, input_file, output_file=output_file

;+
; NAME:
;      CONVERT_ELVLC
;
; PURPOSE:
;      This routine converts an old format ELVLC file (pre-v8) to the
;      new format file described in CHIANTI Technical Report No. 1.
;
; CATEGORY:
;      CHIANTI; file formats; conversion.
;
; CALLING SEQUENCE:
;      CONVERT_ELVLC, Input_File
;
; INPUTS:
;      Input_File:  The name of a CHIANTI .elvlc file.
;
; OPTIONAL INPUTS:
;      Output_File: The name of the output file. If not specified,
;                   then the output file is the same as the input
;                   file, only with '_new' appended.
;
; OUTPUTS:
;      The input file is reformatted and written to a new file in the
;      same directory.
;
; CALLS:
;      READ_ELVLC_OLD
;
; MODIFICATION HISTORY:
;      Ver.1, 31-Jul-2017, Peter Young
;-

chck=file_search(input_file,count=count)
IF count EQ 0 THEN BEGIN
  print,'%CONVERT_ELVLC: the input file was not found. Returning...'
  return
ENDIF 

IF n_elements(output_file) EQ 0 THEN output_file=input_file+'_new'

read_elvlc_old,input_file,l1a,term,confa,ssa,lla,jja,ecma,eryda,ecmtha,erydth,ref, $
               desig=desig, spda=spda


openw,lout,output_file,/get_lun

n=n_elements(l1a)
FOR i=0,n-1 DO BEGIN
  confstr=strpad(trim(desig[i]),30,fill=' ')
  IF ecma[i] EQ 0. AND i NE 0 THEN BEGIN
    en_form='(i15)'
    ecma[i]=-1
  ENDIF ELSE BEGIN
    en_form='(f15.3)'
  ENDELSE 
  printf,lout,format='(i7,a30,5x,i5,a5,f5.1,'+en_form+',f15.3)', $
         l1a[i],confstr,ssa[i],spda[i],jja[i],ecma[i],ecmtha[i]
ENDFOR 

printf,lout,' -1'
printf,lout,'%file: ',input_file
n=n_elements(ref)
FOR i=0,n-1 DO BEGIN
  printf,lout,ref[i]
ENDFOR
printf,lout,''
printf,lout,'File converted to v.8 format, '+systime()

free_lun,lout

END
