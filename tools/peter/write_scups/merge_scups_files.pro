
PRO merge_scups_files, file1, file2, outfile=outfile

;+
; NAME:
;     MERGE_SCUPS_FILES
;
; PURPOSE:
;     Merge two CHIANTI scups files, removing any duplicate
;     transitions. 
;
; CATEGORY:
;     CHIANTI; file handling.
;
; CALLING SEQUENCE:
;     MERGE_SCUPS_FILES, File1, File2
;
; INPUTS:
;     File1:   This is the reference file. The merged file will
;              contain all of the transitions from this file.
;     File2:   This is the secondary file. A transition from this file
;              will only end up in the merged file if the transition
;              is not present in FILE1.
;
; OPTIONAL INPUTS:
;     OutFile:  The name of the output file. If not specified, then
;               the output file will be "FILE1_merged".
;	
; OUTPUTS:
;     Creates a new scups file with the name FILE1_merged, containing
;     all of the transitions from FILE1 and any transitions in FILE2
;     that are not duplicated in FILE1.
;
; CALLS:
;     READ_SCUPS
;
; EXAMPLE:
;     IDL> merge_scups_files,'fe_7.scups','fe_7.scups_extra'
;
; MODIFICATION HISTORY:
;     Ver.1, 16-Sep-2020, Peter Young
;-

IF n_params() LT 2 THEN BEGIN
  print,'Use:  IDL> merge_scups_files, file1, file2'
  return
ENDIF 

read_scups,file1,scup1
read_scups,file2,scup2

n2=scup2.info.ntrans
file2_index=-1
FOR i=0,n2-1 DO BEGIN
  l1=scup2.data[i].lvl1
  l2=scup2.data[i].lvl2
 ;
  k=where(l1 EQ scup1.data.lvl1 AND l2 EQ scup1.data.lvl2,nk)
  IF nk EQ 0 THEN file2_index=[file2_index,i]
ENDFOR

IF n_elements(file2_index) GT 1 THEN BEGIN
  file2_index=file2_index[1:*]
  n=n_elements(file2_index)
  print,'% MERGE_SCUPS_FILES: there are '+trim(n)+' transitions in FILE2 that are not in FILE1.'
ENDIF ELSE BEGIN
  print,'% MERGE_SCUPS_FILES: there are no transitions in FILE2 that are not in FILE1. Returning...'
  return
ENDELSE 


scup1_data=scup1.data
scup2_data=scup2.data[file2_index]

new_data=[scup1_data,scup2_data]
n=n_elements(new_data)

IF n_elements(outfile) EQ 0 THEN outfile=file1+'_merged'
openw,lout,outfile,/get_lun
FOR i=0,n-1 DO BEGIN
    lim=new_data[i].lim
    IF lim EQ -1 THEN lim='-1' ELSE lim=string(format='(e12.3)',lim)
    printf,lout,format='(i7,i7,e12.3,e12.3,a12,i5,i5,f12.5)',new_data[i].lvl1,new_data[i].lvl2, $
           new_data[i].de,new_data[i].gf,lim,new_data[i].nspl,new_data[i].t_type, $
           new_data[i].c_ups
    form='('+trim(new_data[i].nspl)+'e12.3)'
    printf,lout,format=form,new_data[i].stemp
    printf,lout,format=form,new_data[i].spl
ENDFOR
printf,lout,' -1'

comments=scup1.info.comments
nc=n_elements(comments)
FOR i=0,nc-1 DO BEGIN
  printf,lout,comments[i]
ENDFOR 

free_lun,lout

END
