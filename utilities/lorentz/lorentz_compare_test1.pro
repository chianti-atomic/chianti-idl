
PRO lorentz_compare_test1, dir1, dir2, d1=d1, d2=d2


;+
; NAME:
;     LORENTZ_COMPARE_TEST1
;
; PURPOSE:
;     Compares the radiative loss curves stored in the Lorentz Test 1
;     result files.
;
; CATEGORY:
;     CHIANTI; Lorentz tests.
;
; CALLING SEQUENCE:
;     LORENTZ_COMPARE_TEST1, Dir1, Dir2
;
; INPUTS:
;     Dir1:   Directory containing the results from Lorentz Test 2.
;     Dir2:   Directory containing the results from a different run of
;             Lorentz Test 2.
;
; OUTPUTS:
;     For each ion, the routine compares the ion fractions between the
;     two data files. If a ratio has changed by more than 10%, then a
;     message is printed to the screen giving the ion and the three
;     ratio values.
;
; OPTIONAL OUTPUTS:
;     D1:   The data structure for DIR1.
;     D2:   The data structure for DIR2.
;
; EXAMPLE:
;     Here results from CHIANTI Version 9.0.0 are stored in
;     '09_00_00', and from Version 10.0.0 are stored in '10_00_00'. 
;
;     IDL> lorentz_compare_test1,'09_00_00','10_00_00'
;
; MODIFICATION HISTORY:
;     Ver.1, 06-Nov-2020, Peter Young
;-


chck=file_search(dir1,'lorentz_test_1_chianti.txt',count=count)
IF count EQ 1 THEN BEGIN
   file1=chck[0]
ENDIF ELSE BEGIN
   print,'% LORENTZ_COMPARE_TEST1: File not found in DIR1. Returning...'
   return
ENDELSE

str1=''
openr,lin1,file1,/get_lun
FOR i=0,5 DO readf,lin1,str1
print,'***File 1 info***'
FOR i=0,1 DO BEGIN
   readf,lin1,str1
   print,'  ',str1
ENDFOR
str={atom: 0, charge: 0, data: fltarr(3)}
d1=0
WHILE eof(lin1) NE 1 DO BEGIN
   readf,lin1,format='(2i3,3f12.3)',str
   IF n_tags(d1) EQ 0 THEN d1=str ELSE d1=[d1,str]
ENDWHILE 
free_lun,lin1

;---
chck=file_search(dir2,'lorentz_test_1_chianti.txt',count=count)
IF count EQ 1 THEN BEGIN
   file2=chck[0]
ENDIF ELSE BEGIN
   print,'% LORENTZ_COMPARE_TEST1: File not found in DIR2. Returning...'
   return
ENDELSE

str1=''
openr,lin2,file2,/get_lun
FOR i=0,5 DO readf,lin2,str1
print,'***File 2 info***'
FOR i=0,1 DO BEGIN
   readf,lin2,str1
   print,'  ',str1
ENDFOR
str={atom: 0, charge: 0, data: fltarr(3)}
d2=0
WHILE eof(lin2) NE 1 DO BEGIN
   readf,lin2,format='(2i3,3f12.3)',str
   IF n_tags(d2) EQ 0 THEN d2=str ELSE d2=[d2,str]
ENDWHILE 
free_lun,lin2

;---
print,''
n1=n_elements(d1)
n2=n_elements(d2)
IF n1 NE n2 THEN BEGIN
   print,'**Different numbers of elements** '+trim(n1)+' '+trim(n2)
ENDIF 

FOR i=0,n1-1 DO BEGIN
   k=where(d2.atom EQ d1[i].atom AND d2.charge EQ d1[i].charge,nk)
   IF nk NE 0 THEN BEGIN
      frac1=10.^d1[i].data
      frac2=10.^d2[k[0]].data
      ratio=frac2/frac1
      j=where(ratio GT 1.1 OR ratio LT 0.9,nj)
      IF nj GT 0 THEN BEGIN
         print,'   ** Problem, '+trim(d1[i].atom)+' '+trim(d1[i].charge)+': '+ $
               string(format='(3f10.3)',ratio)
      ENDIF 
   ENDIF ELSE BEGIN
      print,'** No match for '+trim(d1[i].atom)+' '+trim(d1[i].charge)
   ENDELSE 
ENDFOR 

END
