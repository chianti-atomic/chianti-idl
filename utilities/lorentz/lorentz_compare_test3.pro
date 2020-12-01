
PRO lorentz_compare_test3, dir1, dir2, temp=temp, d1=d1, d2=d2


;+
; NAME:
;     LORENTZ_COMPARE_TEST3
;
; PURPOSE:
;     Compare results of Lorentz Test 3 between different versions of
;     CHIANTI.
;
; CATEGORY:
;     CHIANTI; Lorentz tests.
;
; CALLING SEQUENCE:
;     LORENTZ_COMPARE_TEST3, Dir1, Dir2
;
; INPUTS:
;     Dir1:   Directory containing the results from Lorentz Test 3.
;     Dir2:   Directory containing the results from a different run of
;             Lorentz Test 3.
;
; OPTIONAL INPUTS:
;     Temp:   The temperature in MK units for which the comparison is
;             done (the data files are only for specific
;             temperatures). The default is 6 MK.
;	
; OUTPUTS:
;     For each of the 100 lines in the DIR2 file, the DIR2 file is
;     searched for the same wavelength (within 10 milli-angstroms). If
;     found, the intensity ratio is printed. Changes of more than 10%
;     are indicated with a '**'. If the wavelength is not found, then
;     a message is printed. It is recommended that the routine is run
;     twice, swapping DIR1 and DIR2.
;
; OPTIONAL OUTPUTS:
;     D1:   The data structure for the DIR1 file.
;     D2:   The data structure for the DIR2 file.
;
; EXAMPLE:
;     Here the results stored in directories 09_00_00 and 10_00_00
;     (i.e., CHIANTI versions 9 and 10) are compared.
;
;     IDL> lorentz_compare_test3, '09_00_00', '10_00_00'
;
; MODIFICATION HISTORY:
;     Ver.1, 05-Nov-2020, Peter Young
;-


IF n_elements(temp) EQ 0 THEN temp=6.0
temp_chck=[1.0, 6.0, 46.0]
getmin=min(abs(temp-temp_chck),imin)

print,'% LORENTZ_COMPARE_TEST3: using temperature = '+trim(string(format='(f10.1)',temp_chck[imin]))+' MK.'

txt=['1MK','6MK','46MK']

chck=file_search(dir1,'lorentz_test_3_chianti_'+txt[imin]+'.txt',count=count)
IF count EQ 1 THEN BEGIN
   file1=chck[0]
ENDIF ELSE BEGIN
   print,'% LORENTZ_COMPARE_TEST3: File not found in DIR1. Returning...'
   return
ENDELSE

str1=''
openr,lin1,file1,/get_lun
FOR i=0,11 DO readf,lin1,str1
;
str={ind: 0, wvl: 0., z: 0, charge: 0, int: 0.}
d1=0
WHILE eof(lin1) NE 1 DO BEGIN
   readf,lin1,format='(i3,f12.3,2i4,f10.3)',str
   IF n_tags(d1) EQ 0 THEN d1=str ELSE d1=[d1,str]
ENDWHILE 
free_lun,lin1

;--
chck=file_search(dir2,'lorentz_test_3_chianti_'+txt[imin]+'.txt',count=count)
IF count EQ 1 THEN BEGIN
   file2=chck[0]
ENDIF ELSE BEGIN
   print,'% LORENTZ_COMPARE_TEST3: File not found in DIR2. Returning...'
   return
ENDELSE

str1=''
openr,lin2,file2,/get_lun
FOR i=0,11 DO readf,lin2,str1
;
str={ind: 0, wvl: 0., z: 0, charge: 0, int: 0.}
d2=0
WHILE eof(lin2) NE 1 DO BEGIN
   readf,lin2,format='(i3,f12.3,2i4,f10.3)',str
   IF n_tags(d2) EQ 0 THEN d2=str ELSE d2=[d2,str]
ENDWHILE 
free_lun,lin2

;--
n1=n_elements(d1)

FOR i=0,n1-1 DO BEGIN
   getmin=min(abs(d2.wvl-d1[i].wvl),imin)
   zion2name,d1[i].z,d1[i].charge+1,name
   IF getmin LE 0.01 THEN BEGIN 
      ratio=10.^d2[imin].int/10.^d1[i].int
      IF (ratio GE 1.1) OR (ratio LT 0.9) THEN extra_txt=' **' ELSE extra_txt=''
      print,format='(i3,a7,f10.3,2f10.3,f10.3,a3)',i+1,name,d1[i].wvl,d1[i].int,d2[imin].int,ratio,extra_txt
   ENDIF ELSE BEGIN
      print,format='(i3,a7,f10.3,"  no match found")',i+1,name,d1[i].wvl
   ENDELSE 
ENDFOR 

END
