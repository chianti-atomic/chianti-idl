
FUNCTION lorentz_compare_test2, dir1, dir2, z=z


;+
; NAME:
;     LORENTZ_COMPARE_TEST2
;
; PURPOSE:
;     Compares the radiative loss curves stored in the Lorentz Test 2
;     result files.
;
; CATEGORY:
;     CHIANTI; Lorentz tests.
;
; CALLING SEQUENCE:
;     Result = LORENTZ_COMPARE_TEST2 ( Dir1, Dir2 )
;
; INPUTS:
;     Dir1:   Directory containing the results from Lorentz Test 2.
;     Dir2:   Directory containing the results from a different run of
;             Lorentz Test 2.
;
; OPTIONAL INPUTS:
;     Z:      Atomic number of element to be plotted. If not set, then
;             Z=1 (hydrogen) is assumed.
;
; OUTPUTS:
;     An IDL plot object showing the radiative loss curve from DIR2 in
;     blue, and the curve from DIR1 in black.
;
; EXAMPLE:
;     Here results from CHIANTI Version 9.0.0 are stored in
;     '09_00_00', and from Version 10.0.0 are stored in '10_00_00'. 
;
;     IDL> p=lorentz_compare_test2( '09_00_00', '10_00_00', z=26)
;
; MODIFICATION HISTORY:
;     Ver.1, 04-Nov-2020, Peter Young
;-


chck=file_search(dir1,'lorentz_test_2_chianti.txt',count=count)
IF count EQ 1 THEN BEGIN
   file1=chck[0]
ENDIF ELSE BEGIN
   print,'% LORENTZ_COMPARE_TEST2: File not found in DIR1. Returning...'
   return,-1
ENDELSE

str1=''
openr,lin1,file1,/get_lun
FOR i=0,3 DO readf,lin1,str1
print,'***File 1 info***'
readf,lin1,str1
print,'  ',str1
;
ltemp=fltarr(51)
readf,lin1,format='(3x,51f10.0)',ltemp
str={z: 0, data: fltarr(51)}
d1=0
WHILE eof(lin1) NE 1 DO BEGIN
   readf,lin1,format='(i3,51e10.0)',str
   IF n_tags(d1) EQ 0 THEN d1=str ELSE d1=[d1,str]
ENDWHILE 
free_lun,lin1

;--
chck=file_search(dir2,'lorentz_test_2_chianti.txt',count=count)
IF count EQ 1 THEN BEGIN
   file2=chck[0]
ENDIF ELSE BEGIN
   print,'% LORENTZ_COMPARE_TEST2: File not found in DIR2. Returning...'
   return,-1
ENDELSE

str1=''
openr,lin2,file2,/get_lun
FOR i=0,3 DO readf,lin2,str1
print,'***File 2 info***'
readf,lin2,str1
print,'  ',str1
;
ltemp=fltarr(51)
readf,lin2,format='(3x,51f10.0)',ltemp
str={z: 0, data: fltarr(51)}
d2=0
WHILE eof(lin2) NE 1 DO BEGIN
   readf,lin2,format='(i3,51e10.0)',str
   IF n_tags(d2) EQ 0 THEN d2=str ELSE d2=[d2,str]
ENDWHILE 
free_lun,lin2


IF n_elements(z) EQ 0 THEN z=1

i=where(d1.z EQ z,ni)
IF ni NE 0 THEN r1=d1[i].data ELSE return,-1

i=where(d1.z EQ z,ni)
IF ni NE 0 THEN r2=d2[i].data ELSE return,-1

z2element,z,elt
p=plot(ltemp,r2,/ylog,color='blue',name=dir2,thick=2, $
       xtitle='Log!d10!n (Temperature/K)', $
       ytitle='Radiative power', $
       title='Z='+trim(z)+', '+trim(elt))
q=plot(/overplot,ltemp,r1,name=dir1)

l=legend(target=[p,q],pos=[0.9,0.8])


return,p

END
