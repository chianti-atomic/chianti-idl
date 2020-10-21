
PRO lorentz_compare_test7, dir1, dir2, d1=d1, d2=d2


;+
; NAME:
;     LORENTZ_COMPARE_TEST7
;
; PURPOSE:
;     Compare results of Lorentz Test 7 between different versions of
;     CHIANTI.
;
; CATEGORY:
;     CHIANTI; Lorentz tests.
;
; CALLING SEQUENCE:
;     LORENTZ_COMPARE_TEST7, Dir1, Dir2
;
; INPUTS:
;     Dir1:   Directory containing the results from Lorentz Test 7.
;     Dir2:   Directory containing the results from a different run of
;             Lorentz Test 7.
;
; OUTPUTS:
;     For each emission line, each temperature and each density, the
;     routine compares the set of rates between the DIR1 and DIR2
;     files. If there is a difference of more than 10%, then a message
;     is printed to the IDL window. The ion name, wavelength,
;     temperature and density are printed, along with the min and max
;     of the ratios.
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


lorentz_test_7_levelpop,/only_list,linelist=linelist

nlines=n_elements(linelist)

chck=file_search(dir1,'lorentz_test_7_chianti.txt',count=count)
IF count EQ 1 THEN BEGIN
   file1=chck[0]
ENDIF ELSE BEGIN
   print,'% LORENTZ_COMPARE_TEST7: File not found in DIR1. Returning...'
   return
ENDELSE


str1=''
openr,lin1,file1,/get_lun
FOR i=0,15 DO readf,lin1,str1
;
ldens=[0.,12.]
temp=[1.,6.,46.4]*1e6
ltemp=alog10(temp)
;
data_line=dblarr(10)
str={index: 0, wvl: 0., ion: '', data: dblarr(3,2,10) }
d1=replicate(str,17)
;
FOR i=0,nlines-1 DO BEGIN
   d1[i].index=i+1
   d1[i].ion=linelist[i].ion
   d1[i].wvl=linelist[i].wvl
   FOR j=0,5 DO BEGIN
      readf,lin1,format='(23x,10e10.0)',data_line
      d1[i].data[j/2,j MOD 2,*]=data_line
   ENDFOR
ENDFOR
free_lun,lin1

;--
chck=file_search(dir2,'lorentz_test_7_chianti.txt',count=count)
IF count EQ 1 THEN BEGIN
   file2=chck[0]
ENDIF ELSE BEGIN
   print,'% LORENTZ_COMPARE_TEST7: File not found in DIR2. Returning...'
   return
ENDELSE


str1=''
openr,lin2,file2,/get_lun
FOR i=0,15 DO readf,lin2,str1
;
ldens=[0.,12.]
temp=[1.,6.,46.4]*1e6
ltemp=alog10(temp)
;
data_line=dblarr(10)
str={index: 0, wvl: 0., ion: '', data: dblarr(3,2,10) }
d2=replicate(str,17)
;
FOR i=0,nlines-1 DO BEGIN
   d2[i].index=i+1
   d2[i].ion=linelist[i].ion
   d2[i].wvl=linelist[i].wvl
   FOR j=0,5 DO BEGIN
      readf,lin2,format='(23x,10e10.0)',data_line
      d2[i].data[j/2,j MOD 2,*]=data_line
   ENDFOR
ENDFOR 
free_lun,lin2

;--
FOR i=0,nlines-1 DO BEGIN
   FOR j=0,2 DO BEGIN
      FOR k=0,1 DO BEGIN
         line1=d1[i].data[j,k,*]
         line2=d2[i].data[j,k,*]
         ii=where(line1 NE 0.,nii)
         IF nii NE 0 THEN BEGIN
            ratio=line2[ii]/line1[ii]
            IF max(ratio) GT 1.1 OR min(ratio) LT 0.9 THEN BEGIN
               print,format='(i2,2x,a5,f8.3,2f7.2,2f10.3)', $
                     d1[i].index,d1[i].ion,d1[i].wvl,ltemp[j],ldens[k],min(ratio),max(ratio)
            ENDIF 
         ENDIF 
      ENDFOR
   ENDFOR
ENDFOR




END
