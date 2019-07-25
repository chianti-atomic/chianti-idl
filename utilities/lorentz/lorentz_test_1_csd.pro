
PRO lorentz_test_1_csd
;+
; NAME:
;     LORENTZ_TEST_1_CSD
;
; PURPOSE:
;     Computes output for Lorentz Test No. 1, charge state
;     distributions.
;
; CATEGORY:
;	Spectral modeling codes; output.
;
; CALLING SEQUENCE:
;       LORENTZ_TEST_1_CSD
;
; INPUTS:
;	None.
;
; OUTPUTS:
;       Creates a text file containing the results.
;
; EXAMPLE:
;       IDL>  lorentz_test_1_csd
;
; MODIFICATION HISTORY:
;       Ver.1, 17-Aug-2016, Peter Young
;-

IF n_elements(ioneq_file) EQ 0 THEN ioneq_file=!ioneq_file

basename='lorentz_test_1_chianti.txt'

temp=[1.0,6.0,46.42]
tempstr=strpad(trim(round(temp)),2,fill='0')
temp=temp*1e6
nt=n_elements(temp)

ioneqfiles=strarr(3)

FOR i=0,nt-1 DO BEGIN
  outname='isothermal_'+tempstr[i]+'MK.ioneq'
  make_ioneq_all,temp[i],outname=outname
  ioneqfiles[i]=outname
ENDFOR


elts=[1,2,6,7,8,10,12,14,16,18,20,26,28]
nelt=n_elements(elts)
elts_txt=strarr(nelt)
FOR i=0,nelt-1 DO BEGIN
  z2element,elts[i],txt,/symbol
  elts_txt[i]=txt
ENDFOR 

openw,lout,basename,/get_lun
printf,lout,'#Column 1: Atomic number of element'
printf,lout,'#Column 2: Ion charge'
printf,lout,'#Column 3: Log10 of ion fraction at temperature 1 MK'
printf,lout,'#Column 4: Log10 of ion fraction at temperature 6 MK'
printf,lout,'#Column 5: Log10 of ion fraction at temperature 46.4 MK (4 keV)'
printf,lout,'#Comment: the minimum ion fraction is forced to be 10^-20'
printf,lout,'#Derived using CHIANTI version '+ch_get_version()
printf,lout,'#File created: '+systime()

FOR i=0,nelt-1 DO BEGIN
;  outname=basename+strlowcase(trim(elts_txt[i]))+'.txt'
  data=fltarr(elts[i]+1,3)
  FOR j=0,nt-1 DO BEGIN
    read_ioneq,ioneqfiles[j],tt,ioneq,ref
    data[*,j]=ioneq[0,elts[i]-1,0:elts[i]] > 1e-20
  ENDFOR 
  FOR k=0,elts[i] DO printf,lout,format='(i3,i3,3f12.3)',elts[i],k,alog10(data[k,*])
ENDFOR

free_lun,lout


END
