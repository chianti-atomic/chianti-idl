;+
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. See www.chiantidatabase.org
;                   
; Name        : CHANGE_ABUND
;     		          
; Purpose     : to change the abundance values editing an abundance file.
;
; Category    : diagnostic analysis
;               
;    
; Inputs      :  abund_name,abund,abund_ref
;              
; Outputs     : the modified abund_name,abund,abund_ref
;
; Calls       : read_abund
;		
; Common      : none
;
; Restrictions: the routine works within the CHIANTI software.
;
;
; Written     : Giulio Del Zanna (GDZ)
;               version 1, 11 Mar 2014
;
; VERSION     :  1, 11 Mar 2014
;-

PRO change_abund ,  abund_name,abund,abund_ref


IF n_params() LT 3 THEN message, 'error '


gz=where(abund gt 0.)
nz=n_elements(gz)



print,' abundances'
for kz=0,nz-1 do begin
   jz=gz(kz)
   z2element,jz+1,ele
   print,jz+1,strpad(ele,12,/after),abund(jz),format='(i5,2x,a12,e10.2)'
endfor

new_abund_name=''
read,'Enter the new core file name (a suffix ".abund" will be added) : ',$
  new_abund_name
new_abund_name=new_abund_name+'.abund'
spawn,'cp '+abund_name+' '+new_abund_name
print,'The file " ',new_abund_name,'"  has been created in the working directory !'
spawn,'chmod u+w '+new_abund_name
editor='emacs'
read,'Which editor do you want to use ? (e.g. emacs) ',editor
spawn,editor+' '+new_abund_name

;overwrite the string definition of the abundance file:
;-----------------------------------------------------
abund_name=new_abund_name

	print,' '
	print,' The new abundances are: ' 
	print,' '

; the abundances are converted from
; logaritmic values....abund(g)=10.^(abund(g)-12.)
;
; if you want the converse: 
; abund= alog10(abund) +12.
;----------------------------------------------------------
read_abund,abund_name,abund,abund_ref
;
gz=where(abund gt 0.)
nz=n_elements(gz)

for kz=0,nz-1 do begin
   jz=gz(kz)
   z2element,jz+1,ele
   print,jz+1,strpad(ele,12,/after),abund(jz),format='(i5,2x,a12,e10.2)'
endfor

print,' '

return

END
