;+
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. See www.chiantidatabase.org
;                   
; Name        : GET_ABUND_DATA
;     		          
; Purpose     : selects and reads a CHIANTI abundance dataset
;
; Category    : diagnostic analysis
;               
;    
; Inputs      : none
;              
; Outputs     : the abundance file and (within a COMMON) the values.
;
;
; Calls       : read_abund,  ch_get_file, change_abund 
;
;		
; Common      : elements
;
; Restrictions: the routine works within the CHIANTI software.
;
;
; Written     : Giulio Del Zanna (GDZ)
;               version 1, 11 Mar 2014
;
; Modified:    v.2, GDZ, 18-Sept-2015  changed the default
;               photospheric abundance file to avoid crash with v.8
;
; VERSION     :  2, 18-Sept-2015
;
;-


PRO get_abund_data, abund_name


COMMON elements, abund,abund_ref, ioneq,ioneq_logt,ioneq_ref


;get the abundance file
;----------------------

;read the  photospheric values:

abund_name =!xuvtop+'/abundance/sun_photospheric_2009_asplund.abund'
read_abund,abund_name,abund_phot,abund_ref


abund_name = ch_get_file(path=!xuvtop+'/abundance', filter='*.abund',  title=' Select an abundance file ')

;read the abundance file and print it
;------------------------------------
; the abundances are converted from  
; logaritmic values....abund(g)=10.^(abund(g)-12.)
;----------------------------------------------------------
read_abund,abund_name,abund,abund_ref
;
gz=where(abund gt 0.)
nz=n_elements(gz)

IF NOT keyword_set(quiet) THEN BEGIN 

   print,' '
   print,' abundances (log values) , and  log differences from photospheric (Asplund et al. 2009): '
   print,' '

   for kz=0,nz-1 do begin
      jz=gz(kz)
      z2element,jz+1,ele
      print,jz+1,strpad(ele,12,/after), alog10(abund(jz))+12. ,$
        alog10(abund(jz))-alog10(abund_phot(jz)), $
        format='(i5,2x,a12,f5.2,1x,f5.2)'
   endfor

ENDIF

;this is necessary to write down the name of the ab file used for the dem:
;--------------------------------------------------------------------------

yes_no,'Do you want to change the abundances?  ',yesno, 'N'
if yesno  THEN $
  change_abund ,  abund_name,abund,abund_ref


END  
