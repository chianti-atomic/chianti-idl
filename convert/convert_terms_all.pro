;+
; PROJECT     :  CHIANTI
;
;       CHIANTI   http://wwwsolar.nrl.navy.mil/chianti.html
;
;       CHIANTI is an atomic database package for the calculation of
;       continuum and emission line spectra from astrophysical plasmas. It is a 
;       collaborative project involving the Naval Research Laboratory
;       (Washington DC, USA), the Arcetri Observatory (Firenze, Italy), and the
;       Cambridge University (United Kingdom).
; 
;                   
; NAME        : CONVERT_TERMS_ALL
;     		          
; PURPOSE     : 
;              to convert level information into readable formats for all the levels
;              of an ion.
;               
;
; CALLING SEQUENCE:
;
;       IDL>convert_terms_all, res_ascii, res_latex
;
; PROCEDURE: 
; 		This function is used to convert the level information into
; 		readable formats, including spaces and calculating the J values.
;		A very useful latex-style output is also created.
;               It has been adapted to manage all the levels of a given ion
;               at once from the original routine CONVERT_TERMS.PRO.
;
;    
; INPUTS      : all the input info is taken from the COMMON.
;
; OPT. INPUTS : none
;
;               
; OUTPUTS     : 
;	       RES_ASCII: a string array with the transition information in ASCII format. 
;		
;              RES_LATEX: a string array with the transition information in LATEX format.
;    
;
; CALLS       : 		
;		CONVERT_CONFIG (included here) and other SolarSoft routines.
;
;		
; COMMON BLOCKS: elvlc,l1a,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref
; 		from here, all the information needed is taken.
;		
;
; RESTRICTIONS: It will not convert every case. 
;
;               
; SIDE EFFECTS: None known yet.
;               
;
; EXAMPLES    : 
;		convert_terms_all, res_ascii, res_latex
;
;
; CATEGORY    : 
;               spectral synthesis.
;
; PREV. HIST. :
;             convert_terms_all has been derived from the existing CHIANTI
;             routine CONVERT_TERMS.
;
; WRITTEN     : 
;             Enrico Landi, 11-Apr-2005
;	      Naval Research Laboratory
;
; MODIFIED    : Version 1, EL 11-Apr-2005
;
;               Version 2, Peter Young, 27-Jul-2012
;               Made the jvalue array larger to accommodate large
;               data-sets. 
;
;               Version 3, Peter Young, 14-Jan-2013
;               convert_config.pro is now a separate routine,
;               so I've removed the subroutine of the same name
;               from this routine.
;
;
; VERSION     : 3, 14-Jan-2013
;
;-
;


PRO  convert_terms_all, res_ascii, res_latex

COMMON elvlc,l1a,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref

IF n_params() lt 1 THEN BEGIN
   print, 'Error !'
   result =''
ENDIF

;I need to add more checks, etc....


spd=['S','P','D','F','G','H','I','K']
;jvalue=['0','1/2','1','3/2','2','5/2','3','7/2','4','9/2','5',  $
;        '11/2','13/2','15/2']

;
; Added the code below on 27-Jul-2012, PRY.
; If you ever need to add more terms, just change the value of nj.
;
nj=100
jvalue=strarr(nj)
ind=indgen(nj/2)
jvalue[2*ind]=trim(ind)
jvalue[2*ind+1]=trim(2*ind+1)+'/2'

;; jvalue=['0','1/2','1','3/2','2','5/2','3','7/2','4','9/2','5',  $
;;         '11/2','6', '13/2','7', '15/2', '8', ']

nlevels=n_elements(term)
result_txt=strarr(nlevels)
result_latex=strarr(nlevels)

FOR i = 0,nlevels-1 DO BEGIN        ; loop for each level

;  get term designation
;------------------------------
    term0=strtrim(term(i),2)
    term1=''
    blank=strpos(term0,' ')
    WHILE blank gt 0 DO BEGIN
       term1=term1+' '+strmid(term0,0,blank)
       term0=strmid(term0,blank,100)
       term0=strtrim(term0,2)
       blank=strpos(term0,' ')
    ENDWHILE
;remchar,term1,'^'


    jinteger=fix(2.*jj(i))
    jstring1=jvalue(jinteger)

    spins=strtrim(string(ss(i),'(i2)'),2)

    desig1 = spins+spd(ll(i))+jstring1
    desig1_latex ='$^'+spins+'$'+spd(ll(i))+'$_{'+jstring1+'}$'


;  get configuration designation
;------------------------------------------

desig2 = convert_config(term1)+' '
desig2_latex = convert_config(term1, /latex)+' '

    result_txt(i) = desig2 + desig1
    
    result_latex(i) = desig2_latex + desig1_latex
  
ENDFOR

res_ascii=result_txt
res_latex=result_latex

END
