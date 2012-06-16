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
; NAME        : CONVERT_TERMS
;     		          
; PURPOSE     : 
;              to convert the transition information into readable formats.
;               
;
; CALLING SEQUENCE:
;
;       IDL>this_design =  convert_terms(l1, l2, result_latex =result_latex)
;
; PROCEDURE: 
; 		This function is used to convert the transition information into
; 		readable formats, including spaces and calulating the J values.
;		The process is quite complex since the notation of the level
;		designations in the CHIANTI files is not standard. 
;		A very useful latex-style output is also created.
;
;    
; INPUTS      : 
;		l1:  the index of the lower level
;               l2:  the index of the upper level
;               
;               the rest of the input info is taken from the COMMON.
;
; OPT. INPUTS : none
;
;               
; OUTPUTS     : 
;	        a string with the transition information. 
;		
;		
; OPT. OUTPUTS:
;	       a string with the transition information in latex format. 
;	
;		
; KEYWORDS    : 
;              RESULT_LATEX: the name of an output string.
;    
;
; CALLS       : 		
;		CONVERT_CONFIG and other SolarSoft routines.
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
;		this_design =  convert_terms(l1, l2, result_latex =result_latex)		
;
;
; CATEGORY    : 
;               spectral synthesis.
;
; PREV. HIST. :
;             parts of the function convert_terms  are derived from the CHIANTI
;             routines.  
;
; WRITTEN     : 
;
;       Giulio Del Zanna (GDZ), 10-Oct-2000 
;	DAMTP  (University of Cambridge, UK) 
;
; MODIFIED    : Version 1, GDZ 10-Oct-2000
;               Version 2, GDZ 10-Oct-2001  
;               Uses standard SolarSoft routines.
;               V. 3, 18-Sep-2002, GDZ 
;                Fixed a bug with the jvalue
;               Ver. 4, Peter Young, 19-Nov-2012
;                 I've taken out the convert_config routine so
;                 that it is a separate routine.
;
;
; VERSION     : 4, 19-Nov-2012
;
;-
;

;-------------------------------------------------------------
;Main routine starts here:
;-------------------------------------------------------------


FUNCTION  convert_terms,  l1, l2, result_latex =result_latex 

COMMON elvlc,l1a,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref

if n_params() lt 2 then begin
   print, 'Error !'
    result_latex =''
   return, ''
ENDIF

;I need to add more checks, etc....


spd=['S','P','D','F','G','H','I','K']
;jvalue=['0','1/2','1','3/2','2','5/2','3','7/2','4','9/2','5',  $
;        '11/2','13/2','15/2']

jvalue=['0','1/2','1','3/2','2','5/2','3','7/2','4','9/2','5',  $
        '11/2','6', '13/2','7', '15/2']

;  get lower level designation
;------------------------------
term0=strtrim(term(l1-1),2)
term1=''
blank=strpos(term0,' ')
while blank gt 0 do begin
   term1=term1+' '+strmid(term0,0,blank)
   term0=strmid(term0,blank,100)
   term0=strtrim(term0,2)
   blank=strpos(term0,' ')
ENDWHILE
;remchar,term1,'^'


jinteger=fix(2.*jj(l1-1))
jstring1=jvalue(jinteger)

spins=strtrim(string(ss(l1-1),'(i2)'),2)

desig1 = spins+spd(ll(l1-1))+jstring1
desig1_latex ='$^'+spins+'$'+spd(ll(l1-1))+'$_{'+jstring1+'}$'


;  get upper level designation
;------------------------------
term0=strtrim(term(l2-1),2)
;
term2=''
blank=strpos(term0,' ')
;           print,term0
while blank gt 0 do begin
   term2=term2+' '+strmid(term0,0,blank)
   term0=strmid(term0,blank,100)
   term0=strtrim(term0,2)
   blank=strpos(term0,' ')
endwhile
;remchar,term2,'^'


jinteger=fix(2.*jj(l2-1))
jstring2=jvalue(jinteger)
;
spins=strtrim(string(ss(l2-1),'(i2)'),2)

desig2=spins+spd(ll(l2-1))+jstring2
desig2_latex= '$^'+spins+'$'+spd(ll(l2-1))+'$_{'+jstring2+'}$'

result_txt = convert_config(term1)+ ' '+desig1+$
  ' - '+ convert_config(term2)+ ' '+desig2

result_latex = convert_config(term1, /latex)+ ' '+desig1_latex+$
  ' - '+ convert_config(term2, /latex)+ ' '+desig2_latex
  
return, result_txt

END


