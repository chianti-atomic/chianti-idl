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
;
;
; VERSION     : 3, 18-Sep-2002
;
;-
;

FUNCTION convert_config, config, latex=latex


;replace '.' with spaces.
pp = repstr(config, '.', ' ')

;insert a space where parentheses are:
pp = repstr(pp, '(', ' (')
pp = repstr(pp, ')', ') ')

;deal with the stuff in parentheses.  Watch that there might be more than one
;(or two) parentheses !!!!

index1 = str_index(pp, '(')
index2 = str_index(pp, ')')

IF n_elements(index1) NE n_elements(index2) THEN BEGIN 
   print,  'wrong number of parentheses !'
;return the original input...
   return, config
END

IF min(index1) GE 0 AND min(index2) GE 0 THEN BEGIN 

;put lowercase what is outside the   parentheses. there is at least one "(". 

   part = strlowcase(strmid(pp, 0, index1(0)+1)) ;get the '(' too.

;do the first parenthesis.
;allow only things like (2D*) or (2P) to become ($^2$D*)

   IF  (index2(0) -index1(0)-1 EQ 2 OR index2(0) -index1(0)-1 EQ 3) AND $
     valid_num(strmid(pp, index1(0)+1, 1),/integer)  THEN BEGIN 
      
;if the first character after the '(' is an integer then proceed
; and uppercase what is inside.

      IF keyword_set(latex) THEN $ 
        part = part+ '$^'+strmid(pp, index1(0)+1, 1)+'$'+$
        STRUPCASE(strmid(pp, index1(0)+2, index2(0)-index1(0)-1)) ELSE $
        part = part+strmid(pp, index1(0)+1, 1)+$
        STRUPCASE(strmid(pp, index1(0)+2, index2(0)-index1(0)-1))

   ENDIF ELSE  part = part+strmid(pp, index1(0)+1, index2(0)-index1(0))


;now I have closed the first ()
   IF n_elements(index1) EQ 1 THEN $
     pp = part+strlowcase(strmid(pp,  index2(0)+1, 100)) ELSE BEGIN 


;in this case there are more than one ()

      FOR p=0, n_elements(index1)-2 DO BEGIN 
         
         part = part+strlowcase(strmid(pp,  index2(p)+1, index1(p+1)-index2(p)))
         
         IF  index2(p+1) -index1(p+1)-1 LE 3 THEN BEGIN 
            
            IF valid_num(strmid(pp, index1(p+1)+1, 1),/integer)  THEN BEGIN 

;if the first character after the '(' is an integer then proceed
               IF keyword_set(latex) THEN $ 
                 part = part+ '$^'+strmid(pp, index1(p+1)+1, 1)+'$'+$
                 STRUPCASE(strmid(pp, index1(p+1)+2, index2(p+1)-index1(p+1)-1)) ELSE $
                 part = part+strmid(pp, index1(p+1)+1, 1)+$
                 STRUPCASE(strmid(pp, index1(p+1)+2, index2(p+1)-index1(p+1)-1))
               
            ENDIF 
         ENDIF ELSE   part = part+strmid(pp, index1(p+1)+1, index2(p+1)-index1(p+1)+1)


      ENDFOR 
;do the last bit 

      pp = part+strlowcase(strmid(pp,  index2(p)+1, 100))


   ENDELSE

;if there are no parentheses,then lowercase everything:
ENDIF     ELSE  pp = strlowcase(pp)  


spdf = ['s', 'p', 'd', 'f', 'g','h','i','k'] 

FOR p=0, 4 DO BEGIN 

;there should be only one occurence of each s,p,d,f,g (could add a check in the
;future....... 

   index1 = strpos(pp, spdf(p) )

   IF index1 NE -1 THEN BEGIN 

      IF valid_num(strmid(pp, index1+1, 1),/integer) THEN BEGIN

;if the next is a number, then proceed, allowing various cases:

;if the next is a number, ;like 2s22p, then insert the dollars and add  a space
         IF  valid_num(strmid(pp, index1+2, 1),/integer)   THEN BEGIN  

            IF keyword_set(latex) THEN $ 
              pp = strmid(pp, 0, index1+1)+'$^'+ $
              strmid(pp, index1+1, 1)+'$'+$
              ' '+ strmid(pp, index1+2, 100)   ELSE $
              pp = strmid(pp, 0, index1+1)+ $
              strmid(pp, index1+1, 1)+$
              ' '+ strmid(pp, index1+2, 100)  


         ENDIF ELSE  IF  strmid(pp, index1+2, 1) EQ ' ' OR  $
           strmid(pp, index1+2, 1) EQ '' THEN BEGIN 

;if there is a space already, then produce 2s$^2$ 2p if no space (at the end) do
;it anyway.
            IF keyword_set(latex) THEN $
              pp = strmid(pp, 0, index1+1)+'$^'+ $
              strmid(pp, index1+1, 1)+'$'+strmid(pp, index1+2, 100) ELSE $
              pp = strmid(pp, 0, index1+1)+ $
              strmid(pp, index1+1, 1)+strmid(pp, index1+2, 100)

         ENDIF ELSE $
;if string e.g. 2s2p then add a space
         pp = strmid(pp, 0, index1+1)+' '+ $
           strmid(pp, index1+1, 100)

      ENDIF 
   ENDIF                        ;found an spdf....

ENDFOR

return, strtrim(pp, 2)

END

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


