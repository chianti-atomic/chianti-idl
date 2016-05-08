
FUNCTION convert_config, config, latex=latex, parity=parity


;+
; NAME
;
;     CONVERT_CONFIG()
;
; PROJECT
;
;     CHIANTI
;
; EXPLANATION
;
;     Takes the configuration descriptor from the CHIANTI energy level
;     and converts it to a standard format. In addition, by giving the
;     /LATEX keyword it can be converted to a latex format.
;
; INPUTS
;
;     CONFIG    A string containing the configuration
;               descriptor. E.g., '3s2.3p3'.
;
; KEYWORDS
;
;     LATEX     If set, then input string is converted to a latex
;               format. 
;
; OPTIONAL OUTPUTS
;
;     PARITY    The parity of the configuration. Either 1 for odd or 0
;               for even.
;
; OUTPUTS
;
;     A string containing the new configuration format.
;
; EXAMPLES
;
;     IDL> print,convert_config('3s2.3p3')
;     3s2 3p3
;
;     IDL> print,convert_config('3s2.3p3',/latex)
;     3s$^2$ 3p$^3$
;
; HISTORY
;
;     Ver. 1, Peter Young, 19-Nov-2012
;        This routine was originally embedded within
;        CONVERT_TERMS.pro. I have extracted it into its own routine.
;        No change has been made to the code.
;     Ver. 2, Peter Young, 20-Dec-2012
;        The method of re-formatting the string has been modified, and
;        the parity can now be output (PARITY= keyword).
;
;     (Code originally written by Giulio Del Zanna, 18-Sep-2002.)
;-

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
n_spdf=n_elements(spdf)

parity=0

;
; Notes
;   1. Inner shell excitation can lead to the case of e.g., '2p5.3p6',
;      in which there are two 'p's in the description. Since strpos
;      can only find one instance of a match, then I've
;      introduced the for loop on i which allows two instances of 'p'
;      to be found. 
;   2. The method will fail if the configuration '2s2 10p' is written
;      as '2s210p', but such a configuration should never be written
;      this way in the first place (!).
;   3. The routine will work for configurations such as '3d10 4s'.
;
save_index1=-1
FOR p=0,n_spdf-1 DO BEGIN 

  FOR i=0,1 DO BEGIN 

  index1 = strpos(pp, spdf(p) , reverse_search=i)

  IF index1 NE -1 AND index1 NE save_index1 THEN BEGIN 

    IF i EQ 0 THEN save_index1=index1 ELSE save_index1=-1

    s1=strmid(pp,index1+1,1)
    s2=strmid(pp,index1+2,1)
    s3=strmid(pp,index1+3,1)

    CASE 1 OF 
     ;
     ; E.g., '3s2 ' or '3s2'
     ;
      valid_num(s1,/integer) AND (s2 EQ ' ' OR s2 EQ ''): BEGIN
        occup_str=s1
        occup=fix(s1)
        j=index1+2
      END 
     ;
     ; E.g., '3s '
     ;
      NOT valid_num(s1,/integer): BEGIN 
        occup_str=''
        occup=1
        j=index1+1
      END 
     ;
     ; E.g., '3d10 4s'
     ;
      valid_num(s1,/integer) AND valid_num(s2,/integer) AND (s3 EQ ' ' OR s3 EQ ''): BEGIN
        occup_str=s1+s2
        occup=fix(s1+s2)
        j=index1+3
      END 
     ;
     ; E.g., '3d104s'
     ;
      valid_num(s1,/integer) AND valid_num(s2,/integer) AND valid_num(s3,/integer): BEGIN
        occup_str=s1+s2
        occup=fix(s1+s2)
        j=index1+2
      END 
     ;
     ; E.g., '3s22p'
     ;
      valid_num(s1,/integer) AND valid_num(s2,/integer) AND NOT valid_num(s3,/integer): BEGIN
        occup_str=s1
        occup=fix(s1)
        j=index1+2
      END 
     ;
     ; E.g., '3s2p'
     ;
      valid_num(s1,/integer) AND NOT valid_num(s2,/integer) AND s2 NE ' ': BEGIN
        occup_str=''
        occup=1
        j=index1+1
      END 
     ;
     ; E.g., '4s'
     ;
      ELSE: BEGIN
        occup_str=''
        occup=1
        j=index1+1
      END 

    ENDCASE 

    IF (p/2)*2 NE p THEN parity=parity+occup

    IF occup EQ 1 THEN occup_str=''
    IF keyword_set(latex) THEN BEGIN
      IF occup_str NE '' THEN add_str='$^{'+occup_str+'}$' ELSE add_str=''
      pp=strmid(pp,0,index1+1)+add_str+strmid(pp,j,100)
    ENDIF ELSE BEGIN
      pp=strmid(pp,0,index1+1)+occup_str+strmid(pp,j,100)
    ENDELSE 
  ENDIF 
ENDFOR 
ENDFOR


IF (parity/2)*2 EQ parity THEN parity=0 ELSE parity=1

return, strtrim(pp, 2)

END
