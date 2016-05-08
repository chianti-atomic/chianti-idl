

PRO nist_energies, infile, outfile, add_empty_theo=add_empty_theo

;+
; NAME
;
;     NIST_ENERGIES
;
; EXPLANATION
;
;     Converts NIST energy table into CHIANTI format. It requires that the
;     NIST table be output in ASCII format.
;
; INPUTS
;
;     INFILE  A file containing the NIST energies obtained directly from
;             NIST. (The columns should be separated by '|'.)
;
;     OUTFILE CHIANTI format energy file.
;
; KEYWORDS
;
;     ADD_EMPTY_THEO  Add 'empty' theoretical energy columns. I.e.,
;                     all energies set to zero.
;
; RESTRICTIONS
;
;     Sometimes the term is given as, e.g., 'aD' which currently gives a
;     format warning but does not crash the routine. Usually such levels
;     won't be in the CHIANTI model, so not a worry.
;
; HISTORY
;
;     Beta 1, 17-Mar-2005, Peter Young
;     Beta 2, 24-May-2005, Peter Young
;        now deals with even ions (fractional J values)
;        also I print an 'x' for the configuration index
;     Beta 3, 21-Jan-2009, Peter Young
;        added /add_empty_theo
;-

IF n_params() LT 1 THEN BEGIN
  print,'Use: IDL> nist_energies, infile, outfile [ /add_empty_theo ]'
  return
ENDIF

openr,lin,infile,/get_lun
openw,lout,outfile,/get_lun

orb=['S','P','D','F','G','H','I']

str1=''
i=1
WHILE eof(lin) NE 1 DO BEGIN
  readf,lin,str1
  bits=str_sep(str1,'|')
 ;
  conf=trim(bits[0])
  IF conf NE '' THEN BEGIN
    config=conf
    config=str_replace(config,'*','')
    config=str_replace(config,'.(','(')
    config=str_replace(config,'.',' ')
  ENDIF
 ;
  te=trim(bits[1])
  IF te NE '' THEN BEGIN
    term=te
    term=str_replace(term,'*','')
    CASE strlen(term) OF
      2: BEGIN
        ll=strmid(term,1,1)
        j=where(ll EQ orb)
        lln=j
        ss=fix(strmid(term,0,1))
      END
      3:
      ELSE: print,'** Problem with term string'
    ENDCASE
  ENDIF
 ;
  jj=trim(bits[2])
  jj_bits=str_sep(jj,'/')
  IF n_elements(jj_bits) EQ 1 THEN BEGIN
    jj=fix(jj_bits[0])
  ENDIF ELSE BEGIN
    jj=float(jj_bits[0])/float(jj_bits[1])
  ENDELSE
 ;
  encm=trim(bits[3])
  IF encm NE '' THEN BEGIN
    IF keyword_set(add_empty_theo) THEN BEGIN
      printf,lout,format='(i3,a6,a15,2i3,a2,f5.1,i3,f15.3,f15.6,f15.3,f15.6)', $
             i,'x',config,ss,lln,ll,jj,2.*jj+1,encm,encm/109737.32, 0.,0.
    ENDIF ELSE BEGIN
      printf,lout,format='(i3,a6,a15,2i3,a2,f5.1,i3,f15.3,f15.6)', $
             i,'x',config,ss,lln,ll,jj,2.*jj+1,encm,encm/109737.32
    ENDELSE
    i=i+1
  ENDIF
ENDWHILE

free_lun,lin,lout

END
