

PRO nist_energies, infile, outfile, add_empty_theo=add_empty_theo

;+
; NAME:
;      NIST_ENERGIES
;
; PURPOSE:
;      Convert a NIST format energy level table into the CHIANTI elvlc
;      format. 
;
; CATEGORY:
;      CHIANTI; NIST; formatting.
;
; CALLING SEQUENCE:
;      NIST_ENERGIES, INFILE, OUTFILE
;
; INPUTS:
;      Infile:   The name of the NIST file.
;      Outfile:  The name of the CHIANTI format output file.
;
; OUTPUTS:
;      Creates a new file in the user's working directory.
;
; RESTRICTIONS:
;      Sometimes the term is given as, e.g., 'aD' which currently gives a
;      format warning but does not crash the routine. Usually such levels
;      won't be in the CHIANTI model, so not a worry.
;
;      The first line in the NIST file should give the energy of the
;      ground level.
;
; EXAMPLE:
;      IDL> nist_energies, 'ar_2_nist_energies.txt','ar_2.elvlc'
;
; MODIFICATION HISTORY:
;      Ver.1, 22-May-2017, Peter Young
;        Updated routine for the new format ELVLC files.
;      Ver.2 23-Jul-2020, Peter Young
;        Now deals with level labels correctly. For example, "a 6D". 
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
        label=''
      END
      4: BEGIN
        chck=str_sep(term,' ')
        IF n_elements(chck) GT 1 THEN BEGIN
          label=trim(chck[0])
          term=chck[1]
          ll=strmid(term,1,1)
          j=where(ll EQ orb)
          lln=j
          ss=fix(strmid(term,0,1))
        ENDIF 
      END
      
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
      printf,lout,format='(i7,a30,a5,i5,a5,f5.1,f15.3)', $
             i,config,label,ss,ll,jj,encm
      ;; printf,lout,format='(i3,a6,a15,2i3,a2,f5.1,i3,f15.3,f15.6,f15.3,f15.6)', $
      ;;        i,'x',config,ss,lln,ll,jj,2.*jj+1,encm,encm/109737.32, 0.,0.
    ENDIF ELSE BEGIN
      printf,lout,format='(i7,a30,a5,i5,a5,f5.1,2f15.3)', $
             i,config,label,ss,ll,jj,encm,0.
      ;; printf,lout,format='(i3,a6,a15,2i3,a2,f5.1,i3,f15.3,f15.6)', $
      ;;        i,'x',config,ss,lln,ll,jj,2.*jj+1,encm,encm/109737.32
    ENDELSE
    i=i+1
  ENDIF
ENDWHILE

printf,lout,' -1'
printf,lout,'%file:  '+trim(outfile)

free_lun,lin,lout

END
