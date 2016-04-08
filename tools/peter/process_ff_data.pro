;+
; Froese Fischer has a website:
;
;  http://www.vuse.vanderbilt.edu/~cff/mchf_collection.html
;
; where her data is stored in a uniform format. Unfortunately it's difficult
; to process this data into the CHIANTI format in an easy way, so this
; routine tries to make the process easier.
;
; HISTORY
;
;    Ver.1, 17-May-2005
;      modified process_conf to deal with configurations of O I.
;    Ver.2, 26-Nov-2007
;      corrected bug related to gf-values
;    Ver.3, 28-Nov-2012
;      I found a case where a tab had been wrongly inserted into the
;      data file. I thus now check the length of the data strings to
;      make sure they are 103 characters.
;-



FUNCTION process_conf, conf
;+
; This routine processes the various formats FF has for configurations.
; Some examples are
; 2p(3)4S3_4S.3s -> 2p3.3s
; 2p(4)3P2       -> 2p4
; 2s.2p_3P.3s    -> 2s.2p.3s
;-

pos1=strpos(conf,'_')
IF pos1 GE 0 THEN BEGIN
  pos2=strpos(conf,'.',pos1)
  sl=strlen(conf)
  str1=strmid(conf,0,pos1)
  str2=strmid(conf,pos2,sl-1)
  conf=str1+str2
ENDIF


bits=str_sep(conf,')')
IF n_elements(bits) EQ 1 THEN return,conf

pos1=strpos(conf,')')
pos2=strpos(conf,'.',pos1)

IF pos2 LT pos1 THEN BEGIN
  conf=bits[0]+')'
ENDIF ELSE BEGIN
  midbit=strmid(conf,pos1,pos2-pos1)
  conf=str_replace(conf,midbit,'')
ENDELSE


; chck=bits[1]
; chck=strmid(chck,0,1)
; IF chck NE '.' THEN conf=bits[0]+')'

conf=str_replace(conf,'(','')
conf=str_replace(conf,')','')

return,conf

END


FUNCTION wgt2string, wgt

wgt=fix(wgt)
IF wgt NE 2*(wgt/2) THEN BEGIN
  str1=trim((wgt-1)/2)
  return,str1
ENDIF ELSE BEGIN
  str1=trim(wgt-1)+'/2'
  return,str1
ENDELSE

END


PRO process_ff_data, fname, trans, levstr=levstr

openr,lun,fname,/get_lun

;
; The following reads the 7 lines of the header
;
str1=''
FOR i=0,6 DO readf,lun,str1

str={conf1: '', lev1: '', conf2: '', lev2: '', gf: 0., aval: 0., $
    chckstr: ''}
trans=0

type=''

count=7
swtch=0
WHILE eof(lun) NE 1 DO BEGIN
  readf,lun,str1
  count=count+1
  IF trim(str1) EQ '' THEN BEGIN
    swtch=0
  ENDIF ELSE BEGIN
    CASE swtch OF 
      0: BEGIN
        bits=str_sep(str1,' - ')
        IF n_elements(bits) LT 5 AND n_elements(bits) GT 1 THEN BEGIN
          conf1=trim(bits[0])
          conf2=trim(bits[1])
          conf1=process_conf(conf1)
          conf2=process_conf(conf2)
          swtch=1
        ENDIF
       ;
       ; the section of code below had to be introduced for O I
       ; where, if the configurations has to changed from the previous
       ; entry, then FF skipped writing them
       ;
        IF n_elements(bits) EQ 1 THEN BEGIN
          str1=trim(str1)
          str1=strcompress(str1)
          bits=str_sep(strcompress(str1),' ')
          term1=str_replace(bits[0],'*','')
          term2=str_replace(bits[1],'*','')
          swtch=2
        ENDIF
      END

      1: BEGIN
        str1=trim(str1)
        str1=strcompress(str1)
        bits=str_sep(strcompress(str1),' ')
        term1=str_replace(bits[0],'*','')
        term2=str_replace(bits[1],'*','')
        swtch=2
      END

      2: BEGIN
        IF strlen(str1) NE 103 THEN BEGIN
          print,'%PROCESS_FF_DATA: warning, string length is not 103 characters!'
          print,'                  may be due to tabs in input string'
          print,'         Line : '+trim(count)
        ENDIF 
        reads,str1,format='(7x,2i3,a4,65x,e10.0,e11.0)', $
             wgt1,wgt2,type,f,aval
        str.conf1=trim(conf1)
        str.conf2=trim(conf2)
        str.lev1=term1+wgt2string(wgt1)
        str.lev2=term2+wgt2string(wgt2)
        str.gf=0.
        chckstr=str.conf1+str.conf2+str.lev1+str.lev2
        str.chckstr=chckstr
        IF trim(type) EQ 'E1' THEN str.gf=f*wgt1
        str.aval=aval
        IF n_tags(trans) EQ 0 THEN BEGIN
          trans=str        
        ENDIF ELSE BEGIN
         ;
         ; for some transitions, two components to the A-value may be given
         ; (E1+M2 or E2+M1). The following checks for this and adds together
         ; the two A-values. For the gf value the maximum of the tabulated
         ; gf values is given
         ;
          ind=where(chckstr EQ trans.chckstr)
          IF ind[0] EQ -1 THEN BEGIN
            trans=[trans,str]
          ENDIF ELSE BEGIN
            trans[ind[0]].gf=max([trans[ind[0]].gf,str.gf])
            trans[ind[0]].aval=trans[ind[0]].aval+str.aval
          ENDELSE
        ENDELSE
      END

    ENDCASE

  ENDELSE
ENDWHILE


free_lun,lun

END
