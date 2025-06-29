
;+
; NAME:
;     FF_READ_DATA
;
; PURPOSE:
;     Read the radiative data stored in the MCHF format files produced
;     by Charlotte Froese Fischer and colleagues. The data are
;     returned in an IDL structure. See Programming Notes below on
;     formatting of the MCHF files.
;
; CATEGORY:
;     CHIANTI; data formatting.
;
; CALLING SEQUENCE:
;     FF_READ_DATA, fname, trans
;
; INPUTS:
;     Fname:   The name of the MCHF data file.
;
; OUTPUTS:
;     Trans:   A structure containing the MCHF data. The tags are, e.g.,
;              CONF1           STRING    '3s2.3p'
;              LEV1            STRING    '2P1/2'
;              CONF2           STRING    '3s2.3p'
;              LEV2            STRING    '2P3/2'
;              GF              FLOAT           0.00000
;              AVAL            FLOAT       5.87923e-08
;              E1_CM           FLOAT           0.00000
;              E2_CM           FLOAT           2181.42
;              DE_CM           FLOAT           2181.42
;              CHCKSTR         STRING    '3s2.3p3s2.3p2P1/22P3/2'
;
; OPTIONAL OUTPUTS:
;	Describe optional outputs here.  If the routine doesn't have any, 
;	just delete this section.
;
; EXAMPLE:
;       IDL> ff_read_data, 'ar_6_froese.txt', trans
;
; PROGRAMMING NOTES:
;       The fortran format of the data lines is assumed to be:
;        7x
;        i3  Weight of lower level
;        i3  Weight of upper level
;        a4  Transition type (E1, M1, E2, etc.)
;        f14 Energy (cm-1) of lower level
;        f14 Energy (cm-1) of upper level
;        f14 Energy separation (cm-1)
;        f13 Wavelength in angstroms (not read)
;        e10 Line strength, S (not read)
;        e10 Oscillator strength
;        e11 A-value
;
;       For some transitions there may be more than one entry (e.g.,
;       forbidden transitions with E2 and M1 components). These are
;       added.
;
;       Since the file gives the oscillator strengths, then these are
;       multiplied by the weight of the lower level to give the
;       weighted oscillator strength (gf) for CHIANTI, but only for
;       the E1 transitions (otherwise gf=0).
;
; MODIFICATION HISTORY:
;       Ver.1, 12-Jun-2017, Peter Young
;         Modified from earlier routine, and renamed.
;       Ver.2, 29-Jun-2025, Peter Young
;         When checking for duplicate transitions, I now check the
;         energies of the lower and upper levels. Previously I was
;         using "chckstr" (this didn't work for Ne IV).
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


;-----------------------
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


;-----------------------------------------------
PRO ff_read_data, fname, trans, alt=alt

openr,lun,fname,/get_lun

;
; The following reads the 7 lines of the header
;
str1=''
FOR i=0,6 DO readf,lun,str1

str={conf1: '', lev1: '', conf2: '', lev2: '', gf: 0., aval: 0., $
     e1_cm: 0., e2_cm: 0., de_cm: 0., chckstr: ''}
trans=0

type=''

IF keyword_set(alt) THEN BEGIN
  format='(6x,i3,i4,a4,2f11.0,21x,e10.0,e10.0)'
  line_length=80
ENDIF ELSE BEGIN 
  format='(7x,2i3,a4,2f14.0,37x,e10.0,e11.0)'
  line_length=103
ENDELSE 


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
        ;
        ; If swtch=0 then expect to read the configuration strings.
        ;
        bits=str_sep(str1,' - ')
        IF n_elements(bits) LT 5 AND n_elements(bits) GT 1 THEN BEGIN
          conf1=trim(bits[0])
          conf2=trim(bits[1])
          conf1=process_conf(conf1)
          conf2=process_conf(conf2)
          swtch=1
        ENDIF
       ;
       ; Sometimes multiple term transitions are given under the same
       ; configuration transition. The additional term transitions are
       ; identified by the lack of a '-'. I then process str1 to get
       ; the two terms.
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
        IF strlen(str1) NE line_length THEN BEGIN
          print,'%PROCESS_FF_DATA: warning, string length is not '+trim(line_length)+' characters!'
          print,'                  may be due to tabs in input string'
          print,'         Line : '+trim(count)
        ENDIF 
        reads,str1,format=format, $
             wgt1,wgt2,type,e1,e2,f,aval
        str.conf1=trim(conf1)
        str.conf2=trim(conf2)
        str.lev1=term1+wgt2string(wgt1)
        str.lev2=term2+wgt2string(wgt2)
        str.gf=0.
        str.e1_cm=e1
        str.e2_cm=e2
        str.de_cm=e2-e1
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
          ind=where(str.e1_cm EQ trans.e1_cm AND str.e2_cm EQ trans.e2_cm)
;          ind=where(chckstr EQ trans.chckstr)
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
