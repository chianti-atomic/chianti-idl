
;+
; NAME:
;    CH_ADF04_CONVERT
;
; PURPOSE:
;    Converts an ADAS adf04 file to the CHIANTI file formats,
;    producing .elvlc, .wgfa and .ups files.
;
; CATEGORY:
;    CHIANTI; format.
;
; CALLING SEQUENCE:
;    CH_ADF04_CONVERT, Filename
;
; INPUTS:
;    Filename:  The name of an ADF04 file.
;
; OUTPUTS:
;    Creates the following files:
;      [ion].ups_adf04
;      [ion].wgfa_adf04
;      [ion].elvlc_adf04
;    where [ion] is the ion name in CHIANTI format (e.g., 'o_6' for O
;    VI). The ion name is automatically taken from the adf04 file.
;
; PROGRAMMING NOTES:
;    The use of the high temperature limit point in the data files
;    seems to vary, so I summarize here.
;
;    For modern adf04 files (since 2010?), there should be a high
;    temperature limit for all transitions, both allowed and
;    forbidden. If a transition is forbidden, then the limit point is
;    given as a *positive* number while for allowed transitions they
;    are given as *negative* numbers.
;
;    For older adf04 files, the limit point was given as
;    zero. I've also seen examples (e.g., Fe V) where the
;    A-value was given as 1e-30.
;
;    There's at least one intermediate case (Ar VI) for which
;    the limit point is non-zero and positive for all transitions. For
;    this case I have code that computes the gf value from the A-value
;    and from the limit point. If they're the same, then I
;    assume a type 1 transition, otherwise it's a type 2 transitions.
;
;    If a limit value of 1e-30 is found, then this implies a limit
;    value of zero. A limit of 0 actually implies the limit is
;    undefined.
;
;    -----
;    If an energy separation is zero, then it is forced to be 1 cm^-1
;    when computing the wavelengths in the .wgfa file.
;
; RESTRICTIONS:
;    Old ADAS files do not have the high temperature limit point and
;    this routine does not work for these. 
;
; MODIFICATION HISTORY:
;    Ver.1, 9-Feb-2017, Peter Young
;      Modified from chianti_adf04_v2 in order to write the new format
;      elvlc and ups files.
;    Ver.2, 14-Jun-2017, Peter Young
;      Changed how gf value and limit point are computed (see
;      Programming Notes above).
;    Ver.3, 15-Jun-20107, Peter Young
;      Renamed routine again and tidied up header and code.
;-


;----------------
FUNCTION conv_num, str
;
; The input is assumed to be an eight character string of the form 
; ' 8.01+02' or '-7.96-06'.
;
str2=trim(str)

a=strmid(str,0,5)
b=strmid(str,5,3)

return,float(a)*10.^(float(b))

;; bits=str_sep(str,'+')
;; IF n_elements(bits) EQ 2 THEN return,float(bits[0])*10.^(float(bits[1]))

;; bits=str_sep(str,'-')
;; IF n_elements(bits) EQ 2 THEN return,float(bits[0])*10.^(-float(bits[1]))

END


;---------------
PRO ch_adf04_convert, fname

IF n_params() LT 1 THEN BEGIN
  print,'Use:  ch_adf04_ups, fname'
  return
ENDIF


lvals=['S','P','D','F','G','H','I','J']


;
; Open adf04 file and extract ion information
;
openr,lin,fname,/get_lun
readf,lin,format='(5x,2i10)',iz,ion
zion2name,iz,ion,ionname

;
; Open the .elvlc file for writing.
;
elvlcname=ionname+'.elvlc_adf04'
openw,lel,elvlcname,/get_lun

energs=0.   ; this will contain array of theor. energies

tst1=0
i=0
conf=''
lev=''
en=0.
str1=''
WHILE tst1 EQ 0 DO BEGIN
  readf,lin,str1
  IF trim(str1) EQ '-1' THEN BEGIN
    tst1=1
  ENDIF ELSE BEGIN
    reads,str1,format='(i5,a17,3x,i1,1x,i1,1x,f4.0,1x,f17.0)',i,conf,spin,orb,jval,en
    conf=trim(conf)
    energs=[energs,en]
   ;
    printf,lel,format='(i7,a30,5x,i5,a5,f5.1,f15.3,f15.3)', $
         i,conf,spin,lvals[orb],jval, $
         0.,en
  ENDELSE
ENDWHILE
energs=energs[1:*]

printf,lel,' -1'
printf,lel,'%file:  '+elvlcname
printf,lel,' -1'

free_lun,lel

;
; Now read the elvlc just created. I'm doing this in order to
; access the statistical weights (elvlc.data.mult) which are used
; later. 
;
read_elvlc,elvlcname,elvlc=elvlc
wgt=elvlc.data.weight

wgfaname=ionname+'.wgfa_adf04'
upsname=ionname+'.ups_adf04'
openw,lwg,wgfaname,/get_lun
openw,lup,upsname,/get_lun

readf,lin,str1
bits=str_sep(trim(strcompress(str1)),' ')
temp=bits[2:*]
nt=n_elements(temp)
t=fltarr(nt)
FOR i=0,nt-1 DO t[i]=conv_num(temp[i])

;
; line1_form_1  -> case when limit value is defined
; line1_form_2  -> case when limit value is not defined (-1)
;
line1_form_1='(2i7,e12.3,e12.3,e12.3,i5)'
line1_form_2='(2i7,e12.3,e12.3,i12,i5)'
line2_form='('+trim(nt)+'e12.3)'   ; format for temps and upsilons


read_form='(2i4,'+trim(nt+2)+'a8)'
line_len=2*4+(nt+2)*8  ; I use this to check the data line length to see if there are problems.

tst1=0
WHILE tst1 EQ 0 DO BEGIN
  readf,lin,str1
  IF trim(str1) EQ '-1' THEN BEGIN
    tst1=1
  ENDIF ELSE BEGIN
    IF strlen(str1) NE line_len THEN print,'%CH_ADF04_CONVERT: data line length is not correct!!',strlen(str1),line_len
    l2=0
    l1=0
    bits=strarr(nt+2)
    reads,str1,format=read_form,l2,l1,bits
   ;
    aval=conv_num(bits[0])
    IF aval EQ 1e-30 THEN aval=0.   ; 1e-30 is equivalent to zero
    ups_str=bits[1:nt]
    ups=fltarr(nt)
    FOR i=0,nt-1 DO ups[i]=conv_num(ups_str[i])
    lim=conv_num(bits[1+nt])
    de_cm=abs(energs[l2-1]-energs[l1-1])
    IF de_cm EQ 0. THEN de_cm=1.   ; force energy to be non-zero
    de_ryd=de_cm/109737.32
   ;
   ; The following works out the limit point and gf value. See
   ; 'Programming Notes' in header for more details of what's
   ; going on here. 
   ;
    CASE 1 OF
      lim EQ 0.: BEGIN
        gf=0.
        lim_print=-1
        form=line1_form_2
      END
      lim LT 0.: BEGIN
        lim_print=-lim
        gf=-lim*de_ryd/4.
        form=line1_form_1
      END
      lim GT 0.: BEGIN
        IF lim EQ 1e-30 THEN BEGIN
          lim_print=0.
          gf=0.
        ENDIF ELSE BEGIN 
          lim_print=lim
          gf_lim=lim*de_ryd/4.
          gf_a=wgt[l2-1]*aval/6.67e-1/de_cm^2
          IF abs(gf_a-gf_lim)/gf_lim LE 0.02 THEN BEGIN
            gf=gf_lim
            form=line1_form_1
          ENDIF ELSE BEGIN
            gf=0.
            form=line1_form_1
          ENDELSE
        ENDELSE 
      END
    ENDCASE 
   ;
    printf,lup,format=form, $
           l1,l2,de_ryd,gf,lim_print,nt
    printf,lup,format=line2_form,t
    printf,lup,format=line2_form,ups
   ;
    wvl=1d8/de_cm
    IF aval NE 0. THEN printf,lwg,format='(2i5,f15.3,2e15.3)',l1,l2,wvl,gf,aval
  ENDELSE
ENDWHILE

printf,lup,' -1'
printf,lup,'%file:  '+upsname
printf,lup,' -1'

printf,lwg,' -1'
printf,lwg,'%file:  '+wgfaname
printf,lwg,' -1'

free_lun,lin,lwg,lup

print,'%CH_ADF04_CONVERT: created the files:'
print,'    '+elvlcname
print,'    '+upsname
print,'    '+wgfaname


END
