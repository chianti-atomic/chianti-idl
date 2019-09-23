
;+
; NAME
;
;    CHIANTI_ADF04_V2
;
; PROJECT
;
;    CHIANTI
;
; EXPLANATION
;
;    Converts an ADAS adf04 file to CHIANTI format.
;
; INPUTS
;
;    FNAME    Name of the ADF04 file.
;
; OUTPUTS
;
;    Creates CHIANTI .elvlc, .wgfa and .upsdat files from the data in
;    the adf04 file. The ion name is automatically obtained from the
;    contents of the ad04 file.
;
; HISTORY
;
;    Beta 1, 2-Feb-2009, Peter Young
;       I was having problems with earlier version of routine
;       (chianti_adf04.pro) so I've written this version.
;-


FUNCTION conv_num, str

str2=trim(str)

bits=str_sep(str,'+')
IF n_elements(bits) EQ 2 THEN return,float(bits[0])*10.^(float(bits[1]))

bits=str_sep(str,'-')
IF n_elements(bits) EQ 2 THEN return,float(bits[0])*10.^(-float(bits[1]))


END


PRO chianti_adf04_v2, fname


IF n_params() LT 1 THEN BEGIN
  print,'Use:  chianti_adf04, fname'
  return
ENDIF


openr,lin,fname,/get_lun

readf,lin,format='(5x,2i10)',iz,ion

zion2name,iz,ion,ionname

lvals=['S','P','D','F','G','H','I','J']


;
; read energies
;
openw,lel,ionname+'.elvlc',/get_lun

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
    printf,lel,format='(i3,6x,a15,2i3,a3,f4.1,i3,f15.3,f15.6,f15.3,f15.6)', $
         i,conf,spin,orb,lvals[orb],jval,fix(jval*2.+1.), $
         0.,0.,en,en/109737.32
  ENDELSE
ENDWHILE
energs=energs[1:*]

printf,lel,' -1'
printf,lel,'%file:  '+ionname+'.elvlc'
printf,lel,' -1'

free_lun,lel

openw,lwg,ionname+'.wgfa',/get_lun
openw,lup,ionname+'.upsdat',/get_lun

readf,lin,str1
bits=str_sep(trim(strcompress(str1)),' ')
temp=bits[2:*]
nt=n_elements(temp)
t=fltarr(nt)
FOR i=0,nt-1 DO t[i]=conv_num(temp[i])

printf,lup,format='(i3)',nt

form1='(2i5,f10.6,'+trim(nt)+'e10.3)'
form2='(2i5,e10.3,'+trim(nt)+'e10.3)'

tst1=0
WHILE tst1 EQ 0 DO BEGIN
  readf,lin,str1
  IF trim(str1) EQ '-1' THEN BEGIN
    tst1=1
  ENDIF ELSE BEGIN
    bits=str_sep(trim(strcompress(str1)),' ')
    l2=fix(bits[0])
    l1=fix(bits[1])
    aval=conv_num(bits[2])
    ups_str=bits[3:3+nt-1]
    ups=fltarr(nt)
    FOR i=0,nt-1 DO ups[i]=conv_num(ups_str[i])
    IF aval NE 0. THEN lim=conv_num(bits[3+nt]) ELSE lim=0.
    de=abs(energs[l2-1]-energs[l1-1])/109737.32
    gf=lim*de/4.
   ;
    printf,lup,format=form1,l1,l2,de,t
    printf,lup,format=form2,l1,l2,gf,ups
   ;
    IF aval NE 0. THEN printf,lwg,format='(2i5,f15.3,2e15.3)',l1,l2,0.,0.,aval
  ENDELSE
ENDWHILE

printf,lup,' -1'
printf,lup,'%file:  '+ionname+'.upsdat'
printf,lup,' -1'

printf,lwg,' -1'
printf,lwg,'%file:  '+ionname+'.wgfa'
printf,lwg,' -1'

free_lun,lin,lwg,lup

END
