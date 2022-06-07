
PRO ch_config_lines, ionname, conf1, conf2, output=output


;+
; NAME:
;     CH_CONFIG_LINES
;
; PURPOSE:
;     Prints the list of transitions between the two specified
;     configurations. 
;
; CATEGORY:
;     CHIANTI; transitions.
;
; CALLING SEQUENCE:
;     CH_CONFIG_LINES, IonName, Conf1, Conf2 
;
; INPUTS:
;     IonName:  The name of an ion in CHIANTI format. For example,
;               'o_6' for O VI.
;     Conf1:    An integer specifying the lower configuration. The
;               ground configuration has index 1.
;     Conf2:    An integer specifying the upper configuration.
;
; OUTPUTS:
;     Prints a list of the transitions between the two
;     configurations. The columns are: wavelength, lower level index,
;     upper level index, A-value and transition information.
;
; OPTIONAL OUTPUTS:
;     Output:  An IDL structure that contains the information printed
;              to the screen.
;
; EXAMPLE:
;     IDL> ch_config_lines, 'fe_13', 1, 3
;
; MODIFICATION HISTORY:
;     Ver.1, 02-Feb-2022, Peter Young
;-


convertname,ionname,iz,ion
zion2filename,iz,ion,fname

read_elvlc,fname+'.elvlc',elvlc=elvlc
read_wgfa_str,fname+'.wgfa',wgfa

ii=where(elvlc.data.conf_index EQ conf1,ni)
jj=where(elvlc.data.conf_index EQ conf2,nj)

l1=elvlc.data[ii].index
l2=elvlc.data[jj].index

str={ lvl1: 0, lvl2: 0, wvl: 0., aval: 0., lvl1_str: '', lvl2_str: ''}
output=0

FOR i=0,ni-1 DO BEGIN
  FOR j=0,nj-1 DO BEGIN
    k=where(wgfa.lvl1 EQ l1[i] AND wgfa.lvl2 EQ l2[j],nk)
    IF nk NE 0 THEN BEGIN
      str.lvl1=l1[i]
      str.lvl2=l2[j]
      str.wvl=wgfa[k[0]].wvl
      str.aval=wgfa[k[0]].aval
      str.lvl1_str=elvlc.data[ii[i]].full_level
      str.lvl2_str=elvlc.data[jj[j]].full_level
     ;
      IF n_tags(output) EQ 0 THEN output=str ELSE output=[output,str]
    ENDIF
  ENDFOR 
ENDFOR 

IF n_tags(output) NE 0 THEN BEGIN 
  k=sort(abs(output.wvl))
  output=output[k]

  n=n_elements(output)
  FOR i=0,n-1 DO BEGIN
    trans=trim(output[i].lvl1_str)+' - '+trim(output[i].lvl2_str)
    trans=strpad(trans,60,/after)
    print,format='(f12.3,2i5,e12.2,3x,a60)',output[i].wvl,output[i].lvl1,output[i].lvl2, $
          output[i].aval,trans
  ENDFOR
ENDIF ELSE BEGIN
  print,'% CH_CONFIG_LINES: no transitions between these configurations.'
ENDELSE 


END
