
;+
; NAME
;
;     WGFA_PLOT_COMP
;
; EXPLANATION
;
;     Takes two .wgfa files for the same ion (and with the same level
;     indices) and makes a plot comparing the A-values.
;
;     FILE1 is treated as the reference and the A-values are plotted
;     on the X-axis. The ratios of the A-values from FILE2 to those
;     from FILE1 are plotted on the Y-axis. The Y-axis is scaled in
;     such as way that it represents 6 orders of magnitude, yet the
;     region between +/- 30 % is shown in "close-up".
;
; INPUTS
;
;     FNAME1   The name of the first file which is used as the
;              reference.
;
;     FNAME2   The name of the comparison file.
;
; OPTIONAL INPUTS
;
;     XRANGE   The X-range for the plot (default: 1e0 to 1e11).
;
;     ENCAP    The name of an encapsulated file to send the output
;              plot to.
;
;     FILE1_NAME  The name associated with FNAME1 for putting on the
;                 plot axes. Default is 'File 1'.
;
;     FILE2_NAME  The name associated with FNAME2 for putting on the
;                 plot axes. Default is 'File 2'.
;
;     XYOUT1      A label that is placed in the top-left corner of the
;                 plot. E.g., one can put the ion name.
;
; OUTPUT
;
;     Creates a plot comparing the A-values from FILE2 to those from
;     FILE1. 
;
; HISTORY
;
;     Ver.1, 5-Oct-2009, Peter Young
;-


FUNCTION ticks,axis,index,value


IF value GE 0. THEN value2=10.^((value/3.)^2*3.) ELSE $
     value2=10.^(-(value/3.)^2*3.)
nval = string(format='(e10.1)',value2)

IF nval EQ '   0.0e+00' THEN BEGIN
  nval2 = '0' 
ENDIF ELSE BEGIN
  bits = str_sep(nval,'e')
  nval2 = '10!a'+strtrim(string(fix(bits[1])),2)+'!n'
ENDELSE

IF (value2 LT 9.) AND (value2 GT 0.11) THEN nval2=string(format='(f3.1)',value2)
IF (nval2 EQ '1.1') OR (nval2 EQ '2.0') OR (nval2 EQ '5.0') THEN nval2=''
IF (nval2 EQ '0.9') OR (nval2 EQ '0.5') OR (nval2 EQ '0.2') THEN nval2=''

return,nval2

END


PRO wgfa_plot_comp, fname1, fname2, encap=encap, xrange=xrange, $
                    file1_name=file1_name,file2_name=file2_name, $
                    xyout1=xyout1

read_wgfa2,fname1,ll_1,ul_1,wvl_1,gf_1,aval_1,ref
read_wgfa2,fname2,ll_2,ul_2,wvl_2,gf_2,aval_2,ref

n=n_elements(aval_1)

IF n_elements(encap) NE 0 THEN BEGIN 
  !p.font=0
  dname=!d.name
  SET_PLOT,'ps'
  file=encap
  DEVICE,/encap,file=file,xsiz=5,ysiz=4,/inches
ENDIF

IF n_elements(file1_name) EQ 0 THEN file1_name='File 1'
IF n_elements(file2_name) EQ 0 THEN file2_name='File 2'

xtitle=file1_name+' A-values / s!u-1!n'
ytitle='A-value ratio ('+file2_name+'/'+file1_name+')'

IF n_elements(xrange) EQ 0 THEN xrange=[1d0,1d11]

plot,/nodata,/xlog,xrange,[-3,3], $
     ytickformat='ticks', $
     xtitle=xtitle,ytitle=ytitle,xth=3,yth=3, $
     ytickv=[-3,-2.449,-1.732,-1.448, -0.95031,-0.68169,-0.3705, $
             0.,0.35239,0.58466,0.95031, 1.448, 1.732,2.449,3], $
     yticks=16,/xsty

oplot,[1d-8,1d12],[0,0],th=3
oplot,[1d-8,1d20],[1,1]*0.58466,line=2,th=3
oplot,[1d-8,1d20],[1,1]*(-0.68169),line=2,th=3

count=0

;
; Go through transitions in FILE1
;
FOR i=0,n-1 DO BEGIN
  ind=where(ll_1[i] EQ ll_2 AND ul_1[i] EQ ul_2)
  IF ind[0] NE -1 THEN BEGIN
    ratio=alog10(aval_2[ind[0]]/aval_1[i])
    IF ratio GT 3 OR ratio LT -3 THEN BEGIN
      count=count+1
    ENDIF ELSE BEGIN 
      jj=where(ratio GE 0)
      IF jj[0] NE -1 THEN BEGIN
        plot_rat=(ratio[jj]/3.)^0.5 *3.
        plots,aval_1[i],plot_rat,psym=2
      ENDIF
      jj=where(ratio LT 0)
      IF jj[0] NE -1 THEN BEGIN
        plot_rat=-(-ratio[ind]/3.)^0.5 *3.
        plots,aval_1[i],plot_rat,psym=2
      ENDIF
    ENDELSE
  ENDIF

ENDFOR

IF count GT 0 THEN BEGIN
  print,'*** '+trim(count)+' transitions lie outside plot range ***'
ENDIF 

IF n_elements(xyout1) NE 0 THEN BEGIN
  yr=!y.crange
  xr=!x.crange
  xyouts,/normal,0.22,0.85,xyout1
ENDIF 

IF keyword_set(encap) THEN BEGIN
  device,/close
  set_plot,dname
  !p.font=-1
ENDIF


END
