
FUNCTION ch_table_summary, reverse=reverse, states=states, reduced=reduced, $
                           all=all, simple=simple



;+
; NAME:
;     CH_TABLE_SUMMARY
;
; PURPOSE:
;     Creates a table showing all of the ions that are in CHIANTI,
;     with a number giving the number of levels for the ion. A color
;     coding is also applied to each box corresponding to the numbers
;     of levels.
;
; CATEGORY:
;     CHIANTI; display.
;
; CALLING SEQUENCE:
;     Result = Ch_Table_Summary()
;
; INPUTS:
;     None.
;
; KEYWORD PARAMETERS:
;     REVERSE:  Uses a reversed color table (i.e., white on
;               black). This also changes the colors of the boxes.
;     STATES:   The word "levels" is replaced with "states".
;     REDUCED:  If set, then only the more abundant elements are
;               shown.
;     ALL:      If set, then all elements are shown, including those
;               for which we don't have any data.
;     SIMPLE:   If set, then the table entries are filled circles
;               rather than the number of levels for the ion, and
;               there are no colors.
;
; OUTPUTS:
;     Returns an IDL plot object containing the table. The table can
;     be saved in a variety of formats. For example:
;      IDL> w=ch_table_summary()
;      IDL> w.save,'ch_table_summary.eps'
;      IDL> w.save,'ch_table_summary.pdf'
;      IDL> w.save,'ch_table_summary.png',resolution=192
;
; CALLS:
;     CH_GET_VERSION, CONVERTNAME, CH_LOOKUP_TABLE_INTERP
;
; MODIFICATION HISTORY:
;     Ver.1, 21-Mar-2019, Peter Young
;     Ver.2, 29-Jun-2023, Peter Young
;       Now uses the lookup tables to speed up the routine; added
;       /simple keyword to give a more simple table; no longer
;       creates a save file containing the ion level numbers.
;-



el_all=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si',$
        'P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co', $
        'Ni','Cu','Zn']

ionstage=['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII', $
          'XIII','XIV','XV','XVI','XVII','XVIII','XIX','XX','XXI','XXII', $
          'XXIII','XXIV','XXV','XXVI','XXVII','XXVIII','XXIX','XXX']

IF keyword_set(states) THEN levstring='states' ELSE levstring='levels'



IF keyword_set(reduced) THEN BEGIN 
  element=['H','He','C','N','O','Ne','Mg', $
           'Al','Si','S','Ar','Ca', $
           'Fe','Ni']
  ionstage=ionstage[0:27]
ENDIF ELSE BEGIN 
  element=['H','He','C','N','O','Ne','Na','Mg', $
           'Al','Si','P','S','Cl','Ar','K','Ca','Ti','Cr', $
           'Mn','Fe','Ni','Zn']
ENDELSE

IF keyword_set(all) THEN element=el_all


IF keyword_set(reverse) THEN BEGIN
  bgcolor='black'
  color='white'
  bcol1='medium sea green'
  bcol2='blue'
  bcol3='red'
ENDIF ELSE BEGIN
  bgcolor='white'
  color='black'
  bcol1='yellow'
  bcol2='aqua'
  bcol3='light salmon'
ENDELSE 

;
; Getting the numbers of levels for each ion takes a while, so I
; output them to a save file that can be restored for the next run. 
;
ver=ch_get_version()
savefile='levels_v'+trim(ver)+'.save'
chck=file_info(savefile)
IF NOT chck.exists THEN BEGIN
  lev_array=intarr(30,30)
ENDIF ELSE BEGIN
  restore,savefile
ENDELSE 



dimen=[800,600]
angle=43.5
IF keyword_set(reduced) THEN BEGIN
  dimen=[800,450]
ENDIF
IF keyword_set(all) THEN BEGIN
  dimen=[800,700]
ENDIF



w=window(dim=dimen,background_color=bgcolor)

th=2
fs=11
fs2=8


nel=n_elements(element)
ni=n_elements(ionstage)

xra=[-1.5,nel+0.5]

p=plot(/nodata,[-1,float(ni)+0.5],xra,xshowtext=0,yshowtext=0,xth=th,yth=th,xsty=1,ysty=1,/current, $
       axis_style=0, xticklayout=1,yticklayout=1,xminor=1,yminor=0,xticklen=0,yticklen=0, $
      pos=[0.05,0.05,0.95,0.95],color=color)


FOR i=0,nel-1 DO BEGIN
  k=where(trim(element[i]) EQ trim(el_all))
  zz=k[0]
  FOR j=0,ni-1 DO BEGIN
    ionname=strlowcase(element[i])+'_'+trim(j+1)
    convertname,ionname,iz,ion
    
    IF ion GT iz THEN continue
    
    p=ch_lookup_table_interp(ionname,1e10,ch_tmax(ionname),/quiet)
    IF n_tags(p) NE 0 THEN BEGIN
      s=size(p.pop,/dim)
      mlev=s[2]
    ENDIF ELSE BEGIN
      mlev=0
    ENDELSE 
    
    IF keyword_set(simple) THEN BEGIN
      q=plot(j+[0.5,1.5,1.5,0.5,0.5],xra[1]-i-0.5+[-0.5,-0.5,0.5,0.5,-0.5], $
             /overplot,thick=th)
      IF mlev GT 0 THEN BEGIN
        xpos=j+1
        ypos=nel-i
        s=plot(/overplot,xpos*[1,1],ypos*[1,1],symbol='o',/sym_filled,sym_size=1)
      ENDIF 
    ENDIF ELSE BEGIN 
    
      IF mlev EQ 0 THEN label='-' ELSE label=trim(mlev)
      IF mlev GT 0 THEN BEGIN
        CASE 1 OF
          mlev GT 100: col=bcol3
          mlev LE 100 AND mlev GT 10: col=bcol2
          mlev LE 10: col=bcol1
        ENDCASE
        q=plot(j+[0.5,1.5,1.5,0.5,0.5],xra[1]-i-0.5+[-0.5,-0.5,0.5,0.5,-0.5], $
               fill_color=col,/overplot,fill_background=1,color=color,thick=th)
      ENDIF ELSE BEGIN
        q=plot(j+[0.5,1.5,1.5,0.5,0.5],xra[1]-i-0.5+[-0.5,-0.5,0.5,0.5,-0.5], $
               /overplot,color=color,thick=th)
      ENDELSE 
      t1=text(j+1,nel-i-0.25,label,align=0.5,font_size=fs2,color=color,/data)
    ENDELSE 
  ENDFOR
ENDFOR


;
; Plots the element symbols on the Y-axis.
;
FOR i=0,nel-1 DO eltxt=text(-0.9,xra[1]-0.5-i-0.2,'!3'+element[i],font_size=fs,/data,color=color)

;
; Plots the roman numerals on the X-axis
;
ni=n_elements(ionstage)
FOR i=0,ni-1 DO BEGIN
  etxt=text(/data,i-0.4,-1.2,ionstage[i],align=0.,orient=angle,font_size=fs,color=color)
  ep=plot([0.+i-1.2,i+.5],[-1.2,0.5],th=th,/overplot,color=color)
ENDFOR



x0=20.5
y0=18
IF keyword_set(reduced) THEN BEGIN 
  x0=18.5
  y0=10
ENDIF
IF keyword_set(all) THEN BEGIN
  y0=20
ENDIF 


t2=text(x0-0.5,y0+1,'Ions in CHIANTI '+ver,font_size=fs+2,/data,color=color)

IF ~ keyword_set(simple) THEN BEGIN 
  b1=plot(x0+[0,1,1,0,0],y0-1+[0,0,1,1,0],fill_background=1,fill_color=bcol1, $
          thick=th,color=color,/overplot)
  b2=plot(x0+[0,1,1,0,0],y0-1-1.2+[0,0,1,1,0],fill_background=1,fill_color=bcol2, $
          thick=th,color=color,/overplot)
  b3=plot(x0+[0,1,1,0,0],y0-1-2.4+[0,0,1,1,0],fill_background=1,fill_color=bcol3, $
          thick=th,color=color,/overplot)
  bt1=text(x0+1.5,y0-0.7,'$\le$ 10 '+levstring,font_size=fs,/data,color=color)
  bt2=text(x0+1.5,y0-1.2-0.7,'< 100 '+levstring,font_size=fs,/data,color=color)
  bt3=text(x0+1.5,y0-2.4-0.7,'$\ge$ 100 '+levstring,font_size=fs,/data,color=color)
ENDIF 


return,w

END


