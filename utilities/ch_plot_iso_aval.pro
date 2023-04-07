

function ch_plot_iso_aval, ionname, upper_level, neutrals=neutrals, quiet=quiet, $
                           outstr=outstr, config_match=config_match, _extra=extra


;+
; NAME:
;	CH_PLOT_ISO_AVAL
;
; PURPOSE:
;	This routine plots the A-values for transitions from the
;       specified upper atomic level along the isoelectronic sequence.
;
; CATEGORY:
;	CHIANTI; A-values.
;
; CALLING SEQUENCE:
;       Result = CH_PLOT_ISO_AVAL( Ionname, Upper_Level )
;
; INPUTS:
;	Ionname: The name of an ion in CHIANTI format. For example,
;	         'fe_13' for Fe XIII.
;       Upper_Level:  The CHIANTI level index for an atomic level belonging
;               to the ion identified by Ionname. Must be a scalar.
;	
; KEYWORD PARAMETERS:
;	NEUTRALS: By default the neutral element on the sequence is
;    	         omitted, so setting this keyword enables it to be shown.
;       CONFIG_MATCH: If set, then a configuration match is required
;               otherwise no match will be returned.
;       QUIET:   If set, then a plot is not produced.
;
; OPTIONAL OUTPUTS:
;       Outstr:  A structure containing the data that has been
;                plotted. 
;
; OUTPUTS:
;	A plot object is created showing how the A-values vary along the
;	isoelectronic sequence. The X-axis shows the atomic
;	number. For each transition, an integer index is assigned and
;       the key to the right of the plot identifies the lower level(s)
;       for the transition(s).
;
; RESTRICTIONS:
;	The routine assumes that the routine Ch_find_iso_level was
;	able to correctly identify the atomic level for each ion on
;	the sequence. 
;
; EXAMPLE:
;       IDL> p=ch_plot_iso_aval('o_4', 8)
;
; MODIFICATION HISTORY:
;       Ver.1, 07-Apr-2023, Peter Young
;-


IF n_params() LT 2 THEN BEGIN
  print,'Use:  IDL> p=ch_plot_iso_aval( ion_name, upper_level [, /neutrals, /quiet, '
  print,'                               outstr=, /config_match ])'
  return,-1
ENDIF 

IF keyword_set(neutrals) THEN i_neut=0 ELSE i_neut=1

convertname,ionname,iz,ion
diff=iz-ion


nlev=n_elements(upper_level)

;
; Read the search ion data into the structure 'str'
;
zion2filename,iz,ion,fname
z2element,iz-ion+1,iso_elt
elvlcname=fname+'.elvlc'
wgfaname=fname+'.wgfa'
;
read_elvlc,elvlcname,elvlcstr=str
latex_str=strarr(nlev)
upper_latex=str.data[upper_level-1].full_level_latex

read_wgfa_str,wgfaname,wgfa
k=where(wgfa.lvl2 EQ upper_level,nk)
avals=wgfa[k].aval
i=where(avals GE max(avals)/1e5)
lower_levels=wgfa[k[i]].lvl1

nlev=n_elements(lower_levels)

mlistname=concat_dir(!xuvtop,'masterlist')
mlistname=concat_dir(mlistname,'masterlist.ions')
read_masterlist,mlistname,mlist


;
; Define output structure.
;
str={ion: '', spect: 0, lvl1: lonarr(nlev), lvl2: lonarr(nlev), $
     aval: fltarr(nlev), $
     lower_latex: strarr(nlev), upper_latex: upper_latex, $
     element: 0}
outstr=0

FOR i=diff+1+i_neut,30 DO BEGIN
  zion2name,i,i-diff,iname
  zion2filename,i,i-diff,fname
 ;
  k=where(trim(iname) EQ mlist,nk)
  IF nk NE 0 THEN BEGIN
    read_wgfa_str,fname+'.wgfa',wgfa
    read_elvlc,fname+'.elvlc',elvlc=elvlc

    str.ion=iname
    str.element=i
    str.spect=i-diff

    
    FOR j=0,nlev-1 DO BEGIN 
      str.lvl2[j]=ch_find_iso_level(iname,ionname,upper_level,/quiet,outlev=outlev,config_match=config_match)
      str.lvl1[j]=ch_find_iso_level(iname,ionname,lower_levels[j],/quiet,outlev=outlev,config_match=config_match)
     ;
      ii=where(wgfa.lvl1 EQ str.lvl1[j] AND wgfa.lvl2 EQ str.lvl2[j],nii)
      IF nii NE 0 THEN BEGIN
        str.aval[j]=wgfa[ii[0]].aval
        str.lower_latex[j]=elvlc.data[str.lvl1[j]-1].full_level_latex
      ENDIF 
    ENDFOR
   ;
    IF str.lvl2[0] NE -1 THEN BEGIN 
      IF n_tags(outstr) EQ 0 THEN outstr=str ELSE outstr=[outstr,str]
    ENDIF 
  ENDIF
ENDFOR 


x0=0.10
x1=0.75
y0=0.10
y1=0.90
pos=[x0,y0,x1,y1]

IF NOT keyword_set(quiet) THEN BEGIN
  w=window(dim=[850,500])
  ss=2.0
  fs=14
  xrange=[min(outstr.element)-1,max(outstr.element)+1]
  yrange=[min(outstr.aval),max(outstr.aval)]

  title=iso_elt+' sequence'
  p=plot(/nodata, xrange, $
         yrange,  $
         ytitle='A-value / s!u-1!n', $
         xtitle='Atomic number',/current,xticklen=0.020,yticklen=0.015,/ylog, $
         font_size=fs,_extra=extra,/xsty, $
         title=title,pos=pos)

  ltxt=text(x1+0.02,0.95,'Decays from:', $
            font_size=fs-2,target=p,vertical_align=1.0)
  ltxt2=text(x1+0.02,0.95-0.06,outstr[0].upper_latex, $
            font_size=fs-2,target=p,vertical_align=1.0)
  ltxt=text(x1+0.02,0.95-0.12,'to the levels:', $
            font_size=fs-2,target=p,vertical_align=1.0)
  
  FOR i=0,nlev-1 DO BEGIN
    pl=plot(outstr.element,outstr.aval[i],th=th,/overplot)
    n_ions=n_elements(outstr)
    FOR j=0,n_ions-1 DO  pt=text(/data,align=0.5,vertical_align=0.5, $
                             outstr[j].element,outstr[j].aval[i], $
                                 trim(i+1),font_size=fs,target=p)
    ltxt=text(x1+0.02,0.95-(i+3)*0.06,trim(i+1)+' - '+outstr[0].lower_latex[i], $
              font_size=fs-2,target=p,vertical_align=1.0)
  ENDFOR

  p.xrange=xrange
  
  return,p
ENDIF ELSE BEGIN
  return,-1
ENDELSE 



END
