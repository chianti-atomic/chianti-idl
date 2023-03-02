

function ch_plot_iso_pops, ionname, level, ldens=ldens, neutrals=neutrals, quiet=quiet, $
                           outstr=outstr, config_match=config_match, _extra=extra, $
                           lookup=lookup


;+
; NAME:
;	CH_PLOT_ISO_POPS
;
; PURPOSE:
;	This routine plots the population of a specific atomic level
;	as a function of the member of the isoelectronic sequence.
;
; CATEGORY:
;	CHIANTI; level populations.
;
; CALLING SEQUENCE:
;       CH_PLOT_ISO_POPS, Ionname, Level
;
; INPUTS:
;	Ionname: The name of an ion in CHIANTI format. For example,
;	        'fe_13' for Fe XIII.
;       Level:  The CHIANTI level index for an atomic level belonging
;               to the ion identified by Ionname. Can be an array of
;               indices (in which case all levels will be shown on the
;               output plot).
;
; OPTIONAL INPUTS:
;	Ldens:	The logarithm of the electron number density (units:
;	        cm^-3). If not specified then a value of 10 is assumed.
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
;	A plot object is created showing how the populations vary along the
;	isoelectronic sequence. The X-axis shows the atomic
;	number. 
;
; RESTRICTIONS:
;	The routine assumes that the routine Ch_find_iso_level was
;	able to correctly identify the atomic level for each ion on
;	the sequence. 
;
; EXAMPLE:
;       IDL> p=ch_plot_iso_pops('o_4', 4)
;       IDL> p=ch_plot_iso_pops('o_4', [1,2,3,4,5])
;       IDL> p=ch_plot_iso_pops('ne_6', 12, ldens=9.0)
;
; MODIFICATION HISTORY:
;       Ver.1, 11-Feb-2014, Peter Young
;       Ver.2, 21-May-2014, Peter Young
;          Added /config_match, /quiet and outstr=.
;       Ver.3, 7-Aug-2017, Peter Young
;          Now creates plot object; added titles; converted from
;          procedure to function (to return object).
;       Ver.4, 29-Jan-2022, Peter Young
;          LEVEL can now be an array.
;       Ver.5, 02-Mar-2023, Peter Young
;          The plot was not showing the input ion (IONNAME), so this has been
;          fixed.
;-


IF n_params() LT 2 THEN BEGIN
  print,'Use:  IDL> p=ch_plot_iso_pops( ion_name, level [, ldens=, /neutrals, /quiet, '
  print,'                               outstr=, /config_match ])'
  return,-1
ENDIF 

IF n_elements(ldens) EQ 0 THEN ldens=10.0
IF keyword_set(neutrals) THEN i_neut=0 ELSE i_neut=1

convertname,ionname,iz,ion
diff=iz-ion

nlev=n_elements(level)

;
; Read the search ion data into the structure 'str'
;
zion2filename,iz,ion,fname
elvlcname=fname+'.elvlc'
;
read_elvlc,elvlcname,elvlcstr=str
latex_str=strarr(nlev)
print,'Selected levels are: '
FOR i=0,nlev-1 DO BEGIN
  k=where(str.data.index EQ level[i],nk)
  print,format='(3x,i5,a20,"  -- plot_index:",i3)',level[i],str.data[k[0]].full_level,i+1
;  print,'    ',strpad(i+1str.data[k[0]].full_level
  latex_str[i]=str.data[k[0]].full_level_latex
ENDFOR 



mlistname=concat_dir(!xuvtop,'masterlist')
mlistname=concat_dir(mlistname,'masterlist.ions')
read_masterlist,mlistname,mlist


;
; Define output structure.
;
str={ion: '', spect: 0, lev: lonarr(nlev), pop: dblarr(nlev), latex: latex_str, element: 0}
outstr=0

FOR i=diff+1+i_neut,30 DO BEGIN
  zion2name,i,i-diff,iname
  k=where(trim(iname) EQ mlist,nk)
  IF nk NE 0 THEN BEGIN

    str.ion=iname
    str.element=i
    str.spect=i-diff

    l=lonarr(nlev)
    FOR j=0,nlev-1 DO BEGIN 
      l[j]=ch_find_iso_level(iname,ionname,level[j],/quiet,outlev=outlev,config_match=config_match)
    ENDFOR
    str.lev=l
    
    tmax=ch_tmax(iname)
    IF keyword_set(lookup) THEN BEGIN
      pop=ch_lookup_pops(iname,ldens=ldens,temp=tmax)
    ENDIF ELSE BEGIN
      pop=ch_pops(iname,ldens=ldens,temp=tmax,/quiet)
    ENDELSE

    ind=where(l NE -1,nind)
    IF nind NE 0 THEN str.pop[ind]=pop.level[l-1].pop
    
;    IF NOT keyword_set(quiet) THEN print,format='(5x,a5,2x,a15,e10.2)',iname,outlev,str.pop
     ;
    IF n_tags(outstr) EQ 0 THEN outstr=str ELSE outstr=[outstr,str]
  ENDIF
ENDFOR 



IF NOT keyword_set(quiet) THEN BEGIN
  w=window(dim=[700,500])
  ss=2.0
  fs=14
  xrange=[min(outstr.element)-1,max(outstr.element)+1]
  yrange=[min(outstr.pop),max(outstr.pop)]

  p=plot(/nodata, xrange, $
         yrange,  $
         ytitle='Level population', $
         xtitle='Atomic number',/current,xticklen=0.020,yticklen=0.015,/ylog, $
         font_size=14,_extra=extra,/xsty)

  FOR i=0,nlev-1 DO BEGIN
    pl=plot(outstr.element,outstr.pop[i],th=th,/overplot)
    n_ions=n_elements(outstr)
    FOR j=0,n_ions-1 DO  pt=text(/data,align=0.5,vertical_align=0.5, $
                             outstr[j].element,outstr[j].pop[i], $
                             trim(i+1),font_size=fs,target=p)
  ENDFOR

  p.xrange=xrange
  
  ;; p=plot(outstr[k].element,outstr[k].pop,symbol='+',sym_size=ss, $
  ;;        xrange=[min(outstr.element)-1,max(outstr.element)+1], $
  ;;        ytitle='Level population', $
  ;;        xtitle='Atomic number',/current,xticklen=0.015,yticklen=0.015,/ylog, $
  ;;        title='Level '+trim(level)+': '+latex_str,linestyle='none', $
  ;;        font_size=14,_extra=extra)
  ;; k=where(outstr.element eq iz) 
  ;; q=plot(/overplot,[1,1]*outstr[k].element,[1,1]*outstr[k].pop,symbol='triangle', $
  ;;        sym_size=ss,_extra=extra)
  return,p
ENDIF ELSE BEGIN
  return,-1
ENDELSE 



END
