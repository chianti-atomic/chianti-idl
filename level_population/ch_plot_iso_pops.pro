

function ch_plot_iso_pops, ionname, level, ldens=ldens, neutrals=neutrals, quiet=quiet, $
                      outstr=outstr, config_match=config_match, _extra=extra


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
;               to the ion identified by Ionname.
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
;	isoelectronic sequence. The X-axis shows the spectroscopic
;	number. For example, 13 would correspond to XIII. In addition,
;	the routine prints information to the IDL command window
;	showing which level has been found for each ion. The user
;	should check to make sure the correct level has been selected.
;
; RESTRICTIONS:
;	The routine assumes that the routine Ch_find_iso_level was
;	able to correctly identify the atomic level for each ion on
;	the sequence. 
;
; EXAMPLE:
;       ch_plot_iso_pops, 'o_4', 4
;       ch_plot_iso_pops, 'ne_6', 12, ldens=9.0
;
; MODIFICATION HISTORY:
;       Ver.1, 11-Feb-2014, Peter Young
;       Ver.2, 21-May-2014, Peter Young
;          Added /config_match, /quiet and outstr=.
;       Ver.3, 7-Aug-2017, Peter Young
;          Now creates plot object; added titles; converted from
;          procedure to function (to return object).
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

;
; Read the search ion data into the structure 'str'
;
zion2filename,iz,ion,fname
elvlcname=fname+'.elvlc'
;
read_elvlc,elvlcname,elvlcstr=str
k=where(str.data.index EQ level,nk)
print,'Selected level is: ',str.data[k].full_level

latex_str=str.data[k].full_level_latex


mlistname=concat_dir(!xuvtop,'masterlist')
mlistname=concat_dir(mlistname,'masterlist.ions')
read_masterlist,mlistname,mlist

str={ion: '', spect: 0, lev: 0, pop: 0., latex: '', element: 0}
outstr=0

FOR i=diff+1+i_neut,30 DO BEGIN
  zion2name,i,i-diff,iname
  k=where(trim(iname) EQ mlist,nk)
  IF nk NE 0 THEN BEGIN
    l=ch_find_iso_level(iname,ionname,level,/quiet,outlev=outlev,config_match=config_match)
    IF l NE -1 THEN BEGIN 
      tmax=get_tmax(iname)
      show_pops,i,i-diff,pstr,dens=ldens,temp=alog10(tmax),/quiet
     ;
      str.ion=iname
      str.lev=l
      str.element=i
      str.spect=i-diff
      str.pop=pstr.level[l-1].pop
      str.latex=latex_str
      print,format='(5x,a5,2x,a15,e10.2)',iname,outlev,str.pop
     ;
      IF n_tags(outstr) EQ 0 THEN outstr=str ELSE outstr=[outstr,str]
    ENDIF
  ENDIF  
ENDFOR 

w=window(dim=[700,500])

ss=2.0
fs=14

IF NOT keyword_set(quiet) THEN BEGIN
  k=where(outstr.element NE iz) 
  p=plot(outstr[k].element,outstr[k].pop,symbol='+',sym_size=ss, $
         xrange=[min(outstr.element)-1,max(outstr.element)+1], $
         ytitle='Level population', $
         xtitle='Atomic number',/current,xticklen=0.015,yticklen=0.015,/ylog, $
         title='Level '+trim(level)+': '+latex_str,linestyle='none', $
         font_size=14,_extra=extra)
  k=where(outstr.element eq iz) 
  q=plot(/overplot,[1,1]*outstr[k].element,[1,1]*outstr[k].pop,symbol='triangle', $
         sym_size=ss,_extra=extra)
  return,p
ENDIF ELSE BEGIN
  return,-1
ENDELSE 



END
