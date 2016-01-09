

PRO ch_plot_iso_pops, ionname, level, ldens=ldens, neutrals=neutrals, quiet=quiet, $
                      outstr=outstr, config_match=config_match


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
;	A plot is created showing how the populations vary along the
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
; 	Written by:	Peter Young, 11-Feb-2014.
;
;       Ver.2, 21-May-2014, Peter Young
;          Added /config_match, /quiet and outstr=.
;-


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

str={ion: '', spect: 0, lev: 0, pop: 0., latex: ''}
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
      str.spect=i-diff
      str.pop=pstr.level[l-1].pop
      str.latex=latex_str
      print,format='(5x,a5,2x,a15,e10.2)',iname,outlev,str.pop
     ;
      IF n_tags(outstr) EQ 0 THEN outstr=str ELSE outstr=[outstr,str]
    ENDIF
  ENDIF  
ENDFOR 

IF NOT keyword_set(quiet) THEN BEGIN 
  plot,outstr.spect,outstr.pop,psym=2,/ylog, $
       xra=[min(outstr.spect)-1,max(outstr.spect)+1],/xsty
ENDIF 

END
