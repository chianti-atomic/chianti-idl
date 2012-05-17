
PRO plot_config_energies, ionname, elvlcfile=elvlcfile, ryd=ryd, ev=ev, yrange=yrange

;+
; NAME
;
;    PLOT_CONFIG_ENERGIES
;
; PROJECT
;
;    CHIANTI
;
; EXPLANATION
;
;    This makes a plot showing the range of energies of configurations in the
;    CHIANTI models.
;
; INPUTS
;
;    IONNAME   E.g., 'fe_13' for Fe XIII.
;
; OPTIONAL INPUTS
;
;    ELVLCFILE  Directly specify a .elvlc file, rather than use the file
;               from the CHIANTI database.
;
; KEYWORDS
;
;    RYD     Use Rydbergs for energy scale.
;
;    EV      Use eV for energy scale.
;
; OUTPUTS
;
;    Creates a plot window showing the configuration
;    energies. Configurations are arranged along the X-axis according
;    to their CHIANTI configuration indices. For each a box is drawn
;    indicating the minimum and maximum values of the level energies
;    belonging to that configuration. The Y-axis shows the energy.
; 
; CALLS
;
;    ION2FILENAME, READ_ELVLC, TRIM
;
; HISTORY
;
;    Ver.1, 13-Feb-2009, Peter Young
;-


IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> plot_config_energies, ionname [,elvlcname=, /eV, /Ryd]'
  return
ENDIF 

IF n_elements(elvlcfile) EQ 0 THEN BEGIN
  ion2filename,ionname,rootname
  elvlcfile=rootname+'.elvlc'
ENDIF

read_elvlc,elvlcfile,l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,ref

str={n: 0, min: 0., max: 0., name: ''}

FOR i=0,100 DO BEGIN
  j=where(conf EQ i)
  IF j[0] NE -1 THEN BEGIN
    str.n=i
    str.max=max(ecm[j])
    str.min=min(ecm[j])
    IF str.min EQ 0. THEN str.min=min(ecmth[j])
    IF str.max EQ 0. THEN str.max=max(ecmth[j])
    name=term[j[0]]
    bits=str_sep(name,' ')
    str.name=''
    FOR k=0,n_elements(bits)-2 DO BEGIN
      str.name=str.name+bits[k]+' '
    ENDFOR

    IF n_tags(struc) EQ 0 THEN struc=str ELSE struc=[struc,str]
  ENDIF
ENDFOR

CASE 1 OF 
  keyword_set(ryd): scale=109737.32
  keyword_set(ev): scale=8065.54
  ELSE: scale=1.0
ENDCASE
;
struc.min=struc.min/scale
struc.max=struc.max/scale

CASE 1 OF
  keyword_set(ev): en_units='eV'
  keyword_set(ryd): en_units='Ryd'
  ELSE: en_units='cm!u-1!n'
ENDCASE 

IF n_elements(yrange) EQ 0 THEN yrange=[0,max(struc.max)*1.1]

n=n_elements(struc)
plot,/nodata,[0,max(struc.n)+1.0],yrange,/xsty,/ysty, $
     ytit='Energy / '+en_units, $
     xtit='Configuration index', $
     xticklen=1e-6
FOR i=0,n-1 DO BEGIN
  xbox=struc[i].n+[-.5,.5,.5,-.5,-.5]
  y0=struc[i].min
  y1=struc[i].max
  ybox=[y0,y0,y1,y1,y0]
  oplot,xbox,ybox
  xyouts,align=0.5,struc[i].n,0.5*(y0+y1),trim(struc[i].name)
ENDFOR


END
