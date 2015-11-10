;+
; Project     : SOHO - CDS     
;                   
; Name        : PLOT_CHIANTI_NE
;               
; Purpose     : Plots a density sensitive ration saved from CHIANTI_NE
;               
; Explanation : The routine CHIANTI_NE allows the calculated ratio to be
;               saved in an IDL save file.  This routine retrieves and
;               plots the ratio.
;               
; Use         : IDL> plot_chianti_ne, file, data
;    
; Inputs      : file - save file name (an extension of .CH_NE will have
;                                      been added, specifying this is optional)
;               
; Opt. Inputs : None
;               
; Outputs     : None
;               
; Opt. Outputs: data - returns the retrieved ratio structure
;               
; Keywords    : log. If set, a log-log plot is produced.
;
; Calls       : print2d_plot 
;
; Common      : None
;               
; Restrictions: None
;               
; Side effects: None
;               
; Category    : Synthetic spectra
;               
; Prev. Hist. : None
;
; Written     : C D Pike, RAL, 7-Oct-96
;               
; Modified    : v.2 Added a few extra things, including possibility to create a
;               postscript file. 
;               Giulio Del Zanna (DAMTP), 10-Oct-2000
;
; Version     : Version 2, GDZ 10-Oct-2000 
;-            


pro plot_chianti_ne, file, data=data, log=log

;
;  if no file given then allow selection
;
if n_params() eq 0 then begin
;   ffile = cdspickfile(/read,filter='*.CH_NE')
   ffile = pickfile(filter='*.CH_NE')
endif else begin
   ffile = file
endelse

if ffile eq '' then return

;
;  check file has the correct extension
;
npos = strpos(ffile,'.CH_NE')
if npos le 0 then begin
   ffile = ffile + '.CH_NE'
endif
if not file_exist(ffile) then begin
   print,'Invalid file name '+ffile
   return
endif

desc = ''
density = 0
ratio   = 0
restore,ffile


IF tag_exist(ne_ratio, 'units') THEN ytit=ne_ratio.units ELSE ytit='Ratio'
IF tag_exist(ne_ratio, 'comment') THEN print, ne_ratio.comment

tit='CHIANTI: '+ne_ratio.desc

IF tag_exist(ne_ratio, 'temperature') THEN $
  tit=tit+' calculated at constant T='+$
  string(ne_ratio.temperature,format='(e10.3)')+' (K)'


window, 0
circle_sym,/fill

x_min=min(ne_ratio.density) 
x_max=max(ne_ratio.density)
y_min=min(ne_ratio.ratio)
y_max=max(ne_ratio.ratio)


pr=''
begin_post:

!p.multi=0


IF keyword_set(log) THEN $
  plot_oo, ne_ratio.density,ne_ratio.ratio,psym=-8,$
  xrange=[x_min,x_max], yrange=[y_min,y_max],$
  xstyle=1,ystyle=1,$
  xtit='Log density (cm-3)',ytit='Log '+ytit,$
  title=tit ELSE $
  plot, ne_ratio.density,ne_ratio.ratio,psym=-8,$
  xrange=[x_min,x_max], yrange=[y_min,y_max],$
  xstyle=1,ystyle=1,$
  xtit='Density (cm-3)',ytit=ytit,$
  title=tit

print2d_plot, pr=pr , x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,$
  go_to_line=go_to_line, out_name=out_name , /ask_name

if go_to_line eq 'y' then goto,begin_post


data = ne_ratio

end


