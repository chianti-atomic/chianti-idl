;+
; Project     : SOHO - CDS     
;                   
; Name        : PLOT_CHIANTI_TE
;               
; Purpose     : Plots a temperature sensitive ratio saved from CHIANTI_TE
;               
; Explanation : The routine CHIANTI_TE allows the calculated ratio to be
;               saved in an IDL save file.  This routine retrieves and
;               plots the ratio.
;               
; Use         : IDL> plot_chianti_te, file, data
;    
; Inputs      : file - save file name (an extension of .CH_TE will have
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
; Calls       : None
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
; Modified    : v.2 Added a few extra things,  including possibility to create a
;               postscript file. 
;               Giulio Del Zanna (DAMTP), 10-Oct-2000
;
; Version     : Version 2, GDZ 10-Oct-2000 
;-            


pro plot_chianti_te, file, data=data, log=log

;
;  if no file given then allow selection
;
if n_params() eq 0 then begin
;   ffile = cdspickfile(/read,filter='*.CH_TE')
   ffile = pickfile(filter='*.CH_TE')
endif else begin
   ffile = file
endelse

if ffile eq '' then return

;
;  check file has the correct extension
;
npos = strpos(ffile,'.CH_TE')
if npos le 0 then begin
   ffile = ffile + '.CH_TE'
endif
if not file_exist(ffile) then begin
   print,'Invalid file name '+ffile
   return
endif

desc = ''
temp = 0
ratio   = 0
restore,ffile

IF tag_exist(te_ratio, 'units') THEN ytit=te_ratio.units ELSE ytit='Ratio'
IF tag_exist(te_ratio, 'comment') THEN print, te_ratio.comment

tit='CHIANTI: '+te_ratio.desc

IF tag_exist(te_ratio, 'density') THEN $
  tit=tit+' calculated at constant Ne='+$
  string(te_ratio.density,format='(e10.3)')+' (cm-3)'

window, 0
circle_sym,/fill

x_min=min(te_ratio.temp) 
x_max=max(te_ratio.temp)
y_min=min(te_ratio.ratio)
y_max=max(te_ratio.ratio)


pr=''
begin_post:

!p.multi=0

IF keyword_set(log) THEN $
  plot_oo, te_ratio.temp,te_ratio.ratio,psym=-8,$
  xrange=[x_min,x_max], yrange=[y_min,y_max],$
  xstyle=1,ystyle=1,$
  xtit='Log temperature (K)',ytit='Log '+ytit,$
  title=tit ELSE $
  plot, te_ratio.temp,te_ratio.ratio,psym=-8,$
  xrange=[x_min,x_max], yrange=[y_min,y_max],$
  xstyle=1,ystyle=1,$
  xtit='Temperature (K)',ytit=ytit,$
  title=tit

print2d_plot, pr=pr , x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,$
  go_to_line=go_to_line, out_name=out_name , /ask_name

if go_to_line eq 'y' then goto,begin_post


data = te_ratio

end


