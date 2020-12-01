;+
; Project     : SOHO - CDS     
;                   
; Name        : PLOT_DMM_DR_FIG
;               
; Purpose     : Plots line intensities in DMM_DR
;               
; Explanation : Internal routine used by CHIANTI_NE for plotting individual
;               line intensities.
;               
; Use         : Called by CHIANTI_NE
;    
; Inputs      : Uses common
;               
; Opt. Inputs : None
;               
; Outputs     : None
;               
; Opt. Outputs: None
;               
; Keywords    : None
;
; Calls       : None
;
; Common      : dmm_lines
;               
; Restrictions: None
;               
; Side effects: None
;               
; Category    : Spectral, plotting
;               
; Prev. Hist. : Originally included in density_ratios by K Dere
;
; Written     : C D Pike, RAL, 15-Jan-96
;               
; Modified    : V.3 Fix bug in searching for peak line (IMAX).  CDP, 27-Jan-96
;               Changed charsize in plots.  CDP, 17-Jul-97
;               V.4 changed yranges in the plots. Removed common calls not
;               needed. Giulio Del Zanna 6-Oct-2000 
;
;
; Version     : Version 4, 6-Oct-2000
;-            

pro plot_dmm_dr_fig

common dmm_lines, list_wvl, list_int, list_descr1, list_descr2, species,$
  density, nden, ratio, description, nlist , savetext,$
  temper, intens

;
intmax=max(list_int)
intmin= min(list_int) > intmax*intens 
;

;
plot_oo,fltarr(2),xr=[density(0),density(nden-1)],yr=[intmin,intmax],$
  tit='Rel. Int. (Pop. of upper level)*A/Density',$
  chars=0.8, xtitle='Density (cm-3)'

;
for ilist=0,nlist-1 do begin
   oplot,density,list_int(ilist,*)
;   oplot,density,list_int(ilist,*), psym=1
endfor
;
rand=randomu(seed,nden)
isrt=sort(rand)                 ; to try to fix problems with labeling curves
;
;
for ilist=0,nlist-1 do begin
   iden=isrt(ilist<(nden-1))
   xyouts,density(iden),list_int(ilist,iden),string(ilist,'(i3)'),chars=1.5,$
     align=0.5,noclip=0
   xyouts,density(1),list_int(ilist,1),string(ilist,'(i3)'),chars=1.5 ,$
     align=0.5,noclip=0
   xyouts,density(nden-2),list_int(ilist,nden-2),string(ilist,'(i3)'),$
     chars=1.5,align=0.5,noclip=0
endfor
;
;
;  get a reference line
;
imax=where(list_int eq max(list_int))
imax = imax(0)
ref=imax-fix(imax/nlist)*nlist

oplot,density,list_int(ref,*),thick=2
;
;
;
refstr=strtrim(string(ref,'(i2)'),2)
wvlstr=strtrim(string(list_wvl(ref),'(f10.4)'),2)

ymin = 1.
for ilist=0,nlist-1 DO ymin = ymin < min(list_int(ilist,*)/list_int(ref,*))

ymin = ymin > intens
ymax = max(list_int(*,*)/list_int(ref,*)) >  2.

;
plot_oo,fltarr(2),xr=[density(0),density(nden-1)],yr=[0.9*ymin, 1.1*ymax],$
   xst=1, yst=1, $
  tit='Intensity Rel. to line '+refstr+' at '+wvlstr,$
  chars=0.8, xtitle='Density (cm-3)'

for ilist=0,nlist-1 do begin
   oplot,density,list_int(ilist,*)/list_int(ref,*)
;   oplot,density,list_int(ilist,*)/list_int(ref,*), psym=1
endfor
;
;
for ilist=0,nlist-1 do begin
   iden=isrt(ilist<(nden-1))
   xyouts,density(iden),list_int(ilist,iden)/list_int(ref,iden),$
     string(ilist,'(i3)'),chars=1.5,align=0.5,noclip=0
   xyouts,density(1),list_int(ilist,1)/list_int(ref,1),$
     string(ilist,'(i3)'),chars=1.5,align=0.5,noclip=0
   xyouts,density(nden-2),list_int(ilist,nden-2)/list_int(ref,nden-1),$
     string(ilist,'(i3)'),chars=1.5,align=0.5,noclip=0
endfor
;
end
