;+
; Project     : SOHO - CDS     
;                   
; Name        : PLOT_DMM_TR_FIG
;               
; Purpose     : Plots line intensities in DMM_TR
;               
; Explanation : Internal routine used by CHIANTI_TE for plotting individual
;               line intensities.
;               
; Use         : Called by CHIANTI_TE
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
; Modified    : Fix bug in searching for peak line (IMAX).  CDP, 27-Jan-96
;               V.1  H E Mason, 03-Oct-96
;               V.2 changed yranges in the plots. Removed common calls not
;               needed. reinstated the two plots. Changed various things.
;               Giulio Del Zanna (DAMTP), 10-Oct-2000
;
; Version     : Version 2, 
;
;
;-            

pro plot_dmm_tr_fig

common dmm_lines_te, list_wvl, list_int, list_descr1, list_descr2, species,$
                  temperature, ntem, ratio, description, nlist, savetext, $
                  dens, intens
;
intmax=max(list_int)
intmin= min(list_int) > intmax*intens

;
 plot_oo,fltarr(2),xr=[temperature(0),temperature(ntem-1)],yr=[intmin,intmax],$
        ytit='Rel. Int. (Pop. of upper level)*A ',$
       xtitle='Temperature (K)',chars=0.8
;

 for ilist=0,nlist-1 do begin
    oplot,temperature,list_int(ilist,*)
;    oplot,temperature,list_int(ilist,*), psym=1
 endfor

;
rand=randomu(seed,ntem)
isrt=sort(rand)   ; to try to fix problems with labeling curves
;
;
 for ilist=0,nlist-1 do begin
 item=isrt(ilist<(ntem-1))
 xyouts,temperature(item),list_int(ilist,item),string(ilist,'(i3)'),chars=1.3,$
            align=0.5,noclip=0
 xyouts,temperature(1),list_int(ilist,1),string(ilist,'(i3)'),chars=1.3 ,$
            align=0.5,noclip=0
 xyouts,temperature(ntem-2),list_int(ilist,ntem-2),string(ilist,'(i3)'),$
            chars=1.3,align=0.5,noclip=0
 endfor

;
;  get a reference line
;
imax=where(list_int eq max(list_int))
imax = imax(0)
ref=imax-fix(imax/nlist)*nlist

 oplot,temperature,list_int(ref,*),thick=2
;
;
refstr=strtrim(string(ref,'(i2)'),2)
wvlstr=strtrim(string(list_wvl(ref),'(f10.4)'),2)

ymin = 1.
for ilist=0,nlist-1 DO ymin = ymin < min(list_int(ilist,*)/list_int(ref,*))

ymin = ymin > intens
ymax = max(list_int(*,*)/list_int(ref,*)) >  2.


;
plot_oo,fltarr(2),xr=[temperature(0),temperature(ntem-1)],yr=[0.9*ymin, 1.1*ymax],$
      xst=1, yst=1, ytit='Intensity Rel. to line '+refstr+' at '+wvlstr,$
       xtitle='Temperature (K)',chars=0.8

for ilist=0,nlist-1 do begin
   oplot,temperature,list_int(ilist,*)/list_int(ref,*)
;   oplot,temperature,list_int(ilist,*)/list_int(ref,*), psym=1
endfor
;
;
for ilist=0,nlist-1 do begin
   iden=isrt(ilist<(ntem-1))
   xyouts,temperature(item),list_int(ilist,iden)/list_int(ref,iden),$
          string(ilist,'(i3)'),charsiz=1.3,align=0.5,noclip=0
   xyouts,temperature(1),list_int(ilist,1)/list_int(ref,1),$
          string(ilist,'(i3)'),charsiz=1.3,align=0.5,noclip=0
   xyouts,temperature(ntem-2),list_int(ilist,ntem-2)/list_int(ref,ntem-1),$
          string(ilist,'(i3)'),charsiz=1.3,align=0.5,noclip=0
endfor
;
end
