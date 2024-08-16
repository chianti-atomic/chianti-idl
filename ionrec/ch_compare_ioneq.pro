;+
;
; PROJECT:  CHIANTI
;
;        This program was developed as part of CHIANTI-VIP. 
;        CHIANTI-VIP (Virtual IDL and Python) is  a member of the CHIANTI family
;        mantained by Giulio Del Zanna, to develop additional features and
;        provide them to the astrophysics community.
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the
;       University of Cambridge, Goddard Space Flight Center, and University of Michigan.
;
; Name        : CH_COMPARE_IONEQ
;               
; Purpose     : Compares the ionisation equilibrium values for an element, present 
;               in different files.
;               NOTE:  the temperature grids of the files do not need to be the same.
;                      For the initial plot, the first file is taken as reference.
;
; Explanation : 
;              The first file name is taken as reference and plotted in black by default,
;              the others are overplotted with different colours.
;              The axes of the plot can be modified and a postscript file can be created.
;
; Restrictions: Only CHIANTI-format ioneq files can be read.
;               If colours are not supplied up to 8 files can be plotted.
;
; Example     : IDL> ch_compare_ioneq, 'C', files=['test.ioneq','chianti.ioneq'], lines=[0,2]
;    
; Inputs      : ELEMENT - the element name (lower or upper case) to be plotted.
;               FILES - the list of ioneq file names
;
; Opt. Inputs : ION_RANGE - specify range of ions to be plotted, e.g., ion=[2,4]
;
;               LABELS - the labels for the legend
;
;               COLORS - the colors numbers:
;                 black=0 & maroon=1 & yellow=5 & red = 2 & green = 7
;                 blue = 10 & orange=4 & purple=13 & magenta=12
;                 pink=3 & olive=6 & dgreen=8 & cyan=9 & dblue=11
;               
;               LINESTYLES - the line styles
;
;               TITLE - the title for the plot
;
;               /INTERPOLATE - to interpolate the data on a fine grid, 0.01 in log T
;
;               CHAR_LEGEND - the charsize of the legend.
; 
;               Some optional inputs can be passed to  AL_LEGEND, the
;               most useful ones are e.g.  /top,/right for the placement of the
;               legend
;
; Opt. Outputs: a postscript file.
;
; Calls:    READ_IONEQ, PRINT2D_PLOT, AL_LEGEND, LINECOLORS
;
; Prev. Hist. :  None
;
; Written     : Giulio Del Zanna (GDZ), 16 Sept 2023
;
; Modified    :
;
; VERSION     : 1, 16 Sept 2023 
;
;- 

PRO ch_compare_ioneq,  element , ion_range=ion , files=files, $
                       labels=labels, colors=colors, linestyles=linestyles,$
                       title=title, interpolate=interpolate,CHAR_legend=CHAR_legend,$
                       psyms=psyms,  _extra=_extra


  if n_params() eq 0 then begin
     print,'Use:  IDL> ch_compare_ioneq, element, files=files  '
     return
  endif

  elements=['H','He','Li','Be','B','C','N','O','F','Ne','Na',$
            'Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti',$
            'V','Cr','Mn','Fe','Co','Ni','Cu','Zn']


  iz = where(strlowcase(elements) eq strlowcase(element))
  if iz[0] eq -1 then begin
     bell
     print,'ERROR, unrecognised element' 
     return
  endif else begin
     iz = iz[0]
     element=elements[iz]       ; to avoid lowercase problems.
  end 
  
  loadct,0
  !p.background=255
  !p.color=0
  
  ionstage=['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII',$
            'XIII','XIV','XV','XVI','XVII','XVIII','XIX','XX','XXI',' XXII',$
            'XXIII','XXIV','XXV','XXVI','XXVII','XXVIII','XXIX']

;
;  which ions to plot
;

  if n_elements(ion) gt 0 then begin
     min_ion = min(ion)-1
     max_ion = max(ion)
  endif else begin
     min_ion = 0
     max_ion = iz+2 
  endelse

  nf=n_elements(files)

  if n_elements(title) eq 0 then title='CHIANTI ionization equilibrium for '+element

  if n_elements(charsize) eq 0 then charsize = 2
  if n_elements(char_legend) eq 0 then char_legend = 1.5

  if n_elements(files) le 1 then $
     message,'ERROR - at least two ioneq files need to be given as input'

  if  n_elements(labels) eq 0 then labels=trim(1+indgen(nf))

  if n_elements(linestyles) eq 0 then linestyles=indgen(nf)

  if n_elements(psyms) eq 0 then psyms=intarr(nf) ; do not plot
  
  linecolors

  black=0 & maroon=1 & yellow=5 & red = 2 & green = 7
  blue = 10 & orange=4 & purple=13 & magenta=12
  pink=3 & olive=6 & dgreen=8 & cyan=9 & dblue=11

; if colors are not defined, define them:  
  these_col=[black, blue, red, purple, olive, orange, green, magenta]
; allow maximum of 8
  if n_elements(colors) ne nf then colors=these_col[0:nf-1]


; read the first one.
  
  ioneq_name=files[0]
  read_ioneq, ioneq_name, logtemp1,data1n,  ref
  nt=n_elements(logtemp1)
  
  if keyword_set(interpolate) then begin 
; now interpolate over a fine grid ? 
     step=0.01
     n_log_t_step=fix( (max(logtemp1)-min(logtemp1))/step) +1
     log_t_step=min(logtemp1)+findgen(n_log_t_step)*step 

     data1=fltarr(n_log_t_step,iz+2)
     for jj=0, iz+1 do data1[*, jj]=interpol(reform(data1n[*,iz, jj]), $
                                             logtemp1, log_t_step,/spline)

     logt1=log_t_step
  endif else begin
     logt1=logtemp1
     data1=reform(data1n[*,iz, 0:iz+1])
  end 

  circle_sym,/fill

  x_min=min(logtemp1) 
  x_max=max(logtemp1)
  y_min=0.
  y_max=1.4


  pr=' '
begin_post:

  
  if pr eq ' ' then begin
     th=1
     !p.background=255
     !p.color=0
  endif else th=3

  
  plot, logt1 , data1[*, min_ion],$
        xrange=[x_min,x_max], yrange=[y_min,y_max],$
        xstyle=1,ystyle=1,$
        title=title, th=th,charth=th,xth=th,yth=th,$
        /nodata, chars=charsize,xtit='log T [K]'

  oplot, logt1 , data1[*, min_ion], line=linestyles[0],col=colors[0], th=th,psym=-psyms[0]

  hit = where(data1[*, min_ion] eq max(data1[*, min_ion]))
  xyouts,  logt1[hit[0]]+0.02,max(data1[*, min_ion])+0.02,'!17'+ionstage[min_ion]+'!3',$
           chars=1.5,charth=th

  if max_ion ne min_ion then begin
     for i = min_ion+1, max_ion-1 do begin

        oplot, logt1, data1[*, i], line=linestyles[0], col=colors[0],th=th,psym=-psyms[0]

        hit = where(data1[*, i] eq max(data1[*,i]))
        if max(data1[*, i]) gt 0 then BEGIN 
           t = logt1(last_item(hit)) + 0.02
           if t lt x_max and t gt x_min then begin
              xyouts,t,max(data1[*,i])+0.02,'!17'+ionstage[i]+'!3',chars=1.5,charth=th
           endif
        endif
     endfor
  endif 

  
  for ii=1,n_elements(files)-1 do begin 

     ioneq_name=files[ii]
     read_ioneq, ioneq_name, logtemp,datan,  ref

     if keyword_set(interpolate) then begin 
; now interpolate over a fine grid 
        step=0.01
        n_log_t_step=fix( (max(logtemp)-min(logtemp))/step) +1
        log_t_step=min(logtemp)+findgen(n_log_t_step)*step 
        logt=log_t_step

        data=fltarr(n_log_t_step,iz+2)
        for jj=0, iz+1 do data[*, jj]=interpol(reform(datan[*,iz, jj]), $
                                               logtemp, log_t_step,/spline)

     endif else begin
        logt=logtemp
        data=reform(datan[*,iz, 0:iz+1])
     end

     oplot, logt , data[*, min_ion], line=linestyles[ii],col=colors[ii], th=th,psym=-psyms[ii]

     if max_ion ne min_ion then begin
        for i = min_ion+1, max_ion-1 do begin

           oplot, logt, data[*, i], line=linestyles[ii], col=colors[ii],th=th,psym=-psyms[ii]
        endfor
     endif  

     
  endfor   
  
  al_legend, CHARSIZE =CHAR_legend,CHARTHICK=th, THICK=th, _extra=_extra,$
             labels, LINESTYLE=LINESTYLES, psym=psyms, colors=colors, textcolors=colors ; /top, /right 


  print2d_plot, pr=pr , x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,$
                go_to_line=go_to_line, out_name=out_name , /ask_name

  if go_to_line eq 'y' then goto,begin_post

  print, 'NOTE: to fix an upside down landscape postscript plot, run e.g.:'
  print," cgFixPS, 'plot.ps','plot2.ps' "
  

END 
