
FUNCTION ch_plot_cit_year, xx, yy, _extra=extra, no_scale=no_scale, reverse=reverse


;+
; NAME:
;     CH_PLOT_CIT_YEAR
;
; PURPOSE:
;     Creates a plot showing the numbers of citations per year the
;     CHIANTI papers have received. This plot is shown on the CHIANTI
;     citations webpage.
;
; CATEGORY:
;     CHIANTI; citations; plot.
;
; CALLING SEQUENCE:
;     Result = CH_PLOT_CIT_YEAR( XX, YY )
;
; INPUTS:
;     Xx:    A 1D array containing a list of years.
;     Yy:    A 1D array containing the citations for the years in XX. 
;
; OPTIONAL INPUTS:
;     _Extra:  Allows additional IDL plot object commands to be passed
;              on to the routine.
;
; KEYWORD PARAMETERS:
;     NO_SCALE: By default, the routine scales upwards the citations
;               for the current year based on the number of days that
;               have passed. Setting this keyword switches off the
;               scaling. 
;     REVERSE:  If set, then the plot will be white-on-black rather
;               than black-on-white.
;
; OUTPUTS:
;     Returns an IDL plot object containing the plot of citations
;     vs. year. 
;
; EXAMPLE:
;     IDL> w=cit_plot_cit_year(xx,yy,/buffer)
;     IDL> w.save,'plot_cit_year.png',resolution=96
;
; MODIFICATION HISTORY:
;     Ver.1, 17-Sep-2019, Peter Young
;-


IF keyword_set(reverse) THEN BEGIN
  bg_color='black'
  color='white'
ENDIF ELSE BEGIN
  bg_color='white'
  color='black'
ENDELSE 

;
; Get rid of any years with zero citations.
; 
k=where(yy NE 0)
yr=xx[k]
ncit=yy[k]

;
; For the current year I scale the number of citations based on the
; number of elapsed days in the year. I check the last modification
; time of the bib fle to get the elapsed number of days.
;
IF NOT keyword_set(no_scale) THEN BEGIN
  t=systime(/jul,/utc)
  caldat,t,mm,dd,year
  k=where(year EQ yr,nk)
  IF nk NE 0 THEN BEGIN
    mjd=t-2400000.5d
    mjd_str={ mjd: floor(mjd), time: (mjd-floor(mjd))*8.64d7 }
    ndays=doy(mjd_str)
    ncit[k]=round(float(ncit[k])*365./float(ndays))
  ENDIF 
ENDIF 



th=2
IF n_elements(font_size) EQ 0 THEN fs=14 ELSE fs=font_size
dimensions=[700,450]
w=window(dim=dimensions,_extra=extra)
w.background_color=bg_color

;
; I wanted to color the current year red, but this didn't work. 
;
color_string=strarr(n_elements(yr))
color_string[*]='blue'
color_string[-1]='red'

yrange=[0,ceil( (max(ncit)+1)/10. )*10.]

p=barplot(yr,ncit, $
          ytit='No. of citations per year', $
          pos=[0.13,0.11,0.96,0.98], $
          thick=th, xthick=th,ythick=th,color=color, $
          xtext_color=color,ytext_color=color, $
          outline=0, $
          xcolor=color,ycolor=color, $
          font_size=fs, $
          xticklen=0.015, yticklen=0.015, $
          xrange=[min(yr)-1,max(yr)+1],xsty=1,  $
          yrange=yrange,ysty=1, $
          xtitle='Year',/current, $
          _extra=extra)

return,w

END
