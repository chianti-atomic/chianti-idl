;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory (USA), the University of Florence (Italy), the
;       University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
; NAME:
;	synthetic_plot
;
; PURPOSE:
;
;       to plot out synthetic spectra calculated with Synthetic
;       and interactively identify spectral lines
;
;
; CATEGORY:
;	
;	spectroscopy
;
; CALLING SEQUENCE:
;
;       SYNTHETIC_PLOT,Wvl,Spectrum,List_wvl,List_ident,fwhm
;
;
; INPUTS:
;
;       Wvl:  wavelength array from synthetic
;       Spectrum:  spectrum intensity array from synthetic
;       List_wvl:  string array of spectral line wavelengths
;       List_ident:  string array of spectral line identifications
;       Fwhm:  when the cursor is clicked, spectral lines with fwhm
;              (Angstroms) of the cursor are printed out
;
;
; KEYWORDS
;
;	xrange:  similar to IDL keyword to determine wavelength range of plot
;
; OUTPUTS:
;
;       None
;
;
; PROCEDURE:
;
;	Click the left mouse button to select a wavelength
;       Click the right mouse button to exit
;
; EXAMPLE:
;
;      > synthetic,100.,200.,.1,1.e+15,wvl,spectrum,list_wvl,list_ident
;      > synthetic_plot,wvl,spectrum,list_wvl,list_ident,0.1
;                     
;      note:  it is not necessary for the two fwhm values to be the same      
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	May 1996:     Version 2.0
;       Dec. 1998:    revised by Ken Dere
;       V.4,  23 Oct 2000 GDZ, added the log keyword, and changed a few things
;       in the plot. 
;
;       Ver.5, 12-Dec-2001, Peter Young
;           Changed style of printing, and made method of extracting the 
;           intensity from list_ident compatible with the new version of 
;           isothermal.pro.
;
; VERSION 5   12 Dec 2001 Peter Young
;-
pro synthetic_plot,lambda,spectrum,list_wvl,list_ident,fwhm_a,xrange=xrange, log=log

;
if n_params(0) lt 5 then begin
;
print,'use IDL> synthetic_plot,lambda,spectrum,list_wvl,list_ident,delta_wvl [,xrange]'
return
endif
;
maxi=max(spectrum)

IF keyword_set(log) THEN mini = maxi/10. ELSE mini=0.  ;maxi/10.

IF n_elements(xrange) EQ 2 THEN xrange = xrange ELSE xrange = [min(lambda), max(lambda)]

IF keyword_set(log) THEN $
   plot_io,lambda,spectrum,xrange=xrange, xst=1, /yno ELSE $
 plot,lambda,spectrum, xrange=xrange, xst=1, /yno

;
delta_lambda=fwhm_a  ;search for lines within delta_lambda
;
coord=2   ; data coordinates

;
;  derived from IDL userlib routine rdpix.pro
;
on_error,2              ;Return to caller if an error occurs
;
print,'Press left or center mouse button to display wvl and identification'
print,'... right mouse button to exit.'
;
!err=0
;
cr = string("15b)	;this codes a newline
;
;
form="($,'x=',f12.3,', y=',e12.2,'    data coordinates',a)"
;
;
while !err ne 4 do begin

  cursor,x,y,2,/data

  if (!err and 3) ne 0 then begin ;New line?
    print,form="($,a)",string("12b)
    print,format='("Selected wavelength: ",f10.3)',x
;
;
    g=where(abs(x-list_wvl) le delta_lambda)
    ng=n_elements(g)
    if max(g) ge 0 then begin
      contr_int=fltarr(n_elements(list_wvl))
      for i=0,ng-1 do BEGIN
        bits=str_sep(list_ident[g[i]],'Int=')
        reads,bits[1],format='(e10.0)',int
        contr_int[g[i]]=int
;         i1=strpos(list_ident(g(i)),'Int=')
;         i2=strpos(list_ident(g(i)),'Tmax')
;         str_int=strmid(list_ident(g(i)),i1+4,i2-i1-5)
;         contr_int(g(i))=float(str_int)
      endfor
      total_int=total(contr_int(g))
      for i=0,ng-1 do begin
        if(contr_int(g(i)) gt 0.01*total_int) then BEGIN
          print_str=string(format='(f10.3)',list_wvl[g[i]])+ $
               '  '+list_ident[g[i]]
          print,print_str
;          print,list_wvl(g(i)),'  '+list_ident(g(i)) ;,string(contr_int(g(i)),'(f15.2)')
        endif
      endfor
    endif

    while (!err ne 0) do begin 
      wait,.1
      cursor,x,y,2,/data
    endwhile
;
  endif

  print,form = form, x,y,cr
;
endwhile
print,form="(/)"
end


