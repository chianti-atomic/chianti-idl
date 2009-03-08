;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory (USA), the University of Florence (Italy), the
;       University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
;
; NAME:
;	DENSITY_RATIOS
;
; PURPOSE:
;
;	to calculate line intensity ratios as a function of electron density
;
; CATEGORY:
;
;	scientific analysis
;
; CALLING SEQUENCE:
;
;       DENSITY_RATIOS,Ion,Wmin,Wmax,Dmin,Dmax,Density,Ratio,Description
;
;
; INPUTS:
;
;       Ion:   the CHIANTI style name of the ion, i.e., 'c_5' for C V
;	wmin:  minimum of the wavelength range of interest in Angstroms
;	wmax:  maximum of the wavelength range of interest in Angstroms
;       dmin:  log10 of the minimum desired density (8. = 10^8 cm^(-3) )
;       dmax:  log10 of the maximum desired density range
;
; INTERACTIVE INPUTS:
;
;	Must select the line for the numerator and denominator 
;       It is possible to select multiple lines to be summed
;	
; KEYWORD PARAMETERS:
;
;        OUTFILE:  the (optional) name of the output ascii file where a 
;                   listing of the line ratio intensity as a function of 
;                   density is saved.  For example, outfile='den_rat.lis'
;
;        PSFILE:  the (optional) name of the output postscript file 
;                  where a plot of the choses density sensitive line
;                  ratio is saved.  For example, psfile='den_rat.ps'
;
;        TEMP:   to specify the temperature, otherwise the temperature at 
;                  the peak 
;                  of the ionization equilibrium is used.  For example, 
;                  temp=1.e+6
;
;        /PHOTONS:  if set, the ratio will be in photon units, as opposed 
;                   to ergs    
;
;	RADTEMP	Specify background radiation temperature (default: 6000 K)
;
;	RPHOT   Distance from the centre of the star in stellar radius units.
;               I.e., RPHOT=1 corresponds to the star's surface. (Default is
;               infinity, i.e., no photoexcitation.)
;
;       NOPROT  If set, then the default setting will be NOT to use 
;               proton rates. This can be changed within the routine.
;
;       VERBOSE To print out information about the lines.
;
; OUTPUTS:
;
;	Density:  an array of the density values for which the selected 
;                   intensity ratio calculated 
;       Ratio:    an array of line intensity ratios
;       Description:  a string describing the transitions selected
;
;       Plots the intensity ratio of the selected line as a function of density
;
;
; COMMON BLOCKS:
;
;       None.
;
; CALLS
;
;       EMISS_CALC, ION2SPECTROSCOPIC, CONVERTNAME, READ_IONEQ,
;       CH_XMENU_SEL
;
; EXAMPLE:
;
;       density_ratios,'o_5',1000.,1500.,8.,13.,den,rat,desc 
;
;       choose the ratio of 1371.294 to 1218.393 
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	March 1996:     Version 2.0
;       May 28, 1996:  Ken Dere added psfile keyword/option
;       Sept 1996:     modified to work with VMS
;                      and added keyword TEMP, Ken Dere
;       Feb. 2000:     Modified for Version 3, K. Dere
;       14-Jul-2000    Peter Young, now calls pop_solver
;       26-Sep-2001    Modified for 9-point splines and proton rates; 
;                      added radtemp and rphot keywords for photoexcitation.
;       20-Nov-2001    Routine now calls emiss_calc to get emissivities.
;
;       V.9, 21-May-2002, Giulio Del Zanna (GDZ) 
;                      generalized directory concatenation to work for
;                               Unix, Windows  and VMS.
;       V.10, 1-Aug-02 GDZ
;          Changed all the formats.
;
;       V.11, 06-Aug-02 GDZ
;              Changed the use of CHIANTI system variables. 
;
;       V.12, 15-Aug-02, GDZ 
;              Major revision:
;              -Added the keyword VERBOSE, to avoid printing out long lists of lines.
;              -Removed the call to ch_xselect_s, that did not work for long lists.
;              -Added a '*' in the line lists, to identify 'unobserved' lines.
;              -Replaced the commands to create PS file, to make it
;               cross-platform compatible.
;              -Added a large number of cosmetics, mainly lables to the axes and
;               titles, that were missing.
;              -Removed plotting in windows already present.
;              -Added the CHIANTI version number in the outputs.
;
;        V.13, 3-Nov-03  GDZ
;           Modified format e8.2 to e9.2 for Windows compatibility.
;
; VERSION     :   13, 3-Nov-03
;
;-
pro density_ratios,ions,wmin,wmax,denmin,denmax,density,ratio,description, $
           outfile=outfile,psfile=psfile,TEMP=temp,PHOTONS=photons, $
           noprot=noprot, radtemp=radtemp, rphot=rphot, VERBOSE=VERBOSE

;
;  wmin=short wavelength limit in Angstroms
;  wmax=long wavelength limit in Angstroms
;
;  dmin= log 10 of minimum density to be considered
;  dmax= log 10 of maximum density to be considered
;
;
;
;
if n_params() lt 6 then begin
   print,'Use: IDL> density_ratios,ion,wmin,wmax,denmin,denmax,density,ratio,description'
   print,'                         [outfile=, psfile= ,temp=, /photons, $'
   print,'                          /noprot, radtemp=radtemp, rphot=rphot, /VERBOSE ]'
   print,'    i.e.> density_ratios,''c_5'',40.,50.,5.,15.,den,rat,desc'
   print,' '
   return
endif

angstrom = ' ('+string(197B)+')'


ion2spectroscopic,ions,species

convertname,ions,iz,ion
locname=strlowcase(ions)
pos=strpos(locname,'_')
l=strlen(pos)
first=strmid(locname,0,pos)
last=strmid(locname,pos+1,l-pos-1)
;
if strpos(last,'d') ge 0 then dielectronic=1 else dielectronic=0

dlogden=0.25                    ;  increment in logNe for calculating ratios
nden=fix((denmax-denmin)/dlogden) +1
density = 10.^(denmin+findgen(nden)*dlogden)


;
if keyword_set(temp) then begin
   tmax=temp
endif else begin

;read default !ioneq_file: 

   print,'Getting T from maximum of default ioniz. fraction'

   ioneq_name=!ioneq_file 
   read_ioneq,ioneq_name,ioneq_t,ioneq,ioneq_ref

   this_ioneq=ioneq(*,iz-1,ion-1+dielectronic)
   hit=where(this_ioneq eq max(this_ioneq))
   tmax=ioneq_t(hit)
   tmax=10.^(tmax(0))
endelse

print,'Calculating emissivities at T (K)='+string(tmax, format='(e10.3)')

IF keyword_set(photons) THEN no_de=1 ELSE no_de=0

em=emiss_calc(iz,ion,diel=dielectronic,temp=alog10(tmax), $
              dens=alog10(density), noprot=noprot, $
              /quiet,rphot=rphot,radtemp=radtemp,no_de=no_de)

FOR i=0,n_elements(em)-1 DO em[i].em=em[i].em/density

ind=where( (em.lambda GE wmin) AND (em.lambda LE wmax) )

list_int=em[ind].em[0,*]

nlist=n_elements(ind)

IF trim(em[0].version) EQ '' THEN version = '' ELSE $
  version = 'V. '+em[0].version+' '

;
xtitle='Electron Density (cm!e-3!n)'
ytitle='Relative Intensity (Population of upper level)*A/Density'
title = 'CHIANTI '+version+species+' emissivities at  T = '+string(tmax,'(e9.2)')+' (K)'


intmax=max(em[ind].em)
intmin=intmax/1d3
;
;device,window_state=ws
;if(ws(0) ne 1) then  window,0,ysize=425,xsize=525 else wset,0
window,0,ysize=425,xsize=525


plot,fltarr(2),xr=[density(0),density(nden-1)],yr=[intmin,intmax], $
  xtitle=xtitle,ytitle=ytitle, title=title, /xlog,/ylog
;
;
FOR ilist=0,nlist-1 DO BEGIN
   IF max(em[ind[ilist]].em) EQ intmax THEN BEGIN
      th=2
      ref=ilist
   ENDIF ELSE BEGIN
      th=1
   ENDELSE
   oplot,density,reform(em[ind[ilist]].em),th=th
ENDFOR
th=1

;
rand1=(nden-1)*randomu(seed,nlist) ; to help scatter the labels about the plot
rand2=(nden-1)*randomu(seed,nlist)
;
;
FOR ilist=0,nlist-1 DO BEGIN
   iden=rand1(ilist)
   xyouts,density(iden),em[ind[ilist]].em[0,iden], $
     string(ilist,'(i4)'), $
     size=charsiz,align=0.5,noclip=0
ENDFOR

;
;device,window_state=ws
;if(ws(2) ne 1) then  window,2,ysize=425,xsize=525 else wset,2
window,1,ysize=425,xsize=525

refstr=strtrim(string(ref,'(i4)'),2)
wvlstr=strtrim(string(em[ind[ref]].lambda,'(f10.4)'),2)

ytitle='Intensity Relative to line '+refstr+' at '+wvlstr+angstrom
;
plot,fltarr(2),xr=[density(0),density(nden-1)],yr=[1.e-2,1.e+1], $
  xtitle=xtitle,ytitle=ytitle, title=title,/xlog,/ylog

for ilist=0,nlist-1 do begin
   oplot,density,reform(em[ind[ilist]].em)/reform(em[ind[ref]].em)
endfor
;
;
FOR ilist=0,nlist-1 DO BEGIN
   iden=rand2(ilist)
   xyouts,density(iden),em[ind[ilist]].em[0,iden]/em[ind[ref]].em[0,iden], $
     string(ilist,'(i4)') $
     ,charsiz=charsiz,align=0.5,noclip=0
ENDFOR
;

options=strarr(nlist)

print,' ------------------------------------------------------'

FOR i=0,nlist-1 DO BEGIN

   IF em[ind[i]].flag EQ -1 THEN text_add = ' *' ELSE text_add ='  '

   options[i]=string(i, format='(i4)')+ $
     string(em[ind[i]].lambda,format='(f10.4)')+text_add+angstrom+$
     '  '+string(em[ind[i]].em[0,nden-1],format='(e10.2)')+'  '+ $
     em[ind[i]].lvl1_desc+' - '+em[ind[i]].lvl2_desc

   IF keyword_set(VERBOSE) THEN   print,options[i]

ENDFOR

numer_index = ch_xmenu_sel(options, tit=' Select lines for numerator', $
                           text=['If more than one line is selected, ',$
                                 ' the emissivities of the lines will be summed.', $
                                 'A  * indicates that the wavelength is theoretical' ])

IF numer_index[0] EQ -1 THEN BEGIN 
   print, '% DENSITY_RATIOS: Routine Aborted'
   return
END 

denom_index = ch_xmenu_sel(options, tit=' Select lines for denominator', $
                           text=['If more than one line is selected, ',$
                                 ' the emissivities of the lines will be summed.'] )

IF denom_index[0] EQ -1 THEN BEGIN 
   print, '% DENSITY_RATIOS: Routine Aborted'
   return
END 

;ch_xselect_s,options,nstatus,abort,title='Select lines for numerator'
;IF (abort EQ 1) OR (total(nstatus) EQ 0) THEN return
;ch_xselect_s,options,dstatus,abort,title='Select lines for denominator'
;IF (abort EQ 1) OR (total(dstatus) EQ 0) THEN return
;i_n=ind[where(nstatus EQ 1)]
;i_d=ind[where(dstatus EQ 1)]

i_n=ind[numer_index]
i_d=ind[denom_index]

numerator=fltarr(nden)

n_n=n_elements(i_n)
n_d=n_elements(i_d)

IF n_n GT 1 THEN numerator=total(reform(em[i_n].em),2) $
ELSE numerator=reform(em[i_n].em)
IF n_d GT 1 THEN denominator=total(reform(em[i_d].em),2) $
ELSE denominator=reform(em[i_d].em)

ratio=numerator/denominator

print,'-------------------------------------------------------'

description='CHIANTI '+version+species+' '
desc_ascii = species+' '

FOR i=0,n_n-1 DO BEGIN
   description=description+strtrim(string(em[i_n[i]].lambda,'(f10.4)'),2)
   IF i NE n_n-1 THEN description=description+'+'
ENDFOR
desc_ascii = desc_ascii+description+' (A) /'
description=description+angstrom+'/'
FOR i=0,n_d-1 DO BEGIN
   desc_ascii = desc_ascii+strtrim(string(em[i_d[i]].lambda,'(f10.4)'),2)
   description=description+strtrim(string(em[i_d[i]].lambda,'(f10.4)'),2)
   IF i NE n_d-1 THEN BEGIN 
      description=description+'+'
      desc_ascii = desc_ascii+'+'
   END 
ENDFOR

desc_ascii = desc_ascii+' (A)'

description=description+angstrom+'  T = '+string(tmax,'(e9.2)')+' (K)'


;device,window_state=ws
;if(ws(1) ne 1) then  window,1,ysize=425,xsize=525 else wset,1
window,2,ysize=425,xsize=525

if keyword_set(photons) then begin
   ytitle='Intensity (Photons) Ratio'
endif else begin
   ytitle='Intensity (Ergs) Ratio'
endelse


plot,density,ratio,xr=[density(0),density(nden-1)],title=description, $
  xtitle=xtitle,ytitle=ytitle,/xlog,/ylog

;
IF  keyword_set(psfile) THEN  BEGIN 
   ps, psfile, /LANDSCAPE
   plot_oo,density,ratio,xr=[density(0),density(nden-1)],title=description, $
     xthick=thick,ythick=thick,thick=thick,xtitle=xtitle, ytitle=ytitle
   psclose
ENDIF 


;
if keyword_set(outfile) then begin
   openw,luo,outfile,/get_lun
;
   printf,luo, outfile
   printf,luo, ytitle+' calculated with CHIANTI '+version
   printf,luo,desc_ascii

   printf,luo,'at constant Temperature='+string(tmax,'(e9.2)')+' (K)'

   printf,luo, ' '
   printf,luo, 'Density (cm-3) ratio '
   printf,luo, ' '

   for i=0,nden-1 do begin
      printf,luo,density(i),ratio(i)
   endfor
;
   close,luo
;
endif

end

