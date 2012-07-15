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
;	TEMPERATURE_RATIOS
;
; PURPOSE:
;
;	calculate and display temperature sensitivity of line intensity ratios
;
; CATEGORY:
;
;	spectral diagnostics
;
; PROCEDURE :
;
;       The intensities (Population of the upper level * A)  of the lines within
;       the selected ion are first  calculated, either at constant pressure or
;       at constant density (however specified in the input). They are plotted
;       in window 0. The intensities  relative to the brightest  reference line
;       are then plotted in window 1. A widget allows the user to select a
;       number of lines (at least one!) for the numerator of the ratio, and a
;       number of lines for the denominator. In case of multiple selections, the
;       line intensities are summed. The ratio values are plotted in window 2, and
;       optionally also saved in a text file. A postscript file can also be
;       created. The ratio values, calculated at twice and half the prescribed
;       density  are also calculated and overplotted, to show how the
;       temperature ratio also depends on the density.
;
; CALLING SEQUENCE:
;
;
;      > temperature_ratios,ion,wmin,wmax,Log10(tempmin),Log10(tempmax),temperature,ratio,description,$
;       [pressure= ,density= , psfile= , outfile= ]
;
;
; EXAMPLE:
;
;      > temperature_ratios,'c_4',100.,1600.,4.,6.,temp,ratio,desc,density=1.e+10,$ 
;        psfile='test.ps', outfile='test.txt'
;
;       then select  the  ratio of (384.17 + 384.19) to 1550.772
;
;
;
; INPUTS:
;
;       Ion:   the CHIANTI style name of the ion, i.e., 'c_5' for C V
;
;       Wmin:  minimum  wavelength  limit in Angstroms
;
;       Wmax:  maximum  wavelength  limit in Angstroms
;
;       Tempmin:  log10 of lowest temperature of interest, i.e. 4 for 10.^4 K
;
;       Tempmax:  log10 of highest temperature of interest
;
; OPTIONAL INPUTS:
;
;	Must specify indices of lines which are to form the ratio
;	
;       RADTEMP   The blackbody radiation field temperature (default 
;                 6000 K).
;
;       RPHOT     Distance from the centre of the star in stellar radius 
;                 units. I.e., RPHOT=1 corresponds to the star's surface. 
;                 (Default is infinity, i.e., no photoexcitation.)
;
;       ABUND_FILE  The name of a CHIANTI abundance file. This is used for 
;                 calculating the proton to electron ratio. Default is 
;                 !abund_file.
;
;       IONEQ_FILE  The name of a CHIANTI ion balance file. This is used for 
;                 calculating the proton to electron ratio and evaluating 
;                 the T_max of the ion. Default is !ioneq_file.
;
; OUTPUTS:
;
;	Temperature:  an array of temperatures spanning Tempmin to Tempmax
;
;       Ratio:  an array of the intensity ratio of the selected lines
;
;       Desc:  a short string description of the selected line ratio
;
;
; OPTIONAL OUTPUTS:
;
;	Ps and/or  text file with the  intensity ratio.
;
;
; KEYWORD PARAMETERS:
;
;		
;	DENSITY:  calculates the intensity ratios for constant density.
;
;                  If neither density or pressure are specified, a constant
;                  density of 1.e+10 cm^-3 is assumed as default.
;
;       OUTFILE:  the (optional) name of the output ascii file where a 
;                   listing of the line ratio intensities as a function of 
;                   temperature is saved.
;
;       PSFILE:  the (optional) name of the output postscript file 
;                  where a plot of the chosen temperature sensitive line
;                  ratio is saved.
;
;       NOPROT    Switches off inclusion of proton rates.
;
;       VERBOSE   prints out information
;
; CALLS: 
;
;        read_ioneq, convertname, ion2spectroscopic,ion2filename,
;        ch_xmenu_sel, emiss_calc
;
; COMMON:
;
;       elvlc,l1a,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref
;       wgfa, wvl,gf,a_value
;       upsilon,t_type,deu,c_ups,splups
;       proton, pstr, pe_ratio
;       radiative, radt, dilute
;
; RESTRICTIONS:
;
;               
; SIDE EFFECTS: None known yet.
;                 
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	May 1996:     Version 2.0, Ken Dere
;       April 2000:   V. 3 Ken Dere modified for V3
;       14-Jul-2000   V. 4  Peter Young, now calls pop_solver
;        2-Oct-2000   V. 5 Giulio Del Zanna, corrected an error in the 
;             creation of the string list of the lines in the ratio. 
;             Also corrected a few minor errors.
;             Removed the device,window_state call, and added a few 
;             other 'cosmetic' adjustments.  
;	Version 6, 21-Dec-2000, William Thompson
;	      Modified for better cross-platform capability.
;
;       Ver.7, 6-Dec-2001, Peter Young
;             Revised to call emiss_calc for the emissivities.
;
;       V.8, 21-May-2002, GDZ 
;                   generalized directory concatenation to work for
;                   Unix, Windows  and VMS. 
;
;       V.9, 1-Aug-02 GDZ
;          Fixed label mistake, and changed all the formats.
;
;       V.10, 06-Aug-02 GDZ
;              Changed the use of CHIANTI system variables. 
;
;       V.11, 15-Aug-02, GDZ 
;              Major revision:
;              -Added the keyword VERBOSE, to avoid printing out long lists of lines.
;              -Removed the call to ch_xselect_s, that did not work for long lists.
;              -Added a '*' in the line lists, to identify 'unobserved' lines.
;              -Replaced the commands to create PS file, to make it
;               cross-platform compatible.
;              -Added a large number of cosmetics, mainly lables to the axes and
;               titles, that were missing.
;              -Removed plotting in windows already present.
;              -Removed the pressure keyword.
;              -Added the CHIANTI version number in the outputs.
;
;        V.12, 3-Nov-03  GDZ
;           Modified format e8.2 to e9.2 for Windows compatibility.
;
;        V.13, 6-Feb-2012, Peter Young
;           Changed all for loops to long integers (due to crashes
;           with big ions in IDL 7); also increased the index format
;           in the line selection widget from i4 to i6.
;
; VERSION     :   13, 6-Feb-12
;
;-
pro temperature_ratios,ions,wmin,wmax,tempmin,tempmax,$
               temperature,ratio,description,$
               density=density,psfile=psfile, $
               outfile=outfile,noprot=noprot, $
               radtemp=radtemp,rphot=rphot,photons=photons, $
               ioneq_file=ioneq_file, abund_file=abund_file, $
               VERBOSE=VERBOSE

common elvlc,l1a,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref
common wgfa, wvl,gf,a_value
common upsilon, splstr
COMMON proton, pstr, pe_ratio
COMMON radiative, radt, dilute
;
if n_params(0) lt 8 then begin
   print,' IDL> temperature_ratios,ion,wmin,wmax,Log10(tempmin),Log10(tempmax),$ '
   print,'           temperature,ratio,description,$ '
   print,'           [pressure= ,density= , psfile= , outfile= , radtemp= ,'
   print,'            rphot= , /noprot, ioneq_file= , abund_file=, /VERBOSE ]'
   print,' '
   print,' i.e.: '
   print,' IDL>  temperature_ratios, "c_5" ,40.,50.,5.,7.,temp,rat,desc '
   return
endif

angstrom = ' ('+string(197B)+')'


IF n_elements(density) EQ   0 THEN  density=1.e+10
dens = alog10(density) 

ion2spectroscopic,ions,species

charsiz=1.5


IF keyword_set(ioneq_file) THEN BEGIN
   ioneq_name=ioneq_file
ENDIF ELSE BEGIN
;ioneqdir=concat_dir(!xuvtop, 'ioneq')
;ioneq_name=concat_dir(ioneqdir, !ioneq_file)
   ioneq_name= !ioneq_file
ENDELSE

IF keyword_set(VERBOSE) THEN $
  print,' using '+ioneq_name+ ' as ionization equilibrium file'

read_ioneq,ioneq_name,ioneq_t,ioneq,ioneq_ref
;
convertname,ions,iz,ion
this_ioneq=ioneq(*,iz-1,ion-1)
hit=where(this_ioneq eq max(this_ioneq))
tmax=ioneq_t(hit)
tmax=10.^(tmax(0))

IF keyword_set(VERBOSE) THEN $
  print,'Maximum Temperature (ionization equilibrium) = ',tmax, $
  format='(a45,e10.3)'

;
dlogtemp=0.10                   ;  increment in logT for calculating ratios
;
ntemp=fix((tempmax-tempmin)/dlogtemp)+1
;


temperature=10.^(findgen(ntemp)*dlogtemp +tempmin)

print,  'Calculating '+species+'  emissivities at Ne = '+$
  string(density,'(e9.2)')+' (cm-3)'


em=emiss_calc(iz,ion,temp=alog10(temperature),dens=dens, $
              pressure=pressure,noprot=noprot,radtemp=radtemp, $
              rphot=rphot,/quiet,no_de=photons, ioneq_file=ioneq_name, $
              abund_file=abund_file)

ind=where( (em.lambda GE wmin) AND (em.lambda LE wmax) )

list_int=em[ind].em[0,*]

nlist=n_elements(ind)

IF trim(em[0].version) EQ '' THEN version = '' ELSE $
  version = 'V. '+em[0].version+' '


;;---------------------------------
;;Plot line intensities for all lines in wavelength range.
;;

xtitle='Electron Temperature (K)'
ytitle='Relative Intensity (Population of upper level)*A/Density'
title = 'CHIANTI '+version+species+' emissivities at Ne = '+$
  string(density,'(e9.2)')+' (cm!e-3!n)'


intmax=max(em[ind].em)
intmin=intmax/1d3
;
;device,window_state=ws
;if(ws(0) ne 1) then  window,0,ysize=425,xsize=525 else wset,0
window,0,ysize=425,xsize=525

plot,fltarr(2),xr=[temperature[0],temperature[ntemp-1]],yr=[intmin,intmax], $
  xtitle=xtitle,ytitle=ytitle,title=title, /xlog,/ylog
;
;
FOR ilist=0l,nlist-1 DO BEGIN
   IF max(em[ind[ilist]].em) EQ intmax THEN BEGIN
      th=2
      ref=ilist
   ENDIF ELSE BEGIN
      th=1
   ENDELSE
   oplot,temperature,reform(em[ind[ilist]].em),th=th
ENDFOR
th=1

;
rand1=(ntemp-1)*randomu(seed,nlist) ; to help scatter the labels about the plot
rand2=(ntemp-1)*randomu(seed,nlist)
;
;
FOR ilist=0l,nlist-1 DO BEGIN
   iden=rand1(ilist)
   xyouts,temperature(iden),em[ind[ilist]].em[iden], $
     string(ilist,'(i4)'), $
     size=charsiz,align=0.5,noclip=0
ENDFOR


;;-----------------------------------
;;Plot ratios relative to strongest line.
;;
;device,window_state=ws
;if(ws(2) ne 1) then  window,2,ysize=425,xsize=525 else wset,2
window,1,ysize=425,xsize=525

;
refstr=strtrim(string(ref,'(i4)'),2)
wvlstr=strtrim(string(em[ind[ref]].lambda,'(f10.4)'),2)
ytitle='Intensity Relative to line '+refstr+' at '+wvlstr+angstrom
;
plot,fltarr(2),xr=[temperature(0),temperature(ntemp-1)],yr=[1.e-3,1.e+1], $
  xtitle=xtitle,ytitle=ytitle,title=title, /xlog,/ylog

for ilist=0l,nlist-1 do begin
   oplot,temperature,reform(em[ind[ilist]].em)/reform(em[ind[ref]].em)
endfor
;
;
FOR ilist=0l,nlist-1 DO BEGIN
   iden=rand2(ilist)
   xyouts,temperature(iden),em[ind[ilist]].em[iden]/em[ind[ref]].em[iden], $
     string(ilist,'(i4)') $
     ,charsiz=charsiz,align=0.5,noclip=0
ENDFOR
;


;;-------------------------------
;;The following allows the user to select a numerator and denominator 
;;for the line ratio using widgets.
;;
options=strarr(nlist)
;
print,' ------------------------------------------------------'


FOR i=0l,nlist-1 DO BEGIN

   IF em[ind[i]].flag EQ -1 THEN text_add = ' *' ELSE text_add ='  '

   options[i]=string(i, format='(i6)')+ $
     string(em[ind[i]].lambda,format='(f10.4)')+text_add+angstrom+$
     '  '+string(em[ind[i]].em[ntemp-1],format='(e10.2)')+'  '+ $
     em[ind[i]].lvl1_desc+' - '+em[ind[i]].lvl2_desc

   IF keyword_set(VERBOSE) THEN   print,options[i]

ENDFOR

numer_index = ch_xmenu_sel(options, tit=' Select lines for numerator', $
                           text=['If more than one line is selected, ',$
                                 ' the emissivities of the lines will be summed.', $
                                 'A  * indicates that the wavelength is theoretical' ])

IF numer_index[0] EQ -1 THEN BEGIN 
   print, '% TEMPERATURE_RATIOS: Routine Aborted'
   return
END 

denom_index = ch_xmenu_sel(options, tit=' Select lines for denominator', $
                           text=['If more than one line is selected, ',$
                                 ' the emissivities of the lines will be summed.'] )

IF denom_index[0] EQ -1 THEN BEGIN 
   print, '% TEMPERATURE_RATIOS: Routine Aborted'
   return
END 

i_n=ind[numer_index]
i_d=ind[denom_index]


numerator=fltarr(ntemp)

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

FOR i=0l,n_n-1 DO BEGIN
   description=description+strtrim(string(em[i_n[i]].lambda,'(f10.4)'),2)
   IF i NE n_n-1 THEN description=description+'+'
ENDFOR
desc_ascii = desc_ascii+description+' (A) /'
description=description+angstrom+'/'
FOR i=0l,n_d-1 DO BEGIN
   desc_ascii = desc_ascii+strtrim(string(em[i_d[i]].lambda,'(f10.4)'),2)
   description=description+strtrim(string(em[i_d[i]].lambda,'(f10.4)'),2)
   IF i NE n_d-1 THEN BEGIN 
      description=description+'+'
      desc_ascii = desc_ascii+'+'
   END 
ENDFOR

desc_ascii = desc_ascii+' (A)'

description=description+angstrom+' Ne = '+$
  string(density,'(e9.2)')+' (cm!e-3!n)'


;;---------------------------------
;; Plot the line ratio

;;device,window_state=ws
;if(ws(1) ne 1) then  window,1,ysize=425,xsize=525 else wset,1
window,2,ysize=425,xsize=525

;
if keyword_set(photons) then begin
   ytitle='Intensity (Photons) Ratio'
endif else begin
   ytitle='Intensity (Ergs) Ratio'
endelse


plot,temperature,ratio,xr=[temperature(0),temperature(ntemp-1)], $
  title=description, $
  xtitle=xtitle,ytitle=ytitle,/xlog,/ylog, xstyle=1


;;--------------------------
;;Create file containing ratio numbers
;;
if keyword_set(outfile) then begin
   openw,luo,outfile,/get_lun
;
   printf,luo, outfile
   printf,luo, ytitle+' calculated with CHIANTI '+version
   printf,luo,desc_ascii

;   if keyword_set(pressure) then begin
;      printf,luo,'Constant Pressure = '+string(pressure,'(e10.2)')
;   endif else if density NE default_density then begin
;      printf,luo,'Constant Density = '+string(density,'(e10.2)')
;   endif else begin
   printf,luo,'Constant Density = '+string(density,'(e10.2)')+' (cm-3)'
;   endelse

   printf,luo, ' '
   printf,luo, 'Temperature (K) ratio '
   printf,luo, ' '

   for i=0l,ntemp-1 do begin
      printf,luo, temperature(i),ratio(i)
   endfor
;
   close,luo
;
ENDIF 


;;----------------------
;;Plot ratio at half the density
;;
newdens=0.5*density
em=emiss_calc(iz,ion,temp=alog10(temperature),dens=alog10(newdens), $
              pressure=pressure,noprot=noprot,radtemp=radtemp, $
              rphot=rphot,/quiet,no_de=photons, ioneq_file=ioneq_name, $
              abund_file=abund_file)

IF n_n GT 1 THEN new_num=total(reform(em[i_n].em),2) $
ELSE new_num=reform(em[i_n].em)
IF n_d GT 1 THEN new_den=total(reform(em[i_d].em),2) $
ELSE new_den=reform(em[i_d].em)

ratiolo=new_num/new_den
oplot,temperature,ratiolo,line=2


;;----------------------
;;Plot ratio at twice the density
;;
newdens=2.0*density
em=emiss_calc(iz,ion,temp=alog10(temperature),dens=alog10(newdens), $
              pressure=pressure,noprot=noprot,radtemp=radtemp, $
              rphot=rphot,/quiet,no_de=photons, ioneq_file=ioneq_name, $
              abund_file=abund_file)

IF n_n GT 1 THEN new_num=total(reform(em[i_n].em),2) $
ELSE new_num=reform(em[i_n].em)
IF n_d GT 1 THEN new_den=total(reform(em[i_d].em),2) $
ELSE new_den=reform(em[i_d].em)

ratiohi=new_num/new_den
oplot,temperature,ratiohi,line=2


print,' '
print,' The dashed lines overplotted are the ratios calculated at '
print,' twice and half the electron density'
print, ' '


;;--------------------------
;; Create ps file containing plot
;;
IF keyword_set(psfile) THEN BEGIN

   ps, psfile, /LANDSCAPE

   plot,temperature,ratio,xr=[temperature(0),temperature(ntemp-1)], $
     title=description, $
     xthick=thick,ythick=thick,thick=thick, $
     xtitle=xtitle,ytitle=ytitle , xstyle=1,/xlog, /ylog

   yplot=0.5*(!y.crange[1]+!y.crange[0])

;   if keyword_set(pressure) then begin
;      xyouts,temperature(1),yplot,'Constant Pressure = '+ $
;        string(pressure,'(e10.2)')+' (cm!e-3!n !eo!nK)',$
;        charsiz=charsiz
;   endif else if density NE default_density then BEGIN

   xyouts,temperature(1),yplot,'Constant Density = '+ $
     string(density,'(e10.2)')+' (cm!e-3!n)',$
     charsiz=charsiz
;   endif else begin
;      xyouts,temperature(1),yplot,'Constant Density = '+ $
;        string(default_density,'(e10.2)')+' (cm!e-3!n)',$
;        charsiz=charsiz
;   endelse
;
   psclose

endif


END

