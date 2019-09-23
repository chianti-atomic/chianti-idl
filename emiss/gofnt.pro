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
;	GOFNT
;
; PURPOSE:
;
;	calculate G(n,T) function (line intensity per unit emission measure)
;       
; PROCEDURE:
;
;	Must specify line to form numerator and denominator
;       Multiple lines can be selected and summed
;       This can now be done interactively or not.
;
; CALLING SEQUENCE: 
;
;	GOFNT,Ion,Wmin,Wmax,Temperature,G,Desc,density=
;
;
; INPUTS:
;
;       Ion:   the CHIANTI style name of the ion, i.e., 'c_5' for C V
;
;       Wmin:  minimum of wavelength wavelength range of interest in Angstroms
;
;       Wmax:  maximum of wavelength wavelength range of interest
;
; OPTIONAL INPUTS:
;
;       Many - see the keywords below.
;
; OUTPUTS:
;
;	Temperature:  an array of temperatures 
;
;       G:  Intensity  per unit emission measure N_e*N_H*dh [cm^-5].
;           The resulting units are therefore erg cm^+3 s^-1 sr-1
;
;           C(T)= 1/(4*!pi)* A_ji*(N_j(X^+m)/N(X^+m))*(N(X^+m)/N(X))*(N(X)/N(H))/N_e
;
;            unless  /NOABUND is set, in which case 
;           C(T)= 1/(4*!pi)* A_ji*(N_j(X^+m)/N(X^+m))*(N(X^+m)/N(X))/N_e
;
;           G(T)=(hc/lambda_ij)*C(T) 
;           G(T)= C(T)   if /PHOTONS is set
;
;       Desc:  a short string description of the selected line
;
;
; OPTIONAL OUTPUTS:
;
;	Postscript file withthe plot of G(T).
;
;       Ascii file with the values of G(T). 
;        
;       VALUE      The array of G(T) values corresponding to logt0.
;
;
; KEYWORD PARAMETERS:
;
;
;	PRESSURE:  specifies the pressure in units of NeT (cm^-3 K).  G is then
;                  calculated at that constant pressure 
;		
;	DENSITY:  specifies the electron density in units of cm^-3.  G is then 
;                 calculated at that value of the electron density.  If neither the 
;                 density or pressure keywords are specified, a constant
;                  density of 1.e+10 cm^-3 is assumed
;
;       PHOTONS:  sets output in photons/s
;
;       RPHOT   Distance from the centre of the star in stellar radius units.
;               I.e., RPHOT=1 corresponds to the star's surface. (Default is
;               infinity, i.e., no photoexcitation.)
;
;       RADTEMP The blackbody radiation field temperature (default 6000 K).
;
;       OUTFILE:  the (optional) name of the output ascii file where a 
;                   listing of the line ratio intensity as a function of 
;                   temperature is saved.
;
;       PSFILE:  the (optional) name of the output postscript file 
;                  where a plot of the chosen G(T) is saved.
;
;       NOABUND: If set, the G(T)'s are not multiplied by the abundance 
;                factor.
;
;       NOPROT   If set, then proton rates are not included.
;
;
;	ABUND_NAME:  Name of the abundance file to use.  If not passed, then
;		     the user is prompted for it.
;
;	IONEQ_NAME:  Name of the ionization equilization name to use.  If not
;		     passed, then the user is prompted for it.
;
;       ALL          If set, all lines are calculated, including
;                    the 'unobserved' ones. 
;
;       LOWER_LEVELS
;       UPPER_LEVELS
;                    Arrays with the indices of the lower and upper levels
;                    pertaining to the transitions you want to get. 
;                    If more than one couple is given, the G(T) of the 
;                    lines are summed.
;                    Obviously, the given indices must correspond to transitions
;                    that are present in the database.
;
;       ARCSECS  
;                 If set, units are photons (ergs) cm^+3 s^-1 arcsecs^-2
;
;       VERBOSE
;
;       LOGT0       An array of log T values for which the G(T) are wanted.
;       VALUE       The array of G(T) values corresponding to logt0.
; 
;                   If logt0 is defined, and within the limits of the 
;                   temperatures for which G(T) NE 0, the array VALUE
;                   is returned with a simple spline interpolation.
;
; CALLS:
;
;       CH_SYNTHETIC, CH_XMENU_SEL
;
; COMMON BLOCKS: None
;
; RESTRICTIONS:
;
; SIDE EFFECTS:
;
; EXAMPLE:
;
;      IDL> gofnt,'o_5',1000.,1500.,temp,goft,desc,density=1.e+16
;
;
; CATEGORY:
;
;	spectral diagnostics
;
;
; MODIFICATION HISTORY:
;
; 	Written by:	Ken Dere
;	October 4, 1996:     Version 1
;       14-Jul-2000     Peter Young, now calls pop_solver
;
;       26-Oct-2000 GDZ, added keyword NOABUND to not multiply for the abundence
;       factor. Corrected header for a wrong description.
;
;	Version 4, 21-Dec-2000, William Thompson, GSFC
;		Modified for better cross-platform graphics capability
;
;  
;       Version 5, 8-Nov-01, Giulio Del Zanna (GDZ). 
;
;       Rewritten as a wrapper routine using the new procedures.
;           Corrected a few inconsistencies in the plots.
;
;       Version 6, 18-Nov-01, Peter Young
;           Added /noprot, rphot and radtemp keywords.
;
;       Version 7, 11-Dec-01, Peter Young
;           Changed call to ch_strpad to strpad.
;
;       V. 8, GDZ, 28-Apr-02 Added abund_name,ioneq_name keywords.
;       v. 9  21-May-2002, GDZ
;             generalized directory concatenation to work for
;             Unix, Windows  and VMS. 
;
;       V.10, 15-Aug-02, GDZ 
;              Major revision:
;              -Removed the call to ch_xselect_s, that did not work for long lists.
;              -Added a '*' in the line lists, to identify 'unobserved' lines.
;              -Replaced the commands to create PS file, to make it
;               cross-platform compatible.
;              -Added a large number of cosmetics, mainly lables to the axes and
;               titles.
;              -Added keyword ALL. If set, all lines are calculated, including
;              the 'unobserved' ones. 
;              -Added the CHIANTI version number in the outputs.
;
;       V. 11, 19-Sep-02, GDZ
;              Clarified output units.
;
;       V.12, 25-Jun-03, GDZ, 
;              Added many new keywords. Now is possible to use the routine 
;              with background jobs, in not-interactive mode.
;              Rounded the wavelengths. 
;
;       V.13, 24-Sept-2003, GDZ 
;              Corrected a bug when logt0 is not defined.
;
;       V.14, 3-Nov-03  GDZ
;             Modified format e8.2 to e9.2 for Windows compatibility.
;
;       V.15, 22-Aug-2008, Peter Young
;             Corrected code where temperature range is reduced based
;             on zero G(T) values.
;	V.16, 22-Jun-2011 Corrected bug that arose when upper_levels or 
;             lower_levels is an array with one element, Terry Kucera
;
;       V.17, 21-Feb-2014, Peter Young
;             Modified behavior so that if logt0 is specified and
;             some of the logt0 values are outside of the range at
;             which G(T) is defined, then Value is set to zero for
;             these temperatures. (Previously the routine just exited
;             if any temperatures were outside the range.)
;
; VERSION     :   17, 21-Feb-2014
;
;-
pro gofnt, ions, wmin, wmax, temperature, gof, description,$
           pressure=pressure, density=density, $
           psfile=psfile, outfile=outfile, photons=photons,  $
           noabund=noabund, noprot=noprot, rphot=rphot, radtemp=radtemp, $
           abund_name=abund_name,ioneq_name=ioneq_name, all=all, $
           lower_levels=lower_levels, upper_levels=upper_levels, $
           arcsecs=arcsecs, VERBOSE=VERBOSE, logt0=logt0, value=value 

on_error, 2

;
if n_params(0) lt 5 then begin
   print,' '
   print,' IDL> gofnt, ion, wmin, wmax, temperature, g, desc'
   print,'        [pressure= , density= , psfile= , outfile= , /photons, '
   print,'         /noabund, abund_name= , ioneq_name= , rphot= , radtemp= ,/noprot, /all, '
   print,'       lower_levels= ,upper_levels= ,/arcsecs, /VERBOSE, logt0= ,value=   ]'
   print,' '
   print,' i.e.>  gofnt, "c_5",40.,50.,temp,g,desc    for C V'
   print,' '
   return
ENDIF

angstrom = ' ('+string(197B)+')'

IF  (keyword_set(pressure) ne 1) and (keyword_set(density) ne 1) then BEGIN
   density=1.e+10
   print,' using constant default density = ',density
ENDIF

IF n_elements(lower_levels) GT 0 AND  n_elements(upper_levels) GT 0 THEN BEGIN 
   IF  n_elements(lower_levels) EQ n_elements(upper_levels) THEN $
     automatic = 1  ELSE BEGIN 
      automatic = 0 
      print,  ' A problem - lower and upper levels do not have the same number of elements'
   END 
ENDIF ELSE automatic = 0 


ch_synthetic, wmin, wmax, output=TRANSITIONS,  err_msg=err_msg, $
  ioneq_name=ioneq_name, $
  pressure=pressure, density=density, all=all,sngl_ion=ions, $
  photons=photons, verbose=verbose, /goft, noprot=noprot, $
  rphot=rphot, radtemp=radtemp


IF err_msg NE '' THEN BEGIN 
   print, err_msg
   return 
ENDIF   

IF  NOT keyword_set(noabund) THEN BEGIN 

   IF n_elements(abund_name) EQ 0 THEN BEGIN 

      dir=concat_dir(!xuvtop,'abundance')

      abund_name=ch_get_file(path=dir,filter='*.abund',title='Select Abundance File')

      ff = findfile(abund_name)
      IF   ff(0) NE '' THEN $
        read_abund,abund_name,abund,abund_ref ELSE BEGIN 
         print, 'Error,  no abundance  file found !'
         return
      END 
   ENDIF ELSE BEGIN 
      ff = findfile(abund_name)
      IF   ff(0) NE '' THEN $
        read_abund,abund_name,abund,abund_ref ELSE BEGIN 
         print, 'Error,  no abundance  file found !'
         return
      END 
   END 

   read_abund,abund_name,abund,abund_ref

END  


nlist = n_elements(TRANSITIONS.lines)


; remove the temperatures where the G(T)'s are zero:

good_t = where(TRANSITIONS.lines[0].goft GT 0, nn)

IF nn LE 1 THEN BEGIN 
   print, 'problems, G(T)=0 !?'
   return
END 

;;
;; PRY, 22-Aug-2008
;; I've added nn=mm which changes the value of nn if mm is smaller.
;;
FOR  i=0,nlist-1 DO  BEGIN
   index = where(TRANSITIONS.lines[i].goft GT 0, mm)
   IF mm LT nn THEN BEGIN
     good_t = index
     nn=mm
   ENDIF
ENDFOR 


temperature = 10.^TRANSITIONS.ioneq_logt[good_t]
n_ioneq_t = n_elements(temperature)

list_wvl = TRANSITIONS.lines[*].wvl

; SORT IN WAVELENGTH 

isort = sort(TRANSITIONS.lines[*].wvl)

list_wvl = TRANSITIONS.lines[isort].wvl
list_descr = TRANSITIONS.lines[isort].ident
list_flag = TRANSITIONS.lines[isort].flag

list_lvl1 = TRANSITIONS.lines[isort].lvl1
list_lvl2 = TRANSITIONS.lines[isort].lvl2


list_int=dblarr(nlist,n_ioneq_t)

FOR  i=0,nlist-1 DO  BEGIN

   ilist = isort[i]

   IF  keyword_set(noabund) THEN $  
     list_int(i,*) = TRANSITIONS.lines[ilist].goft[good_t] ELSE $
     list_int(i,*) = TRANSITIONS.lines[ilist].goft[good_t] *abund[transitions.lines[ilist].iz-1]
ENDFOR 

IF trim(TRANSITIONS.version) EQ '' THEN version = '' ELSE $
  version = 'V. '+TRANSITIONS.version+' '


ion2spectroscopic,ions,species

gof=fltarr(n_ioneq_t)

;description='CHIANTI '+version+species+' '
description=species+' '


IF  keyword_set(photons) THEN  units = 'photons' ELSE units = 'ergs'

IF keyword_set(arcsecs) THEN BEGIN 

   list_int = list_int*(!PI/(180.*60.^2))^2
   units = units+' cm^+3 s^-1 arcsecs^-2 '

ENDIF ELSE  units = units+' cm^+3 s^-1 sr^-1 '

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

IF automatic THEN BEGIN 

   nnum = n_elements(lower_levels)

   IF  nnum EQ  1 THEN  BEGIN 
      index = where(list_lvl1 EQ  lower_levels[0]  AND $
                    list_lvl2 EQ upper_levels[0], nn)
      IF nn NE 1 THEN message, '% GOFNT: Error in the level numbering - EXIT '
      description=description+trim(round(list_wvl(index[0])*100.)/100.)+angstrom
      gof=reform(list_int(index[0],*))

   ENDIF   ELSE BEGIN 

      FOR  j=0,nnum-1 DO  BEGIN 
         index = where(list_lvl1 EQ  lower_levels[j]  AND $
                       list_lvl2 EQ upper_levels[j], nn)
         IF nn NE 1 THEN message, '% GOFNT: Error in the level numbering - EXIT '

         gof=gof+list_int(index[0],*)
         description=description+trim(round(list_wvl(index[0])*100.)/100.)
         IF j NE nnum-1 THEN description=description+'+'
      ENDFOR 

      description=description+angstrom

   END 

ENDIF  ELSE BEGIN 

;
   charsiz=1.5
   xtitle='Temperature (K)'
   title='CHIANTI '+version+species+' '


   ytitle='Relative Intensity (Population of upper level)*A'
;
   intmax=max(list_int)
   intmin=min(list_int)         ;intmax/1000.
;

;device,window_state=ws
;if(ws(0) ne 1) then  
   window,0,ysize=425,xsize=525 

   plot_oo,fltarr(2),xr=[temperature(0),temperature(n_ioneq_t-1)],$
     yr=[intmin,intmax], $
     xtitle=xtitle,ytitle=ytitle,title=title
;
;
   for ilist=0,nlist-1 do begin
      oplot,temperature,list_int(ilist,*)
   endfor
;
   rand=randomu(seed,n_ioneq_t)
   isrt=sort(rand)              ; to try to help with problem of labeling lines
;

   FOR  ilist=0,nlist-1 DO  BEGIN 
      itemp=isrt(ilist<(n_ioneq_t-1))
      xyouts,temperature(itemp),list_int(ilist,itemp),string(ilist,'(i4)'),size=charsiz,align=0.5,noclip=0
      xyouts,temperature(1),list_int(ilist,1),string(ilist,'(i4)'),charsiz=charsiz ,align=0.5,noclip=0
      xyouts,temperature(n_ioneq_t-2),list_int(ilist,n_ioneq_t-2),string(ilist,'(i4)'),charsiz=charsiz,align=0.5,noclip=0
   ENDFOR

;
;  get a reference line - the brightest one.
;
   imax=where(list_int eq max(list_int))
   ref=imax-fix(imax/nlist)*nlist
   oplot,temperature,list_int(ref,*),thick=2
;
   refstr=strtrim(string(ref,'(i4)'),2)
   wvlstr=strtrim(string(list_wvl(ref),'(f10.4)'),2)

   ytitle='Intensity relative to line '+refstr+' at '+wvlstr+angstrom


   window,1,ysize=425,xsize=525 

   y_min = min(list_int/list_int(ref,*))

;  plot ratio to reference line 
;
   plot_oo,fltarr(2),xr=[temperature(0),temperature(n_ioneq_t-1)],$
     yr=[y_min ,1.3], $
     xtitle=xtitle,ytitle=ytitle,title=title

   for ilist=0,nlist-1 do begin
      oplot,temperature,list_int(ilist,*)/list_int(ref,*)
   endfor
;
;
   for ilist=0,nlist-1 do begin
      itemp=isrt(ilist<(n_ioneq_t-1))
      xyouts,temperature(itemp),list_int(ilist,itemp)/list_int(ref,itemp),string(ilist,'(i4)'),charsiz=charsiz,align=0.5,noclip=0
      xyouts,temperature(1),list_int(ilist,1)/list_int(ref,1),string(ilist,'(i4)'),charsiz=charsiz,align=0.5,noclip=0
      xyouts,temperature(n_ioneq_t-2),list_int(ilist,n_ioneq_t-2)/list_int(ref,n_ioneq_t-2),string(ilist,'(i4)'),charsiz=charsiz,align=0.5,noclip=0
   endfor
;

   print,' ------------------------------------------------------'

;
   options=strarr(nlist)

   FOR  i=0,nlist-1 DO  BEGIN

      IF list_flag[i] EQ -1 THEN text_add = ' *' ELSE text_add ='  '

      options[i]=string(i, format='(i4)')+ $
        string(list_wvl(i),format='(f10.4)')+text_add+angstrom+$
        '  '+string(max(list_int(i,*))/max(list_int),format='(e10.2)')+'  '+ $
        list_descr(i)

      print, options(i)

   ENDFOR 

   index = ch_xmenu_sel(options, tit=' Select lines ', $
                        text=['If more than one line is selected, ',$
                              ' the emissivities of the lines will be summed.', $
                              'A  * indicates that the wavelength is theoretical' ])

   IF index[0] EQ -1 THEN BEGIN 
      print, '% GOFNT: Routine Aborted'
      return
   END 

   nnum = n_elements(index)

   IF  nnum EQ  1 THEN  BEGIN 
      description=description+trim(round(list_wvl(index[0])*100.)/100.)+angstrom
      gof=reform(list_int(index[0],*))
   ENDIF  ELSE BEGIN 
      FOR  j=0,nnum-1 DO  BEGIN 
         gof=gof+list_int(index[j],*)
         description=description+ trim(round(list_wvl(index[j])*100.)/100.)
         IF j NE nnum-1 THEN description=description+'+'
      ENDFOR 

      description=description+angstrom

   END 

;endif else if nnum eq 2 then begin
;   description=description+strtrim(string(list_wvl(i(0)),'(f10.3)'),2)+'+'
;   description=description+strtrim(string(list_wvl(i(1)),'(f10.3)'),2)
;endif else begin
;   for j=0,nnum-1 do begin
;      description=description+strtrim(string(list_wvl(i(j)),'(f10.3)'),2)+'+'
;   endfor
;   description=description+strtrim(string(list_wvl(i(nnum-1)),'(f10.3)'),2)
;endelse
;


   print,'-------------------------------------------------------'

   window,2,ysize=425,xsize=525 

   ytitle='Intensity per Emission Measure'
   plot_oo,temperature,gof,xr=[temperature(0),temperature(n_ioneq_t-1)],$
     yr=[min(gof),max(gof)], $
     title=description,xtitle=xtitle,ytitle=ytitle

;
   if keyword_set(pressure) then begin
      xyouts,temperature(1),max(gof)/2.,'Constant Pressure = '+string(pressure,'(e10.2)'),charsiz=charsiz
   endif else if keyword_set(density) then begin
      xyouts,temperature(1),max(gof)/2.,'Constant Density = '+string(density,'(e10.2)'),charsiz=charsiz
   endif else begin
      xyouts,temperature(1),max(gof)/2.,'Constant Density = '+string(density,'(e10.2)'),charsiz=charsiz
   endelse
;

ENDELSE                         ;automatic 


logt = alog10(temperature)      ;TRANSITIONS.ioneq_logt

; check that we have all positive numbers ! 

good = where(gof  GT   0., ngood)

IF ngood NE n_elements(logt) THEN BEGIN 
   print, ' removing zero G(T)'
   logt = logt(good)
   temperature = temperature[good]
   gof = gof(good)
END 

;
; 21-Feb-2014, PRY
; I've modified the section below so that the routine no longer
; exits if logt0 is outside of the range of logt. Instead, G(T) is
; just set to zero for these points.
;
; I've also modified the spline interpolation to use spl_init
; and spl_interp as I think these are more reliable, and I also do the
; interpolation on log(G) rather than just G.
;
IF n_elements(logt0) EQ 1 THEN BEGIN 
   IF  logt0 GT  max(logt) OR logt0 LE min(logt) THEN BEGIN 
      ;; print, 'Error, Log T value outside ranges: '+$
      ;;   trim(min(logt))+' - '+ trim(max(logt))
      ;; return
     value=0d0
   ENDIF ELSE BEGIN 
     x=logt
     y=alog10(gof)
     y2=spl_init(x,y)
     yi=spl_interp(x,y,y2,logt0)
     value=10.^yi
   ENDELSE 
ENDIF ELSE IF  n_elements(logt0) GT  1 THEN BEGIN 

  nt=n_elements(logt0)
  value=dblarr(nt)

  k=where(logt0 GE min(logt) AND logt0 LE max(logt),nk)
  IF nk NE 0 THEN BEGIN
    x=logt
    y=alog10(gof)
    y2=spl_init(x,y)
    yi=spl_interp(x,y,y2,logt0[k])
    value[k]=10.^yi
  ENDIF 

   ;; IF max(logt0) GT  max(logt) OR min(logt0) LE min(logt) THEN BEGIN 
   ;;    print, 'Error, Log T values outside ranges : '+$
   ;;      trim(min(logt))+' - '+ trim(max(logt))
   ;;    return
   ;; END 

   ;; value = fltarr( n_elements(logt0))
   ;; FOR i=0,  n_elements(logt0)-1 DO $
   ;;   value[i] = spline(logt, gof, logt0[i])

END


if keyword_set(psfile) then begin

   ps, psfile, /LANDSCAPE

   plot_oo,temperature,gof,xr=[temperature(0),temperature(n_ioneq_t-1)], $
     yr=[min(gof),max(gof)],title=description, $
     xthick=thick,ythick=thick,thick=thick,xtitle=xtitle, $
     ytitle=ytitle

   psclose

endif

;
; determine temperature of maximum emission
;
atemp=alog10(temperature)
agof=alog10(gof)
dt=0.01
nt1=fix(min(atemp)/dt)+1
nt2=fix(max(atemp)/dt)
itemp=nt1*dt+dt*findgen(nt2-nt1+1)
igof=spline(atemp,agof,itemp)
;  plot,itemp,igof
mt=where(igof eq max(igof))

IF keyword_set(VERBOSE) THEN BEGIN 
   print,'  '+description
   print,' Tmax = '+trim(string(10.^itemp(mt),format='(e9.2)')) +' (K) - '+$
     trim(itemp(mt))+' (log)'
ENDIF 

;
if keyword_set(outfile) then begin
   openw,luo,outfile,/get_lun
;
   printf,luo, outfile
   printf,luo,description+'  intensity per emission measure'
   printf, luo, ' Tmax = '+trim(string(10.^itemp(mt),format='(e9.2)')) +' (K) - '+$
     trim(itemp(mt))+' (log)'
   printf,luo,'  calculated with CHIANTI '+version+' with the gofnt function'

   printf,luo,units 

   IF  NOT keyword_set(noabund) THEN $
     printf,luo,' Abundance file name:  ',abund_name

   printf,luo,' Ionization equilibrium file name:  ',ioneq_name
   if keyword_set(pressure) then begin
      printf,luo,'Constant Pressure = '+string(pressure,'(e10.2)')
   endif else if keyword_set(density) then begin
      printf,luo,'Constant Density = '+string(density,'(e10.2)')
   endif else begin
      printf,luo,'Constant Density = '+string(density,'(e10.2)')
   endelse
   printf,luo, ' '
;
   for i=0,n_ioneq_t-1 do begin
      printf,luo,temperature(i),gof(i),format='(2e12.3)'
   endfor
;
   free_lun,luo
;
endif                           ; outfile
;,format='(a8,f5.2,3x,e8.2)'
;
END 
