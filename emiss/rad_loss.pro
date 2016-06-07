
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
;	RAD_LOSS
;
; PURPOSE:
;
;       Calculates energy loss rate by free-free (ff), radiative 
;       recombination (fb) and by line (bound-bound) radiation.
;
; CATEGORY:
;	
;	synthetic spectra
;
; CALLING SEQUENCE:
;
;       RAD_LOSS,Temperature,Loss_rate
;
;
; INPUTS:
;	None:  The user will select an elemental abundance file and a 
;              ionization equilibrium file through widgets.
;
; OPTIONAL INPUTS:
;	Abund_File:  Specifies the element abundance file to be
;                    used. If not set, then the file is selected
;                    through a widget interface.
;       Ioneq_File:  Specifies the ionization equilibrium file. If not
;                    set, then the default file (!ioneq_file) is used.
;
;       Pressure:  Electron pressure in emitting region (cm^-3 K)
;                  density=pressure/temperature(K)
;
;       Density:   The electron number density (cm^-3), constant for all temperatures.
;                  If neither density or pressure is set, then a default 
;                  constant density of 10^10 cm^-3 is used.
;
;       Psfile:    to create a postscript plot of the radiative loss in 
;                  the file specified by the name assigned to 'Psfile'
;
;       Outfile:   name of an ascii file where the radiative loss rate 
;                  as a function of temperature is output
;                   
;	RADTEMP	Specify background radiation temperature (default: 6000 K)
;
;	RPHOT   Distance from the centre of the star in stellar radius units.
;               I.e., RPHOT=1 corresponds to the star's surface. (Default is
;               infinity, i.e., no photoexcitation.)
;
;  KEYWORD PARAMETERS:
;       NOPROT  If set, then the default setting will be NOT to use 
;               proton rates. This can be changed within the routine.
;
; OUTPUTS:
;
;       Temperature:  array of temperatures (K)
;       Loss_rate:    energy loss rate in erg cm^3 s^-1
;
;
; PROCEDURE:
;
;
;  if keyword pressure is set then calculations performed at constant pressure
;  if keyword density is set then calculations performed at constant density
;  otherwise, density is set to 1.e+10
;  
;  pressure = density * temperature  (cm^-3 K)
;
;	the user will be asked to select an abundance file and a 
;       ionization balance file.
;
; EXAMPLE:
;
;       > rad_loss,t,r
;
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	January 1999:  version 1, adopted from synthetic.pro
;	Version 2, 21-Dec-2000, William Thompson
;		Modified for better cross-platform capability.
;
;       Ver.3, 6-Dec-2001, Peter Young
;               Added /noprot, radtemp and dilute keywords.
;               Removed elvlc, wgfa and upsilon common blocks.
;               Removed calls to read_ip and read_masterlist (not needed).
;
;       Ver.4, 8-Jul-2003, Peter Young
;               Updated routine header (no changes to code).
;
;       V 5, 25-May-2005, GDZ 
;                  corrected routine header.
;
;       Ver.6, 12-Jul-2007, EL
;                  corrected the output units in the y-axis of the display
;
;       Ver.7, 4-Apr-2016, Peter Young
;               Updated header to state that the pressure and density
;               inputs are for electrons.
;
;       Ver.8, 8-Aug-2017, Peter Young
;               Added abund_file and ioneq_file optional inputs.
;
; VERSION     :   v.8, 8-Aug-2017, Peter Young
;-

pro RAD_LOSS,t,loss_rate,pressure=pressure,density=density, $
             psfile=psfile,outfile=outfile,noprot=noprot, $
             radtemp=radtemp,rphot=rphot, abund_file=abund_file, ioneq_file=ioneq_file 
;
;
;  pressure= electron pressure (Ne x T)  cm-3 K
;
;
;
;
common elements,abund,abund_ref,ioneq,ioneq_t,ioneq_ref
;
if n_params() lt 2 then begin
   print,' type> rad_loss,temperature,loss_rate [,pressure= , density= ,psfile= '
   print,'                   radtemp= , rphot= , /noprot ]'
   return
endif
;
;
default_density=1.e+10
if (keyword_set(pressure) ne 1) and (keyword_set(density) ne 1) then begin
      print,' using constant default density = ',default_density
endif
;
;
;  Read elemental abundances, ionization equilibrium and ionization 
;  potentials
;
setup_elements, abund_file=abund_file, ioneq_file=ioneq_file

;
; Call the three individual components of the rad loss. No that
; /no_setup is used.
;
ff_rad_loss,t,ff,/no_setup
fb_rad_loss,t,fb,/no_setup
if keyword_set(pressure) then begin
  bb_rad_loss,t,bb,/no_setup,pressure=pressure,noprot=noprot, $
       radtemp=radtemp,rphot=rphot
endif else bb_rad_loss,t,bb,/no_setup,density=density,noprot=noprot, $
     radtemp=radtemp,rphot=rphot
;
;
loss_rate=ff+fb+bb
;
;
   thick=1
   xtitle='Temperature (K)'
   ytitle='Radiative energy losses (ergs cm!u3!n s!u-1!n)'
   title='CHIANTI:  radiative loss rate'
   plot_oo,t,loss_rate,xr=minmax(t),yr=minmax(loss_rate),title=title, $
     xthick=thick,ythick=thick,thick=thick,xtitle=xtitle,ytitle=ytitle
;
if keyword_set(psfile) then begin
   dname = !d.name
   set_plot,'ps'
   device,filename=psfile
   landscape
   device,/times,/isolatin1,font_size=16
   !p.font=0
   thick=3
   plot_oo,t,loss_rate,xr=minmax(t),yr=minmax(loss_rate),title=title, $
     xthick=thick,ythick=thick,thick=thick,xtitle=xtitle,ytitle=ytitle
   !p.font=-1
   device,/close
   set_plot,dname
endif
;
;
if keyword_set(outfile) then begin
   openw,luo,outfile,/get_lun
;
   printf,luo, outfile
   printf,luo,'  Radiative loss rate calculated by CHIANTI'
   if keyword_set(pressure) then begin
      printf,luo,'Constant Pressure = '+string(pressure,'(e10.2)')
   endif else if keyword_set(density) then begin
      printf,luo,'Constant Density = '+string(density,'(e10.2)')
   endif else begin
      printf,luo,'Default Density used'
   endelse
   printf,luo,' '
   printf,luo,'  Elemental abundances'
   nref=n_elements(abund_ref)
   for i=0,nref-1 do printf,luo,abund_ref(i)
   printf,luo,' '
   printf,luo,'  Ionization equilibrium'
   nref=n_elements(ioneq_ref)
   for i=0,nref-1 do printf,luo,ioneq_ref(i)
   printf,luo, ' '
;
   ntemp=n_elements(t)
   for i=0,ntemp-1 do begin
       printf,luo,t(i),loss_rate(i),format='(2e10.2)'
   endfor
;
   close,luo
endif
;
;
;
end
