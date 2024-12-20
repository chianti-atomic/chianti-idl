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
;	LATEX_WVL_DEM
;
; PURPOSE:
;
;	create a latex file of predicted spectral line intensities and
;       wavelengths corresponding to a selected abundance and differential
;       emission measure (DEM)
;
;
; PROCEDURE:
;
;                 Calculations are done assuming either constant density or
;                 constant pressure. See CH_SYNTHETIC for details.
;
;
; CALLING SEQUENCE:
;
;       LATEX_WVL_DEM, Wmin, Wmax, Pressure= , [density= ], $ 
;              [outfile= , mini= , sngl_ion=, /photons, /all, /masterlist]
;
;
; INPUTS:
;
;	Wmin:   lower limit of the wavelength/energy range of interest (Angstroms)
;               if kev keyword set, then wmin is in kev	
;
;	Wmax:   upper limit of the wavelength/energy range of interest (Angstroms)
;               if kev keyword set, then wmax is in kev	
;
;       Pressure:  pressure in emitting region (cm^-3 K), or 
;       Density:   density in emitting region (cm^-3).
;
;
; OPTIONAL INPUTS:
;
;       OUTFILE: the name of the output latex file to be written.
;	
;	MINI:	Minimum intensity for a line to be included in the output.
;
;	SNGL_ION:  specifies  a single ion (e.g. SNGL_ION='Fe_10' to include
;                 only Fe X lines) or an array (e.g. SNGL_ION=['Fe_10','Fe_11']
;                 to include only Fe X and Fe XI lines) of ions to be used
;                 instead of the complete set of ions specified in
;                 !xuvtop/masterlist/masterlist.ions 
;
;       MASTERLIST: string of a specific masterlist file (full path). 
;                   If defined as a keyword (i.e. MASTERLIST=1 or /MASTERLIST)
;                   then a widget allows the user to select a  user-defined
;                   masterlist file. Shortcut for SNGL_ION.   
;
;       RPHOT   Distance from the centre of the star in stellar radius units.
;               I.e., RPHOT=1 corresponds to the star's surface. (Default is
;               infinity, i.e., no photoexcitation.)
;
;       RADTEMP The blackbody radiation field temperature (default 6000 K).
;
;       ABUND_NAME:  The name of a CHIANTI abundance file. If not specified, then
;                    the user will be asked to select a file from a drop-down
;                    list.
; 
;       IONEQ_NAME:  The name of a CHIANTI ioneq file. If not specified, then
;                    the user will be asked to select a file from a drop-down
;                    list.
;
; OUTPUTS:
;       A Latex file containing the table of line intensities. The default name
;       is 'linelist.text' but this can be over-ridden with OUTFILE=.
;
; KEYWORD PARAMETERS:
;
;       PHOTONS:  units will be in photons rather than ergs
;
;       KEV:  wavelengths will be given in kev rather than Angstroms
;
;       ALL:  if set, then all lines are included.  This means that lines for which
;             only an approximate wavelength is known (only theoretical energy
;             levels are known) are included.
;
;       NOPROT   If set, then proton rates are not included.
;
;       LOOKUP:  If set, then lookup tables are used, greatly speeding up the
;                calculation. 
;
;       ADVANCED_MODEL: Set this to zero to switch off the advanced models. 
;
; CALLS:
;       
;       CH_SYNTHETIC, CH_LINE_LIST
;
;
; EXAMPLE:
;
;       > latex_wvl_dem, 400.,800., mini=1, pressure=1.e+15,sngl_ion='o_4'
;
;
; CATEGORY:
;
;	spectral synthesis.
;
; WRITTEN: 
;
;       Version 1, 8-Nov-01, Giulio Del Zanna (GDZ). 
;
;       Compared to the previous LATEX_WVL_DEM, these are the main changes:
;
;       1-Rewritten as a wrapper routine using the new procedures.
;       2-Now the PRESSURE value is a keyword.
;       3-The calculations can be done at constant DENSITY.
;       4-MASTERLIST can now be used both as an input string or as a keyword.
;
;
; MODIFICATION HISTORY:
;
;       Version 2, 18-Nov-01, Peter Young
;           Added /noprot, rphot and radtemp keywords.
;
;       V. 3, 22-May-2002 GDZ.  Removed const_net definitions.
;
;       Ver.4, 08-Jun-2022, Peter Young
;         Added abund_name= optional input.
;
;       Ver.5, 02-Dec-2024, Peter Young
;         Added ioneq_name= optional input and advanced_model keyword.
;             
; VERSION     : Version 5, 02-Dec-2024
;
;-

PRO  latex_wvl_dem,wmin,wmax,pressure=pressure, density=density,$
                   minI=minI,photons=photons,kev=kev,all=all,$
                   sngl_ion=sngl_ion,masterlist=masterlist, $
                   outfile=outfile, noprot=noprot, rphot=rphot, radtemp=radtemp, $
                   dem_name=dem_name, abund_name=abund_name, lookup=lookup, $
                   ioneq_name=ioneq_name, advanced_model=advanced_model

;
if n_params(0) lt 2 then begin
   print,'  IDL>  latex_wvl_dem,wmin,wmax,pressure=, [density= ], $ '
   print,'     [outfile= , mini= , sngl_ion=, /photons, /all, /masterlist, '
   print,'      radtemp= , rphot= , /noprot ]'
   return
endif
;

IF n_elements(outfile) EQ 0 THEN outfile ='linelist.tex'

IF n_elements(ioneq_name) EQ 0 THEN ioneq_name=!ioneq_file

IF  keyword_set(pressure)  AND  keyword_set(density)  THEN  begin
   print,' You  have  to decide if you want the intensities calculated '+$
     'for constant pressure or not ! '
   print,'(e.g.  pressure=1.e16) or '+$
     'for constant electron density (e.g.  density=1.e9) !!!'
   err_msg ='Error '
   return
endif

IF n_elements(density) NE 0 THEN BEGIN 
   print,' using constant Density = ',density
ENDIF ELSE IF  n_elements(pressure) NE 0 THEN BEGIN 
 print,' using constant Pressure =', pressure 
ENDIF ELSE message, 'No density or pressure defined -- EXIT '


ch_synthetic, wmin, wmax, output=TRANSITIONS, err_msg=err_msg, $
              pressure=pressure, density=density,all=all,sngl_ion=sngl_ion, $
              photons=photons,  masterlist=masterlist,  $
              verbose=verbose,$
              logt_isothermal=logt_isothermal,  logem_isothermal=logem_isothermal,$
              noprot=noprot, rphot=rphot, radtemp=ratemp, $
              dem_name=dem_name,lookup=lookup, ioneq_name=ioneq_name, $
              advanced_model=advanced_model


IF err_msg NE '' THEN BEGIN 
   print, err_msg
   return 
ENDIF ELSE BEGIN 

   ch_line_list, TRANSITIONS, OUTFILE, /latex, $
                 wmin=wmin,wmax=wmax,$
                 minI=minI,photons=photons,kev=kev, $
                 all=all,no_sort=no_sort, sngl_ion=sngl_ion, $
                 abundfile=abund_name, max_int=max_int

   print,''
   print,format='("NOTE: maximum intensity displayed in table is: ",e10.2)',max_int

END 

;
END 
;


