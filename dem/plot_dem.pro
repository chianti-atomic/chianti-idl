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
;	PLOT_DEM
;
; PURPOSE:
;
;	To plot differential emission measure (DEM) values
;
; CATEGORY:
;
;	Widgets.
;
; CALLING SEQUENCE:
;
;       PLOT_DEM,filename
;
;
; INPUTS:
;
;	filename:  the name of the DEM file to be plotted.  The file must b
;                  in the standard CHIANTI format for DEM files.  If filename 
;                  is a blank string ('') then an interactive window will come 
;                  up to allow the user to select a DEM file from the CHIANTI 
;                  DEM directory or some other directory.	
;
;	
; KEYWORD PARAMETERS:
;
;	PSFILE:	If set, the a postscript plot will be place in the 
;               file 'psfile' specified by the user
;
;
; OUTPUTS:
;
;       None, other than a plot
;
;
;
;
; EXAMPLE:
;
;             > plot_dem,'ademfile.dem'
;         or
;             > plot_dem,''
;
; MODIFICATION HISTORY:
;
; 	Written by:	Ken Dere
;	June 1998:     Version 1.0
;	Version 2, 21-Dec-2000, William Thompson, GSFC
;		Modified for better cross-platform graphics capability
;
;       V.   3, 21-May-2002, Giulio Del Zanna (GDZ) 
;                   generalized directory concatenation to work for
;                   Unix, Windows  and VMS. 
;
; VERSION     :   3, 21-May-2002 
;
;
;-
pro plot_dem,name,psfile=psfile
;
if n_params(0) lt 1 then begin
   print,''
   print,' > plot_dem, filename, [psfile=] '
   print,'      or'
   print,' > plot_dem,''','''
   print,''
   return
endif
;
if strtrim(name,2) eq '' then begin
;
    dir=concat_dir(!xuvtop,'dem')
    name=dialog_pickfile(path=dir,filter='*.dem',title='Select DEM File')
endif
;
read_dem,name,logt,dem,demref
;
plot,logt,dem,xr=[4.,7.5],yr=[20.,24.],title=name
;
if keyword_set(psfile) then begin
;
   dname = !d.name
   set_plot,'ps'
   device,filename=psfile
   device,/landscape
   plot,logt,dem,xr=[4.,7.5],yr=[20.,24.],title=name
   device,/close
   set_plot, dname
endif
;
END


