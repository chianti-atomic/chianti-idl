;+
; PROJECT
;
;       CHIANTI   http://wwwsolar.nrl.navy.mil/chianti.html
;
;       CHIANTI is an atomic database package for the calculation of
;       continuum and emission line spectra from astrophysical plasmas. It is a 
;       collaborative project involving the Naval Research Laboratory
;       (Washington DC, USA), the Arcetri Observatory (Firenze, Italy), and the
;       Cambridge University (United Kingdom).
;
; NAME
;        temp_plotter
;
; PURPOSE:
;        A widget-based routine to allow the analysis of temperature sensitive
;        ratios. **** See RATIO_PLOTTER for details. *****
;
; CALLING SEQUENCE:
;
;       IDL>  temp_plotter,  name,$ 
;                  EM, PATH=PATH, NOPROT=NOPROT, $
;                  IONEQ_FILE=IONEQ_FILE, ABUND_FILE=ABUND_FILE
;
;
; INPUTS:
;        The ion name (e.g. 'si_3' for Si III)
;
; OPTIONAL INPUTS : none
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
;	EM:	Save the displayed emissivities to structure EM.
;               **** See RATIO_PLOTTER for details. *****
;
; KEYWORDS:
;
;
;	PATH:	Data in the CHIANTI format that is not in the CHIANTI 
;		database can be read by specifying the directory in which 
;		it lies through PATH.
;
;       ABUND_FILE  The name of a CHIANTI abundance file. This is used for 
;               calculating the proton to electron ratio. Default is 
;               !abund_file.
;
;       IONEQ_FILE  The name of a CHIANTI ion balance file. This is used for 
;               calculating the proton to electron ratio and evaluating 
;               the T_max of the ion. Default is !ioneq_file.
;
;       NOPROT  If set, then the default setting will be NOT to use 
;               proton rates. This can be changed within the routine.
;
;       LOOKUP: If set, then routine will attempt to use population
;               lookup tables to compute the ratios.
;
;
; CALLS:   CONVERTNAME RATIO_PLOTTER
;      
;
; COMMON BLOCKS: none
;
;
; RESTRICTIONS:
;
; SIDE EFFECTS:
;
; CATEGORY:
;	spectral synthesis.
;	
; EXAMPLE:
;             IDL> temp_plotter, 'fe_13'
;
; WRITTEN     : 
;              Ver.1, 18-Apr-2002, Giulio Del Zanna (GDZ) written as a wrapper
;              routine to call RATIO_PLOTTER (written by P. Young).
;
; MODIFIED:   V.2,  2-Aug-2005, GDZ
;              Now the routine handles the dielectronic case
;
;             V.3, 18-Dec-2019, Peter Young
;
;-


PRO temp_plotter,  name,$ 
                  EM, PATH=PATH, NOPROT=NOPROT, $
                   IONEQ_FILE=IONEQ_FILE, ABUND_FILE=ABUND_FILE, $
                   lookup=lookup

;
if n_params(0) lt 1 then begin
   print,' '
   print,' IDL>  temp_plotter,  name,$ '
   print,'       EM, PATH=PATH, NOPROT=NOPROT, $ '
   print,'       IONEQ_FILE=IONEQ_FILE, ABUND_FILE=ABUND_FILE, /LOOKUP'
   print,'       '
   print,' i.e.> temp_plotter, "o_5" '
   print,' '
   return
ENDIF

convertname,name,ION_Z, ION_SP,dielectronic=diel

RATIO_PLOTTER, ION_Z, ION_SP, diel=diel, /TEMP, $
     EM, PATH=PATH, NOPROT=NOPROT, $
               IONEQ_FILE=IONEQ_FILE, ABUND_FILE=ABUND_FILE, $
               lookup=lookup

END

