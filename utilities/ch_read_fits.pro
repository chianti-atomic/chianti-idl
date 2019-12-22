;+
; PROJECT:  CHIANTI
;
;       http://wwwsolar.nrl.navy.mil/chianti.html
;
;       CHIANTI is an atomic database package for the calculation of
;       continuum and emission line spectra from astrophysical plasmas. It is a 
;       collaborative project involving the Naval Research Laboratory
;       (Washington DC, USA), the Arcetri Observatory (Firenze, Italy), and the
;       Cambridge University (United Kingdom).
;
;
; NAME:
;       CH_READ_FITS
;
; PURPOSE:
;
;       Read  standard CHIANTI FITS binary table data containing the output from
;       CH_SYNTHETIC and output a TRANSITIONS structure.
;
; CALLING SEQUENCE:
;
;       CH_READ_FITS,  Filename, TRANSITIONS
;
; INPUTS:
;
;       Filename = String containing the name of the CHIANTI FITS file written
;       by CH_WRITE_FITS.
;
; OUTPUTS:
;
;       TRANSITIONS = Structure to be written.
;
; OPTIONAL INPUTS: none
;
; KEYWORDS: none
;
;
; NOTES:
; 
;
; CALLS:
;        MRDFITS, ADD_TAG
;	
; COMMON BLOCKS: none.
;
; RESTRICTIONS:
;
;       (3)     The input FITS file must have been written by  CH_WRITE_FITS
;
; PREV. HIST. :
;
;
; EXAMPLE:
;
;         ch_read_fits,  'file.fits', transitions 
;
; WRITTEN     : 
;
;       Ver.1, 8-Apr-02 Giulio Del Zanna (GDZ)
;       V.2 GDZ 31 May 2002 added more checks.
;
; MODIfICATION HISTORY:
;
; VERSION     : 2, 31 May 2002 
;
;-

PRO  ch_read_fits,  file,  transitions, err_msg=err_msg

if n_params() NE 2  THEN  BEGIN 
   print, 'CH_READ_FITS: Usage:'
   print, '    CH_READ_FITS, struct_name, file'
   return 
endif

;on_ioerror, open_error

err_msg = ''

found= findfile(file, count=count)

;check if it exists:

IF count EQ 1 THEN BEGIN 

;check that is a fits file:

test = valid_fits(file)
IF NOT  test THEN BEGIN 
err_msg = file+ ' is not a FITS file !!!!'
return
END 

   index1=mrdfits(file,1)
   index2=mrdfits(file,2)
   transitions=add_tag(index2,  index1, 'lines') 

;check that it is a 'TRANSITIONS' structure at least by having the minimal tags:

test = ch_check_str (transitions, error=error)

IF NOT  test THEN BEGIN 
err_msg =  file+ 'is a wrong FITS file !!!!'
return
END


ENDIF ELSE BEGIN 
err_msg = 'No FITS file found!!!'
return
END 

END 
