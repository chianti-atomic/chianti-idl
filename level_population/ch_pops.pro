

FUNCTION ch_pops, ionname, dens=dens, temp=temp, ldens=ldens, ltemp=ltemp, $
                  _extra=extra

;+
; NAME:
;      CH_POPS()
;
; PURPOSE:
;      Compute level populations for the specified ion. Note that this
;      routine is a wrapper for show_pops, and it accepts most of the
;      same keywords.
;
; CATEGORY:
;      CHIANTI; level populations.
;
; CALLING SEQUENCE:
;      Result = CH_POPS( IonName )
;
; INPUTS:
;      Ionname:  The name of an ion in CHIANTI format, e.g., 'fe_13'
;                for Fe XIII.
;
; OPTIONAL INPUTS:
;      Temp:  A temperature in kelvin. If not specified then the
;             temperature of maximum ionization is used. 
;      Ltemp: The logarithm of temperature (in K).
;      Dens:  An electron number density in cm^-3. If not specified,
;             then 10^10 is assumed.
;      Ldens: The logarithm of electron number density (in cm^-3).
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;      Returns a structure with the tags:
;         DENS            DOUBLE       1.0000000e+10
;         TEMP            DOUBLE           281838.17
;         LEVEL           STRUCT    -> <Anonymous> Array[86]
;         RADTEMP         FLOAT           10000.0
;         RPHOT           FLOAT           0.00000
;         PROTON          STRING    'yes'
;         VERSION         STRING    'CHIANTI 8.0.2'
;         DATE            STRING    'Tue May 17 11:01:38 2016'
;         SUM_MWL         INT              0
;         SUM_MWL_COEFFS  FLOAT          -1.00000
;
; EXAMPLE:
;      IDL> p=ch_pops('o_6')
;      IDL> p=ch_pops('fe_13',ldens=10,ltemp=6.1)
;      IDL> p=ch_pops('mg_5',temp=1e5)
;
; MODIFICATION HISTORY:
;      Ver.1, 17-May-2016, Peter Young
;-


IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> pop = ch_pops ( ion_name [, dens=, ldens=, temp=, '
  print,'                           ltemp=, rphot=, radtemp=, /all, /noprot'
  print,'                           path=, n_levels=, /diel, ioneq_file=, '
  print,'                           abund_file=, sum_mwl_coeffs=, radfunc=, '
  print,'                           level=, /quiet )'
  print,''
  print,"e.g.,  IDL> p=ch_pops('o_6')"
ENDIF 

convertname,ionname,iz,ion

;
; Handle the temperature and density inputs. Note that if both temp
; and ltemp are input, then temp is ignored.
;
IF n_elements(temp) NE 0 THEN temp_in=alog10(temp)
IF n_elements(dens) NE 0 THEN dens_in=alog10(dens)
;
IF n_elements(ltemp) NE 0 THEN temp_in=ltemp
IF n_elements(ldens) NE 0 THEN dens_in=ldens

show_pops,iz,ion,popstr,dens=dens_in, temp=temp_in, _extra=extra

return,popstr

END
