

FUNCTION ch_ip, ion_name

;+
; NAME:
;      CH_IP
;
; PURPOSE:
;      Returns the ionization potential of the specified ion, i.e.,
;      the minimum energy required to remove an electron from the ion.
;
; CATEGORY:
;      CHIANTI; ionization potential.
;
; CALLING SEQUENCE:
;	Result = CH_IP( Ion_Name )
;
; INPUTS:
;      Ion_name:  The name of an ion in CHIANTI format, e.g., 'o_6'
;                 for O VI.
;
; OUTPUTS:
;      The ionization potential in eV. If a problem is found then -1
;      is returned.
;
; CALLS:
;      READ_IP, CONVERTNAME
;
; EXAMPLE:
;      IDL> print,ch_ip('h_1')
;              109678.77
;
; MODIFICATION HISTORY:
;      Ver.1, 14-Jul-2016, Peter Young
;-


IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> ip=ch_ip(ion_name)'
  print,''
  print,'  Ionization potential returned in eV units.'
  return,-1
ENDIF 

ipfile=concat_dir(!xuvtop,'ip/chianti.ip')
read_ip,ipfile,ip,ref

convertname,ion_name,iz,ion

return,ip[iz-1,ion-1]

END
