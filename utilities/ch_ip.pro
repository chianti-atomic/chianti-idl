

FUNCTION ch_ip, ion_name, cm=cm, ryd=ryd, units=units

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
; KEYWORD PARAMETERS:
;      CM:    If set, then the ionization potential is returned in
;             cm^-1 units.
;      RYD:   If set, then the ionization potential is returned in
;             Rydberg units.
;
; OUTPUTS:
;      The ionization potential in eV. If a problem is found then -1
;      is returned.
;
; OPTIONAL OUTPUTS:
;      Units:  A string containing the units of the output number.
;
; CALLS:
;      READ_IP, CONVERTNAME
;
; EXAMPLE:
;      IDL> print,ch_ip('h_1')
;              13.598434
;      IDL> print,ch_ip('h_1',/cm)
;              109678.77
;      IDL> print,ch_ip('h_1',/ryd)
;              0.99946648
;
; MODIFICATION HISTORY:
;      Ver.1, 14-Jul-2016, Peter Young
;      Ver.2, 4-Jun-2018, Peter Young
;         Added /cm and /ryd keywords, and units= optional output.
;-


IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> ip=ch_ip(ion_name)'
  print,''
  print,'  Ionization potential returned in eV units. Use /cm or /ryd for alternative units.'
  return,-1
ENDIF 

ipfile=concat_dir(!xuvtop,'ip/chianti.ip')
read_ip,ipfile,ip,ref

convertname,ion_name,iz,ion


output=ip[iz-1,ion-1]

CASE 1 OF
  keyword_set(ryd): BEGIN
    output=output/109737.32
    units='Rydberg'
  END 
  keyword_set(cm): units='cm^-1'
  ELSE: BEGIN
    output=output/8065.5446
    units='eV'
  END 
ENDCASE 


return,output

END
