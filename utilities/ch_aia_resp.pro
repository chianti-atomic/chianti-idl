;+
; NAME:
;
;      CH_AIA_RESP()
;       
; PURPOSE:
;
;      To compute the SDO AIA temperature responses of the six coronal
;      channels for a specific date.
;
; CATEGORY:
;
;      CHIANTI; emissivity.
;
; EXPLANATION:
;
;      This routine calculates  the SDO AIA temperature responses of the six coronal
;      channels. The routine first obtained  the time-dependent AIA
;      effective areas from the AIA SSW, then calculates  CHIANTI
;      isothermal spectra between 10^4 and 10^8 K and between 25 and
;      900 Angstroms, at a  0.1  resolution, either at constant denity
;      or constant pressure.   Note that due to the
;      off-band AIA sensitivity, significant emission in some bands
;      comes from wavelengths above 450 Angstroms, which were not
;      originally included. Also note that the CHIANTI model for these
;      longer wavelengths is not very accurate. 
;      The  isothermal spectra are then multiplied by the effective
;      areas and integrated over wavelength to produce the temperature
;      responses. 
;
; CALLING SEQUENCE:
;
;      AIA_RESP=CH_AIA_RESP(time)
;
; EXAMPLE:
;
;      To reproduce responses close to the AIA ones which were calculated
;      with CHIANTI v.7.1:
;
;      AIA_RESP= ch_aia_resp('20100522', pressure=1e15, $
;      abund_name=!xuvtop+'/abundance/sun_coronal_1992_feldman_ext.abund',$ 
;      ioneq_name=!xuvtop+'/ioneq/chianti.ioneq')
;
; INPUTS:
;
;       TIME 
;       The time of the observations, passed to the AIA routines.
; 
;	PRESSURE
;       The pressure, assumed to be in 	K * cm^-3.
;
;       DENSITY
;       The density,  assumed to be in 	 cm^-3.
;
;
; OPTIONAL INPUTS:
;
;       ABUND_FILE  The name of a CHIANTI abundance file.
;
;       IONEQ_FILE  The name of a CHIANTI ion balance file.
;
;
; OUTPUT:
;       An IDL structure, containing:
;
;       TEMPERATURES
;          the array of temperatures
;       COMMENT
;          a comment on what parameters were used for the calculation.
;       RESP_X
;          the response for the X channel 
;
; CALLS:
;       AIA_GET_RESPONSE (AIA), ISOTHERMAL (CHIANTI), INTERPOL (IDL)
;
; HISTORY:
;	Version 1, Giulio Del Zanna, 2018 Dec 5
;
;
; VERSION     :  1
;
;-


function ch_aia_resp, time,  pressure=pressure, density=density, $ 
           abund_name=abund_name, ioneq_name=ioneq_name, verbose=verbose

if n_elements(density) gt 0 and  n_elements(pressure) gt 0 then begin 
print,' % CH_AIA_RESP: Error, you need to define either a constant density or constant pressure' 
return, -1
end 

if n_elements(density) eq 0 and  n_elements(pressure) eq 0 then begin 
print,' % CH_AIA_RESP: Error, you need to define either a constant density or constant pressure' 
return, -1
end 


ssw_path,/aia

; Solid angle of one AIA pixel:
sterad_aia_pix=8.4d-12

; Get the time-dependent AIA effective areas:
aia_resp = aia_get_response(/dn,/area, evenorm=0, timedepend_date=time)

verstr=ch_get_version()

comment='AIA responses calculated in the 25-900 Angstrom range with CHIANTI version '+$
   verstr+' data '


if n_elements(density) gt 0 then $
comment=[comment, ' at a fixed density= '+ trim(density)] else $
if n_elements(pressure) gt 0 then $
comment=[comment, ' at a fixed pressure= '+ trim(pressure)]

;logt=indgen(81)*0.05+4.0
logt=indgen(26)*0.1+5.0
temp = 10.^logt

isothermal, 25, 900, 0.1, temp, lambda,spectrum,list_wvl,list_ident,$
                pressure=pressure, edensity=density,/photons,/cont , $
                  abund_name=abund_name, ioneq_name=ioneq_name, verbose=verbose


comment=[comment, ' with the elemental abundance file: '+abund_name]
comment=[comment, ' with the ionisation equilibrium file: '+ioneq_name]



eff_94=interpol(aia_resp.a94.ea, aia_resp.a94.wave,  lambda)
eff_131=interpol(aia_resp.a131.ea, aia_resp.a131.wave,  lambda)
eff_171=interpol(aia_resp.a171.ea, aia_resp.a171.wave,  lambda)
eff_193=interpol(aia_resp.a193.ea, aia_resp.a193.wave,  lambda)
eff_211=interpol(aia_resp.a211.ea, aia_resp.a211.wave,  lambda)
eff_335=interpol(aia_resp.a335.ea, aia_resp.a335.wave,  lambda)


sp_conv= spectrum & sp_conv[*,*]=0.
for i=0,n_elements(temp)-1 do $
sp_conv[*, i]= sterad_aia_pix*spectrum[*, i] * eff_94
resp_94=total(sp_conv,1)

sp_conv= spectrum & sp_conv[*,*]=0.
for i=0,n_elements(temp)-1 do $
sp_conv[*, i]= sterad_aia_pix*spectrum[*, i] * eff_131
resp_131=total(sp_conv,1)

sp_conv= spectrum & sp_conv[*,*]=0.
for i=0,n_elements(temp)-1 do $
sp_conv[*, i]= sterad_aia_pix*spectrum[*, i] * eff_171
resp_171=total(sp_conv,1)

sp_conv= spectrum & sp_conv[*,*]=0.
for i=0,n_elements(temp)-1 do $
sp_conv[*, i]= sterad_aia_pix*spectrum[*, i] * eff_193
resp_193=total(sp_conv,1)

sp_conv= spectrum & sp_conv[*,*]=0.
for i=0,n_elements(temp)-1 do $
sp_conv[*, i]= sterad_aia_pix*spectrum[*, i] * eff_211
resp_211=total(sp_conv,1)

sp_conv= spectrum & sp_conv[*,*]=0.
for i=0,n_elements(temp)-1 do $
sp_conv[*, i]= sterad_aia_pix*spectrum[*, i] * eff_335
resp_335=total(sp_conv,1)

resp={temperatures:temp,comment: comment, $
     resp_94:resp_94,resp_131:resp_131,resp_171:resp_171, $
    resp_193:resp_193,resp_211:resp_211, resp_335:resp_335}

return, resp


end
