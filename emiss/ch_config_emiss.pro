
PRO ch_config_emiss, ion, output, dens=dens, noplot=noplot, joules=joules

;+
; NAME
;
;      CH_CONFIG_EMISS
;
; PROJECT
;
;      CHIANTI
;
; EXPLANATION
;
;      Takes the output from emiss_calc and averages emissivities over
;      configurations.
;
;      The output emissivities are calculated as 
;
;      SUM ( E_ji * n_j * A_ji / N_e )
;
;      where E_ji is the energy of the transition, n_j is the
;      population of the upper emitting level (where SUM(n_j)=1), A_ji
;      is the radiative decay rate for the transition, and N_e is the
;      electron number density. The sum is performed over all
;      transitions emanating from the configuration.
;
; INPUTS
;
;      ION    The name of the ion in CHIANTI format. E.g., 'fe_13' for
;             Fe XIII.
;
; OPTIONAL INPUTS
;
;      DENS   The density at which the emissivities are calculated.
;
; KEYWORDS
;
;      NOPLOT If set, then the emissivities are not plotted to the
;             screen. 
;
;      JOULES If set, then the emissivities will be written out in
;             joules rather than ergs.
;
; OUTPUTS
;
;      OUTPUT A structure array with n elements where n is the
;             number of configurations in the CHIANTI ion model. The
;             tags of the structure are:
;              .n   The index of the configuration
;              .temp The temperatures at which the emissivities are
;                    calculated (units: K).
;              .dens The density at which the emissivities are
;                    calculated (units: cm^-3).
;              .emiss The emissivities (units: erg s^-1 cm^-3).
;
; HISTORY
;
;      Ver. 1, 6-Jan-2009, Peter Young
;-

IF n_params() LT 2 THEN BEGIN
  print,'Use:  IDL> ch_config_emiss, ionname, output [, dens=, /noplot, /joules]'
  return
ENDIF

IF n_elements(dens) EQ 0 THEN dens=1e10

convertname,ion,iz,ion
zion2filename,iz,ion,filename
elvlcname=filename+'.elvlc'

read_elvlc,elvlcname,l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,ref

temp=findgen(41)/10.+4.
nt=n_elements(temp)
emiss=emiss_calc(iz,ion,dens=alog10(dens),temp=temp,/quiet)

;
; conf_em is an array of same size as emiss giving the configuration
; index of each transition
;
conf_em=conf[emiss.level2-1]

n=max(conf)

str={n: 0, temp: 10.^temp, dens: dens, emiss: dblarr(nt)}
output=replicate(str,n)

IF keyword_set(joules) THEN conv=1e-7 ELSE conv=1.0

;
; Don't include the ground configuration in loop
;
; Note that the emissivities are divided by the density. This is
; consistent with the radiative loss curve.
;
FOR i=1,n DO BEGIN
  output[i-1].n=i
  k=where(conf_em EQ i,nk)
  IF nk NE 0 THEN BEGIN
    IF nk GT 1 THEN em=total(emiss[k].em,2) ELSE em=emiss[k].em
    output[i-1].emiss=em/dens*conv
  ENDIF 
ENDFOR

IF NOT keyword_set(noplot) THEN BEGIN
  plot,/nodata,[min(temp),max(temp)],[0,max(output.emiss)*1.10],/xsty,/ysty
  FOR i=0,n-1 DO BEGIN
    oplot,temp,output[i].emiss
    getmax=max(output[i].emiss,imax)
    xyouts,temp[imax],getmax,trim(output[i].n),charsiz=1.5,align=0.5
  ENDFOR
ENDIF

END
