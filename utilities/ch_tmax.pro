FUNCTION ch_tmax, ionname, ioneqname=ioneqname, log=log, interp=interp, $
                  advanced_model=advanced_model, pressure=pressure, $
                  density=density, quiet=quiet

;+
; NAME:
;     CH_TMAX()
;
; PURPOSE:
;     Returns the temperature of maximum ionization (T_max) for the
;     specified ion.
;
; CATEGORY:
;     CHIANTI; ionization fractions.
;
; CALLING SEQUENCE:
;     Result = CH_TMAX( IonName )
;
; INPUTS:
;     Ionname:  The name of the ion in CHIANTI format. E.g., 'fe_13'
;               for Fe XIII.
;
; OPTIONAL INPUTS:
;     Ioneqname: The routine reads the default CHIANTI ion balance
;                file (!ioneq_file). To use a different ion balance
;                file, specify the full pathname with this keyword.
;     Pressure:  Specify the pressure (units: K cm^-3) for which the
;                ion balance should be calculated. Only valid for
;                advanced models.
;     Density:   Specify the density (units: cm^-3) for which the
;                ion balance should be calculated. Only valid for
;                advanced models.
;     Advanced_Model:  The advanced models are switched on by default.
;                      To not use the advanced models, set
;                      advanced_model=0
;
; KEYWORD PARAMETERS:
;     LOG:    If set, then the logarithm (base 10) of T_max is
;             returned.
;     INTERP: If set, then the ion fraction curve is interpolated onto
;             a logT scale with 0.01 dex steps. Note that the
;             interpolation is done on the logarithm of the ion
;             fraction values.
;     QUIET:  If set, then information messages are suppressed.
;
; OUTPUTS:
;     Returns the temperature of maximum ionization in K. If either the
;     pressure or the density are input, then the CHIANTI advanced
;     models are run. If neither are specified, then the default
;     zero-density ion balance file is used (chianti.ioneq).
;
;     If the ion balance file is not found, then a value of -1 is
;     returned.
;
; CALLS:
;     READ_IONEQ, CONVERTNAME, CH_CALC_IONEQ
;
; EXAMPLE:
;     IDL> tmax=ch_tmax('o_4',pressure=1e15)
;     IDL> tmax=ch_tmax('o_4',density=1e9)
;     IDL> tmax=ch_tmax('fe_13',advanced=0)
;     IDL> log_tmax=ch_tmax('fe_13',/log)
;     IDL> tmax=ch_tmax('fe_13',ioneqname='myioneqfile.ioneq')
;
; MODIFICATION HISTORY:
;     Ver.1, 9-May-2013, Peter Young
;         This routine is the same as a previous routine named
;         get_tmax.pro, only I've added the /LOG keyword.
;     Ver.2, 9-Aug-2019, Peter Young
;         Changed behavior so that the one-higher ionization stage is
;         used for the dielectronic "d" ions; updated header format.
;     Ver.3, 21-Jan-2021, Peter Young
;         Added /interp keyword.
;     Ver.4, 04-Nov-2024, Peter Young
;         Now uses the advanced models by default for computing Tmax.
;-

IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> tmax=ch_tmax( ion_name [, /log, ioneqname=, /interp, '
  print,'                            advanced_model=advanced_model ] )'
  print,''
  print,'  Note: the advanced models are on by default; advanced_model=0 switches them off'
  return,-1.
ENDIF 

IF n_elements(advanced_model) EQ 0 THEN adv=1b ELSE adv=advanced_model


convertname,ionname,iz,ion,diel=diel

nd=n_elements(density)
np=n_elements(pressure)

;
; If neither pressure or density are specified, then the advanced models
; are switched off.
;
IF nd EQ 0 AND np EQ 0 THEN BEGIN
  adv=0b
ENDIF

;
; Do some checks on the density and pressure inputs.
;
IF nd GT 1 THEN BEGIN
  message,/info,/CONTINUE,'Density must be a scalar. Returning...'
  return,-1.
ENDIF 
IF np GT 1 THEN BEGIN
  message,/info,/CONTINUE,'Pressure must be a scalar. Returning...'
  return,-1.
ENDIF 
IF np EQ 1 AND nd EQ 1 THEN BEGIN
  message,/info,/CONTINUE,'Please specify either pressure or density, not both. Returning...'
  return,-1.
ENDIF 

IF n_elements(ioneqname) NE 0 THEN BEGIN
  chck=file_info(ioneqname)
  IF chck.exists EQ 0 THEN BEGIN 
    message,/info,/cont,'The specified ionization balance file does not exist. Returning...'
    return,-1.
  ENDIF
ENDIF 



;
; Read the ioneq file and extract the ion fractions. Note that the advanced
; model ioneq file is the default, and the coronal approximation file
; (!ioneq_file) is only used if advanced_model=0
;
IF NOT keyword_set(adv) THEN BEGIN
  IF n_elements(ioneqname) EQ 0 THEN ioneqname=!ioneq_file
  read_ioneq,ioneqname,ioneq_t,ioneq_data,ref
ENDIF ELSE BEGIN
  ioneq_data=ch_calc_ioneq(temp, dens=density, press=pressure, ele=iz, $
                           outname=outname, quiet=quiet)
  ioneq_t=alog10(temp)
  chck=file_info(outname)
  IF chck.exists EQ 1 THEN file_delete,outname
ENDELSE 


ii=reform(ioneq_data[*,iz-1,ion-1+diel])

getmax=max(ii,index)
logtmax=ioneq_t[index]

IF keyword_set(interp) THEN BEGIN
   k=where(ii NE 0.)
   y=alog10(ii[k])
   x=ioneq_t[k]
   y2=spl_init(x,y)
  ;
   xi=findgen(41)/100.+logtmax-0.2
   k=where(xi GE min(ioneq_t) AND xi LE max(ioneq_t))
   xi=xi[k]
  ;
   yi=spl_interp(x,y,y2,xi)
   getmax=max(yi,index)
   logtmax=xi[index]
ENDIF 


IF keyword_set(log) THEN return,logtmax ELSE return,10.^logtmax

END
