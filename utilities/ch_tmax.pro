FUNCTION ch_tmax, ionname, ioneqname=ioneqname, log=log, interp=interp

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
;
; KEYWORD PARAMETERS:
;     LOG:    If set, then the logarithm (base 10) of T_max is
;             returned.
;     INTERP: If set, then the ion fraction curve is interpolated onto
;             a logT scale with 0.01 dex steps. Note that the
;             interpolation is done on the logarithm of the ion
;             fraction values.
;
; OUTPUTS:
;     Returns the temperature of maximum ionization in K. If the ion
;     balance file is not found, then a value of -1 is returned.
;
; CALLS:
;     READ_IONEQ, CONVERTNAME
;
; EXAMPLE:
;     IDL> tmax=ch_tmax('fe_13')
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
;-

IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> tmax=ch_tmax( ion_name [, /log, ioneqname=, /interp ] )'
  return,-1.
ENDIF 

IF n_elements(ioneqname) EQ 0 THEN ioneqname=!ioneq_file

IF NOT file_exist(ioneqname) THEN BEGIN
  print,'%GET_TMAX: the specified ionization balance file does not exist. Returning...'
  return,-1.
ENDIF 

read_ioneq,ioneqname,tt,ii,ref

convertname,ionname,iz,ion,diel=diel

ii=reform(ii[*,iz-1,ion-1+diel])

getmax=max(ii,index)
logtmax=tt[index]

IF keyword_set(interp) THEN BEGIN
   k=where(ii NE 0.)
   y=alog10(ii[k])
   x=tt[k]
   y2=spl_init(x,y)
  ;
   xi=findgen(41)/100.+logtmax-0.2
   k=where(xi GE min(tt) AND xi LE max(tt))
   xi=xi[k]
  ;
   yi=spl_interp(x,y,y2,xi)
   getmax=max(yi,index)
   logtmax=xi[index]
ENDIF 


IF keyword_set(log) THEN return,logtmax ELSE return,10.^logtmax

END
