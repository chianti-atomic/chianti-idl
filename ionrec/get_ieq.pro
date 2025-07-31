

FUNCTION get_ieq, temp, iz, ion, ioneq_logt=ioneq_logt, ioneq_frac=ioneq_frac, $
                  ioneq_name=ioneq_name, dielectronic=dielectronic

;+
; NAME
;      GET_IEQ()
;
; PURPOSE:
;      For a specified ion (IZ, ION) and set of temperatures (TEMP) this 
;      routine takes the ion fraction values tabulated in one of the CHIANTI 
;      .IONEQ files, interpolates and extracts the values of the ion 
;      fraction at the input temperatures.
;      
; CATEGORY:
;      CHIANTI; ionization fractions; interpolation.
;
; INPUTS:
;      Temp:  The temperature(s) at which the ion fractions are
;             required.
;      Iz:    The name of an ion in CHIANTI format (e.g., 'fe_13'). Or
;             it can be the atomic number of the element (integer). If
;             the latter, then ION must be specified.
;
; OPTIONAL INPUTS:
;      Ion:   The spectroscopic number of the ion (e.g., 13 =
;             XIII). This is only used if IZ has been specified as an
;             integer. 
;      Ioneq_LogT: The temperature output from the READ_IONEQ routine.
;      Ioneq_Frac: The ion fractions from the READ_IONEQ routine.
;      Ioneq_Name: The name of an ionization balance file in CHIANTI format.
;
; OUTPUT:
;      A vector of same length as the input TEMP containing the ion 
;      fractions at these temperatures. If the requested temperatures lie
;      outside of the temperature range of the ioneq file, then ion
;      fractions at those temperatures will be zero.
;
;      If a problem is found, then a value of -1 is returned.
;
; CALLS:
;      READ_IONEQ, CONVERTNAME
;
; EXAMPLES:
;      IDL> ltemp=findgen(101)/100+5.2
;      IDL> y=get_ieq(10.^ltemp,'o_6')
;
;      IDL> y=get_ieq(10.^ltemp,8,6)
;
;      IDL> read_ioneq,!ioneq_file,ioneq_logt,ioneq_frac,ref
;      IDL> y=get_ieq(10.^ltemp,'o_6',ioneq_logt=ioneq_logt,ioneq_frac=ioneq_frac)
;
; MODIFICATION HISTORY:
;      Ver.1, 24-Jul-2002, Peter Young
;      Ver.2, 14-Sept-2005, GDZ
;         Modified error message.
;      Ver.3, 24-Jul-2018, Peter Young
;         Now uses !ioneq_file by default; added ioneq_name= optional
;         input; now accepts the ion name as an input.
;      Ver.4, 31-Jul-2025, Peter Young
;         Now handles the case when the ioneq file contains only one
;         temperature. This became a problem in CHIANTI 11. Also swapped datatype
;         for isa. Added /dielectronic keyword.
;-

IF n_params() LT 2 THEN BEGIN
  print,'Use:'
  print,'   Either: IDL> result=get_ieq(temp,ion_name [,ioneq_logt= , ioneq_frac=, ioneq_name=, dielectronic= ]'
  print,''
  print,'       Or: IDL> result=get_ieq(temp,iz,ion [,ioneq_logt= , ioneq_frac=, ioneq_name=, dielectronic= ]'
  return,-1.
ENDIF

;
; Checks if IZ is a string, in which case it is interpreted as the ion
; name. 
;
IF isa(iz,/string) THEN BEGIN
  convertname,iz,ion_z,ion_n,dielectronic=diel
  IF diel EQ 1 THEN dielectronic=1b
ENDIF ELSE BEGIN
  ion_z=iz
  ion_n=ion
ENDELSE 

IF NOT keyword_set(dielectronic) THEN dielectronic=0b

IF n_elements(ioneq_name) EQ 0 THEN BEGIN 
  ioneq_name=!ioneq_file
ENDIF ELSE BEGIN
  chck=file_search(ioneq_name,count=count)
  IF count EQ 0 THEN BEGIN
    message,/info,/cont,'the specified ioneq file was not found. Returning...'
    return,-1.
  ENDIF 
ENDELSE 
  
IF n_elements(ioneq_logt) EQ 0 OR n_elements(ioneq_frac) EQ 0 THEN BEGIN 
  read_ioneq,ioneq_name,ioneq_logt,ioneq_frac,ioneq_ref 
END 

nt=n_elements(temp)
ltemp=alog10(temp)

answer=dblarr(nt)
this_ioneq=ioneq_frac[*,ion_z-1,ion_n-1+dielectronic]

;
; The following checks if the input temperatures exactly match those in
; ioneq_logt. If they do, then we simply return the ion fractions for the
; ioneq_logt values without doing an interpolation. This is to get round
; a problem in CHIANTI 11 whereby it's possible for the ioneq file to
; contain only one temperature.
;
; Note: "exactly match" means within 0.001 dex of the requested
; temperature.
;
indx=intarr(nt)-1
FOR i=0,nt-1 DO BEGIN
  getmin=min(abs(ltemp[i]-ioneq_logt),imin)
  IF getmin LT 0.001 THEN BEGIN
    indx[i]=imin
    answer[i]=this_ioneq[imin]
  ENDIF 
ENDFOR
;
; If all of the temperatures match, then return.
;
IF min(indx) GE 0 THEN return,answer

;
; If the ioneq file only has one temperature, then exit since the spline
; interpolation will not work.
;
IF n_elements(ioneq_logt) EQ 1 THEN return,answer

;
; Exit if the ioneq array is empty.
;
i=where(this_ioneq NE 0.)
IF i[0] EQ -1 THEN return,dblarr(nt)

x=ioneq_logt[i]
y=alog10(this_ioneq[i])

ind=where(ltemp GE min(x) AND ltemp LE max(x))
IF ind[0] EQ -1 THEN return,dblarr(nt)

;
; Perform spline fit to get the ioneq values.
;
xi=ltemp[ind]
y2=spl_init(x,y)
yi=spl_interp(x,y,y2,xi)
answer[ind]=10.^yi

return,answer

END
