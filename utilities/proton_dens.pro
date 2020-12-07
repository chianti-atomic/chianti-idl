;+
; NAME:
;    PROTON_DENS()
;
; PURPOSE:
;    Calculates the ratio of the proton density to electron density using 
;    abundance and ion balance files.
;
; CATEGORY:
;    CHIANTI; protons.
;
; CALLING SEQUENCE:
;    Result = PROTON_DENS( Log_Temp )
;
; INPUTS:
;    Temp:   The logarithm (base 10) of the temperature(s) for which the 
;            ratio is required. Can be an array.
;
; OPTIONAL INPUTS:
;    Abund_file:   The name of an element abundance file, in CHIANTI
;                  format. This over-rides (and replaces) any
;                  abundances in the common block.
;    Ioneq_file:   The name of an ionization fraction file, in CHIANTI
;                  format. This over-rides (and replaces) any
;                  ion fractions in the common block.
;	
; KEYWORD PARAMETERS:
;    HYDROGEN: If set then the routine computes the ratio of hydrogen to 
;              free electrons.
;    QUIET:  If set, then information messages will not be printed to
;            the IDL window.
;
; OUTPUTS:
;    An array of same size as TEMP containing the proton-to-electron
;    ratio. If /hydrogen is set, then the hydrogen-to-electron ratio
;    is returned. If a problem is found, then a value of -1 is
;    returned. 
;
; CALLS:
;    READ_IONEQ, READ_ABUND
;
; PROGRAMMING NOTES:
;    The proton/electron ratio is computed using the ionization
;    fractions (ioneq) and the element abundance (abund) file. If the
;    names of these files are not specified, then !ioneq_file and
;    !abund_file are used.
;
;    It is possible to input the ioneq and abund information through
;    the common block, but this is not recommended.
;
;    If a temperature is requested that is outside the range of
;    validity of the ioneq file, then the values at the ends of the
;    range will be returned.
;
;    I've added a check to see if the temperatures are tabulated at 0.1 
;    dex intervals (e.g., 4.0, 4.1, etc.). If they are, then a quicker 
;    algorithm is used to calculate p/e ratios. This is useful for 
;    synthetic.pro.
;
; EXAMPLE:
;    IDL> log_temp=findgen(41)/10.+4.0
;    IDL> np_ne=proton_dens(log_temp)
;    IDL> nh_ne=proton_dens(log_temp,/hydrogen)
;
; MODIFICATION HISTORY:
;    Ver.1, 5-Dec-2001, Peter Young
;
;    Ver.2, 3-Dec-2001, Peter Young
;        Added /hydrogen keyword.
;
;    V. 3, 22-May-2002, GDZ: 
;                  generalize directory concatenation to work for Unix, Windows
;                  and VMS.
;
;    V.4, 06-Aug-02 GDZ
;              Changed the use of CHIANTI system variables.
;
;    Ver.5, 16-Aug-2016, Peter Young
;              Introduced IONEQ_FILE and ABUND_FILE optional inputs.
;    Ver.6, 07-Dec-2020, Peter Young
;       Added check on input parameters; removed the common block;
;       updated header; added /quiet keyword.
;-


FUNCTION  proton_dens, temp, hydrogen=hydrogen, abund_file=abund_file, $
                       ioneq_file=ioneq_file, quiet=quiet


IF n_params() LT 1 THEN BEGIN
   print,'Use:  IDL> p_e=proton_dens( log_temp [, /hydrogen, abund_file=, ioneq_file=, /quiet ] )'
   print,''
   print,'   Returns the proton-to-electron ratio for specified temperatures.'
   print,'     /hydrogen  - used to return hydrogen-to-electon ratio'
   return,-1
ENDIF 

;
; The logic below is as follows:
;    - if abund_file has been specified and it exists, then use this
;    - otherwise, if the common block has been defined, then use this
;    - if neither of above, then read the default abundance file
;
IF n_elements(abund_file) NE 0 THEN BEGIN
  chck=file_search(abund_file,count=count)
ENDIF ELSE BEGIN
  count=0
ENDELSE 
;
IF count GT 0 THEN BEGIN
  read_abund, abund_file, abund,abund_ref
ENDIF ELSE BEGIN
  IF n_elements(abund) EQ 0 THEN BEGIN
    read_abund, !abund_file, abund,abund_ref
    IF NOT keyword_set(quiet) THEN print,'% PROTON_DENS: Using the default abundance file '+!abund_file +' to calculate the p/e ratio'
  ENDIF 
ENDELSE 


;
; The logic below is as follows:
;    - if ioneq_file has been specified and it exists, then use this
;    - otherwise, if the common block has been defined, then use this
;    - if neither of above, then read the default ion fraction file
;
IF n_elements(ioneq_file) NE 0 THEN BEGIN
  chck=file_search(ioneq_file,count=count)
ENDIF ELSE BEGIN
  count=0
ENDELSE 
;
IF count GT 0 THEN BEGIN
  read_ioneq,ioneq_file, ioneq_logt,ioneq,ioneq_ref
ENDIF ELSE BEGIN
  IF n_elements(ioneq) EQ 0 THEN BEGIN 
    read_ioneq,!ioneq_file, ioneq_logt,ioneq,ioneq_ref
    IF NOT keyword_set(quiet) THEN print,'% PROTON_DENS: Using the default ion fraction file '+!ioneq_file+' to calculate the p/e ratio'
  ENDIF
ENDELSE 


temp=double(temp)
nt=n_elements(temp)

siz=size(ioneq)
n_ion=siz[3]

;
; tst1 is set to 1 if the temperatures are found NOT to be on (roughly) 0.1
; dex intervals.
;
tst1=0
index=0
FOR i=0,nt-1 DO BEGIN
   getmin=min( abs(temp[i]-ioneq_logt), ind)
   index=[index,ind]
   IF getmin GT 0.001 THEN tst1=1
ENDFOR
index=index[1:*]


IF tst1 EQ 0 THEN BEGIN
                                ;
   nprot=abund[0]*ioneq[index,0,1]
   nh=nprot+abund[0]*ioneq[index,0,0]
                                ;
   i_a=where(abund NE 0.,nia)
   nelec=0d0
   IF n_elements(index) EQ 1 THEN BEGIN
      FOR i=0,nia-1 DO BEGIN
         nelec=nelec+abund[i_a[i]]* $
           (transpose(reform(ioneq[index,i_a[i],*])) # findgen(n_ion))
      ENDFOR
   ENDIF ELSE BEGIN
      FOR i=0,nia-1 DO BEGIN
         nelec=nelec+abund[i_a[i]]* $
           (reform(ioneq[index,i_a[i],*]) # findgen(n_ion))
      ENDFOR
   ENDELSE

   IF keyword_set(hydrogen) THEN return,nh/nelec ELSE return,nprot/nelec
                                ;
ENDIF ELSE BEGIN
                                ;
   edens=dblarr(nt)
   pdens=dblarr(nt)
   FOR i=0,nt-1 DO BEGIN        ; loop through each temperature
      t=temp[i]
      getmin=min( abs(ioneq_logt-t), it )
      FOR j=0,siz[2]-1 DO BEGIN ; loop through each element
         ind=where(reform(ioneq[it,j,*]) NE 0.)
         IF ind[0] NE -1 THEN BEGIN
            ni=n_elements(ind)
            FOR k=0,ni-1 DO BEGIN ; loop through each ion
               y=reform(ioneq[*,j,ind[k]])
               ind2=where(y NE 0.)
                                ;
                                ; fit spline to log of ion balance data
                                ;
               y=alog10(y[ind2])
               x=ioneq_logt[ind2]
               y2=spl_init(x,y)
               yi=spl_interp(x,y,y2,t)
                                ;
                                ; for electron density, need elt. abundance, ion fraction, and 
                                ; charge of ion
                                ;
               edens[i]=edens[i]+abund[j]*10.^yi*double(ind[k])
                                ;
                                ; for proton density, just need number of protons (ionized H); if 
                                ; /hydrogen set, though, then also add number of neutral hydrogen 
                                ; atoms.
                                ;
               IF (j EQ 0) AND (ind[k] EQ 1) THEN pdens[i]=pdens[i]+abund[j]*10.^yi
               IF keyword_set(hydrogen) AND (j EQ 0) AND (ind[k] EQ 0) THEN $
                 pdens[i]=pdens[i]+abund[j]*10.^yi
            ENDFOR
         ENDIF
      ENDFOR
   ENDFOR
   return,pdens/edens
ENDELSE

END
