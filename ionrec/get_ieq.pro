

FUNCTION get_ieq, temp, iz, ion, ioneq_logt=ioneq_logt, ioneq_frac=ioneq_frac

;+
; PROJECT:  CHIANTI
;
;      CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;      Astrophysical Plasmas. It is a collaborative project involving the Naval
;      Research Laboratory (USA), the University of Florence (Italy), the
;      University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
; NAME
;
;      GET_IEQ()
;
; EXPLANATION
;
;      For a specified ion (IZ, ION) and set of temperatures (TEMP) this 
;      routine takes the ion fraction values tabulated in one of the CHIANTI 
;      .IONEQ files, interpolates and extracts the values of the ion 
;      fraction at the input temperatures.
;
; INPUTS
;
;      TEMP   The temperature(s) at which the ion fractions are required.
;
;      IZ     The atomic number of the element (e.g., 26 = iron).
;
;      ION    The spectroscopic number of the ion (e.g., 13 = XIII).
;
; OPTIONAL INPUTS
;
;      IONEQ_LOGT  The temperature output from the READ_IONEQ routine.
;
;      IONEQ_FRAC  The ion fractions from the READ_IONEQ routine.
;
; OUTPUT
;
;      A vector of same length as the input TEMP containing the ion 
;      fractions at these temperatures.
;
; CALLS
;
;      READ_IONEQ
;
; HISTORY
;
;      Ver.1, 24-Jul-2002, Peter Young
;      v.2, 14-Sept-2005, GDZ
;      modified error message.
;
;
; VERSION 2 14-Sept-2005
;-

IF n_params() LT 3 THEN BEGIN
  print,'Use:  IDL> result=get_ieq(temp,iz,ion [,ioneq_logt= , ioneq_frac= ]'
  return,-1.
ENDIF

IF n_elements(ioneq_logt) EQ 0 OR n_elements(ioneq_frac) EQ 0 THEN BEGIN 
  dir=concat_dir(!xuvtop,'ioneq')
  ioneq_name=ch_get_file(path=dir,filter='*.ioneq', $
                         title='Select Ionization Equilibrium File')
  read_ioneq,ioneq_name,ioneq_logt,ioneq_frac,ioneq_ref 
END 

nt=n_elements(temp)
ltemp=alog10(temp)

answer=dblarr(nt)

this_ioneq=ioneq_frac[*,iz-1,ion-1]
i=where(this_ioneq NE 0.)
IF i[0] EQ -1 THEN return,dblarr(nt)

x=ioneq_logt[i]
y=alog10(this_ioneq[i])

ind=where(ltemp GE min(x) AND ltemp LE max(x))
IF ind[0] EQ -1 THEN return,dblarr(nt)

xi=ltemp[ind]

y2=spl_init(x,y)
yi=spl_interp(x,y,y2,xi)
answer[ind]=10.^yi

return,answer

END
