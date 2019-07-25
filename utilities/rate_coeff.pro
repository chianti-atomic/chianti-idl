
function RATE_COEFF, IONNAME, TEMP, TRANS=TRANS, PATH=PATH, $
                     QUIET=QUIET

;+
; NAME:
;       RATE_COEFF()
;
; PURPOSE:
;       Returns the electron rate coefficient in units of cm^3 s^-1.
;       Either the entire excitation/de-excitation array for the specified
;       ion or, if the TRANS= keyword is used, the excitation rate for
;       the specified transition.
;
; CATEGORY:
;       CHIANTI; rate coefficient.
;
; CALLING SEQUENCE:
;       IDL> OUTPUT=RATE_COEFF( IONNAME, TEMP )
;
; INPUTS:
;	IonName: The CHIANTI ion identifier, e.g., 'o_6' for O VI.
;
;	Temp: 	 Specify temperature(s) at which the rate
;	         coefficient(s) is required. Note that this can be an
;	         array. Units: K. 
;
; KEYWORDS PARAMETERS:
;       QUIET   Prevents printing of information to the screen.
;
; OPTIONAL INPUTS:
;       Trans:  Level indices for transition, e.g., [1,2] for transition 
;               1-2.
;
;       Path:   If the data files are in a different directory from the
;               standard CHIANTI distribution, then this keyword should be
;               set to the directory. E.g., PATH='/home/mydata/o/o_6'
;
; OUTPUTS:
;       Returns the 2D array containing excitation and de-excitation
;       coefficients. If E(i) < E(j) then CC(i,j) will be the excitation
;       coefficient and CC(j,i) will be the de-excitation coefficent.
;
;       If the keyword TRANS= is set, then the routine only returns the
;       excitation coefficient for the specified transition.
;
; EXAMPLES:
;       IDL> help,rate_coeff('o_6',3e5)
;       <Expression>    DOUBLE    = Array[40, 40]
;
;       IDL> result=rate_coeff('o_6',3e5,trans=[1,3])
;        Wavelength:            1031.91
;        Exc. rate coeff:    1.913e-008
;        De-exc. rate coeff: 1.522e-008
;       IDL> print,result
;         1.9130775e-008
;
; CALLS:
;       CONVERTNAME, ZION2FILENAME, ZION2NAME, CH_SETUP_ION, POP_SOLVER
;
; MODIFICATION HISTORY:
;       Ver.1, 22-Jun-2004, Peter Young
;       Ver.2, 5-Jul-2005, Peter Young
;           updated for v.5 of CHIANTI
;       Ver.3, 30-Jan-2018, Peter Young
;           added call to ch_setup_ion and removed common blocks;
;           tidied up code and header.
;-


IF n_params() LT 2 THEN BEGIN
  print,'Use:  IDL> result=rate_coeff(ionname, temp, [trans= , path= ,'+ $
       '/quiet]'
  return,-1.
ENDIF

convertname,ionname,iz,ion
zion2filename,iz,ion,filename,name=name,diel=diel
IF N_ELEMENTS(path) NE 0 THEN BEGIN
  zion2name,iz,ion,name
  filename=concat_dir(path,name)
ENDIF


;
; Load the atomic data into INPUT.
;
input=ch_setup_ion(name,rphot=rphot,radtemp=radtemp,noprot=noprot, $
                   ioneq_file=ioneq_file,abund_file=abund_file,path=path)
wvl=input.wvl


;
; the density does not affect the rate coefficient, and so I give an
; arbitrary value of 10^9.
;
dens=1d9
pop_solver,input, temp[0],dens,pop, data_str=data_str

cc=data_str.cc/dens

IF n_elements(trans) EQ 0 THEN return,cc

i=trans[0]-1
j=trans[1]-1
exc=cc[i,j]
deexc=cc[j,i]

IF NOT keyword_set(quiet) THEN BEGIN
  PRINT,FORMAT='(" Wavelength:         ",5f10.2)',wvl[i,j]
  PRINT,FORMAT='(" Exc. rate coeff:    ",5e10.3)',exc
  PRINT,FORMAT='(" De-exc. rate coeff: ",5e10.3)',deexc
ENDIF

return,exc

END

