
function RATE_COEFF, IONNAME, TEMP, TRANS=TRANS, PATH=PATH, $
                     QUIET=QUIET

;+
; NAME
;
;       RATE_COEFF()
;
; PROJECT
;
;       CHIANTI
;
; EXPLANATION
;
;       Returns the electron rate coefficient in units of cm^3 s^-1.
;       Either the entire excitation/de-excitation array for the specified
;       ion or, if the TRANS= keyword is used, the excitation rate for
;       the specified transition.
;
; INPUTS
;
;	IONNAME	The CHIANTI ion identifier, e.g., 'o_6' for O VI.
;
;	TEMP	Specify temperature(s) at which upsilon(s) are required. 
;		Note that this can be an array. Units: K.
;
; KEYWORDS
;
;       QUIET   Prevents printing of information to the screen.
;
; OPTIONAL INPUTS
;
;       TRANS   Level indices for transition, e.g., [1,2] for transition 
;               1-2.
;
;       PATH    If the data files are in a different directory from the
;               standard CHIANTI distribution, then this keyword should be
;               set to the directory. E.g., PATH='/home/mydata/o/o_6'
;
; OUTPUTS
;
;       Returns the 2D array containing excitation and de-excitation
;       coefficients. If E(i) < E(j) then CC(i,j) will be the excitation
;       coefficient and CC(j,i) will be the de-excitation coefficent.
;
;       If the keyword TRANS= is set, then the routine only returns the
;       excitation coefficient for the specified transition.
;
; EXAMPLES
;
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
; CALLS
;
;       CONVERTNAME, ZION2FILENAME, ZION2NAME, SETUP_ION, POP_SOLVER,
;       READ_IONEQ, READ_ABUND, PROTON_DENS
;
; HISTORY
;
;       Ver.1, 22-Jun-2004, Peter Young
;       Ver.2, 5-Jul-2005, Peter Young
;           updated for v.5 of CHIANTI
;-


COMMON elvlc,l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref
COMMON wgfa, wvl,gf,a_value
COMMON upsilon,splstr
COMMON radiative, radt,dilute
COMMON proton, pstr, pe_ratio
COMMON elements,abund,abund_ref,ioneq,temp_all,ioneq_ref
COMMON ionrec,rec_rate,ci_rate,temp_ionrec,luprec,lupci,status

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

setup_ion,name,-1,-1,wvltst,lvl1,lvl2,wvl1,gf1,a_value1,path=path, $
     noprot=noprot

IF NOT keyword_set(no_setup) OR n_elements(abund) EQ 0 OR $
     n_elements(ioneq) EQ 0 THEN BEGIN
  IF n_elements(ioneq_file) EQ 0 THEN ioneq_file=!ioneq_file
  read_ioneq,ioneq_file,temp_all,ioneq,ioneq_ref
  IF n_elements(abund_file) EQ 0 THEN abund_file=!abund_file
  read_abund,abund_file,abund,abund_ref
ENDIF

;
; ionization/recombination aren't needed for this routine, so switch them off
; 
ionrec = {rec:0., ci:0., temp:0., lev_up_rec:0., lev_up_ci:0., $
          status:-1, ioneq:0.}


pe_ratio=proton_dens(alog10(temp[0]))

input = {gname:name, jj:jj, ecm:ecm,ecmth:ecmth, $
 wvl:wvl, a_value:a_value, splstr:splstr, $
 pe_ratio:pe_ratio,prot_struc:pstr, dilute:0.0, radtemp:6000., $
        ionrec: ionrec}

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

