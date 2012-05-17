
;+
; NAME
;
;       SPL2UPS()
;
; PROJECT
;
;       CHIANTI
;
; EXPLANATION
;
;       Collision strength data are stored in CHIANTI in the form of spline
;       fits to the Maxwellian-averaged collision strength (or upsilon). This
;       routine will retrieve the upsilons for the transition TRANS of a
;       specified ion IONNAME, at the specified temperatures TEMP.
;
; INPUTS
;
;	IONNAME The name of the ion in CHIANTI format. E.g., 'fe_13'.
;
;       TRANS   Level indices for transition, e.g., [1,2] for transition 
;               1-2.
;
;	TEMP	Specify temperature(s) at which upsilon(s) are required. 
;		Note that this can be an array. The routine will
;		accept temperatures (in K) and log temperatures. It
;		does this by checking if the max value of TEMP is
;		greater or less than 12.
;		
;
; OPTIONAL INPUTS
;
;	PATH	Directly specify the directory containing the .splups 
;		file
;
; OPTIONAL OUTPUTS
;
;       DE      Energy for transition, Rydbergs
;
; KEYWORDS
;
;	QUIET	If this is set, then the routine will not print anything 
;		to screen.
;
;       PROT    If set will return the proton rate coefficients.
;
; EXAMPLES
;
;
; CALLS
;
;       READ_SPLUPS, ZION2FILENAME, ZION2NAME, DESCALE_ALL, READ_SPLUPS,
;       CONVERTNAME
;
; HISTORY
;
;	Ver 1, 14-Aug-2006, Peter Young
;           Re-written as function rather than procedure.
;
; CONTACT
;
;	Peter Young, NRL/GMU
;-


FUNCTION SPL2UPS, IONNAME, TRANS, TEMP, PATH=PATH, QUIET=QUIET, $
             PROT=PROT, DE=DE

IF N_PARAMS() LT 3 THEN BEGIN
  PRINT,'Use:  IDL> ups = spl2ups (ionname, trans, temp, $'
  PRINT,'                     path=, /quiet, /prot, de=)'
  RETURN,-1.
ENDIF

convertname,ionname,iz,ion

IF KEYWORD_SET(path) THEN BEGIN
  zion2name,iz,ion,name
  filename=path+'/'+name
ENDIF ELSE BEGIN
  zion2filename,iz,ion,filename
ENDELSE

IF keyword_set(prot) THEN BEGIN
  splupsfile=filename+'.psplups'
  result=findfile(splupsfile)
  IF result[0] EQ '' THEN BEGIN
    IF NOT keyword_set(quiet) THEN print,'** The .psplups file does not exist for this ion. **'
    return,-1.
  ENDIF
ENDIF ELSE BEGIN
  splupsfile=filename+'.splups'
  result=findfile(splupsfile)
  IF result[0] EQ '' THEN BEGIN
    IF NOT keyword_set(quiet) THEN print,'** The .splups file does not exist for this ion. **'
    return,-1.
  ENDIF
ENDELSE

ll=trans[0]
ul=trans[1]

read_splups,splupsfile,splstr,splref,prot=prot

ind=where((splstr.lvl1 EQ ll) AND (splstr.lvl2 EQ ul))
IF ind[0] EQ -1 THEN BEGIN
  IF NOT keyword_set(quiet) THEN print,'** No data exists for this transition **'
  data=-1.
  return,-1.
ENDIF

ind=reverse(ind)
ind=ind[0]

eij=splstr[ind].de
k=splstr[ind].t_type
c=splstr[ind].c_ups
nspl=splstr[ind].nspl

n=n_elements(temp)
IF max(temp) LT 12 THEN temp_ten = 10.^temp ELSE temp_ten = temp

descale_all,temp_ten,splstr,ind,ups

de=eij

return,ups

END
