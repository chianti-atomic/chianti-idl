
FUNCTION SPL2UPS, IONNAME, TRANS, TEMP, PATH=PATH, QUIET=QUIET, $
                  PROT=PROT, DE=DE


;+
; NAME:
;       SPL2UPS()
;
; PURPOSE:
;       Collision strength data are stored in CHIANTI in the form of spline
;       fits to the Maxwellian-averaged collision strength (or upsilon). This
;       routine will retrieve the upsilons for the transition TRANS of a
;       specified ion IONNAME, at the specified temperatures TEMP.
;
; CATEGORY:
;       CHIANTI; data access.
;
; CALLING SEQUENCE:
;	Result = SPL2UPS( IonName, Trans, Temp )
;
; INPUTS:
;	IonName:The name of the ion in CHIANTI format. E.g., 'fe_13'.
;
;       Trans:  Level indices for transition, e.g., [1,2] for transition 
;               1-2.
;
;	Temp:   Specify temperature(s) at which upsilon(s) are required. 
;		Note that this can be an array. The routine will
;		accept temperatures (in K) and log temperatures. It
;		does this by checking if the max value of TEMP is
;		greater or less than 12.
;
; OPTIONAL INPUTS:
;	Path:	Directly specify the directory containing the .scups 
;		file
;
; KEYWORD PARAMETERS:
;	QUIET:  If this is set, then the routine will not print anything 
;		to screen.
;
;       PROT:   If set will return the proton rate coefficients.
;
; OUTPUTS:
;       Returns the Maxwellian-averaged collision strengths (upsilons)
;       at the requested temperatures. If a problem is found, then a
;       value of -1 is returned.
;
; OPTIONAL OUTPUTS:
;       De:     Energy for transition, Rydbergs
;
; CALLS:
;       READ_SPLUPS, ZION2FILENAME, ZION2NAME, DESCALE_ALL, READ_SCUPS,
;       CONVERTNAME, DESCALE_SCUPS
;
; EXAMPLE:
;       IDL> ltemp=findgen(21)/10.+5.2
;       IDL> ups=spl2ups('fe_13',[1,20],10.^ltemp)
;
; MODIFICATION HISTORY:
;	Ver 1, 14-Aug-2006, Peter Young
;           Re-written as function rather than procedure.
;       Ver 2, 20-Jun-2014, Peter Young
;           Modified to call the new scups routines for electron
;           data. For proton data the old procedure is retained.
;       Ver. 3, 12-Nov-2020, Peter Young
;           Now handles the dielectronic files; updated header format.
;-





IF N_PARAMS() LT 3 THEN BEGIN
  PRINT,'Use:  IDL> ups = spl2ups (ionname, trans, temp, $'
  PRINT,'                     path=, /quiet, /prot, de=)'
  RETURN,-1.
ENDIF

convertname,ionname,iz,ion,diel=diel

IF KEYWORD_SET(path) THEN BEGIN
  zion2name,iz,ion,name,diel=diel
  filename=path+'/'+name
ENDIF ELSE BEGIN
  zion2filename,iz,ion,filename,diel=diel
ENDELSE

ll=trans[0]
ul=trans[1]

IF keyword_set(prot) THEN BEGIN
 ;
 ; Descale the proton data
 ;
  splupsfile=filename+'.psplups'
  result=findfile(splupsfile)
  IF result[0] EQ '' THEN BEGIN
    IF NOT keyword_set(quiet) THEN print,'** The .psplups file does not exist for this ion. **'
    return,-1.
  ENDIF
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

ENDIF ELSE BEGIN
 ;
 ; Descale the electron data
 ;
  splupsfile=filename+'.scups'
  result=findfile(splupsfile)
  IF result[0] EQ '' THEN BEGIN
    IF NOT keyword_set(quiet) THEN print,'** The .splups file does not exist for this ion. **'
    return,-1.
  ENDIF
  read_scups,splupsfile,splstr

  ind=where((splstr.data.lvl1 EQ ll) AND (splstr.data.lvl2 EQ ul))
  IF ind[0] EQ -1 THEN BEGIN
    IF NOT keyword_set(quiet) THEN print,'** No data exists for this transition **'
    data=-1.
    return,-1.
  ENDIF

  ind=reverse(ind)
  ind=ind[0]

  eij=splstr.data[ind].de
  k=splstr.data[ind].t_type
  c=splstr.data[ind].c_ups
  nspl=splstr.data[ind].nspl

  n=n_elements(temp)
  IF max(temp) LT 12 THEN temp_ten = 10.^temp ELSE temp_ten = temp

  descale_scups,temp_ten,splstr,ind,ups

  de=eij

  return,ups

ENDELSE


END
