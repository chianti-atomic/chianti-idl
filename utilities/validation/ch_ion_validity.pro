
FUNCTION ch_ion_validity, ionname, density, meta_val=meta_val

;+
; NAME:
;     CH_ION_VALIDITY
;
; PURPOSE:
;     Checks if the CHIANTI model for the specified ion is valid for
;     the requested density.
;
; CATEGORY:
;     CHIANTI; validation.
;
; CALLING SEQUENCE:
;     Result = CH_ION_VALIDITY( IonName, Density )
;
; INPUTS:
;     IonName:  Name of an ion in CHIANTI format (e.g., 'o_6' for O
;               VI). 
;     Density:  Electron number density in units cm^-3.
;
; OPTIONAL INPUTS:
;     Meta_Val: The population cutoff value that defines a metastable
;               level. That is, levels with populations of META_VAL or
;               higher are considered metastable. Default is 0.001.
;
; KEYWORD PARAMETERS:
;     QUIET:  If set, then information messages are not printed to the
;             screen. 
;
; OUTPUTS:
;     Returns 1 if the CHIANTI model is considered valid, and 0
;     otherwise. Information messages are also printed to the
;     screen. If a problem is found, then -1 is returned.
;
; CALLS:
;     CH_POPS, READ_SCUPS, CONVERTNAME, ZION2FILENAME
;
; PROGRAMMING NOTES:
;     The CHIANTI model is considered valid if it contains
;     electron excitations rates for all transitions out of all
;     of the ion's metastable levels. A metastable level is
;     defined to be a level with a population greater or equal to
;     META_VAL (set to 0.1%).
;
; EXAMPLE:
;     IDL> r=ch_ion_validity('fe_13',1e13)
;
; MODIFICATION HISTORY:
;     Ver.1, 09-Dec-2020, Peter Young
;-


IF n_params() LT 2 THEN BEGIN
   print,'Use:  IDL> result=ch_ion_validity(ion_name,density)'
   return,-1
ENDIF 

;
; For a level to be considered metastable, it must have a population
; that is greater or equal to META_VAL.
;
IF n_elements(meta_val) EQ 0 THEN meta_val=0.001


p=ch_pops(ionname,dens=density,/quiet)

IF n_tags(p) EQ 0 THEN BEGIN
   print,'% CH_ION_VALIDITY: this ion does not exist in CHIANTI. Returning...'
   return,-1
ENDIF


k=where(p.level.pop GE meta_val,nk)

nlvl=n_elements(p.level)
IF NOT keyword_set(quiet) THEN print,ionname+' has '+trim(nlvl)+' levels, and there are '+trim(nk)+' metastable levels for this ion.'

convertname,ionname,iz,ion
zion2filename,iz,ion,fname
read_scups,fname+'.scups',scup

swtch=bytarr(nk)

FOR i=0,nk-1 DO BEGIN
   lvl=p.level[k[i]].index
   j=where(scup.data.lvl1 EQ lvl OR scup.data.lvl2 EQ lvl,nj)
   IF nj EQ nlvl-1 THEN BEGIN
      swtch[i]=1b
   ENDIF ELSE BEGIN 
      IF NOT keyword_set(quiet) THEN print,'  Level '+trim(lvl)+': there are '+trim(nj)+' of '+trim(nlvl)+' transitions in the CHIANTI model.'
   ENDELSE 
ENDFOR

IF total(swtch) NE nk THEN BEGIN
   IF NOT keyword_set(quiet) THEN print,'The CHIANTI model is not accurate for the specified density.'
   return,0b
ENDIF ELSE BEGIN
   IF NOT keyword_set(quiet) THEN print,'The CHIANTI model should be accurate for the specified density.'
   return,1b
ENDELSE    

END
