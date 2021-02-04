
FUNCTION ch_ion_validity, ionname, density, meta_val=meta_val, path=path, $
                          quiet=quiet, verbose=verbose

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
;               higher are considered metastable. Default is 0.0001.
;     Path:     Specify the path to the ion's data
;               files. Over-rides the default CHIANTI location. 
;
; KEYWORD PARAMETERS:
;     QUIET:  If set, then information messages are not printed to the
;             screen.
;     VERBOSE: If set, then all missing transitions are printed to the
;              screen. 
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
;     Ver.2, 29-Dec-2020, Peter Young
;       Modified so that the routine only checks if electron
;       excitation data is available for transitions to levels below
;       the ionization threshold. This is specifically for dealing
;       with the two-ion models where there are many levels above
;       threshold. 
;     Ver.3, 04-Feb-2021, Peter Young
;        Changed default meta_val to 10^-4; added /verbose option. 
;-


IF n_params() LT 2 THEN BEGIN
   print,'Use:  IDL> result=ch_ion_validity(ion_name,density)'
   return,-1
ENDIF 

;
; For a level to be considered metastable, it must have a population
; that is greater or equal to META_VAL.
;
IF n_elements(meta_val) EQ 0 THEN meta_val=0.0001


p=ch_pops(ionname,dens=density,/quiet,path=path)

IF n_tags(p) EQ 0 THEN BEGIN
   print,'% CH_ION_VALIDITY: this ion does not exist in CHIANTI. Returning...'
   return,-1
ENDIF


k=where(p.level.pop GE meta_val,nk)

nlvl=n_elements(p.level)
IF NOT keyword_set(quiet) THEN print,ionname+' has '+trim(nlvl)+' levels, and there are '+trim(nk)+' metastable levels for this ion.'

IF n_elements(path) EQ 0 THEN BEGIN 
   convertname,ionname,iz,ion
   zion2filename,iz,ion,fname
   scups_name=fname+'.scups'
   elvlc_name=fname+'.elvlc'
ENDIF ELSE BEGIN
   scups_name=concat_dir(path,trim(ionname)+'.scups')
   elvlc_name=concat_dir(path,trim(ionname)+'.elvlc')
ENDELSE
read_elvlc,elvlc_name,elvlc=elvlc
read_scups,scups_name,scup
   
;
; Get ionization potential in cm^-1.
;
ip=ch_ip(ionname,/cm)

;
; Extract the levels that are below the ionization threshold (ip).
;
en=elvlc.data.energy
en_data=elvlc.data
k_ip=where(en_data.energy LT ip,n_ip)
en_data=en_data[k_ip]

IF NOT keyword_set(quiet) THEN print,ionname+' has '+trim(n_ip)+' levels below the ionization threshold.'

;
; This creates a string array that I use for checking missing
; transitions later.
;
lvl=[[scup.data.lvl1],[scup.data.lvl2]]
lvl=minmax(lvl,dimen=2)
lvl_check=trim(reform(lvl[0,*]))+'-'+trim(reform(lvl[1,*]))

;
; For each metastable level (i-loop) check all of the levels below the
; ionization threshold (j-loop), and see if scups data are available
; for each of these levels. If yes, then swtch=1 for this metastable. 
; 
swtch=bytarr(nk)
FOR i=0,nk-1 DO BEGIN
   lvl=p.level[k[i]].index
   lvl_chck=bytarr(n_ip)
   FOR j=0,n_ip-1 DO BEGIN
      lvl_ip=en_data[j].index
      IF lvl_ip NE lvl THEN BEGIN 
         str1=trim(min([lvl,lvl_ip]))+'-'+trim(max([lvl,lvl_ip]))
         chck=where(str1 EQ lvl_check,nchck)
         IF nchck GT 0 THEN BEGIN
            lvl_chck[j]=1b
         ENDIF ELSE BEGIN
            IF keyword_set(verbose) THEN print,format='("  Missing transition: ",i5," -",i5)',lvl,lvl_ip
         ENDELSE 
      ENDIF
   ENDFOR
   n_lvl_chck=total(lvl_chck)
   IF n_lvl_chck EQ n_ip-1 THEN BEGIN
      swtch[i]=1b
   ENDIF ELSE BEGIN 
      IF NOT keyword_set(quiet) THEN print,'  Level '+trim(lvl)+': there are '+trim(n_lvl_chck)+' of '+trim(n_ip-1)+' transitions in the CHIANTI model.'
   ENDELSE 
ENDFOR 
   


IF total(swtch) NE nk THEN BEGIN
   IF NOT keyword_set(quiet) THEN BEGIN
      print,''
      print,'The CHIANTI model may not be accurate for the specified density.'
      print,'Depending on how many transitions are missing, and for which levels, the model'
      print,'may still be accurate. Use the /verbose keyword to print the list of missing'
      print,'transitions.'
   ENDIF 
   return,0b
ENDIF ELSE BEGIN
   IF NOT keyword_set(quiet) THEN print,'The CHIANTI model should be accurate for the specified density.'
   return,1b
ENDELSE    


END
