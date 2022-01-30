
FUNCTION ch_lookup_pops, ionname, dens=dens, temp=temp, ldens=ldens, ltemp=ltemp

;+
; NAME:
;     CH_LOOKUP_POPS
;
; PURPOSE:
;     Returns level populations for the specified ion in the same
;     structure format as CH_POPS, but calculated with the population
;     lookup tables.
;
; CATEGORY:
;     CHIANTI; level populations.
;
; CALLING SEQUENCE:
;     Result = CH_LOOKUP_POPS( IonName )
;
; INPUTS:
;     IonName:  The name of an ion in CHIANTI format. For example,
;               'o_5' for O V.
;
; OPTIONAL INPUTS:
;     Temp:   The temperature at which populations should be
;             calculated. If not set, then the Tmax of the ion is used
;             (calculated with ch_tmax.pro).
;     Ltemp:  Alternatively the log of the temperature can be
;             specified with this input.
;     Dens:   The density at which populations should be
;             calculated. If not set, then a value of 10^10 cm-3 is
;             used. 
;     Ldens:  Alternatively the log of the density can be specified
;             with this input.
;	
; OUTPUTS:
;      Returns a structure with the tags:
;         DENS            DOUBLE       1.0000000e+10
;         TEMP            DOUBLE           281838.17
;         LEVEL           STRUCT    -> <Anonymous> Array[86]
;         RADTEMP         FLOAT           10000.0
;         RPHOT           FLOAT           0.00000
;         PROTON          STRING    'yes'
;         VERSION         STRING    'CHIANTI 8.0.2'
;         DATE            STRING    'Tue May 17 11:01:38 2016'
;         SUM_MWL         INT              0
;         SUM_MWL_COEFFS  FLOAT          -1.00000
;
;      If a problem is found then a value of -1 is returned.
;
; EXAMPLE:
;     IDL> p=ch_lookup_pops('o_5')
;     IDL> p=ch_lookup_pops('fe_13',ltemp=6.4,ldens=8.5)
;
; MODIFICATION HISTORY:
;     Ver.1, 29-Jan-2022, Peter Young
;-


IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> output = ch_lookup_pops( IonName [, dens=, ldens=, temp=, ltemp= ] )'
  return,-1
ENDIF 

dir=getenv('CHIANTI_LOOKUP')
chck=file_info(dir)
IF chck.exists EQ 0 THEN BEGIN
  print,'% CH_LOOKUP_POPS: The environment variable $CHIANTI_LOOKUP is not defined. Returning...'
  return,-1
ENDIF

nlt=n_elements(ltemp)
nt=n_elements(temp)
nld=n_elements(ldens)
nd=n_elements(dens)

swtch=0
IF nlt GT 1 OR nt GT 1 THEN swtch=1
IF nld GT 1 OR nd GT 1 THEN swtch=1
IF swtch EQ 1 THEN BEGIN
  print,'% CH_LOOKUP_POPS: density and temperature must be specified as scalars. Returning...'
  return,-1
ENDIF 

;
; Set default parameters.
;
IF nlt EQ 0 AND nt EQ 0 THEN temp=ch_tmax(ionname)
IF nld EQ 0 AND nd EQ 0 THEN dens=1e10

IF n_elements(temp) NE 0 THEN temp_in=temp
IF n_elements(dens) NE 0 THEN dens_in=dens
;
IF n_elements(ltemp) NE 0 THEN temp_in=10.^ltemp
IF n_elements(ldens) NE 0 THEN dens_in=10.^ldens


p=ch_lookup_table_interp(ionname,dens_in,temp_in,/quiet)


nlev=max(p.levels)

str={index: 0, term: '', pop: 0d0}

lev=replicate(str,nlev)
lev.index=indgen(nlev)+1
lev[p.levels-1].pop=reform(p.pop[*,*,p.levels-1])

IF total(lev.pop) EQ 0. THEN BEGIN
  print,'% CH_LOOKUP_POPS: The specified parameters are outside the range of validity of the lookup tables. '
  print,'                  Try using ch_pops.pro instead.'
  return,-1
ENDIF 

;
; Note that radtemp, rphot and sum_mwl are not implemented through the
; lookup tables, but proton rates are.
;
popstr={dens: dens_in, temp: temp_in, level: lev, radtemp: -1., $
        rphot: -1., proton: 'yes', version: ch_get_version(), date: systime(), $
        sum_mwl: 0, sum_mwl_coeffs: -1.}


return,popstr

END
