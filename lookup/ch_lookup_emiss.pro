
FUNCTION ch_lookup_emiss, iz, ion, temp=temp, ltemp=ltemp, dens=dens, ldens=ldens, $
                          dir_lookup=dir_lookup, no_de=no_de, quiet=quiet, $
                          pressure=pressure, verbose=verbose

;+
; NAME:
;     CH_LOOKUP_EMISS
;
; PURPOSE:
;     Reproduces the behavior of EMISS_CALC, but using lookup tables.
;     The emissivity is defined as E_ij*A_ji*n_j, where E_ij is the
;     transition energy, A_ji is the radiative decay rate, and n_j is
;     the level population relative to the ion population.
;
; CATEGORY:
;     CHIANTI; emissivity.
;
; CALLING SEQUENCE:
;	Result = EMISS_CALC( IZ )
;
; INPUTS:
;     IZ:    Either, the name of an ion in CHIANTI format (e.g., o_6
;            for O VI), or, the atomic number of an element. If the
;            latter, then the input ION must also be specified.
;
; OPTIONAL INPUTS:
;     Ion:   The spectroscopic number of the ion. For example, 13 for
;            XIII. Should be specified only if IZ is the atomic
;            number.
;     Ltemp:  A 1D array of Log10 temperatures for which emissivities
;             are required. If not specified then a 3-element array is
;             used that is logTmax+[-0.15,0,0.15].
;     Temp:   A 1D array of temperatures (K).
;     Ldens:  A 1D array of Log10 electron number densities (cm^-3)
;             for which emissivities are required. If not specified,
;             then a 9-element array giving log densities from 8 to 12
;             at 0.5 dex intervals is assumed.
;     Dens:   A 1D array of electron number densities (cm^-3).
;     Dir_Lookup: If set, then the routine looks for the lookup table
;                 in this directory.
;     Pressure: If set, then the density is set to pressure/temp and
;               the emissivity array will be a 1D array. The pressure
;               has units K cm^-3.
;
; KEYWORD PARAMETERS:
;     NO_DE:  If set, then the emissivity is does not include the
;             energy factor.
;     VERBOSE:If set, then information is printed to the IDL window. 
;     QUIET:  No longer used, but kept for backwards compatibility.
;
; OUTPUTS:
;     An IDL structure with the following tags:
;      ION_NAME        STRING    'fe_13'
;      LAMBDA          FLOAT           10749.1
;      LEVEL1          INT              1
;      LVL1_DESC       STRING    '3s2 3p2 3P0'
;      LEVEL2          INT              2
;      LVL2_DESC       STRING    '3s2 3p2 3P1'
;      FLAG            INT              0
;      EM              DOUBLE    Array[3, 9]
;      VERSION         STRING    '9.0.1'
;
;      Flag is set to -1 for transitions with theoretical wavelengths,
;      and 0 otherwise.
;
; RESTRICTIONS:
;     If the densities and/or temperatures are outside
;     the ranges in the lookup table, then the emissivity will be
;     returned as zero for these points.
;
; EXAMPLE:
;     IDL> em=ch_lookup_emiss('fe_13')
;     IDL> em=ch_lookup_emiss('fe_13',ldens=[9,10],temp=[1e6,2e6])
;     IDL> em=ch_lookup_emiss('fe_13',dir_lookup='~/my_lookup_dir')
;
; MODIFICATION HISTORY:
;     Ver.1, 7-Aug-2019, Peter Young
;     Ver.2, 18-Dec-2019, Peter Young
;       Changed processing of lookup filename.
;     Ver.3, 16-Feb-2021, Peter Young
;       Fixed bug when only one density and one temperature.
;     Ver.4, 04-Nov-2021, Peter Young
;       Transitions with negative wavelengths in wgfa file are set to
;       positive wavelength in the output structure, but flag is set
;       to -1. This reproduces behavior of emiss_calc. Disabled /quiet
;       and introduced /verbose. 
;-



IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> em=ch_lookup_emiss( ion_name [, temp=, ltemp=, dens=, ldens=, dir_lookup= '
  print,'                                  pressure=, /no_de, /verbose )'
  print,''
  print,"   ion_name is 'fe_9', for example"
  print,'   Routine requires the CHIANTI lookup tables to exist.'
  return,-1
ENDIF 

;
; If IZ is a string then it is assumed that the ion is being specified
; as, e.g., 'fe_13'.
;
IF datatype(iz) EQ 'STR' THEN BEGIN
  ionname=iz
ENDIF ELSE BEGIN
  zion2name,iz,ion,ionname,diel=diel
  IF keyword_set(diel) THEN ionname=ionname+'.d'
ENDELSE 

;
; Check temperature and density inputs and set default values if
; necessary. 
;
IF n_elements(ltemp) NE 0 THEN temp=10.^ltemp
IF n_elements(ldens) NE 0 THEN dens=10.^ldens
;
IF n_elements(temp) EQ 0 THEN BEGIN
  ltmax=ch_tmax(ionname,/log)
  ltemp=ltmax+[-0.15,0.,0.15]
  temp=10.^ltemp
ENDIF 
;
IF n_elements(pressure) NE 0 THEN BEGIN
  dens=pressure/temp
ENDIF 
;
IF n_elements(dens) EQ 0 THEN BEGIN
  ldens=findgen(9)/2.+8.
  dens=10.^ldens
ENDIF 
nd=n_elements(dens)
nt=n_elements(temp)

IF KEYWORD_SET(verbose) THEN BEGIN
  PRINT,''
  PRINT,'Log_10 temperatures...'
  PRINT,FORMAT='("   ",9f6.2)',alog10(temp)
  PRINT,''
  PRINT,'Log_10 densities...'
  PRINT,FORMAT='("   ",9f6.2)',alog10(dens)
ENDIF


p=ch_lookup_table_interp(ionname, dens, temp, /pad)
IF n_tags(p) EQ 0 THEN BEGIN
  message,' the lookup file was not found. Returning...',/continue
  return,-1
ENDIF
;
IF NOT keyword_set(p.all_levels) THEN BEGIN
  print,"% CH_LOOKUP_EMISS: Warning. The lookup table does not contain populations for all of the ion's levels."
  print,'                   It is recommended that you create a new lookup table.'
ENDIF 

;
; Need to get the list of transitions from the wgfa file.
;
convertname,ionname,iz,ion,diel=diel
zion2filename,iz,ion,fname,diel=diel
read_wgfa_str,fname+'.wgfa',wgfa

IF n_elements(pressure) NE 0 THEN BEGIN
  em_arr=dblarr(nt)
ENDIF ELSE BEGIN
  em_arr=dblarr(nt,nd)
ENDELSE 

str={ion_name: ionname, $
     lambda: 0., $
     level1: 0, $
     lvl1_desc: '', $
     level2: 0, $
     lvl2_desc: '', $
     flag: 0, $
     em: em_arr, $
     version: ch_get_version() }

;
; Filter out transitions with zero A-value (autoionization
; transitions), zero wavelength (2-photon transitions), and for which
; the level indices are higher than those in the lookup file. 
; 
k=where(wgfa.aval NE 0. AND wgfa.wvl NE 0. AND wgfa.lvl1 LE max(p.levels) AND wgfa.lvl2 LE max(p.levels),nk)
wgfa=wgfa[k]
emstr=replicate(str,nk)


;
; Populate emissivity tag.
;
id_t=make_array(nt,value=1.)
IF keyword_set(no_de) THEN BEGIN
  mult_array=(id_t) # (wgfa.aval)
ENDIF ELSE BEGIN
  mult_array=(id_t) # (wgfa.aval*1.986e-8/abs(wgfa.wvl))
ENDELSE
;
IF n_elements(pressure) NE 0 THEN BEGIN
  FOR i=0,nt-1 DO emstr[*].em[i]=reform(p.pop[i,i,wgfa.lvl2-1])*reform(mult_array[i,*])
ENDIF ELSE BEGIN
   CASE 1 OF
      nd EQ 1 AND nt EQ 1: emstr[*].em[0]=reform(p.pop[*,*,wgfa.lvl2-1])*reform(mult_array)
      nd EQ 1 AND nt GT 1: emstr.em=reform(p.pop[*,*,wgfa.lvl2-1])*reform(mult_array)
      ELSE: BEGIN
         IF nt EQ 1 THEN BEGIN
            FOR i=0,nd-1 DO emstr[*].em[0,i]=reform(p.pop[i,*,wgfa.lvl2-1])*reform(mult_array)
         ENDIF ELSE BEGIN 
            FOR i=0,nd-1 DO emstr[*].em[*,i]=reform(p.pop[i,*,wgfa.lvl2-1])*reform(mult_array)
         ENDELSE 
      END 
   ENDCASE 
ENDELSE 

k=where(wgfa.wvl LT 0,nk)
IF nk NE 0 THEN emstr[k].flag=-1
emstr.lambda=abs(wgfa.wvl)
emstr.level1=wgfa.lvl1
emstr.level2=wgfa.lvl2

read_elvlc,fname+'.elvlc',elvlc=elvlc
emstr.lvl1_desc=elvlc.data[emstr.level1-1].full_level
emstr.lvl2_desc=elvlc.data[emstr.level2-1].full_level

return,emstr

END
