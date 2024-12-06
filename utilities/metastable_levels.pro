
PRO metastable_levels, ionname, meta, cutoff=cutoff, path=path, $
                       quiet=quiet, density=density

;+
; NAME:
;    METASTABLE_LEVELS
;
; PURPOSE:
;    This routine calculates which of an ion's levels are expected
;    to be metastables, based on data from the elvlc and wgfa files
;    (i.e., without having to solve the level balance equations).
;    The calculation comes from the coronal approximation, whereby
;    the excited level is assumed to be populated only from the
;    ground (with population 1) by a transition with collision
;    strength 1. 
;
; CATEGORY:
;    CHIANTI; utility.
;
; CALLING SEQUENCE:
;    METASTABLE_LEVELS, Ionname, Meta
;
; INPUTS:
;    Ionname:   Ion name in CHIANTI format. E.g., 'fe_13' for Fe XIII.
;
; OPTIONAL INPUTS:
;    Cutoff:   Cutoff for determining which levels are metastable. The
;              default is 5e4. Decreasing the cutoff will give more
;              metastable levels.
;    Path:     By default the routine looks in the user's
;              CHIANTI distribution for the .wgfa file. Setting PATH
;              allows you to directly specify a directory containing
;              the .wgfa file.
;    Density:  Density (cm^-3) for which the metastable calculation is
;              performed. If an array is specified, then the maximum
;              density in the array is used.
;	
; KEYWORD PARAMETERS:
;    QUIET:    If set, then no information is printed to the screen.
;
; OUTPUTS:
;    Meta:     A byte array with N elements, where N is the number of
;              levels in the ion model (as judged from the .wgfa file).
;              Metastable levels are identified in this array with the
;              value 1.
;
; CALLS:
;    READ_WGFA_STR, READ_ELVLC, CONVERTNAME, ZION2FILENAME, CH_TMAX
;
; EXAMPLE:
;    IDL> metastable_levels, 'fe_13', meta
;
; MODIFICATION HISTORY:
;    Ver.1, 4-Feb-2009, Peter Young
;    Ver.2, 13-Feb-2009, Peter Young
;      changed input to ionname; added print out of levels; added
;      /quiet keyword.
;    Ver.3, 26-Jul-2019, Peter Young
;      fixed error when a dielectronic ion is specified; updated
;      header format.
;    Ver.4, 2-Aug-2019, Peter Young
;      moved read_elvlc inside if statement; now check whether the sum
;      of A-value + autoionization rate is greater than cutoff.
;    Ver.5, 06-Nov-2023, Peter Young
;      Major revision to routine so that it more accurately works
;      out potential metastables using the coronal approximation.
;    Ver.6, 04-Nov-2024, Peter Young
;      Modified call to ch_tmax for CHIANTI 11; updated header.
;-


IF n_params() LT 1 THEN BEGIN
  print,'Use: IDL> metastable_levels, ionname, meta [, density=, cutoff=, path=, /quiet]'
  return
ENDIF

;
; If an array of densities is given, then the maximum density is used.
;
IF n_elements(density) EQ 0 THEN dens=1e10 ELSE dens=max(density)

;
; This defines the cutoff, above which a level is considered to be a
; metastable. Decreasing cutoff will lead to more metastables.
; I chose the value below by comparing the actual metastables with the
; estimated metastables for densities of 1e8, 1e10 and 1e13. All
; metastables with populations with populations within a factor 100
; of the most populated level were correctly found. 
;
IF n_elements(cutoff) EQ 0 THEN cutoff=5e4

convertname,ionname,iz,ion,diel=diel


IF n_elements(path) EQ 0 THEN BEGIN 
  zion2filename,iz,ion,filename,diel=diel
  wgfaname=filename+'.wgfa'
  elvlcname=filename+'.elvlc'
ENDIF ELSE BEGIN
  wgfaname=concat_dir(path,ionname+'.wgfa')
  elvlcname=concat_dir(path,ionname+'.elvlc')
ENDELSE 

chck=file_info(wgfaname)
IF chck.exists EQ 0 THEN BEGIN
  meta=-1
  message,/info,/cont,'The .wgfa file can not be found. Returning...'
  return
ENDIF
read_wgfa_str,wgfaname,str,ref


n=max(str.lvl2)
meta=bytarr(n)

;
; For the dielectronic models there are no metastables.
;
IF diel THEN BEGIN
  meta[0]=1b
  return
ENDIF 


read_elvlc,elvlcname,l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,ref,elvlc=elvlcstr
e_ev=elvlcstr.data.energy/8065.54
wgt=elvlcstr.data.weight

;
; The electron excitation rates are temperature sensitive, so need to
; specify a temperature. Using Tmax calculated from zero-density ion
; balance.
;
kt=ch_tmax(ionname,advanced_model=0,/quiet)*8.617e-5
exp_de_kt=exp(-e_ev/kt)
kt_sqrt=sqrt(13.606/kt)

;
; PRY, 2-Aug-2019: I now sum aval and auto to fix a problem for fe_11d
; whereby levels 11-20 had no A-values, only auto rates.
;
FOR i=1,n DO BEGIN
  k=where(str.lvl2 EQ i)
  IF k[0] EQ -1 THEN BEGIN
    meta[i-1]=1b
  ENDIF ELSE BEGIN
    amax=max(str[k].aval+str[k].auto)
    IF dens/amax*exp_de_kt[i-1]*kt_sqrt/wgt[i-1] GE cutoff THEN meta[i-1]=1b
  ENDELSE
ENDFOR


IF NOT keyword_set(quiet) THEN BEGIN 
  k=where(meta EQ 1,nk)
  print,'Metastable levels are: '
  FOR i=0,nk-1 DO BEGIN
    print,format='(i3,a25)',l1[k[i]],term[k[i]]
  ENDFOR
ENDIF 

END
