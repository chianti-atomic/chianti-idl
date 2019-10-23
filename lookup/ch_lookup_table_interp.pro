
function ch_lookup_table_interp, ion_name, rdens, rtemp, cubic=cubic, $
                                 quiet=quiet,popdir=popdir, $
                                 pad=pad, dir_lookup=dir_lookup

;+
; NAME:
;      CH_LOOKUP_TABLE_INTERP
;
; PURPOSE:
;      This takes the output from ch_read_pop_lookup_table and performs
;      an interpolation to retrieve the level population for the
;      specified densities and temperatures.
;
;      The level population, n_j, is defined such that sum(n_j)=1 for
;      all levels j in the CHIANTI model.
;
; CATEGORY:
;      CHIANTI; populations; interpolation.
;
; CALLING SEQUENCE:
;      Result = CH_LOOKUP_TABLE_INTERP( Ion_Name, Rdens, Rtemp )
;
; INPUTS:
;      Ion_Name: Name of the ion for which the population lookup table
;              is required. Should be given in CHIANTI format, e.g.,
;              "o_6" for O VI. By default the routine will look for
;              the table in $CHIANTI_LOOKUP. See also the FILENAME and
;              POPDIR optional inputs.
;      Rdens:  Densities (units: cm^-3) for which interpolated
;              populations are required. See Programming Notes below
;              for how to specify RDENS. If values are outside of the
;              range within the lookup table, then the populations at
;              these densities will be returned as zeros.
;      Rtemp:  Temperatures (units: K) for which interpolated
;              populations are required. See Programming Notes below
;              for how to specify RTEMP. If values are outside of the
;              range within the lookup table, then the populations at
;              these temperatures will be returned as zeros.
;
; OPTIONAL INPUTS:
;     Dir_Lookup:  Directory containing the lookup table file. The
;                  file name is assumed to be of the form
;                  "pop_lookup_[ION_NAME].txt". This input over-rides
;                  the default directory ($CHIANTI_LOOKUP). 
;     PopDir:  **Obsolete**. Performs same function as
;              DIR_LOOKUP. Retained for backwards compatibility. 
;
; KEYWORD PARAMETERS:
;      PAD:    The lookup table can be created for only a
;              subset of ion's levels (see
;              ch_write_pop_lookup_table.pro). By default,
;              ch_lookup_table_interp will return an array containing
;              only those levels. The array can be expanded to all
;              levels by setting /PAD. Note that the extra levels will
;              all have zero population.
;      QUIET:  If set, then no information messages will be printed to
;              the IDL window.
;      CUBIC:  If set then cubic interpolation is performed with
;              parameter set to -0.5. **Obsolete**
;
; OUTPUTS:
;      An IDL structure with the following tags:
;      POP             DOUBLE    Array[9, 11, 912]
;      LDENS           FLOAT     Array[9]
;      LTEMP           FLOAT     Array[11]
;      LEVELS          INT       Array[912]
;      ALL_LEVELS      INT              1
;      TEMP            FLOAT     Array[11]
;      DENS            FLOAT     Array[9]
;      LDENS_RANGE     FLOAT     Array[2]
;      LTEMP_RANGE     FLOAT     Array[2]
;      CHIANTI_VERSION STRING    '9.0.1'
;
;      LDENS and LTEMP are log10 of DENS and TEMP. LEVELS gives the
;      CHIANTI indices for the levels. POP contains the level
;      populations tabulated for the specified densities and
;      temperatures. LDENS_RANGE and LTEMP_RANGE give the min-max
;      values of the Log-Dens and Log-Temp ranges. 
;
;      ALL_LEVELS indicates whether ch_write_pop_lookup_table
;      attempted to write out all of the levels. Sometimes level
;      populations are zero and these do not get written to the lookup
;      table, so this keyword lets the software know that these levels
;      should be assigned a zero population.
;
;      If a problem is found, then a value of -1 is returned.
;
; CALLS:
;      CH_READ_POP_LOOKUP_TABLE, CH_LOOKUP_FILENAME
;
; EXAMPLE:
;      IDL> rtemp=10.^(findgen(11)/10.+5.0)
;      IDL> rdens=10.^(findgen(9)/2.+8.0)
;      IDL> output=ch_lookup_table_interp('fe_12',rdens,rtemp)
;
; PROGRAMMING NOTES:
;      This routine uses the builtin IDL routine BILINEAR to perform
;      the interpolation. There are two options for how to specify
;      RDENS and RTEMP:
;
;      Option 1 - 2D arrays
;      RDENS and RTEMP are 2D arrays of the same dimensions, such that
;      a pixel (i,j) has a density RDENS(i,j) and temperature
;      RTEMP(i,j). The output population array will have dimensions
;      [NX,NY,NL] where NX is the size of the 1st dimension of
;      RDENS, NY is the size of the 2nd dimension, and NL is the
;      number of levels.
;
;      Option 2 - scalar or 1D arrays, different dimensions
;      RDENS and RTEMP are scalars or 1D arrays of different dimension. The
;      output population array will have dimensions [ND,NT,NL], where
;      ND is the size of RDENS, NT is the size of RTEMP and NL is the
;      number of levels.
;      -------
;      I compared the routines BILINEAR and INTERPOLATE for this
;      routine and decided to use BILINEAR as the cubic interpolation
;      with INTERPOLATE was not always reliable. For BILINEAR
;      it's necessary to use a fine mesh for density (specified
;      by the ch_write_lookup_table routine).
;
; MODIFICATION HISTORY:
;      Ver.1, 2-Oct-2019, Peter Young
;      Ver.2, 18-Dec-2019, Peter Young
;         Now calls to ch_lookup_filename to get the lookup table
;         filename. Removed option to directly specify the filename. 
;-


IF n_params() LT 3 THEN BEGIN
  print,'Use:  IDL> output=ch_lookup_table_interp(ion_name, dens, temp [, /quiet, dir_lookup='
  print,'                                         /pad, /quiet ] )'
  return,-1
ENDIF

;
; Get the lookup table filename
;
filename=ch_lookup_filename(ion_name,dir_lookup=dir_lookup,status=status)
;
IF status EQ 0 THEN BEGIN
  print,'% CH_LOOKUP_TABLE_INTERP: could not find the lookup table '+filename+'.'
  print,'                          Returning...'
  return,-1
ENDIF 

;
; Read the lookup table. 
;
data=ch_read_pop_lookup_table(filename)
IF n_tags(data) EQ 0 THEN return,-1


;
; Do check on the CHIANTI version number.
;
ch_ver=ch_get_version()
IF ch_ver NE data.chianti_version THEN BEGIN
  print,'% CH_LOOKUP_TABLE_INTERP: the lookup table was generated with an old version of CHIANTI. Please '
  print,'                          create a new version. Returning...'
  return,-1
ENDIF 

ldens=data.ldens
dens=data.dens
ltemp=data.ltemp
temp=data.temp
meta=data.meta
popx=data.pop
levels=data.levels
alpha=data.alpha



nd=n_elements(rdens)
nt=n_elements(rtemp)
nl=n_elements(levels)

sd=size(rdens)
st=size(rtemp)


;
; If rtemp and/or rdens are scalars, then to get bilinear to work I
; need to convert them to 1D arrays. This is done by simply creating
; 2-element arrays with duplicate entries. The "swtch" parameters are
; used to restore the original scalar values at the end of the
; routine. 
;
rdens_swtch=0
IF sd[0] EQ 0 THEN BEGIN
  rdens_swtch=1
  rdens=[rdens,rdens]
  nd=2
ENDIF
;
rtemp_swtch=0
IF st[0] EQ 0 THEN BEGIN
  rtemp_swtch=1
  rtemp=[rtemp,rtemp]
  nt=2
ENDIF


;
; Convert RDENS and RTEMP to indices that will be input to BILINEAR.
;
log_rdens=alog10(rdens)
ddens=data.ldens_step
ind_rdens=(log_rdens-ldens[0])/ddens
;
log_rtemp=alog10(rtemp)
dtemp=data.ltemp_step
ind_rtemp=(log_rtemp-ltemp[0])/dtemp



;
; See 'Programming Notes' in header for descriptions of Options 1 and
; 2. 
;
IF sd[0] LE 1 AND sd[0] LE 1 THEN BEGIN
 ;
 ; This is option 2
 ; ----------------
  id_temp=make_array(nt,value=1.0)
  id_dens=make_array(nd,value=1.0)
  ind_rdens=ind_rdens#id_temp
  ind_rtemp=id_dens#ind_rtemp
  ldens_arr=log_rdens#id_temp
  pop=make_array(nd,nt,nl,/double)
ENDIF ELSE BEGIN
 ;
 ; This is option 1
 ; ----------------
  IF sd[0] EQ 2 AND st[0] EQ 2 THEN BEGIN
    ldens_arr=log_rdens
    pop=make_array(sd[1],sd[2],nl,/double)
  ENDIF ELSE BEGIN
    print,'%CH_LOOKUP_TABLE_INTERP: the dimensions of RDENS and RTEMP are not compatible with this routine.'
    print,'                         Please check the Programming Notes in the routine header. Returning...'
    return,-1
  ENDELSE 
ENDELSE 

IF nl EQ 1 THEN BEGIN
  alpha_t=-1./rtemp*alpha[0]
  alpha_t_arr=dblarr(nd,nt)
  FOR i=0,nd-1 DO alpha_t_arr[i,*]=alpha_t
ENDIF ELSE BEGIN
  alpha_t=- ( (1./rtemp) # alpha ) ; * alog10(exp(1.))
  alpha_t_arr=dblarr(nd,nt,nl)
  FOR i=0,nd-1 DO alpha_t_arr[i,*,*]=alpha_t
ENDELSE 

;
; I have to check whether there is only one level, in which case the
; arrays are 2D rather than 3D.
;
; If any of the input density and/or temperatures are outside of the
; ranges of the lookup tables, then the populations are set to
; zero. (Note I'm using 100 as a missing value, since the
; lookup table values can never by 100.)
;
missing_val=100.
IF nl EQ 1 THEN BEGIN
  output=bilinear(popx[*,*],ind_rdens,ind_rtemp,missing=missing_val)
  IF meta EQ 1 THEN BEGIN
    pop=10.^(output-5.0+alpha_t_arr)
  ENDIF ELSE BEGIN
    pop=10.^(output+ldens_arr-25.0+alpha_t_arr)
  ENDELSE
 ;
  k=where(output EQ missing_val,nk)
  IF nk NE 0 THEN BEGIN
    pop[k]=0.
    miss_flag=1b
  ENDIF 
 ;
  IF rdens_swtch EQ 1 THEN pop=pop[0,*]
  IF rtemp_swtch EQ 1 THEN pop=pop[*,0]
 ;  
ENDIF ELSE BEGIN 
  FOR i=0,nl-1 DO BEGIN
   ;
   ; Note: I tried 'interpolate' to do cubic spline interpolation, but
   ; the results were worse than bilinear.
   ;
    output=bilinear(reform(popx[*,*,i]),ind_rdens,ind_rtemp,missing=missing_val)
      ;; output=interpolate(reform(popx[*,*,i]),ind_rdens,ind_rtemp, $
      ;;                    /grid,cubic=0.5)
    IF meta[i] EQ 1 THEN BEGIN
      pop_2d=10.^(output-5.0+alpha_t_arr[*,*,i])
    ENDIF ELSE BEGIN
      pop_2d=10.^(output+ldens_arr-25.0+alpha_t_arr[*,*,i])
    ENDELSE 
   ;
    k=where(output EQ missing_val,nk)
    IF nk NE 0 THEN BEGIN
      pop_2d[k]=0.
      miss_flag=1b
    ENDIF
    pop[*,*,i]=temporary(pop_2d)
  ENDFOR 
 ;
  IF rdens_swtch EQ 1 THEN pop=pop[0,*,*]
  IF rtemp_swtch EQ 1 THEN pop=pop[*,0,*]
 ;
ENDELSE



;
; Reset rdens and rtemp to their original values (if they were scalars
; in the first place).
;
IF rdens_swtch EQ 1 THEN BEGIN
  rdens=rdens[0]
  log_rdens=log_rdens[0]
  nd=1
ENDIF 
IF rtemp_swtch EQ 1 THEN BEGIN
  rtemp=rtemp[0]
  log_rtemp=log_rtemp[0]
  nt=1
ENDIF


;
; Due to fact that populations are approximations of the real
; populations, then I need to renormalize the sum to 1, but only if
; the file contains all the levels. Note that this correction should
; be tiny, though.
;
IF keyword_set(data.all_levels) THEN BEGIN 
  FOR i=0,nd-1 DO BEGIN
    FOR j=0,nt-1 DO BEGIN
      p=reform(pop[i,j,*])
      IF total(p) NE 0. THEN pop[i,j,*]=p/total(p)
    ENDFOR
  ENDFOR
ENDIF 

IF keyword_set(miss_flag) AND NOT keyword_set(quiet) THEN BEGIN
  print,'% CH_LOOKUP_TABLE_INTERP: some of the specified temperature and/or density values are outside of the range'
  print,'                          in the lookup table. Populations for these values are set to zero.'
  status=1
ENDIF

;
; The keyword /PAD creates a population array that includes all levels
; (NLEV_ALL). 
;
; Note that some ions (e.g., o_6) have levels with zero
; population. These aren't stored in the lookup table, but
; using /pad can restore the levels to the population array. 
; 
IF keyword_set(pad) THEN BEGIN
  pop_all=dblarr(nd,nt,data.nlev_all)
  pop_all[*,*,levels-1]=pop
  pop=temporary(pop_all)
  levels=indgen(data.nlev_all)+1
ENDIF 

output={pop: pop, $
        ldens: log_rdens, $
        ltemp: log_rtemp, $
        levels: levels, $
        all_levels: data.all_levels, $
        temp: rtemp, $
        dens: rdens, $
        ldens_range: [min(data.ldens),max(data.ldens)], $
        ltemp_range: [min(data.ltemp),max(data.ltemp)], $
        chianti_version: data.chianti_version }

return,output

END
