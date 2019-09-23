
PRO ch_write_pop_lookup_table, ion_name, $
                               ldens_start=ldens_start, ldens_end=ldens_end, $
                               ldens_step=ldens_step, $
                               ltemp_start=ltemp_start, ltemp_end=ltemp_end, $
                               ltemp_step=ltemp_step, $
                               levels=levels, $
                               outfile=outfile, no_execute=no_execute, $
                               ioneq_file=ioneq_file, $
                               dir_lookup=dir_lookup, $
                               overwrite=overwrite


;+
; NAME:
;      CH_WRITE_POP_LOOKUP_TABLE
;
; PURPOSE:
;      This routine creates a lookup table of scaled level
;      populations. See CHIANTI Technical Report No. 16 for more
;      details. 
;
; CATEGORY:
;      CHIANTI; input/output; level populations.
;
; CALLING SEQUENCE:
;      CH_WRITE_POP_LOOKUP_TABLE, Ion_Name
;
; INPUTS:
;      Ion_Name:  The name of the ion in CHIANTI format for which the
;                 population table is required.
;
; OPTIONAL INPUTS:
;      Ldens_Start:  The starting value of the log density (units:
;                    cm^-3) range for which populations are
;                    required. Default is 8.0.
;      Ldens_End:    The ending value of the log density (units:
;                    cm^-3) range for which populations are
;                    required. Default is 12.0.
;      Ldens_Step:   The step size for log density values. Default is
;                    0.2. 
;      Ltemp_Start:  The starting value of log temperature (units:
;                    K) range for which populations are required. If
;                    not set, then the minimum value in the ioneq file
;                    is used.
;      Ltemp_End:    The ending value of log temperature (units:
;                    K) range for which populations are required. If
;                    not set, then the maximum value in the ioneq file
;                    is used.
;      Ltemp_Step:   The step size for log density values. Default is
;                    0.05.
;      Outfile:      The name of the output file. If not specified,
;                    then it is of the form 'pop_lookup_IONNAME.txt'.
;      Levels:       An integer, or integer array, specifying the
;                    CHIANTI level indices of the levels to be written
;                    to the file. If not set, then all levels are
;                    written.
;      Ioneq_File:   The name of an ionization equilibrium file. The
;                    default is to use the CHIANTI file (!ioneq_file).
;      Dir_Lookup:   The name of a directory to which the lookup file
;                    will be sent. If OUTFILE is specified then
;                    DIR_LOOKUP is ignored.
;
; KEYWORD PARAMETERS:
;      NO_EXECUTE:   If set, then the routine will not calculate or
;                    write the lookup table, but it will populate the
;                    temperature and density range parameters. This is
;                    useful for the user to find out what the default
;                    temperature and density ranges are without having
;                    to do a calculation.
;      OVERWRITE:    By default the routine checks if the lookup table
;                    already exists and whether it's consistent with
;                    the one that would be written. If yes, then the
;                    lookup table is not written. You can force it to
;                    be written with /OVERWRITE.
; 
; OUTPUTS:
;      Creates a text file in the working directory with the name
;      'pop_lookup_IONNAME.txt' where IONNAME is the ion name in
;      CHIANTI format. The name can be overwritten with the keyword
;      OUTFILE. The format of the file is as follows:
;
;      Line 1: CHIANTI version number.
;      Line 2: States whether all levels were calculated (0/1).
;      Line 3: Specifies sizes of arrays.
;      Line 4: Log temperatures.
;      Line 5: Log densities.
;      Line 6 onwards: modified populations as function of density,
;              tabulated for each temperature.
;
;      See CHIANTI Technical Report No. 16 for information on the
;      modified populations.
;
; EXAMPLE:
;      IDL> ch_write_pop_lookup_table, 'c_3'
;      IDL> ch_write_pop_lookup_table, 'c_3', dir_lookup='~/mydata'
;      IDL> ch_write_pop_lookup_table, 'c_3', outfile='my_lookup_file.txt'
;
; MODIFICATION HISTORY:
;      Ver.1, 2-Oct-2019, Peter Young
;-


IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> ch_write_pop_lookup_table, ion_name [, ldens_start=, ldens_end= '
  print,'                       ldens_step=, ltemp_start=, ltemp_end=, ltemp_step='
  print,'                       levels=, outfile=, /no_execute, ioneq_file=, dir_lookup= ]'
  return
ENDIF 

IF n_elements(ioneq_file) EQ 0 THEN ioneq_file=!ioneq_file

IF n_elements(ltemp_step) EQ 0 THEN ltemp_step=0.05
;
; If the ltemp inputs are not specified, then I determine the
; temperature range from the ioneq file. I include temperatures for
; which the ion fraction is > 10^-8 of the max ion fraction.
;
IF n_elements(ltemp_start) EQ 0 OR n_elements(ltemp_end) EQ 0 THEN BEGIN
  read_ioneq,ioneq_file,tt,ii,ref
  convertname,ion_name,iz,ion,diel=diel
  ioneq=ii[*,iz-1,ion-1+diel]
  k=where(ioneq GE max(ioneq)*1e-8,nk)
  ltemp_start=tt[k[0]]
  ltemp_end=tt[k[nk-1]]
ENDIF 
lt0=ltemp_start
lt1=ltemp_end
dt=ltemp_step
;
ntemp=round((lt1-lt0)/dt) +1
ltemp=findgen(ntemp)*dt+lt0
temp=10.^ltemp

;
; Set default density parameters.
;
IF n_elements(ldens_step) EQ 0 THEN ldens_step=0.2
IF n_elements(ldens_start) EQ 0 OR n_elements(ldens_end) EQ 0 THEN BEGIN
  ldens_start=8.0
  ldens_end=12.0
ENDIF
ld0=ldens_start
ld1=ldens_end
dd=ldens_step
;
ndens=round((ld1-ld0)/dd) +1
ldens=findgen(ndens)*dd+ld0
dens=10.^ldens

IF keyword_set(no_execute) THEN return

;
; all_levels is written to the output file later
;
IF n_elements(levels) NE 0 THEN BEGIN
  all_levels=0
  txt='(no) '
ENDIF ELSE BEGIN
  all_levels=1
  txt='(yes)'
ENDELSE 



;
; Here I create the filename of the lookup table. The priority of the
; inputs is:
;   - use FILENAME, if specified.
;   - create DIR_LOOKUP/pop_lookup_[ion_name].txt, if dir_lookup
;     specified.
;   - create $CHIANTI_LOOKUP/pop_lookup_[ion_name].txt
;   - otherwise create pop_lookup_[ion_name].txt in current working
;     directory 
;
IF n_elements(outfile) EQ 0 THEN BEGIN
  outfile='pop_lookup_'+ion_name+'.txt'
  IF n_elements(dir_lookup) NE 0 THEN BEGIN
    chck=file_info(dir_lookup)
    IF chck.directory EQ 0 THEN file_mkdir,dir_lookup
    outfile=concat_dir(dir_lookup,outfile)
  ENDIF ELSE BEGIN
    chck=getenv('CHIANTI_LOOKUP')
    IF chck NE '' THEN outfile=concat_dir(chck,outfile)
  ENDELSE 
ENDIF

;
; Check if the lookup table already exists and is consistent with the
; input parameters (i.e., ltemp and ldens arrays are the same). Also
; check if CHIANTI version number is the same. The /overwrite keyword
; can be used to bypass this check.
;
IF NOT keyword_set(overwrite) THEN BEGIN 
  chck=file_info(outfile)
  IF chck.exists EQ 1 THEN BEGIN
    swtch=0
    pop_table=ch_read_pop_lookup_table(outfile)
    IF pop_table.chianti_version EQ ch_get_version() THEN swtch=swtch+1
    IF n_elements(ldens) EQ n_elements(pop_table.ldens) THEN swtch=swtch+1
    IF n_elements(ltemp) EQ n_elements(pop_table.ltemp) THEN swtch=swtch+1
    IF max(abs(ldens-pop_table.ldens)/pop_table.ldens) LE 0.001 THEN swtch=swtch+1
    IF max(abs(ltemp-pop_table.ltemp)/pop_table.ldens) LE 0.001 THEN swtch=swtch+1
   ;
    IF swtch EQ 5 THEN BEGIN
      print,'% CH_WRITE_POP_LOOKUP_TABLE: lookup table already exists, and no need to update. An over-write'
      print,'                             can be forced with the /OVERWRITE keyword.'
      return
    ENDIF 
  ENDIF 
ENDIF 


;
; Identify metastable levels.
;
metastable_levels,ion_name,meta,/quiet

;
; Load atomic data for input to pop_solver
;
input=ch_setup_ion(ion_name)

;
; Get level populations.
;
pop_solver,input, temp, dens ,pop

s=size(pop,/dim)
nlev_all=s[2]
nlev=nlev_all

;
; The following computes the population scaling factor "alpha" (see
; CHIANTI Technical Report #16). Note that the definition depends on
; the requested temperature and density ranges.
;
; alpha_t is alpha/T.
;
; If the populations of a level are zero, then alpha is not
; defined. These levels don't get written to the output file, though. 
;
alpha=fltarr(nlev)
tfac=(temp[ntemp-1]*temp[0])/(temp[ntemp-1]-temp[0])
FOR i=0,nlev-1 DO BEGIN
  mp0=mean(pop[0,*,i])
  mp1=mean(pop[ntemp-1,*,i])
  IF mp0 NE 0. AND mp1 NE 0. THEN alpha[i]=(alog10(mp1)-alog10(mp0))*tfac
ENDFOR
id_t=make_array(ntemp,value=1.)
id_lev=make_array(nlev,value=1.)
alpha_t=(id_t # alpha) / (temp#id_lev)


IF n_elements(levels) EQ 0 THEN levels=findgen(nlev)+1 ELSE nlev=n_elements(levels)


;
; This fixes a problem I found with fe_11d. The populations of 3
; levels were zero because they were not being excited by any
; mechanism. This is OK for the normal CHIANTI software, but because I
; take the log here, then I get infinities. The lines below remove the
; zero-population levels from the table.
;
; There's a potential problem here in that the population may be
; zero for the lowest temperature, but non-zero for other
; temperatures. The check is only done on the lowest temperature. 
;
pchck=reform(pop[0,0,levels-1])
k=where(pchck GT 0.,nk)
IF nk NE nlev THEN BEGIN
  levels=levels[k]
  print,'% CH_WRITE_POP_LOOKUP_TABLE: Warning - '+trim(nlev-nk)+' levels of this ion have zero population. '
  print,'                             They are not written to the lookup table.'
  nlev=nk
ENDIF 


openw,lout,outfile,/get_lun

;
; Prints lines 1-5 to the output file.
;
ch_ver=ch_get_version()
ch_ver=strpad(ch_ver,10,fill=' ')
printf,lout,format='("CHIANTI version: ",a10)',ch_ver
printf,lout,format='("ION: ",a5,"   NO. of LEVELS: ",i7)',ion_name,nlev_all
printf,lout,format='("ALL_LEVELS? ",i3,"   ",a5)',all_levels,txt
printf,lout,format='("SIZES:   ",2i4,i7,"   [temp, dens, lev]")',ntemp,ndens,nlev
printf,lout,format='("TEMP STEP: ",f10.4,"   DENS STEP: ",f10.4)',ltemp_step,ldens_step
printf,lout,format='("LOG TEMP:",'+trim(ntemp)+'e11.3)',ltemp
printf,lout,format='("LOG DENS:",'+trim(ndens)+'e11.3)',ldens

;
; Now write out "alpha" for the levels
;
printf,lout,'POPULATION SCALE FACTORS (ALPHA):'
printf,lout,format='(10e12.3)',alpha[levels-1]


;
; This is the format for the main data. The 2nd format is in case the
; pop parameter goes below -100.
;
format1='(i4,i7,i3,'+trim(ndens)+'f11.6)'
format2='(i4,i7,i3,'+trim(ndens)+'f11.5)'

;
; Now print the populations.
;    i4  - temperature index
;    i7  - CHIANTI level index
;    i3  - metastable flag (1-yes, 0-no)
;    Ne11.3 - scaled populations
;
; Note that scaling varies on whether the level is metastable or not.
;
printf,lout,'POPULATION LOOKUP DATA:'
FOR i=0,ntemp-1 DO BEGIN
  FOR j=0,nlev-1 DO BEGIN
    mt=meta[levels[j]-1]
    lpop=alog10(pop[i,*,levels[j]-1])+alpha_t[i,levels[j]-1]
    IF mt EQ 1 THEN pop_print=lpop+5. $
    ELSE pop_print=lpop-alog10(dens)+25.
    IF min(lpop) LE -100. THEN format=format2 ELSE format=format1
    printf,lout,format=format,i,levels[j],mt,pop_print
  ENDFOR 
ENDFOR


free_lun,lout

print,'% CH_WRITE_POP_LOOKUP_TABLE: the lookup table has been written to '+outfile+'.'

END
