

FUNCTION ch_read_pop_lookup_table, filename

;+
; NAME:
;      CH_READ_POP_LOOKUP_TABLE
;
; PURPOSE:
;      Read a CHIANTI level population lookup table that has been
;      created with ch_write_pop_lookup_table.pro.
;
; CATEGORY:
;      CHIANTI; input/output; level populations.
;
; CALLING SEQUENCE:
;      Result = CH_READ_POP_LOOKUP_TABLE( Filename )
;
; INPUTS:
;      Filename:  Name of the population lookup table.
;
; OUTPUTS:
;      A structure with the following tags
;       pop  3D array (ndens*ntemp*nlev) containing scaled level
;            populations 
;       dens Electron number densities (cm^-3) at which populations
;            tabulated.
;       temp Temperatures (K) at which populations tabulated.
;       ldens Log densities.
;       ltemp Log temperatures.
;       meta Flag (0/1) indicating if level is metastable (1).
;       levels CHIANTI level indices of levels in pop.
;       all_levels Flag (0/1) indicating if all of the levels were
;                  calculated.
;       alpha Scaling parameter used for deriving populations. 
;       chianti_version  String containing CHIANTI version number.
;
; EXAMPLE:
;      IDL> pop=ch_read_pop_lookup_table('pop_lookup_si_12.txt')
;
; MODIFICATION HISTORY:
;      Ver.1, 2-Oct-2019, Peter Young
;-


IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> output=ch_read_pop_lookup_table(filename)'
  return,-1
ENDIF 

chck=file_search(filename,count=count)
IF count EQ 0 THEN BEGIN
  print,'%CH_READ_POP_LOOKUP_TABLE: input file not found. Returning...'
  return,-1
ENDIF ELSE BEGIN 
  openr,lin,filename,/get_lun
ENDELSE

ch_ver=''
ionname=''
nlev_all=0l & nt=0l & nd=0l & nlev=0l   ; make these long integers
all_levels=-1
readf,lin,format='(17x,a10)',ch_ver
readf,lin,format='(5x,a5,18x,i7)',ionname,nlev_all
readf,lin,format='(12x,a3)',all_levels
readf,lin,format='(9x,2i4,i7)',nt,nd,nlev
readf,lin,format='(11x,f10.0,14x,f10.0)',ltemp_step,ldens_step
ltemp=fltarr(nt)
ldens=fltarr(nd)
readf,lin,format='(9x,'+trim(nt)+'e11.0)',ltemp
readf,lin,format='(9x,'+trim(nd)+'e11.0)',ldens
alpha=fltarr(nlev)
str1=''
readf,lin,str1
readf,lin,format='(10e12.0)',alpha

str={ t_index: 0, level: 0, meta: 0b, pop: dblarr(nd) }
data=replicate(str,nt*nlev)
readf,lin,str1
readf,lin,data

free_lun,lin


;
; Extract the lookup table into a 3D array (pop)
;
pop=dblarr(nd,nt,nlev)
ind=findgen(nt)*nlev
FOR i=0,nlev-1 DO BEGIN
  pop[*,*,i]=data[ind+i].pop
ENDFOR

;
; Construct the output structure.
;
dens=10.^ldens
temp=10.^ltemp
meta=data[0:nlev-1].meta
levels=data[0:nlev-1].level
;
output={ion: trim(ionname), $
        pop: pop, $
        dens: dens, $
        temp: temp, $
        meta: meta, $
        levels: levels, $
        nlev_all: nlev_all, $
        all_levels: all_levels, $
        ldens: ldens, $
        ldens_step: ldens_step, $
        ltemp: ltemp, $
        ltemp_step: ltemp_step, $
        alpha: alpha, $
        chianti_version: trim(ch_ver)}

return,output

END

