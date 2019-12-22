

PRO read_elvlc_str, filename, elvlcstr, ref, time_taken=time_taken

;+
; NAME
;
;     READ_ELVLC_STR
;
; PROJECT
;
;     CHIANTI
;
; EXPLANATION 
;
;     The CHIANTI .elvlc file format was updated in version 8 and this
;     routine reads the new format.
;
;
; INPUTS
;
;     FILENAME  Name of a CHIANTI .elvlc file.
;
; OUTPUTS
;
;     ELVLCSTR  A structure array containing data for each energy
;               level. The tags are (using a Fe X level as an
;               example): 
;
;               ION_NAME        STRING    'fe_10'
;               ION_Z           INT             26
;               ION_N           INT             10
;               ION_ROMAN       STRING    'Fe X'
;               ION_LATEX       STRING    '\ion{Fe}{x}'
;               ION_LATEX_ALT   STRING    '\ion{Fe}{10}'
;               COMMENTS        STRING    Array[10]
;               CHIANTI_VER     STRING    '7.1'
;               TIME_STAMP      STRING    'Tue Dec  4 13:07:39 2012'
;               DATA            STRUCT    -> <Anonymous> Array[825]
;
;               The individual level data are stored in ELVLCSTR.DATA
;               which has the following tags:
;
;               INDEX           INT              8
;               CONF            STRING    '3s2.3p4(3P).3d'
;               CONF_LATEX      STRING    '3s$^2$ 3p$^4$ ($^3$P) 3d'
;               CONF_INDEX      INT              3
;               TERM            STRING    '4F'
;               TERM_LATEX      STRING    '$^4$F'
;               LEVEL           STRING    '4F9/2'
;               LEVEL_LATEX     STRING    '$^4$F$_{9/2}$'
;               FULL_LEVEL      STRING    '3s2.3p4(3P).3d 4F9/2'
;               FULL_LEVEL_LATEX
;                               STRING    '3s$^2$ 3p$^4$ ($^3$P) 3d $^4$F$_{9/2}$'
;               LABEL           STRING    ''
;               MULT            INT              4
;               S               FLOAT           1.50000
;               L               INT              3
;               L_SYM           STRING    'F'
;               J               FLOAT           4.50000
;               J_STR           STRING    '9/2'
;               PARITY          INT              0
;               PARITY_STR      STRING    'e'
;               WEIGHT          FLOAT           10.0000
;               OBS_ENERGY      DOUBLE           417652.00
;               THEORY_ENERGY   DOUBLE          -1.0000000
;               ENERGY          DOUBLE           417652.00                                   
;               ENERGY_UNITS    STRING    'cm^-1'
;               EXTRA_INFO      STRING    ''
;
;     TIME_TAKEN Returns the amount of time (in seconds) that it took
;                to read and process the file.
;
; CALLS
;
;     ZION2SPECTROSCOPIC, CONVERTNAME, CONVERT_CONFIG
;
; HISTORY
;
;     Ver.1, 19-Nov-2012, Peter Young
;     Ver.2, 29-Nov-2012, Peter Young
;        Now sets 'missing' values of the energy arrays to 0 instead
;        of -1 to replicate the original read_elvlc routine. Missing
;        values in the output structure are still set to -1. If data
;        lines contain additional information (beyond the 87
;        characters of required information) then this is stored in
;        the tag EXTRA_INFO. Added ENERGY_UNITS tag.
;     Ver.3, 4-Dec-2012, Peter Young
;        modified output structure.
;     Ver.4, 3-Jan-2012, Peter Young
;        added parity to output structure
;     Ver.5, 14-April-2013, Ken Dere
;		this was originally read_elvlc but has been modified to be read_elvlc_str with just structure as output
;-

t0=systime(1)

IF n_params() LT 1 THEN BEGIN
  print,' '
  print,'   IDL> read_elvlc_str, filename, elvlcstr, ref'
  print,' '
  return 
ENDIF

chck=file_search(filename)
IF chck[0] EQ '' THEN BEGIN
  print,'%READ_ELVLC: file not found. Returning...'
  return
ENDIF

;
; This routine returns the number of lines in the ascii elvlc
; file. This is useful for speeding up the routine (see later). 
;
result=query_ascii(filename,info)
nlines=info.lines

;
; Extract information about the ion that will go into the output structure.
;
fname=file_basename(filename)
bits=strsplit(fname,'.',/extract)
ionname=bits[0]
convertname,ionname,iz,ion
IF iz NE 0 AND ion NE 0 THEN BEGIN
  zion2spectroscopic,iz,ion,ion_roman
  iz=fix(iz)
  ion_roman=strcompress(ion_roman)
  bits=strsplit(ion_roman,' ',/extract)
  ion_latex='\ion{'+bits[0]+'}{'+strlowcase(bits[1])+'}'
  ion_latex_alt='\ion{'+bits[0]+'}{'+trim(ion)+'}'
ENDIF ELSE BEGIN
  print,'%READ_ELVLC:  ion name can not be extracted from filename. Ion name parameters will not be set in output structure.'
  ionname=''
  iz=0
  ion=0
  ion_roman=''
  ion_latex=''
  ion_latex_alt=''
ENDELSE 


spd=['S','P','D','F','G','H','I','K','L','M','N','O','Q','R','T','U','V','W','X','Y']

;
; Read the entire file into a string array
;
file_string=strarr(nlines)
openr,lin,filename,/get_lun
readf,lin,file_string
free_lun,lin

;
; By working out where the -1 is, the number of data lines can be found.
;
k=where(trim(file_string) EQ '-1')
ndata=k[0]
data_string=file_string[0:ndata-1]

;
; The main data will be in the first 87 characters of the string. We
; allow for there being extra characters that may carry additional
; information. The additional information ends up going into the
; 'extra_info' tag.
;
data_string_trunc=data_string
data_string_trunc=strmid(data_string,0,87)
data_string_extra=strmid(data_string,87)

;
; Define structure for reading data, and read data in one go.
;
str={i: 0, conf: '', lbl: '',ss: 0, llstr: '', jj: 0., obs_en: 0d0, th_en: 0d0}
readstr=replicate(str,ndata)

reads,data_string_trunc,format='(i7,a30,a5,i5,a5,f5.1,f15.3,f15.3)',readstr


;
; Read reference lines
; 
nref=nlines-ndata-1
IF nref NE 0 THEN BEGIN 
  ref_arr=file_string[ndata+1:nlines-1]
ENDIF 
ref=ref_arr

vfile=concat_dir(!xuvtop,'VERSION')
chck=file_search(vfile)
vnum=''
IF chck[0] NE '' THEN BEGIN
  openr,lin,vfile,/get_lun
  readf,lin,vnum
  free_lun,lin
  vnum=strtrim(vnum,2)
ENDIF


n=n_elements(readstr)

newstr={ index: 0, $
         conf: '', $
         conf_latex: '', $
         conf_index: 0, $
         term: '', $
         term_latex: '', $
         level: '', $
         level_latex: '', $
         full_level: '', $
         full_level_latex: '', $
         label: '', $
         mult: 0, $
         s: 0.0, $
         l: 0, $
         l_sym: '', $
         j: 0., $
         j_str: '', $
         parity: 0, $
         parity_str: '', $
         weight: 0., $
         obs_energy: 0d0, $
         theory_energy: 0d0, $
         energy: 0d0, $
         energy_units: 'cm^-1', $
         extra_info: ''}
levstr=replicate(newstr,n)

elvlcstr={ ion_name: ionname, $
           ion_z: iz, $
           ion_n: ion, $
           ion_roman: strcompress(ion_roman), $
           ion_latex: ion_latex, $
           ion_latex_alt: ion_latex_alt, $
           comments: ref, $
           chianti_ver: vnum, $
           time_stamp: systime(), $
           data: levstr}
           



elvlcstr.data.index=readstr.i
elvlcstr.data.conf=trim(readstr.conf)
elvlcstr.data.label=trim(readstr.lbl)
elvlcstr.data.mult=readstr.ss
elvlcstr.data.l_sym=trim(readstr.llstr)
elvlcstr.data.j=readstr.jj
elvlcstr.data.obs_energy=readstr.obs_en
elvlcstr.data.theory_energy=readstr.th_en
elvlcstr.data.extra_info=strtrim(data_string_extra,2)

;
; 'energy' contains the best energy for a level.
;
elvlcstr.data.energy=elvlcstr.data.obs_energy
k=where(elvlcstr.data.obs_energy EQ -1,nk)
IF nk GT 0 THEN elvlcstr.data[k].energy=elvlcstr.data[k].theory_energy



;
; Get configuration index. I order them according to their 1st
; occurrence within the level ordering.
;
; I also to the configuration latex conversion here.
;
uniq_conf=get_uniq(elvlcstr.data.conf,count=nconf)
sort_array=intarr(nconf)
FOR i=0,nconf-1 DO begin
  k=where(elvlcstr.data.conf EQ uniq_conf[i])
  sort_array[i]=min(k)
ENDFOR
uniq_conf=uniq_conf[sort(sort_array)]
par_str=['e','o']
FOR i=0,nconf-1 DO begin
  k=where(elvlcstr.data.conf EQ uniq_conf[i])
  elvlcstr.data[k].conf_index=i+1
 ;
  elvlcstr.data[k].conf_latex=strcompress(convert_config(uniq_conf[i],/latex,parity=parity))
  elvlcstr.data[k].parity=parity
  elvlcstr.data[k].parity_str=par_str[parity]
ENDFOR



nspd=n_elements(spd)
FOR i=0,nspd-1 DO BEGIN
  k=where(elvlcstr.data.l_sym EQ spd[i],nk)
  IF nk GT 0 THEN elvlcstr.data[k].l=i
ENDFOR 


elvlcstr.data.s=(float(elvlcstr.data.mult)-1)/2.

elvlcstr.data.weight=2.*elvlcstr.data.j+1.0

elvlcstr.data.term=trim(elvlcstr.data.mult)+trim(elvlcstr.data.l_sym)
elvlcstr.data.term_latex='$^'+trim(elvlcstr.data.mult)+'$'+trim(elvlcstr.data.l_sym)

k=where(elvlcstr.data.j EQ fix(elvlcstr.data.j),nk)
IF nk GT 0 THEN elvlcstr.data.j_str=trim(elvlcstr.data.j)
;
k=where(elvlcstr.data.j NE fix(elvlcstr.data.j),nk)
IF nk GT 0 THEN elvlcstr.data.j_str=trim(elvlcstr.data.j*2)+'/2'


elvlcstr.data.level=elvlcstr.data.term+elvlcstr.data.j_str

elvlcstr.data.level_latex=elvlcstr.data.term_latex+'$_{'+elvlcstr.data.j_str+'}$'
elvlcstr.data.full_level=trim(elvlcstr.data.conf)+' '+trim(elvlcstr.data.level)
elvlcstr.data.full_level_latex=trim(elvlcstr.data.conf_latex)+' '+trim(elvlcstr.data.level_latex)

;
; Create the old-style arrays
;
l1=elvlcstr.data.index
term=elvlcstr.data.full_level
conf=elvlcstr.data.conf_index
ss=elvlcstr.data.mult
ll=elvlcstr.data.l
jj=elvlcstr.data.j
;
; Note: levels with 'missing' energies are set to zero (not -1) in the
; ecm and eryd arrays in order to replicate the behavior of the
; original read_elvlc routine.
;
ecm=dblarr(n_elements(l1))
eryd=ecm
k=where(elvlcstr.data.obs_energy NE -1.,nk)
IF nk NE 0 THEN BEGIN
  ecm[k]=elvlcstr.data[k].obs_energy
  eryd[k]=elvlcstr.data[k].obs_energy/109737.32
ENDIF 
;
ecmth=dblarr(n_elements(l1))
erydth=ecmth
k=where(elvlcstr.data.theory_energy NE -1.,nk)
IF nk NE 0 THEN BEGIN
  ecmth[k]=elvlcstr.data[k].theory_energy
  erydth[k]=elvlcstr.data[k].theory_energy/109737.32
ENDIF 


t1=systime(1)
time_taken=t1-t0

END

