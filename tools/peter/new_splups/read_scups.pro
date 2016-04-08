

PRO read_scups, splfile, splstr, add_empty=add_empty


;+
; NAME
;
;     READ_SCUPS
;
; PROJECT
;
;     CHIANTI
;
; EXPLANATION
;
;     Reads a file in the CHIANTI SCUPS format.
;
; INPUTS
;
;     SPLFILE  The name of the SCUPS file to read.
;
; KEYWORDS
;
;     ADD_EMPTY  If set, then an additional line of data is added to
;                SPLSTR.DATA. This extra line is empty, and is
;                intended for use by WRITE_SCUPS_GUI.
;
; OUTPUTS
;
;     SPLSTR   A structure containing the spline data. There are two
;              tags called 'info' and 'data' that are each structures.
;
;              UPSSTR.INFO has the following tags:
;
;              ION_NAME        STRING    'ca_2'
;              ION_Z           INT             20
;              ION_N           INT              2
;              ION_ROMAN       STRING    'Ca II'
;              ION_LATEX       STRING    '\ion{Ca}{ii}'
;              ION_LATEX_ALT   STRING    '\ion{Ca}{2}'
;              COMMENTS        STRING    Array[12]
;              CHIANTI_VER     STRING    '7.1'
;              TIME_STAMP      STRING    'Fri May 24 16:54:46 2013'
;              FILENAME        STRING    'ca_2.scups_auto'
;              MISSING         FLOAT          -1.00000
;              NTRANS          LONG               766
;
;              UPSSTR.DATA has the following tags:
;
;              LVL1            INT              1
;              LVL2            INT              2
;              T_TYPE          INT              2
;              DE              FLOAT          0.124400
;              GF              FLOAT           0.00000
;              LIM             FLOAT          -1.00000
;              C_UPS           FLOAT          0.560000
;              NSPL            INT             13
;              STEMP           FLOAT     Array[18]
;              SPL             FLOAT     Array[18]
;
; HISTORY
;
;     Ver. 1, 24-May-2013, Peter Young
;     Ver. 2, 28-May-2013, Peter Young
;         added /ADD_EMPTY keyword.
;     Ver. 3, 31-May-2013, Peter Young
;         the arrays i1,i2,i3 are now long integers.
;-

IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL>  read_scups, splfile, splstr'
  return 
ENDIF 


missing_val=-1


;
; This routine returns the number of lines in the ascii elvlc
; file. This is useful for speeding up the routine (see later). 
;
result=query_ascii(splfile,info)
nlines=info.lines

;
; Read the entire file into a string array
;
file_string=strarr(nlines)
openr,lin,splfile,/get_lun
readf,lin,file_string
free_lun,lin

; By working out where the -1 is, the number of data lines can be found.
;
k=where(trim(file_string) EQ '-1')
ndata=k[0]
data_string=file_string[0:ndata-1]
ref=file_string[ndata+1:nlines-1]

ntrans=ndata/3
IF 3*ntrans NE ndata THEN BEGIN
  print,'%READ_SCUPS: Each transition should have 3 lines of data!!  Returning...'
  return
ENDIF 

;
; Each transition has 3 lines of data, so put these into separate
; arrays (data1, data2, data3).
;
i1=lindgen(ntrans)*3
i2=lindgen(ntrans)*3+1
i3=lindgen(ntrans)*3+2

data1=data_string(i1)
data2=data_string(i2)
data3=data_string(i3)

;
; Read data1 into a structure.
;
str={lvl1: 0, lvl2: 0, de: 0., gf: 0., lim: 0., nspl: 0, t_type: 0, c_ups: 0. }
readstr=replicate(str,ntrans)
;
reads,data1,format='(2i7,3e12.3,2i5,e12.3)',readstr


; nt_max (maximum number of temperatures for a transition) sets the
; size of the temperature and upsilon arrays 
;
nt_max=max(readstr.nspl)

;
; Pad out data2 so that each line has the correct size for reading.
;
data2=strpad(data2,nt_max*12,/after)

;
; Read data2 into a temperature array. Note that empty columns will
; end up as zeros.
;
temps=fltarr(nt_max,ntrans)
format_str='('+trim(nt_max)+'e12.3)'
reads,data2,format=format_str,temps


;
; Set the empty columns to be missing data. This method is quicker
; than going through the transitions one-by-one.
;
FOR i=nt_max-1,1,-1 DO BEGIN
  j=where(readstr.nspl EQ i,nj)
  IF nj NE 0 THEN BEGIN
    temps[i:*,j[0]]=missing_val
  ENDIF 
ENDFOR 


;
; Now follow same procedure for the upsilon data
;
data3=strpad(data3,nt_max*12,/after)
ups=fltarr(nt_max,ntrans)
reads,data3,format=format_str,ups
;
k=where(temps EQ missing_val,nk)
IF nk NE 0 THEN ups[k]=missing_val



;
; Create the data structure.
;
str2={lvl1: 0, lvl2: 0, t_type: 0, de: 0., gf: 0., lim: 0., c_ups: 0., nspl: 0, $
     stemp: fltarr(nt_max), spl: fltarr(nt_max)}
IF keyword_set(add_empty) THEN BEGIN
  datastr=replicate(str2,ntrans+1)
ENDIF ELSE BEGIN
  datastr=replicate(str2,ntrans)
ENDELSE 
;
datastr[0:ntrans-1].lvl1=readstr.lvl1
datastr[0:ntrans-1].lvl2=readstr.lvl2
datastr[0:ntrans-1].de=readstr.de
datastr[0:ntrans-1].gf=readstr.gf
datastr[0:ntrans-1].lim=readstr.lim
datastr[0:ntrans-1].nspl=readstr.nspl
datastr[0:ntrans-1].c_ups=readstr.c_ups
datastr[0:ntrans-1].t_type=readstr.t_type
;
datastr[0:ntrans-1].stemp=temps
datastr[0:ntrans-1].spl=ups



; Get CHIANTI version number
;
vfile=concat_dir(!xuvtop,'VERSION')
chck=file_search(vfile)
vnum=''
IF chck[0] NE '' THEN BEGIN
  openr,lin,vfile,/get_lun
  readf,lin,vnum
  free_lun,lin
  vnum=strtrim(vnum,2)
ENDIF


;
; Create the various ion identifiers from the filename.
;
basename=file_basename(splfile)
bits=str_sep(basename,'.')
ion_name=bits[0]
convertname,ion_name,ion_z,ion_n
ion_z=fix(ion_z)
zion2spectroscopic,ion_z,ion_n,ion_roman
ion_roman=strcompress(ion_roman)
bits=strsplit(ion_roman,' ',/extract)
ion_latex='\ion{'+bits[0]+'}{'+strlowcase(bits[1])+'}'
ion_latex_alt='\ion{'+bits[0]+'}{'+trim(ion_n)+'}'


IF keyword_set(add_empty) THEN ntrans=ntrans+1

;
; Create the info structure
;
info={ion_name: ion_name, ion_z: ion_z, ion_n: ion_n, ion_roman: ion_roman, $
      ion_latex: ion_latex, ion_latex_alt: ion_latex_alt, comments: ref, $
      chianti_ver: vnum, time_stamp: systime(), $
      filename: splfile, missing: float(missing_val), $
      ntrans: ntrans}

;
; Combine the info and data structures into the final output
; structure. 
;
splstr={info: info, data: datastr}




END
