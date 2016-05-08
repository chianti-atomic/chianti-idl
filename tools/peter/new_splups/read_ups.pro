
PRO read_ups, upsname, upsstr


;+
; NAME
;
;     READ_UPS
;
; PROJECT
;
;     CHIANTI
;
; EXPLANATION
;
;     Reads a file in the CHIANTI UPS format that contains
;     effective collision strengths (upsilons).
;
; INPUTS
;
;     UPSFILE  The name of the CHIANTI UPS format file to read.
;
; OUTPUTS
;
;     UPSSTR   A structure containing the upsilon data. There are two
;              tags called 'info' and 'data' that are each structures.
;
;              UPSSTR.INFO has the following tags:
;
;              ION_NAME        STRING    'mg_5'
;              ION_Z           INT             12
;              ION_N           INT              5
;              ION_ROMAN       STRING    'Mg V'
;              ION_LATEX       STRING    '\ion{Mg}{v}'
;              ION_LATEX_ALT   STRING    '\ion{Mg}{5}'
;              COMMENTS        STRING    Array[2]
;              CHIANTI_VER     STRING    '7.1'
;              TIME_STAMP      STRING    'Fri May 17 18:01:09 2013'
;              FILENAME        STRING    'mg_5.upsdatx'
;              MISSING         FLOAT          -1.00000
;              NTRANS          LONG               666
;
;              UPSSTR.DATA has the following tags:
;
;              LVL1            INT              1
;              LVL2            INT              2
;              DE              FLOAT         0.0155100
;              GF              FLOAT           0.00000
;              LIM             FLOAT          -1.00000
;              NT              INT             41
;              TEMP            FLOAT     Array[41]
;              UPS             FLOAT     Array[41]
;
;
; HISTORY
;
;     Ver.1, 31-May-2013, Peter Young
;     Ver.2, 31-May-2013, Peter Young
;       - the arrays i1,i2,i3 are now long integers.
;-


missing_val=-1


;
; This routine returns the number of lines in the ascii elvlc
; file. This is useful for speeding up the routine (see later). 
;
result=query_ascii(upsname,info)
nlines=info.lines

;
; Read the entire file into a string array
;
file_string=strarr(nlines)
openr,lin,upsname,/get_lun
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
  print,'%READ_UPSDATX: Each transition should have 3 lines of data!!  Returning...'
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
str={lvl1: 0, lvl2: 0, de: 0., gf: 0., lim: 0., nt: 0}
readstr=replicate(str,ntrans)
;
reads,data1,format='(2i7,3e12.3,i5)',readstr

;
; nt_max (maximum number of temperatures for a transition) sets the
; size of the temperature and upsilon arrays 
;
nt_max=max(readstr.nt)

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
; Set the empty columns to be missing data.
;
k=where(temps EQ 0.,nk)
IF nk NE 0 THEN temps[k]=missing_val


;
; Now follow same procedure for the upsilon data
;
data3=strpad(data3,nt_max*12,/after)
ups=fltarr(nt_max,ntrans)
reads,data3,format=format_str,ups
;
k=where(ups EQ 0.,nk)
IF nk NE 0 THEN ups[k]=missing_val



;
; Create the data structure.
;
str2={lvl1: 0, lvl2: 0, de: 0., gf: 0., lim: 0., nt: 0, $
     temp: fltarr(nt_max), ups: fltarr(nt_max)}
datastr=replicate(str2,ntrans)
;
datastr.lvl1=readstr.lvl1
datastr.lvl2=readstr.lvl2
datastr.de=readstr.de
datastr.gf=readstr.gf
datastr.lim=readstr.lim
datastr.nt=readstr.nt
;
datastr.temp=temps
datastr.ups=ups

;
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
basename=file_basename(upsname)
bits=str_sep(basename,'.')
ion_name=bits[0]
convertname,ion_name,ion_z,ion_n
ion_z=fix(ion_z)
zion2spectroscopic,ion_z,ion_n,ion_roman
ion_roman=strcompress(ion_roman)
bits=strsplit(ion_roman,' ',/extract)
ion_latex='\ion{'+bits[0]+'}{'+strlowcase(bits[1])+'}'
ion_latex_alt='\ion{'+bits[0]+'}{'+trim(ion_n)+'}'

;
; Create the info structure
;
info={ion_name: ion_name, ion_z: ion_z, ion_n: ion_n, ion_roman: ion_roman, $
      ion_latex: ion_latex, ion_latex_alt: ion_latex_alt, comments: ref, $
      chianti_ver: vnum, time_stamp: systime(), $
      filename: upsname, missing: float(missing_val), $
      ntrans: ntrans}

;
; Combine the info and data structures into the final output
; structure. 
;
upsstr={info: info, data: datastr}


;
; Tidy up
;
datastr=0 & file_string=0 & readstr=0

END
