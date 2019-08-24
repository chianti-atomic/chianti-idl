

PRO read_scups, splfile, splstr, add_empty=add_empty, verbose=verbose


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
;     Ver. 4, 28 Apr 2014, Giulio Del Zanna 
;         rewritten to speed up by almost a factor of two.
;
; 
;-

IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL>  read_scups, splfile, splstr'
  return 
ENDIF 

t1=systime(/sec)

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

; now read the comments in the end. Start from the end going
; backwards. Assume that there are only two lines with '-1'.
lines_comment=0
comment=0
ref=''

WHILE  comment lt  2    DO  BEGIN  
   if trim(file_string[nlines-1-lines_comment]) EQ '-1' then $
          comment=comment+1 else ref=[ref,trim(file_string[nlines-1-lines_comment])]
   lines_comment=lines_comment+1
endwhile 

ndata=nlines-lines_comment

;k=where(trim(file_string) EQ '-1')
;ndata=k[0]

; data_string=temporary(file_string[0:ndata-1])

;ref=file_string[ndata+1:nlines-1]



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
;

bigarr=fltarr(8, ntrans)
reads, file_string[i1], bigarr

; number of temperatures for each transition:  
nspl=reform(bigarr[5,*]) 


; nt_max (maximum number of temperatures for a transition) sets the
; size of the temperature and upsilon arrays 
;
nt_max=max(nspl)

;
; Create the output data structure.
;
str2={lvl1: 0, lvl2: 0, t_type: 0, de: 0., gf: 0., lim: 0., c_ups: 0., nspl: 0, $
      stemp: fltarr(nt_max)-1, spl: fltarr(nt_max)-1}

IF keyword_set(add_empty) THEN BEGIN
  datastr=replicate(str2,ntrans+1)
ENDIF ELSE BEGIN
  datastr=replicate(str2,ntrans)
ENDELSE 
;
datastr[0:ntrans-1].lvl1=reform(bigarr[0,*])
datastr[0:ntrans-1].lvl2=reform(bigarr[1,*])
datastr[0:ntrans-1].de=reform(bigarr[2,*])
datastr[0:ntrans-1].gf=reform(bigarr[3,*])
datastr[0:ntrans-1].lim=reform(bigarr[4,*])
datastr[0:ntrans-1].nspl=reform(bigarr[5,*]) 
datastr[0:ntrans-1].t_type=reform(bigarr[6,*])
datastr[0:ntrans-1].c_ups=reform(bigarr[7,*])
;

if min(nspl) eq   nt_max then num_temp=nt_max else  begin 
; find the unique numbers of temperatures
   
   num_temp=-1
   
   for ii=1,nt_max do begin  
      indt=where(nspl eq ii, nnt)
      if nnt gt 0 then num_temp=[num_temp,ii]
   endfor 
   num_temp=num_temp[1:*]
   
endelse 

nvt=n_elements(num_temp)

; 
; Fill in each temperature array at a time
;

for ii=0,nvt-1 do begin   
   
   ind=where(nspl eq  num_temp[ii], nn)
   bigarr=fltarr( num_temp[ii], nn)-1
   
   reads,file_string[i2[ind]], bigarr
   datastr[ind].stemp[0:num_temp[ii]-1]=bigarr

   reads,file_string[i3[ind]], bigarr
   datastr[ind].spl[0:num_temp[ii]-1]=bigarr
   
endfor 


; Get the CHIANTI version number
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

if keyword_set(verbose) then begin 
t2=systime(/sec)
print,splfile+', seconds:', t2-t1
endif 



END
