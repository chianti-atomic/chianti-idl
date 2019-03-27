
PRO read_auto, autoname, lvl1, lvl2, auto, autostr=autostr, ref

;+
; NAME:
;      READ_AUTO
;
; PURPOSE:
;      Reads the CHIANTI autoionization files (extension: .auto).
;
; CATEGORY:
;      CHIANTI; autoionization;read.
;
; CALLING SEQUENCE:
;      READ_AUTO, Lvl2, Lvl2, Auto
;
; INPUTS:
;      Autoname:  The name of the file to be read.
;
; OUTPUTS:
;      Lvl1:  A 1D array containing the indices of the final states
;             within the ionized ion.
;      Lvl2:  A 1D array containing the indices of the autoionizing
;             levels in the recombined ion.
;      Auto:  A 1D array containing the autoionization rates (units:
;             s^-1). 
;
; OPTIONAL OUTPUTS:
;      Autostr:  A structure with the following tags:
;                 .lvl1   Same as the LVL1 output.
;                 .lvl2   Same as the LVL2 output.
;                 .auto   Same as the AUTO output.
;
; EXAMPLE:
;      IDL> zion2filename,20,18,fname
;      IDL> autoname=fname+'.auto'
;      IDL> read_auto, autoname, lvl1, lvl2, auto, ref
;         or
;      IDL> read_auto, autoname, autostr=autostr
;
; MODIFICATION HISTORY:
;      Ver.1, 1-Jun-2018, Peter Young
;         Modified version of read_wgfa_str.pro.
;      Ver.2, 6-Jun-2018, Peter Young
;         Modified header; added extra information print-out. 
;-


IF n_params() LT 1 THEN BEGIN
  print,'Use:  read_auto, autoname, autostr=autostr, ref'
  print,' or:  read_auto, autoname, lvl1, lvl2, auto, ref'
  print,''
  print,'    lvl1 - level indices within the ionized ion'
  print,'    lvl2 - level indices within the original ion'
  return
ENDIF 

chck=file_search(autoname)
IF chck EQ '' THEN BEGIN
  print,'%READ_AUTO: file not found. Returning...'
  return
ENDIF

;
; This routine returns the number of lines in the ascii elvlc
; file. This is useful for speeding up the routine (see later). 
;
result=query_ascii(autoname,info)
nlines=info.lines

;
; Read the entire file into a string array
;
file_string=strarr(nlines)
openr,lin,autoname,/get_lun
readf,lin,file_string
free_lun,lin

;
; By working out where the -1 is, the number of data lines can be found.
;
k=where(trim(file_string) EQ '-1')
ndata=k[0]
data_string=file_string[0:ndata-1]
ref=file_string[ndata+1:nlines-2]

;
; The main data will be in the first 55 characters of the string. 
;
data_string_trunc=data_string
data_string_trunc=strmid(data_string,0,55)


;
; Define structure for reading data, and read data in one go.
;
str={lvl1: 0, lvl2: 0, auto: 0. }
autostr=replicate(str,ndata)

reads,data_string_trunc,format='(2i7,e12.2)',autostr

lvl1=autostr.lvl1
lvl2=autostr.lvl2
auto=autostr.auto

END
