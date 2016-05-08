
PRO read_ff_map, infile, levstr

;+
; NAME
;
;     READ_FF_MAP
;
; EXPLANATION
;
;     Reads a text file containing a mapping from Froese
;     Fischer's level format to CHIANTI level indices. The
;     mapping is stored in an IDL structure.
;
;     The three columns in the file must be separated by white space. 
;
; INPUTS
;
;     INFILE  The name of the map file.
;
; OUTPUTS
;
;     LEVSTR  An IDL structure containing the mapping. The tags are:
;             .conf   FF configuration name
;             .lev    FF level name
;             .chck   String concatenation .conf and .lev
;             .ind    CHIANTI level index
;
; HISTORY
;
;     Ver.1, 4-Oct-2009, Peter Young
;
;     Ver.2, 5-Oct-2009, Peter Young
;        Removed fixed format for reading, so columns just need to be
;        separated by white space now.
;-

str={conf: '', lev: '', chck: '', ind: 0}
levstr=0

openr,lin,infile,/get_lun


str1=''
str2=''
WHILE eof(lin) NE 1 DO BEGIN
  readf,lin,str1
  str1=trim(strcompress(str1))
  bits=str_sep(str1,' ')
  str.conf=bits[0]
  str.lev=bits[1]
  str.ind=fix(bits[2])
  str.chck=trim(bits[0])+' '+trim(bits[1])
  IF n_tags(levstr) EQ 0 THEN levstr=str ELSE levstr=[levstr,str]
ENDWHILE 

free_lun,lin

END
