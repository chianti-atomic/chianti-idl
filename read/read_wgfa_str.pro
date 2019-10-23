
PRO read_wgfa_str, wgfaname, wgfastr, ref

;+
; NAME
;
;    READ_WGFA_STR
;
; PROJECT
;
;    CHIANTI
;
; EXPLANATION
;
;    Reads the CHIANTI .wgfa file into a structure. For the most part
;    CHIANTI wgfa files contain radiative decay rates. I.e., rates for
;    transitions that yield photons of a definite wavelength. However
;    for some ions the .wgfa file is also used to contain
;    autoionization rates and/or two photon transitions. These are
;    needed for accurately modelling the level populations within
;    ions, but do not yield photons of a definite wavelength
;    (autoionization yields no photons; two photon transitions yield a
;    continuum which is separately modelled in CHIANTI).
;
;    These "radiationless" transitions are denoted in the .wgfa file
;    with a zero wavelength. 
;
;    In the output structure all radiative decay rates are stored in
;    the .AVAL tag, while autoionization and two photon rates are
;    stored in the .AUTO tag.
;
; INPUTS
;
;    WGFANAME  The name of the file to be read.
;
; OUTPUTS
;
;    STR   A structure with the following tags:
;          .lvl1  The index of the lower level of the transition. 
;          .lvl2  The index of the upper level of the transition.
;          .wvl   The wavelength of the transition (angstroms).
;          .gf    The weighted oscillator strength.
;          .aval  The radiative decay rate (s^-1).
;          .auto  The autoionization rate (s^-1).
;
;    REF   A string array containing the file references.
;
; HISTORY
;
;    Ver.1, 13-Feb-2009, Peter Young
;    Ver.2, 24-Apr-2013, Peter Young
;        This routine was very slow for large ions so I've
;        completely changed the way the data are read.
;-

chck=file_search(wgfaname)
IF chck EQ '' THEN BEGIN
  print,'%READ_WGFA_STR: file not found. Returning...'
  return
ENDIF

;
; This routine returns the number of lines in the ascii elvlc
; file. This is useful for speeding up the routine (see later). 
;
result=query_ascii(wgfaname,info)
nlines=info.lines

;
; Read the entire file into a string array
;
file_string=strarr(nlines)
openr,lin,wgfaname,/get_lun
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
str={lvl1: 0, lvl2: 0, wvl: 0d, gf: 0d, aval: 0d }
readstr=replicate(str,ndata)

reads,data_string_trunc,format='(2i5,f15.0,2e15.3)',readstr

;
; This defines the output structure which contains an extra tag for
; the autoionization rate.
;
str2={lvl1: 0, lvl2: 0, wvl: 0d, gf: 0d, aval: 0d, auto: 0d}
wgfastr=replicate(str2,ndata)

wgfastr.lvl1=readstr.lvl1
wgfastr.lvl2=readstr.lvl2
wgfastr.wvl=readstr.wvl
wgfastr.gf=readstr.gf

;
; Only load A-values into wgfastr for those transitions with
; non-negative wavelengths.
;
k=where(readstr.wvl NE 0.)
wgfastr[k].aval=readstr[k].aval

;
; The following loads the autoionization data into the output
; structure. The 'flag' array is used to flag lines of data that need
; to be removed from the output structure. The removed lines are those
; for which both a A-value and a autoionization rate are available. 
;
flag=bytarr(ndata)+1b
k=where(readstr.wvl EQ 0.,nk)
IF nk NE 0 THEN BEGIN
  FOR i=0,nk-1 DO BEGIN
    l1=readstr[k[i]].lvl1
    l2=readstr[k[i]].lvl2
   ;
    ii=where(readstr.lvl1 EQ l1 AND readstr.lvl2 EQ l2,nii)
    wgfastr[ii].auto=readstr[k[i]].aval
   ;
    IF nii EQ 2 THEN flag[k[i]]=0b
  ENDFOR 
ENDIF 

k=where(flag EQ 1b)
wgfastr=wgfastr[k]

END
