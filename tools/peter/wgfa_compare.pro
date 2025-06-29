
PRO wgfa_compare, file1, file2, list_file1=list_file1, list_file2=list_file2, $
                  limit=limit, lev_cutoff=lev_cutoff, aval_cutoff=aval_cutoff

;+
; NAME
;
;    wgfa_compare
;
; PROJECT
;
;    CHIANTI
;
; EXPLANATION
;
;    Compares two .wgfa files for the same ion, displaying
;    differences. The two files must have the same level index scheme.
;
; INPUTS
;
;    FILE1   The name of the reference .wgfa file.
;
;    FILE2   The name of the comparison .wgfa file.
;
; KEYWORDS
;
;    LIST_FILE1  Lists transitions in FILE1 that are not in FILE2.
;
;    LIST_FILE2  Lists transitions in FILE2 that are not in FILE1.
;
; OPTIONAL INPUTS
;
;    LIMIT   By default the routine shows transitions for which there is 
;            more than a 10% difference in A-values. This keyword allows 
;            you to change this. E.g., LIMIT=0.5 corresponds to 50%
;
;    LEV_CUTOFF If set, then only transitions with lower level
;               indices less than or equal to LEV_CUTOFF will be
;               considered.
;
;    AVAL_CUTOFF: Only A-values above this value will be displayed. The
;                 cutoff is applied to the stronger of the transitions
;                 from the two models.
;
; HISTORY
;
;    Ver.1, 3-Feb-03, Peter Young
;    Ver.2, 17-Feb-2009, Peter Young
;       added lev_cutoff= keyword
;    Ver.3, 6-Oct-2009, Peter Young
;       renamed as wgfa_compare
;    Ver.4, 26-Jun-2025, Peter Young
;       added aval_cutoff= optional input.
;-

IF n_params() LT 2 THEN BEGIN
  print,'Use: IDL> wgfa_compare, file1, file2, /list_file1, /list_file2,'
  print,'                          limit=, lev_cutoff=, aval_cutoff='
  return
ENDIF
 

IF n_elements(limit) EQ 0 THEN limit=0.1 ELSE limit=float(limit)

IF n_elements(aval_cutoff) EQ 0 THEN aval_cutoff=1e-10

read_wgfa_pry,file1,str1,ref
read_wgfa_pry,file2,str2,ref

IF keyword_set(lev_cutoff) THEN BEGIN
  i1=where(str1.lvl1 LE lev_cutoff)
  str1=str1[i1]
 ;
  i2=where(str2.lvl1 LE lev_cutoff)
  str2=str2[i2]
ENDIF 

n1=n_elements(str1)
n2=n_elements(str2)

flag1=bytarr(n1)
flag2=bytarr(n2)

print,format='(2a3,3a10,a8)','i','j','wvl','A_1','A_2','%diff'

FOR i=0,n1-1 DO BEGIN
  l1=str1[i].lvl1
  l2=str1[i].lvl2
  ind=where(str2.lvl1 EQ l1 AND str2.lvl2 EQ l2)
  IF ind[0] NE -1 THEN BEGIN
    IF n_elements(ind) GT 1 THEN print,'WARNING: transition duplicated!'
    aval1=str1[i].aval
    aval2=str2[ind[0]].aval
    tst=(aval2-aval1)/aval1
    IF abs(tst) GT limit AND max([aval1,aval2]) GE aval_cutoff THEN BEGIN
      print,format='(2i3,f10.3,2e10.2,f8.2)',l1,l2,str1[i].wvl,aval1,aval2,tst*100.
    ENDIF
  ENDIF ELSE BEGIN
    flag1[i]=1
  ENDELSE
ENDFOR

FOR i=0,n2-1 DO BEGIN
  l1=str2[i].lvl1
  l2=str2[i].lvl2
  ind=where(str1.lvl1 EQ l1 AND str1.lvl2 EQ l2)
  IF ind[0] EQ -1 THEN flag2[i]=1
END

ind=where(flag1 EQ 1)
n=n_elements(ind)
IF ind[0] EQ -1 THEN n=0
print,format='("Transitions in file1 not in file2: ",i3)',n
IF keyword_set(list_file1) AND n NE 0 THEN BEGIN
  FOR i=0,n-1 DO BEGIN
    print,format='(2i3,f15.3,e12.2)',str1[ind[i]].lvl1, $
         str1[ind[i]].lvl2,str1[ind[i]].wvl,str1[ind[i]].aval
  ENDFOR
ENDIF

ind=where(flag2 EQ 1)
n=n_elements(ind)
IF ind[0] EQ -1 THEN n=0
print,format='("Transitions in file2 not in file1: ",i3)',n
IF keyword_set(list_file2) AND n NE 0 THEN BEGIN
  FOR i=0,n-1 DO BEGIN
    print,format='(2i3,f15.3,e12.2)',str2[ind[i]].lvl1, $
         str2[ind[i]].lvl2,str2[ind[i]].wvl,str2[ind[i]].aval
  ENDFOR
ENDIF

END
