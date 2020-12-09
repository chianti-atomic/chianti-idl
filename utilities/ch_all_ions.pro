
FUNCTION ch_all_ions, count=count, sequence=sequence, element=element, $
                      diel_only=diel_only, no_diel=no_diel

;+
; NAME:
;     CH_ALL_IONS
;
; PURPOSE:
;     Returns a list of all the ions in the CHIANTI database. 
;
; CATEGORY:
;     CHIANTI; information.
;
; CALLING SEQUENCE:
;     Result = CH_ALL_IONS()
;
; INPUTS:
;     None.
;
; OPTIONAL INPUTS:
;     Sequence:  An integer or string specifying an isoelectronic
;                sequence. Only those ions belonging to the sequence
;                wil be returned. For example, to return the oxygen
;                isoelectronic sequence, input 8 or 'o'.
;     Element:   An integer or string specifying an element. Only
;                those ions belonging to the element will be
;                returned. For example, to return only oxygen ions,
;                input 8 or 'o'.
;
; KEYWORD PARAMETERS:
;     DIEL_ONLY: If set, then only return the dielectronic ions (i.e.,
;                those with 'd' appended to the ion name.
;     NO_DIEL:   If set, then do not return any dielectronic ions.
;
; OUTPUTS:
;     A string array containing the list of ions. The ion names are in
;     the CHIANTI format. For example, 'o_6' for O VI.
;
; OPTIONAL OUTPUTS:
;     Count:  The number of ions.
;
; EXAMPLE:
;     IDL> ions=ch_all_ions(count=count)
;     IDL> ions=ch_all_ions(element='fe')
;     IDL> ions=ch_all_ions(sequence=5)
;     IDL> ions=ch_all_ions(/diel_only)
;     IDL> ions=ch_all_ions(/no_diel)
;
; MODIFICATION HISTORY:
;     Ver.1, 09-Dec-2020, Peter Young
;-

count=0

IF n_elements(sequence) NE 0 AND n_elements(element) NE 0 THEN BEGIN
   print,'% CH_ALL_IONS: please specify ELEMENT= or SEQUENCE=, but not both. Returning...'
   return,-1
ENDIF 

IF n_elements(diel_only) NE 0 AND n_elements(no_diel) NE 0 THEN BEGIN
   print,'% CH_ALL_IONS: please specify /DIEL_ONLY or /NO_DIEL, but not both. Returning...'
   return,-1
ENDIF 

mlist_dir=concat_dir(!xuvtop,'masterlist')
mlist_file=concat_dir(mlist_dir,'masterlist.ions')

read_masterlist,mlist_file,mlist

;
; Here I trim the ion names in case there is any extra space.
;
n=n_elements(mlist)
z=intarr(n)
ion=intarr(n)
FOR i=0,n-1 DO BEGIN
   mlist[i]=trim(mlist[i])
   convertname,mlist[i],a,b
   z[i]=a
   ion[i]=b
ENDFOR 

;
; Handle the input 'sequence'.
;
seq=-1
IF n_elements(sequence) NE 0 THEN BEGIN
   IF datatype(sequence) EQ 'STR' THEN BEGIN
    z2element,indgen(30)+1,elt,/symbol,/lower_CASE
    k=where(strlowcase(sequence) EQ elt,nk)
    IF nk NE 0 THEN seq=k[0]+1
  ENDIF ELSE BEGIN 
    seq=sequence[0]
  ENDELSE
  k=where(z-ion+1 EQ seq,nk)
  IF nk NE 0 THEN mlist=mlist[k]
ENDIF 

;
; Handle the input 'element'.
;
elt_iz=-1
IF n_elements(element) NE 0 THEN BEGIN
  IF datatype(element) EQ 'STR' THEN BEGIN
    z2element,indgen(30)+1,elt,/symbol,/lower_CASE
    k=where(strlowcase(element) EQ elt,nk)
    IF nk NE 0 THEN elt_iz=k[0]+1
  ENDIF ELSE BEGIN 
    elt_iz=element[0]
  ENDELSE
  k=where(z EQ elt_iz,nk)
  IF nk NE 0 THEN mlist=mlist[k]
ENDIF

;
; Deal with dielectronic keywords.
;
chck=strpos(mlist,'d')
i_diel=where(chck GE 3,n_diel)
i_nodiel=where(chck LT 3,n_nodiel)
IF keyword_set(no_diel) THEN BEGIN
   IF n_nodiel NE 0 THEN BEGIN 
      mlist=mlist[i_nodiel]
   ENDIF ELSE BEGIN
      return,-1
   ENDELSE
ENDIF
;
IF keyword_set(diel_only) THEN  BEGIN
   IF n_diel NE 0 THEN BEGIN 
      mlist=mlist[i_diel]
   ENDIF ELSE BEGIN
      return,-1
   ENDELSE
ENDIF

count=n_elements(mlist)


return,mlist

END
