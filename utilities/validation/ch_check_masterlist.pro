

PRO ch_check_masterlist

;+
; NAME:
;     CH_CHECK_MASTERLIST
;
; PURPOSE:
;     Compares the 'masterlist' ions against the contents of the
;     database to make sure they are consistent.
;
; CATEGORY:
;     CHIANTI; validation.
;
; CALLING SEQUENCE:
;     CH_CHECK_MASTERLIST
;
; INPUTS:
;     None.
;
; OUTPUTS:
;     Prints messages to the IDL window.
;
; MODIFICATION HISTORY:
;     Ver.1, 09-Dec-2020, Peter Young
;-



ml=ch_all_ions(/no_diel)
mld=ch_all_ions(/diel_only)

ion_add=''
ion_miss=''

;
; I go through every ion of every element up to zinc and check that
; the elvlc, scups and wgfa files exist. I then compare against the
; masterlist. I do the dielectronic ions separately.
;
FOR i=1,30 DO BEGIN
   FOR j=1,i DO BEGIN
      zion2name,i,j,name
      k=where(name EQ ml,nk)
      IF nk EQ 0 THEN swtch=0 ELSE swtch=1
     ;
      zion2filename,i,j,fname
      chck1=file_search(fname+'.elvlc',count=c1)
      chck2=file_search(fname+'.scups',count=c2)
      chck3=file_search(fname+'.wgfa',count=c3)
      IF c1+c2+c3 EQ 3 THEN BEGIN
         IF swtch EQ 0 THEN ion_add=[ion_add,name]
      ENDIF ELSE BEGIN
         IF swtch EQ 1 THEN ion_miss=[ion_miss,name]
      ENDELSE 
     ;
      zion2name,i,j,name,/diel
      k=where(name EQ mld,nk)
      IF nk EQ 0 THEN swtch=0 ELSE swtch=1
     ;
      zion2filename,i,j,fname,/diel
      chck1=file_search(fname+'.elvlc',count=c1)
      chck2=file_search(fname+'.scups',count=c2)
      chck3=file_search(fname+'.wgfa',count=c3)
      IF c1+c2+c3 EQ 3 THEN BEGIN
         IF swtch EQ 0 THEN ion_add=[ion_add,name]
      ENDIF ELSE BEGIN
         IF swtch EQ 1 THEN ion_miss=[ion_miss,name]
      ENDELSE 
   ENDFOR
ENDFOR

swtch1=0
IF n_elements(ion_add) GT 1 THEN BEGIN
   print,'The following ions are missing from masterlist:'
   FOR i=1,n_elements(ion_add)-1 DO print,format='("    ",a6)',ion_add[i]
   swtch1=1
ENDIF 

swtch2=0
IF n_elements(ion_miss) GT 1 THEN BEGIN
   print,'The following ions are in masterlist but missing from the database:'
   FOR i=1,n_elements(ion_miss)-1 DO print,format='("    ",a6)',ion_miss[i]
   swtch2=1
ENDIF

IF swtch1+swtch2 EQ 0 THEN BEGIN
   print,'% CH_CHECK_MASTERLIST: the masterlist and database are consistent with each other.'
ENDIF 

END
