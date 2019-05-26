;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory (USA), the University of Florence (Italy), the
;       University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
;
; NAME:
;	READ_FBLVL
;
; PURPOSE:
;
;	To read files containing observed free-bound energy levels
;
; CATEGORY:
;
;	Science.
;
; CALLING SEQUENCE:
;
;       READ_FBLVL, File, L1, Conf, pqn, ll, spd, mult, Ecm, Ecmth, Ref
;
;
; INPUTS:
;
;	File:	the name of the file 
;               e.g., !xuvtop+'/si/si_12/si_12.fblvl' for Si XII
;
; OPTIONAL INPUTS:
;
;	None:
;	
;
; OUTPUTS:       L1    - level index
;                CONF  - configuration description
;                PQN   - principal quantum number of electron being ionized
;                LL    - L-value of ionized electron
;                SPD   - 'S', 'P', etc to denote L value
;                MULT  - multiplicity  2J+1
;                ECM   - energy (cm^-1)
;                ECMTH - energy (cm^-1)
;                REF   - reference
;               
; EXAMPLE:
;
;             > file = !xuvtop+'/si/si_12/si_12.elvlc'
;             > read_fblvl,file,l1,conf,pqn,ll,mult,ecm,ecmth,ref
;             > 
;
; HISTORY
;
;      Ver.1, 31-Jul-2002, Peter Young
;
;-

pro read_fblvl,filename,l1,conf,pqn,ll,spd,mult,ecm,ecmth,ref
;
;
if n_params(0) lt 9 then begin
     print,' '
     print,' IDL> read_fblvl,filename,l1,conf,pqn,ll,spd,mult,ecm,ecmth,ref'
     print,' '
     return
endif
;
;
openr,lue,filename,/get_lun


tst1=0
str1=''
datastring=''
ref=''
WHILE tst1 LT 2 DO BEGIN
  readf,lue,str1
  CASE 1 OF 
    strtrim(str1,2) NE '-1' AND tst1 EQ 0: datastring=[datastring,str1]

    strtrim(str1,2) NE '-1' AND tst1 EQ 1: ref=[ref,str1]

    strtrim(str1,2) EQ '-1' AND tst1 EQ 0: tst1=1

    strtrim(str1,2) EQ '-1' AND tst1 EQ 1: tst1=2
  ENDCASE
ENDWHILE

free_lun,lue

ref=ref[1:*]
datastring=datastring[1:*]
n=n_elements(datastring)

;
; using a structure the complete set of data can be read with one line of code
; without using a for loop
; 
str={l1: 0, conf: '', pqn: 0, ll: 0, spd: '', mult: 0, ecm: 0., ecmth: 0.}
struc=replicate(str,n)
reads,datastring,format='(i5,a20,2i5,a3,i5,2f20.3)',struc

l1=struc.l1
conf=struc.conf
pqn=struc.pqn
ll=struc.ll
mult=struc.mult
spd=struc.spd
ecm=struc.ecm
ecmth=struc.ecmth

END
