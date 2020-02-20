

PRO ch_ptable_pos, iz, ion, row=row, col=col

;+
; NAME:
;      CH_PTABLE_POS
;
; PURPOSE:
;      Given an isoelectronic sequence number, this routine returns
;      the position within the periodic table as a row and a column.
;
; CATEGORY:
;      CHIANTI; atomic; periodic table.
;
; CALLING SEQUENCE:
;      CH_PTABLE_POS, IZ, ION
;
; INPUTS:
;      Iz:   The atomic number of the species.
;
; OPTIONAL INPUTS:
;      Ion:  If specified, then it is interpreted as the spectroscopic
;            number of an ion (for example 13=XIII, hence a charge of
;            +12). This modifies the isoelectronic sequence number
;            such that iso_num=IZ-ION+1.
;
; OUTPUTS:
;      None. See optional outputs.
;
; OPTIONAL OUTPUTS:
;      Row:    The row that the species belongs to in the periodic table.
;      Col:    The column that the species belongs to in the periodic table.
;
; RESTRICTIONS:
;      Only works for sequences up to Zn (atomic number=30).
;
; EXAMPLE:
;      IDL> ch_ptable_pos, 14, row=row, col=col   ; silicon
;      IDL> print,row,col
;             3      14
;      IDL> ch_ptable_pos, 26, 3, row=row, col=col   ; Fe III (Cr-like)
;      IDL> print,row,col
;             4       6
;
; MODIFICATION HISTORY:
;      Ver.1, 25-Apr-2017, Peter Young
;-

IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> ch_ptable_pos, iz [, ion, row=row, col=col]'
  return
ENDIF 


IF n_elements(ion) NE 0 THEN pos=iz-ion+1 ELSE pos=iz

IF pos GT 30 THEN BEGIN
  print,'%CH_PTABLE_POS: This routine only works for sequences up to N=30. Returning...'
  return
ENDIF 

;
; Get row
;
CASE 1 OF
  pos LE 2: row=1
  pos GT 2 AND pos LE 10: row=2
  pos GT 10 AND pos LE 18: row=3
  pos GT 18 AND pos LE 30: row=4
ENDCASE


;
; Get column
;
all_cols=intarr(30)-1
i=indgen(3)*8
FOR j=0,1 DO all_cols[i+j+2]=1+j
i=indgen(2)*8
FOR j=0,5 DO all_cols[i+j+4]=13+j
all_cols[0:1]=[1,18]
all_cols[20:29]=indgen(10)+3
;
col=all_cols[pos-1]

END
