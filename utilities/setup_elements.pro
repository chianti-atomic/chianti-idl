;+
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory (USA), the University of Florence (Italy), the
;       University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
; NAME:
;	SETUP_ELEMENTS
;
; PURPOSE:
;       Loads element abundances and ionization fractions into a
;       common block that is used in other CHIANTI routines. The
;       abundance and ion fraction files are selected using widgets
;       unless the abund_file or ioneq_file inputs are used.
;
; CALLING SEQUENCE:
;       SETUP_ELEMENTS
;
; INPUTS:
;	None.
;
; OPTIONAL INPUTS:
;       Abund_file:  The name of a file containing an element
;                    abundance table in CHIANTI format.
;       Ioneq_file:  The name of a file containing an ionization
;                    fraction table in CHIANTI format.
;
; OUTPUTS:
;       Reads data into a common block.
;
; COMMON BLOCKS:
;	common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	Feb.  2000:     Version 1.0 for CHIANTI version 3
;       V. 2, 21-May-2002, Giulio Del Zanna (GDZ):
;                   generalized directory concatenation to work for
;                   Unix, Windows  and VMS. 
;
;       V. 3, 16-Apr-2003, Peter Young
;                   I had to use expand_path for DIR in order for the paths to be 
;                   recognised in Windows.
;
;       Ver.4, 15-Aug-2016, Peter Young
;          Added optional inputs ABUND_FILE and IONEQ_FILE.
;
; VERSION     :   4, 15-Aug-2016, Peter Young
;
;-
pro setup_elements, ioneq_file=ioneq_file, abund_file=abund_file
;
common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref


;
; Read element abundance file.
;
IF n_elements(abund_file) NE 0 THEN BEGIN
  abund_name=file_search(abund_file,count=n)
  IF n EQ 0 THEN BEGIN
    print,'%SETUP_ELEMENTS:  input abund_file not found.'
  ENDIF
ENDIF ELSE BEGIN
  abund_name=''
ENDELSE 
;
IF abund_name EQ '' THEN BEGIN
  dir=concat_dir(!xuvtop,'abundance')
  dir=expand_path(dir)
  WHILE abund_name EQ '' DO BEGIN
    abund_name=dialog_pickfile(path=dir,filter='*.abund',title='Select Abundance File')
    if abund_name EQ '' THEN print,' Please select abundance file'
  ENDWHILE
ENDIF 
read_abund,abund_name,abund,abund_ref

  
;
; Read ionization equilibrium file
;
IF n_elements(ioneq_file) NE 0 THEN BEGIN
  ioneq_name=file_search(ioneq_file,count=n)
  IF n EQ 0 THEN BEGIN
    print,'%SETUP_ELEMENTS:  input ioneq_file not found.'
  ENDIF
ENDIF ELSE BEGIN
  ioneq_name=''
ENDELSE 
;
IF ioneq_name EQ '' THEN BEGIN
  dir=concat_dir(!xuvtop,'ioneq')
  dir=expand_path(dir)
 ;
  WHILE ioneq_name EQ '' DO BEGIN
    ioneq_name=dialog_pickfile(path=dir,filter='*.ioneq', $
                               title='Select Ionization Equilibrium File')
    IF ioneq_name EQ '' THEN print,' Please select ionization equilibrium file'
  ENDWHILE
ENDIF 
read_ioneq,ioneq_name,ioneq_logt,ioneq,ioneq_ref

END 
