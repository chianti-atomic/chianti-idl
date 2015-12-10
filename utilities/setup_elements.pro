;+
;
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
;
; PURPOSE:  read in ionization equilibrium and abundances
;           not for general user use
;
;
; CALLING SEQUENCE:
;
;       SETUP_ELEMENTS
;
;
; INPUTS:
;
;	None	
;
; OUTPUTS:  reads data into common block
;
;
; COMMON BLOCKS:
;
;	common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref
;
;
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
; VERSION     :   3, 16-Apr-2003 
;
;-
pro setup_elements
;
common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref
;
;   read elemental abundances
;
dir=concat_dir(!xuvtop,'abundance')
dir=expand_path(dir)
abund_name=''
while abund_name eq '' do begin
   abund_name=dialog_pickfile(path=dir,filter='*.abund',title='Select Abundance File')
   if abund_name eq '' then print,' Please select abundance file'
endwhile
read_abund,abund_name,abund,abund_ref
;
;  read ionization equilibrium
;
dir=concat_dir(!xuvtop,'ioneq')
dir=expand_path(dir)
;
ioneq_name=''
while ioneq_name eq '' do begin
   ioneq_name=dialog_pickfile(path=dir,filter='*.ioneq', $
     title='Select Ionization Equilibrium File')
   if ioneq_name eq '' then print,' Please select ionization equilibrium file'
endwhile
read_ioneq,ioneq_name,ioneq_logt,ioneq,ioneq_ref
;
;
end
