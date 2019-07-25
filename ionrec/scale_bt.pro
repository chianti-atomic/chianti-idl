 ;+
; PROJECT     :  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory, George Mason University (USA), the University of
;	Florence (Italy), and the University of Cambridge(UK).
;       
;
; NAME
;
;     SCALE_BT
;
;  PURPOSE:
;
;     This routine scales ionization energies and cross sections using
;     a Burgess-Tully (1992) type scaling developed for ionization
;     
;
; INPUTS
;
;
;    ENERGY    Incident electron energy.
;    CROSS_SECTION  The ionization cross section
;    CROSS_SECTION_ERROR  Error values for the cross section (not necessary)
;    F         The scaling factor - arbitrary but must be greater than 1.0
;    EV1       The ionization potential in the same units as the energy
;    
;
; OUTPUTS
;
;    BTE       The scaled energy
;    BTCROSS   The scaled cross section
;    BTERR     The scaled cross section error
;
; OPTIONAL INPUTS
;
;    NONE 
;
; KEYWORDS
;
;    NONE 
;
; CALLS
;
;    NONE
;
; COMMON BLOCKS
;
;    NONE
;
; PROGRAMMING NOTES
;
;    NONE
;
; MODIFICATION HISTORY
;
;    Ver.1, 17-Nov-2006, Ken Dere
;
; VERSION     :  1, 17-Nov-2006
;
;-
PRO scale_bt,energy,cross,crosserr,f,ev1,bte,btcross,bterr
;
;
u=double(energy/ev1)
bte=1.-alog(f)/alog(u-1.+f)
btcross=u*cross*ev1^2/(alog(u)+1.)
;
if n_elements(crosserr) gt 0 then begin
   bterr=u*crosserr*ev1^2/(alog(u)+1.)
endif
;
end
