;+
; PROJECT     :  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory, George Mason University (USA), the University of
;	Florence (Italy), and the University of Cambridge (UK).
;       
;
; NAME
;
;     DESCALE_BT
;
;  PURPOSE:
;
;     This routine returns ionization energies and cross sections from scaled ionization energies 
;     and cross sections using a Burgess-Tully (1992) type scaling developed for ionization
;     
;     
;
; INPUTS
;
;
;    BTE       The scaled energy
;    BTCROSS   The scaled cross section
;    BTERR     The scaled cross section error
;    F         The scaling factor - arbitrary but must be greater than 1.0
;    EV1       The ionization potential in the same units as the original scaling
;    
;
; OUTPUTS
;
;    ENERGY    Incident electron energy in the same units as the orginal scaling.
;    CROSS_SECTION  The ionization cross section
;    CROSS_SECTION_ERROR  Error values for the cross section (not necessary)
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
PRO descale_bt,bte,btcross,bterr,f,ev1,energy,cross,crosserr
;
;
u=1.d -f +exp(alog(f)/(1.-bte)) ;  spaces are necessary in this line
energy=double(u*ev1)
cross=(alog(u)+1.d)*btcross/(u*ev1^2)
if n_elements(bterr) gt 0 then begin
   crosserr=(alog(u)+1.)*bterr/(u*ev1^2)
endif
;
end
