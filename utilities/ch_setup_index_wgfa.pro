
FUNCTION ch_setup_index_wgfa, wgfastr, wvlmin=wvlmin, wvlmax=wvlmax, levmax=levmax, $
                              count=count, obs_only=obs_only

;+
; NAME:
;     CH_SETUP_INDEX_WGFA
;
; PURPOSE:
;     Derive the index_wgfa array that is used by ch_synthetic.
;
; CATEGORY:
;     CHIANTI.
;
; CALLING SEQUENCE:
;     Result = CH_SETUP_INDEX_WGFA( Wgfastr )
;
; INPUTS:
;     Wgfastr:  The structure returned by read_wgfa_str.pro.
;
; OPTIONAL INPUTS:
;     Wvlmin:  A wavelength in angstroms. If set, then the routine
;	       checks if the ion has any wavelengths above
;	       WVLMIN. If not, then the routine exits, and an empty
;	       output is returned.
;     Wvlmax:  A wavelength in angstroms. If set, then the routine
;	       checks if the ion has any wavelengths below
;	       WVLMAX. If not, then the routine exits, and an empty
;	       output is returned.
;     Levmax:  An integer specifying the highest atomic level to
;              include in the output.
;	
; KEYWORD PARAMETERS:
;     OBS_ONLY: If a wavelength check is performed (WVLMIN and/or
;               WVLMAX), then the default is to consider observed
;               and theoretical wavelengths. With /OBS_ONLY, only
;               the observed wavelengths are considered.
;
; OUTPUTS:
;     An index of the WGFA structure that picks out those transitions
;     that satisfy the WVLMIN and/or WVLMAX conditions, the A_VALUE NE 0
;     condition, and the LEVMAX condition. This index is used by the
;     routine CH_SYNTHETIC routine (see "anylines" in this  routine).
;
;     If no entries are found, then -1 is returned.
;
; OPTIONAL OUTPUTS:
;     Count:  The number of entries in the output index.
;
; EXAMPLE:
;     IDL> ion2filename,'o_6',fname
;     IDL> read_wgfa_str,fname+'.wgfa',wgfa
;     IDL> ind=ch_setup_index_wgfa(wgfa,wvlmin=1030,wvlmax=1040)
;     IDL> ind=ch_setup_index_wgfa(wgfa,levmax=3,count=count)
;
; MODIFICATION HISTORY:
;     Ver.1, 12-May-2023, Peter Young
;-


IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> ind=ch_setup_index_wgfa(wgfa [, wvlmin=, wvlmax=, levmax=, count= '
  print,'                          /obs_only ] )'
  return,-1
ENDIF 

IF n_elements(levmax) EQ 0 THEN levmax=max(wgfastr.lvl2)

IF keyword_set(obs_only) THEN wvlchck=wgfastr.wvl ELSE wvlchck=abs(wgfastr.wvl)

CASE 1 OF
  n_elements(wvlmin) NE 0 AND n_elements(wvlmax) EQ 0: $
     k=where(wvlchck GE wvlmin AND wgfastr.aval NE 0. AND wgfastr.lvl2 LE levmax,nk)
  n_elements(wvlmin) EQ 0 AND n_elements(wvlmax) NE 0: $
     k=where(wvlchck LE wvlmax AND wgfastr.aval NE 0. AND wgfastr.lvl2 LE levmax,nk)
  n_elements(wvlmin) NE 0 AND n_elements(wvlmax) NE 0: $
     k=where(wvlchck LE wvlmax AND wvlchck GE wvlmin AND wgfastr.aval NE 0. AND wgfastr.lvl2 LE levmax,nk)
  ELSE: BEGIN
    k=where(wgfastr.lvl2 LE levmax,nk)
  ENDELSE 
ENDCASE
index_wgfa=k

count=nk

return,index_wgfa


END
