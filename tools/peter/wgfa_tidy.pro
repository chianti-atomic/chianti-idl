
PRO WGFA_TIDY, FILEIN, FILEOUT, ENFILE=ENFILE, ADD_INFO=ADD_INFO, $
               USE_RYD=USE_RYD, wvl_decimal=wvl_decimal

;+
; NAME
;
;       WGFA_TIDY
;
; PURPOSE
;
;	To reduce the size of Chianti .wgfa files
;
; EXPLANATION
;
;	Some of the Chianti .wgfa files contain transitions with very 
;	small A-values that do not affect the level balance of the ion 
;	and so can safely be removed from the file. In addition quite 
;	a few of the files give the level IDs for each of the transitions 
;	- information which is not read by the read_wgfa routine and so 
;	is unnecessary.
;
;	This routine looks at all the upper levels in the .wgfa file, 
;	finds the largest A-value from that level and retains only 
;	those transitions with A-values greater than or equal to 1/10^7 
;	of this value.
;
;       By specifying the .elvlc file (ENFILE), the routine will use 
;       these energies to compute the wavelengths in the .wgfa file. 
;       This serves to make the wavelengths consistent with the .elvlc 
;       file.
;
;       The wavelength calculation takes account of the third energy
;       column if it is present in the .elvlc file. The priorities for
;       working out the wavelength is as follows:
;
;        1. if observed energies exist for both levels, then use
;        these.
;        2. if an observed energy is not available for one or both
;        levels, then use the third energy column to compute
;        wavelengths. 
;        3. if only a theoretical energy (2nd energy column) is
;        available for either of the two levels, then compute the
;        wavelength using the theoretical energies for *both levels*. 
;
; INPUTS
;
;	FILEIN	The original Chianti .wgfa file
;	FILEOUT	The name of the revised Chianti file
;
;
; OPTIONAL INPUTS
;
;       ENFILE  Allows the input of a specified .elvlc file.
;
;       WVL_DECIMAL  By default the routine prints wavelengths to 3
;                    decimal places. Specifying wvl_decimal to be an
;                    integer, this can be over-ridden.
;
; KEYWORDS
;
;       ADD_INFO This adds the level information for each transition (Ken 
;                and Enrico prefer to have this information in each file).
;
; CALLS
;
;       PRY:  READ_BEST_ENERGS
;       Chianti:  READ_WGFA2, READ_ELVLC
;
; EXAMPLE
;
;	WGFA_TIDY, !xuvtop+'/fe/fe_14/fe_14.wgfa', 'fe_14.wgfa_new'
;
; HISTORY
;
;	Ver.1, 9-Feb-99, PRY
;       Ver.2, 24-Aug-00, PRY
;       Ver.3, 11-Dec-00, PRY
;                Added /USE_RYD keyword.
;       Ver.4, 5-Nov-02, PRY
;                Disabled /USE_RYD option following change to 
;                read_best_energs.pro 
;       Ver.5, 28-Jan-2008, PRY
;                Changed A-value cutoff from 1e-5 to 1e-7.
;       Ver.6, 12-Feb-2009, PRY
;                Added /get_lun in call to openw, and removed direct
;                lun specification.
;       Ver.7, 5-Sep-2012, PRY
;                Added WVL_DECIMAL= keyword; updated header.
;
; CONTACT
;
;	Peter Young, NRL, pyoung@ssd5.nrl.navy.mil
;-

IF n_params() LT 2 THEN BEGIN
  print,'Use:  IDL> wgfa_tidy, filein, fileout [, enfile=enfile, /add_info, /use_ryd, '
  print,'                          wvl_decimal= ]'
  return
ENDIF

;
; Print format for output file
;
IF n_elements(wvl_decimal) EQ 0 THEN wvl_decimal=3
format='(2i5,f15.'+trim(wvl_decimal)+',2e15.3)'

read_wgfa2,filein,l1,l2,wvl,gf,a,ref

IF N_ELEMENTS(enfile) NE 0 THEN BEGIN
  read_best_energs,enfile,ec,w,en_flag,ethc,ethr, enryd=enryd
  read_elvlc,enfile,lev,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth, $
       refen
ENDIF

IF keyword_set(add_info) AND (N_ELEMENTS(enfile) EQ 0) THEN BEGIN
  print,'** please give the name of a .elvlc file through the ENFILE keyword'
  return
ENDIF


OPENW,lout,fileout,/get_lun

mind=MAX([l1,l2])

FOR i=1,mind DO BEGIN
  ind=WHERE(l2 EQ i)
  IF ind[0] NE -1 THEN BEGIN
    avals=a[ind]
    m_aval=MAX(avals)
    ind2=WHERE(avals GE m_aval*10.^(-7))
    n=N_ELEMENTS(ind2)
    FOR j=1,n DO BEGIN
      neg = +1.
      index = ind[ind2[j-1]]
      wavel = wvl[index]
     ;
      IF n_elements(enfile) NE 0 THEN BEGIN
        ll=l1[index]
        en_flag1=en_flag[ll-1]
        en_flag2=en_flag[i-1]
        IF (en_flag1 NE 2) AND (en_flag2 NE 2) THEN BEGIN
          IF keyword_set(use_ryd) THEN BEGIN
            wavel=1d8/(109737.32d0*abs(er[ll-1]-er[i-1]))
          ENDIF ELSE BEGIN
            wavel=1d8/(abs(ec[ll-1]-ec[i-1]))
          ENDELSE
        ENDIF ELSE BEGIN
          IF keyword_set(use_ryd) THEN BEGIN
            wavel=1d8/(109737.32d0*abs(ethr[ll-1]-ethr[i-1]))
          ENDIF ELSE BEGIN
            wavel=1d8/(abs(ethc[ll-1]-ethc[i-1]))
          ENDELSE
        ENDELSE
        IF (en_flag1 NE 0) OR (en_flag2 NE 0) THEN wavel=-wavel
      ENDIF
     ;
      str_print = STRING(FORMAT=format,l1[index], $
               l2[index],wavel,gf[index],a[index])
      IF keyword_set(add_info) THEN BEGIN
         string1 = STRTRIM(term[l1[index]-1],2)+' - '+ $
                   STRTRIM(term[l2[index]-1],2)
         str_print = str_print+'   '+string1
      ENDIF
      printf,lout,str_print
     ;
    ENDFOR
  ENDIF
ENDFOR

PRINTF,lout,FORMAT='(i3)',-1
n=N_ELEMENTS(ref)
FOR i=0,n-1 DO PRINTF,lout,ref(i)
PRINTF,lout,FORMAT='(i3)',-1

CLOSE,lout

END
