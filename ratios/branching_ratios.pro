
PRO branching_ratios, ion, $
                      photons=photons, wrange=wrange, all=all, $
                      extinct=extinct, path=path, level=level, $
                      maxlev=maxlev

;+
; NAME
;
;     BRANCHING_RATIOS
;
; PROJECT
;
;     CHIANTI
;
; EXPLANATION
;
;     For the specified ion this routine computes branching ratios
;     (these are ratios between two lines that emitted from the same
;     upper level; such ratios are insensitive to plasma conditions). 
;
; INPUTS
;
;     ION   The name of the ion in CHIANTI format. E.g., 'fe_13'.
;
; OPTIONAL INPUTS
;
;     PATH      By default the routine uses ion data from the
;               user's CHIANTI distribution. By setting PATH=
;               the user can directly specify a path to the data
;               files. 
;
;     WRANGE    Only branching ratios in the specified wavelength range are
;               printed. 
;
;     EXTINCT   Specify the interstellar extinction, E(B-V), to modify the
;               branching ratio value. The IDL astronomy routine FM_UNRED
;               is used to calculate the extinction.
;
;     LEVEL     Only show branching ratios from the specified level.
;
;     MAXLEV    Only show data for levels up to and including MAXLEV.
;
; KEYWORDS
;
;     PHOTONS   By default the branching ratios are given in energy units.
;               Setting this keyword gives them in photon units.
;
;     ALL       By default the routine only shows ratios > 0.001. By
;               setting /ALL all branching ratios are shown.
;
; CALLS
;
;     READ_WGFA_STR, READ_ELVLC, ION2FILENAME, FM_UNRED
;
; EXAMPLE
;
;     IDL> branching_ratios,'fe_13',lev=20     
;     
;     Wavelength     Ratio    Lower level
;     
;     Upper level  20  3s2.3p3d 3P1
;        202.044     1.000     1   3s2.3p2 3P0
;        205.915     0.006     2   3s2.3p2 3P1
;        209.919     0.150     3   3s2.3p2 3P2
;        223.777     0.015     4   3s2.3p2 1D2
;     
; HISTORY
;
;     Ver.1, 4-Nov-2003, Peter Young
;     Ver.2, 30-Apr-2007, Peter Young
;        introduced /NEG keyword, and corrected bug
;     Ver.3, 13-Feb-2009, Peter Young
;        switched to using read_wgfa_str; modified the printed output;
;        added LEVEL= and MAXLEV= keywords.
;     Ver.4, 20-Apr-2012, Peter Young
;        modified print behavior when /all is set.
;-

IF n_params() LT 1 THEN BEGIN
  print,'Use: IDL> branching_ratios, ion, /photons, wrange=, path=, /all'
  print,'                 extinct=, level=, maxlev='
  return
ENDIF

;
; Get the names of the .elvlc and .wgfa files
;
IF n_elements(path) EQ 0 THEN BEGIN
  ion2filename,ion,rootname
  wgfaname=rootname+'.wgfa'
  elvlcname=rootname+'.elvlc'
ENDIF ELSE BEGIN
  wgfaname=concat_dir(path,ion+'.wgfa')
  elvlcname=concat_dir(path,ion+'.elvlc')
ENDELSE 

IF n_elements(wrange) EQ 0 THEN wrange=[0,1d30]

IF n_elements(extinct) NE 0 THEN BEGIN
  print,''
  print,format='(" E(B-V)=",f7.3)',extinct
ENDIF

read_elvlc,elvlcname,l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,ref
read_wgfa_str,wgfaname,str

;
; Use optional input LEVEL to restrict the following for loop.
;
IF n_elements(level) GT 0 THEN BEGIN
  k=where(l1 EQ level)
  istart=k[0]
  iend=k[0]
ENDIF ELSE BEGIN
  istart=0
  iend=n_elements(l1)-1
 ;
  IF n_elements(maxlev) NE 0 THEN iend=maxlev-1
ENDELSE 

print,''
print,'Wavelength     Ratio    Lower level'

n=n_elements(l1)
FOR i=istart,iend DO BEGIN
 ;
 ; Find all lines that decay from the upper level and lie within WRANGE
 ;
  ind=where(str.lvl2 EQ l1[i] AND abs(str.wvl) GT wrange[0] AND abs(str.wvl) LE wrange[1],n_ind)
 ;
  IF n_ind GE 2 THEN BEGIN  ; more than one transition from upper level
    print,''
    print,'Upper level'+string(format='(i4)',l1[i])+'  '+trim(term[i])
    lev1=str[ind].lvl1
    lev2=str[ind].lvl2
   ;
    wvls=str[ind].wvl
    IF keyword_set(photons) THEN BEGIN
      avals=str[ind].aval
    ENDIF ELSE BEGIN
      avals=str[ind].aval/abs(wvls)
    ENDELSE
   ;
    k=where(avals EQ max(avals))
    aval_ref=avals[k[0]]
    wvl_ref=wvls[k[0]]
    FOR j=0,n_ind-1 DO BEGIN
      ratio=avals[j]/aval_ref
     ;
     ; Modify ratio for interstellar extinction
     ;
      IF n_elements(extinct) NE 0 THEN BEGIN
        fm_unred,wvls[j],1.0,extinct,funred
        ratio=ratio/funred
        fm_unred,wvl_ref,1.0,extinct,funred
        ratio=ratio*funred
      ENDIF
     ;
      CASE 1 OF
        ratio LT 0.001 AND NOT keyword_set(all): 
        ratio LT 0.001 AND  keyword_set(all): BEGIN
          print_str=string(format='(f10.3)',wvls[j])+ $
                    string(format='(e10.2)',ratio)+ $
                    string(format='(i6)',lev1[j])+ $
                    '   '+trim(term[lev1[j]-1])
          print,print_str
        END
        ELSE: BEGIN
         print_str=string(format='(f10.3)',wvls[j])+ $
                  string(format='(f10.3)',ratio)+ $
                  string(format='(i6)',lev1[j])+ $
                  '   '+trim(term[lev1[j]-1])
          print,print_str
       END
      ENDCASE 
         
    ENDFOR
  ENDIF
ENDFOR

IF NOT keyword_set(all) THEN BEGIN
  print,''
  print,'** Only ratios greater than 0.001 shown **'
ENDIF

END
