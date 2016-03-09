
FUNCTION ch_find_iso_level, ionname, ref_ionname, ref_level, quiet=quiet, $
                            outlev_name=outlev_name, config_match=config_match

;+
; NAME:
;       CH_FIND_ISO_LEVEL
;
; PURPOSE: 
;       This routine tries to find a CHIANTI level for an ion, given the
;       ion and index of another ion on the sequence.
;
; CATEGORY:
;	CHIANTI.
;
; CALLING SEQUENCE:
;       Result = CH_FIND_ISO_LEVEL(Ionname, Ref_ionname, Ref_level)
;
; INPUTS:
;	Ionname: The name of the ion for which the level index is
;         	 needed. Specified in CHIANTI format, e.g., 'fe_13'
;         	 for Fe XIII.
;       Ref_ionname: The name of the reference ion for which the level
;                index is known. Specified in CHIANTI format, e.g.,
;                'fe_13' for Fe XIII.
;       Ref_level: The CHIANTI index for the level of Ref_ionname.
;
; KEYWORDS:
;       QUIET:  If set, then no information messages will be printed. 
;       CONFIG_MATCH: If set, then a configuration match is required
;               otherwise no match will be returned.
; 
; OPTIONAL OUTPUTS:
;       Outlev_name:  A string giving the name of the output level.
;
; OUTPUTS:
;	Returns the level index for the corresponding level in the
;	ion. If a level is not found, then a value of -1 is returned. 
;
; RESTRICTIONS:
;	The method may not work if there are many levels with the same
;	name in the same configuration.
;
; PROCEDURE:
;	The level is found by first identifying all levels with the
;	same, S, L, J and parity as the specified level. Amongst these
;	levels the routine simply searches for the one with the level
;	index nearest to the specified index.
;
; EXAMPLE:
;       IDL> lev = ch_find_iso_level('cl_16','ar_17',3
;       Level 3 of Ar XVII is 1s.2p 3P0
;       The corresponding level of Cl XVI is:
;            4  1s.2p 3P0
;
; MODIFICATION HISTORY:
; 	Written by:	Peter Young, 10-Feb-2014
;
;       Ver.1, 21-May-2014, Peter Young
;          Added /config_match.
;-

;
; Read the search ion data into the structure 'str'
;
convertname,ionname,iz,ion
diff=iz-ion
zion2filename,iz,ion,fname
elvlcname=fname+'.elvlc'
;
read_elvlc,elvlcname,elvlcstr=str


;
; Read the reference ion data into the structure 'refstr'. Also check
; that two ions belong to same sequence.
;
convertname,ref_ionname,iz,ion
IF iz-ion NE diff THEN BEGIN
  print,'% CH_FIND_ISO_LEVEL: The two ions do not belong to the same isoelectronic sequence. Returning...'
  return,-1
ENDIF 
zion2filename,iz,ion,fname
elvlcname=fname+'.elvlc'
;
read_elvlc,elvlcname,elvlcstr=refstr


;
; Extract data for the reference level
;
k=where(refstr.data.index EQ ref_level,nk)
IF nk EQ 0 THEN BEGIN
  print,'% CH_FIND_ISO_LEVEL: the reference level does not exist in the reference ion. Returning...'
  return,-1
ENDIF 
jref=refstr.data[k].j
sref=refstr.data[k].s
lref=refstr.data[k].l
pref=refstr.data[k].parity
eref=refstr.data[k].energy
cref=refstr.data[k].conf
IF NOT keyword_set(quiet) THEN print,'Level '+trim(ref_level)+' of '+refstr.info.ion_roman+' is '+refstr.data[k].full_level


;
; Find corresponding level in the search ion.
;
k=where(str.data.j EQ jref AND str.data.s EQ sref AND str.data.l EQ lref AND str.data.parity EQ pref,nk)
IF nk EQ 0 THEN BEGIN
  print,'% CH_FIND_ISO_LEVEL: the level can not be found. Returning...'
  return,-1
ENDIF ELSE BEGIN
  kc=where(str.data[k].conf EQ cref,nkc)
  IF keyword_set(config_match) AND nkc EQ 0 THEN return,-1
  IF nkc GE 1 THEN k=k[kc]
 ;
  getmin=min(abs(str.data[k].index-ref_level),imin)
  i=k[imin]
  IF NOT keyword_set(quiet) THEN BEGIN 
    print,'The corresponding level of '+str.info.ion_roman+' is:'
    print,'     '+trim(str.data[i].index)+'  '+trim(str.data[i].full_level)
  ENDIF 
  outlev_name=str.data[k[imin]].full_level
  return,str.data[k[imin]].index
ENDELSE

END
