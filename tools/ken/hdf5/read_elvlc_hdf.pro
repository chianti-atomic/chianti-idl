;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving 
;       George Mason University USA), the University of Michigan (USA), 
;       and Cambridge University (UK).
;
;
; NAME:
;  READ_ELVLC_HDF
;
; PURPOSE:
;
;  to read HDF file containing the energy level data
;
; CATEGORY:
;
;  science.
;
; CALLING SEQUENCE:
;
;       ELVLCSTR = READ_ELVLC_HDF(Filename, Verbose)
;
;
; INPUTS:
;
;  Filename: the name of the input HDF file, i.e. !xuvtop/si/si_4/si_4.hdf
;
;
; OUTPUTS:
;
;     ELVLCSTR  A structure array containing data for each energy
;               level. The structure has two tags called INFO and
;               DATA. The tags of ELVLCSTR.INFO are:
;
;               ION_NAME        STRING    'fe_10'
;               ION_Z           INT             26
;               ION_N           INT             10
;               ION_ROMAN       STRING    'Fe X'
;               ION_LATEX       STRING    '\ion{Fe}{x}'
;               ION_LATEX_ALT   STRING    '\ion{Fe}{10}'
;               COMMENTS        STRING    Array[10]
;               CHIANTI_VER     STRING    '7.1'
;               TIME_STAMP      STRING    'Tue Dec  4 13:07:39 2012'
;               FILENAME        STRING    'fe_10.elvlc'
;
;               ELVLCSTR.DATA is a structure array (with an entry for
;               each level), and the tags are:
;
;               INDEX           INT              8
;               CONF            STRING    '3s2.3p4(3P).3d'
;               CONF_LATEX      STRING    '3s$^2$ 3p$^4$ ($^3$P) 3d'
;               CONF_INDEX      INT              3
;               TERM            STRING    '4F'
;               TERM_LATEX      STRING    '$^4$F'
;               LEVEL           STRING    '4F9/2'
;               LEVEL_LATEX     STRING    '$^4$F$_{9/2}$'
;               FULL_LEVEL      STRING    '3s2.3p4(3P).3d 4F9/2'
;               FULL_LEVEL_LATEX
;                               STRING    '3s$^2$ 3p$^4$ ($^3$P) 3d $^4$F$_{9/2}$'
;               LABEL           STRING    ''
;               MULT            INT              4
;               S               FLOAT           1.50000
;               L               INT              3
;               L_SYM           STRING    'F'
;               J               FLOAT           4.50000
;               J_STR           STRING    '9/2'
;               PARITY          INT              0
;               PARITY_STR      STRING    'e'
;               WEIGHT          FLOAT           10.0000
;               OBS_ENERGY      DOUBLE           417652.00
;               THEORY_ENERGY   DOUBLE          -1.0000000
;               ENERGY          DOUBLE           417652.00                                   
;               ENERGY_UNITS    STRING    'cm^-1'
;               no longer included:   EXTRA_INFO      STRING    ''
;
;
;
; KEYWORDS
;
;  Verbose: a keyword, if set, some information will be printed
;
;
; PROCEDURE:
;
;  see Burgess and Tully, 1992, Astronomy and Astrophysics, 254, 436.
;
; EXAMPLE:
;
;       > elvlcstr = read_elvlc_hdf(!xuvtop+'/si/si_4/si_4.hdf')
;
; PROGRAMMING NOTES
;
;       This routine is to read the CHIANTI HDF files containing the
;       energy level information.  The goal is to spead up read time.
;
; HISTORY
;
;     Ver.1, 22-Jan-2017, Ken Dere
;        copied a lot of code from read_elvlc_str.pro written by Peter Young and Giulio Del Zanna
;        extra_info tag has been removed
;        the comments/reference strings may need some work
;
;-
pro read_elvlc_hdf, filename, elvlcstr, elvlcref, verbose=verbose
   hdf_id = h5f_open(filename)
   elvlc_group = h5g_open(hdf_id,'/elvlc')
   if keyword_set(verbose) then begin
      nmembers = h5g_get_nmembers(hdf_id,'/elvlc')
      for im = 0,nmembers-1 do begin
         amember = h5g_get_member_name(hdf_id,'/elvlc',im)
         print,' amember = ',amember
      endfor
   endif
   
   lvl_dataset = h5d_open(hdf_id,'/elvlc/lvl')
   i = h5d_read(lvl_dataset)   
   h5d_close, lvl_dataset
      
   term_dataset = h5d_open(hdf_id,'/elvlc/term')
   conf = h5d_read(term_dataset)
   h5d_close, term_dataset

   label_dataset = h5d_open(hdf_id,'/elvlc/label')
   lbl = h5d_read(label_dataset)
   h5d_close, label_dataset

   spin_dataset = h5d_open(hdf_id,'/elvlc/spin')
   ss = h5d_read(spin_dataset)
   h5d_close, spin_dataset

   spd_dataset = h5d_open(hdf_id,'/elvlc/spd')
   llstr = h5d_read(spd_dataset)
   h5d_close, spd_dataset

;   l_dataset = h5d_open(hdf_id,'/elvlc/l')
;   l = h5d_read(l_dataset)
;   h5d_close, l_dataset

   j_dataset = h5d_open(hdf_id,'/elvlc/j')
   jj = h5d_read(j_dataset)
   h5d_close, j_dataset

   ecm_dataset = h5d_open(hdf_id,'/elvlc/ecm')
   obs_en = h5d_read(ecm_dataset)
   h5d_close, ecm_dataset

   ecmth_dataset = h5d_open(hdf_id,'/elvlc/ecmth')
   th_en = h5d_read(ecmth_dataset)
   h5d_close, ecmth_dataset
   
   ref_dataset = h5d_open(hdf_id,'/elvlc/ref')
   elvlcref = h5d_read(ref_dataset)
   h5d_close, ref_dataset

   h5f_close,hdf_id
   
   ndata = n_elements(i)
   if keyword_set(verbose) then begin
      print,' ndata = ',ndata
   endif
   
   str = {i: 0, conf: '', lbl: '',ss: 0, llstr: '', jj: 0., obs_en: 0d0, th_en: 0d0}
   readstr = replicate(str,ndata)
   
   readstr.i = i
   readstr.conf = conf
   readstr.lbl = lbl
   readstr.ss = ss
   readstr.llstr = llstr
   readstr.jj = jj
   readstr.obs_en = obs_en
   readstr.th_en = th_en

   ;
   ; Extract information about the ion that will go into the output structure.
   ;
   fname=file_basename(filename)
   bits=strsplit(fname,'.',/extract)
   ionname=bits[0]
   convertname,ionname,iz,ion
   IF iz NE 0 AND ion NE 0 THEN BEGIN
      zion2spectroscopic,iz,ion,ion_roman
      iz=fix(iz)
      ion_roman=strcompress(ion_roman)
      bits=strsplit(ion_roman,' ',/extract)
      ion_latex='\ion{'+bits[0]+'}{'+strlowcase(bits[1])+'}'
      ion_latex_alt='\ion{'+bits[0]+'}{'+trim(ion)+'}'
   ENDIF ELSE BEGIN
      print,'%READ_ELVLC:  ion name can not be extracted from filename. Ion name parameters will not be set in output structure.'
      ionname=''
      iz=0
      ion=0
      ion_roman=''
      ion_latex=''
      ion_latex_alt=''
   ENDELSE


   spd=['S','P','D','F','G','H','I','K','L','M','N','O','Q','R','T','U','V','W','X','Y']
   
   vfile=concat_dir(!xuvtop,'VERSION')
   chck=file_search(vfile)
   vnum=''
   IF chck[0] NE '' THEN BEGIN
      openr,lin,vfile,/get_lun
      readf,lin,vnum
      free_lun,lin
      vnum=strtrim(vnum,2)
   ENDIF

   n=n_elements(readstr)

   newstr={ index: 0, $
      conf: '', $
      conf_latex: '', $
      conf_index: 0, $
      term: '', $
      term_latex: '', $
      level: '', $
      level_latex: '', $
      full_level: '', $
      full_level_latex: '', $
      label: '', $
      mult: 0, $
      s: 0.0, $
      l: 0, $
      l_sym: '', $
      j: 0., $
      j_str: '', $
      parity: 0, $
      parity_str: '', $
      weight: 0., $
      obs_energy: 0d0, $
      theory_energy: 0d0, $
      energy: 0d0, $
      energy_units: 'cm^-1'} ;, $
;      extra_info: '' kpd - removed
   levstr=replicate(newstr,n)
   
;   ref = ' ref'

   elvlcstr={ ion_name: ionname, $
      ion_z: iz, $
      ion_n: ion, $
      ion_roman: strcompress(ion_roman), $
      ion_latex: ion_latex, $
      ion_latex_alt: ion_latex_alt, $
      comments: elvlcref, $
      chianti_ver: vnum, $
      time_stamp: systime(), $
      data: levstr}

   elvlcstr.data.index=readstr.i
   elvlcstr.data.conf=trim(readstr.conf)
   elvlcstr.data.label=trim(readstr.lbl)
   elvlcstr.data.mult=readstr.ss
   elvlcstr.data.l_sym=trim(readstr.llstr)
   elvlcstr.data.j=readstr.jj
   elvlcstr.data.obs_energy=readstr.obs_en
   elvlcstr.data.theory_energy=readstr.th_en
;   elvlcstr.data.extra_info=strtrim(data_string_extra,2)
   
   end