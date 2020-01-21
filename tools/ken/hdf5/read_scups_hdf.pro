;+
; NAME
;
;     READ_SCUPS_HDF
;
; PROJECT
;
;     CHIANTI
;
; EXPLANATION
;
;     Reads the CHIANTI HDF file containing the  SCUPS data.
;   
;
; INPUTS
;
;     FILENAME   name of the HDF file to read.
;
; KEYWORDS
;
;     None
;
; OUTPUTS
;
;     SPLSTR   A structure containing the spline data. There are two
;              tags called 'info' and 'data' that are each structures.
;
;              UPSSTR.INFO has the following tags:
;
;              ION_NAME        STRING    'ca_2'
;              ION_Z           INT             20
;              ION_N           INT              2
;              ION_ROMAN       STRING    'Ca II'
;              ION_LATEX       STRING    '\ion{Ca}{ii}'
;              ION_LATEX_ALT   STRING    '\ion{Ca}{2}'
;              COMMENTS        STRING    Array[12]
;              CHIANTI_VER     STRING    '7.1'
;              TIME_STAMP      STRING    'Fri May 24 16:54:46 2013'
;              FILENAME        STRING    'ca_2.scups_auto'
;              MISSING         FLOAT          -1.00000
;              NTRANS          LONG               766
;
;              UPSSTR.DATA has the following tags:
;
;              LVL1            INT              1
;              LVL2            INT              2
;              T_TYPE          INT              2
;              DE              FLOAT          0.124400
;              GF              FLOAT           0.00000
;              LIM             FLOAT          -1.00000
;              C_UPS           FLOAT          0.560000
;              NSPL            INT             13
;              STEMP           FLOAT     Array[18]
;              SPL             FLOAT     Array[18]
;
; HISTORY
;
;     Ver. 1, 31-jan-2017, Ken Dere
;        copied a lot of code from read_scups written by Peter Young and Giulio Del Zanna
;
; 
;-
PRO read_scups_hdf, filename, scupstr, scupsref, verbose=verbose


   IF n_params() LT 1 THEN BEGIN
      print,'Use:  IDL>  read_scups_hdf, splfile, splstr'
      return
   ENDIF
   
   t1=systime(/sec)

   chck = file_search(filename)
   IF chck EQ '' THEN BEGIN
      print,'%READ_scups_STR: file not found. Returning...'
      return
   ENDIF

   hdf_id = h5f_open(filename)
   scups_group = h5g_open(hdf_id,'/scups')
   
   if keyword_set(verbose) then begin
      nmembers = h5g_get_nmembers(hdf_id,'/scups')
      for im = 0,nmembers-1 do begin
         amember = h5g_get_member_name(hdf_id,'/scups',im)
         print,' amember = ',amember
      endfor
   endif

   lvl1_dataset = h5d_open(hdf_id,'/scups/lvl1')
   lvl1 = h5d_read(lvl1_dataset)
   h5d_close, lvl1_dataset

   lvl2_dataset = h5d_open(hdf_id,'/scups/lvl2')
   lvl2 = h5d_read(lvl2_dataset)
   h5d_close, lvl2_dataset

   ttype_dataset = h5d_open(hdf_id,'/scups/ttype')
   ttype = h5d_read(ttype_dataset)
   h5d_close, ttype_dataset

   de_dataset = h5d_open(hdf_id,'/scups/de')
   de = h5d_read(de_dataset)
   h5d_close, de_dataset

   gf_dataset = h5d_open(hdf_id,'/scups/gf')
   gf = h5d_read(gf_dataset)
   h5d_close, gf_dataset

   lim_dataset = h5d_open(hdf_id,'/scups/lim')
   lim = h5d_read(lim_dataset)
   h5d_close, lim_dataset

   ntemp_dataset = h5d_open(hdf_id,'/scups/ntemp')
   ntemp = h5d_read(ntemp_dataset)
   h5d_close, ntemp_dataset

   cups_dataset = h5d_open(hdf_id,'/scups/cups')
   cups = h5d_read(cups_dataset)
   h5d_close, cups_dataset

   btemp_dataset = h5d_open(hdf_id,'/scups/btemp')
   btemp = h5d_read(btemp_dataset)
   h5d_close, btemp_dataset

   bscups_dataset = h5d_open(hdf_id,'/scups/bscups')
   bscups = h5d_read(bscups_dataset)
   h5d_close, bscups_dataset

   ref_dataset = h5d_open(hdf_id,'/scups/ref')
   scupsref = h5d_read(ref_dataset)
   h5d_close, ref_dataset
   
   h5f_close,hdf_id
; nt_max (maximum number of temperatures for a transition) sets the
; size of the temperature and upsilon arrays 
;
   nt_max = max(ntemp)

;
; Create the output data structure.
;
   str2 = {lvl1: 0, lvl2: 0, t_type: 0, de: 0., gf: 0., lim: 0., c_ups: 0., nspl: 0, $
      stemp: fltarr(nt_max)-1, spl: fltarr(nt_max)-1}
   ntrans = n_elements(lvl1)
   datastr = replicate(str2,ntrans)
;
   datastr.lvl1 = lvl1
   datastr.lvl2 = lvl1
   datastr.de = de
   datastr.gf = gf
   datastr.lim = lim
   datastr.c_ups = cups
   datastr.nspl = ntemp
   datastr.t_type = ttype
   datastr.stemp = btemp
   datastr.spl = bscups
   ;
   ; Get the CHIANTI version number
   ;
   vfile=concat_dir(!xuvtop,'VERSION')
   chck=file_search(vfile)
   vnum=''
   IF chck[0] NE '' THEN BEGIN
      openr,lin,vfile,/get_lun
      readf,lin,vnum
      free_lun,lin
      vnum=strtrim(vnum,2)
   ENDIF
   ;
   missing_val=-1
   ;
   ; Create the various ion identifiers from the filename.
   ;
   basename=file_basename(filename)
   bits=str_sep(basename,'.')
   ion_name=bits[0]
   convertname,ion_name,ion_z,ion_n
   ion_z=fix(ion_z)
   zion2spectroscopic,ion_z,ion_n,ion_roman
   ion_roman=strcompress(ion_roman)
   bits=strsplit(ion_roman,' ',/extract)
   ion_latex='\ion{'+bits[0]+'}{'+strlowcase(bits[1])+'}'
   ion_latex_alt='\ion{'+bits[0]+'}{'+trim(ion_n)+'}'

   info={ion_name: ion_name, ion_z: ion_z, ion_n: ion_n, ion_roman: ion_roman, $
      ion_latex: ion_latex, ion_latex_alt: ion_latex_alt, comments: scupsref, $
      chianti_ver: vnum, time_stamp: systime(), $
      filename: filename, missing: float(missing_val), $
      ntrans: ntrans}
      
   ; Combine the info and data structures into the final output
   ; structure.
   ;
   scupstr={info: info, data: datastr}

   if keyword_set(verbose) then begin
      t2=systime(/sec)
      print,filename+', seconds:', t2-t1
   endif

   end
