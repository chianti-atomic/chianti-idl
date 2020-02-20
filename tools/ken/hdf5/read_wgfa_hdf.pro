   ;+
   ; NAME
   ;
   ;    READ_WGFA_HDF
   ;
   ; PROJECT
   ;
   ;    CHIANTI
   ;
   ; EXPLANATION
   ;
   ;    Reads the CHIANTI HDF file into a structure. For the most part
   ;    CHIANTI wgfa files contain radiative decay rates. I.e., rates for
   ;    transitions that yield photons of a definite wavelength. However
   ;    for some ions the .wgfa file is also used to contain
   ;    autoionization rates and/or two photon transitions. These are
   ;    needed for accurately modelling the level populations within
   ;    ions, but do not yield photons of a definite wavelength
   ;    (autoionization yields no photons; two photon transitions yield a
   ;    continuum which is separately modelled in CHIANTI).
   ;
   ;    These "radiationless" transitions are denoted in the .wgfa file
   ;    with a zero wavelength.
   ;
   ;    In the output structure all radiative decay rates are stored in
   ;    the .AVAL tag, while autoionization and two photon rates are
   ;    stored in the .AUTO tag.
   ;
   ; INPUTS
   ;
   ;    FILENAME  The name of the HDF5 file to be read.
   ;
   ; OUTPUTS
   ;
   ;    WGFASTR   A structure with the following tags:
   ;          .lvl1  The index of the lower level of the transition.
   ;          .lvl2  The index of the upper level of the transition.
   ;          .wvl   The wavelength of the transition (angstroms).
   ;          .gf    The weighted oscillator strength.
   ;          .aval  The radiative decay rate (s^-1).
   ;          .auto  The autoionization rate (s^-1).
   ;          
   ;   WGFAREF
   ;          Description of the atomic data
   ;
   ;
   ; HISTORY
   ;
   ;    Ver.1, 30-Jan-2017, Ken Dere
   ;        copied a lot of code from read_wgfa_str, written by Peter Young
   ;-
PRO read_wgfa_hdf, filename, wgfastr, wgfaref, verbose=verbose

   chck = file_search(filename)
   IF chck EQ '' THEN BEGIN
      print,'%READ_WGFA_STR: file not found. Returning...'
      return
   ENDIF

   hdf_id = h5f_open(filename)
   wgfa_group = h5g_open(hdf_id,'/wgfa')
   if keyword_set(verbose) then begin
      nmembers = h5g_get_nmembers(hdf_id,'/wgfa')
      for im = 0,nmembers-1 do begin
         amember = h5g_get_member_name(hdf_id,'/wgfa',im)
         print,' amember = ',amember
      endfor
   endif

   lvl1_dataset = h5d_open(hdf_id,'/wgfa/lvl1')
   lvl1 = h5d_read(lvl1_dataset)
   h5d_close, lvl1_dataset

   lvl2_dataset = h5d_open(hdf_id,'/wgfa/lvl2')
   lvl2 = h5d_read(lvl2_dataset)
   h5d_close, lvl2_dataset

   wvl_dataset = h5d_open(hdf_id,'/wgfa/wvl')
   wvl = h5d_read(wvl_dataset)
   h5d_close, wvl_dataset

   gf_dataset = h5d_open(hdf_id,'/wgfa/gf')
   gf = h5d_read(gf_dataset)
   h5d_close, gf_dataset

   pretty1_dataset = h5d_open(hdf_id,'/wgfa/pretty1')
   pretty1 = h5d_read(pretty1_dataset)
   h5d_close, pretty1_dataset

   pretty2_dataset = h5d_open(hdf_id,'/wgfa/pretty2')
   pretty2 = h5d_read(pretty2_dataset)
   h5d_close, pretty2_dataset

   avalue_dataset = h5d_open(hdf_id,'/wgfa/avalue')
   aval = h5d_read(avalue_dataset)
   h5d_close, avalue_dataset

   ref_dataset = h5d_open(hdf_id,'/wgfa/ref')
   wgfaref = h5d_read(ref_dataset)
   h5d_close, ref_dataset
   
   h5f_close,hdf_id
   ;
   ; Define structure for reading data, and read data in one go.
   ;
   str = {lvl1: 0, lvl2: 0, wvl: 0d, gf: 0d, aval: 0d, pretty1: '', pretty2: '' }
   ndata = n_elements(lvl1)
;   readstr = replicate(str,ndata)
;   
;   readstr.lvl1 = lvl1
;   readstr.lvl2 = lvl2
;   readstr.wvl = wvl
;   readstr.gf = gf
;   readstr.aval = aval
;   readstr.pretty1 = pretty1
;   readstr.pretty2 = pretty2

   str2={lvl1: 0, lvl2: 0, wvl: 0d, gf: 0d, aval: 0d, auto: 0d}
   wgfastr = replicate(str2,ndata)

;   wgfastr.lvl1 = lvl1
;   wgfastr.lvl2 = readstr.lvl2
;   wgfastr.wvl = readstr.wvl
;   wgfastr.gf = readstr.gf

   wgfastr.lvl1 = lvl1
   wgfastr.lvl2 = lvl2
   wgfastr.wvl = wvl
   wgfastr.gf = gf

   ;
   ; Only load A-values into wgfastr for those transitions with
   ; non-negative wavelengths.
   ;
   k=where(wvl NE 0.)
   wgfastr[k].aval = aval[k]

   ;
   ; The following loads the autoionization data into the output
   ; structure. The 'flag' array is used to flag lines of data that need
   ; to be removed from the output structure. The removed lines are those
   ; for which both a A-value and a autoionization rate are available.
   ;
   flag = bytarr(ndata)+1b
   k = where(wvl EQ 0., nk)
   IF nk NE 0 THEN BEGIN
      FOR i = 0,nk-1 DO BEGIN
;         l1 = readstr[k[i]].lvl1
;         l2 = readstr[k[i]].lvl2
         l1 = lvl1[k[i]]
         l2 = lvl2[k[i]]
         ;
         ii = where(lvl1 EQ l1 AND lvl2 EQ l2,nii)
         wgfastr[ii].auto = aval[k[i]]
         ;
         IF nii EQ 2 THEN flag[k[i]] = 0b
      ENDFOR
   ENDIF

   k = where(flag EQ 1b)
   wgfastr = wgfastr[k]

   end
