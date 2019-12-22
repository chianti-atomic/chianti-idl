;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an atomic database package for the calculation of
;       continuum and emission line spectra from astrophysical plasmas. It is a
;       collaborative project involving George Mason University, the University
;       of Michigan, the Naval Research Laboratory (USA), Cambridge University (UK),
;       the Arcetri Observatory (Italy).
;
;
; NAME:  CHIANTI_ION__DEFINE
;
; HISTORY
;
;       Ver.1, 6-jan-2011, Ken Dere
;           first release - defines the chianti_ion object
;
;
;       V.1, 6-jan-2011, Ken Dere
;-
;
;
pro chianti_ion::cleanup
  ptr_free,self.wgfa
  ptr_free,self.auto
  ptr_free,self.rrlvl
  ptr_free,self.elvlc
  ptr_free,self.splups
  ptr_free,self.psplups
  ptr_free,self.cisplups
  ptr_free,self.population
  ptr_free,self.emiss
  ptr_free,self.abund
  ptr_free,self.ioneq
end

pro chianti_ion::printit
  ; a junk but temporarily useful method
  print, ' ions = ',self.ions, self.z, self.stage
  mywgfa = *self.wgfa
  mysplups = *self.splups
  for i=0,4 do begin
    print, ' l1, ecm = ', (*self.elvlc).data[i].index, (*self.elvlc).data[i].energy
    print, ' lvl1, lvl2, avalue = ', (*self.wgfa).lvl1[i], (*self.wgfa).lvl2[i], (*self.wgfa).a_value[i]
    print ,' lvl1, lvl2,  splups.de = ', (*self.splups)[i].lvl1, (*self.splups)[i].lvl2, (*self.splups)[i].de
  endfor
  IF ptr_valid(ptr_psplups) THEN BEGIN
    npsplups = n_elements(*self.psplups)
    print, ' npsplups = ',npsplups
    for i=0,npsplups-1 do begin
      print ,' lvl1  psplups.de, tt = ',(*self.psplups)[i].lvl1,(*self.psplups)[i].de, (*self.psplups)[i].t_type
    endfor
  endif else begin
    print ,' psplups not valid/ does not exist '
  endelse
  return
end

function chianti_ion::getElvlc
  elvlc = *self.elvlc
  return, elvlc
end

function chianti_ion::getWgfa
  wgfa = *self.wgfa
  return, wgfa
end

function chianti_ion::getSplups
  splups = *self.splups
  return, splups
end

function chianti_ion::getReclvl
  reclvl = *self.reclvl
  return, reclvl
end

function chianti_ion::getRrlvl
  rrlvl = *self.rrlvl
  return, rrlvl
end

pro chianti_ion::ionizRate,temperature
  ioniz = ioniz_rate(self.ionS,temperature)
  ionizStr = {temperature:temperature, rate:ioniz}
  self.ioniz = ptr_new(ionizStr)
  return
end

function chianti_ion::getIonizRate
  ioniz = *self.ioniz
  return, ioniz
end

pro chianti_ion::recombRate,temperature
  recomb = recomb_rate(self.ionS,temperature)
  recombStr = {temperature:temperature, rate:recomb}
  self.recomb = ptr_new(recombStr)
  return
end

function chianti_ion::getRecombRate
  recomb = *self.recomb
  return, recomb
end


function chianti_ion::getPopulation
  return, *self.population
end

function chianti_ion::getEmiss
  return, *self.emiss
end

function chianti_ion::init,ionS
  constants = self.constants()
  self.ionS = ionS
  convertname,ionS, z, stage, diel=diel
  self.z = z
  self.stage = stage
  self.dielectronic = diel
  print, ' ion, z, stage, diel = ',ionS, z, stage, diel
  print, ' ion, z, stage, diel = ',self.ionS, self.z, self.stage, self.dielectronic
  ion2filename,ionS,fname
  print, ' fname = ',fname
  
  if file_test(fname+'.elvlc') then begin
  
    read_elvlc_str,fname+'.elvlc', elvlcstr,ref
    ptr_elvlc = ptr_new(elvlcstr)
    if ptr_valid(ptr_elvlc) then begin
      self.elvlc = ptr_elvlc
    endif else print, ' elvlc pointer not valid'
    ;
    read_wgfa2_str,fname+'.wgfa', wgfastr,ref
    ptr_wgfa = ptr_new(wgfastr)
    if ptr_valid(ptr_wgfa) then begin
      self.wgfa = ptr_wgfa
    endif else print, ' wgfa pointer not valid'
    ;
    ;
    if file_test(fname+'.splups') then begin
      print,' splups file found'
    endif else begin
      print,' splups file not found'
    endelse
    read_splups,fname+'.splups', splupsstr,ref
    ptr_splups = ptr_new(splupsstr)
    if ptr_valid(ptr_splups) then begin
      self.splups = ptr_splups
    endif else print, ' splups pointer not valid'
    ;
    if file_test(fname+'.psplups') then begin
      print, ' file exists - ',fname+'.psplups'
      read_splups,fname+'.psplups', psplupsstr,ref,/prot
      ptr_psplups = ptr_new(psplupsstr)
      if ptr_valid(ptr_psplups) then begin
        self.psplups = ptr_psplups
        self.pstatus = 1
      endif else print, ' psplups pointer not valid'
    endif else begin
      self.pstatus = 0
    endelse
    ;
    
    read_cireclvl, ionS, reclvlstr, /rec
    ptr_reclvl = ptr_new(reclvlstr)
    if ptr_valid(ptr_reclvl) then begin
      self.reclvl = ptr_reclvl
    endif else print, ' reclvl pointer not valid'
    
    
    read_cireclvl, ionS, rrlvlstr, /rr
    ptr_rrlvl = ptr_new(rrlvlstr)
    if ptr_valid(ptr_rrlvl) then begin
      self.rrlvl = ptr_rrlvl
    endif else print, ' rrlvl pointer not valid'
    
    autoName = fname+'.auto'
    if file_test(autoName) then begin
      print, ' file exists - ',fname+'.psplups'
      read_wgfa_str,autoName , autostr,ref
      ptr_auto = ptr_new(autostr)
      if ptr_valid(ptr_auto) then begin
        self.auto = ptr_auto
        self.autostatus = 1
      endif else print, ' auto pointer not valid - ',fname+'.auto
    endif else begin
      self.autostatus = 0
    endelse
    
    read_cireclvl, ionS, cilvl, /ci
    ptr_cilvl = ptr_new(cilvl)
    if ptr_valid(ptr_cilvl) then begin
      self.cilvl = ptr_cilvl
    endif else print,' cilvl pointer not valid'
    
    ioneqname = getenv('IONEQ')
    ioneqfilename = concat_dir(!xuvtop, 'ioneq') + '/'+ ioneqname
    read_ioneq, ioneqfilename, t, ioneq1, ioneqref
    ioneq = {temperature:t, ioneq:ioneq1, ref:ioneqref}
    ptr_ioneq = ptr_new(ioneq)
    if ptr_valid(ptr_ioneq) then begin
      self.ioneq = ptr_ioneq
    endif else print, ' ioneq pointer not valid'
    ;
    abundname = getenv('ABUND')
    abundfilename = concat_dir(!xuvtop ,'abundance') + '/' + abundname
    read_abund, abundfilename, abund1, abundref
    abund = {abund:abund1, ref:ioneqref}
    ptr_abund = ptr_new(abund)
    self.abund = ptr_abund
    
  endif ; end elvlc file exists
  
  
  return, 1
end
;
pro chianti_ion__define
  struct = {chianti_ion, z:0L, ions:' ', $
    stage:0L, $
    dielectronic:0L, $
    pstatus:0L, $
    autostatus:0L, $
    elvlc:ptr_new(/allocate_heap), $
    wgfa:ptr_new(/allocate_heap), $
    auto:ptr_new(/allocate_heap), $
    splups:ptr_new(/allocate_heap), $
    psplups:ptr_new(/allocate_heap), $
    cisplups:ptr_new(/allocate_heap), $
    ioniz:ptr_new(/allocate_heap), $
    recomb:ptr_new(/allocate_heap), $
    cilvl:ptr_new(/allocate_heap), $
    reclvl:ptr_new(/allocate_heap), $
    rrlvl:ptr_new(/allocate_heap), $
    population:ptr_new(/allocate_heap), $
    emiss:ptr_new(/allocate_heap), $
    abund:ptr_new(/allocate_heap), $
    ioneq:ptr_new(/allocate_heap)}
  print,' struct defined'
end

