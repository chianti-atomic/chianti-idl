;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an atomic database package for the calculation of
;       continuum and emission line spectra from astrophysical plasmas. It is a 
;       collaborative project involving George Mason University, the University
;		of Michigan, the Naval Research Laboratory (USA), Cambridge University (UK),
;		the Arcetri Observatory (Italy).
;
;
; NAME:  CHIANTI_ION__DEFINE
;
; HISTORY
;
;       Ver.1, 6-jan-2011, Ken Dere
;			first release - defines the chianti_ion object
;
;
;		V.1, 6-jan-2011, Ken Dere
;-
;
pro chianti_ion::cleanup
    ptr_free,self.wgfa
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
        print, ' l1, ecm = ',(*self.elvlc).l1[i],(*self.elvlc).ecm[i]
        print, ' lvl1, avalue = ',(*self.wgfa).lvl1[i],(*self.wgfa).a_value[i]
        print ,' lvl1  splups.de = ',(*self.splups)[i].lvl1,(*self.splups)[i].de
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
    self.ionS = ionS
    convertname,ionS, z, stage
    self.z = z
    self.stage = stage
    print, ' ion, z, stage = ',ionS,z,stage
    print, ' ion, z, stage = ',self.ionS,self.z,self.stage
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
        read_splups,fname+'.splups', splupsstr,ref
        ptr_splups = ptr_new(splupsstr)
        if ptr_valid(ptr_splups) then begin
            self.splups = ptr_splups
        endif else print, ' splups pointer not valid'
    ;
        read_splups,fname+'.psplups', psplupsstr,ref,/prot
        print,' psplups.status = ',psplupsstr.status[0]
        if psplupsstr[0].status gt 0 then begin
			ptr_psplups = ptr_new(psplupsstr)
			if ptr_valid(ptr_psplups) then begin
				self.psplups = ptr_psplups
			endif else print, ' psplups pointer not valid'
		endif
        
    endif
    
    read_cireclvl,ionS, reclvl, /rec
    print,' reclvl.status = ',reclvl.status
    if reclvl.status > 0 then begin
		ptr_reclvl = ptr_new(reclvl)
		if ptr_valid(ptr_reclvl) then begin
			self.reclvl = ptr_reclvl
		endif else print, ' reclvl pointer not valid'
	endif
    
    read_cireclvl, ionS, cilvl, /ci
    print,' cilvl.status = ',cilvl.status
    if cilvl.status > 0 then begin
		ptr_cilvl = ptr_new(cilvl)
		if ptr_valid(ptr_cilvl) then begin
			self.cilvl = ptr_cilvl
		endif else print,' cilvl pointer not valid'
	endif
	
; 	ioneqname = getenv('IONEQ')
; 	ioneqfilename = concat_dir(!xuvtop, 'ioneq') + '/'+ ioneqname
; 	read_ioneq, ioneqfilename, t, ioneq1, ioneqref
; 	ioneq = {temperature:t, ioneq:ioneq1, ref:ioneqref}
; 	ptr_ioneq = ptr_new(ioneq)
; 	if ptr_valid(ptr_ioneq) then begin
; 		self.ioneq = ptr_ioneq
; 	endif else print, ' ioneq pointer not valid'
; 	
; 	abundname = getenv('ABUND')
; 	abundfilename = concat_dir(!xuvtop ,'abundance') + '/' + abundname
; 	read_abund, abundfilename, abund1, abundref
; 	abund = {abund:abund1, ref:ioneqref}
; 	ptr_abund = ptr_new(abund)
; 	self.abund = ptr_abund
	
    
    return,1
end

pro chianti_ion__define
    struct = {chianti_ion, z:1L, ions:' ', $
    stage:0L, $
    elvlc:ptr_new(/allocate_heap), $
    wgfa:ptr_new(/allocate_heap), $
    splups:ptr_new(/allocate_heap), $
    psplups:ptr_new(/allocate_heap), $
    cisplups:ptr_new(/allocate_heap), $
    ioniz:ptr_new(/allocate_heap), $
    recomb:ptr_new(/allocate_heap), $
    cilvl:ptr_new(/allocate_heap), $
    reclvl:ptr_new(/allocate_heap), $
    population:ptr_new(/allocate_heap), $
    emiss:ptr_new(/allocate_heap), $
    abund:ptr_new(/allocate_heap), $
    ioneq:ptr_new(/allocate_heap) $
    }
end