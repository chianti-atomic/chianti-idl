pro chianti_ion::cleanup
    ptr_free,self.wgfa
    ptr_free,self.elvlc
    ptr_free,self.splups
    ptr_free,self.psplups
    ptr_free,self.population
end

pro chianti_ion::printit
    print, ' ions = ',self.ions, self.z, self.stage
    mywgfa = *self.wgfa
    mysplups = *self.splups
    for i=0,4 do begin
        print, ' l1, ecm = ',(*self.elvlc).l1[i],(*self.elvlc).ecm[i]
        print, ' lvl1, avalue = ',(*self.wgfa)[i].lvl1,(*self.wgfa)[i].aval
        print ,' lvl1  splups.de = ',(*self.splups)[i].lvl1,(*self.splups)[i].de
    endfor
    IF n_tags(*(self.psplups)) NE 0 THEN BEGIN
        npsplups = n_elements(*self.psplups)
        print, ' npsplups = ',npsplups
        for i=0,npsplups-1 do begin
            print ,' lvl1  psplups.de = ',(*self.psplups)[i].lvl1,(*self.psplups)[i].de
        endfor
    endif
    return
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

pro chianti_ion::populate,temperature,density
    dilute = 0
    radt = 0
    ionrec = {rec:0., ci:0., temp:0., lev_up_rec:0., lev_up_ci:0., status:-1, ioneq:0., temp_ioneq:0.}
;     IF n_tags(*(self.psplups)) NE 0 THEN begin
;         pstr = (*self.psplups)
;     else 
    avalue = (*self.wgfa).a_value
;
    input = {gname:self.ionS, $
    jj:(*self.elvlc).jj, $
    ecm:(*self.elvlc).ecm, $
    ecmth:(*self.elvlc).ecmth, $
    wvl:(*self.wgfa).wvl, $
    a_value:avalue, $
    splstr:(*self.splups), $
    pe_ratio:proton_dens(alog10(temperature)), $
    prot_struc:*self.psplups, $
    ionrec:ionrec $
    }
;     dilute:dilute, $
;     radtemp:radt, $
;
    self->pop_solverx,input, temperature, density, pop
;     for i=0,n_elements(pop)-1 do begin
;         print, i, pop[i]
;     endfor
    population = {temperature:temperature, population:pop, density:density, ion:self.ionS, input:input}
;     population = input
    self.population = ptr_new(population)
    return
end

function chianti_ion::getPopulation
    return, *self.population
end

function chianti_ion::init,ionS
    self.ionS = ionS
    convertname,ionS, z, stage
    self.z = z
    self.stage = stage
    print, ' ion, z, stage = ',ionS,z,stage
    print, ' ion, z, stage = ',self.ionS,self.z,self.stage
    zion2filename,self.z,self.stage,fname
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
        ptr_psplups = ptr_new(psplupsstr)
        if ptr_valid(ptr_psplups) then begin
            self.psplups = ptr_psplups
        endif else print, ' psplups pointer not valid'
        
    endif
    
    read_cireclvl,ionS, reclvl, /rec
    ptr_reclvl = ptr_new(reclvl)
    if ptr_valid(ptr_reclvl) then begin
        self.reclvl = ptr_reclvl
    endif else print, ' reclvl pointer not valid'
    
    read_cireclvl, ionS, cilvl, /ci
    ptr_cilvl = ptr_new(cilvl)
    if ptr_valid(ptr_cilvl) then begin
        self.cilvl = ptr_cilvl
    endif else print,' cilvl pointer not valid'
    
    return,1
end

pro chianti_ion__define
    struct = {chianti_ion, z:1L, ions:' ', $
    stage:0L, $
    elvlc:ptr_new(/allocate_heap), $
    wgfa:ptr_new(/allocate_heap), $
    splups:ptr_new(/allocate_heap), $
    psplups:ptr_new(/allocate_heap), $
    ioniz:ptr_new(/allocate_heap), $
    recomb:ptr_new(/allocate_heap), $
    cilvl:ptr_new(/allocate_heap), $
    reclvl:ptr_new(/allocate_heap), $
    population:ptr_new(/allocate_heap) $
    }
end