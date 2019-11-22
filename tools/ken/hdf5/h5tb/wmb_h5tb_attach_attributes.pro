

; wmb_h5tb_attach_attributes
; 
; Purpose: Private function that creates the conforming table attributes;
;          Used by wmb_h5tb_combine_tables.
; 


pro wmb_h5tb_attach_attributes, table_title, $
                                loc_id, $
                                dset_name, $
                                nfields, $
                                tid

    compile_opt idl2, strictarrsubs
    
    const_TABLE_CLASS = 'TABLE'
    const_TABLE_VER = 3.0
    const_HLTB_MAX_FIELD_LEN = 255
    
    ; attach the CLASS, VERSION, and TITLE attributes
    
    wmb_h5lt_set_attribute_string, loc_id, dset_name, 'CLASS', const_TABLE_CLASS
    wmb_h5lt_set_attribute_string, loc_id, dset_name, 'VERSION', const_TABLE_VER
    wmb_h5lt_set_attribute_string, loc_id, dset_name, 'TITLE', table_title
    
    
    ; attach the FIELD_ name attribute
    
    for i = 0, nfields-1 do begin
    
        member_name = h5t_get_member_name(tid, i)
        tmpstr = 'FIELD_' + strtrim(string(i),2) + '_NAME'
        wmb_h5lt_set_attribute_string, loc_id, dset_name, tmpstr, member_name

    endfor

end

