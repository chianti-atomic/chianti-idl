
;
; wmb_h5tb_read_table
; 
; Purpose: Reads a table. 
;
; Description: wmb_h5tb_read_table reads a table named table_name attached to 
;              the object specified by the identifier loc_id.  
;              
; Parameters:
; 
; loc_id
;     IN: Identifier of the file or group in which the table is located. 
; dset_name
;     IN: The name of the table. 
; databuffer
;     OUT: Output data buffer. 
;

pro wmb_h5tb_read_table, loc_id, $
                         dset_name, $
                         databuffer

    compile_opt idl2, strictarrsubs

    ; open the dataset
    
    did = h5d_open(loc_id, dset_name)
    
    ; read the table
    
    databuffer = h5d_read(did)

    ; close
    
    h5d_close, did
    
end

