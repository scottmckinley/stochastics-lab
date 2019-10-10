%% Cargo velocity in a state with cargo driven by np plus and nm minus motors

function vc = vc(np,nm) 

global Fsp Fsm vbm vfm vbp vfp Fext

if np*Fsp>nm*Fsm
    vc = (np*Fsp-nm*Fsm-Fext)/((np*Fsp/vfp)+(nm*Fsm/vbm));
    
else
    vc = (np*Fsp-nm*Fsm-Fext)/((np*Fsp/vbp)+(nm*Fsm/vfm));
    
end