%% Balance force on cargo pulled by np plus and nm minus motors

function Fc = Fc(np,nm)

global Fsp Fsm vbm vfm vbp vfp

if np*Fsp>nm*Fsm
    lambda = nm*Fsm*vfp/(nm*Fsm*vfp+np*Fsp*vbm);
    Fc = lambda*np*Fsp+(1-lambda)*nm*Fsm;
    
else
    
    lambda = nm*Fsm*vbp/(nm*Fsm*vbp+np*Fsp*vfm);
    Fc = lambda*np*Fsp+(1-lambda)*nm*Fsm;
    
end

end