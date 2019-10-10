%% Unbinding rate for plus motors for np plus and nm minus motors

function epsp = epsp(np,nm) 

global eps0p Fdp 

Fcargo = Fc(np,nm); 

epsp   = (np*eps0p*exp(Fcargo/(np*Fdp)));

end