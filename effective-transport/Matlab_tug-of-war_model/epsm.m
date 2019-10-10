%% Unbinding rate for minus motors for np plus and nm minus motors

function epsm = epsm(np,nm) 

global eps0m Fdm 

Fcargo = Fc(np,nm);

epsm    = (nm*eps0m*exp(Fcargo/(nm*Fdm)));

end

