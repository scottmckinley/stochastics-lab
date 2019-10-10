%% Set up vectors that go on the main or off-diagonals in the transition
% probability matrix

function [v1,v2,v3,v4] = vectors()

v1 = [epsm(0,1)/(epsm(0,1)+pip(2,0)+pim(2,1)) epsm(1,1)/(epsm(1,1)+epsp(1,1)+pip(2,1)+pim(2,1))...
    epsm(2,1)/(epsm(2,1)+epsp(2,1)+pim(2,1)) epsm(0,2)/(epsm(0,2)+pip(2,0))...
    epsm(1,2)/(epsm(1,2)+epsp(1,2)+pip(2,1)) epsm(2,2)/(epsm(2,2)+epsp(2,2))
    ];

v2 = [pim(2,0)/(pim(2,0)+pip(2,0))  pim(2,0)/(epsp(1,0)+pim(2,0)+pip(2,1))...
    pim(2,0)/(pim(2,0)+epsp(2,0))  pim(2,1)/(epsm(0,1)+pim(2,1)+pip(2,0))...
    pim(2,1)/(epsm(1,1)+epsp(1,1)+pim(2,1)+pip(2,1))  pim(2,1)/(epsm(2,1)+epsp(2,1)+pim(2,1))
    ];

v3 = [epsp(1,0)/(epsp(1,0)+pip(2,1)+pim(2,0)) epsp(2,0)/(epsp(2,0)+pim(2,0)) 0 ...
    epsp(1,1)/(epsp(1,1)+epsm(1,1)+pip(2,1)+pim(2,1)) epsp(2,1)/(epsp(2,1)+epsm(2,1)+pim(2,1)) 0 ...
    epsp(1,2)/(epsp(1,2)+epsm(1,2)+pip(2,1)) epsp(2,2)/(epsp(2,2)+epsm(2,2))
    ];

v4 = [pip(2,0)/(pip(2,0)+pim(2,0)) pip(2,1)/(epsp(1,0)+pip(2,1)+pim(2,0)) 0 ...
    pip(2,0)/(epsm(0,1)+pip(2,0)+pim(2,1)) pip(2,1)/(epsp(1,1)+epsm(1,1)+pip(2,1)+pim(2,1)) 0 ...
    pip(2,0)/(epsm(0,2)+pip(2,0)) pip(2,1)/(epsp(1,2)+epsm(1,2)+pip(2,1))
    ];


end