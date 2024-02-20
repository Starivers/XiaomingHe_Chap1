function [w,pt] = gaussValues_1d(n)    
    switch (n)
        case 1
            w=2; pt=0;
        case 2
            w=[1,1]; pt=[-1/sqrt(3), 1/sqrt(3)];
        case 3
            w=[5/9, 8/9, 5/9]; pt=[-sqrt(3/5), 0, sqrt(3/5)];
        case 4
            w=[(18+sqrt(30))/36, (18+sqrt(30))/36, (18-sqrt(30))/36, (18-sqrt(30))/36]; 
            pt=[-sqrt(3/7-2/7*(sqrt(6/5))), sqrt(3/7-2/7*(sqrt(6/5))),...
                -sqrt(3/7+2/7*(sqrt(6/5))), sqrt(3/7+2/7*(sqrt(6/5)))];
        case 5
            w=[(322+13*sqrt(70))/900, (322+13*sqrt(70))/900, 128/225,...
                (322-13*sqrt(70))/900,(322-13*sqrt(70))/900]; 
            pt=[-1/3*sqrt(5-2*sqrt(10/7)), 1/3*sqrt(5-2*sqrt(10/7)),0,...
                -1/3*sqrt(5+2*sqrt(10/7)), 1/3*sqrt(5+2*sqrt(10/7))];
        otherwise
            error('No data are defined for this value');
    end
end