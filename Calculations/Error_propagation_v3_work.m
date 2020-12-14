syms theta kx ky e1x e1y e3x e3z real

syms ssu_theta ssu_kx ssu_ky real
syms ssu_e1x ssu_e1y real
syms ssu_e3x ssu_e3z real

k = [kx, ky, sqrt(1-kx^2-ky^2)];
e1 = [e1x, e1y, sqrt(1-e1x^2-e1y^2)];
assume(e1, 'real')
e3 = [e3x, sqrt(1-e3x^2-e3z^2), e3z];
assume(e3, 'real')
e2 = cross(e3,e1)/sqrt(sum(cross(e3,e1).^2,2));


% Theta1 -----------------------------------------------------------------
ftheta1 = dot(cross(e2,e3),k)*theta;
gtheta1 = dot(cross(e1,e2),e3);
%
dftheta1dtheta = gradient(ftheta1,theta);
dftheta1dkx = gradient(ftheta1,kx);
dftheta1dky = gradient(ftheta1,ky);
dftheta1de1x = gradient(ftheta1,e1x);
dftheta1de1y = gradient(ftheta1,e1y); 
dftheta1de3x = gradient(ftheta1,e3x);
dftheta1de3z = gradient(ftheta1,e3z);
%
dgtheta1dtheta = gradient(gtheta1,theta); % 0
dgtheta1dkx = gradient(gtheta1,kx); % 0
dgtheta1dky = gradient(gtheta1,ky); % 0
dgtheta1de1x = gradient(gtheta1,e1x);
dgtheta1de1y = gradient(gtheta1,e1y); 
dgtheta1de3x = gradient(gtheta1,e3x);
dgtheta1de3z = gradient(gtheta1,e3z);

ssu_theta1 = (1/gtheta1)^4*...
        ((gtheta1*dftheta1dtheta - ftheta1*dgtheta1dtheta)^2*ssu_theta + ...
        (gtheta1*dftheta1dkx - ftheta1*dgtheta1dkx)^2*ssu_kx + ...
        (gtheta1*dftheta1dky - ftheta1*dgtheta1dky)^2*ssu_ky + ...
        (gtheta1*dftheta1de1x - ftheta1*dgtheta1de1x)^2*ssu_e1x + ...
        (gtheta1*dftheta1de1y - ftheta1*dgtheta1de1y)^2*ssu_e1y + ...
        (gtheta1*dftheta1de3x - ftheta1*dgtheta1de3x)^2*ssu_e3x + ...
        (gtheta1*dftheta1de3z - ftheta1*dgtheta1de3z)^2*ssu_e3z);

