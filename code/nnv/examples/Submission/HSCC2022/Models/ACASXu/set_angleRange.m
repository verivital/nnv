function ang = set_angleRange(angle)
% Compute angle in interval [-pi,pi]. For the discontinuity pi/-pi it
% choses always pi.
%ang = angle - ceil(angle/(2*pi) - 0.5)*2*pi;
ang = mod(angle+pi,2*pi)-pi;
end