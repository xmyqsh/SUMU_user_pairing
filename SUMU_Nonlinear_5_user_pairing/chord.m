function d=chord(v,q)
% chordal distance
t=v*v'-q*q';
d=sqrt(trace(t*t'))/sqrt(2);