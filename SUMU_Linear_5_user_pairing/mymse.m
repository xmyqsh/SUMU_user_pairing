function se=mymse(a,b)
a=sqrt(numel(a)).*a./sqrt(trace(a'*a));
b=sqrt(numel(a)).*b./sqrt(trace(a'*a));
e=a-b;
se=trace(e'*e)/numel(e);
se=10*log10(se);