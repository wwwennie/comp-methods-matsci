% find thermodynamic capacities, i.e. fluctuations
function[vau,vak,vae,vat,vap] = variance(un,kn,en,tn,pn,xmin,xmax)

[au,ak,ae,at,ap] = averaging(un,kn,en,tn,pn,xmin,xmax);
[ausq,aksq,aesq,atsq,apsq] = averaging(un.^2,kn.^2,en.^2,tn.^2,pn.^2,xmin,xmax);

vau = au^2 - ausq; 
vak = ak^2 - aksq;
vae = ae^2 - aesq;
vat = at^2 - atsq;
vap = ap^2 - apsq;


