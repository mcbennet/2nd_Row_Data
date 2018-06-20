***,Calculation for Be atom, singlet and triplet
memory,512,m
gthresh,twoint=1.0E-12
geometry={
1
Mg
Mg  0.0 0.0 0.0
}

basis={
ecp,Mg,10,2,0
5
 1, 0.0,  2
 1, 0.0, -10.0
 1  ,  11.99962000  ,   10.00000000 
 3  ,  12.00778300  ,  119.99620000    
 2  ,  11.98802200  ,  -76.19471400 
2
 2  ,  90.16677400  ,   54.58518300 
 2  ,  23.95808900  ,  117.10260900 
1 
 2, 1.0, 0.0
include,aug-cc-pCV5Z.basis
}

include,states.proc

do i=1,3
if (i.eq.1) then
    GS
else if (i.eq.2) then
    IP1
else if (i.eq.3) then
    IP2
endif
scf(i)=energy
_CC_NORM_MAX=2.0
{uccsd(t),maxit=400;core}
ccsd(i)=energy
enddo

table,scf,ccsd
save
type,csv
