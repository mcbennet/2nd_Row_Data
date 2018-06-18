***,Calculation for Be atom, singlet and triplet
memory,512,m
gthresh,twoint=1.0E-12
geometry={
1
Al
Al  0.0 0.0 0.0
}

basis={
ecp,Al,10,2,0;
3;
1,2.978096,3.0;
3,9.96813,8.934288;
2,2.721415,-9.402794;
1;
2,1.910508,17.571802;
1;
2,1.698219,12.148665;
include,/remote/cmelton/data/corrECP/2nd_row/aug-cc-pVTZ.basis
}

include,states.proc

do i=1,4
if (i.eq.1) then
    GS
else if (i.eq.2) then
    IP1
else if (i.eq.3) then
    IP2
else if (i.eq.4) then
    EA
endif
scf(i)=energy
_CC_NORM_MAX=2.0
{uccsd(t),maxit=400;core}
ccsd(i)=energy
enddo

table,scf,ccsd
save
type,csv
