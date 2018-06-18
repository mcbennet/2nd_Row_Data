***,Calculation for Be atom, singlet and triplet
memory,512,m
gthresh,twoint=1.0E-12
geometry={
1
Al
Al  0.0 0.0 0.0
}

basis={
ecp,Al,10,2,0
5
  1, 0.0, 3
  1, 0.0,-11
  1  ,  10.99788500  ,   11.00000000
  3  ,  11.22045000  ,  120.97673500 
  2  ,  11.08325000  ,  -80.39326600 
2
  2  ,  80.97577200  ,   25.00635800
  2  ,  24.39824000  ,  112.36079900
1
  2, 1.0, 0.0
include,/remote/cmelton/data/corrECP/2nd_row/aug-cc-pCVTZ.basis
}

include,states.proc

do i=1,5
if (i.eq.1) then
    GS
else if (i.eq.2) then
    IP1
else if (i.eq.3) then
    IP2
else if (i.eq.4) then
    IP3
else if (i.eq.5) then
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
