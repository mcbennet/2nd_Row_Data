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
3
  1  ,  1.91050800 ,   3.00000000 
  3  ,  1.69821900 ,   8.93428800 
  2  ,  2.97809600 ,  -9.40279400 
1
  2  ,  9.96813000 ,  17.57180200 
1
  2  ,  2.72141500 ,  12.14866500 
include,/remote/cmelton/data/corrECP/2nd_row/aug-cc-pVDZ.basis
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
