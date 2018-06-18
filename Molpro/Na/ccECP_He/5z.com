***,Calculation for Be atom, singlet and triplet
memory,512,m
gthresh,twoint=1.0E-12

geometry={
1
Na
Na  0.0 0.0 0.0
}

basis={
ecp,Na,10,2,0
5
  1,0, 1
  1,0,-9
  1  ,  11.99983500  ,    9.00000000 
  3  ,  11.94985400  ,  107.99851500 
  2  ,  12.04004600  ,  -69.25461900 
1
  2  ,  29.94712200  ,  221.73468900 
1
  1  , 1.0 , 0
include,/remote/cmelton/data/corrECP/2nd_row/aug-cc-pCV5Z.basis
}

include,states.proc

do i=1,3
if (i.eq.1) then
    GS
else if (i.eq.2) then
    IP
else if (i.eq.3) then 
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
