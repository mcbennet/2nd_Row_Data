***,Calculation for Be atom, singlet and triplet
memory,512,m
gthresh,twoint=1.0E-12
geometry={
1
Si
Si  0.0 0.0 0.0
}

basis={
ecp,Si,10,2,0
5
  1, 0, 4
  1, 0 ,-12
  1  ,  16.14603300  ,   12.00000000 
  3  ,  17.07755600  ,  193.75239600 
  2  ,  16.42528400  , -104.82393500 
1
  2  ,  22.54933800  ,   96.90843500 
1
  2 , 1.0, 0.0
include,/remote/cmelton/data/corrECP/2nd_row/aug-cc-pCVDZ.basis
}

include,states.proc

do i=1,6
if (i.eq.1) then
    GS
else if (i.eq.2) then
    IP1
else if (i.eq.3) then
    IP2
else if (i.eq.4) then
    IP3
else if (i.eq.5) then
    IP4
else if (i.eq.6) then
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
