***,Calculation for Be atom, singlet and triplet
memory,512,m


! s  0 x^2,y^2,z^2
! p  1 x
! p -1 y
! d -2 xy
! p  0 z
! d  1 xz
! d -1 yz
! f -2 xyz

print,orbitals=2

gparam,NOBUFF
gthresh,twoint=1.e-12
geometry={
1
P
P  0.0 0.0 0.0
}

basis={
ecp,P,10,2,0;
5;
1,0.0,-13.0;
1,0.0,5.0;
1,11.695975118527,13.0;
3,12.232718710797,152.047676540851;
2,12.489786639170001,-103.251368355488;
2;
2,54.280555383252,35.012435877038996;
2,34.449772213685996,155.677422007114;
1;
2,1.0,0.0;

include,aug-cc-pCVTZ.basis
}

include,ppstates_sc.proc

do i=0,6
    if (i.eq.0) then
        Nn
    else if (i.eq.1) then
        N1
    else if (i.eq.2) then
        N2
    else if (i.eq.3) then
        N3
    else if (i.eq.4) then
        N4
    else if (i.eq.5) then
        N5
    else if (i.eq.6) then
        N6
    else if (i.eq.7) then
        N7
    else if (i.eq.8) then
        N8
    endif
    scf(i+1)=energy
    _CC_NORM_MAX=2.0
    {uccsd(t),maxit=100;core}
    ccsd(i+1)=energy
enddo

table,scf,ccsd
save
type,csv
