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
4;
1,5.69505216,5.0;
3,13.37645903,28.4752608;
2,4.91050582,-19.88716289;
2,9.90588392,15.0693342;
2;
2,15.72747985,6.16882039;
2,3.13007778,33.92271491;
2;
2,2.16392613,4.07177501;
2,2.81692052,14.80148965;

include,aug-cc-pCV5Z.basis
}

include,ppstates.proc

do i=1,6
    if (i.eq.1) then
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
    scf(i)=energy
    _CC_NORM_MAX=2.0
    {uccsd(t),maxit=100;core}
    ccsd(i)=energy
enddo

table,scf,ccsd
save
type,csv
