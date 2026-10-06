// TODO: when they fix type bug in magma, change signature of this to ModFrmBianchiElt!

intrinsic ModPEigenvalues(f::ModFrmHilElt, p::RngIntElt, prime_list::SeqEnum) -> .
{Return list of all mod p reductions of eigenvalues builtin eigenform f.
Output is of the form [pp, A] where pp is a prime above p in the Hecke field of f,
and A is an associative array encoding the eigenvalues mod pp}
    E := Parent(HeckeEigenvalue(f, prime_list[1]));
    output := [* *];
    for pp in Factorization(p*Integers(E)) do
	_, phi := ResidueClassField(pp[1]);
	A := AssociativeArray();
	for qq in prime_list do
	    A[qq] := phi(HeckeEigenvalue(f, qq));
	end for;
	Append(~output,[* pp[1], A *]);
    end for;
    return output;
end intrinsic;

// TODO: when they fix type bug in magma, change signature of this to ModFrmBianchiElt!
intrinsic PrintModPEigensystems(f::ModFrmHilElt, p::RngIntElt)
    {Print mod p eigenvalue systems associated to Bianchi eigenform f}
    primes := SortByLMFDBLabel(PrimesUpTo(100, Parent(f)`Field));
    primes := [qq : qq in primes | IsCoprime(qq, Level(Parent(f)))];
    evs := ModPEigenvalues(f, p, primes);
    for ev in evs do
	print "Mod", LMFDBLabel(ev[1]), "eigenvalues:", [ev[2][qq] : qq in primes];
    end for;
end intrinsic;

intrinsic PrintModPEigensystems(C::ModFrmBianchi, p::RngIntElt)
{Print all mod p eigenvalue systems associated to newforms in Bianchi space C}
    for f in NewformDecomposition(NewSubspace(C)) do
	PrintModPEigensystems(Eigenform(f),p); 
    end for;
end intrinsic;
