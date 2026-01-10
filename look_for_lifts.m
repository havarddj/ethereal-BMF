function LoadForm(d,p,lvl_label, evals : HeckeBd := 100)
    // Helper function to load form from value of d, p, the level label, and eigenvalues evals
    F := QuadFld(d);
    lvl := LMFDBIdeal(F, lvl_label);
    W := BianchiWeight(F, 0,0 : char := p);
    B := BianchiCohomologySpace(lvl, W);
    form := 0;
    good_primes := GoodHeckePrimes(B, HeckeBd);
    SetHeckeBound(B, HeckeBd);
    len := Min(#good_primes, #evals);
    // if Dimension(B) eq 1 then
    // 	form := Basis(B)[1];
    // 	form`isEigenform := true;
    // 	return form;
    // end if;
    for f in Eigenforms(B) do
	found := false;
	for i -> ev in evals[1..len] do
	    if ev ne 0 and ev ne EigenvalueList(f)[i] then 
		break;
	    end if;
	    form := f;
	    break;
	end for;
    end for;
    
    if Type(form) eq Type(0) then
	return "Eigenform not found";
    else
	return form;
    end if;
    
end function;


function FindLifts(f : lvlBd := 100)
    // LR_nums := LevelRaiseFactors(f, lvlBd);
    LR_nums := LevelRaisingPrimes(f);
    B := Parent(f);
    level1 := Level(B);
    NormBd := GetHeckeBound(B);
    K := NumberField(B);
    ZK := Integers(K);
    Wp := Weight(B);
    W0 := BianchiWeight(K, 0,0);

    // TODO: also look for primes dividing level!
    for I in [I : I in LR_nums | IsCoprime(Level(B), I)] do
	printf "Looking for lifts of level %o \n", LMFDBLabel(I);
	level2 := I*level1;

	B0 := BianchiCohomologySpace(level2, W0);
	Bp := BianchiCohomologySpace(level2, Wp);
	red_map := ReductionModPMap(B0,Bp);
	reduction := [Bp`down(red_map((Inverse(B0`down)(B0`forms.i)))) : i in [1..Dimension(B0`forms)]];
	oldspace := DegenerateSubspace(f,Bp);
        inter := sub<Bp`forms | reduction> meet oldspace;

	if Dimension(inter) gt 0 then
	    printf "Found characteristic zero lift of level %o\n", LMFDBLabel(level2);
	    print "Computing builtin Hecke eigenforms";

	    C := BianchiCuspForms(K,level2);
	    for f in NewformDecomposition(NewSubspace(C)) do
		primes := GoodHeckePrimes(Bp, NormBd);
		char0_evals := [HeckeEigenvalue(Eigenform(f), pp) : pp in primes];
		E := Parent(char0_evals[1]);
		print "Eigenvalues of char 0 lift lie in", E, "of discriminant", Discriminant(Integers(E));
		p_primes := [m[1] : m in Factorization(p*Integers(E))];
		for pp in Factorization(p*Integers(E)) do
		    _,phi := ResidueClassField(pp[1]);
		    print [phi(app) : app in char0_evals];
		end for;
		print "";
	    end for;
	    print "---\n";
	    // return inter;
	end if;

    end for;
    return "None found";
end function;
