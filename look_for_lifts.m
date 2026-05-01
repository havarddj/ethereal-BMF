function LoadForm(d,p,lvl_label, evals : HeckeBd := 30)
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
	end for;
	form := f;
	break;
    end for;
    
    if Type(form) eq Type(0) then
	return "Eigenform not found";
    else
	return form;
    end if;
end function;


function FindLifts(f : lvlBd := 100)
    LR_nums := LevelRaiseFactors(f, lvlBd);
    // LR_nums := LevelRaisingPrimes(f);
    if #LR_nums eq 0 then
	print "No level raising primes available";
	return 0;
    end if;
    
    B := Parent(f);
    level1 := Level(B);
    NormBd := GetHeckeBound(B);
    K := NumberField(B);
    ZK := Integers(K);
    Wp := Weight(B);
    W0 := BianchiWeight(K, 0,0);
    p := Characteristic(B);

    // TODO: also look for primes dividing level!
    for I in [I : I in LR_nums | IsCoprime(Level(B), I)] do
	printf "Looking for lifts of level %o \n", LMFDBLabel(I);
	level2 := I*level1;

	B0 := BianchiCohomologySpace(level2, W0);
	Bp := BianchiCohomologySpace(level2, Wp);
	red_map := ReductionModPMap(B0,Bp);
	reduction := [Bp`down(red_map((Inverse(B0`down)(B0`forms.i)))) : i in [1..Dimension(B0`forms)]];
	oldspace := DegenerateSubspace(f, Bp);
        inter := sub<Bp`forms | reduction> meet oldspace;

	if Dimension(inter) gt 0 then
	    print "Found intersection of oldspace and reduction mod p in level", LMFDBLabel(level2);
	end if;
	// vprint User1: "Computing builtin Hecke eigenforms";

	C := BianchiCuspForms(K, level2);
	for j -> F in NewformDecomposition(NewSubspace(C)) do
	    print "Testing eigenform", j;
	    wrong_ctr := 0;
	    is_wrong := false;
	    primes := GoodHeckePrimes(B, NormBd);
	    E := Parent(HeckeEigenvalue(Eigenform(F), primes[1]));
	    print "Eigenvalues of char 0 lift lie in", E, "of discriminant", Discriminant(Integers(E));
	    p_primes := [m[1] : m in Factorization(p*Integers(E))];
	    for pp in Factorization(p*Integers(E)) do
		_, phi := ResidueClassField(pp[1]);
		for qq in primes do 
		    if Eigenvalue(f,qq) ne phi(HeckeEigenvalue(Eigenform(F), qq)) then
			// this is not the right one, move on to another level raising prime
			// but allow a couple of mismatches just in case
			wrong_ctr +:= 1;
			if wrong_ctr gt 2 then
			    is_wrong := true;
			    break;
			end if;
		    end if;
		    printf "%o \t %o \t %o\n", LMFDBLabel(qq), phi(HeckeEigenvalue(Eigenform(F), qq)), Eigenvalue(f,qq);
		end for;
		if not is_wrong then 
		    print "Found correct form!";
		    return F;
		end if;
	    end for;
	end for;
	print "---\n";
	// return inter;

    end for;
    return "None found";
end function;



/*
Attempts to find lifts of Bianchi cohomology class with irrational coefficients which is given as a string in the format

127.1;7;[ 1, 4 ];[ 1, 2 ];[ 5, 2 ];[ 5, 2 ];[ 0, 3 ];x^2 + 6*x + 3

where 
- the first entry is the level
- the second is the prime p,
- the last entry is the characteristic polynomial of the extension in which the Hecke eigenvalues lie
- the rest of the entries correspond to Hecke eigenvalues
*/

function BatchFindIrrationalLifts(d : lvlBd := 100)
    F := QuadFld(d);
    filename := "data/nonEis_d" cat Sprint(d) cat ".csv";
    lines := Split(Read(filename), "\n");
    header := Split(lines[1], ";");
    primeLabels := header[3..#header];
    primeList := [LMFDBIdeal(F,label) : label in primeLabels];
    for line in lines[2..#lines] do
	entries := Split(line, ";");
	if "x" notin entries[#entries] or "[" notin entries[3] then
	    continue;
	end if;

	lvl := LMFDBIdeal(F, entries[1]);
	p := StringToInteger(entries[2]);
	B := BianchiCohomologySpace(lvl, BianchiWeight(F,0,0 : char :=p));
	f := ReadClass(line,B, primeList);
	print "Finding lifts; this may take time";
	print FindLifts(f : lvlBd := lvlBd);
    end for;
    return 0;
end function;

function FindIrrationalLifts(f)
    LRPrimes := LevelRaisingPrimes(f);
end function;
