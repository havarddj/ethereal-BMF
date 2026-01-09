function LoadForm(d, ideal, evals)
    
end function;


function FindLifts(f : normBd := 100)
    LR_ideals := LevelRaiseFactors(f, normBd);
    B := Parent(f);
    // TODO: also look for primes dividing level!
    level1 := Level(B);
    Wp := Weight(B);
    W0 := BianchiWeight(K, 0,0);


    for I in [I : I in LR_ideals | IsCoprime(Level(B), I)] do
	level2 := I*level1;

	B0 := BianchiCohomologySpace(level2, W0);
	Bp := BianchiCohomologySpace(level2, Wp);
	
	reduction := [B2`down(ReductionModPMap(B0,Bp)(Inverse(B0`down)(B0`forms.i))) : i in [1..Dimension(B0`forms)]];
	oldspace := DegenerateSubspace(f,Bp);
        inter := sub<Bp`forms | reduction> meet oldspace;
	if Dimension(inter) gt 0 then
	    printf "Found characteristic zero lift of level %o\n", LMFDBLabel(level2);
	    return inter;
	end if;

    end for;
    return "None found"
end function;
