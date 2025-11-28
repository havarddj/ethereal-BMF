AttachSpec("../spec");
load "homology_torsion.m";


// Compute nonliftable eigenvalue systems with level at least levelLowerBd
function ComputeNonliftable(d : bd := 100, levelLowerBd := 0)
    print "Computing H_1 to find nonliftable forms";
    levels_and_primes := ComputeLevelsAndPrimes(d);
    print "Finished computing H_1";
    
    F:=QuadFld(d);
    ZF:=MaximalOrder(F);
    filename := "data/nonEis_d" cat Sprint(d) cat ".csv";
    // bd := 100;
    labels := [LMFDBLabel(pp) : pp in SortByLMFDBLabel(PrimesUpTo(bd, F)) ];

    // need to mess around a bit to make sure we don't print the final semicolon
    fprintf filename, "level;p;%o%o\n", &cat[l cat ";" : l in labels | l ne labels[#labels]], labels[#labels];

    for L in levels_and_primes do
	level_gen := L[1];
	level := (ZF!level_gen)*ZF;
	// if startBd is non-zero
	if Norm(level) lt levelLowerBd then
	    continue;
	end if;
	wt1 := BianchiWeight(F, 0, 0);
	B1 := BianchiCohomologySpace(level,wt1);
	for p in L[2] do
	    wt2:=BianchiWeight(F, 0,0 : char:= Integers()!p);
	    B2 := BianchiCohomologySpace(level,wt2);
	    SetHeckeBound(B2, bd);
	    BG := GenuineSubspace(B1,B2);
	    for f in Eigenforms(B2) do
		if not IsEisenstein(f) and f in BG and f`eigenspaceDim eq 1 then
		    WriteClass(f, filename : labels := labels);
		    print "Wrote class to file!";
		end if;
	    end for;
	end for;
    end for;
    return "";
end function;
