AttachSpec("../spec");
load "homology_torsion.m";


// Compute nonliftable eigenvalue systems with level at least levelLowerBd
function ComputeNonliftable(d : heckeBd := 100, levelUpperBd := 1000, levelLowerBd := 0)
    print "Computing H_1 to find nonliftable forms";
    levels_and_primes := ComputeLevelsAndPrimes(d : lowerBound:= levelLowerBd, upperBound := levelUpperBd);
    print "Finished computing H_1";
    
    F:=QuadFld(d);
    ZF:=MaximalOrder(F);
    filename := "data/nonEis_d" cat Sprint(d) cat ".csv";
    // bd := 100;
    labels := [LMFDBLabel(pp) : pp in SortByLMFDBLabel(PrimesUpTo(heckeBd, F)) ];

    // need to mess around a bit to make sure we don't print the final semicolon
    fprintf filename, "level;p;%o;Coeff_minpoly\n", Join(labels, ";");

    for L in levels_and_primes do
	level_gen := L[1];
	level := (ZF!level_gen)*ZF;
	// if startBd is non-zero
	if Norm(level) lt levelLowerBd then
	    continue;
	end if;
	wt1 := BianchiWeight(F, 0, 0);
	B1 := BianchiCohomologySpace(level, wt1);
	for p in L[2] do
	    try 
		wt2 := BianchiWeight(F, 0,0 : char:= Integers()!p);
		B2 := BianchiCohomologySpace(level,wt2);
		SetHeckeBound(B2, heckeBd);
		BG := EtherealSubspace(B1,B2);
		print BG;
		for f in Eigenforms(B2) do
		    if not IsEisenstein(f) and f in ChangeRing(BG,BaseRing(f)) and f`eigenspaceDim eq 1 then
			WriteClass(f, filename : labels := labels);
			print "Wrote class to file!";
		    end if;
		end for;
	    catch err;
		print "error:", err;
		fprintf filename, "%o; %o; ERROR\n", LMFDBLabel(level), p;
	    end try;
	end for;
    end for;
    return "";
end function;

/*
Recompute non-liftable forms (useful after changing formatting)
and append them to the end of the file

So we have to order them manually, or using this bash command (on unix systems):

`sort -k1 -n -t";" -o filename filename`
(here:
-k1 sorts by first col,
-n is numerical,
-t";" sets delimiter,
-o filename sets output file (and the second filename specifies input file))

(It might be helpful to sort in-place: perl -i -ne 'print if ! $x{$_}++' filename)
The condition for recomputing is set manually in the function,
see the comment labeled "(*)". 
*/
function RecomputeIrrational(d : heckeBd := 100, lvlLowerBd := 0, lvlUpperBd := 100000)
    F := QuadFld(d);
    filename := "data/nonEis_d" cat Sprint(d) cat ".csv";
    lines := Split(Read(filename), "\n");
    header := Split(lines[1], ";");
    primeLabels := header[3..#header];
    primeList := [LMFDBIdeal(F,label) : label in primeLabels];
    recomputedPairs := [];
    for line in lines[2..#lines] do
	entries := Split(line, ";");
	// (*) current miscomputed ones:
	if not exists{e : e in entries[3..#entries-1] | "w" in e or "x" in e} then
	    continue;
	end if;
	level := LMFDBIdeal(F, entries[1]);
	if Norm(level) lt lvlLowerBd or Norm(level) gt lvlUpperBd then
	    continue;
	end if;
	print line;
	p := StringToInteger(entries[2]);

	// Keep track of levels and primes computed so we don't double count Hecke conjugates
	if [*level,p*] in recomputedPairs then
	    print "Already computed lifts for this system";
	    continue;
	end if;
	Append(~recomputedPairs, [*level,p*]);

	wt1 := BianchiWeight(F, 0, 0);
	B1 := BianchiCohomologySpace(level,wt1);
	    try 
		wt2 := BianchiWeight(F, 0,0 : char:= Integers()!p);
		B2 := BianchiCohomologySpace(level,wt2);
		SetHeckeBound(B2, heckeBd);
		print "Computing ethereal subspace";
	
		BG := EtherealSubspace(B1,B2);
		print BG;
		for f in Eigenforms(B2) do
		    
		    if f`eigenspaceDim gt 1 then
			print "Found eigenspace of dimension", f`eigenspaceDim;
		    end if;
		    if not IsEisenstein(f) and f in ChangeRing(BG, BaseRing(f)) and f`eigenspaceDim eq 1 then
			WriteClass(f, filename : labels := primeLabels);
			print "Wrote class to file!";
		    end if;
		end for;
	    catch err;
		print "error computing forms:", err;
		fprintf filename, "%o; %o; ERROR\n", LMFDBLabel(level), p;
	    end try;
    end for;
    return "";
end function;


// Look for instances of ethereal spaces where the level raise at p has increased multiplicity
function LookForMultiplicity(d : heckeBd := 100, levelUpperBd := 1000, levelLowerBd := 0)
    print "Computing H_1 to find ethereal characteristics.";
    levels_and_primes := ComputeLevelsAndPrimes(d : lowerBound:= levelLowerBd, upperBound := levelUpperBd);
    print "Finished computing ethereal characteristics.";
    
    F := QuadFld(d);
    ZF := MaximalOrder(F);

    labels := [LMFDBLabel(pp) : pp in SortByLMFDBLabel(PrimesUpTo(heckeBd, F)) ];

    for L in levels_and_primes do
	level := (ZF!L[1])*ZF;
	B0 := BianchiCohomologySpace(level, [0,0]);
	for p in L[2] do
	    if p notin [3,5] or not IsCoprime(level, p*ZF) then
		continue;
	    end if;
	    Wp := BianchiWeight(F, 0, 0 : char := Integers()!p);
	    Bp := BianchiCohomologySpace(level, Wp);
	    SetHeckeBound(Bp, heckeBd);
	    BEth := EtherealSubspace(B0, Bp);
	    print "Ethereal subspace:", BEth;
	    for f in Eigenforms(Bp) do
		if f in ChangeRing(BEth, BaseRing(f)) and IsCoprime(p*Integers(F), Level(f)) then
		    BLevelRaise := BianchiCohomologySpace(level*p, Wp);
		    SetHeckeBound(BLevelRaise, heckeBd);
		    if Dimension(BLevelRaise) in [1, Dimension(Bp)] then
			continue;
		    end if;

		    print "Computing eigenforms with p in the level; dimension is", Dimension(BLevelRaise);

		    for g in Eigenforms(BLevelRaise) do
			flag := g`eigenspaceDim gt 1;
			flag and:= forall{1 : pp in PrimesUpTo(heckeBd, F) | Eigenvalue(g,pp) eq Eigenvalue(f,pp)};
			if flag then
			    printf "Located multiplicity %o space of level %o mod %o\n", g`eigenspaceDim, LMFDBLabel(Level(g)), p;
			    return g;
			else
			    print "No multiplicity one";
			    _ := exists(pp){pp : pp in PrimesUpTo(heckeBd, F) | Eigenvalue(g,pp) ne Eigenvalue(f,pp)};
			    print "They differ at", LMFDBLabel(pp);
			end if;
		    end for;
		end if;
		// catch err
		//     print "Failed with error", err;
		//     print f;
		// end try;
	    end for;
	end for;
    end for;
    return "";
end function;

