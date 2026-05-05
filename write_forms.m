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
	B1 := BianchiCohomologySpace(level,wt1);
	for p in L[2] do
	    // try 
		wt2 := BianchiWeight(F, 0,0 : char:= Integers()!p);
		B2 := BianchiCohomologySpace(level,wt2);
		SetHeckeBound(B2, heckeBd);
		BG := GenuineSubspace(B1,B2);
		print BG;
		for f in Eigenforms(B2) do
		    if not IsEisenstein(f) and f in ChangeRing(BG,BaseRing(f)) and f`eigenspaceDim eq 1 then
			WriteClass(f, filename : labels := labels);
			print "Wrote class to file!";
		    end if;
		end for;
	    // catch err;
	    // 	fprintf filename, "%o; %o; ERROR\n", LMFDBLabel(level), p;
	    // end try;
	end for;
    end for;
    return "";
end function;

/*
Recompute non-liftable forms (useful after changing formatting)
and append them to the end of the file

So we have to order them manually, or using this bash command (on unix systems):

`sort -k1 -n -t; -u -o filename filename`
(here:
-k1 sorts by first col,
-n is numerical,
-t; sets delimiter,
-o filename sets output file (and the second filename specifies input file))

(It might be helpful to sort in-place: perl -i -ne 'print if ! $x{$_}++' filename)
The condition for recomputing is set manually in the function,
see the comment labeled "(*)". 
*/
function RecomputeIrrational(d : heckeBd := 100)
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
	p := StringToInteger(entries[2]);
	if [*level,p*] in recomputedPairs then
	    continue;
	end if;
	Append(~recomputedPairs,[*level,p*]);

	wt1 := BianchiWeight(F, 0, 0);
	B1 := BianchiCohomologySpace(level,wt1);
	    // try 
		wt2 := BianchiWeight(F, 0,0 : char:= Integers()!p);
		B2 := BianchiCohomologySpace(level,wt2);
		SetHeckeBound(B2, heckeBd);
		BG := GenuineSubspace(B1,B2);
		print BG;
		for f in Eigenforms(B2) do
		    if not IsEisenstein(f) and f in ChangeRing(BG,BaseRing(f)) and f`eigenspaceDim eq 1 then
			WriteClass(f, filename : labels := primeLabels);
			print "Wrote class to file!";
		    end if;
		end for;
	    // catch err;
	    // 	fprintf filename, "%o; %o; ERROR\n", LMFDBLabel(level), p;
	    // end try;
    end for;
    return "";
end function;
