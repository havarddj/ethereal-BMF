
// Helper function to load form from value of d, p, the level label, and eigenvalues evals
// WARNING: this recomputes eigenforms, so it's expensive - it's better to use something else
function LoadForm(d,p,lvl_label, evals : HeckeBd := 30)
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

// turns a line from d*_liftable.csv into its raw data so we can feed it into VerifyLMFDBLift
function UnpackLiftableString(lmfdbstr)
	ss := Split(lmfdbstr,";");
	levelLabel := ss[1];
	p := StringToInteger(ss[2]);
	eigs := [StringToInteger(u) : u in Split(ss[3][2..#ss[3]-1],",")];

	// little utility function to extract the level from a BMF label 
	function BMFLabelToLevel(foo)
		bar := foo[Index(foo,"-")+1..#foo-1];
		return bar[1..Index(bar,"-")-1];
	end function;

	liftlabels := [BMFLabelToLevel(u) : u in Split(ss[4][2..#ss[4]-1],",")];

	return levelLabel,p,eigs,liftlabels;
end function;

// turns a line from d*_liftable.csv into its BMFCohomClass
function LoadForm_liftable(d,str)
	K := QuadFld(d);
	levelLabel,p,eigs,liftLabels := UnpackLiftableString(str);
	level := LMFDBIdeal(K, levelLabel);
	W := BianchiWeight(K, 0, 0 : char := p);
	B := BianchiCohomologySpace(level, W);

	s := Split(Read("../data/nonEis_d" cat Sprint(d) cat ".csv"),"\n")[1];
	// this cuts out just the part of the string with the labels 
	pp := s[9..#s-14];
	PP := [LMFDBIdeal(K,u) : u in Split(pp,";")];

	// this feels quite silly, having to re-egineer the string.
	fstr := levelLabel cat ";" cat Sprint(p) cat ";" cat &cat[Sprint(u) cat ";" : u in eigs];
	f := ReadClass(fstr,B,PP);
	return f;
end function;


// turns a line from nonEis_d*.csv into its BMFCohomClass
function LoadForm_nonEis(d,str)
	K := QuadFld(d);
	ss := Split(str,";");

	// gathering prime labels 
	s := Split(Read("../data/nonEis_d" cat Sprint(d) cat ".csv"),"\n")[1];
	pp := s[9..#s-34];
	PP := [LMFDBIdeal(K,u) : u in Split(pp,";")];

	B := BianchiCohomologySpace(LMFDBIdeal(K,ss[1]),BianchiWeight(K,0,0 : char := StringToInteger(ss[2])));	
	return ReadClass(str,B,PP);
end function;


