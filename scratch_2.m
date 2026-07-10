load "write_forms.m";


procedure RecomputeSpecific(d,label,p)
	F:=QuadFld(d);
    ZF:=MaximalOrder(F);
    heckeBd := 100;
    filename := "data/nonEis_d" cat Sprint(d) cat ".csv";
    labels := [LMFDBLabel(pp) : pp in SortByLMFDBLabel(PrimesUpTo(heckeBd, F)) ];
    // bd := 100;
    level := LMFDBIdeal(F,label);
	wt1 := BianchiWeight(F, 0, 0);
	B1 := BianchiCohomologySpace(level, wt1);
	try 
		wt2 := BianchiWeight(F, 0,0 : char:= Integers()!p);
		B2 := BianchiCohomologySpace(level,wt2);
		SetHeckeBound(B2, heckeBd);
		//BG := EtherealSubspace(B1,B2);
		//print BG;
		print B2;
		for f in EigenformGaloisRepresentatives(B2) do
		    if not IsEisenstein(f) and HasEtherealEigenvalues(B1,f) and f`eigenspaceDim eq 1 then
			WriteClass(f, filename : labels := labels);
			print "Wrote class to file!\n";
		    end if;
		end for;
	    catch err;
		print "error:", err;
		fprintf filename, "%o; %o; ERROR\n", LMFDBLabel(level), p;
	end try;
end procedure;



d := 7;
levels := [
<"648.2",3>,
<"694.1",3>,
<"709.1",3>,
<"737.1",3>,
<"763.2",3>,
<"961.1",2>
];


for l in levels do 
RecomputeSpecific(d,l[1],l[2]);
end for;

