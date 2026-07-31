load "write_forms.m";


procedure RecomputeSpecific(d,label,p)
	F:=QuadFld(d);
    ZF:=MaximalOrder(F);
    heckeBd := 100;
    filename := "data/nonEis_d" cat Sprint(d) cat "_recomputed.csv";
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



d := 11;
levels := [ < "125.4",  2 >,
< "269.1",  2 >,
< "284.1",  2 >,
< "335.2",  2 >,
< "355.1",  2 >,
< "397.1",  2 >,
< "445.2",  19 >,
< "515.4",  2 >,
< "555.5",  3 >,
< "564.1",  2 >,
< "583.2",  3 >,
< "636.3",  2 >,
< "639.6",  3 >,
< "669.3",  2 >,
< "675.12",  3 >,
< "675.3",  3 >,
< "716.2",  3 >,
< "729.2",  7 >,
< "764.1",  2 >,
< "775.2",  3 >,
< "775.2",  23 >,
< "775.2",  29 >,
< "775.6",  3 >,
< "795.4",  3 >,
< "807.2",  3 >,
< "837.5",  2 >,
< "841.1",  2 >,
< "892.1",  2 >,
< "925.1",  3 >,
< "925.3",  3 >,
< "925.2",  17 >,
< "961.2",  2 >,
< "972.6",  2 > ];

for l in levels do 
RecomputeSpecific(d,l[1],l[2]);
end for;




/*

1: [ <"666.2", 3>, <"941.2", 2>, <"797.2", 2>, <"754.1",3>, <"970.3",11> ]
2: [ <"328.1", 5>, <"267.4", 3>, <"347.1", 2>, <"921.1", 3>, <"801.4", 7>, <"841.1", 113>, <"625.1", 3>, <"801.1", 17>, <"801.1", 23>, <"722.3", 2>, <"683.2", 7>, <"979.1", 17>, <"979.1", 653> ]
3: [ <"919.2",2>]
7: [ <"961.1", 2>, <"709.1", 3>, <"763.2", 3>, <"599.2", 2>, <"648.2", 3>, <"737.1", 3>, <"694.1", 3> ]
11: [ < "125.4",  2 >,
< "269.1",  2 >,
< "284.1",  2 >,
< "335.2",  2 >,
< "355.1",  2 >,
< "397.1",  2 >,
< "445.2",  19 >,
< "515.4",  2 >,
< "555.5",  3 >,
< "564.1",  2 >,
< "583.2",  3 >,
< "636.3",  2 >,
< "639.6",  3 >,
< "669.3",  2 >,
< "675.12",  3 >,
< "675.3",  3 >,
< "716.2",  3 >,
< "729.2",  7 >,
< "764.1",  2 >,
< "775.2",  3 >,
< "775.2",  23 >,
< "775.2",  29 >,
< "775.6",  3 >,
< "795.4",  3 >,
< "807.2",  3 >,
< "837.5",  2 >,
< "841.1",  2 >,
< "892.1",  2 >,
< "925.1",  3 >,
< "925.3",  3 >,
< "925.2",  17 >,
< "961.2",  2 >,
< "972.6",  2 > ]

*/